// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniaPhotonProducer
//
// Implementation of Onia + Photon candidate producer
// Based on chi-analysis-miniaod by Alberto Sanchez-Hernandez

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/OniaPhotonProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include <TMath.h>
#include <TVector3.h>
#include <Math/VectorUtil.h>

OniaPhotonProducer::OniaPhotonProducer(const edm::ParameterSet& iConfig)
    : oniaToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("oniaSrc"))),
      conversionToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("conversionSrc"))),
      pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotSrc"))),
      ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      ttBuilder_(nullptr),
      minConvPt_(iConfig.getParameter<double>("minConvPt")),
      minConvRho_(iConfig.getParameter<double>("minConvRho")),
      rejectPi0_(iConfig.getParameter<bool>("rejectPi0")),
      minOniaPt_(iConfig.getParameter<double>("minOniaPt")),
      maxOniaAbsY_(iConfig.getParameter<double>("maxOniaAbsY")),
      requireTriggerMatch_(iConfig.getParameter<bool>("requireTriggerMatch")),
      dzMax_(iConfig.getParameter<double>("dzMax")),
      deltaMassRange_(iConfig.getParameter<std::vector<double>>("deltaMassRange")),
      minCandPt_(iConfig.getParameter<double>("minCandPt")),
      maxCandAbsY_(iConfig.getParameter<double>("maxCandAbsY")),
      doKinematicRefit_(iConfig.getParameter<bool>("doKinematicRefit")),
      constraintMass_(iConfig.getParameter<double>("constraintMass")),
      parentType_(iConfig.getParameter<std::string>("parentType")),
      nCandidates_(0),
      nDeltaMassFail_(0),
      nDzFail_(0),
      nPi0Fail_(0) {
  if (deltaMassRange_.size() < 2) {
    throw cms::Exception("Configuration")
        << "OniaPhotonProducer: deltaMassRange must contain at least 2 elements [min, max].";
  }
  produces<pat::CompositeCandidateCollection>();
  if (doKinematicRefit_) {
    produces<pat::CompositeCandidateCollection>("refit");
  }
}

void OniaPhotonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections
  desc.add<edm::InputTag>("oniaSrc", edm::InputTag("onia2MuMuPAT"));
  desc.add<edm::InputTag>("conversionSrc", edm::InputTag("oniaPhotonCandidates", "conversions"));
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("beamSpotSrc", edm::InputTag("offlineBeamSpot"));

  // Conversion cuts
  desc.add<double>("minConvPt", 0.0);
  desc.add<double>("minConvRho", 0.0);
  desc.add<bool>("rejectPi0", false);

  // Onia cuts
  desc.add<double>("minOniaPt", 0.0);
  desc.add<double>("maxOniaAbsY", 2.4);
  desc.add<bool>("requireTriggerMatch", false);

  // Matching cuts
  desc.add<double>("dzMax", 1.0);
  desc.add<std::vector<double>>("deltaMassRange", {0.0, 2.0});

  // Combined candidate cuts
  desc.add<double>("minCandPt", 0.0);
  desc.add<double>("maxCandAbsY", 2.4);

  // Kinematic refit
  desc.add<bool>("doKinematicRefit", false);
  desc.add<double>("constraintMass", 3.0969);

  // Parent type: "jpsi", "upsilon1S", "upsilon2S", "upsilon3S", "d0"
  desc.add<std::string>("parentType", "jpsi");

  descriptions.add("oniaPhotonProducer", desc);
}

void OniaPhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto chiCandColl = std::make_unique<pat::CompositeCandidateCollection>();
  auto refitCandColl = std::make_unique<pat::CompositeCandidateCollection>();

  if (doKinematicRefit_) {
    ttBuilder_ = &iSetup.getData(ttbToken_);
  }

  edm::Handle<pat::CompositeCandidateCollection> oniaHandle;
  iEvent.getByToken(oniaToken_, oniaHandle);

  edm::Handle<pat::CompositeCandidateCollection> convHandle;
  iEvent.getByToken(conversionToken_, convHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);

  if (!oniaHandle.isValid() || !convHandle.isValid()) {
    iEvent.put(std::move(chiCandColl));
    if (doKinematicRefit_) {
      iEvent.put(std::move(refitCandColl), "refit");
    }
    return;
  }

  reco::Vertex pv;
  if (pvHandle.isValid() && !pvHandle->empty() && !pvHandle->front().isFake()) {
    pv = pvHandle->front();
  } else if (bsHandle.isValid()) {
    pv = reco::Vertex(bsHandle->position(), bsHandle->covariance3D());
  }

  for (size_t oniaIdx = 0; oniaIdx < oniaHandle->size(); ++oniaIdx) {
    const auto& oniaCand = (*oniaHandle)[oniaIdx];

    if (oniaCand.pt() < minOniaPt_)
      continue;
    if (std::abs(oniaCand.rapidity()) > maxOniaAbsY_)
      continue;

    if (requireTriggerMatch_) {
      if (!oniaCand.hasUserInt("isTriggerMatched") || !oniaCand.userInt("isTriggerMatched")) {
        continue;
      }
    }

    reco::Candidate::Point oniaVtxPos = oniaCand.vertex();
    if (oniaCand.hasUserData("commonVertex")) {
      const reco::Vertex* oniaVtx = oniaCand.userData<reco::Vertex>("commonVertex");
      if (oniaVtx != nullptr) {
        oniaVtxPos = oniaVtx->position();
      }
    }

    const int resolvedOniaIdx =
        oniaCand.hasUserInt("oniaStateIndex") ? oniaCand.userInt("oniaStateIndex") : static_cast<int>(oniaIdx);
    const float nomMass =
        oniaCand.hasUserFloat("nominalMass") ? oniaCand.userFloat("nominalMass") : static_cast<float>(constraintMass_);

    for (size_t convIdx = 0; convIdx < convHandle->size(); ++convIdx) {
      const auto& convCand = (*convHandle)[convIdx];

      if (!passConversionCuts(convCand))
        continue;

      pat::CompositeCandidate chiCand =
          makeChiCandidate(oniaCand, convCand, resolvedOniaIdx, static_cast<int>(convIdx));
      if (oniaCand.hasUserInt("oniaPdgId")) {
        chiCand.addUserInt("oniaPdgId", oniaCand.userInt("oniaPdgId"));
      }
      chiCand.addUserFloat("nominalMass", nomMass);

      if (!passDeltaMassCut(chiCand, oniaCand)) {
        nDeltaMassFail_++;
        continue;
      }

      const float dz = computeDz(convCand, oniaVtxPos);
      chiCand.addUserFloat("dz", dz);

      if (!std::isfinite(dz) || std::abs(dz) > dzMax_) {
        nDzFail_++;
        continue;
      }

      if (rejectPi0_ && convCand.hasUserInt("flags")) {
        const int flags = convCand.userInt("flags") % 32;
        const bool pi0Fail = (flags & 8) != 0;
        if (pi0Fail) {
          nPi0Fail_++;
          continue;
        }
      }

      if (chiCand.pt() < minCandPt_)
        continue;
      if (std::abs(chiCand.rapidity()) > maxCandAbsY_)
        continue;

      chiCand.addUserFloat("convPt", convCand.pt());
      chiCand.addUserFloat("convEta", convCand.eta());
      chiCand.addUserFloat("convPhi", convCand.phi());
      chiCand.addUserFloat("convMass", convCand.mass());

      chiCand.addUserFloat("convVtxX", convCand.vertex().x());
      chiCand.addUserFloat("convVtxY", convCand.vertex().y());
      chiCand.addUserFloat("convVtxZ", convCand.vertex().z());
      chiCand.addUserFloat("convVtxRho", convCand.vertex().rho());

      if (convCand.hasUserInt("flags")) {
        chiCand.addUserInt("convFlags", convCand.userInt("flags"));
      }

      const reco::Track* convTk0 =
          convCand.hasUserData("track0") ? convCand.userData<reco::Track>("track0") : nullptr;
      const reco::Track* convTk1 =
          convCand.hasUserData("track1") ? convCand.userData<reco::Track>("track1") : nullptr;

      if (convTk0 != nullptr && convTk1 != nullptr) {
        const reco::Track* eleHigher = (convTk0->pt() > convTk1->pt()) ? convTk0 : convTk1;
        const reco::Track* eleLower = (convTk0->pt() > convTk1->pt()) ? convTk1 : convTk0;

        chiCand.addUserFloat("ele1_pt", eleHigher->pt());
        chiCand.addUserFloat("ele1_eta", eleHigher->eta());
        chiCand.addUserFloat("ele1_phi", eleHigher->phi());
        chiCand.addUserInt("ele1_charge", eleHigher->charge());
        chiCand.addUserFloat("ele1_chi2", eleHigher->normalizedChi2());
        chiCand.addUserInt("ele1_nHits", eleHigher->numberOfValidHits());
        chiCand.addUserFloat("ele1_dxy", eleHigher->dxy(oniaVtxPos));
        chiCand.addUserFloat("ele1_dz", eleHigher->dz(oniaVtxPos));

        chiCand.addUserFloat("ele2_pt", eleLower->pt());
        chiCand.addUserFloat("ele2_eta", eleLower->eta());
        chiCand.addUserFloat("ele2_phi", eleLower->phi());
        chiCand.addUserInt("ele2_charge", eleLower->charge());
        chiCand.addUserFloat("ele2_chi2", eleLower->normalizedChi2());
        chiCand.addUserInt("ele2_nHits", eleLower->numberOfValidHits());
        chiCand.addUserFloat("ele2_dxy", eleLower->dxy(oniaVtxPos));
        chiCand.addUserFloat("ele2_dz", eleLower->dz(oniaVtxPos));
      }

      chiCand.addUserFloat("oniaMass", oniaCand.mass());
      chiCand.addUserFloat("oniaPt", oniaCand.pt());
      chiCand.addUserFloat("oniaEta", oniaCand.eta());
      chiCand.addUserFloat("oniaPhi", oniaCand.phi());
      chiCand.addUserFloat("oniaY", oniaCand.rapidity());
      if (oniaCand.hasUserFloat("vProb")) {
        chiCand.addUserFloat("oniaVProb", oniaCand.userFloat("vProb"));
      }

      const double deltaM = chiCand.mass() - oniaCand.mass();
      chiCand.addUserFloat("deltaM", deltaM);

      const double correctedMass = deltaM + static_cast<double>(nomMass);
      chiCand.addUserFloat("correctedMass", correctedMass);

      const int rawIdx = static_cast<int>(chiCandColl->size());
      chiCandColl->push_back(chiCand);
      nCandidates_++;

      if (doKinematicRefit_) {
        pat::CompositeCandidate refitCand = doKinematicRefit(chiCand, oniaCand, convCand, pv, resolvedOniaIdx);
        if (refitCand.mass() > 0) {
          refitCand.addUserInt("rawIndex", rawIdx);
          refitCand.addUserInt("Index", rawIdx);
          refitCandColl->push_back(refitCand);
        }
      }
    }
  }

  iEvent.put(std::move(chiCandColl));
  if (doKinematicRefit_) {
    iEvent.put(std::move(refitCandColl), "refit");
  }
}

void OniaPhotonProducer::endStream() {
  edm::LogVerbatim("OniaPhotonProducer")
      << "OniaPhotonProducer summary (parentType=" << parentType_ << "): "
      << "candidates=" << nCandidates_.load() << ", deltaMassFail=" << nDeltaMassFail_.load()
      << ", dzFail=" << nDzFail_.load() << ", pi0Fail=" << nPi0Fail_.load();
}

pat::CompositeCandidate OniaPhotonProducer::makeChiCandidate(
    const pat::CompositeCandidate& onia,
    const pat::CompositeCandidate& photon,
    int resolvedOniaIdx,
    int convIdx) const {
  pat::CompositeCandidate chiCand;
  chiCand.addDaughter(onia, "dimuon");
  chiCand.addDaughter(photon, "photon");

  reco::Candidate::Point vtxPos = onia.vertex();
  if (onia.hasUserData("commonVertex")) {
    const reco::Vertex* oniaVtx = onia.userData<reco::Vertex>("commonVertex");
    if (oniaVtx != nullptr) {
      vtxPos = oniaVtx->position();
    }
  }
  chiCand.setVertex(vtxPos);

  reco::Candidate::LorentzVector chiP4 = onia.p4() + photon.p4();
  chiCand.setP4(chiP4);

  chiCand.addUserInt("oniaIdx", resolvedOniaIdx);
  chiCand.addUserInt("oniaIndex", resolvedOniaIdx);
  chiCand.addUserInt("convIdx", convIdx);
  chiCand.addUserInt("photonSourceIndex", convIdx);

  return chiCand;
}

float OniaPhotonProducer::computeDz(const pat::CompositeCandidate& conv,
                                    const reco::Candidate::Point& vtxPos) const {
  const reco::Candidate::LorentzVector& mom = conv.p4();
  const double rho = mom.Rho();
  if (rho <= 0.) {
    return std::numeric_limits<float>::quiet_NaN();
  }
  const reco::Candidate::Point& convVtx = conv.vertex();
  const double dz = (convVtx.Z() - vtxPos.Z()) -
                    ((convVtx.X() - vtxPos.X()) * mom.X() + (convVtx.Y() - vtxPos.Y()) * mom.Y()) /
                        rho * mom.Z() / rho;
  return static_cast<float>(dz);
}

bool OniaPhotonProducer::passDeltaMassCut(const pat::CompositeCandidate& chiCand,
                                          const pat::CompositeCandidate& oniaCand) const {
  const double deltaM = chiCand.mass() - oniaCand.mass();
  return (deltaM > deltaMassRange_[0] && deltaM < deltaMassRange_[1]);
}

bool OniaPhotonProducer::passConversionCuts(const pat::CompositeCandidate& conv) const {
  if (conv.pt() < minConvPt_)
    return false;
  if (conv.vertex().rho() < minConvRho_)
    return false;
  return true;
}

pat::CompositeCandidate OniaPhotonProducer::doKinematicRefit(
    const pat::CompositeCandidate& chiCand,
    const pat::CompositeCandidate& oniaCand,
    const pat::CompositeCandidate& photonCand,
    const reco::Vertex& pv,
    int oniaIdx) const {
  (void)chiCand;
  pat::CompositeCandidate invalidCand;
  invalidCand.setP4(reco::Candidate::LorentzVector(0, 0, 0, -1));

  if (ttBuilder_ == nullptr)
    return invalidCand;

  const reco::Candidate* dau1 = oniaCand.daughter("muon1");
  const reco::Candidate* dau2 = oniaCand.daughter("muon2");
  if (!dau1 || !dau2) {
    if (oniaCand.numberOfDaughters() >= 2) {
      dau1 = oniaCand.daughter(0);
      dau2 = oniaCand.daughter(1);
    }
  }
  if (!dau1 || !dau2)
    return invalidCand;

  const pat::Muon* mu1 = dynamic_cast<const pat::Muon*>(dau1);
  const pat::Muon* mu2 = dynamic_cast<const pat::Muon*>(dau2);
  if (!mu1 || !mu2)
    return invalidCand;

  const reco::TrackRef muTk1 = mu1->innerTrack();
  const reco::TrackRef muTk2 = mu2->innerTrack();
  if (muTk1.isNull() || muTk2.isNull())
    return invalidCand;

  const reco::Track* convTk0 =
      photonCand.hasUserData("track0") ? photonCand.userData<reco::Track>("track0") : nullptr;
  const reco::Track* convTk1 =
      photonCand.hasUserData("track1") ? photonCand.userData<reco::Track>("track1") : nullptr;
  if (!convTk0 || !convTk1)
    return invalidCand;

  std::vector<reco::TransientTrack> mumuTT;
  mumuTT.push_back(ttBuilder_->build(&muTk1));
  mumuTT.push_back(ttBuilder_->build(&muTk2));

  std::vector<reco::TransientTrack> eeTT;
  eeTT.push_back(ttBuilder_->build(*convTk0));
  eeTT.push_back(ttBuilder_->build(*convTk1));

  const ParticleMass zeroMass(0);
  const float zeroSigma = 1e-6f;
  const ParticleMass eleMass(OniaPhoton::elMass);
  const float eleSigma = 1e-6f;

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> photonParticles;
  photonParticles.push_back(pFactory.particle(eeTT[0], eleMass, 0.f, 0.f, eleSigma));
  photonParticles.push_back(pFactory.particle(eeTT[1], eleMass, 0.f, 0.f, eleSigma));

  KinematicParticleVertexFitter fitter;
  RefCountedKinematicTree photonTree = fitter.fit(photonParticles);

  if (photonTree->isEmpty() || !photonTree->isValid()) {
    edm::ParameterSet pSet;
    pSet.addParameter<double>("maxDistance", 3.0);
    pSet.addParameter<int>("maxNbrOfIterations", 10000);
    KinematicParticleVertexFitter fitter2(pSet);
    photonTree = fitter2.fit(photonParticles);
  }
  if (photonTree->isEmpty() || !photonTree->isValid())
    return invalidCand;

  KinematicParticleFitter csFitter;
  MassKinematicConstraint photonConstraint(zeroMass, zeroSigma);
  photonTree->movePointerToTheTop();
  photonTree = csFitter.fit(&photonConstraint, photonTree);
  if (photonTree->isEmpty() || !photonTree->isValid())
    return invalidCand;

  photonTree->movePointerToTheTop();
  RefCountedKinematicParticle fitPhoton = photonTree->currentParticle();
  if (!fitPhoton || !fitPhoton->currentState().isValid())
    return invalidCand;

  const ParticleMass muonMass(OniaPhoton::muMass);
  const float muonSigma = muonMass * 1e-6f;

  std::vector<RefCountedKinematicParticle> chiDaughters;
  chiDaughters.push_back(pFactory.particle(mumuTT[0], muonMass, 0.f, 0.f, muonSigma));
  chiDaughters.push_back(pFactory.particle(mumuTT[1], muonMass, 0.f, 0.f, muonSigma));
  chiDaughters.push_back(fitPhoton);

  const double fitConstraintMass =
      oniaCand.hasUserFloat("nominalMass") ? static_cast<double>(oniaCand.userFloat("nominalMass")) : constraintMass_;
  KinematicConstrainedVertexFitter constFitter;
  TwoTrackMassKinematicConstraint oniaMassConstraint(fitConstraintMass);
  RefCountedKinematicTree chiTree = constFitter.fit(chiDaughters, &oniaMassConstraint);
  if (chiTree->isEmpty() || !chiTree->isValid())
    return invalidCand;

  chiTree->movePointerToTheTop();
  RefCountedKinematicParticle fitChi = chiTree->currentParticle();
  RefCountedKinematicVertex chiVertex = chiTree->currentDecayVertex();
  if (!fitChi || !fitChi->currentState().isValid() || !chiVertex || !chiVertex->vertexIsValid())
    return invalidCand;

  const float chiM = fitChi->currentState().mass();
  const float chiPx = fitChi->currentState().kinematicParameters().momentum().x();
  const float chiPy = fitChi->currentState().kinematicParameters().momentum().y();
  const float chiPz = fitChi->currentState().kinematicParameters().momentum().z();
  const float chiVtxX = chiVertex->position().x();
  const float chiVtxY = chiVertex->position().y();
  const float chiVtxZ = chiVertex->position().z();
  const float chiVProb = ChiSquaredProbability(chiVertex->chiSquared(), chiVertex->degreesOfFreedom());

  bool child = chiTree->movePointerToTheFirstChild();
  if (!child)
    return invalidCand;
  RefCountedKinematicParticle fitMu1 = chiTree->currentParticle();
  if (!fitMu1 || !fitMu1->currentState().isValid())
    return invalidCand;

  const float mu1M = fitMu1->currentState().mass();
  const float mu1Q = fitMu1->currentState().particleCharge();
  const float mu1Px = fitMu1->currentState().kinematicParameters().momentum().x();
  const float mu1Py = fitMu1->currentState().kinematicParameters().momentum().y();
  const float mu1Pz = fitMu1->currentState().kinematicParameters().momentum().z();
  const float mu1E = std::sqrt(mu1M * mu1M + mu1Px * mu1Px + mu1Py * mu1Py + mu1Pz * mu1Pz);
  reco::CompositeCandidate recoMu1(
      mu1Q, math::XYZTLorentzVector(mu1Px, mu1Py, mu1Pz, mu1E), math::XYZPoint(chiVtxX, chiVtxY, chiVtxZ), 13);
  pat::CompositeCandidate patMu1(recoMu1);

  child = chiTree->movePointerToTheNextChild();
  if (!child)
    return invalidCand;
  RefCountedKinematicParticle fitMu2 = chiTree->currentParticle();
  if (!fitMu2 || !fitMu2->currentState().isValid())
    return invalidCand;

  const float mu2M = fitMu2->currentState().mass();
  const float mu2Q = fitMu2->currentState().particleCharge();
  const float mu2Px = fitMu2->currentState().kinematicParameters().momentum().x();
  const float mu2Py = fitMu2->currentState().kinematicParameters().momentum().y();
  const float mu2Pz = fitMu2->currentState().kinematicParameters().momentum().z();
  const float mu2E = std::sqrt(mu2M * mu2M + mu2Px * mu2Px + mu2Py * mu2Py + mu2Pz * mu2Pz);
  reco::CompositeCandidate recoMu2(
      mu2Q, math::XYZTLorentzVector(mu2Px, mu2Py, mu2Pz, mu2E), math::XYZPoint(chiVtxX, chiVtxY, chiVtxZ), 13);
  pat::CompositeCandidate patMu2(recoMu2);

  child = chiTree->movePointerToTheNextChild();
  if (!child)
    return invalidCand;
  RefCountedKinematicParticle fitGamma = chiTree->currentParticle();
  if (!fitGamma || !fitGamma->currentState().isValid())
    return invalidCand;

  const float gammaM = fitGamma->currentState().mass();
  const float gammaPx = fitGamma->currentState().kinematicParameters().momentum().x();
  const float gammaPy = fitGamma->currentState().kinematicParameters().momentum().y();
  const float gammaPz = fitGamma->currentState().kinematicParameters().momentum().z();
  const float gammaE = std::sqrt(gammaM * gammaM + gammaPx * gammaPx + gammaPy * gammaPy + gammaPz * gammaPz);
  reco::CompositeCandidate recoGamma(
      0, math::XYZTLorentzVector(gammaPx, gammaPy, gammaPz, gammaE), math::XYZPoint(chiVtxX, chiVtxY, chiVtxZ), 22);
  pat::CompositeCandidate patGamma(recoGamma);

  const float chiE = std::sqrt(chiM * chiM + chiPx * chiPx + chiPy * chiPy + chiPz * chiPz);
  reco::CompositeCandidate recoRefitCand(
      0, math::XYZTLorentzVector(chiPx, chiPy, chiPz, chiE), math::XYZPoint(chiVtxX, chiVtxY, chiVtxZ), 50551);

  pat::CompositeCandidate refitCand(recoRefitCand);
  refitCand.addUserFloat("vProb", chiVProb);
  refitCand.addUserInt("oniaIdx", oniaIdx);

  TVector3 vtx(chiVtxX, chiVtxY, 0);
  TVector3 pvtx(pv.x(), pv.y(), 0);
  TVector3 pperp(chiPx, chiPy, 0);
  TVector3 vdiff = vtx - pvtx;

  const double vdiffPerp = vdiff.Perp();
  const double pPerp = pperp.Perp();
  const double cosAlpha = (vdiffPerp > 0. && pPerp > 0.) ? vdiff.Dot(pperp) / (vdiffPerp * pPerp) : -999.;

  VertexDistanceXY vdistXY;
  reco::Vertex myVertex(*chiVertex);
  Measurement1D distXY = vdistXY.distance(myVertex, pv);

  const double ctauPV = (pPerp > 0. && cosAlpha > -900.) ? distXY.value() * cosAlpha * chiM / pPerp : -999.;

  AlgebraicVector3 vpperp(pperp.x(), pperp.y(), 0);
  GlobalError v1e = myVertex.error();
  GlobalError v2e = pv.error();
  AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
  const double sim = ROOT::Math::Similarity(vpperp, vXYe);
  const double ctauErrPV = (pperp.Perp2() > 0. && sim >= 0.) ? std::sqrt(sim) * chiM / pperp.Perp2() : -999.;

  refitCand.addUserFloat("ctauPV", ctauPV);
  refitCand.addUserFloat("ctauErrPV", ctauErrPV);
  refitCand.addUserFloat("cosAlpha", cosAlpha);

  pat::CompositeCandidate dimuon;
  dimuon.addDaughter(patMu1, "muon1");
  dimuon.addDaughter(patMu2, "muon2");
  dimuon.setP4(patMu1.p4() + patMu2.p4());

  refitCand.addDaughter(dimuon, "dimuon");
  refitCand.addDaughter(patGamma, "photon");

  return refitCand;
}

DEFINE_FWK_MODULE(OniaPhotonProducer);
