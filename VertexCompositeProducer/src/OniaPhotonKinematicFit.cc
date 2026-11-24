// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniaPhotonKinematicFit
//
/**
 Description: Kinematic Fit for Onia + Photon

 Implementation:
     Original work from Stefano Argiro, and the Torino group
     Adapted for MINIAOD by Alberto Sanchez-Hernandez
     Unified for VertexCompositeAnalysis
*/

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

class OniaPhotonKinematicFit : public edm::stream::EDProducer<> {
public:
  explicit OniaPhotonKinematicFit(const edm::ParameterSet&);
  ~OniaPhotonKinematicFit() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> trackBuilderToken_;
  double upsilon_mass_;
  std::string product_name_;
  std::vector<int> oniaPDGs_;
  std::vector<double> nominalMasses_;
};

OniaPhotonKinematicFit::OniaPhotonKinematicFit(const edm::ParameterSet& iConfig)
    : chi_Label_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("chi_cand"))),
      trackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
      product_name_(iConfig.getParameter<std::string>("product_name")),
      oniaPDGs_(iConfig.getParameter<std::vector<int>>("oniaPDGs")),
      nominalMasses_(iConfig.getParameter<std::vector<double>>("nominalMasses")) {
  if (!oniaPDGs_.empty() && oniaPDGs_.size() != nominalMasses_.size()) {
    throw cms::Exception("Configuration") << "Each Onia hypothesis requires a nominal mass";
  }
  produces<pat::CompositeCandidateCollection>(product_name_);
}

void OniaPhotonKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto chicCompCandRefitColl = std::make_unique<pat::CompositeCandidateCollection>();

  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByToken(chi_Label_, chiCandHandle);
  if (!chiCandHandle.isValid()) {
    iEvent.put(std::move(chicCompCandRefitColl), product_name_);
    return;
  }

  const auto& theB = iSetup.getData(trackBuilderToken_);

  int indexConversion = -1;
  for (auto chiCand = chiCandHandle->begin(); chiCand != chiCandHandle->end(); ++chiCand) {
    ++indexConversion;

    double constraintMass = upsilon_mass_;
    if (!oniaPDGs_.empty()) {
      if (!chiCand->hasUserInt("oniaPdgId")) {
        throw cms::Exception("Configuration") << "Missing oniaPdgId on chi candidate for multi-hypothesis fit";
      }
      const auto hypothesis = std::find(oniaPDGs_.begin(), oniaPDGs_.end(), chiCand->userInt("oniaPdgId"));
      if (hypothesis == oniaPDGs_.end()) {
        throw cms::Exception("Configuration") << "Unknown Onia fit hypothesis";
      }
      constraintMass = nominalMasses_[std::distance(oniaPDGs_.begin(), hypothesis)];
    }

    const auto* dimuon = dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("dimuon"));
    const auto* photon = dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("photon"));
    if (dimuon == nullptr || photon == nullptr) {
      continue;
    }

    const auto* muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
    const auto* muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));
    if (muon1 == nullptr || muon2 == nullptr) {
      continue;
    }

    const reco::TrackRef jpsiTk0 = muon1->innerTrack();
    const reco::TrackRef jpsiTk1 = muon2->innerTrack();
    if (jpsiTk0.isNull() || jpsiTk1.isNull()) {
      continue;
    }

    const reco::Track* track0 = photon->hasUserData("track0") ? photon->userData<reco::Track>("track0") : nullptr;
    const reco::Track* track1 = photon->hasUserData("track1") ? photon->userData<reco::Track>("track1") : nullptr;
    if (track0 == nullptr || track1 == nullptr) {
      continue;
    }

    const reco::Vertex* pvPtr =
        dimuon->hasUserData("PVwithmuons") ? dimuon->userData<reco::Vertex>("PVwithmuons") : nullptr;
    const reco::Vertex thePrimaryV =
        (pvPtr != nullptr) ? *pvPtr : reco::Vertex(dimuon->vertex(), reco::Vertex::Error());

    std::vector<reco::TransientTrack> MuMuTT;
    MuMuTT.push_back(theB.build(&jpsiTk0));
    MuMuTT.push_back(theB.build(&jpsiTk1));

    std::vector<reco::TransientTrack> EETT;
    EETT.push_back(theB.build(*track0));
    EETT.push_back(theB.build(*track1));

    const ParticleMass zero_mass(0);
    const float zero_sigma = 1E-6f;
    const ParticleMass eleMass(0.000511);
    const float eleSigma = 1E-6f;

    KinematicParticleFactoryFromTransientTrack pFactory;
    std::vector<RefCountedKinematicParticle> PhotonParticles;
    PhotonParticles.push_back(pFactory.particle(EETT[0], eleMass, 0.f, 0.f, eleSigma));
    PhotonParticles.push_back(pFactory.particle(EETT[1], eleMass, 0.f, 0.f, eleSigma));

    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree photonVertexFitTree = fitter.fit(PhotonParticles);

    if (photonVertexFitTree->isEmpty() || !photonVertexFitTree->isValid()) {
      edm::ParameterSet pSet;
      pSet.addParameter<double>("maxDistance", 3.0);
      pSet.addParameter<int>("maxNbrOfIterations", 10000);
      KinematicParticleVertexFitter fitter2(pSet);
      photonVertexFitTree = fitter2.fit(PhotonParticles);
    }

    if (photonVertexFitTree->isEmpty() || !photonVertexFitTree->isValid()) {
      continue;
    }

    KinematicParticleFitter csFitterPhoton;
    MassKinematicConstraint pho_c(zero_mass, zero_sigma);
    photonVertexFitTree->movePointerToTheTop();
    photonVertexFitTree = csFitterPhoton.fit(&pho_c, photonVertexFitTree);

    if (photonVertexFitTree->isEmpty() || !photonVertexFitTree->isValid()) {
      continue;
    }

    const ParticleMass muonMass(0.1056584);
    const float muonSigma = muonMass * 1E-6f;

    photonVertexFitTree->movePointerToTheTop();
    RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
    if (!fitPhoton || !fitPhoton->currentState().isValid()) {
      continue;
    }

    std::vector<RefCountedKinematicParticle> allChiBDaughters;
    allChiBDaughters.push_back(pFactory.particle(MuMuTT[0], muonMass, 0.f, 0.f, muonSigma));
    allChiBDaughters.push_back(pFactory.particle(MuMuTT[1], muonMass, 0.f, 0.f, muonSigma));
    allChiBDaughters.push_back(fitPhoton);

    KinematicConstrainedVertexFitter constVertexFitter;
    TwoTrackMassKinematicConstraint upsilon_mtc(constraintMass);
    RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters, &upsilon_mtc);

    if (ChiBTree->isEmpty() || !ChiBTree->isValid()) {
      continue;
    }

    ChiBTree->movePointerToTheTop();
    RefCountedKinematicParticle fitChiB = ChiBTree->currentParticle();
    RefCountedKinematicVertex ChiBDecayVertex = ChiBTree->currentDecayVertex();

    if (!fitChiB || !fitChiB->currentState().isValid() || !ChiBDecayVertex || !ChiBDecayVertex->vertexIsValid()) {
      continue;
    }

    const float ChiBM_fit = fitChiB->currentState().mass();
    const float ChiBPx_fit = fitChiB->currentState().kinematicParameters().momentum().x();
    const float ChiBPy_fit = fitChiB->currentState().kinematicParameters().momentum().y();
    const float ChiBPz_fit = fitChiB->currentState().kinematicParameters().momentum().z();
    const float ChiBVtxX_fit = ChiBDecayVertex->position().x();
    const float ChiBVtxY_fit = ChiBDecayVertex->position().y();
    const float ChiBVtxZ_fit = ChiBDecayVertex->position().z();
    const float ChiBVtxP_fit = ChiSquaredProbability(static_cast<double>(ChiBDecayVertex->chiSquared()),
                                                     static_cast<double>(ChiBDecayVertex->degreesOfFreedom()));

    reco::CompositeCandidate recoChib(
        0,
        math::XYZTLorentzVector(
            ChiBPx_fit,
            ChiBPy_fit,
            ChiBPz_fit,
            std::sqrt(ChiBM_fit * ChiBM_fit + ChiBPx_fit * ChiBPx_fit + ChiBPy_fit * ChiBPy_fit + ChiBPz_fit * ChiBPz_fit)),
        math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit),
        50551);

    pat::CompositeCandidate patChib(recoChib);
    patChib.addUserFloat("vProb", ChiBVtxP_fit);
    patChib.addUserInt("Index", indexConversion);
    patChib.addUserInt("rawIndex", indexConversion);

    VertexDistanceXY vdistXY;
    reco::Vertex myVertex = *ChiBDecayVertex;

    TVector3 vtx(ChiBVtxX_fit, ChiBVtxY_fit, 0);
    TVector3 pvtx(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
    TVector3 pperp(ChiBPx_fit, ChiBPy_fit, 0);
    AlgebraicVector3 vpperp(pperp.x(), pperp.y(), 0);

    TVector3 vdiff = vtx - pvtx;
    const double vdiffPerp = vdiff.Perp();
    const double pPerp = pperp.Perp();
    const double cosAlpha = (vdiffPerp > 0. && pPerp > 0.) ? vdiff.Dot(pperp) / (vdiffPerp * pPerp) : -999.;
    Measurement1D distXY = vdistXY.distance(myVertex, thePrimaryV);
    const double ctauPV = (pPerp > 0. && cosAlpha > -900.) ? distXY.value() * cosAlpha * ChiBM_fit / pPerp : -999.;
    GlobalError v1e = myVertex.error();
    GlobalError v2e = thePrimaryV.error();
    AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
    const double sim = ROOT::Math::Similarity(vpperp, vXYe);
    const double ctauErrPV = (pperp.Perp2() > 0. && sim >= 0.) ? std::sqrt(sim) * ChiBM_fit / pperp.Perp2() : -999.;

    patChib.addUserFloat("ctauPV", ctauPV);
    patChib.addUserFloat("ctauErrPV", ctauErrPV);
    patChib.addUserFloat("cosAlpha", cosAlpha);

    bool child = ChiBTree->movePointerToTheFirstChild();
    if (!child)
      continue;
    RefCountedKinematicParticle fitMu1 = ChiBTree->currentParticle();
    if (!fitMu1 || !fitMu1->currentState().isValid())
      continue;

    const float mu1M_fit = fitMu1->currentState().mass();
    const float mu1Q_fit = fitMu1->currentState().particleCharge();
    const float mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
    const float mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
    const float mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
    reco::CompositeCandidate recoMu1(
        mu1Q_fit,
        math::XYZTLorentzVector(
            mu1Px_fit,
            mu1Py_fit,
            mu1Pz_fit,
            std::sqrt(mu1M_fit * mu1M_fit + mu1Px_fit * mu1Px_fit + mu1Py_fit * mu1Py_fit + mu1Pz_fit * mu1Pz_fit)),
        math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit),
        13);
    pat::CompositeCandidate patMu1(recoMu1);

    child = ChiBTree->movePointerToTheNextChild();
    if (!child)
      continue;
    RefCountedKinematicParticle fitMu2 = ChiBTree->currentParticle();
    if (!fitMu2 || !fitMu2->currentState().isValid())
      continue;

    const float mu2M_fit = fitMu2->currentState().mass();
    const float mu2Q_fit = fitMu2->currentState().particleCharge();
    const float mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
    const float mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
    const float mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
    reco::CompositeCandidate recoMu2(
        mu2Q_fit,
        math::XYZTLorentzVector(
            mu2Px_fit,
            mu2Py_fit,
            mu2Pz_fit,
            std::sqrt(mu2M_fit * mu2M_fit + mu2Px_fit * mu2Px_fit + mu2Py_fit * mu2Py_fit + mu2Pz_fit * mu2Pz_fit)),
        math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit),
        13);
    pat::CompositeCandidate patMu2(recoMu2);

    pat::CompositeCandidate ups;
    ups.addDaughter(patMu1, "muon1");
    ups.addDaughter(patMu2, "muon2");
    ups.setP4(patMu1.p4() + patMu2.p4());

    child = ChiBTree->movePointerToTheNextChild();
    if (!child)
      continue;
    RefCountedKinematicParticle fitGamma = ChiBTree->currentParticle();
    if (!fitGamma || !fitGamma->currentState().isValid())
      continue;

    const float gammaM_fit = fitGamma->currentState().mass();
    const float gammaPx_fit = fitGamma->currentState().kinematicParameters().momentum().x();
    const float gammaPy_fit = fitGamma->currentState().kinematicParameters().momentum().y();
    const float gammaPz_fit = fitGamma->currentState().kinematicParameters().momentum().z();
    reco::CompositeCandidate recoGamma(
        0,
        math::XYZTLorentzVector(
            gammaPx_fit,
            gammaPy_fit,
            gammaPz_fit,
            std::sqrt(gammaM_fit * gammaM_fit + gammaPx_fit * gammaPx_fit + gammaPy_fit * gammaPy_fit +
                      gammaPz_fit * gammaPz_fit)),
        math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit),
        22);
    pat::CompositeCandidate patGamma(recoGamma);

    patChib.addDaughter(ups, "dimuon");
    patChib.addDaughter(patGamma, "photon");

    chicCompCandRefitColl->push_back(patChib);
  }

  iEvent.put(std::move(chicCompCandRefitColl), product_name_);
}

void OniaPhotonKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("chi_cand", edm::InputTag("chiCConventionalPhotonStates"));
  desc.add<double>("upsilon_mass", 3.0969);
  desc.add<std::string>("product_name", "y1S");
  desc.add<std::vector<int>>("oniaPDGs", {});
  desc.add<std::vector<double>>("nominalMasses", {});
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(OniaPhotonKinematicFit);
