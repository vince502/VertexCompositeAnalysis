// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCTrackPairProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCTrackPairProducer.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

using reco::TransientTrack;

#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>

namespace {
reco::RecoChargedCandidate makeRecoCandidate(const reco::TrackRef& trackRef, double mass) {
  math::XYZTLorentzVector p4;
  const auto& track = *trackRef;
  const double momentum = track.p();
  const double energy = std::sqrt(momentum * momentum + mass * mass);
  p4.SetPxPyPzE(track.px(), track.py(), track.pz(), energy);
  reco::RecoChargedCandidate cand(trackRef->charge(), p4, track.vertex());
  cand.setTrack(trackRef);
  return cand;
}
}

ChiCTrackPairProducer::ChiCTrackPairProducer(const edm::ParameterSet& cfg)
  : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackCollection"))),
    vertexToken_(cfg.exists("vertexRecoAlgorithm") ? consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertexRecoAlgorithm")) : edm::EDGetTokenT<reco::VertexCollection>()),
    beamSpotToken_(cfg.exists("beamSpot") ? consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot")) : edm::EDGetTokenT<reco::BeamSpot>()),
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    daughterMass_(cfg.getParameter<double>("daughterMass")),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    maxTrackChi2_(cfg.getParameter<double>("maxTrackNormalizedChi2")),
    minTrackNHits_(cfg.getParameter<int>("minTrackNHits")),
    minTrackNPix_(cfg.exists("minTrackNPix") ? cfg.getParameter<int>("minTrackNPix") : 0),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    minPairPt_(cfg.getParameter<double>("minPairPt")),
    requiredChargeProduct_(cfg.existsAs<int>("requiredChargeProduct") ? cfg.getParameter<int>("requiredChargeProduct") : -1),
    useVertexFitting_(cfg.getParameter<bool>("useVertexFitting"))
{
  // Set default mass sigma if not provided (1% uncertainty)
  if (cfg.exists("daughterMassSigma")) {
    daughterMassSigma_ = cfg.getParameter<double>("daughterMassSigma");
  } else {
    daughterMassSigma_ = daughterMass_ * 0.01;
  }

  const auto& statePsets = cfg.getParameter<std::vector<edm::ParameterSet> >("states");
  states_.reserve(statePsets.size());
  for (const auto& ps : statePsets) {
    ChiStateConfig state;
    state.name = ps.getParameter<std::string>("name");
    state.pdgId = ps.getParameter<int>("pdgId");
    state.mass = ps.getParameter<double>("mass");
    state.massWindow = ps.getParameter<double>("massWindow");
    states_.push_back(state);
    produces<ChiCollection>(state.name);
  }
}

ChiCTrackPairProducer::~ChiCTrackPairProducer() = default;

void ChiCTrackPairProducer::beginJob() {}

void ChiCTrackPairProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracks;
  event.getByToken(trackToken_, tracks);
  if (!tracks.isValid())
    return;

  // Get primary vertices and beamspot if vertex fitting is enabled
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::BeamSpot> beamSpot;
  const MagneticField* magField = nullptr;

  if (useVertexFitting_) {
    if (!vertexToken_.isUninitialized()) {
      event.getByToken(vertexToken_, vertices);
    }
    if (!beamSpotToken_.isUninitialized()) {
      event.getByToken(beamSpotToken_, beamSpot);
    }
    magField = &iSetup.getData(bFieldToken_);

    if (!vertices.isValid() || !beamSpot.isValid()) {
      throw cms::Exception("InvalidInput") << "Vertex or BeamSpot not available but useVertexFitting is enabled";
    }
  }

  // Get best primary vertex for decay length calculations
  math::XYZPoint primaryVtx(0, 0, 0);
  double xVtxError = 0.0, yVtxError = 0.0, zVtxError = 0.0;
  const reco::Vertex* vtxPrimary = nullptr;
  bool isVtxPV = false;
  unsigned int bestVtxIdx = 0;

  if (useVertexFitting_ && vertices.isValid() && beamSpot.isValid()) {
    using namespace VertexCompositeProducerCommonTools;
    const reco::VertexCollection vtxCollection = *(vertices.product());
    auto [bestvtx, vtxIdx] = getBestVertex(vtxCollection, *beamSpot, 2);
    primaryVtx = bestvtx;
    bestVtxIdx = vtxIdx;
    isVtxPV = (vtxCollection.size() > 0 && vtxIdx < vtxCollection.size());

    if (isVtxPV) {
      vtxPrimary = &(vtxCollection[vtxIdx]);
      xVtxError = vtxPrimary->xError();
      yVtxError = vtxPrimary->yError();
      zVtxError = vtxPrimary->zError();
    } else {
      xVtxError = beamSpot->BeamWidthX();
      yVtxError = beamSpot->BeamWidthY();
      zVtxError = 0.0;
    }
  }

  std::vector<reco::TrackRef> selectedTracks;
  selectedTracks.reserve(tracks->size());

  for (std::size_t idx = 0; idx < tracks->size(); ++idx) {
    reco::TrackRef trackRef(tracks, idx);
    const auto& track = *trackRef;
    if (track.pt() < minTrackPt_)
      continue;
    if (std::abs(track.eta()) > maxTrackEta_)
      continue;
    if (track.normalizedChi2() > maxTrackChi2_)
      continue;
    if (track.numberOfValidHits() < minTrackNHits_)
      continue;
    if (minTrackNPix_ > 0 && track.hitPattern().numberOfValidPixelHits() < minTrackNPix_)
      continue;
    selectedTracks.push_back(trackRef);
  }

  if (selectedTracks.size() < 2)
    return;

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  AddFourMomenta addP4;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  for (std::size_t i = 0; i < selectedTracks.size(); ++i) {
    const auto& trackRef1 = selectedTracks[i];
    const auto& track1 = *trackRef1;

    for (std::size_t j = i + 1; j < selectedTracks.size(); ++j) {
      const auto& trackRef2 = selectedTracks[j];
      const auto& track2 = *trackRef2;

      const int chargeProduct = track1.charge() * track2.charge();
      if (requiredChargeProduct_ == -1 && chargeProduct >= 0)
        continue;  // Require opposite charge
      if (requiredChargeProduct_ == +1 && chargeProduct <= 0)
        continue;  // Require same charge
      // If requiredChargeProduct_ == 0, accept any charge combination

      if ((track1.pt() + track2.pt()) < minPairPt_)
        continue;

      auto dau1 = makeRecoCandidate(trackRef1, daughterMass_);
      auto dau2 = makeRecoCandidate(trackRef2, daughterMass_);

      for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
        const auto& state = states_[stateIdx];
        auto chi = std::make_unique<pat::CompositeCandidate>();
        chi->setPdgId(state.pdgId);
        chi->setCharge(dau1.charge() + dau2.charge());

        // Vertex fitting or simple average
        reco::Candidate::Point chiVtx(0, 0, 0);
        double vtxChi2 = -1.0;
        double vtxNdof = -1.0;
        double vtxProb = -1.0;
        bool vertexFitValid = false;
        SMatrixSym3D vtxCovMatrix;

        if (useVertexFitting_ && magField) {
          // Create TransientTracks
          std::vector<TransientTrack> transTracks;
          transTracks.reserve(2);
          transTracks.emplace_back(*trackRef1, magField);
          transTracks.emplace_back(*trackRef2, magField);

          // Create kinematic particles
          KinematicParticleFactoryFromTransientTrack pFactory;
          std::vector<RefCountedKinematicParticle> chiParticles;

          float massSigmaFloat = static_cast<float>(daughterMassSigma_);
          chiParticles.push_back(pFactory.particle(transTracks[0], daughterMass_, 0.0f, 0.0f, massSigmaFloat));
          chiParticles.push_back(pFactory.particle(transTracks[1], daughterMass_, 0.0f, 0.0f, massSigmaFloat));

          // Fit vertex
          KinematicParticleVertexFitter chiFitter;
          RefCountedKinematicTree chiVertex = chiFitter.fit(chiParticles);

          if (chiVertex->isValid()) {
            chiVertex->movePointerToTheTop();
            RefCountedKinematicParticle chiCand = chiVertex->currentParticle();
            RefCountedKinematicVertex chiDecayVertex = chiVertex->currentDecayVertex();

            if (chiCand->currentState().isValid() && chiDecayVertex->vertexIsValid()) {
              vertexFitValid = true;

              // Get vertex position
              GlobalPoint vtxPos = chiDecayVertex->position();
              chiVtx = reco::Candidate::Point(vtxPos.x(), vtxPos.y(), vtxPos.z());

              // Get vertex errors
              vtxChi2 = chiDecayVertex->chiSquared();
              vtxNdof = chiDecayVertex->degreesOfFreedom();
              vtxProb = TMath::Prob(vtxChi2, vtxNdof);

              // Get covariance matrix
              std::vector<double> vtxEVec;
              vtxEVec.push_back(chiDecayVertex->error().cxx());
              vtxEVec.push_back(chiDecayVertex->error().cyx());
              vtxEVec.push_back(chiDecayVertex->error().cyy());
              vtxEVec.push_back(chiDecayVertex->error().czx());
              vtxEVec.push_back(chiDecayVertex->error().czy());
              vtxEVec.push_back(chiDecayVertex->error().czz());
              vtxCovMatrix = SMatrixSym3D(vtxEVec.begin(), vtxEVec.end());

              // Update 4-momentum from fitted state
              GlobalVector chiMom = chiCand->currentState().globalMomentum();
              double chiE = chiCand->currentState().kinematicParameters().energy();
              chi->setP4(reco::Particle::LorentzVector(chiMom.x(), chiMom.y(), chiMom.z(), chiE));

              // Store vertex object
              reco::Vertex chiVtxObj(chiVtx, vtxCovMatrix, vtxChi2, vtxNdof, 2);
              chi->addUserData("Vtx", chiVtxObj);
            }
          }
        }

        // Fallback to simple average if vertex fitting failed or disabled
        if (!vertexFitValid) {
          chiVtx = reco::Candidate::Point((track1.vx() + track2.vx()) * 0.5,
                                          (track1.vy() + track2.vy()) * 0.5,
                                          (track1.vz() + track2.vz()) * 0.5);
        }

        chi->setVertex(chiVtx);
        chi->addDaughter(dau1, "Track1");
        chi->addDaughter(dau2, "Track2");

        // Set 4-momentum if not already set by vertex fit
        if (!vertexFitValid) {
          addP4.set(*chi);
        }

        const double mass = chi->mass();
        if (applyMassWindow_) {
          if (mass < state.mass - state.massWindow || mass > state.mass + state.massWindow)
            continue;
        }

        // Store vertex quality metrics
        if (vertexFitValid) {
          chi->addUserFloat("VtxChi2", vtxChi2);
          chi->addUserFloat("VtxNdof", vtxNdof);
          chi->addUserFloat("VtxProb", vtxProb);
          chi->addUserFloat("VtxNormalizedChi2", vtxNdof > 0 ? vtxChi2 / vtxNdof : -1.0);
        }

        // Calculate and store decay length information
        if (useVertexFitting_ && vertexFitValid && isVtxPV) {
          GlobalVector chiLineOfFlight(
            chiVtx.x() - primaryVtx.x(),
            chiVtx.y() - primaryVtx.y(),
            chiVtx.z() - primaryVtx.z()
          );

          double decayLength3D = chiLineOfFlight.mag();
          double decayLength2D = chiLineOfFlight.perp();

          // Calculate decay length significance
          SMatrixSym3D totalCov = vtxCovMatrix;
          if (vtxPrimary) {
            totalCov += vtxPrimary->covariance();
          } else {
            totalCov += beamSpot->rotatedCovariance3D();
          }

          SVector3 distanceVector3D(chiLineOfFlight.x(), chiLineOfFlight.y(), chiLineOfFlight.z());
          SVector3 distanceVector2D(chiLineOfFlight.x(), chiLineOfFlight.y(), 0.0);

          double sigmaDecayLength3D = (decayLength3D > 0) ?
            sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / decayLength3D : 999.0;
          double sigmaDecayLength2D = (decayLength2D > 0) ?
            sqrt(ROOT::Math::Similarity(totalCov, distanceVector2D)) / decayLength2D : 999.0;

          double decayLengthSig3D = (sigmaDecayLength3D > 0) ? decayLength3D / sigmaDecayLength3D : 0.0;
          double decayLengthSig2D = (sigmaDecayLength2D > 0) ? decayLength2D / sigmaDecayLength2D : 0.0;

          // Calculate pointing angles
          GlobalVector chiMom(chi->px(), chi->py(), chi->pz());
          double pointingAngle3D = angle(
            chiLineOfFlight.x(), chiLineOfFlight.y(), chiLineOfFlight.z(),
            chiMom.x(), chiMom.y(), chiMom.z()
          );
          double pointingAngle2D = angle(
            static_cast<double>(chiLineOfFlight.x()), static_cast<double>(chiLineOfFlight.y()), 0.0,
            static_cast<double>(chiMom.x()), static_cast<double>(chiMom.y()), 0.0
          );

          chi->addUserFloat("decaylength2D", decayLength2D);
          chi->addUserFloat("decaylength3D", decayLength3D);
          chi->addUserFloat("decaylengthsignif2D", decayLengthSig2D);
          chi->addUserFloat("decaylengthsignif3D", decayLengthSig3D);
          chi->addUserFloat("alpha2D", pointingAngle2D);
          chi->addUserFloat("alpha3D", pointingAngle3D);
          chi->addUserFloat("cosAlpha2D", std::cos(pointingAngle2D));
          chi->addUserFloat("cosAlpha3D", std::cos(pointingAngle3D));
        }

        // Calculate DCA between tracks
        if (useVertexFitting_ && magField) {
          std::vector<TransientTrack> transTracks;
          transTracks.emplace_back(*trackRef1, magField);
          transTracks.emplace_back(*trackRef2, magField);

          if (transTracks[0].impactPointTSCP().isValid() &&
              transTracks[1].impactPointTSCP().isValid()) {
            FreeTrajectoryState state1 = transTracks[0].impactPointTSCP().theState();
            FreeTrajectoryState state2 = transTracks[1].impactPointTSCP().theState();

            TwoTrackMinimumDistance minDist;
            minDist.calculate(state1, state2);
            double dca = std::abs(minDist.distance());

            chi->addUserFloat("trackDCA", dca);
          }
        }

        // Calculate impact parameter significances for each track
        if (useVertexFitting_ && isVtxPV) {
          std::array<reco::TrackRef, 2> trackRefs{{trackRef1, trackRef2}};
          for (size_t k = 0; k < 2; ++k) {
            const auto& trackRef = trackRefs[k];
            double dzvtx = trackRef->dz(primaryVtx);
            double dxyvtx = trackRef->dxy(primaryVtx);
            double dzerror = sqrt(trackRef->dzError() * trackRef->dzError() + zVtxError * zVtxError);
            double dxyerror = sqrt(trackRef->d0Error() * trackRef->d0Error() + xVtxError * yVtxError);

            double dauLongImpactSig = (dzerror > 0) ? dzvtx / dzerror : 0.0;
            double dauTransImpactSig = (dxyerror > 0) ? dxyvtx / dxyerror : 0.0;

            std::string suffix = std::to_string(k + 1);
            chi->addUserFloat("dauLongImpactSig_" + suffix, dauLongImpactSig);
            chi->addUserFloat("dauTransImpactSig_" + suffix, dauTransImpactSig);
          }
        }

        // Store primary vertex information
        if (useVertexFitting_ && isVtxPV && vtxPrimary) {
          chi->addUserInt("assocVtxIndex", static_cast<int>(bestVtxIdx));
          chi->addUserFloat("pvX", primaryVtx.x());
          chi->addUserFloat("pvY", primaryVtx.y());
          chi->addUserFloat("pvZ", primaryVtx.z());
          chi->addUserFloat("pvXError", xVtxError);
          chi->addUserFloat("pvYError", yVtxError);
          chi->addUserFloat("pvZError", zVtxError);
        }

        outputs[stateIdx]->push_back(*chi);
      }
    }
  }

  for (std::size_t idx = 0; idx < outputs.size(); ++idx) {
    event.put(std::move(outputs[idx]), states_[idx].name);
  }
}

void ChiCTrackPairProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCTrackPairProducer);
