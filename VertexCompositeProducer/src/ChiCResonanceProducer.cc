// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCResonanceProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCResonanceProducer.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
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

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <vector>

ChiCResonanceProducer::ChiCResonanceProducer(const edm::ParameterSet& cfg)
  : resonanceToken_(consumes<ResonanceCollection>(cfg.getParameter<edm::InputTag>("resonanceCollection"))),
    vertexToken_(cfg.exists("vertexRecoAlgorithm") ? consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertexRecoAlgorithm")) : edm::EDGetTokenT<reco::VertexCollection>()),
    beamSpotToken_(cfg.exists("beamSpot") ? consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot")) : edm::EDGetTokenT<reco::BeamSpot>()),
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    requireUniqueTracks_(cfg.getParameter<bool>("requireUniqueTracks")),
    useVertexFitting_(cfg.getParameter<bool>("useVertexFitting"))
{

  // Set default resonance mass sigmas if not provided (for Ks: ~0.497 GeV with ~1% uncertainty)
  if (cfg.exists("resonanceMassSigmas")) {
    resonanceMassSigmas_ = cfg.getParameter<std::vector<double> >("resonanceMassSigmas");
  } else {
    // Default: 1% uncertainty for Ks mass (0.497614 GeV)
    resonanceMassSigmas_ = {0.00497614, 0.00497614};
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

ChiCResonanceProducer::~ChiCResonanceProducer() = default;

void ChiCResonanceProducer::beginJob() {}

void ChiCResonanceProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {
  edm::Handle<ResonanceCollection> resonances;
  event.getByToken(resonanceToken_, resonances);

  if (!resonances.isValid() || resonances->size() < 2)
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

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  const auto& coll = *resonances;
  AddFourMomenta addP4;

  // Helper function to extract tracks from a resonance (Ks)
  auto extractTracks = [](const reco::VertexCompositeCandidate& resonance) -> std::vector<reco::TrackRef> {
    std::vector<reco::TrackRef> tracks;
    for (size_t idx = 0; idx < resonance.numberOfDaughters(); ++idx) {
      const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(resonance.daughter(idx));
      if (dau && dau->track().isNonnull()) {
        tracks.push_back(dau->track());
      }
    }
    return tracks;
  };

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  for (std::size_t i = 0; i < coll.size(); ++i) {
    const auto& first = coll[i];

    for (std::size_t j = i + 1; j < coll.size(); ++j) {
      const auto& second = coll[j];

      if (requireUniqueTracks_ && shareTracks(first, second))
        continue;

      // Extract tracks from both resonances
      std::vector<reco::TrackRef> tracks1 = extractTracks(first);
      std::vector<reco::TrackRef> tracks2 = extractTracks(second);

      if (tracks1.size() != 2 || tracks2.size() != 2) {
        continue; // Skip if we don't have exactly 2 tracks per resonance
      }

      // Combine all tracks (total: 4 tracks)
      std::vector<reco::TrackRef> allTracks;
      allTracks.reserve(4);
      allTracks.insert(allTracks.end(), tracks1.begin(), tracks1.end());
      allTracks.insert(allTracks.end(), tracks2.begin(), tracks2.end());

      for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
        const auto& state = states_[stateIdx];
        auto chi = std::make_unique<pat::CompositeCandidate>();
        chi->setPdgId(state.pdgId);
        chi->setCharge(first.charge() + second.charge());

        // Vertex fitting or simple average
        reco::Candidate::Point chiVtx(0, 0, 0);
        double vtxChi2 = -1.0;
        double vtxNdof = -1.0;
        double vtxProb = -1.0;
        bool vertexFitValid = false;
        SMatrixSym3D vtxCovMatrix;

        if (useVertexFitting_ && magField) {
          // Create TransientTracks from all 4 tracks
          std::vector<TransientTrack> transTracks;
          transTracks.reserve(4);
          for (const auto& trackRef : allTracks) {
            TransientTrack transTrack(*trackRef, magField);
            transTracks.push_back(transTrack);
          }

          // Create kinematic particles
          // For ChiC → Ks Ks, we treat each Ks as a single particle with its mass
          KinematicParticleFactoryFromTransientTrack pFactory;
          std::vector<RefCountedKinematicParticle> chiParticles;

          // First, fit the two Ks separately to get their 4-momenta, then combine
          // Actually, we need to fit all 4 tracks together, but treat them as 2 Ks
          // For simplicity, we'll fit all 4 tracks as pions, then constrain to Ks masses
          // But kinematic fitter expects individual tracks, so we fit all 4 tracks together
          // and use pion masses for the tracks (since Ks decays to π+π-)
          const double pionMass = 0.13957018;

          for (size_t k = 0; k < 4; ++k) {
            double massSigma = (k < resonanceMassSigmas_.size()) ? resonanceMassSigmas_[k] : resonanceMassSigmas_.back();
            float massSigmaFloat = static_cast<float>(massSigma);
            chiParticles.push_back(pFactory.particle(transTracks[k], pionMass, 0.0f, 0.0f, massSigmaFloat));
          }

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
              reco::Vertex chiVtxObj(chiVtx, vtxCovMatrix, vtxChi2, vtxNdof, 4);
              chi->addUserData("Vtx", chiVtxObj);
            }
          }
        }

        // Fallback to simple average if vertex fitting failed or disabled
        if (!vertexFitValid) {
          const double vx = 0.5 * (first.vx() + second.vx());
          const double vy = 0.5 * (first.vy() + second.vy());
          const double vz = 0.5 * (first.vz() + second.vz());
          chiVtx = reco::Candidate::Point(vx, vy, vz);
        }

        chi->setVertex(chiVtx);
        chi->addDaughter(first, "Resonance1");
        chi->addDaughter(second, "Resonance2");

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
          for (const auto& trackRef : allTracks) {
            transTracks.emplace_back(*trackRef, magField);
          }

          // Calculate minimum DCA among all track pairs
          double minDCA = std::numeric_limits<double>::max();
          double maxDCA = -1.0;
          double sumDCA = 0.0;
          int nPairs = 0;

          for (size_t k = 0; k < 4; ++k) {
            for (size_t l = k + 1; l < 4; ++l) {
              if (!transTracks[k].impactPointTSCP().isValid() ||
                  !transTracks[l].impactPointTSCP().isValid()) continue;

              FreeTrajectoryState state1 = transTracks[k].impactPointTSCP().theState();
              FreeTrajectoryState state2 = transTracks[l].impactPointTSCP().theState();

              TwoTrackMinimumDistance minDist;
              minDist.calculate(state1, state2);
              double dca = std::abs(minDist.distance());

              minDCA = std::min(minDCA, dca);
              maxDCA = std::max(maxDCA, dca);
              sumDCA += dca;
              nPairs++;
            }
          }

          if (nPairs > 0) {
            chi->addUserFloat("trackDCA_min", minDCA);
            chi->addUserFloat("trackDCA_max", maxDCA);
            chi->addUserFloat("trackDCA_avg", sumDCA / nPairs);
          }
        }

        // Calculate impact parameter significances for each track
        if (useVertexFitting_ && isVtxPV) {
          for (size_t k = 0; k < 4; ++k) {
            const auto& trackRef = allTracks[k];
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

void ChiCResonanceProducer::endJob() {}

bool ChiCResonanceProducer::shareTracks(const reco::VertexCompositeCandidate& first,
                                        const reco::VertexCompositeCandidate& second) const {
  std::set<reco::TrackRef> tracks;
  for (size_t idx = 0; idx < first.numberOfDaughters(); ++idx) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(first.daughter(idx));
    if (dau && dau->track().isNonnull()) {
      tracks.insert(dau->track());
    }
  }

  for (size_t idx = 0; idx < second.numberOfDaughters(); ++idx) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(second.daughter(idx));
    if (dau && dau->track().isNonnull() && tracks.count(dau->track())) {
      return true;
    }
  }

  return false;
}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCResonanceProducer);
