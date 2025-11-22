// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCFourTrackProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCFourTrackProducer.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

using reco::TransientTrack;

#include <Eigen/Dense>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>

namespace {
reco::RecoChargedCandidate makeRecoCandidate(const reco::TrackRef& trackRef, double mass) {
  const auto& track = *trackRef;
  const double momentum2 = track.momentum().mag2();
  const double energy = std::sqrt(momentum2 + mass * mass);
  math::XYZTLorentzVector p4(track.px(), track.py(), track.pz(), energy);
  reco::RecoChargedCandidate cand(trackRef->charge(), p4, track.vertex());
  cand.setTrack(trackRef);
  return cand;
}

bool hasTwoPositiveTwoNegative(const std::array<int, 4>& charges) {
  int nPos = 0;
  int nNeg = 0;
  for (auto q : charges) {
    if (q > 0)
      ++nPos;
    else if (q < 0)
      ++nNeg;
  }
  return (nPos == 2 && nNeg == 2);
}

struct EventShapeResult {
  double sphericity{0.0};
  std::array<double, 3> eigenvalues{{0.0, 0.0, 0.0}};
};

EventShapeResult computeEventShape(const std::array<reco::RecoChargedCandidate, 4>& daughters) {
  EventShapeResult result;
  Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
  double sumP2 = 0.0;
  for (const auto& dau : daughters) {
    Eigen::Vector3d p(dau.px(), dau.py(), dau.pz());
    tensor += p * p.transpose();
    sumP2 += p.squaredNorm();
  }

  if (sumP2 <= 0.0)
    return result;

  tensor /= sumP2;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(tensor);
  if (eigenSolver.info() != Eigen::Success)
    return result;

  result.eigenvalues = {{eigenSolver.eigenvalues()(0), eigenSolver.eigenvalues()(1), eigenSolver.eigenvalues()(2)}};
  std::sort(result.eigenvalues.begin(), result.eigenvalues.end(), std::greater<double>());
  result.sphericity = 1.5 * (result.eigenvalues[1] + result.eigenvalues[2]);
  return result;
}

double computeAcoplanarity(const std::array<reco::RecoChargedCandidate, 4>& daughters) {
  std::array<std::pair<double, std::size_t>, 4> ptIndex;
  for (std::size_t idx = 0; idx < daughters.size(); ++idx) {
    ptIndex[idx] = std::make_pair(daughters[idx].pt(), idx);
  }
  std::sort(ptIndex.begin(), ptIndex.end(), [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

  const auto& lead = daughters[ptIndex[0].second];
  const auto& sublead = daughters[ptIndex[1].second];
  const double deltaPhi = reco::deltaPhi(lead.phi(), sublead.phi());
  return 1.0 - std::abs(deltaPhi) / M_PI;
}
}

ChiCFourTrackProducer::ChiCFourTrackProducer(const edm::ParameterSet& cfg)
  : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackCollection"))),
    vertexToken_(cfg.exists("vertexRecoAlgorithm") ? consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertexRecoAlgorithm")) : edm::EDGetTokenT<reco::VertexCollection>()),
    beamSpotToken_(cfg.exists("beamSpot") ? consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot")) : edm::EDGetTokenT<reco::BeamSpot>()),
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    daughterMasses_(cfg.getParameter<std::vector<double> >("daughterMasses")),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    maxTrackChi2_(cfg.getParameter<double>("maxTrackNormalizedChi2")),
    minTrackNHits_(cfg.getParameter<int>("minTrackNHits")),
    minTrackNPix_(cfg.exists("minTrackNPix") ? cfg.getParameter<int>("minTrackNPix") : 0),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    minCandidatePt_(cfg.getParameter<double>("minCandidatePt")),
    minAcoplanarity_(cfg.getParameter<double>("minAcoplanarity")),
    maxSphericity_(cfg.getParameter<double>("maxSphericity")),
    maxCandidateAbsEta_(cfg.getParameter<double>("maxCandidateAbsEta")),
    storeEventShape_(cfg.getParameter<bool>("storeEventShape")),
    useVertexFitting_(cfg.getParameter<bool>("useVertexFitting"))
{
  if (daughterMasses_.empty()) {
    throw cms::Exception("InvalidConfiguration") << "Parameter 'daughterMasses' must contain at least one value.";
  }

  // Set default mass sigmas if not provided
  if (cfg.exists("daughterMassSigmas")) {
    daughterMassSigmas_ = cfg.getParameter<std::vector<double> >("daughterMassSigmas");
  } else {
    // Default: 1% uncertainty for pion mass
    daughterMassSigmas_.resize(daughterMasses_.size());
    for (size_t i = 0; i < daughterMasses_.size(); ++i) {
      daughterMassSigmas_[i] = daughterMasses_[i] * 0.01;
    }
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

ChiCFourTrackProducer::~ChiCFourTrackProducer() = default;

void ChiCFourTrackProducer::beginJob() {}

void ChiCFourTrackProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracks;
  event.getByToken(trackToken_, tracks);
  if (!tracks.isValid())
    return;

  // Get primary vertices and beamspot
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::BeamSpot> beamSpot;
  const MagneticField* magField = nullptr;
  
  if (useVertexFitting_) {
    event.getByToken(vertexToken_, vertices);
    event.getByToken(beamSpotToken_, beamSpot);
    magField = &iSetup.getData(bFieldToken_);
    
    if (!vertices.isValid() || !beamSpot.isValid()) {
      throw cms::Exception("InvalidInput") << "Vertex or BeamSpot not available but useVertexFitting is enabled";
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

  if (selectedTracks.size() < 4)
    return;

  // Get best primary vertex for decay length calculations
  math::XYZPoint primaryVtx(0, 0, 0);
  double xVtxError = 0.0, yVtxError = 0.0, zVtxError = 0.0;
  const reco::Vertex* vtxPrimary = nullptr;
  bool isVtxPV = false;
  unsigned int bestVtxIdx = 0;
  
  if (useVertexFitting_) {
    if (!vertexToken_.isUninitialized()) {
      event.getByToken(vertexToken_, vertices);
    }
    if (!beamSpotToken_.isUninitialized()) {
      event.getByToken(beamSpotToken_, beamSpot);
    }
    if (vertices.isValid() && beamSpot.isValid()) {
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
  }

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  AddFourMomenta addP4;
  const bool applyAcoplanarityCut = (minAcoplanarity_ > 0.0);
  const bool applySphericityCut = (maxSphericity_ >= 0.0 && maxSphericity_ < std::numeric_limits<double>::infinity());

  const auto massForIndex = [this](unsigned int index) {
    if (daughterMasses_.size() == 1)
      return daughterMasses_.front();
    if (index < daughterMasses_.size())
      return daughterMasses_[index];
    return daughterMasses_.back();
  };

  const std::size_t nTracks = selectedTracks.size();
  for (std::size_t i = 0; i < nTracks - 3; ++i) {
    const auto& ref1 = selectedTracks[i];
    for (std::size_t j = i + 1; j < nTracks - 2; ++j) {
      const auto& ref2 = selectedTracks[j];
      for (std::size_t k = j + 1; k < nTracks - 1; ++k) {
        const auto& ref3 = selectedTracks[k];
        for (std::size_t l = k + 1; l < nTracks; ++l) {
          const auto& ref4 = selectedTracks[l];

          const std::array<int, 4> charges{{ref1->charge(), ref2->charge(), ref3->charge(), ref4->charge()}};
          if (!hasTwoPositiveTwoNegative(charges))
            continue;

          const std::array<reco::RecoChargedCandidate, 4> daughters{{
              makeRecoCandidate(ref1, massForIndex(0)),
              makeRecoCandidate(ref2, massForIndex(1)),
              makeRecoCandidate(ref3, massForIndex(2)),
              makeRecoCandidate(ref4, massForIndex(3))}};

          for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
            const auto& state = states_[stateIdx];
            auto chi = std::make_unique<pat::CompositeCandidate>();
            chi->setPdgId(state.pdgId);
            chi->setCharge(0);

            // Vertex fitting or simple average
            reco::Candidate::Point chiVtx(0, 0, 0);
            double vtxChi2 = -1.0;
            double vtxNdof = -1.0;
            double vtxProb = -1.0;
            bool vertexFitValid = false;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            SMatrixSym3D vtxCovMatrix;
            
            if (useVertexFitting_ && magField) {
              // Create TransientTracks
              std::array<reco::TrackRef, 4> trackRefs{{ref1, ref2, ref3, ref4}};
              std::vector<TransientTrack> transTracks;
              transTracks.reserve(4);
              
              for (const auto& trackRef : trackRefs) {
                TransientTrack transTrack(*trackRef, magField);
                transTracks.push_back(transTrack);
              }
              
              // Create kinematic particles
              KinematicParticleFactoryFromTransientTrack pFactory;
              std::vector<RefCountedKinematicParticle> chiParticles;
              
              for (size_t i = 0; i < 4; ++i) {
                double mass = massForIndex(i);
                double massSigma = (i < daughterMassSigmas_.size()) ? daughterMassSigmas_[i] : daughterMassSigmas_.back();
                float massSigmaFloat = static_cast<float>(massSigma);
                chiParticles.push_back(pFactory.particle(transTracks[i], mass, 0.0f, 0.0f, massSigmaFloat));
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
              const auto& v1 = daughters[0].vertex();
              const auto& v2 = daughters[1].vertex();
              const auto& v3 = daughters[2].vertex();
              const auto& v4 = daughters[3].vertex();
              chiVtx = reco::Candidate::Point(
                0.25 * (v1.x() + v2.x() + v3.x() + v4.x()),
                0.25 * (v1.y() + v2.y() + v3.y() + v4.y()),
                0.25 * (v1.z() + v2.z() + v3.z() + v4.z())
              );
            }
            
            chi->setVertex(chiVtx);

            chi->addDaughter(daughters[0], "Track1");
            chi->addDaughter(daughters[1], "Track2");
            chi->addDaughter(daughters[2], "Track3");
            chi->addDaughter(daughters[3], "Track4");

            // Set 4-momentum if not already set by vertex fit
            if (!vertexFitValid) {
              addP4.set(*chi);
            }

            if (chi->pt() < minCandidatePt_)
              continue;

            const double mass = chi->mass();
            if (applyMassWindow_) {
              if (mass < state.mass - state.massWindow || mass > state.mass + state.massWindow)
                continue;
            }

            const double absEta = std::abs(chi->eta());
            if (maxCandidateAbsEta_ >= 0.0 && absEta > maxCandidateAbsEta_)
              continue;

            double acoplanarity = 0.0;
            EventShapeResult eventShape;

            if (applyAcoplanarityCut || applySphericityCut || storeEventShape_) {
              acoplanarity = computeAcoplanarity(daughters);
              eventShape = computeEventShape(daughters);

              if (applyAcoplanarityCut && acoplanarity < minAcoplanarity_)
                continue;

              if (applySphericityCut && eventShape.sphericity > maxSphericity_)
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
              std::array<reco::TrackRef, 4> trackRefs{{ref1, ref2, ref3, ref4}};
              std::vector<TransientTrack> transTracks;
              for (const auto& trackRef : trackRefs) {
                transTracks.emplace_back(*trackRef, magField);
              }
              
              // Calculate minimum DCA among all track pairs
              double minDCA = std::numeric_limits<double>::max();
              double maxDCA = -1.0;
              double sumDCA = 0.0;
              int nPairs = 0;
              
              for (size_t i = 0; i < 4; ++i) {
                for (size_t j = i + 1; j < 4; ++j) {
                  if (!transTracks[i].impactPointTSCP().isValid() || 
                      !transTracks[j].impactPointTSCP().isValid()) continue;
                  
                  FreeTrajectoryState state1 = transTracks[i].impactPointTSCP().theState();
                  FreeTrajectoryState state2 = transTracks[j].impactPointTSCP().theState();
                  
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
              std::array<reco::TrackRef, 4> trackRefs{{ref1, ref2, ref3, ref4}};
              for (size_t i = 0; i < 4; ++i) {
                const auto& trackRef = trackRefs[i];
                double dzvtx = trackRef->dz(primaryVtx);
                double dxyvtx = trackRef->dxy(primaryVtx);
                double dzerror = sqrt(trackRef->dzError() * trackRef->dzError() + zVtxError * zVtxError);
                double dxyerror = sqrt(trackRef->d0Error() * trackRef->d0Error() + xVtxError * yVtxError);
                
                double dauLongImpactSig = (dzerror > 0) ? dzvtx / dzerror : 0.0;
                double dauTransImpactSig = (dxyerror > 0) ? dxyvtx / dxyerror : 0.0;
                
                std::string suffix = std::to_string(i + 1);
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

            if (storeEventShape_) {
              chi->addUserFloat("acoplanarity", acoplanarity);
              chi->addUserFloat("sphericity", eventShape.sphericity);
              chi->addUserFloat("pca_lambda1", eventShape.eigenvalues[0]);
              chi->addUserFloat("pca_lambda2", eventShape.eigenvalues[1]);
              chi->addUserFloat("pca_lambda3", eventShape.eigenvalues[2]);
            }

            outputs[stateIdx]->push_back(*chi);
          }
        }
      }
    }
  }

  for (std::size_t idx = 0; idx < outputs.size(); ++idx) {
    event.put(std::move(outputs[idx]), states_[idx].name);
  }
}

void ChiCFourTrackProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCFourTrackProducer);
