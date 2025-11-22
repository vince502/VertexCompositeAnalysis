// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiMuonFromTracksFitter
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuonFromTracksFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <cmath>
#include <algorithm>
#include <vector>

using reco::TransientTrack;

const double muonMass = 0.1056583745;
const double muonMassSquared = muonMass * muonMass;
const double jpsiMass = 3.096916;  // J/ψ mass

DiMuonFromTracksFitter::DiMuonFromTracksFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC)
  : token_tracks_(iC.consumes<reco::TrackCollection>(params.getParameter<edm::InputTag>("trackRecoAlgorithm"))),
    token_vertices_(iC.consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("vertexRecoAlgorithm"))),
    token_beamSpot_(iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
    bFieldToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    tkChi2Cut_(params.getParameter<double>("tkChi2Cut")),
    tkNhitsCut_(params.getParameter<int>("tkNhitsCut")),
    tkPtCut_(params.getParameter<double>("tkPtCut")),
    tkEtaCut_(params.getParameter<double>("tkEtaCut")),
    tkDCACut_(params.getParameter<double>("tkDCACut")),
    mllCutMin_(params.getParameter<double>("mllCutMin")),
    mllCutMax_(params.getParameter<double>("mllCutMax")),
    jpsiMassCut_(params.getParameter<double>("jpsiMassCut")),
    chi2Cut_(params.getParameter<double>("vtxChi2Cut")),
    dauTransImpactSigCut_(params.getParameter<double>("dauTransImpactSigCut")),
    dauLongImpactSigCut_(params.getParameter<double>("dauLongImpactSigCut")),
    doVertexFit_(params.getParameter<bool>("doVertexFit")),
    vtxFitter_(params.getParameter<std::string>("vertexFitter"))
{
}

DiMuonFromTracksFitter::~DiMuonFromTracksFitter() = default;

void DiMuonFromTracksFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  theDiMuons_.clear();

  edm::Handle<reco::TrackCollection> trackHandle;
  edm::Handle<reco::VertexCollection> vertexHandle;
  edm::Handle<reco::BeamSpot> beamSpotHandle;

  iEvent.getByToken(token_tracks_, trackHandle);
  iEvent.getByToken(token_vertices_, vertexHandle);
  iEvent.getByToken(token_beamSpot_, beamSpotHandle);

  if (!trackHandle.isValid() || trackHandle->empty()) return;

  const MagneticField* magField = &iSetup.getData(bFieldToken_);

  // Get primary vertex
  bool isVtxPV = false;
  double xVtx = beamSpotHandle->position().x();
  double yVtx = beamSpotHandle->position().y();
  double zVtx = 0.0;
  double xVtxError = beamSpotHandle->BeamWidthX();
  double yVtxError = beamSpotHandle->BeamWidthY();
  double zVtxError = 0.0;

  if (vertexHandle.isValid() && !vertexHandle->empty()) {
    const auto& vtxPrimary = vertexHandle->front();
    if (!vtxPrimary.isFake() && vtxPrimary.tracksSize() >= 2) {
      isVtxPV = true;
      xVtx = vtxPrimary.x();
      yVtx = vtxPrimary.y();
      zVtx = vtxPrimary.z();
      xVtxError = vtxPrimary.xError();
      yVtxError = vtxPrimary.yError();
      zVtxError = vtxPrimary.zError();
    }
  }

  // Select tracks
  std::vector<reco::TrackRef> trackRefs;
  std::vector<reco::TransientTrack> transTracks;

  for (size_t i = 0; i < trackHandle->size(); ++i) {
    reco::TrackRef trackRef(trackHandle, i);
    if (trackRef->normalizedChi2() >= tkChi2Cut_) continue;
    if (trackRef->numberOfValidHits() < tkNhitsCut_) continue;
    if (trackRef->pt() <= tkPtCut_) continue;
    if (std::abs(trackRef->eta()) > tkEtaCut_) continue;

    reco::TransientTrack transTrack(*trackRef, magField);
    math::XYZPoint bestvtx(xVtx, yVtx, zVtx);
    double dzvtx = trackRef->dz(bestvtx);
    double dxyvtx = trackRef->dxy(bestvtx);
    double dzerror = std::sqrt(trackRef->dzError() * trackRef->dzError() + zVtxError * zVtxError);
    double dxyerror = std::sqrt(trackRef->d0Error() * trackRef->d0Error() + xVtxError * yVtxError);

    double dauLongImpactSig = dzvtx / dzerror;
    double dauTransImpactSig = dxyvtx / dxyerror;

    if (std::abs(dauTransImpactSig) > dauTransImpactSigCut_ && std::abs(dauLongImpactSig) > dauLongImpactSigCut_) {
      trackRefs.push_back(trackRef);
      transTracks.push_back(transTrack);
    }
  }

  // Loop over track pairs
  for (size_t i = 0; i < trackRefs.size(); ++i) {
    for (size_t j = i + 1; j < trackRefs.size(); ++j) {
      // Require opposite charge
      if (trackRefs[i]->charge() * trackRefs[j]->charge() >= 0) continue;

      reco::TrackRef posTrack, negTrack;
      reco::TransientTrack* posTrans = nullptr;
      reco::TransientTrack* negTrans = nullptr;

      if (trackRefs[i]->charge() > 0) {
        posTrack = trackRefs[i];
        negTrack = trackRefs[j];
        posTrans = &transTracks[i];
        negTrans = &transTracks[j];
      } else {
        posTrack = trackRefs[j];
        negTrack = trackRefs[i];
        posTrans = &transTracks[j];
        negTrans = &transTracks[i];
      }

      // Calculate DCA
      FreeTrajectoryState posState = posTrans->impactPointTSCP().theState();
      FreeTrajectoryState negState = negTrans->impactPointTSCP().theState();

      if (!posTrans->impactPointTSCP().isValid() || !negTrans->impactPointTSCP().isValid()) continue;

      ClosestApproachInRPhi cApp;
      cApp.calculate(posState, negState);
      if (!cApp.status()) continue;

      float dca = std::abs(cApp.distance());
      GlobalPoint cxPt = cApp.crossingPoint();

      if (dca < 0.0 || dca > tkDCACut_) continue;
      if (std::sqrt(cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y()) > 120.0 || std::abs(cxPt.z()) > 300.0) continue;

      // Get trajectory states at POCA
      TrajectoryStateClosestToPoint posTSCP = posTrans->trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint negTSCP = negTrans->trajectoryStateClosestToPoint(cxPt);

      if (!posTSCP.isValid() || !negTSCP.isValid()) continue;

      // Calculate invariant mass with muon hypothesis
      double totalE = std::sqrt(posTSCP.momentum().mag2() + muonMassSquared) +
                      std::sqrt(negTSCP.momentum().mag2() + muonMassSquared);
      double totalESq = totalE * totalE;
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
      double mass = std::sqrt(totalESq - totalPSq);

      if (mass < mllCutMin_ || mass > mllCutMax_) continue;

      // Vertex fit
      std::vector<reco::TransientTrack> fitTracksVec;
      fitTracksVec.push_back(*posTrans);
      fitTracksVec.push_back(*negTrans);

      TransientVertex recoVertex;
      if (doVertexFit_ && vtxFitter_ == "KalmanVertexFitter") {
        KalmanVertexFitter fitter(false);
        recoVertex = fitter.vertex(fitTracksVec);
      } else {
        continue;  // Require vertex fit
      }

      if (!recoVertex.isValid() || recoVertex.totalChiSquared() < 0.0) continue;
      if (recoVertex.totalChiSquared() / recoVertex.degreesOfFreedom() > chi2Cut_) continue;

      // Create candidate
      reco::Vertex vtx = recoVertex;
      GlobalVector momentum = posTSCP.momentum() + negTSCP.momentum();
      double energy = std::sqrt(momentum.mag2() + jpsiMass * jpsiMass);
      reco::Particle::LorentzVector p4(momentum.x(), momentum.y(), momentum.z(), energy);

      reco::VertexCompositeCandidate* jpsi = new reco::VertexCompositeCandidate(0, p4, vtx.position(), vtx.covariance(), vtx.chi2(), vtx.ndof());

      reco::RecoChargedCandidate posCand(posTrack->charge(), reco::Particle::LorentzVector(posTSCP.momentum().x(), posTSCP.momentum().y(), posTSCP.momentum().z(), std::sqrt(posTSCP.momentum().mag2() + muonMassSquared)), vtx.position());
      posCand.setTrack(posTrack);
      reco::RecoChargedCandidate negCand(negTrack->charge(), reco::Particle::LorentzVector(negTSCP.momentum().x(), negTSCP.momentum().y(), negTSCP.momentum().z(), std::sqrt(negTSCP.momentum().mag2() + muonMassSquared)), vtx.position());
      negCand.setTrack(negTrack);

      jpsi->addDaughter(posCand);
      jpsi->addDaughter(negCand);
      jpsi->setPdgId(443);  // J/ψ

      AddFourMomenta addP4;
      addP4.set(*jpsi);

      // Apply J/ψ mass cut
      if (std::abs(jpsi->mass() - jpsiMass) < jpsiMassCut_) {
        theDiMuons_.push_back(*jpsi);
      }
      delete jpsi;
    }
  }
}

void DiMuonFromTracksFitter::resetAll() {
  theDiMuons_.clear();
}
