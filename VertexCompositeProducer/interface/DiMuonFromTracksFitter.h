// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiMuonFromTracksFitter
//
/**\class DiMuonFromTracksFitter
 *
 *  Fitter for J/ψ(μ+μ-) reconstruction from general tracks with muon mass hypothesis
 *  Similar to V0Fitter but specifically for dimuons
 */

#ifndef VertexCompositeAnalysis__DI_MUON_FROM_TRACKS_FITTER_H
#define VertexCompositeAnalysis__DI_MUON_FROM_TRACKS_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include <vector>

class MagneticField;

class DiMuonFromTracksFitter {
 public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  DiMuonFromTracksFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC);
  ~DiMuonFromTracksFitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const CandidateCollection& getDiMuons() const { return theDiMuons_; }

  void resetAll();

 private:
  CandidateCollection theDiMuons_;

  edm::EDGetTokenT<reco::TrackCollection> token_tracks_;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices_;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  // Cuts
  double tkChi2Cut_;
  int tkNhitsCut_;
  double tkPtCut_;
  double tkEtaCut_;
  double tkDCACut_;
  double mllCutMin_;
  double mllCutMax_;
  double jpsiMassCut_;
  double chi2Cut_;
  double dauTransImpactSigCut_;
  double dauLongImpactSigCut_;
  bool doVertexFit_;
  std::string vtxFitter_;
};

#endif
