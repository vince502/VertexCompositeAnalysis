// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniapipiFitterPAT
// 
/**\class OniapipiFitterPAT OniapipiFitterPAT.h VertexCompositeAnalysis/VertexCompositeProducer/interface/OniapipiFitterPAT.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

#ifndef VertexCompositeAnalysis__OTT_FITTER_H
#define VertexCompositeAnalysis__OTT_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"


#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/GBRForest/interface/GBRForest.h"

#include <string>
#include <fstream>
#include <typeinfo>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>

class OniapipiFitterPAT {
 public:
  OniapipiFitterPAT(const edm::ParameterSet& theParams, edm::ConsumesCollector && iC);
  ~OniapipiFitterPAT();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const pat::CompositeCandidateCollection& getB() const;
  const std::vector<int>& getOniaIdx() const;
//  const std::vector<float>& getMVAVals() const; 

  void resetAll();

 private:
  const float piMassB = 0.13957018;
  const float piMassBSquared = piMassB*piMassB;
  const float kaonMassB = 0.493677;
  const float kaonMassBSquared = kaonMassB*kaonMassB;
  const float oniaMass = 3.094;
  float piMassB_sigma = 3.5E-4f;
  float kaonMassB_sigma = 1.6E-4f;
  pat::CompositeCandidateCollection theBs;
  std::vector<int> oniaIdx;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bField_esToken_;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::InputTag d0Alg;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> token_onias;
  edm::EDGetTokenT<reco::TrackCollection> token_pixeltracks;
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;

 // Track option
  bool   usePixelTracks;
  // Cuts
  double mPiDCutMin;
  double mPiDCutMax;
  double batTkDCACut;
  double batTkChi2Cut;
  int    batTkNhitsCut;
  double batTkPtErrCut;
  double batTkPtCut;
  double batTkEtaCut;
  double batTrkEtaDiffCut;
  double batTkPtSumCut;
  double batTkEtaDiffCut;
  double bVtxChi2Cut;
  double bRVtxCut;
  double bRVtxSigCut;
  double bLVtxCut;
  double bLVtxSigCut;
  double bCollinCut2D;
  double bCollinCut3D;
  double bMassCut;
  double bQMassCut;
  std::vector<double> bOniaMass;
  std::vector<double> bOniaWindow;
  double batDauTransImpactSigCut;
  double batDauLongImpactSigCut;
  double bVtxChiProbCut;
  double bPtCut;
  double bAlphaCut;
  double bAlpha2DCut;
  bool   isWrongSignB;

  double trk1cosPhiCut;
  double trk2cosPhiCut;
  double trk12cosPhiCut;
  double trk1dRCut;
  double trk2dRCut;
  double trk12dRCut;
  double oniapTCut;
  double trk1pTCut;
  double trk2pTCut;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  //setup mva selector
/*
  bool useAnyMVA_;
  std::vector<bool> useMVA_;
  std::vector<double> min_MVA_;
  std::string mvaType_;
  std::string forestLabel_;
  GBRForest * forest_;
  bool useForestFromDB_;

  std::vector<float> mvaVals_;
*/
//  auto_ptr<edm::ValueMap<float> >mvaValValueMap;
//  MVACollection mvas; 

//  std::string dbFileName_;

  bool match_OniaPiPi(const reco::Track* a, const reco::Track* b){
          if( a->charge() != b->charge() ) return false;
          if(
                  fabs(a->pt()- b->pt()) < 0.03 &&
                  reco::deltaR(*a, *b) < 0.03
          ) return true;
          return false;
  };
  double Dphi(double phi1, double phi2) {
      double dphi = phi1 - phi2;
      return std::atan2(std::sin(dphi), std::cos(dphi)); // Ensures result is in [-π, π]
  };

  };

#endif
