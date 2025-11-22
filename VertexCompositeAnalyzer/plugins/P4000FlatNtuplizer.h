// -*- C++ -*-
//
// Package:    VertexCompositeAnalyzer
// Class:      P4000FlatNtuplizer
//
// Flat ntuplizer for P(4000) -> J/ψ(μ+μ-) + φ(K+K-)

#ifndef VertexCompositeAnalysis__P4000_FLAT_NTUPLIZER_H
#define VertexCompositeAnalysis__P4000_FLAT_NTUPLIZER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TTree.h"
#include "TFile.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>
#include <array>

class P4000FlatNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit P4000FlatNtuplizer(const edm::ParameterSet&);
  ~P4000FlatNtuplizer() override;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void resetBranches();
  void fillCandidate(const pat::CompositeCandidate& cand, const reco::Vertex* primaryVertex);

  edm::EDGetTokenT<pat::CompositeCandidateCollection> p4000Token_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;

  edm::Service<TFileService> fileService_;
  TTree* tree_;

  // Event-level
  unsigned int run_;
  unsigned int lumi_;
  unsigned long long event_;

  // P(4000) candidate
  float p4000Mass_;
  float p4000Pt_;
  float p4000Eta_;
  float p4000Phi_;
  float p4000Rapidity_;
  float p4000Vx_;
  float p4000Vy_;
  float p4000Vz_;
  float p4000VtxChi2_;
  float p4000VtxNdof_;
  float p4000VtxProb_;

  // J/ψ from P(4000)
  float jpsiMass_;
  float jpsiPt_;
  float jpsiEta_;
  float jpsiPhi_;
  float jpsiVx_;
  float jpsiVy_;
  float jpsiVz_;

  // φ from P(4000)
  float phiMass_;
  float phiPt_;
  float phiEta_;
  float phiPhi_;
  float phiVx_;
  float phiVy_;
  float phiVz_;

  // J/ψ daughters (muons)
  std::array<float, 2> muPt_;
  std::array<float, 2> muEta_;
  std::array<float, 2> muPhi_;
  std::array<int, 2> muCharge_;
  std::array<float, 2> muDxy_;
  std::array<float, 2> muDz_;

  // φ daughters (kaons)
  std::array<float, 2> kaonPt_;
  std::array<float, 2> kaonEta_;
  std::array<float, 2> kaonPhi_;
  std::array<int, 2> kaonCharge_;
  std::array<float, 2> kaonDxy_;
  std::array<float, 2> kaonDz_;

  // Primary vertex
  int nPV_;
  float pvX_;
  float pvY_;
  float pvZ_;
  int pvNTracks_;
};

#endif
