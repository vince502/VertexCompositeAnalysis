// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCTrackPairProducer

#ifndef VertexCompositeAnalysis__CHIC_TRACK_PAIR_PRODUCER_H
#define VertexCompositeAnalysis__CHIC_TRACK_PAIR_PRODUCER_H

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class ChiCTrackPairProducer : public edm::one::EDProducer<> {
public:
  using ChiCollection = pat::CompositeCandidateCollection;

  explicit ChiCTrackPairProducer(const edm::ParameterSet&);
  ~ChiCTrackPairProducer() override;

private:
  struct ChiStateConfig {
    std::string name;
    int pdgId;
    double mass;
    double massWindow;
  };

  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  double daughterMass_;
  double minTrackPt_;
  double maxTrackEta_;
  double maxTrackChi2_;
  int minTrackNHits_;
  bool applyMassWindow_;
  double minPairPt_;
  int requiredChargeProduct_;  // -1 for opposite charge, +1 for same charge, 0 for any

  std::vector<ChiStateConfig> states_;
};

#endif
