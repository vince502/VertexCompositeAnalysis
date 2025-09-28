// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCFourTrackProducer

#ifndef VertexCompositeAnalysis__CHIC_FOUR_TRACK_PRODUCER_H
#define VertexCompositeAnalysis__CHIC_FOUR_TRACK_PRODUCER_H

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

class ChiCFourTrackProducer : public edm::one::EDProducer<> {
public:
  using ChiCollection = pat::CompositeCandidateCollection;

  explicit ChiCFourTrackProducer(const edm::ParameterSet&);
  ~ChiCFourTrackProducer() override;

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
  std::vector<double> daughterMasses_;
  double minTrackPt_;
  double maxTrackEta_;
  double maxTrackChi2_;
  int minTrackNHits_;
  bool applyMassWindow_;
  double minCandidatePt_;
  double minAcoplanarity_;
  double maxSphericity_;

  std::vector<ChiStateConfig> states_;
};

#endif
