// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiMuonFromTracksProducer
//

#ifndef VertexCompositeAnalysis__DI_MUON_FROM_TRACKS_PRODUCER_H
#define VertexCompositeAnalysis__DI_MUON_FROM_TRACKS_PRODUCER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuonFromTracksFitter.h"

class DiMuonFromTracksProducer : public edm::one::EDProducer<> {
public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  explicit DiMuonFromTracksProducer(const edm::ParameterSet&);
  ~DiMuonFromTracksProducer() override;

private:
  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  DiMuonFromTracksFitter fitter_;
};

#endif
