// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCResonanceProducer

#ifndef VertexCompositeAnalysis__CHIC_RESONANCE_PRODUCER_H
#define VertexCompositeAnalysis__CHIC_RESONANCE_PRODUCER_H

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

class ChiCResonanceProducer : public edm::one::EDProducer<> {
public:
  using ResonanceCollection = reco::VertexCompositeCandidateCollection;
  using ChiCollection = pat::CompositeCandidateCollection;

  explicit ChiCResonanceProducer(const edm::ParameterSet&);
  ~ChiCResonanceProducer() override;

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

  bool shareTracks(const reco::VertexCompositeCandidate& first,
                   const reco::VertexCompositeCandidate& second) const;

  edm::EDGetTokenT<ResonanceCollection> resonanceToken_;
  std::vector<ChiStateConfig> states_;
  bool applyMassWindow_;
  bool requireUniqueTracks_;
};

#endif
