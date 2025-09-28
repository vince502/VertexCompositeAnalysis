// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiKaonProducer

#ifndef VertexCompositeAnalysis__DI_KAON_PRODUCER_H
#define VertexCompositeAnalysis__DI_KAON_PRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiKaonFitter.h"

class DiKaonProducer : public edm::one::EDProducer<> {
public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  explicit DiKaonProducer(const edm::ParameterSet&);
  ~DiKaonProducer() override;

private:
  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  DiKaonFitter fitter_;
};

#endif
