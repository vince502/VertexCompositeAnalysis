// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      KshortProducer

#ifndef VertexCompositeAnalysis__KSHORT_PRODUCER_H
#define VertexCompositeAnalysis__KSHORT_PRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/KshortFitter.h"

class KshortProducer : public edm::one::EDProducer<> {
public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  explicit KshortProducer(const edm::ParameterSet&);
  ~KshortProducer() override;

private:
  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  KshortFitter fitter_;
};

#endif
