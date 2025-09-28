// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiC2Producer

#ifndef VertexCompositeAnalysis__CHIC2_PRODUCER_H
#define VertexCompositeAnalysis__CHIC2_PRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiC2Fitter.h"

class ChiC2Producer : public edm::one::EDProducer<> {
public:
  using ChiCollection = pat::CompositeCandidateCollection;

  explicit ChiC2Producer(const edm::ParameterSet&);
  ~ChiC2Producer() override;

private:
  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  ChiC2Fitter fitter_;
};

#endif
