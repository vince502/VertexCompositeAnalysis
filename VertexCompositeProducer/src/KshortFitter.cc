// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      KshortFitter

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/KshortFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

KshortFitter::KshortFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC)
  : baseFitter_(params, std::move(iC))
{
  if (!params.getParameter<bool>("selectKshorts")) {
    edm::LogWarning("KshortFitter")
        << "Parameter 'selectKshorts' is false. Enable it to populate Kshort candidates.";
  }
}

KshortFitter::~KshortFitter() = default;

void KshortFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  baseFitter_.fitAll(iEvent, iSetup);
}

const KshortFitter::CandidateCollection& KshortFitter::getKshorts() const {
  return baseFitter_.getKshorts();
}

void KshortFitter::resetAll() {
  baseFitter_.resetAll();
}
