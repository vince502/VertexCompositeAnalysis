// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiKaonFitter
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiKaonFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

DiKaonFitter::DiKaonFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC)
  : baseFitter_(params, std::move(iC))
{
  if (!params.getParameter<bool>("selectPhis")) {
    edm::LogWarning("DiKaonFitter")
        << "Parameter 'selectPhis' is false. The di-kaon producer relies on V0Fitter"
           " phi selection; please enable it in the configuration.";
  }
}

DiKaonFitter::~DiKaonFitter() = default;

void DiKaonFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::LogInfo("DiKaonFitter") << "DiKaonFitter::fitAll - Calling baseFitter_.fitAll() for event " 
                                << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event();
  baseFitter_.fitAll(iEvent, iSetup);
  const auto& phis = baseFitter_.getPhis();
  edm::LogInfo("DiKaonFitter") << "DiKaonFitter::fitAll - After baseFitter_.fitAll(), found " << phis.size() << " phi candidates";
}

const DiKaonFitter::CandidateCollection& DiKaonFitter::getDiKaons() const {
  return baseFitter_.getPhis();
}

void DiKaonFitter::resetAll() {
  baseFitter_.resetAll();
}
