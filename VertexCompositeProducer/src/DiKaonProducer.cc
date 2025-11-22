// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiKaonProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiKaonProducer.h"

#include <algorithm>
#include <iterator>
#include <string>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

DiKaonProducer::DiKaonProducer(const edm::ParameterSet& cfg)
  : fitter_(cfg, consumesCollector())
{
  produces<CandidateCollection>("DiKaon");
}

DiKaonProducer::~DiKaonProducer() = default;

void DiKaonProducer::beginJob() {}

void DiKaonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::LogInfo("DiKaonProducer") << "Processing event " << event.id().run() << ":" << event.id().luminosityBlock() << ":" << event.id().event();
  
  fitter_.fitAll(event, setup);

  auto output = std::make_unique<CandidateCollection>();
  const auto& diKaons = fitter_.getDiKaons();
  
  edm::LogInfo("DiKaonProducer") << "Found " << diKaons.size() << " di-kaon candidates";
  
  output->reserve(diKaons.size());
  std::copy(diKaons.begin(), diKaons.end(), std::back_inserter(*output));
  
  if (diKaons.size() > 0) {
    edm::LogInfo("DiKaonProducer") << "Storing " << output->size() << " di-kaon candidates to event";
    for (size_t i = 0; i < std::min(diKaons.size(), size_t(3)); ++i) {
      const auto& dk = diKaons[i];
      edm::LogInfo("DiKaonProducer") << "  DiKaon[" << i << "]: mass=" << dk.mass() 
                                     << " GeV, pt=" << dk.pt() << " GeV, eta=" << dk.eta()
                                     << ", daughters=" << dk.numberOfDaughters();
    }
  } else {
    edm::LogWarning("DiKaonProducer") << "No di-kaon candidates found in event " << event.id().run() << ":" << event.id().luminosityBlock() << ":" << event.id().event();
  }

  event.put(std::move(output), std::string("DiKaon"));

  fitter_.resetAll();
}

void DiKaonProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(DiKaonProducer);
