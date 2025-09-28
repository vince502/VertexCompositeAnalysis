// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      KshortProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/KshortProducer.h"

#include <algorithm>
#include <iterator>
#include <string>

KshortProducer::KshortProducer(const edm::ParameterSet& cfg)
  : fitter_(cfg, consumesCollector())
{
  produces<CandidateCollection>("Kshort");
}

KshortProducer::~KshortProducer() = default;

void KshortProducer::beginJob() {}

void KshortProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  fitter_.fitAll(event, setup);

  auto output = std::make_unique<CandidateCollection>();
  const auto& kshorts = fitter_.getKshorts();
  output->reserve(kshorts.size());
  std::copy(kshorts.begin(), kshorts.end(), std::back_inserter(*output));

  event.put(std::move(output), std::string("Kshort"));

  fitter_.resetAll();
}

void KshortProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(KshortProducer);
