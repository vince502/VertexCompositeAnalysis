// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiC2Producer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiC2Producer.h"

#include <algorithm>
#include <iterator>
#include <string>

ChiC2Producer::ChiC2Producer(const edm::ParameterSet& cfg)
  : fitter_(cfg, consumesCollector())
{
  produces<ChiCollection>("ChiC2");
}

ChiC2Producer::~ChiC2Producer() = default;

void ChiC2Producer::beginJob() {}

void ChiC2Producer::produce(edm::Event& event, const edm::EventSetup& setup) {
  fitter_.fitAll(event, setup);

  auto output = std::make_unique<ChiCollection>();
  const auto& chiCandidates = fitter_.getChiC2();
  output->reserve(chiCandidates.size());
  std::copy(chiCandidates.begin(), chiCandidates.end(), std::back_inserter(*output));

  event.put(std::move(output), std::string("ChiC2"));

  fitter_.resetAll();
}

void ChiC2Producer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiC2Producer);
