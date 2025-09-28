// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiKaonProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiKaonProducer.h"

#include <algorithm>
#include <iterator>
#include <string>

DiKaonProducer::DiKaonProducer(const edm::ParameterSet& cfg)
  : fitter_(cfg, consumesCollector())
{
  produces<CandidateCollection>("DiKaon");
}

DiKaonProducer::~DiKaonProducer() = default;

void DiKaonProducer::beginJob() {}

void DiKaonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  fitter_.fitAll(event, setup);

  auto output = std::make_unique<CandidateCollection>();
  const auto& diKaons = fitter_.getDiKaons();
  output->reserve(diKaons.size());
  std::copy(diKaons.begin(), diKaons.end(), std::back_inserter(*output));

  event.put(std::move(output), std::string("DiKaon"));

  fitter_.resetAll();
}

void DiKaonProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(DiKaonProducer);
