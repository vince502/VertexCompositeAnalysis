// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiMuonFromTracksProducer
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuonFromTracksProducer.h"

#include <algorithm>
#include <iterator>
#include <string>

DiMuonFromTracksProducer::DiMuonFromTracksProducer(const edm::ParameterSet& cfg)
  : fitter_(cfg, consumesCollector())
{
  produces<CandidateCollection>("Jpsi");
}

DiMuonFromTracksProducer::~DiMuonFromTracksProducer() = default;

void DiMuonFromTracksProducer::beginJob() {}

void DiMuonFromTracksProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  fitter_.fitAll(event, setup);

  auto output = std::make_unique<CandidateCollection>();
  const auto& diMuons = fitter_.getDiMuons();
  output->reserve(diMuons.size());
  std::copy(diMuons.begin(), diMuons.end(), std::back_inserter(*output));

  event.put(std::move(output), std::string("Jpsi"));

  fitter_.resetAll();
}

void DiMuonFromTracksProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(DiMuonFromTracksProducer);
