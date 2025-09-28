// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCResonanceProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCResonanceProducer.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include <algorithm>
#include <iterator>
#include <memory>
#include <set>

ChiCResonanceProducer::ChiCResonanceProducer(const edm::ParameterSet& cfg)
  : resonanceToken_(consumes<ResonanceCollection>(cfg.getParameter<edm::InputTag>("resonanceCollection"))),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    requireUniqueTracks_(cfg.getParameter<bool>("requireUniqueTracks"))
{
  const auto& statePsets = cfg.getParameter<std::vector<edm::ParameterSet> >("states");
  states_.reserve(statePsets.size());
  for (const auto& ps : statePsets) {
    ChiStateConfig state;
    state.name = ps.getParameter<std::string>("name");
    state.pdgId = ps.getParameter<int>("pdgId");
    state.mass = ps.getParameter<double>("mass");
    state.massWindow = ps.getParameter<double>("massWindow");
    states_.push_back(state);
    produces<ChiCollection>(state.name);
  }
}

ChiCResonanceProducer::~ChiCResonanceProducer() = default;

void ChiCResonanceProducer::beginJob() {}

void ChiCResonanceProducer::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<ResonanceCollection> resonances;
  event.getByToken(resonanceToken_, resonances);

  if (!resonances.isValid() || resonances->size() < 2)
    return;

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  const auto& coll = *resonances;
  AddFourMomenta addP4;

  for (std::size_t i = 0; i < coll.size(); ++i) {
    const auto& first = coll[i];

    for (std::size_t j = i + 1; j < coll.size(); ++j) {
      const auto& second = coll[j];

      if (requireUniqueTracks_ && shareTracks(first, second))
        continue;

      const double vx = 0.5 * (first.vx() + second.vx());
      const double vy = 0.5 * (first.vy() + second.vy());
      const double vz = 0.5 * (first.vz() + second.vz());
      const reco::Candidate::Point vertex(vx, vy, vz);

      for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
        const auto& state = states_[stateIdx];
        auto chi = std::make_unique<pat::CompositeCandidate>();
        chi->setPdgId(state.pdgId);
        chi->setCharge(first.charge() + second.charge());
        chi->setVertex(vertex);

        chi->addDaughter(first, "Resonance1");
        chi->addDaughter(second, "Resonance2");

        addP4.set(*chi);
        const double mass = chi->mass();
        if (applyMassWindow_) {
          if (mass < state.mass - state.massWindow || mass > state.mass + state.massWindow)
            continue;
        }

        outputs[stateIdx]->push_back(*chi);
      }
    }
  }

  for (std::size_t idx = 0; idx < outputs.size(); ++idx) {
    event.put(std::move(outputs[idx]), states_[idx].name);
  }
}

void ChiCResonanceProducer::endJob() {}

bool ChiCResonanceProducer::shareTracks(const reco::VertexCompositeCandidate& first,
                                        const reco::VertexCompositeCandidate& second) const {
  std::set<reco::TrackRef> tracks;
  for (size_t idx = 0; idx < first.numberOfDaughters(); ++idx) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(first.daughter(idx));
    if (dau && dau->track().isNonnull()) {
      tracks.insert(dau->track());
    }
  }

  for (size_t idx = 0; idx < second.numberOfDaughters(); ++idx) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(second.daughter(idx));
    if (dau && dau->track().isNonnull() && tracks.count(dau->track())) {
      return true;
    }
  }

  return false;
}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCResonanceProducer);
