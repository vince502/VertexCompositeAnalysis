// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCTrackPairProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCTrackPairProducer.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <memory>

namespace {
reco::RecoChargedCandidate makeRecoCandidate(const reco::TrackRef& trackRef, double mass) {
  math::XYZTLorentzVector p4;
  const auto& track = *trackRef;
  const double momentum = track.p();
  const double energy = std::sqrt(momentum * momentum + mass * mass);
  p4.SetPxPyPzE(track.px(), track.py(), track.pz(), energy);
  reco::RecoChargedCandidate cand(trackRef->charge(), p4, track.vertex());
  cand.setTrack(trackRef);
  return cand;
}
}

ChiCTrackPairProducer::ChiCTrackPairProducer(const edm::ParameterSet& cfg)
  : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackCollection"))),
    daughterMass_(cfg.getParameter<double>("daughterMass")),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    maxTrackChi2_(cfg.getParameter<double>("maxTrackNormalizedChi2")),
    minTrackNHits_(cfg.getParameter<int>("minTrackNHits")),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    minPairPt_(cfg.getParameter<double>("minPairPt")),
    requiredChargeProduct_(cfg.existsAs<int>("requiredChargeProduct") ? cfg.getParameter<int>("requiredChargeProduct") : -1)
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

ChiCTrackPairProducer::~ChiCTrackPairProducer() = default;

void ChiCTrackPairProducer::beginJob() {}

void ChiCTrackPairProducer::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<reco::TrackCollection> tracks;
  event.getByToken(trackToken_, tracks);
  if (!tracks.isValid())
    return;

  std::vector<reco::TrackRef> selectedTracks;
  selectedTracks.reserve(tracks->size());

  for (std::size_t idx = 0; idx < tracks->size(); ++idx) {
    reco::TrackRef trackRef(tracks, idx);
    const auto& track = *trackRef;
    if (track.pt() < minTrackPt_)
      continue;
    if (std::abs(track.eta()) > maxTrackEta_)
      continue;
    if (track.normalizedChi2() > maxTrackChi2_)
      continue;
    if (track.numberOfValidHits() < minTrackNHits_)
      continue;
    selectedTracks.push_back(trackRef);
  }

  if (selectedTracks.size() < 2)
    return;

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  AddFourMomenta addP4;

  for (std::size_t i = 0; i < selectedTracks.size(); ++i) {
    const auto& trackRef1 = selectedTracks[i];
    const auto& track1 = *trackRef1;

    for (std::size_t j = i + 1; j < selectedTracks.size(); ++j) {
      const auto& trackRef2 = selectedTracks[j];
      const auto& track2 = *trackRef2;

      const int chargeProduct = track1.charge() * track2.charge();
      if (requiredChargeProduct_ == -1 && chargeProduct >= 0)
        continue;  // Require opposite charge
      if (requiredChargeProduct_ == +1 && chargeProduct <= 0)
        continue;  // Require same charge
      // If requiredChargeProduct_ == 0, accept any charge combination

      if ((track1.pt() + track2.pt()) < minPairPt_)
        continue;

      auto dau1 = makeRecoCandidate(trackRef1, daughterMass_);
      auto dau2 = makeRecoCandidate(trackRef2, daughterMass_);

      for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
        const auto& state = states_[stateIdx];
        auto chi = std::make_unique<pat::CompositeCandidate>();
        chi->setPdgId(state.pdgId);
        chi->setCharge(dau1.charge() + dau2.charge());
        chi->setVertex(reco::Candidate::Point((track1.vx() + track2.vx()) * 0.5,
                                              (track1.vy() + track2.vy()) * 0.5,
                                              (track1.vz() + track2.vz()) * 0.5));

        chi->addDaughter(dau1, "Track1");
        chi->addDaughter(dau2, "Track2");

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

void ChiCTrackPairProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCTrackPairProducer);
