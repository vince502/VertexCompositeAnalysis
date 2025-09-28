// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiC2Fitter

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiC2Fitter.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <iterator>
#include <set>

ChiC2Fitter::ChiC2Fitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC)
  : diKaonToken_(iC.consumes<DiKaonCollection>(params.getParameter<edm::InputTag>("diKaonCollection"))),
    chiMass_(params.getParameter<double>("chiMass")),
    chiMassWindow_(params.getParameter<double>("chiMassWindow")),
    applyMassWindow_(params.getParameter<bool>("applyMassWindow"))
{}

ChiC2Fitter::~ChiC2Fitter() = default;

void ChiC2Fitter::fitAll(const edm::Event& event, const edm::EventSetup&) {
  edm::Handle<DiKaonCollection> diKaons;
  event.getByToken(diKaonToken_, diKaons);

  if (!diKaons.isValid() || diKaons->size() < 2)
    return;

  const auto& diKaonColl = *diKaons;
  const double minMass = chiMass_ - chiMassWindow_;
  const double maxMass = chiMass_ + chiMassWindow_;

  AddFourMomenta addP4;

  for (std::size_t i = 0; i < diKaonColl.size(); ++i) {
    const auto& first = diKaonColl[i];

    for (std::size_t j = i + 1; j < diKaonColl.size(); ++j) {
      const auto& second = diKaonColl[j];

      if (shareTracks(first, second))
        continue;

      pat::CompositeCandidate chiCandidate;
      chiCandidate.setPdgId(445); // chi_c2(1P)
      chiCandidate.setCharge(0);

      chiCandidate.addDaughter(first, "DiKaon1");
      chiCandidate.addDaughter(second, "DiKaon2");

      const double vx = 0.5 * (first.vx() + second.vx());
      const double vy = 0.5 * (first.vy() + second.vy());
      const double vz = 0.5 * (first.vz() + second.vz());
      chiCandidate.setVertex(reco::Candidate::Point(vx, vy, vz));

      addP4.set(chiCandidate);

      const double mass = chiCandidate.mass();
      if (applyMassWindow_ && (mass < minMass || mass > maxMass))
        continue;

      chiCandidates_.push_back(chiCandidate);
    }
  }
}

const ChiC2Fitter::ChiCollection& ChiC2Fitter::getChiC2() const {
  return chiCandidates_;
}

void ChiC2Fitter::resetAll() {
  chiCandidates_.clear();
}

bool ChiC2Fitter::shareTracks(const reco::VertexCompositeCandidate& first,
                              const reco::VertexCompositeCandidate& second) const {
  std::set<reco::TrackRef> tracks;
  for (std::size_t i = 0; i < first.numberOfDaughters(); ++i) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(first.daughter(i));
    if (dau && dau->track().isNonnull()) {
      tracks.insert(dau->track());
    }
  }

  for (std::size_t i = 0; i < second.numberOfDaughters(); ++i) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(second.daughter(i));
    if (dau && dau->track().isNonnull() && tracks.count(dau->track())) {
      return true;
    }
  }

  return false;
}
