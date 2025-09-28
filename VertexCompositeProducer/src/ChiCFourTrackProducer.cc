// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiCFourTrackProducer

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ChiCFourTrackProducer.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>

namespace {
reco::RecoChargedCandidate makeRecoCandidate(const reco::TrackRef& trackRef, double mass) {
  const auto& track = *trackRef;
  const double momentum2 = track.momentum().mag2();
  const double energy = std::sqrt(momentum2 + mass * mass);
  math::XYZTLorentzVector p4(track.px(), track.py(), track.pz(), energy);
  reco::RecoChargedCandidate cand(trackRef->charge(), p4, track.vertex());
  cand.setTrack(trackRef);
  return cand;
}

bool hasTwoPositiveTwoNegative(const std::array<int, 4>& charges) {
  int nPos = 0;
  int nNeg = 0;
  for (auto q : charges) {
    if (q > 0)
      ++nPos;
    else if (q < 0)
      ++nNeg;
  }
  return (nPos == 2 && nNeg == 2);
}

struct EventShapeResult {
  double sphericity{0.0};
  std::array<double, 3> eigenvalues{{0.0, 0.0, 0.0}};
};

EventShapeResult computeEventShape(const std::array<reco::RecoChargedCandidate, 4>& daughters) {
  EventShapeResult result;
  Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
  double sumP2 = 0.0;
  for (const auto& dau : daughters) {
    Eigen::Vector3d p(dau.px(), dau.py(), dau.pz());
    tensor += p * p.transpose();
    sumP2 += p.squaredNorm();
  }

  if (sumP2 <= 0.0)
    return result;

  tensor /= sumP2;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(tensor);
  if (eigenSolver.info() != Eigen::Success)
    return result;

  result.eigenvalues = {{eigenSolver.eigenvalues()(0), eigenSolver.eigenvalues()(1), eigenSolver.eigenvalues()(2)}};
  std::sort(result.eigenvalues.begin(), result.eigenvalues.end(), std::greater<double>());
  result.sphericity = 1.5 * (result.eigenvalues[1] + result.eigenvalues[2]);
  return result;
}

double computeAcoplanarity(const std::array<reco::RecoChargedCandidate, 4>& daughters) {
  std::array<std::pair<double, std::size_t>, 4> ptIndex;
  for (std::size_t idx = 0; idx < daughters.size(); ++idx) {
    ptIndex[idx] = std::make_pair(daughters[idx].pt(), idx);
  }
  std::sort(ptIndex.begin(), ptIndex.end(), [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

  const auto& lead = daughters[ptIndex[0].second];
  const auto& sublead = daughters[ptIndex[1].second];
  const double deltaPhi = reco::deltaPhi(lead.phi(), sublead.phi());
  return 1.0 - std::abs(deltaPhi) / M_PI;
}
}

ChiCFourTrackProducer::ChiCFourTrackProducer(const edm::ParameterSet& cfg)
  : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackCollection"))),
    daughterMasses_(cfg.getParameter<std::vector<double> >("daughterMasses")),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    maxTrackChi2_(cfg.getParameter<double>("maxTrackNormalizedChi2")),
    minTrackNHits_(cfg.getParameter<int>("minTrackNHits")),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    minCandidatePt_(cfg.getParameter<double>("minCandidatePt")),
    minAcoplanarity_(cfg.getParameter<double>("minAcoplanarity")),
    maxSphericity_(cfg.getParameter<double>("maxSphericity")),
    maxCandidateAbsEta_(cfg.getParameter<double>("maxCandidateAbsEta")),
    storeEventShape_(cfg.getParameter<bool>("storeEventShape"))
{
  if (daughterMasses_.empty()) {
    throw cms::Exception("InvalidConfiguration") << "Parameter 'daughterMasses' must contain at least one value.";
  }

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

ChiCFourTrackProducer::~ChiCFourTrackProducer() = default;

void ChiCFourTrackProducer::beginJob() {}

void ChiCFourTrackProducer::produce(edm::Event& event, const edm::EventSetup&) {
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

  if (selectedTracks.size() < 4)
    return;

  std::vector<std::unique_ptr<ChiCollection> > outputs;
  outputs.reserve(states_.size());
  for (std::size_t idx = 0; idx < states_.size(); ++idx) {
    outputs.push_back(std::make_unique<ChiCollection>());
  }

  AddFourMomenta addP4;
  const bool applyAcoplanarityCut = (minAcoplanarity_ > 0.0);
  const bool applySphericityCut = (maxSphericity_ >= 0.0 && maxSphericity_ < std::numeric_limits<double>::infinity());

  const auto massForIndex = [this](unsigned int index) {
    if (daughterMasses_.size() == 1)
      return daughterMasses_.front();
    if (index < daughterMasses_.size())
      return daughterMasses_[index];
    return daughterMasses_.back();
  };

  const std::size_t nTracks = selectedTracks.size();
  for (std::size_t i = 0; i < nTracks - 3; ++i) {
    const auto& ref1 = selectedTracks[i];
    for (std::size_t j = i + 1; j < nTracks - 2; ++j) {
      const auto& ref2 = selectedTracks[j];
      for (std::size_t k = j + 1; k < nTracks - 1; ++k) {
        const auto& ref3 = selectedTracks[k];
        for (std::size_t l = k + 1; l < nTracks; ++l) {
          const auto& ref4 = selectedTracks[l];

          const std::array<int, 4> charges{{ref1->charge(), ref2->charge(), ref3->charge(), ref4->charge()}};
          if (!hasTwoPositiveTwoNegative(charges))
            continue;

          const std::array<reco::RecoChargedCandidate, 4> daughters{{
              makeRecoCandidate(ref1, massForIndex(0)),
              makeRecoCandidate(ref2, massForIndex(1)),
              makeRecoCandidate(ref3, massForIndex(2)),
              makeRecoCandidate(ref4, massForIndex(3))}};

          for (std::size_t stateIdx = 0; stateIdx < states_.size(); ++stateIdx) {
            const auto& state = states_[stateIdx];
            auto chi = std::make_unique<pat::CompositeCandidate>();
            chi->setPdgId(state.pdgId);
            chi->setCharge(0);

            const auto& v1 = daughters[0].vertex();
            const auto& v2 = daughters[1].vertex();
            const auto& v3 = daughters[2].vertex();
            const auto& v4 = daughters[3].vertex();
            const double vx = 0.25 * (v1.x() + v2.x() + v3.x() + v4.x());
            const double vy = 0.25 * (v1.y() + v2.y() + v3.y() + v4.y());
            const double vz = 0.25 * (v1.z() + v2.z() + v3.z() + v4.z());
            chi->setVertex(reco::Candidate::Point(vx, vy, vz));

            chi->addDaughter(daughters[0], "Track1");
            chi->addDaughter(daughters[1], "Track2");
            chi->addDaughter(daughters[2], "Track3");
            chi->addDaughter(daughters[3], "Track4");

            addP4.set(*chi);

            if (chi->pt() < minCandidatePt_)
              continue;

            const double mass = chi->mass();
            if (applyMassWindow_) {
              if (mass < state.mass - state.massWindow || mass > state.mass + state.massWindow)
                continue;
            }

            const double absEta = std::abs(chi->eta());
            if (maxCandidateAbsEta_ >= 0.0 && absEta > maxCandidateAbsEta_)
              continue;

            double acoplanarity = 0.0;
            EventShapeResult eventShape;

            if (applyAcoplanarityCut || applySphericityCut || storeEventShape_) {
              acoplanarity = computeAcoplanarity(daughters);
              eventShape = computeEventShape(daughters);

              if (applyAcoplanarityCut && acoplanarity < minAcoplanarity_)
                continue;

              if (applySphericityCut && eventShape.sphericity > maxSphericity_)
                continue;
            }

            if (storeEventShape_) {
              chi->addUserFloat("acoplanarity", acoplanarity);
              chi->addUserFloat("sphericity", eventShape.sphericity);
              chi->addUserFloat("pca_lambda1", eventShape.eigenvalues[0]);
              chi->addUserFloat("pca_lambda2", eventShape.eigenvalues[1]);
              chi->addUserFloat("pca_lambda3", eventShape.eigenvalues[2]);
            }

            outputs[stateIdx]->push_back(*chi);
          }
        }
      }
    }
  }

  for (std::size_t idx = 0; idx < outputs.size(); ++idx) {
    event.put(std::move(outputs[idx]), states_[idx].name);
  }
}

void ChiCFourTrackProducer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ChiCFourTrackProducer);
