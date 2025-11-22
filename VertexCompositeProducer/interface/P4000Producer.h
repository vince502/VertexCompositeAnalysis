// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      P4000Producer
//
// Producer for P(4000) -> J/ψ(μ+μ-) + φ(K+K-)

#ifndef VertexCompositeAnalysis__P4000_PRODUCER_H
#define VertexCompositeAnalysis__P4000_PRODUCER_H

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

class MagneticField;

class P4000Producer : public edm::one::EDProducer<> {
public:
  using ResonanceCollection = reco::VertexCompositeCandidateCollection;
  using P4000Collection = pat::CompositeCandidateCollection;

  explicit P4000Producer(const edm::ParameterSet&);
  ~P4000Producer() override;

private:
  struct P4000StateConfig {
    std::string name;
    int pdgId;
    double mass;
    double massWindow;
  };

  void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  bool shareTracks(const reco::VertexCompositeCandidate& jpsi,
                   const reco::VertexCompositeCandidate& phi) const;

  edm::EDGetTokenT<ResonanceCollection> jpsiToken_;
  edm::EDGetTokenT<ResonanceCollection> phiToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  std::vector<double> resonanceMassSigmas_;
  bool applyMassWindow_;
  bool requireUniqueTracks_;
  bool useVertexFitting_;

  std::vector<P4000StateConfig> states_;
};

#endif
