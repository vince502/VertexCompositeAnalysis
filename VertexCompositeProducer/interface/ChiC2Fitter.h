// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      ChiC2Fitter

#ifndef VertexCompositeAnalysis__CHIC2_FITTER_H
#define VertexCompositeAnalysis__CHIC2_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

class ChiC2Fitter {
 public:
  using DiKaonCollection = reco::VertexCompositeCandidateCollection;
  using ChiCollection = pat::CompositeCandidateCollection;

  ChiC2Fitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC);
  ~ChiC2Fitter();

  void fitAll(const edm::Event& event, const edm::EventSetup& setup);

  const ChiCollection& getChiC2() const;

  void resetAll();

 private:
  bool shareTracks(const reco::VertexCompositeCandidate& first,
                   const reco::VertexCompositeCandidate& second) const;

  edm::EDGetTokenT<DiKaonCollection> diKaonToken_;
  double chiMass_;
  double chiMassWindow_;
  bool applyMassWindow_;

  ChiCollection chiCandidates_;
};

#endif
