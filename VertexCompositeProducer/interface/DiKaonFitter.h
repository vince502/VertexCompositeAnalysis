// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DiKaonFitter
//
/**\class DiKaonFitter
 *
 *  Lightweight wrapper around V0Fitter to provide a dedicated
 *  di-kaon (phi-like) candidate builder. This module keeps the
 *  existing selection and vertexing logic from V0Fitter while
 *  exposing only the di-kaon collection needed for higher-level
 *  combinations.
 */

#ifndef VertexCompositeAnalysis__DI_KAON_FITTER_H
#define VertexCompositeAnalysis__DI_KAON_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/V0Fitter.h"

class DiKaonFitter {
 public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  DiKaonFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC);
  ~DiKaonFitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const CandidateCollection& getDiKaons() const;

  void resetAll();

 private:
  V0Fitter baseFitter_;
};

#endif
