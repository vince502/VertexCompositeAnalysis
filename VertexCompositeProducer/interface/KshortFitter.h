// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      KshortFitter
//
/**\class KshortFitter
 *  Lightweight access to V0Fitter restricted to Kshort reconstruction.
 */

#ifndef VertexCompositeAnalysis__KSHORT_FITTER_H
#define VertexCompositeAnalysis__KSHORT_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/V0Fitter.h"

class KshortFitter {
 public:
  using CandidateCollection = reco::VertexCompositeCandidateCollection;

  KshortFitter(const edm::ParameterSet& params, edm::ConsumesCollector&& iC);
  ~KshortFitter();

  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  const CandidateCollection& getKshorts() const;

  void resetAll();

 private:
  V0Fitter baseFitter_;
};

#endif
