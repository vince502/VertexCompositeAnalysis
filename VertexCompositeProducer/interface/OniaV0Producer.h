// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniaV0Producer
// 
/**\class OniaV0Producer OniaV0Producer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/OniaV0Producer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis_OniaV0_PRODUCER_H
#define VertexCompositeAnalysis_OniaV0_PRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/OniaV0Fitter.h"

class OniaV0Producer : public edm::one::EDProducer<> {
public:
//  using MVACollection = std::vector<float>;

  explicit OniaV0Producer(const edm::ParameterSet&);
  ~OniaV0Producer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//  bool useAnyMVA_;

  OniaV0Fitter theVees; 
//  edm::ParameterSet theParams;
};

#endif
