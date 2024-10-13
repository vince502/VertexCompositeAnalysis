// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniaPiPiProducer
// 
/**\class OniaPiPiProducer OniaPiPiProducer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/OniaPiPiProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis_OniaPiPi_PRODUCER_H
#define VertexCompositeAnalysis_OniaPiPi_PRODUCER_H

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/OniapipiFitter.h"

class OniaPiPiProducer : public edm::one::EDProducer<> {
public:
//  using MVACollection = std::vector<float>;

  explicit OniaPiPiProducer(const edm::ParameterSet&);
  ~OniaPiPiProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//  bool useAnyMVA_;

  OniapipiFitter theVees; 
//  edm::ParameterSet theParams;
};

#endif
