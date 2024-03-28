// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DDProducer
// 
/**\class DDProducer DDProducer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/DDProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis__DD_PRODUCER_H
#define VertexCompositeAnalysis__DD_PRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DDFitter.h"

class DDProducer : public edm::EDProducer {
public:
  using MVAPairCollection = std::vector<std::pair<float,float>>;

  explicit DDProducer(const edm::ParameterSet&);
  ~DDProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  bool useAnyMVA_;

  DDFitter theVees; 
//  edm::ParameterSet theParams;
};

#endif
