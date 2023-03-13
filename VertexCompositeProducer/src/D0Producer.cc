// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      D0Producer
// 
/**\class D0Producer D0Producer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0Producer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//


// system include files
#include <memory>

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Producer.h"

// Constructor
D0Producer::D0Producer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("D0");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesD0");
  produces<std::vector<float > >("DCAValuesD0");
  produces<std::vector<float > >("DCAErrorsD0");
}

// (Empty) Destructor
D0Producer::~D0Producer() {
}


//
// Methods
//

// Producer Method
void D0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create D0Fitter object which reconstructs the vertices and creates
//   D0Fitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     d0Candidates( new reco::VertexCompositeCandidateCollection );
//
   auto d0Candidates = std::make_unique<reco::VertexCompositeCandidateCollection>();
   d0Candidates->reserve( theVees.getD0().size() );

   std::copy( theVees.getD0().begin(),
              theVees.getD0().end(),
              std::back_inserter(*d0Candidates) );

   // Write the collections to the Event
   iEvent.put( std::move(d0Candidates), std::string("D0") );
    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesD0"));
   }
   auto dcaVals = std::make_unique<std::vector<float > >(theVees.getDCAVals().begin(), theVees.getDCAVals().end());
   iEvent.put(std::move(dcaVals), std::string("DCAValuesD0"));
   auto dcaErrs = std::make_unique<std::vector<float > >(theVees.getDCAErrs().begin(), theVees.getDCAErrs().end());
   iEvent.put(std::move(dcaErrs), std::string("DCAErrorsD0"));

   theVees.resetAll();
}


//void D0Producer::beginJob() {
void D0Producer::beginJob() {
}


void D0Producer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(D0Producer);
