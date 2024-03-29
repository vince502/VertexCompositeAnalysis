// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      DDroducer
// 
/**\class DDProducer DDProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/DDProducer.cc

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DDProducer.h"

// Constructor
DDProducer::DDProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("DD");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesDD1");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesDD2");
  produces<std::vector<float > >("DCAValuesDD");
  produces<std::vector<float > >("DCAErrorsDD");
}

// (Empty) Destructor
DDProducer::~DDProducer() {
}


//
// Methods
//

// Producer Method
void DDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create DDFitter object which reconstructs the vertices and creates
//   DDFitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     d0Candidates( new reco::VertexCompositeCandidateCollection );
//
   auto d0Candidates = std::make_unique<reco::VertexCompositeCandidateCollection>();
   d0Candidates->reserve( theVees.getDD().size() );

   std::copy( theVees.getDD().begin(),
              theVees.getDD().end(),
              std::back_inserter(*d0Candidates) );

   // Write the collections to the Event
   iEvent.put( std::move(d0Candidates), std::string("DD") );
    
   if(useAnyMVA_) 
   {
     auto mvas1 = std::make_unique<MVACollection>(theVees.getMVAVals1().begin(),theVees.getMVAVals1().end());
     iEvent.put(std::move(mvas1), std::string("MVAValuesDD1"));
     auto mvas2 = std::make_unique<MVACollection>(theVees.getMVAVals2().begin(),theVees.getMVAVals2().end());
     iEvent.put(std::move(mvas2), std::string("MVAValuesDD2"));
   }
   auto dcaVals = std::make_unique<std::vector<float > >(theVees.getDCAVals().begin(), theVees.getDCAVals().end());
   iEvent.put(std::move(dcaVals), std::string("DCAValuesDD"));
   auto dcaErrs = std::make_unique<std::vector<float > >(theVees.getDCAErrs().begin(), theVees.getDCAErrs().end());
   iEvent.put(std::move(dcaErrs), std::string("DCAErrorsDD"));

   theVees.resetAll();
}


//void DDProducer::beginJob() {
void DDProducer::beginJob() {
}


void DDProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(DDProducer);
