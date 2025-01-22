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
  produces< std::vector<float> >("DCAValCollection1");
  produces< std::vector<float> >("DCAErrCollection1");
  produces< std::vector<float> >("Angle2D1");
  produces< std::vector<float> >("Angle3D1");
  produces< std::vector<float> >("DCAValCollection2");
  produces< std::vector<float> >("DCAErrCollection2");
  produces< std::vector<float> >("Angle2D2");
  produces< std::vector<float> >("Angle3D2");
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
   auto dcaval1 = std::make_unique<std::vector<float>>();
   auto dcaerr1 = std::make_unique<std::vector<float>>();
   auto angle2D1 = std::make_unique<std::vector<float>>();
   auto angle3D1 = std::make_unique<std::vector<float>>();
   auto dcaval2 = std::make_unique<std::vector<float>>();
   auto dcaerr2 = std::make_unique<std::vector<float>>();
   auto angle2D2 = std::make_unique<std::vector<float>>();
   auto angle3D2 = std::make_unique<std::vector<float>>();
   d0Candidates->reserve( theVees.getDD().size() );
   dcaval1->reserve( theVees.getDcaVal1().size() );
   dcaerr1->reserve( theVees.getDcaErr1().size() );
   angle2D1->reserve( theVees.getAngle2D1().size() );
   angle3D1->reserve( theVees.getAngle3D1().size() );
   dcaval2->reserve( theVees.getDcaVal2().size() );
   dcaerr2->reserve( theVees.getDcaErr2().size() );
   angle2D2->reserve( theVees.getAngle2D2().size() );
   angle3D2->reserve( theVees.getAngle3D2().size() );

   std::copy( theVees.getDD().begin(), theVees.getDD().end(), std::back_inserter(*d0Candidates) );
   std::copy( theVees.getDcaVal1().begin(), theVees.getDcaVal1().end(), std::back_inserter(*dcaval1));
   std::copy( theVees.getDcaErr1().begin(), theVees.getDcaErr1().end(), std::back_inserter(*dcaerr1));
   std::copy( theVees.getAngle2D1().begin(), theVees.getAngle2D1().end(), std::back_inserter(*angle2D1));
   std::copy( theVees.getAngle3D1().begin(), theVees.getAngle3D1().end(), std::back_inserter(*angle3D1));
   std::copy( theVees.getDcaVal2().begin(), theVees.getDcaVal2().end(), std::back_inserter(*dcaval2));
   std::copy( theVees.getDcaErr2().begin(), theVees.getDcaErr2().end(), std::back_inserter(*dcaerr2));
   std::copy( theVees.getAngle2D2().begin(), theVees.getAngle2D2().end(), std::back_inserter(*angle2D2));
   std::copy( theVees.getAngle3D2().begin(), theVees.getAngle3D2().end(), std::back_inserter(*angle3D2));

   // Write the collections to the Event
   iEvent.put( std::move(d0Candidates), std::string("DD") );
   iEvent.put( std::move(dcaval1), std::string("DCAValCollection1") );
   iEvent.put( std::move(dcaerr1), std::string("DCAErrCollection1") );
   iEvent.put( std::move(angle2D1), std::string("Angle2D1") );
   iEvent.put( std::move(angle3D1), std::string("Angle3D1") );

   iEvent.put( std::move(dcaval2), std::string("DCAValCollection2") );
   iEvent.put( std::move(dcaerr2), std::string("DCAErrCollection2") );
   iEvent.put( std::move(angle2D2), std::string("Angle2D2") );
   iEvent.put( std::move(angle3D2), std::string("Angle3D2") );
    
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
