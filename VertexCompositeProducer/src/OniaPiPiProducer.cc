#include <memory>

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/OniaPiPiProducer.h"

// Constructor
OniaPiPiProducer::OniaPiPiProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
//  useAnyMVA_ = false;
//  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("B");
//   produces< pat::CompositeCandidateCollection >("B");
  produces< std::vector<int> >("OniaIdx");
//  if(useAnyMVA_) produces<MVACollection>("MVAValuesB");
}

// (Empty) Destructor
OniaPiPiProducer::~OniaPiPiProducer() {
}


//
// Methods
//

// Producer Method
void OniaPiPiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create BFitter object which reconstructs the vertices and creates
//   BFitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     bCandidates( new reco::VertexCompositeCandidateCollection );
//
   auto bCandidates = std::make_unique<reco::VertexCompositeCandidateCollection>();
   // auto bCandidates = std::make_unique<pat::CompositeCandidateCollection>();
   auto oniaIndices = std::make_unique<std::vector<int> >(theVees.getOniaIdx().begin(), theVees.getOniaIdx().end());
   bCandidates->reserve( theVees.getB().size() );

   std::copy( theVees.getB().begin(),
              theVees.getB().end(),
              std::back_inserter(*bCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(bCandidates), std::string("B") );
   iEvent.put( std::move(oniaIndices), std::string("OniaIdx") );
/*    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesB"));
   }
*/
   theVees.resetAll();
}


//void OniaPiPiProducer::beginJob() {
void OniaPiPiProducer::beginJob() {
}


void OniaPiPiProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(OniaPiPiProducer);
