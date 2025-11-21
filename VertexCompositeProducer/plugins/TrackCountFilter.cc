// -*- C++ -*-
//
// Package:    VertexCompositeAnalysis/VertexCompositeProducer
// Class:      TrackCountFilter
// 
/**\class TrackCountFilter TrackCountFilter.cc VertexCompositeAnalysis/VertexCompositeProducer/plugins/TrackCountFilter.cc

 Description: Filter events based on number of selected tracks matching ChiC producer criteria

 Implementation:
     Applies track selection cuts (pt, eta, chi2, nhits) and filters based on count
*/
//

// system include files
#include <memory>
#include <limits>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

//
// class declaration
//

class TrackCountFilter : public edm::stream::EDFilter<> {
   public:
      explicit TrackCountFilter(const edm::ParameterSet&);
      ~TrackCountFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
      double minTrackPt_;
      double maxTrackEta_;
      double maxTrackChi2_;
      int minTrackNHits_;
      int minNSelectedTracks_;
      int maxNSelectedTracks_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackCountFilter::TrackCountFilter(const edm::ParameterSet& iConfig) :
  trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"))),
  minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
  maxTrackEta_(iConfig.getParameter<double>("maxTrackEta")),
  maxTrackChi2_(iConfig.getParameter<double>("maxTrackNormalizedChi2")),
  minTrackNHits_(iConfig.getParameter<int>("minTrackNHits")),
  minNSelectedTracks_(iConfig.getParameter<int>("minNSelectedTracks")),
  maxNSelectedTracks_(iConfig.getParameter<int>("maxNSelectedTracks"))
{
   //now do what ever initialization is needed
}


TrackCountFilter::~TrackCountFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrackCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  
  if (!tracks.isValid()) {
    LogDebug("TrackCountFilter") << "Track collection is invalid" << std::endl;
    return false;
  }

  int nSelected = 0;
  
  for (const auto& track : *tracks) {
    if (track.pt() < minTrackPt_)
      continue;
    if (std::abs(track.eta()) > maxTrackEta_)
      continue;
    if (track.normalizedChi2() > maxTrackChi2_)
      continue;
    if (track.numberOfValidHits() < minTrackNHits_)
      continue;
    nSelected++;
  }

  bool pass = (nSelected >= minNSelectedTracks_ && nSelected <= maxNSelectedTracks_);
  
  LogDebug("TrackCountFilter") << "N selected tracks: " << nSelected 
    << ", min: " << minNSelectedTracks_ << ", max: " << maxNSelectedTracks_
    << ", pass: " << (pass ? "true" : "false") << std::endl;

  return pass;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TrackCountFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TrackCountFilter::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackCollection", edm::InputTag("generalTracks"))->setComment("Track collection");
  desc.add<double>("minTrackPt", 0.4)->setComment("Minimum track pT");
  desc.add<double>("maxTrackEta", 2.4)->setComment("Maximum track |eta|");
  desc.add<double>("maxTrackNormalizedChi2", 10.0)->setComment("Maximum track normalized chi2");
  desc.add<int>("minTrackNHits", 6)->setComment("Minimum number of valid hits");
  desc.add<int>("minNSelectedTracks", 4)->setComment("Minimum number of selected tracks (>=)");
  desc.add<int>("maxNSelectedTracks", std::numeric_limits<int>::max())->setComment("Maximum number of selected tracks (<=)");
  descriptions.add("trackCountFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackCountFilter);
