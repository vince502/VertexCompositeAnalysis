#ifndef __PVANAMULTIFS__
#define __PVANAMULTIFS__

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <math.h>

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


class PrimaryVertexAnalyzerMultiFS : public edm::EDAnalyzer {

public:
    explicit PrimaryVertexAnalyzerMultiFS(const edm::ParameterSet &); 
    ~PrimaryVertexAnalyzerMultiFS();

private:
    virtual void beginJob(){};
    virtual void analyze(const edm::Event &, const edm::EventSetup &);
    virtual void endJob(){};

    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;

    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;

    bool doGen;
    bool doReco;
};

PrimaryVertexAnalyzerMultiFS::PrimaryVertexAnalyzerMultiFS(const edm::ParameterSet &iConfig){
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    doGen = iConfig.getParameter<bool>("doGen");
    doReco = iConfig.getParameter<bool>("doReco");
};

PrimaryVertexAnalyzerMultiFS::~PrimaryVertexAnalyzerMultiFS(){

};

DEFINE_FWK_MODULE(PrimaryVertexAnalyzerMultiFS);
#endif