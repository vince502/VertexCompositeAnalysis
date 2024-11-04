// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm> 



class GenParicleSimpleAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns>{
public:
    explicit GenParicleSimpleAnalyzer(const edm::ParameterSet &);
    ~GenParicleSimpleAnalyzer(){};
private:
    virtual void analyze(const edm::Event &, const edm::EventSetup &);
    virtual void beginJob(){};
    virtual void beginRun(const edm::Run&, const edm::EventSetup&){};
    virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
    virtual void endJob() {};


    std::vector<int> wanted;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    
};

GenParicleSimpleAnalyzer::GenParicleSimpleAnalyzer(const edm::ParameterSet& iConfig):
    wanted(iConfig.getUntrackedParameter<std::vector<int> >("pdgIDs")),
    tok_genParticle_ (consumes<reco::GenParticleCollection>( edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection"))))
//    tok_pgenParticle_ (consumes<pat::PackedGenParticleCollection>( edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("PackedGenParticleCollection"))))
{
};


//void GenParicleSimpleAnalyzer::beginRun(){};
//void GenParicleSimpleAnalyzer::endRun(){};
//void GenParicleSimpleAnalyzer::beginJob(){};
//void GenParicleSimpleAnalyzer::endJob(){};

void GenParicleSimpleAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSteup){
    using namespace std;
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_, genpars);
    if(!genpars.isValid()){
        cout << "Faulty gen collection.. I'm out" << endl;
        return;
    }

    auto const nGens = genpars->size();
    cout << "nGens : " << nGens << endl;
    for(unsigned int igen=0; igen < nGens; igen++){
        const reco::GenParticle &trk = (*genpars)[igen];
        auto hasP = std::find(wanted.begin(), wanted.end(), trk.pdgId() );
        if(hasP == wanted.end()) continue;
        cout << trk.pdgId();
        cout << Form(" -> %zu (",trk.numberOfDaughters());
        string dauInfo = "    Dau : ";
        for(unsigned int idau=0; idau <trk.numberOfDaughters(); ++idau){
            auto const& dau = trk.daughter(idau);
            cout << dau->pdgId() << " ";
            dauInfo += Form("(%.3f, %.3f) ", dau->pt(), dau->eta());
            string dauInfo2 = "        Grand Dau : ";
	    
	    for(unsigned int idau2=0; idau2 < dau->numberOfDaughters(); ++idau2){
                auto const& dau2 = dau->daughter(idau2);
                cout << dau2->pdgId() << " ";
                dauInfo2 += Form("(%.3f, %.3f) ", dau2->pt(), dau2->eta());
	    }
            cout << "    )" << endl;
            cout << dauInfo2 << endl;
        }
        cout << ")" << endl;
        cout << dauInfo << endl;
    }
};

DEFINE_FWK_MODULE(GenParicleSimpleAnalyzer);
