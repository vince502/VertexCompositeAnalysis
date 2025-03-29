#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/VertexCompositeTreeProducer2.h"

//#define DEBUG false

#define DEBUG true
#define PI 3.1416
#define MAXCAN 1000000

VertexCompositeTreeProducer2::VertexCompositeTreeProducer2(const edm::ParameterSet& iConfig)
{
    //options
    doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
    doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
    twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
    threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
    doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
    doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
    hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
    decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
    doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
    doMuonFull_ = iConfig.getUntrackedParameter<bool>("doMuonFull");
    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    if(threeProngDecay_) PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");
    
    saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
    saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");
    saveAllHistogram_ = iConfig.getUntrackedParameter<bool>("saveAllHistogram");
    massHistPeak_ = iConfig.getUntrackedParameter<double>("massHistPeak");
    massHistWidth_ = iConfig.getUntrackedParameter<double>("massHistWidth");
    massHistBins_ = iConfig.getUntrackedParameter<int>("massHistBins");

    useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
    isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA"); 

    //cut variables
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);

    pTBins_ = iConfig.getUntrackedParameter< std::vector<double> >("pTBins");
    yBins_  = iConfig.getUntrackedParameter< std::vector<double> >("yBins");

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
    tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

    isCentrality_ = false;
    if(iConfig.exists("isCentrality")) isCentrality_ = iConfig.getParameter<bool>("isCentrality");
    if(isCentrality_)
    {
      tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
      tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
    }

    isEventPlane_ = false;
    if(iConfig.exists("isEventPlane")) isEventPlane_ = iConfig.getParameter<bool>("isEventPlane");
    if(isEventPlane_)
    {
      tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));
    }

    if(useAnyMVA_ && iConfig.exists("MVACollection"))
      MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
}


VertexCompositeTreeProducer2::~VertexCompositeTreeProducer2()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeTreeProducer2::analyze(const edm::Event& iEvent, const edm::EventSetup&
iSetup)
{
    using std::vector;
    using namespace edm;
    using namespace reco;

    if(doGenNtuple_) fillGEN(iEvent,iSetup);
    if(doRecoNtuple_) fillRECO(iEvent,iSetup);

    if(saveTree_) VertexCompositeNtuple->Fill();
    //clear vector;
    clearVectors();

}

void
VertexCompositeTreeProducer2::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef DEBUG
    using std::cout;
    using std::endl;
#endif
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    const reco::VertexCompositeCandidateCollection * v0candidates_ = v0candidates.product();
    
    edm::Handle<MVACollection> mvavalues;
    if(useAnyMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
    }

    edm::Handle<reco::GenParticleCollection> genpars;
    if(doGenMatching_ || doGenMatchingTOF_) iEvent.getByToken(tok_genParticle_,genpars);

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    centrality=-1;
    if(isCentrality_)
    {
      edm::Handle<reco::Centrality> cent;
      iEvent.getByToken(tok_centSrc_, cent);

      iEvent.getByToken(tok_centBinLabel_,cbin_);
      centrality = *cbin_;  

      HFsumETPlus = cent->EtHFtowerSumPlus();
      HFsumETMinus = cent->EtHFtowerSumMinus();
      Npixel = cent->multiplicityPixel();
      ZDCPlus = cent->zdcSumPlus();
      ZDCMinus = cent->zdcSumMinus();
    }

    if(isEventPlane_)
    {
      edm::Handle<reco::EvtPlaneCollection> eventplanes;
      iEvent.getByToken(tok_eventplaneSrc_,eventplanes);

      const reco::EvtPlane & ephfp1 = (*eventplanes)[0];
      const reco::EvtPlane & ephfm1 = (*eventplanes)[1];
      const reco::EvtPlane & ephfp2 = (*eventplanes)[6];
      const reco::EvtPlane & ephfm2 = (*eventplanes)[7];
      const reco::EvtPlane & ephfp3 = (*eventplanes)[13];
      const reco::EvtPlane & ephfm3 = (*eventplanes)[14];
     
      ephfpAngle[0] = ephfp1.angle(2);
      ephfpAngle[1] = ephfp2.angle(2);
      ephfpAngle[2] = ephfp3.angle(2);

      ephfmAngle[0] = ephfm1.angle(2);
      ephfmAngle[1] = ephfm2.angle(2);
      ephfmAngle[2] = ephfm3.angle(2);

      ephfpQ[0] = ephfp1.q(2);
      ephfpQ[1] = ephfp2.q(2);
      ephfpQ[2] = ephfp3.q(2);

      ephfmQ[0] = ephfm1.q(2);
      ephfmQ[1] = ephfm2.q(2);
      ephfmQ[2] = ephfm3.q(2);

      ephfpSumW = ephfp2.sumw();
      ephfmSumW = ephfm2.sumw();
    }

    //best vertex
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    //Ntrkoffline
    Ntrkoffline = 0;
    if(multMax_!=-1 && multMin_!=-1)
    {
      for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        Ntrkoffline++;
      }
    }

    std::vector<reco::GenParticleRef> genRefs;
    if(doGenMatching_){
        if(!genpars.isValid())
        { cout<<"Gen matching cannot be done without Gen collection!!"<<endl; return; }
        for(unsigned int it=0; it<genpars->size(); ++it){
            const reco::GenParticle & trk = (*genpars)[it];
            int id = trk.pdgId();
            if(fabs(id)!=PID_) continue; //check is target
            if(decayInGen_ && trk.numberOfDaughters()!=2 && !threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
            if(decayInGen_ && trk.numberOfDaughters()!=3 && threeProngDecay_) continue; //check 2-pron decay if target decays in Gen

            // wrong when considering two layer decay
            int nDau = threeProngDecay_ ? 3 : 2;
            std::vector<unsigned int> idxs;
            std::vector<unsigned int> permutations(nDau);
            std::iota(permutations.begin(), permutations.end(), 0);
            std::sort(permutations.begin(), permutations.end());
            if (!threeProngDecay_) {
              do {
                auto Dd1 = trk.daughter( permutations.at(0) );
                // if(fabs(Dd1->pdgId())!=421){cout << "wrongmatching: " << Dd1->pdgId() << endl;}
                auto Dd2 = trk.daughter( permutations.at(1) );
                // if(fabs(Dd2->pdgId())==421){cout << "wrongmatching: " << Dd2->pdgId() << " Dd1? " << Dd1->pdgId() << endl;}
                if (abs(Dd1->pdgId()) == PID_dau1_ && abs(Dd2->pdgId()) == PID_dau2_) {
                  if(twoLayerDecay_){
                    // Magic numbers, _permutations -> number of D0 daughters;
                    std::vector<unsigned int> _permutations(2);
                    std::iota(_permutations.begin(), _permutations.end(), 0);
                    std::sort(_permutations.begin(), _permutations.end());
                    do {
                      auto Ddd1 = Dd1->daughter( _permutations.at(0) );
                      auto Ddd2 = Dd1->daughter( _permutations.at(1) );
                      if (abs(Ddd1->pdgId()) == 211 && abs(Ddd2->pdgId()) == 321) {
                        idxs = permutations;
                        break;
                      }
                    } while (std::next_permutation(_permutations.begin(), _permutations.end()));
                    if(!idxs.empty()) break;
                  } else {
                    if (abs(Dd1->pdgId()) == PID_dau1_
                        && abs(Dd2->pdgId()) == PID_dau2_
                        ) {
                      idxs = permutations;
                      break;
                    }// while (std::next_permutation(permutations.begin(), permutations.end()));
                    if(!idxs.empty()) break;
                  }
                }
              } while (std::next_permutation(permutations.begin(), permutations.end()));
            } else {
              do {
                auto Dd1 = trk.daughter( permutations.at(0) );
                auto Dd2 = trk.daughter( permutations.at(1) );
                auto Dd3 = trk.daughter( permutations.at(2) );

                if (abs(Dd1->pdgId()) == PID_dau1_
                    && abs(Dd2->pdgId()) == PID_dau2_
                    && abs(Dd3->pdgId() == PID_dau3_)) {
                  idxs = permutations;
                  break;
                }
              } while (std::next_permutation(permutations.begin(), permutations.end()));
            }
            if (decayInGen_ && idxs.empty()) continue;
            genRefs.push_back(reco::GenParticleRef(genpars, it));
        }
        //if (genRefs.size()>1) std::cout << "More than one target of generated particles\n";
    }

    //RECO Candidate info
    candSize = v0candidates_->size();
    #ifdef DEBUG
    cout << "candSize : " << candSize << endl;
    #endif

    for(unsigned it=0; it<v0candidates_->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();

        eta.push_back(trk.eta());
        y.push_back(trk.rapidity());
        pt.push_back(trk.pt());
        #ifdef DEBUG
        #endif
        phi.push_back(trk.phi());
        flavor.push_back(trk.pdgId()/abs(trk.pdgId()));
        mva.push_back(0.0);
        if(useAnyMVA_) mva.push_back((*mvavalues)[it]);

        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass.push_back(trk.mass());
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * gd1;
        const reco::Candidate * gd2;
        if(twoLayerDecay_){
          gd1 = d1->daughter(0);
          gd2 = d1->daughter(1);
        }
        const reco::Candidate * d2 = trk.daughter(1);
        const reco::Candidate * d3 = 0;        
        if(threeProngDecay_) d3 = trk.daughter(2);

        //Gen match
        if(doGenMatching_ )
        {
          if( twoLayerDecay_ ){
            bool _matchGEN=false;
            unsigned int nGen = genRefs.size();
            bool _isSwap=false;
            int _idmom_reco=-77;
            int _idBAnc_reco=-77;

            for( unsigned int igen=0; igen<nGen; igen++){
              auto const theGenDStar = genRefs.at(igen);
              unsigned int idxD0 = 1;
              if( abs(theGenDStar->daughter(0)->pdgId()) == 421 ) idxD0 = 0;
              auto const* theGenD0 = genRefs.at(igen)->daughter(idxD0);
              auto const* theGenPion = genRefs.at(igen)->daughter(1- idxD0);
              // Only works for 2 body two layer decay
              reco::Candidate const* recoD1;
              reco::Candidate const* recoPi;
              unsigned int idxRecoD0 = 1;
              if (abs(trk.daughter(0)->pdgId())== 421) idxRecoD0 = 0;
              recoD1 = trk.daughter(idxRecoD0);
              recoPi = trk.daughter(1-idxRecoD0);
              const auto nGenDau = theGenD0->numberOfDaughters();

              matchGEN.push_back(_matchGEN || (matchHadron(recoD1, *theGenD0,true) && matchHadron(recoPi, *theGenPion,false)));
              #ifdef DEBUG
              cout << "matchGEN[igen]: " << matchGEN[igen] << endl;
              #endif
              if(matchGEN.back()){
                  isSwap.push_back(checkSwap(recoD1, *theGenD0));
                  auto mom_ref = findMother(theGenDStar);
                  if (mom_ref.isNonnull()) idmom_reco.push_back(mom_ref->pdgId());
                  else idmom_reco.push_back(_idmom_reco);
                  int __count_anc__ = 0;
                  auto __ref_anc__ = mom_ref;
                  if (!__ref_anc__.isNonnull()) idBAnc_reco.push_back(_idBAnc_reco);
                  while ( __ref_anc__.isNonnull() && __count_anc__ < 50 ){
                      __ref_anc__ = findMother(__ref_anc__);
                      if( __ref_anc__.isNonnull()){
                          if( ((int) abs(__ref_anc__->pdgId())) % 1000 / 100 == 5){ 
                              idBAnc_reco.push_back(__ref_anc__->pdgId());
                          } else {
                              idBAnc_reco.push_back(_idBAnc_reco);
                          }
                      }
                  }
              } else {
              isSwap.push_back(_isSwap);
              idmom_reco.push_back(_idmom_reco);
              idBAnc_reco.push_back(_idBAnc_reco);
              }

              matchGen_D0pT_.push_back(theGenD0->pt());
              matchGen_D0eta_.push_back(theGenD0->eta());
              matchGen_D0phi_.push_back(theGenD0->phi());
              matchGen_D0mass_.push_back(theGenD0->mass());
              matchGen_D0y_.push_back(theGenD0->rapidity());
              matchGen_D0charge_.push_back(theGenD0->charge());
              matchGen_D0pdgId_.push_back(theGenD0->pdgId());

              genDecayLength(*theGenD0, matchGen_D1decayLength2D_, matchGen_D1decayLength3D_, matchGen_D1angle2D_, matchGen_D1angle3D_ );
              getAncestorId(*theGenD0, matchGen_D1ancestorId_, matchGen_D1ancestorFlavor_ );

              const auto* genDau0 = theGenD0->daughter(0);
              const auto* genDau1 = theGenD0->daughter(1);

              matchGen_D0Dau1_pT_.push_back(genDau0->pt());
              matchGen_D0Dau1_eta_.push_back(genDau0->eta());
              matchGen_D0Dau1_phi_.push_back(genDau0->phi());
              matchGen_D0Dau1_mass_.push_back(genDau0->mass());
              matchGen_D0Dau1_y_.push_back(genDau0->rapidity());
              matchGen_D0Dau1_charge_.push_back(genDau0->charge());
              matchGen_D0Dau1_pdgId_.push_back(genDau0->pdgId());

              matchGen_D0Dau2_pT_.push_back(genDau1->pt());
              matchGen_D0Dau2_eta_.push_back(genDau1->eta());
              matchGen_D0Dau2_phi_.push_back(genDau1->phi());
              matchGen_D0Dau2_mass_.push_back(genDau1->mass());
              matchGen_D0Dau2_y_.push_back(genDau1->rapidity());
              matchGen_D0Dau2_charge_.push_back(genDau1->charge());
              matchGen_D0Dau2_pdgId_.push_back(genDau1->pdgId());

              matchGen_D1pT_.push_back(theGenPion->pt());
              matchGen_D1eta_.push_back(theGenPion->eta());
              matchGen_D1phi_.push_back(theGenPion->phi());
              matchGen_D1mass_.push_back(theGenPion->mass());
              matchGen_D1y_.push_back(theGenPion->rapidity());
              matchGen_D1charge_.push_back(theGenPion->charge());
              matchGen_D1pdgId_.push_back(theGenPion->pdgId());
              }
            } // END for nGen
         
            else {
            matchGEN.push_back(false);
            unsigned int nGen = genRefs.size();
            isSwap.push_back(false);
            idmom_reco.push_back(-77);
            idBAnc_reco.push_back(-77);

            for( unsigned int igen=0; igen<nGen; igen++){
              auto const theGenP = genRefs.at(igen);
              #ifdef DEBUG 
              cout << "theGenP pdgId : " << theGenP->pdgId() << endl;
              #endif
              matchGEN.push_back(matchGEN.at(it) || matchHadron(&trk, *theGenP,true));
              if(matchGEN.at(it)){
              isSwap.push_back(checkSwap(&trk, *theGenP));
              auto mom_ref = findMother(theGenP);
              if (mom_ref.isNonnull()) idmom_reco.push_back(mom_ref->pdgId());
              int __count_anc__ = 0;
              auto __ref_anc__ = mom_ref;
              while ( __ref_anc__.isNonnull() && __count_anc__ < 50 ){
                __ref_anc__ = findMother(__ref_anc__);
                if( __ref_anc__.isNonnull()){
                if( ((int) abs(__ref_anc__->pdgId())) % 1000 / 100 == 5){ 
                  idBAnc_reco.push_back(__ref_anc__->pdgId());
                } 
                } 
              }

              matchGen_D0pT_.push_back(theGenP->pt());
              matchGen_D0eta_.push_back(theGenP->eta());
              matchGen_D0phi_.push_back(theGenP->phi());
              matchGen_D0mass_.push_back(theGenP->mass());
              matchGen_D0y_.push_back(theGenP->rapidity());
              matchGen_D0charge_.push_back(theGenP->charge());
              matchGen_D0pdgId_.push_back(theGenP->pdgId());

              genDecayLength(*theGenP, matchGen_D1decayLength2D_.at(it), matchGen_D1decayLength3D_.at(it), matchGen_D1angle2D_.at(it), matchGen_D1angle3D_.at(it) );
              getAncestorId(*theGenP, matchGen_D1ancestorId_.at(it), matchGen_D1ancestorFlavor_.at(it) );

              const auto* genDau0 = theGenP->daughter(0);
              const auto* genDau1 = theGenP->daughter(1);

              matchGen_D0Dau1_pT_.push_back(genDau0->pt());
              matchGen_D0Dau1_eta_.push_back(genDau0->eta());
              matchGen_D0Dau1_phi_.push_back(genDau0->phi());
              matchGen_D0Dau1_mass_.push_back(genDau0->mass());
              matchGen_D0Dau1_y_.push_back(genDau0->rapidity());
              matchGen_D0Dau1_charge_.push_back(genDau0->charge());
              matchGen_D0Dau1_pdgId_.push_back(genDau0->pdgId());

              matchGen_D0Dau2_pT_.push_back(genDau1->pt());
              matchGen_D0Dau2_eta_.push_back(genDau1->eta());
              matchGen_D0Dau2_phi_.push_back(genDau1->phi());
              matchGen_D0Dau2_mass_.push_back(genDau1->mass());
              matchGen_D0Dau2_y_.push_back(genDau1->rapidity());
              matchGen_D0Dau2_charge_.push_back(genDau1->charge());
              matchGen_D0Dau2_pdgId_.push_back(genDau1->pdgId());
                }
              } // END for nGen
            }
          }
          
          double pxd1 = d1->px();
          double pyd1 = d1->py();
          double pzd1 = d1->pz();
          double pxd2 = d2->px();
          double pyd2 = d2->py();
          double pzd2 = d2->pz();
          
          TVector3 dauvec1(pxd1,pyd1,pzd1);
          TVector3 dauvec2(pxd2,pyd2,pzd2);
          
          //pt
          pt1.push_back(d1->pt());
          pt2.push_back(d2->pt());
          
          //momentum
          p1.push_back(d1->p());
          p2.push_back(d2->p());
          
          //eta
          eta1.push_back(d1->eta());
          eta2.push_back(d2->eta());
          
          //phi
          phi1.push_back(d1->phi());
          phi2.push_back(d2->phi());
          
          //charge
          charge1.push_back(d1->charge());
          charge2.push_back(d2->charge());
          
          double pxd3 = -999.9;
          double pyd3 = -999.9;
          double pzd3 = -999.9;
          if(threeProngDecay_ && d3)
          {
            pxd3 = d3->px();
            pyd3 = d3->py();
            pzd3 = d3->pz();
            pt3.push_back(d3->pt());
            p3.push_back(d3->p());
            eta3.push_back(d3->eta());
            phi3.push_back(d3->phi());
            charge3.push_back(d3->charge());
          }
          TVector3 dauvec3(pxd3,pyd3,pzd3);

          pid1.push_back(-99999);
          pid2.push_back(-99999);
          if(doGenMatchingTOF_)
          {
            for(unsigned it=0; it<genpars->size(); ++it){

                const reco::GenParticle & trk = (*genpars)[it];

                if(trk.pt()<0.001) continue;

                int id = trk.pdgId();
                TVector3 trkvect(trk.px(),trk.py(),trk.pz());

                if(fabs(id)!=PID_ && trk.charge())
                {
                  // matching daughter 1
                  double deltaR = trkvect.DeltaR(dauvec1);
                  if(deltaR < deltaR_ && fabs((trk.pt()-pt1.at(it))/pt1.at(it)) < 0.5 && trk.charge()==charge1.at(it) && pid1.push_back(=-99999)
                  {
                    pid1.push_back(id));
                  } 

                  // matching daughter 2
                  deltaR = trkvect.DeltaR(dauvec2);
                  if(deltaR < deltaR_ && fabs((trk.pt()-pt2.at(it))/pt2.at(it)) < 0.5 && trk.charge()==charge2.at(it) && pid2.push_back(=-99999)
                  {
                    pid2.push_back(id));
                  }
                }

                if(fabs(id)==PID_ && trk.numberOfDaughters()==2)
                {
                  const reco::Candidate * Dd1 = trk.daughter(0);
                  const reco::Candidate * Dd2 = trk.daughter(1);
                  TVector3 d1vect(Dd1->px(),Dd1->py(),Dd1->pz());
                  TVector3 d2vect(Dd2->px(),Dd2->py(),Dd2->pz());
                  int id1 = Dd1->pdgId();
                  int id2 = Dd2->pdgId();
                
                  double deltaR = d1vect.DeltaR(dauvec1);
                  if(deltaR < deltaR_ && fabs((Dd1->pt()-pt1.at(it))/pt1.at(it)) < 0.5 && Dd1->charge()==charge1.at(it) && pid1.push_back(=-99999)
                  {
                    pid1.push_back(id1));
                  }
                  deltaR = d2vect.DeltaR(dauvec1);
                  if(deltaR < deltaR_ && fabs((Dd2->pt()-pt1.at(it))/pt1.at(it)) < 0.5 && Dd2->charge()==charge1.at(it) && pid1.push_back(=-99999)
                  {
                    pid1.push_back(id1));
                  }

                  deltaR = d1vect.DeltaR(dauvec2);
                  if(deltaR < deltaR_ && fabs((Dd1->pt()-pt2.at(it))/pt2.at(it)) < 0.5 && Dd1->charge()==charge2.at(it) && pid2.push_back(=-99999)
                  {
                    pid2.push_back(id2));
                  }
                  deltaR = d2vect.DeltaR(dauvec2);
                  if(deltaR < deltaR_ && fabs((Dd2->pt()-pt2.at(it))/pt2.at(it)) < 0.5 && Dd2->charge()==charge2.at(it) && pid2.push_back(=-99999)
                  {
                    pid2.push_back(id2));
                  }
                }

                if(pid1.at(it)!=-99999 && pid2.at(it)!=-99999) break;
            }
          }

          //vtxChi2
          vtxChi2.push_back(trk.vertexChi2());
          ndf.push_back(trk.vertexNdof());
          VtxProb.push_back(TMath::Prob(vtxChi2.at(it),ndf.at(it)));
          
          //PAngle
          TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
          TVector3 secvec(px,py,pz);
          
          TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
          TVector3 secvec2D(px,py,0);
          
          agl.push_back(cos(secvec.Angle(ptosvec)));
          agl_abs.push_back(secvec.Angle(ptosvec));
          
          agl2D.push_back(cos(secvec2D.Angle(ptosvec2D)));
          agl2D_abs.push_back(secvec2D.Angle(ptosvec2D));
          
          //Decay length 3D
          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;
          typedef ROOT::Math::SVector<double, 6> SVector6;
          
          SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
          SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
          
          dl.push_back(ROOT::Math::Mag(distanceVector));
          dlerror.push_back(sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl.at(it));
          
          dlos.push_back(dl.at(it)/dlerror.at(it));
          
          //Decay length 2D
          SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
          SVector6 v2(trk.vertexCovariance(0,0), trk.vertexCovariance(0,1),trk.vertexCovariance(1,1),0,0,0);
          
          SMatrixSym3D sv1(v1);
          SMatrixSym3D sv2(v2);
          
          SMatrixSym3D totalCov2D = sv1 + sv2;
          SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
          
          dl2D.push_back(ROOT::Math::Mag(distanceVector2D));
          double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D.at(it);
          
          dlos2D.push_back(dl2D.at(it)/dl2Derror);

          //trk info
          auto dau1 = d1->get<reco::TrackRef>();
          if(!twoLayerDecay_)
          {
              //trk quality
              trkquality1.push_back(dau1->quality(reco::TrackBase::highPurity));
              
              //trk dEdx
              H2dedx1.push_back(-999.9);
              
              if(dEdxHandle1.isValid()){
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                  H2dedx1.push_back(dEdxTrack[dau1].dEdx());
              }
              
              T4dedx1.push_back(-999.9);
              
              if(dEdxHandle2.isValid()){
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                  T4dedx1.push_back(dEdxTrack[dau1].dEdx());
              }
              
              //track Chi2
              trkChi1.push_back(dau1->normalizedChi2());
              
              //track pT error
              ptErr1.push_back(dau1->ptError());
              
              //vertexCovariance 00-xError 11-y 22-z
              secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
              
              //trkNHits
              nhit1.push_back(dau1->numberOfValidHits());
              
              //DCA
              math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
              
              double dzbest1 = dau1->dz(bestvtx);
              double dxybest1 = dau1->dxy(bestvtx);
              double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
              double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
              
              dzos1.push_back(dzbest1/dzerror1);
              dxyos1.push_back(dxybest1/dxyerror1);
          }
          
          auto dau2 = d2->get<reco::TrackRef>();
          
          //trk quality
          trkquality2.push_back(dau2->quality(reco::TrackBase::highPurity));
          
          //trk dEdx
          H2dedx2.push_back(-999.9);
          
          if(dEdxHandle1.isValid()){
              const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
              H2dedx2.push_back(dEdxTrack[dau2].dEdx());
          }
          
          T4dedx2.push_back(-999.9);
          
          if(dEdxHandle2.isValid()){
              const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
              T4dedx2.push_back(dEdxTrack[dau2].dEdx());
          }
          
          //track Chi2
          trkChi2.push_back(dau2->normalizedChi2());
          
          //track pT error
          ptErr2.push_back(dau2->ptError());
          
          //vertexCovariance 00-xError 11-y 22-z
          secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
          
          //trkNHits
          nhit2.push_back(dau2->numberOfValidHits());
          
          //DCA
          math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
          
          double dzbest2 = dau2->dz(bestvtx);
          double dxybest2 = dau2->dxy(bestvtx);
          double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
          double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
          
          dzos2.push_back(dzbest2/dzerror2);
          dxyos2.push_back(dxybest2/dxyerror2);
          
          if(doMuon_)
          {
            edm::Handle<reco::MuonCollection> theMuonHandle;
            iEvent.getByToken(tok_muon_, theMuonHandle);
              
            // nmatchedch1.push_back(-1);
            // nmatchedst1.push_back(-1);
            // matchedenergy1.push_back(-1);
            // nmatchedch2.push_back(-1);
            // nmatchedst2.push_back(-1);
            // matchedenergy2.push_back(-1);
              
            // double x_exp = -999.;
            // double y_exp = -999.;
            // double xerr_exp = -999.;
            // double yerr_exp = -999.;
            // double dxdz_exp = -999.;
            // double dydz_exp = -999.;
            // double dxdzerr_exp = -999.;
            // double dydzerr_exp = -999.;
              
            // double x_seg = -999.;
            // double y_seg = -999.;
            // double xerr_seg = -999.;
            // double yerr_seg = -999.;
            // double dxdz_seg = -999.;
            // double dydz_seg = -999.;
            // double dxdzerr_seg = -999.;
            // double dydzerr_seg = -999.;
              
            // double dx_seg = 999.;
            // double dy_seg = 999.;
            // double dxerr_seg = 999.;
            // double dyerr_seg = 999.;
            // double dxSig_seg = 999.;
            // double dySig_seg = 999.;
            // double ddxdz_seg = 999.;
            // double ddydz_seg = 999.;
            // double ddxdzerr_seg = 999.;
            // double ddydzerr_seg = 999.;
            // double ddxdzSig_seg = 999.;
            // double ddydzSig_seg = 999.;
              
            // onestmuon1.push_back(false);
            // pfmuon1.push_back(false);
            // glbmuon1.push_back(false);
            // trkmuon1.push_back(false);
            // calomuon1.push_back(false); 
            // softmuon1.push_back(false);
            // onestmuon2.push_back(false);
            // pfmuon2.push_back(false);
            // glbmuon2.push_back(false);
            // trkmuon2.push_back(false);
            // calomuon2.push_back(false);
            // softmuon2.push_back(false);

            const int muId1 = muAssocToTrack( dau1, theMuonHandle );
            const int muId2 = muAssocToTrack( dau2, theMuonHandle );

            if( muId1 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId1];

              onestmuon1.push_back(muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight")));
              pfmuon1.push_back( cand.isPFMuon());
              glbmuon1.push_back( cand.isGlobalMuon());
              trkmuon1.push_back( cand.isTrackerMuon());
              calomuon1.push_back( cand.isCaloMuon());

              if( 
                  //glbmuon1.at(it) && 
                  trkmuon1.at(it) &&
                  cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && 
                  cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 && 
                  fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
                  fabs(cand.innerTrack()->dz(vtx.position())) < 20.
                ) softmuon1.push_back(true);
            }

            if( muId2 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId2];

              onestmuon2.push_back(muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight")));
              pfmuon2.push_back( cand.isPFMuon());
              glbmuon2.push_back( cand.isGlobalMuon());
              trkmuon2.push_back( cand.isTrackerMuon());
              calomuon2.push_back( cand.isCaloMuon());

              if(
                  //glbmuon2.at(it) && 
                  trkmuon2.at(it) &&
                  cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                  cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
                  fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
                  fabs(cand.innerTrack()->dz(vtx.position())) < 20.
                ) softmuon2.push_back(true);
            }

            if(doMuonFull_)
            {

            if( muId1 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId1];

              nmatchedch1.push_back(cand.numberOfMatches());
              nmatchedst1.push_back(cand.numberOfMatchedStations());
                    
              reco::MuonEnergy muenergy = cand.calEnergy();
              matchedenergy1.push_back(muenergy.hadMax);
                      
              const std::vector<reco::MuonChamberMatch>& muchmatches = cand.matches();
                    
              for(unsigned int ich=0;ich<muchmatches.size();ich++)
              {
                x_exp = muchmatches[ich].x;
                y_exp = muchmatches[ich].y;
                xerr_exp = muchmatches[ich].xErr;
                yerr_exp = muchmatches[ich].yErr;
                dxdz_exp = muchmatches[ich].dXdZ;
                dydz_exp = muchmatches[ich].dYdZ;
                dxdzerr_exp = muchmatches[ich].dXdZErr;
                dydzerr_exp = muchmatches[ich].dYdZErr;
                          
                std::vector<reco::MuonSegmentMatch> musegmatches = muchmatches[ich].segmentMatches;
                          
                if(!musegmatches.size()) continue;
                for(unsigned int jseg=0;jseg<musegmatches.size();jseg++)
                {
                  x_seg = musegmatches[jseg].x;
                  y_seg = musegmatches[jseg].y;
                  xerr_seg = musegmatches[jseg].xErr;
                  yerr_seg = musegmatches[jseg].yErr;
                  dxdz_seg = musegmatches[jseg].dXdZ;
                  dydz_seg = musegmatches[jseg].dYdZ;
                  dxdzerr_seg = musegmatches[jseg].dXdZErr;
                  dydzerr_seg = musegmatches[jseg].dYdZErr;
                              
                  if(sqrt((x_seg-x_exp)*(x_seg-x_exp)+(y_seg-y_exp)*(y_seg-y_exp))<sqrt(dx_seg*dx_seg+dy_seg*dy_seg))
                  {
                    dx_seg = x_seg - x_exp;
                    dy_seg = y_seg - y_exp;
                    dxerr_seg = sqrt(xerr_seg*xerr_seg+xerr_exp*xerr_exp);
                    dyerr_seg = sqrt(yerr_seg*yerr_seg+yerr_exp*yerr_exp);
                    dxSig_seg = dx_seg / dxerr_seg;
                    dySig_seg = dy_seg / dyerr_seg;
                    ddxdz_seg = dxdz_seg - dxdz_exp;
                    ddydz_seg = dydz_seg - dydz_exp;
                    ddxdzerr_seg = sqrt(dxdzerr_seg*dxdzerr_seg+dxdzerr_exp*dxdzerr_exp);
                    ddydzerr_seg = sqrt(dydzerr_seg*dydzerr_seg+dydzerr_exp*dydzerr_exp);
                    ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                    ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                  }
                }
                        
                dx1_seg_.push_back(dx_seg);
                dy1_seg_.push_back(dy_seg);
                dxSig1_seg_.push_back(dxSig_seg);
                dySig1_seg_.push_back(dySig_seg);
                ddxdz1_seg_.push_back(ddxdz_seg);
                ddydz1_seg_.push_back(ddydz_seg);
                ddxdzSig1_seg_.push_back(ddxdzSig_seg);
                ddydzSig1_seg_.push_back(ddydzSig_seg);
              }
            } 

            if( muId2 != -1 )
            {
              const reco::Muon& cand = (*theMuonHandle)[muId2];

              nmatchedch2.push_back(cand.numberOfMatches());
              nmatchedst2.push_back(cand.numberOfMatchedStations());
                      
              reco::MuonEnergy muenergy = cand.calEnergy();
              matchedenergy2.push_back(muenergy.hadMax);
                      
              const std::vector<reco::MuonChamberMatch>& muchmatches = cand.matches();
              for(unsigned int ich=0;ich<muchmatches.size();ich++)
                          //                        for(unsigned int ich=0;ich<1;ich++)
              {
                x_exp = muchmatches[ich].x;
                y_exp = muchmatches[ich].y;
                xerr_exp = muchmatches[ich].xErr;
                yerr_exp = muchmatches[ich].yErr;
                dxdz_exp = muchmatches[ich].dXdZ;
                dydz_exp = muchmatches[ich].dYdZ;
                dxdzerr_exp = muchmatches[ich].dXdZErr;
                dydzerr_exp = muchmatches[ich].dYdZErr;
                          
                std::vector<reco::MuonSegmentMatch> musegmatches = muchmatches[ich].segmentMatches;
                          
                if(!musegmatches.size()) continue;
                for(unsigned int jseg=0;jseg<musegmatches.size();jseg++)
                {
                  x_seg = musegmatches[jseg].x;
                  y_seg = musegmatches[jseg].y;
                  xerr_seg = musegmatches[jseg].xErr;
                  yerr_seg = musegmatches[jseg].yErr;
                  dxdz_seg = musegmatches[jseg].dXdZ;
                  dydz_seg = musegmatches[jseg].dYdZ;
                  dxdzerr_seg = musegmatches[jseg].dXdZErr;
                  dydzerr_seg = musegmatches[jseg].dYdZErr;
                              
                  if(sqrt((x_seg-x_exp)*(x_seg-x_exp)+(y_seg-y_exp)*(y_seg-y_exp))<sqrt(dx_seg*dx_seg+dy_seg*dy_seg))
                  {
                    dx_seg = x_seg - x_exp;
                    dy_seg = y_seg - y_exp;
                    dxerr_seg = sqrt(xerr_seg*xerr_seg+xerr_exp*xerr_exp);
                    dyerr_seg = sqrt(yerr_seg*yerr_seg+yerr_exp*yerr_exp);
                    dxSig_seg = dx_seg / dxerr_seg;
                    dySig_seg = dy_seg / dyerr_seg;
                    ddxdz_seg = dxdz_seg - dxdz_exp;
                    ddydz_seg = dydz_seg - dydz_exp;
                    ddxdzerr_seg = sqrt(dxdzerr_seg*dxdzerr_seg+dxdzerr_exp*dxdzerr_exp);
                    ddydzerr_seg = sqrt(dydzerr_seg*dydzerr_seg+dydzerr_exp*dydzerr_exp);
                    ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                    ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                  }
                }
                          
                dx2_seg_.push_back(dx_seg);
                dy2_seg_.push_back(dy_seg);
                dxSig2_seg_.push_back(dxSig_seg);
                dySig2_seg_.push_back(dySig_seg);
                ddxdz2_seg_.push_back(ddxdz_seg);
                ddydz2_seg_.push_back(ddydz_seg);
                ddxdzSig2_seg_.push_back(ddxdzSig_seg);
                ddydzSig2_seg_.push_back(ddydzSig_seg);
              }
            }
            } // doMuonFull
          }
          
          if(twoLayerDecay_)
          {
              grand_mass.push_back(d1->mass());
              
              // const reco::Candidate * gd1 = d1->daughter(0);
              // const reco::Candidate * gd2 = d1->daughter(1);
              
              double gpxd1 = gd1->px();
              double gpyd1 = gd1->py();
              double gpzd1 = gd1->pz();
              double gpxd2 = gd2->px();
              double gpyd2 = gd2->py();
              double gpzd2 = gd2->pz();
              
              TVector3 gdauvec1(gpxd1,gpyd1,gpzd1);
              TVector3 gdauvec2(gpxd2,gpyd2,gpzd2);
              
              auto gdau1 = gd1->get<reco::TrackRef>();
              auto gdau2 = gd2->get<reco::TrackRef>();
              
              //trk quality
              
              grand_trkquality1.push_back(gdau1->quality(reco::TrackBase::highPurity));
              grand_trkquality2.push_back(gdau2->quality(reco::TrackBase::highPurity));
              
              //trk dEdx
              grand_H2dedx1.push_back(-999.9);
              grand_H2dedx2.push_back(-999.9);
              
              if(dEdxHandle1.isValid()){
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                  grand_H2dedx1.push_back(dEdxTrack[gdau1].dEdx());
                  grand_H2dedx2.push_back(dEdxTrack[gdau2].dEdx());
              }
              
              grand_T4dedx1.push_back(-999.9);
              grand_T4dedx2.push_back(-999.9);
              
              if(dEdxHandle2.isValid()){
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                  grand_T4dedx1.push_back(dEdxTrack[gdau1].dEdx());
                  grand_T4dedx2.push_back(dEdxTrack[gdau2].dEdx());
              }
              
              //track pt
              grand_pt1.push_back(gd1->pt());
              grand_pt2.push_back(gd2->pt());
              
              //track momentum
              grand_p1.push_back(gd1->p());
              grand_p2.push_back(gd2->p());
              
              //track eta
              grand_eta1.push_back(gd1->eta());
              grand_eta2.push_back(gd2->eta());
              
              //track charge
              grand_charge1.push_back(gd1->charge());
              grand_charge2.push_back(gd2->charge());
              
              //track Chi2
              grand_trkChi1.push_back(gdau1->normalizedChi2());
              grand_trkChi2.push_back(gdau2->normalizedChi2());
              
              //track pT error
              grand_ptErr1.push_back(gdau1->ptError());
              grand_ptErr2.push_back(gdau2->ptError());
              
              //vertexCovariance 00-xError 11-y 22-z
              secvz = d1->vz(); secvx = d1->vx(); secvy = d1->vy();
              
              //trkNHits
              grand_nhit1.push_back(gdau1->numberOfValidHits());
              grand_nhit2.push_back(gdau2->numberOfValidHits());
              
              //DCA
              math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
              
              double gdzbest1 = gdau1->dz(bestvtx);
              double gdxybest1 = gdau1->dxy(bestvtx);
              double gdzerror1 = sqrt(gdau1->dzError()*gdau1->dzError()+bestvzError*bestvzError);
              double gdxyerror1 = sqrt(gdau1->d0Error()*gdau1->d0Error()+bestvxError*bestvyError);
              
              grand_dzos1.push_back(gdzbest1/gdzerror1);
              grand_dxyos1.push_back(gdxybest1/gdxyerror1);
              
              double gdzbest2 = gdau2->dz(bestvtx);
              double gdxybest2 = gdau2->dxy(bestvtx);
              double gdzerror2 = sqrt(gdau2->dzError()*gdau2->dzError()+bestvzError*bestvzError);
              double gdxyerror2 = sqrt(gdau2->d0Error()*gdau2->d0Error()+bestvxError*bestvyError);
              
              grand_dzos2.push_back(gdzbest2/gdzerror2);
              grand_dxyos2.push_back(gdxybest2/gdxyerror2);
              
              //vtxChi2
              grand_vtxChi2.push_back(d1->vertexChi2());
              grand_ndf.push_back(d1->vertexNdof());
              grand_VtxProb.push_back(TMath::Prob(grand_vtxChi2.at(it),grand_ndf.at(it)));
              
              //PAngle
              TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
              TVector3 secvec(d1->px(),d1->py(),d1->pz());
              
              TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
              TVector3 secvec2D(d1->px(),d1->py(),0);
              
              grand_agl.push_back(cos(secvec.Angle(ptosvec)));
              grand_agl_abs.push_back(secvec.Angle(ptosvec));
              
              grand_agl2D.push_back(cos(secvec2D.Angle(ptosvec2D)));
              grand_agl2D_abs.push_back(secvec2D.Angle(ptosvec2D));
              
              //Decay length 3D
              typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
              typedef ROOT::Math::SVector<double, 3> SVector3;
              typedef ROOT::Math::SVector<double, 6> SVector6;
              
              SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
              SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
              
              grand_dl.push_back(ROOT::Math::Mag(distanceVector));
              grand_dlerror.push_back(sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl.at(it));
              
              grand_dlos.push_back(grand_dl.at(it)/grand_dlerror.at(it));
              
              //Decay length 2D
              SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
              SVector6 v2(d1->vertexCovariance(0,0), d1->vertexCovariance(0,1),d1->vertexCovariance(1,1),0,0,0);
              
              SMatrixSym3D sv1(v1);
              SMatrixSym3D sv2(v2);
              
              SMatrixSym3D totalCov2D = sv1 + sv2;
              SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
              
              double gdl2D = ROOT::Math::Mag(distanceVector2D);
              double gdl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
              
              grand_dlos2D.push_back(gdl2D/gdl2Derror);
          }
  #ifdef DEBUG
  cout << "Done reco single iter" << endl;
  #endif

          if(saveHistogram_)
          {
            for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
              for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
              {
                if(pt.at(it)<pTBins_[ipt+1] && pt.at(it)>pTBins_[ipt] && y.at(it)<yBins_[iy+1] && y.at(it)>yBins_[iy])
                {
                  hMassVsMVA[iy][ipt]->Fill(mva.at(it),mass.at(it));
  //                h3DDCAVsMVA[iy][ipt]->Fill(mva.at(it),dl.at(it)*sin(agl_abs.at(it)));
  //                h2DDCAVsMVA[iy][ipt]->Fill(mva.at(it),dl2D.at(it)*sin(agl2D_abs.at(it)));

                  if(saveAllHistogram_)
                  {
                  hpTVsMVA[iy][ipt]->Fill(mva.at(it),pt.at(it));
                  hetaVsMVA[iy][ipt]->Fill(mva.at(it),eta.at(it));
                  hyVsMVA[iy][ipt]->Fill(mva.at(it),y.at(it));
                  hVtxProbVsMVA[iy][ipt]->Fill(mva.at(it),VtxProb.at(it));
                  h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva.at(it),agl.at(it));
                  h3DPointingAngleVsMVA[iy][ipt]->Fill(mva.at(it),agl_abs.at(it));
                  h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva.at(it),agl2D.at(it));
                  h2DPointingAngleVsMVA[iy][ipt]->Fill(mva.at(it),agl2D_abs.at(it));
                  h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva.at(it),dlos.at(it));
                  h3DDecayLengthVsMVA[iy][ipt]->Fill(mva.at(it),dl.at(it));
                  h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva.at(it),dlos2D.at(it));
                  h2DDecayLengthVsMVA[iy][ipt]->Fill(mva.at(it),dl2D.at(it));
                  hzDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva.at(it),dzos1.at(it));
                  hxyDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva.at(it),dxyos1.at(it));
                  hNHitD1VsMVA[iy][ipt]->Fill(mva.at(it),nhit1.at(it));
                  hpTD1VsMVA[iy][ipt]->Fill(mva.at(it),pt1.at(it));
                  hpTerrD1VsMVA[iy][ipt]->Fill(mva.at(it),ptErr1.at(it)/pt1.at(it));
                  hEtaD1VsMVA[iy][ipt]->Fill(mva.at(it),eta1.at(it));
                  hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva.at(it),H2dedx1.at(it));
                  hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1.at(it),H2dedx1.at(it));
                  hzDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva.at(it),dzos2.at(it));
                  hxyDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva.at(it),dxyos2.at(it));
                  hNHitD2VsMVA[iy][ipt]->Fill(mva.at(it),nhit2.at(it));
                  hpTD2VsMVA[iy][ipt]->Fill(mva.at(it),pt2.at(it));
                  hpTerrD2VsMVA[iy][ipt]->Fill(mva.at(it),ptErr2.at(it)/pt2.at(it));
                  hEtaD2VsMVA[iy][ipt]->Fill(mva.at(it),eta2.at(it));
                  hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva.at(it),H2dedx2.at(it));
                  hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2.at(it),H2dedx2.at(it));
                  if(threeProngDecay_)
                  {
                    hzDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva.at(it),dzos3.at(it));
                    hxyDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva.at(it),dxyos3.at(it));
                    hNHitD3VsMVA[iy][ipt]->Fill(mva.at(it),nhit3.at(it));
                    hpTD3VsMVA[iy][ipt]->Fill(mva.at(it),pt3.at(it));
                    hpTerrD3VsMVA[iy][ipt]->Fill(mva.at(it),ptErr3.at(it)/pt3.at(it));
                    hEtaD3VsMVA[iy][ipt]->Fill(mva.at(it),eta3.at(it));
                    hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva.at(it),H2dedx3.at(it));
                    hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1.at(it),H2dedx3.at(it));
                  }

                  }
                }
              }
          }

      }
  #ifdef DEBUG
  cout << "Fill reco done" << endl;
  #endif
  }

  void
  VertexCompositeTreeProducer2::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
      #ifdef DEBUG
      cout << "Fill GEN Start" << endl;
      #endif
      edm::Handle<reco::GenParticleCollection> genpars;
      iEvent.getByToken(tok_genParticle_,genpars);
      std::vector<reco::GenParticleRef> genRefs;
      for(unsigned it=0; it<genpars->size(); ++it){

          const reco::GenParticle & trk = (*genpars)[it];

          int id = trk.pdgId();

          if(fabs(id)!=PID_) continue; //check is target
          if(decayInGen_ && trk.numberOfDaughters()!=2 && !threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
          if(decayInGen_ && trk.numberOfDaughters()!=3 && threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
          int nDau = threeProngDecay_ ? 3 : 2;
          std::vector<unsigned int> idxs;
          std::vector<unsigned int> permutations(nDau);
          std::iota(permutations.begin(), permutations.end(), 0);
          std::sort(permutations.begin(), permutations.end());
          if (!threeProngDecay_) {
            do {
              auto Dd1 = trk.daughter( permutations.at(0) );
              auto Dd2 = trk.daughter( permutations.at(1) );
              if (abs(Dd1->pdgId()) == PID_dau1_ && abs(Dd2->pdgId()) == PID_dau2_) {
                if(twoLayerDecay_){
                  // Magic numbers, _permutations -> number of D0 daughters;
                  std::vector<unsigned int> _permutations(2);
                  std::iota(_permutations.begin(), _permutations.end(), 0);
                  std::sort(_permutations.begin(), _permutations.end());
                  do {
                    auto Ddd1 = Dd1->daughter( _permutations.at(0) );
                    auto Ddd2 = Dd1->daughter( _permutations.at(1) );
                    if (abs(Ddd1->pdgId()) == 211 && abs(Ddd2->pdgId()) == 321) {
                      idxs = permutations;
                      break;
                    }
                  } while (std::next_permutation(_permutations.begin(), _permutations.end()));
                  if(!idxs.empty()) break;
                } else {
                  if (abs(Dd1->pdgId()) == PID_dau1_
                      && abs(Dd2->pdgId()) == PID_dau2_
                      ) {
                    idxs = permutations;
                    break;
                  } 
                  if(!idxs.empty()) break;
                }
              }
            } while (std::next_permutation(permutations.begin(), permutations.end()));
          } else {
            do {
              auto Dd1 = trk.daughter( permutations.at(0) );
              auto Dd2 = trk.daughter( permutations.at(1) );
              auto Dd3 = trk.daughter( permutations.at(2) );

              if (abs(Dd1->pdgId()) == PID_dau1_
                  && abs(Dd2->pdgId()) == PID_dau2_
                  && abs(Dd3->pdgId() == PID_dau3_)) {
                idxs = permutations;
                break;
              }
            } while (std::next_permutation(permutations.begin(), permutations.end()));
          }
              if (decayInGen_ && idxs.empty()) continue;
              genRefs.push_back(reco::GenParticleRef(genpars, it));
      }
      if(twoLayerDecay_){
        unsigned int nGen = genRefs.size();
        candSize_gen = nGen;
        for( unsigned int igen=0; igen<nGen; igen++){
          #ifdef DEBUG
          cout <<"nGen : "<< nGen << endl;
          #endif
          auto const theGenDStar = genRefs.at(igen);
          if(abs(theGenDStar->pdgId())!=413) cout << "id : " << theGenDStar->pdgId() << endl;
          unsigned int idxD0 = 1;
          cout << idxD0 << endl;
          cout <<  genRefs.at(igen)->daughter(idxD0)->pdgId() << endl;
          cout <<  genRefs.at(igen)->daughter(1-idxD0)->pdgId() << endl;
          if( fabs(theGenDStar->daughter(0)->pdgId()) == 421 ) idxD0 = 0;
          auto const* theGenD0 = genRefs.at(igen)->daughter(idxD0);
          auto const* theGenPion = genRefs.at(igen)->daughter(1- idxD0);
          mass_gen.at(igen) = theGenDStar->mass();
          pt_gen.at(igen) = theGenDStar->pt();
          eta_gen.at(igen) = theGenDStar->eta(); 
          phi_gen.at(igen) = theGenDStar->phi();
          y_gen.at(igen) = theGenDStar->rapidity();
          status_gen.at(igen) = theGenDStar->status();
          idmom.at(igen) = -77;
          if(theGenDStar->numberOfMothers()!=0){
            const reco::Candidate * mom = theGenDStar->mother();
            idmom.at(igen) = mom->pdgId();
          }


          gen_D0mass_.at(igen) = theGenD0->mass();
          gen_D0pT_.at(igen) = theGenD0->pt();
          gen_D0eta_.at(igen) = theGenD0->eta();
          gen_D0phi_.at(igen) = theGenD0->phi();
          gen_D0y_.at(igen) = theGenD0->rapidity();
          gen_D0pdgId_.at(igen) = theGenD0->pdgId();
          gen_D1mass_.at(igen) = theGenPion->mass(); 
          gen_D1pT_.at(igen) = theGenPion->pt();
          gen_D1eta_.at(igen) = theGenPion->eta();
          gen_D1phi_.at(igen) = theGenPion->phi();
          gen_D1y_.at(igen) = theGenPion->rapidity();
          gen_D1pdgId_.at(igen) = theGenPion->pdgId();

          const auto* genDau0 = theGenD0->daughter(0);
          const auto* genDau1 = theGenD0->daughter(1);
          gen_D0Dau1_pT_.at(igen) = genDau0->pt();
          gen_D0Dau1_eta_.at(igen) = genDau0->eta();
          gen_D0Dau1_phi_.at(igen) = genDau0->phi();
          gen_D0Dau1_y_.at(igen) = genDau0->rapidity();
          gen_D0Dau1_pdgId_.at(igen) = genDau0->pdgId();
          #ifdef DEBUG
          cout << "D0 dau1 pdgId : " << genDau0->pdgId() << endl;
          #endif
          // cout << "D0 dau1 pdgId : " <<  gen_D0Dau1_pdgId_.at(igen) << endl;
          gen_D0Dau2_pT_.at(igen) = genDau1->pt();

          gen_D0Dau2_eta_.at(igen) = genDau1->eta();
          gen_D0Dau2_phi_.at(igen) = genDau1->phi();
          gen_D0Dau2_y_.at(igen) = genDau1->rapidity();
          gen_D0Dau2_pdgId_.at(igen) = genDau1->pdgId();
    }
    
  }

          //if (genRefs.size()>1) std::cout << "More than one target of generated particles\n";
      
      //     if(trk.numberOfMothers()!=0)
      //     {
      //         const reco::Candidate * mom = trk.mother();
      //         idmom[candSize_gen-1] = mom->pdgId();
      //     }
              
      //         const reco::Candidate * Dd1 = trk.daughter(0);
      //         const reco::Candidate * Dd2 = trk.daughter(1);
          



      //     pt_gen[candSize_gen-1] = trk.pt();
      //     eta_gen[candSize_gen-1] = trk.eta();
      //     status_gen[candSize_gen-1] = trk.status();
      //     idmom[candSize_gen-1] = -77;
      //     y_gen[candSize_gen-1] = trk.rapidity();

      //     if(trk.numberOfMothers()!=0)
      //     {
      //         const reco::Candidate * mom = trk.mother();
      //         idmom[candSize_gen-1] = mom->pdgId();
      //     }

      //     if(!decayInGen_) continue;

      //     const reco::Candidate * Dd1 = trk.daughter(0);
      //     const reco::Candidate * Dd2 = trk.daughter(1);
      //     const reco::Candidate * Dd3 = trk.daughter(2);

      //     iddau1[candSize_gen-1] = fabs(Dd1->pdgId());
      //     iddau2[candSize_gen-1] = fabs(Dd2->pdgId());
      //     if(Dd3) iddau3[candSize_gen-1] = fabs(Dd3->pdgId());
      // }
  }

  // ------------ method called once each job just before starting event
  //loop  ------------
  void
  VertexCompositeTreeProducer2::beginJob()
  {
      TH1D::SetDefaultSumw2();
      
      if(!doRecoNtuple_ && !doGenNtuple_)
      {
          cout<<"No output for either RECO or GEN!! Fix config!!"<<endl; return;
      }

      if(twoLayerDecay_ && doMuon_)
      {
          cout<<"Muons cannot be coming from two layer decay!! Fix config!!"<<endl; return;
      }
      
      if(saveHistogram_) initHistogram();
      if(saveTree_) initTree();
  }

  void
  VertexCompositeTreeProducer2::initHistogram()
  {
    for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
    {
      for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
    {
    hMassVsMVA[iy][ipt] = fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d",iy,ipt),";mva;mass(GeV)",100,-1.,1.,massHistBins_,massHistPeak_-massHistWidth_,massHistPeak_+massHistWidth_);
  //   h3DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;3D DCA;",100,-1.,1.,1000,0,10);
  //   h2DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;2D DCA;",100,-1.,1.,1000,0,10);

    if(saveAllHistogram_)
    {
    hpTVsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTVsMVA_y%d_pt%d",iy,ipt),";mva;pT;",100,-1,1,100,0,10);
    hetaVsMVA[iy][ipt] = fs->make<TH2F>(Form("hetaVsMVA_y%d_pt%d",iy,ipt),";mva;eta;",100,-1.,1.,40,-4,4);
    hyVsMVA[iy][ipt] = fs->make<TH2F>(Form("hyVsMVA_y%d_pt%d",iy,ipt),";mva;y;",100,-1.,1.,40,-4,4);
    hVtxProbVsMVA[iy][ipt] = fs->make<TH2F>(Form("hVtxProbVsMVA_y%d_pt%d",iy,ipt),";mva;VtxProb;",100,-1.,1.,100,0,1);
    h3DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DCosPointingAngle;",100,-1.,1.,100,-1,1);
    h3DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
    h2DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DCosPointingAngle;",100,-1.,1.,100,-1,1);
    h2DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
    h3DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLengthSignificance;",100,-1.,1.,300,0,30);
    h2DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLengthSignificance;",100,-1.,1.,300,0,30);
    h3DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLength;",100,-1.,1.,300,0,30);
    h2DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLength;",100,-1.,1.,300,0,30);
    hzDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther1VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther1;",100,-1.,1.,100,-10,10);
    hxyDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther1VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther1;",100,-1.,1.,100,-10,10);
    hNHitD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD1VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD1;",100,-1.,1.,100,0,100);
    hpTD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTD1;",100,-1.,1.,100,0,10);
    hpTerrD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD1;",100,-1.,1.,50,0,0.5);
    hEtaD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD1VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD1;",100,-1.,1.,40,-4,4);
    hdedxHarmonic2D1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D1;",100,-1.,1.,100,0,10);
    hdedxHarmonic2D1VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D1",100,0,10,100,0,10);
    hzDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther2VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther2;",100,-1.,1.,100,-10,10);
    hxyDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther2VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther2;",100,-1.,1.,100,-10,10);
    hNHitD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD2VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD2;",100,-1.,1.,100,0,100);
    hpTD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTD2;",100,-1.,1.,100,0,10);
    hpTerrD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD2;",100,-1.,1.,50,0,0.5);
    hEtaD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD2VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD2;",100,-1.,1.,40,-4,4);
    hdedxHarmonic2D2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D2;",100,-1.,1.,100,0,10);
    hdedxHarmonic2D2VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D2",100,0,10,100,0,10);

    if(threeProngDecay_)
    {
      hzDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther3VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificanceDaugther3;",100,-1.,1.,100,-10,10);
      hxyDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther3VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificanceDaugther3;",100,-1.,1.,100,-10,10);
      hNHitD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD3VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD3;",100,-1.,1.,100,0,100);
      hpTD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD3VsMVA_y%d_pt%d",iy,ipt),";mva;pTD3;",100,-1.,1.,100,0,10);
      hpTerrD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD3VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD3;",100,-1.,1.,50,0,0.5);
      hEtaD3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD3VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD3;",100,-1.,1.,40,-4,4);
      hdedxHarmonic2D3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D3;",100,-1.,1.,100,0,10);
      hdedxHarmonic2D3VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D3",100,0,10,100,0,10);
    }

    }
    }
  }
  }

  void 
  VertexCompositeTreeProducer2::initTree()
  { 
      VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");
      
      if(doRecoNtuple_) 
      { 
    
      // Event info
      VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline);
      VertexCompositeNtuple->Branch("Npixel",&Npixel);
      VertexCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus);
      VertexCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus);
      VertexCompositeNtuple->Branch("ZDCPlus",&ZDCPlus);
      VertexCompositeNtuple->Branch("ZDCMinus",&ZDCMinus);
      VertexCompositeNtuple->Branch("bestvtxX",&bestvx);
      VertexCompositeNtuple->Branch("bestvtxY",&bestvy);
      VertexCompositeNtuple->Branch("bestvtxZ",&bestvz);
      VertexCompositeNtuple->Branch("candSize",&candSize);
      if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality);
      if(isEventPlane_) 
      {
        VertexCompositeNtuple->Branch("ephfpAngle",&ephfpAngle);
        VertexCompositeNtuple->Branch("ephfmAngle",&ephfmAngle);
        VertexCompositeNtuple->Branch("ephfpQ",&ephfpQ);
        VertexCompositeNtuple->Branch("ephfmQ",&ephfmQ);
        VertexCompositeNtuple->Branch("ephfpSumW",&ephfpSumW);
        VertexCompositeNtuple->Branch("ephfmSumW",&ephfmSumW);
      }

      // particle info
      VertexCompositeNtuple->Branch("pT",&pt);
      VertexCompositeNtuple->Branch("y",&y);
      VertexCompositeNtuple->Branch("eta",&eta);
      VertexCompositeNtuple->Branch("phi",&phi);
      VertexCompositeNtuple->Branch("mass",&mass);
      if(useAnyMVA_) VertexCompositeNtuple->Branch("mva",&mva);

      if(!isSkimMVA_)  
      {
          //Composite candidate info RECO
          VertexCompositeNtuple->Branch("flavor",&flavor);
          VertexCompositeNtuple->Branch("VtxProb",&VtxProb);
          VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2);
          VertexCompositeNtuple->Branch("VtxNDF",&ndf);
          VertexCompositeNtuple->Branch("3DCosPointingAngle",&agl);
          VertexCompositeNtuple->Branch("3DPointingAngle",&agl_abs);
          VertexCompositeNtuple->Branch("2DCosPointingAngle",&agl2D);
          VertexCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs);
          VertexCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos);
          VertexCompositeNtuple->Branch("3DDecayLength",&dl);
          VertexCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D);
          VertexCompositeNtuple->Branch("2DDecayLength",&dl2D);
      
          if(doGenMatching_)
          {
              VertexCompositeNtuple->Branch("isSwap",&isSwap);
              VertexCompositeNtuple->Branch("idmom_reco",&idmom_reco);
              VertexCompositeNtuple->Branch("idBAnc_reco",&idBAnc_reco);
              VertexCompositeNtuple->Branch("matchGEN",&matchGEN);
              VertexCompositeNtuple->Branch("matchGen3DPointingAngle",&gen_agl_abs);
              VertexCompositeNtuple->Branch("matchGen2DPointingAngle",&gen_agl2D_abs);
              VertexCompositeNtuple->Branch("matchGen3DDecayLength",&gen_dl);
              VertexCompositeNtuple->Branch("matchGen2DDecayLength",&gen_dl2D);
              // VertexCompositeNtuple->Branch("matchgen_D0pT",&gen_D0pT_);
              // VertexCompositeNtuple->Branch("matchgen_D0eta",&gen_D0eta_);
              // VertexCompositeNtuple->Branch("matchgen_D0phi",&gen_D0phi_);
              // VertexCompositeNtuple->Branch("matchgen_D0mass",&gen_D0mass_);
              // VertexCompositeNtuple->Branch("matchgen_D0y",&gen_D0y_);
              // VertexCompositeNtuple->Branch("matchgen_D0charge",&gen_D0charge_);
              // VertexCompositeNtuple->Branch("matchgen_D0pdgId",&gen_D0pdgId_);
              if(twoLayerDecay_){
                VertexCompositeNtuple->Branch("matchGen_D0pT",&matchGen_D0pT_);
                VertexCompositeNtuple->Branch("matchGen_D0eta",&matchGen_D0eta_);
                VertexCompositeNtuple->Branch("matchGen_D0phi",&matchGen_D0phi_);
                VertexCompositeNtuple->Branch("matchGen_D0mass",&matchGen_D0mass_);
                VertexCompositeNtuple->Branch("matchGen_D0y",&matchGen_D0y_);
                VertexCompositeNtuple->Branch("matchGen_D0charge",&matchGen_D0charge_);
                VertexCompositeNtuple->Branch("matchGen_D0pdgId",&matchGen_D0pdgId_);

                VertexCompositeNtuple->Branch("matchGen_D0Dau1_pT",&matchGen_D0Dau1_pT_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_eta",&matchGen_D0Dau1_eta_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_phi",&matchGen_D0Dau1_phi_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_mass",&matchGen_D0Dau1_mass_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_y",&matchGen_D0Dau1_y_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_charge",&matchGen_D0Dau1_charge_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau1_pdgId",&matchGen_D0Dau1_pdgId_);

                VertexCompositeNtuple->Branch("matchGen_D0Dau2_pT",&matchGen_D0Dau2_pT_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_eta",&matchGen_D0Dau2_eta_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_phi",&matchGen_D0Dau2_phi_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_mass",&matchGen_D0Dau2_mass_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_y",&matchGen_D0Dau2_y_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_charge",&matchGen_D0Dau2_charge_);
                VertexCompositeNtuple->Branch("matchGen_D0Dau2_pdgId",&matchGen_D0Dau2_pdgId_);

                VertexCompositeNtuple->Branch("matchGen_D1pT",&matchGen_D1pT_);
                VertexCompositeNtuple->Branch("matchGen_D1eta",&matchGen_D1eta_);
                VertexCompositeNtuple->Branch("matchGen_D1phi",&matchGen_D1phi_);
                VertexCompositeNtuple->Branch("matchGen_D1mass",&matchGen_D1mass_);
                VertexCompositeNtuple->Branch("matchGen_D1y",&matchGen_D1y_);
                VertexCompositeNtuple->Branch("matchGen_D1charge",&matchGen_D1charge_);
                VertexCompositeNtuple->Branch("matchGen_D1pdgId",&matchGen_D1pdgId_);
                VertexCompositeNtuple->Branch("matchGen_D1decayLength2D_",&matchGen_D1decayLength2D_);
                VertexCompositeNtuple->Branch("matchGen_D1decayLength3D_",&matchGen_D1decayLength3D_);
                VertexCompositeNtuple->Branch("matchGen_D1angle2D_",&matchGen_D1angle2D_);
                VertexCompositeNtuple->Branch("matchGen_D1angle3D_",&matchGen_D1angle3D_);
                VertexCompositeNtuple->Branch("matchGen_D1ancestorId_",&matchGen_D1ancestorId_);
                VertexCompositeNtuple->Branch("matchGen_D1ancestorFlavor_",&matchGen_D1ancestorFlavor_);
              }
          }
          
          if(doGenMatchingTOF_)
          {
            VertexCompositeNtuple->Branch("PIDD1",&pid1);
            VertexCompositeNtuple->Branch("PIDD2",&pid1);
            VertexCompositeNtuple->Branch("TOFD1",&tof1);
            VertexCompositeNtuple->Branch("TOFD2",&tof1);
          }

          //daugther & grand daugther info
          if(twoLayerDecay_)
          {
              VertexCompositeNtuple->Branch("massDaugther1",&grand_mass);
              VertexCompositeNtuple->Branch("pTD1",&pt1);
              VertexCompositeNtuple->Branch("EtaD1",&eta1);
              VertexCompositeNtuple->Branch("PhiD1",&phi1);
              VertexCompositeNtuple->Branch("VtxProbDaugther1",&grand_VtxProb);
              VertexCompositeNtuple->Branch("VtxChi2Daugther1",&grand_vtxChi2);
              VertexCompositeNtuple->Branch("VtxNDFDaugther1",&grand_ndf);
              VertexCompositeNtuple->Branch("3DCosPointingAngleDaugther1",&grand_agl);
              VertexCompositeNtuple->Branch("3DPointingAngleDaugther1",&grand_agl_abs);
              VertexCompositeNtuple->Branch("2DCosPointingAngleDaugther1",&grand_agl2D);
              VertexCompositeNtuple->Branch("2DPointingAngleDaugther1",&grand_agl2D_abs);
              VertexCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",&grand_dlos);
              VertexCompositeNtuple->Branch("3DDecayLengthDaugther1",&grand_dl);
              VertexCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",&grand_dlerror);
              VertexCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",&grand_dlos2D);
              VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2);
              VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2);
              VertexCompositeNtuple->Branch("NHitD2",&nhit2);
              VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2);
              VertexCompositeNtuple->Branch("pTD2",&pt2);
              VertexCompositeNtuple->Branch("EtaD2",&eta2);
              VertexCompositeNtuple->Branch("PhiD2",&phi2);
              VertexCompositeNtuple->Branch("pTerrD1",&ptErr2);
              VertexCompositeNtuple->Branch("pTerrD2",&ptErr2);
              VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2);
  //            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2);
              VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther1",&grand_dzos1);
              VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther2",&grand_dzos2);
              VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther1",&grand_dxyos1);
              VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther2",&grand_dxyos2);
              VertexCompositeNtuple->Branch("NHitGrandD1",&grand_nhit1);
              VertexCompositeNtuple->Branch("NHitGrandD2",&grand_nhit2);
              VertexCompositeNtuple->Branch("HighPurityGrandDaugther1",&grand_trkquality1);
              VertexCompositeNtuple->Branch("HighPurityGrandDaugther2",&grand_trkquality2);
              VertexCompositeNtuple->Branch("pTGrandD1",&grand_pt1);
              VertexCompositeNtuple->Branch("pTGrandD2",&grand_pt2);
              VertexCompositeNtuple->Branch("pTerrGrandD1",&grand_ptErr1);
              VertexCompositeNtuple->Branch("pTerrGrandD2",&grand_ptErr2);
              VertexCompositeNtuple->Branch("EtaGrandD1",&grand_eta1);
              VertexCompositeNtuple->Branch("EtaGrandD2",&grand_eta2);
              VertexCompositeNtuple->Branch("dedxHarmonic2GrandD1",&grand_H2dedx1);
              VertexCompositeNtuple->Branch("dedxHarmonic2GrandD2",&grand_H2dedx2);
          }
          else
          {
              VertexCompositeNtuple->Branch("zDCASignificanceDaugther1",&dzos1);
              VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1",&dxyos1);
              VertexCompositeNtuple->Branch("NHitD1",&nhit1);
              VertexCompositeNtuple->Branch("HighPurityDaugther1",&trkquality1);
              VertexCompositeNtuple->Branch("pTD1",&pt1);
              VertexCompositeNtuple->Branch("pTerrD1",&ptErr1);
  //            VertexCompositeNtuple->Branch("pD1",&p1);
              VertexCompositeNtuple->Branch("EtaD1",&eta1);
              VertexCompositeNtuple->Branch("PhiD1",&phi1);
  //            VertexCompositeNtuple->Branch("chargeD1",&charge1);
              VertexCompositeNtuple->Branch("dedxHarmonic2D1",&H2dedx1);
  //            VertexCompositeNtuple->Branch("dedxTruncated40Daugther1",&T4dedx1);
  //            VertexCompositeNtuple->Branch("normalizedChi2Daugther1",&trkChi1);
              VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2);
              VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2);
              VertexCompositeNtuple->Branch("NHitD2",&nhit2);
              VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2);
              VertexCompositeNtuple->Branch("pTD2",&pt2);
              VertexCompositeNtuple->Branch("pTerrD2",&ptErr2);
  //            VertexCompositeNtuple->Branch("pD2",&p2);
              VertexCompositeNtuple->Branch("EtaD2",&eta2);
              VertexCompositeNtuple->Branch("PhiD2",&phi2);
  //            VertexCompositeNtuple->Branch("chargeD2",&charge2);
              VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2);
  //            VertexCompositeNtuple->Branch("dedxTruncated40Daugther2",&T4dedx2);
  //            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2);
              if(threeProngDecay_)
              {
                VertexCompositeNtuple->Branch("zDCASignificanceDaugther3",&dzos3);
                VertexCompositeNtuple->Branch("xyDCASignificanceDaugther3",&dxyos3);
                VertexCompositeNtuple->Branch("NHitD3",&nhit3);
                VertexCompositeNtuple->Branch("HighPurityDaugther3",&trkquality3);
                VertexCompositeNtuple->Branch("pTD3",&pt1);
                VertexCompositeNtuple->Branch("pTerrD3",&ptErr3);
                VertexCompositeNtuple->Branch("EtaD3",&eta1);
                VertexCompositeNtuple->Branch("dedxHarmonic2D3",&H2dedx1);
              }
          }
          
          if(doMuon_)
          {
              VertexCompositeNtuple->Branch("OneStMuon1",&onestmuon1);
              VertexCompositeNtuple->Branch("OneStMuon2",&onestmuon2);
              VertexCompositeNtuple->Branch("PFMuon1",&pfmuon1);
              VertexCompositeNtuple->Branch("PFMuon2",&pfmuon2);
              VertexCompositeNtuple->Branch("GlbMuon1",&glbmuon1);
              VertexCompositeNtuple->Branch("GlbMuon2",&glbmuon2);
              VertexCompositeNtuple->Branch("trkMuon1",&trkmuon1);
              VertexCompositeNtuple->Branch("trkMuon2",&trkmuon2);
              VertexCompositeNtuple->Branch("caloMuon1",&calomuon1);
              VertexCompositeNtuple->Branch("caloMuon2",&calomuon2);
              VertexCompositeNtuple->Branch("SoftMuon1",&softmuon1);
              VertexCompositeNtuple->Branch("SoftMuon2",&softmuon2);

              if(doMuonFull_)
            {
              VertexCompositeNtuple->Branch("nMatchedChamberD1",&nmatchedch1);
              VertexCompositeNtuple->Branch("nMatchedStationD1",&nmatchedst1);
              VertexCompositeNtuple->Branch("EnergyDepositionD1",&matchedenergy1);
              VertexCompositeNtuple->Branch("nMatchedChamberD2",&nmatchedch2);
              VertexCompositeNtuple->Branch("nMatchedStationD2",&nmatchedst2);
              VertexCompositeNtuple->Branch("EnergyDepositionD2",&matchedenergy2);
              VertexCompositeNtuple->Branch("dx1_seg",        &dx1_seg_);
              VertexCompositeNtuple->Branch("dy1_seg",        &dy1_seg_);
              VertexCompositeNtuple->Branch("dxSig1_seg",     &dxSig1_seg_);
              VertexCompositeNtuple->Branch("dySig1_seg",     &dySig1_seg_);
              VertexCompositeNtuple->Branch("ddxdz1_seg",     &ddxdz1_seg_);
              VertexCompositeNtuple->Branch("ddydz1_seg",     &ddydz1_seg_);
              VertexCompositeNtuple->Branch("ddxdzSig1_seg",  &ddxdzSig1_seg_);
              VertexCompositeNtuple->Branch("ddydzSig1_seg",  &ddydzSig1_seg_);
              VertexCompositeNtuple->Branch("dx2_seg",        &dx2_seg_);
              VertexCompositeNtuple->Branch("dy2_seg",        &dy2_seg_);
              VertexCompositeNtuple->Branch("dxSig2_seg",     &dxSig2_seg_);
              VertexCompositeNtuple->Branch("dySig2_seg",     &dySig2_seg_);
              VertexCompositeNtuple->Branch("ddxdz2_seg",     &ddxdz2_seg_);
              VertexCompositeNtuple->Branch("ddydz2_seg",     &ddydz2_seg_);
              VertexCompositeNtuple->Branch("ddxdzSig2_seg",  &ddxdzSig2_seg_);
              VertexCompositeNtuple->Branch("ddydzSig2_seg",  &ddydzSig2_seg_);
           }
        }
    }

    } // doRecoNtuple_

    if(doGenNtuple_)
    {
        VertexCompositeNtuple->Branch("candSize_gen",&candSize_gen);
        VertexCompositeNtuple->Branch("gen_mass",&mass_gen);
        VertexCompositeNtuple->Branch("gen_pT",&pt_gen);
        VertexCompositeNtuple->Branch("gen_eta",&eta_gen);
        VertexCompositeNtuple->Branch("gen_phi",&phi_gen);
        VertexCompositeNtuple->Branch("gen_y",&y_gen);
        VertexCompositeNtuple->Branch("gen_status",&status_gen);
        VertexCompositeNtuple->Branch("gen_MotherID",&idmom);

        if(decayInGen_)
        {

            VertexCompositeNtuple->Branch("gen_DauID1",&iddau1);
            VertexCompositeNtuple->Branch("gen_DauID2",&iddau2);
            VertexCompositeNtuple->Branch("gen_DauID3",&iddau3);
        }
        if(twoLayerDecay_){
          VertexCompositeNtuple->Branch("gen_D0pT",&gen_D0pT_);
          VertexCompositeNtuple->Branch("gen_D0eta",&gen_D0eta_);
          VertexCompositeNtuple->Branch("gen_D0phi",&gen_D0phi_);
          VertexCompositeNtuple->Branch("gen_D0mass",&gen_D0mass_);
          VertexCompositeNtuple->Branch("gen_D0y",&gen_D0y_);
          VertexCompositeNtuple->Branch("gen_D0charge",&gen_D0charge_);
          VertexCompositeNtuple->Branch("gen_D0pdgId",&gen_D0pdgId_);

          VertexCompositeNtuple->Branch("gen_D0Dau1_pT",&gen_D0Dau1_pT_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_eta",&gen_D0Dau1_eta_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_phi",&gen_D0Dau1_phi_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_mass",&gen_D0Dau1_mass_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_y",&gen_D0Dau1_y_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_charge",&gen_D0Dau1_charge_);
          VertexCompositeNtuple->Branch("gen_D0Dau1_pdgId",&gen_D0Dau1_pdgId_);

          VertexCompositeNtuple->Branch("gen_D0Dau2_pT",&gen_D0Dau2_pT_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_eta",&gen_D0Dau2_eta_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_phi",&gen_D0Dau2_phi_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_mass",&gen_D0Dau2_mass_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_y",&gen_D0Dau2_y_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_charge",&gen_D0Dau2_charge_);
          VertexCompositeNtuple->Branch("gen_D0Dau2_pdgId",&gen_D0Dau2_pdgId_);

          VertexCompositeNtuple->Branch("gen_D1pT",&gen_D1pT_);
          VertexCompositeNtuple->Branch("gen_D1eta",&gen_D1eta_);
          VertexCompositeNtuple->Branch("gen_D1phi",&gen_D1phi_);
          VertexCompositeNtuple->Branch("gen_D1mass",&gen_D1mass_);
          VertexCompositeNtuple->Branch("gen_D1y",&gen_D1y_);
          VertexCompositeNtuple->Branch("gen_D1charge",&gen_D1charge_);
          VertexCompositeNtuple->Branch("gen_D1pdgId",&gen_D1pdgId_);
        }
    }
}

int VertexCompositeTreeProducer2::
muAssocToTrack( const reco::TrackRef& trackref,
                const edm::Handle<reco::MuonCollection>& muonh) const {
  auto muon = std::find_if(muonh->cbegin(),muonh->cend(),
                           [&](const reco::Muon& m) {
                             return ( m.track().isNonnull() &&
                                      m.track() == trackref    );
                           });
  return ( muon != muonh->cend() ? std::distance(muonh->cbegin(),muon) : -1 );
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
VertexCompositeTreeProducer2::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducer2);
