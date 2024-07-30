// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416

using namespace std;

class VertexCompositeNtupleProducer2 : public edm::EDAnalyzer {
public:
  explicit VertexCompositeNtupleProducer2(const edm::ParameterSet&);
  ~VertexCompositeNtupleProducer2();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();

  int muAssocToTrack( const reco::TrackRef& trackref, const edm::Handle<reco::MuonCollection>& muonh) const;


  reco::GenParticleRef findMother(const reco::GenParticleRef&);
  void genDecayLength(const reco::GenParticle&);
  void getAncestorId(const reco::GenParticle&);

  // ----------member data ---------------------------

    edm::Service<TFileService> fs;

    TTree* VertexCompositeNtuple;
    TTree* genCandidateNtuple;
    TH2F*  hMassVsMVA[6][10];
    TH2F*  hpTVsMVA[6][10];
    TH2F*  hetaVsMVA[6][10];
    TH2F*  hyVsMVA[6][10];
    TH2F*  hVtxProbVsMVA[6][10];
    TH2F*  h3DCosPointingAngleVsMVA[6][10];
    TH2F*  h3DPointingAngleVsMVA[6][10];
    TH2F*  h2DCosPointingAngleVsMVA[6][10];
    TH2F*  h2DPointingAngleVsMVA[6][10];
    TH2F*  h3DDecayLengthSignificanceVsMVA[6][10];
    TH2F*  h3DDecayLengthVsMVA[6][10];
    TH2F*  h2DDecayLengthSignificanceVsMVA[6][10];
    TH2F*  h2DDecayLengthVsMVA[6][10];
    TH2F*  h3DDCAVsMVA[6][10];
    TH2F*  h2DDCAVsMVA[6][10];
    TH2F*  hzDCASignificanceDaugther1VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther1VsMVA[6][10];
    TH2F*  hNHitD1VsMVA[6][10];
    TH2F*  hpTD1VsMVA[6][10];
    TH2F*  hpTerrD1VsMVA[6][10];
    TH2F*  hEtaD1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsP[6][10];
    TH2F*  hzDCASignificanceDaugther2VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther2VsMVA[6][10];
    TH2F*  hNHitD2VsMVA[6][10];
    TH2F*  hpTD2VsMVA[6][10];
    TH2F*  hpTerrD2VsMVA[6][10];
    TH2F*  hEtaD2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsP[6][10];
    TH2F*  hzDCASignificanceDaugther3VsMVA[6][10];
    TH2F*  hxyDCASignificanceDaugther3VsMVA[6][10];
    TH2F*  hNHitD3VsMVA[6][10];
    TH2F*  hpTD3VsMVA[6][10];
    TH2F*  hpTerrD3VsMVA[6][10];
    TH2F*  hEtaD3VsMVA[6][10];
    TH2F*  hdedxHarmonic2D3VsMVA[6][10];
    TH2F*  hdedxHarmonic2D3VsP[6][10];

    bool   saveTree_;
    bool   saveHistogram_;
    bool   saveAllHistogram_;
    double massHistPeak_;
    double massHistWidth_;
    int    massHistBins_;

    //options
    bool doRecoNtuple_;
    bool doGenNtuple_;
    bool doGenMatching_;
    bool doGenMatchingTOF_;
    bool hasSwap_;
    bool decayInGen_;
    bool twoLayerDecay_;
    bool threeProngDecay_;
    bool doMuon_;
    bool doMuonFull_;
    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    int PID_dau3_;

    //cut variables
    double multMax_;
    double multMin_;
    double deltaR_; //deltaR for Gen matching

    vector<double> pTBins_;
    vector<double> yBins_;

    //tree branches
    //event info
    int centrality;
    int Ntrkoffline;
    int Npixel;
    float HFsumETPlus;
    float HFsumETMinus;
    float ZDCPlus;
    float ZDCMinus;
    float bestvx;
    float bestvy;
    float bestvz;

    //Composite candidate info
    float mva;
    float pt;
    float eta;
    float phi;
    float flavor;
    float y;
    float mass;
    float VtxProb;
    float dlos;
    float dl;
    float dlerror;
    float agl;
    float dca3D;
    float dcaErr3D;
    float vtxChi2;
    float ndf;
    float agl_abs;
    float agl2D;
    float agl2D_abs;
    float dlos2D;
    float dl2D;
    bool isSwap;
    bool matchGEN;
    int idmom_reco;

    //dau candidate info
    float grand_mass;
    float grand_VtxProb;
    float grand_dlos;
    float grand_dl;
    float grand_dlerror;
    float grand_agl;
    float grand_vtxChi2;
    float grand_ndf;
    float grand_agl_abs;
    float grand_agl2D;
    float grand_agl2D_abs;
    float grand_dlos2D;

    //dau info
    float dzos1;
    float dzos2;
    float dzos3;
    float dxyos1;
    float dxyos2;
    float dxyos3;
    float nhit1;
    float nhit2;
    float nhit3;
    bool trkquality1;
    bool trkquality2;
    bool trkquality3;
    float pt1;
    float pt2;
    float pt3;
    float ptErr1;
    float ptErr2;
    float ptErr3;
    float p1;
    float p2;
    float p3;
    float eta1;
    float eta2;
    float eta3;
    float phi1;
    float phi2;
    float phi3;
    int charge1;
    int charge2;
    int charge3;
    int pid1;
    int pid2;
    int pid3;
    float tof1;
    float tof2;
    float tof3;
    float H2dedx1;
    float H2dedx2;
    float H2dedx3;
    float T4dedx1;
    float T4dedx2;
    float T4dedx3;
    float trkChi1;
    float trkChi2;
    float trkChi3;

    //grand-dau info
    float grand_dzos1;
    float grand_dzos2;
    float grand_dxyos1;
    float grand_dxyos2;
    float grand_nhit1;
    float grand_nhit2;
    bool grand_trkquality1;
    bool grand_trkquality2;
    float grand_pt1;
    float grand_pt2;
    float grand_ptErr1;
    float grand_ptErr2;
    float grand_p1;
    float grand_p2;
    float grand_eta1;
    float grand_eta2;
    int grand_charge1;
    int grand_charge2;
    float grand_H2dedx1;
    float grand_H2dedx2;
    float grand_T4dedx1;
    float grand_T4dedx2;
    float grand_trkChi1;
    float grand_trkChi2;

    //dau muon info
    bool  onestmuon1;
    bool  onestmuon2;
    bool  pfmuon1;
    bool  pfmuon2;
    bool  glbmuon1;
    bool  glbmuon2;
    bool  trkmuon1;
    bool  trkmuon2;
    bool  calomuon1;
    bool  calomuon2;
    bool  softmuon1;
    bool  softmuon2;
    float nmatchedst1;
    float nmatchedch1;
    float ntrackerlayer1;
    float npixellayer1;
    float matchedenergy1;
    float nmatchedst2;
    float nmatchedch2;
    float ntrackerlayer2;
    float npixellayer2;
    float matchedenergy2;
    float dx1_seg_;
    float dy1_seg_;
    float dxSig1_seg_;
    float dySig1_seg_;
    float ddxdz1_seg_;
    float ddydz1_seg_;
    float ddxdzSig1_seg_;
    float ddydzSig1_seg_;
    float dx2_seg_;
    float dy2_seg_;
    float dxSig2_seg_;
    float dySig2_seg_;
    float ddxdz2_seg_;
    float ddydz2_seg_;
    float ddxdzSig2_seg_;
    float ddydzSig2_seg_;

    // gen info
    float pt_gen;
    float eta_gen;
    float phi_gen;
    int status_gen;
    int idmom;
    float ptmom;
    float ymom;
    float etamom;
    float phimom;
    int statusmom;
    float y_gen;
    int iddau1;
    int iddau2;
    int iddau3;

  // gen information for # of daughters == 2
    int   gen_ancestorFlavor_;
    int   gen_ancestorId_;
  float gen_PVx_;
  float gen_PVy_;
  float gen_PVz_;

  float gen_pT_;
  float gen_eta_;
  float gen_phi_;
  float gen_mass_;
  float gen_y_;
  float gen_charge_;
  float gen_pdgId_;

  float gen_decayLength3D_;
  float gen_decayLength2D_;
  float gen_angle3D_;
  float gen_angle2D_;

  float gen_pTD1_;
  float gen_etaD1_;
  float gen_phiD1_;
  float gen_massD1_;
  float gen_yD1_;
  float gen_chargeD1_;
  float gen_pdgIdD1_;

  float gen_pTD2_;
  float gen_etaD2_;
  float gen_phiD2_;
  float gen_massD2_;
  float gen_yD2_;
  float gen_chargeD2_;
  float gen_pdgIdD2_;

  reco::Particle::Point genVertex_;

  // options
    bool useAnyMVA_;
    bool useDCA_;
    bool isSkimMVA_;
    bool isCentrality_;

    edm::Handle<int> cbin_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;

    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

    edm::EDGetTokenT<int> tok_centBinLabel_;
    edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

    // for DCA
    edm::EDGetTokenT<std::vector<float > > tok_DCAVal_;
    edm::EDGetTokenT<std::vector<float > > tok_DCAErr_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// static data member definitions
//

//
// constructors and destructor
//

VertexCompositeNtupleProducer2::VertexCompositeNtupleProducer2(const edm::ParameterSet& iConfig){
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
    if(isCentrality_) {
      tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
      tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
    }

    if(useAnyMVA_ && iConfig.exists("MVACollection"))
      MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
    if(iConfig.exists("DCAValCollection") && iConfig.exists("DCAErrCollection")) {
      useDCA_ = true;
      tok_DCAVal_ = consumes<std::vector<float > >(iConfig.getParameter<edm::InputTag>("DCAValCollection"));
      tok_DCAErr_ = consumes<std::vector<float > >(iConfig.getParameter<edm::InputTag>("DCAErrCollection"));
    }
};


VertexCompositeNtupleProducer2::~VertexCompositeNtupleProducer2() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

};

void VertexCompositeNtupleProducer2::beginJob(){
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
};

void VertexCompositeNtupleProducer2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using std::vector;
    using namespace edm;
    using namespace reco;

    if(doGenNtuple_) fillGEN(iEvent,iSetup);
    if(doRecoNtuple_) fillRECO(iEvent,iSetup);
};

void VertexCompositeNtupleProducer2::endJob() {

};

reco::GenParticleRef VertexCompositeNtupleProducer2::findMother(const reco::GenParticleRef& genParRef) {
  if(genParRef.isNull()) return genParRef;
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId(); const int pdg_OLD = pdg;
  while(pdg==pdg_OLD && genMomRef->numberOfMothers()>0)
  {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  if(pdg==pdg_OLD) genMomRef = reco::GenParticleRef();
  return genMomRef;
};

void VertexCompositeNtupleProducer2::genDecayLength(const reco::GenParticle& gCand) {
  gen_decayLength2D_ = -99.;
  gen_decayLength3D_ = -99.;
  gen_angle2D_ = -99;
  gen_angle3D_ = -99;

  if(gCand.numberOfDaughters()==0 || !gCand.daughter(0)) return;
  const auto& dauVtx = gCand.daughter(0)->vertex();
  TVector3 ptosvec(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), dauVtx.Z() - genVertex_.z());
  TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
  gen_angle3D_ = secvec.Angle(ptosvec);
  gen_decayLength3D_ = ptosvec.Mag();
  TVector3 ptosvec2D(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), 0.0);
  TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
  gen_angle2D_ = secvec2D.Angle(ptosvec2D);
  gen_decayLength2D_ = ptosvec2D.Mag();
};

void VertexCompositeNtupleProducer2::getAncestorId(const reco::GenParticle& gCand) {
  gen_ancestorId_ = 0;
  gen_ancestorFlavor_ = 0;
  for (auto mothers = gCand.motherRefVector();
      !mothers.empty(); ) {
    auto mom = mothers.at(0);
    mothers = mom->motherRefVector();
    gen_ancestorId_ = mom->pdgId();
    const auto idstr = std::to_string(std::abs(gen_ancestorId_));
    gen_ancestorFlavor_ = std::stoi(std::string{idstr.begin(), idstr.begin()+1});
    if (idstr[0] == '5') {
      break;
    }
    if (std::abs(gen_ancestorId_) <= 40) break;
  }
};

int VertexCompositeNtupleProducer2::muAssocToTrack( const reco::TrackRef& trackref, const edm::Handle<reco::MuonCollection>& muonh) const {
  auto muon = std::find_if(muonh->cbegin(),muonh->cend(),
                           [&](const reco::Muon& m) {
                             return ( m.track().isNonnull() &&
                                      m.track() == trackref    );
                           });
  return ( muon != muonh->cend() ? std::distance(muonh->cbegin(),muon) : -1 );
};

void VertexCompositeNtupleProducer2::initTree() {
    VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");
    genCandidateNtuple = fs->make< TTree>("genCandidateNtuple","genCandidateNtuple");

    VertexCompositeNtuple->Branch("pT",&pt,"pT/F");
    VertexCompositeNtuple->Branch("y",&y,"y/F");
    VertexCompositeNtuple->Branch("eta",&eta,"eta/F");
    VertexCompositeNtuple->Branch("phi",&phi,"phi/F");
    VertexCompositeNtuple->Branch("mass",&mass,"mass/F");

    if(useAnyMVA_) VertexCompositeNtuple->Branch("mva",&mva,"mva/F");
    if(useDCA_) {
      VertexCompositeNtuple->Branch("dca3D", &dca3D, "dca3D/F");
      VertexCompositeNtuple->Branch("dcaErr3D", &dcaErr3D, "dcaErr3D/F");
    }
    if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality,"centrality/I");
    if(!isSkimMVA_) {
        //Event info
        VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
        VertexCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
        VertexCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
        VertexCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
        VertexCompositeNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
        VertexCompositeNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
        VertexCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
        VertexCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
        VertexCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");

        //Composite candidate info RECO
        VertexCompositeNtuple->Branch("flavor",&flavor,"flavor/F");
        VertexCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb/F");
      //  VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2/F");
      //  VertexCompositeNtuple->Branch("VtxNDF",&ndf,"VtxNDF/F");
        VertexCompositeNtuple->Branch("3DCosPointingAngle",&agl,"3DCosPointingAngle/F");
        VertexCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle/F");
        VertexCompositeNtuple->Branch("2DCosPointingAngle",&agl2D,"2DCosPointingAngle/F");
        VertexCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle/F");
        VertexCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance/F");
        VertexCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength/F");
      //  VertexCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError/F");
        VertexCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance/F");
        VertexCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength/F");

        if(doGenMatching_) {
            VertexCompositeNtuple->Branch("isSwap",&isSwap,"isSwap/O");
            VertexCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco/I");
            VertexCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN/O");

            if (!twoLayerDecay_ && !threeProngDecay_) {
              VertexCompositeNtuple->Branch("gen_ancestorFlavor", &gen_ancestorFlavor_);
              VertexCompositeNtuple->Branch("gen_ancestorId", &gen_ancestorId_);

              VertexCompositeNtuple->Branch("gen_PVx_", &gen_PVx_);
              VertexCompositeNtuple->Branch("gen_PVy", &gen_PVy_);
              VertexCompositeNtuple->Branch("gen_PVz", &gen_PVz_);

              VertexCompositeNtuple->Branch("gen_pT", &gen_pT_);
              VertexCompositeNtuple->Branch("gen_eta", &gen_eta_);
              VertexCompositeNtuple->Branch("gen_phi", &gen_phi_);
              VertexCompositeNtuple->Branch("gen_mass", &gen_mass_);
              VertexCompositeNtuple->Branch("gen_y", &gen_y_);

              VertexCompositeNtuple->Branch("gen_decayLength3D", &gen_decayLength3D_);
              VertexCompositeNtuple->Branch("gen_decayLength2D", &gen_decayLength2D_);
              VertexCompositeNtuple->Branch("gen_angle3D", &gen_angle3D_);
              VertexCompositeNtuple->Branch("gen_angle2D", &gen_angle2D_);

              VertexCompositeNtuple->Branch("gen_pTD1", &gen_pTD1_);
              VertexCompositeNtuple->Branch("gen_etaD1", &gen_etaD1_);
              VertexCompositeNtuple->Branch("gen_phiD1", &gen_phiD1_);
              VertexCompositeNtuple->Branch("gen_massD1", &gen_massD1_);
              VertexCompositeNtuple->Branch("gen_yD1", &gen_yD1_);
              VertexCompositeNtuple->Branch("gen_chargeD1", &gen_chargeD1_);
              VertexCompositeNtuple->Branch("gen_pdgIdD1", &gen_pdgIdD1_);

              VertexCompositeNtuple->Branch("gen_pTD2", &gen_pTD2_);
              VertexCompositeNtuple->Branch("gen_etaD2", &gen_etaD2_);
              VertexCompositeNtuple->Branch("gen_phiD2", &gen_phiD2_);
              VertexCompositeNtuple->Branch("gen_massD2", &gen_massD2_);
              VertexCompositeNtuple->Branch("gen_yD2", &gen_yD2_);
              VertexCompositeNtuple->Branch("gen_chargeD2", &gen_chargeD2_);
              VertexCompositeNtuple->Branch("gen_pdgIdD2", &gen_pdgIdD2_);
            }
        }

        if(doGenMatchingTOF_) {
          VertexCompositeNtuple->Branch("PIDD1",&pid1,"PIDD1/I");
          VertexCompositeNtuple->Branch("PIDD2",&pid1,"PIDD2/I");
          VertexCompositeNtuple->Branch("TOFD1",&tof1,"TOFD1/F");
          VertexCompositeNtuple->Branch("TOFD2",&tof1,"TOFD2/F");
        }

        //daugther & grand daugther info
        if(twoLayerDecay_) {
            VertexCompositeNtuple->Branch("massDaugther1",&grand_mass,"massDaugther1/F");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1/F");
            VertexCompositeNtuple->Branch("PhiD1",&phi1,"PhiD1/F");
            VertexCompositeNtuple->Branch("VtxProbDaugther1",&grand_VtxProb,"VtxProbDaugther1/F");
            VertexCompositeNtuple->Branch("VtxChi2Daugther1",&grand_vtxChi2,"VtxChi2Daugther1/F");
            VertexCompositeNtuple->Branch("VtxNDFDaugther1",&grand_ndf,"VtxNDFDaugther1/F");
            VertexCompositeNtuple->Branch("3DCosPointingAngleDaugther1",&grand_agl,"3DCosPointingAngleDaugther1/F");
            VertexCompositeNtuple->Branch("3DPointingAngleDaugther1",&grand_agl_abs,"3DPointingAngleDaugther1/F");
            VertexCompositeNtuple->Branch("2DCosPointingAngleDaugther1",&grand_agl2D,"2DCosPointingAngleDaugther1/F");
            VertexCompositeNtuple->Branch("2DPointingAngleDaugther1",&grand_agl2D_abs,"2DPointingAngleDaugther1/F");
            VertexCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",&grand_dlos,"3DDecayLengthSignificanceDaugther1/F");
            VertexCompositeNtuple->Branch("3DDecayLengthDaugther1",&grand_dl,"3DDecayLengthDaugther1/F");
            VertexCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",&grand_dlerror,"3DDecayLengthErrorDaugther1/F");
            VertexCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",&grand_dlos2D,"2DDecayLengthSignificanceDaugther1/F");
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2,"zDCASignificanceDaugther2/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2,"xyDCASignificanceDaugther2/F");
            VertexCompositeNtuple->Branch("NHitD2",&nhit2,"NHitD2/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2,"HighPurityDaugther2/O");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2/F");
            VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2/F");
            VertexCompositeNtuple->Branch("pD2",&p2,"pD2/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2/F");
            VertexCompositeNtuple->Branch("PhiD2",&phi2,"PhiD2/F");
            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2/F");
            VertexCompositeNtuple->Branch("dedxTruncated40Daugther2",&T4dedx2,"dedxTruncated40Daugther2/F");
            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2,"normalizedChi2Daugther2/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther1",&grand_dzos1,"zDCASignificanceGrandDaugther1/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther2",&grand_dzos2,"zDCASignificanceGrandDaugther2/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther1",&grand_dxyos1,"xyDCASignificanceGrandDaugther1/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther2",&grand_dxyos2,"xyDCASignificanceGrandDaugther2/F");
            VertexCompositeNtuple->Branch("NHitGrandD1",&grand_nhit1,"NHitGrandD1/F");
            VertexCompositeNtuple->Branch("NHitGrandD2",&grand_nhit2,"NHitGrandD2/F");
            VertexCompositeNtuple->Branch("HighPurityGrandDaugther1",&grand_trkquality1,"HighPurityGrandDaugther1/O");
            VertexCompositeNtuple->Branch("HighPurityGrandDaugther2",&grand_trkquality2,"HighPurityGrandDaugther2/O");
            VertexCompositeNtuple->Branch("pTGrandD1",&grand_pt1,"pTGrandD1/F");
            VertexCompositeNtuple->Branch("pTGrandD2",&grand_pt2,"pTGrandD2/F");
            VertexCompositeNtuple->Branch("pTerrGrandD1",&grand_ptErr1,"pTerrGrandD1/F");
            VertexCompositeNtuple->Branch("pTerrGrandD2",&grand_ptErr2,"pTerrGrandD2/F");
            VertexCompositeNtuple->Branch("pGrandD1",&grand_p1,"pGrandD1/F");
            VertexCompositeNtuple->Branch("pGrandD2",&grand_p2,"pGrandD2/F");
            VertexCompositeNtuple->Branch("EtaGrandD1",&grand_eta1,"EtaGrandD1/F");
            VertexCompositeNtuple->Branch("EtaGrandD2",&grand_eta2,"EtaGrandD2/F");
            VertexCompositeNtuple->Branch("chargeGrandD1",&grand_charge1,"chargeGrandD1/I");
            VertexCompositeNtuple->Branch("chargeGrandD2",&grand_charge2,"chargeGrandD2/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD1",&grand_H2dedx1,"dedxHarmonic2GrandD1/F");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD2",&grand_H2dedx2,"dedxHarmonic2GrandD2/F");
            VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther1",&grand_T4dedx1,"dedxTruncated40GrandDaugther1/F");
            VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther2",&grand_T4dedx2,"dedxTruncated40GrandDaugther2/F");
            VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther1",&grand_trkChi1,"normalizedChi2GrandDaugther1/F");
            VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther2",&grand_trkChi2,"normalizedChi2GrandDaugther2/F");
        }
        else {
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther1",&dzos1,"zDCASignificanceDaugther1/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1",&dxyos1,"xyDCASignificanceDaugther1/F");
            VertexCompositeNtuple->Branch("NHitD1",&nhit1,"NHitD1/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther1",&trkquality1,"HighPurityDaugther1/O");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1/F");
            VertexCompositeNtuple->Branch("pTerrD1",&ptErr1,"pTerrD1/F");
            VertexCompositeNtuple->Branch("pD1",&p1,"pD1/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1/F");
            VertexCompositeNtuple->Branch("PhiD1",&eta1,"PhiD1/F");
            VertexCompositeNtuple->Branch("chargeD1",&charge1,"chargeD1/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D1",&H2dedx1,"dedxHarmonic2D1/F");
            VertexCompositeNtuple->Branch("dedxTruncated40Daugther1",&T4dedx1,"dedxTruncated40Daugther1/F");
            VertexCompositeNtuple->Branch("normalizedChi2Daugther1",&trkChi1,"normalizedChi2Daugther1/F");
            VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2,"zDCASignificanceDaugther2/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2,"xyDCASignificanceDaugther2/F");
            VertexCompositeNtuple->Branch("NHitD2",&nhit2,"NHitD2/F");
            VertexCompositeNtuple->Branch("HighPurityDaugther2",&trkquality2,"HighPurityDaugther2/O");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2/F");
            VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2/F");
            VertexCompositeNtuple->Branch("pD2",&p2,"pD2/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2/F");
            VertexCompositeNtuple->Branch("PhiD2",&eta2,"PhiD2/F");
            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2/F");
            VertexCompositeNtuple->Branch("dedxTruncated40Daugther2",&T4dedx2,"dedxTruncated40Daugther2/F");
            VertexCompositeNtuple->Branch("normalizedChi2Daugther2",&trkChi2,"normalizedChi2Daugther2/F");
          if(threeProngDecay_)
            {
              VertexCompositeNtuple->Branch("zDCASignificanceDaugther3",&dzos3,"zDCASignificanceDaugther3/F");
              VertexCompositeNtuple->Branch("xyDCASignificanceDaugther3",&dxyos3,"xyDCASignificanceDaugther3/F");
              VertexCompositeNtuple->Branch("NHitD3",&nhit3,"NHitD3/F");
              VertexCompositeNtuple->Branch("HighPurityDaugther3",&trkquality3,"HighPurityDaugther3/O");
              VertexCompositeNtuple->Branch("pTD3",&pt1,"pTD3/F");
              VertexCompositeNtuple->Branch("pTerrD3",&ptErr3,"pTerrD3/F");
              VertexCompositeNtuple->Branch("EtaD3",&eta1,"EtaD3/F");
              VertexCompositeNtuple->Branch("dedxHarmonic2D3",&H2dedx1,"dedxHarmonic2D3/F");
            }
        }

        if(doMuon_) {
            VertexCompositeNtuple->Branch("OneStMuon1",&onestmuon1,"OneStMuon1/O");
            VertexCompositeNtuple->Branch("OneStMuon2",&onestmuon2,"OneStMuon2/O");
            VertexCompositeNtuple->Branch("PFMuon1",&pfmuon1,"PFMuon1/O");
            VertexCompositeNtuple->Branch("PFMuon2",&pfmuon2,"PFMuon2/O");
            VertexCompositeNtuple->Branch("GlbMuon1",&glbmuon1,"GlbMuon1/O");
            VertexCompositeNtuple->Branch("GlbMuon2",&glbmuon2,"GlbMuon2/O");
            VertexCompositeNtuple->Branch("trkMuon1",&trkmuon1,"trkMuon1/O");
            VertexCompositeNtuple->Branch("trkMuon2",&trkmuon2,"trkMuon2/O");
            VertexCompositeNtuple->Branch("caloMuon1",&calomuon1,"caloMuon1/O");
            VertexCompositeNtuple->Branch("caloMuon2",&calomuon2,"caloMuon2/O");
            VertexCompositeNtuple->Branch("SoftMuon1",&softmuon1,"SoftMuon1/O");
            VertexCompositeNtuple->Branch("SoftMuon2",&softmuon2,"SoftMuon2/O");
            if(doMuonFull_)
            {
              VertexCompositeNtuple->Branch("nMatchedChamberD1",&nmatchedch1,"nMatchedChamberD1/F");
              VertexCompositeNtuple->Branch("nMatchedStationD1",&nmatchedst1,"nMatchedStationD1/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD1",&matchedenergy1,"EnergyDepositionD1/F");
              VertexCompositeNtuple->Branch("nMatchedChamberD2",&nmatchedch2,"nMatchedChamberD2/F");
              VertexCompositeNtuple->Branch("nMatchedStationD2",&nmatchedst2,"nMatchedStationD2/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD2",&matchedenergy2,"EnergyDepositionD2/F");
              VertexCompositeNtuple->Branch("dx1_seg",        &dx1_seg_, "dx1_seg/F");
              VertexCompositeNtuple->Branch("dy1_seg",        &dy1_seg_, "dy1_seg/F");
              VertexCompositeNtuple->Branch("dxSig1_seg",     &dxSig1_seg_, "dxSig1_seg/F");
              VertexCompositeNtuple->Branch("dySig1_seg",     &dySig1_seg_, "dySig1_seg/F");
              VertexCompositeNtuple->Branch("ddxdz1_seg",     &ddxdz1_seg_, "ddxdz1_seg/F");
              VertexCompositeNtuple->Branch("ddydz1_seg",     &ddydz1_seg_, "ddydz1_seg/F");
              VertexCompositeNtuple->Branch("ddxdzSig1_seg",  &ddxdzSig1_seg_, "ddxdzSig1_seg/F");
              VertexCompositeNtuple->Branch("ddydzSig1_seg",  &ddydzSig1_seg_, "ddydzSig1_seg/F");
              VertexCompositeNtuple->Branch("dx2_seg",        &dx2_seg_, "dx2_seg/F");
              VertexCompositeNtuple->Branch("dy2_seg",        &dy2_seg_, "dy2_seg/F");
              VertexCompositeNtuple->Branch("dxSig2_seg",     &dxSig2_seg_, "dxSig2_seg/F");
              VertexCompositeNtuple->Branch("dySig2_seg",     &dySig2_seg_, "dySig2_seg/F");
              VertexCompositeNtuple->Branch("ddxdz2_seg",     &ddxdz2_seg_, "ddxdz2_seg/F");
              VertexCompositeNtuple->Branch("ddydz2_seg",     &ddydz2_seg_, "ddydz2_seg/F");
              VertexCompositeNtuple->Branch("ddxdzSig2_seg",  &ddxdzSig2_seg_, "ddxdzSig2_seg/F");
              VertexCompositeNtuple->Branch("ddydzSig2_seg",  &ddydzSig2_seg_, "ddydzSig2_seg/F");
           }
        }
    }

    if(doGenNtuple_) {
        genCandidateNtuple->Branch("pT_gen",&pt_gen,"pT_gen/F");
        genCandidateNtuple->Branch("eta_gen",&eta_gen,"eta_gen/F");
        genCandidateNtuple->Branch("y_gen",&y_gen,"y_gen/F");
        genCandidateNtuple->Branch("status_gen",&status_gen,"status_gen/I");
        genCandidateNtuple->Branch("MotherID_gen",&idmom,"MotherID_gen/I");
        genCandidateNtuple->Branch("MotherPt_gen",&ptmom,"MotherPt_gen/F");
        genCandidateNtuple->Branch("MotherEta_gen",&etamom,"MotherEta_gen/F");
        genCandidateNtuple->Branch("MotherPhi_gen",&phimom,"MotherPhi_gen/F");
        genCandidateNtuple->Branch("MotherY_gen",&ymom,"MotherY_gen/F");
        genCandidateNtuple->Branch("MotherStatus_gen",&statusmom,"MotherStatus_gen/I");

        if(decayInGen_) {
            genCandidateNtuple->Branch("DauID1_gen",&iddau1,"DauID1_gen/I");
            genCandidateNtuple->Branch("DauID2_gen",&iddau2,"DauID2_gen/I");
            genCandidateNtuple->Branch("DauID3_gen",&iddau3,"DauID3_gen/I");
        }

        genCandidateNtuple->Branch("gen_ancestorFlavor", &gen_ancestorFlavor_);
        genCandidateNtuple->Branch("gen_ancestorId", &gen_ancestorId_);

        genCandidateNtuple->Branch("gen_PVx_", &gen_PVx_);
        genCandidateNtuple->Branch("gen_PVy", &gen_PVy_);
        genCandidateNtuple->Branch("gen_PVz", &gen_PVz_);

        genCandidateNtuple->Branch("gen_pT", &gen_pT_);
        genCandidateNtuple->Branch("gen_eta", &gen_eta_);
        genCandidateNtuple->Branch("gen_phi", &gen_phi_);
        genCandidateNtuple->Branch("gen_mass", &gen_mass_);
        genCandidateNtuple->Branch("gen_y", &gen_y_);

        genCandidateNtuple->Branch("gen_decayLength3D", &gen_decayLength3D_);
        genCandidateNtuple->Branch("gen_decayLength2D", &gen_decayLength2D_);
        genCandidateNtuple->Branch("gen_angle3D", &gen_angle3D_);
        genCandidateNtuple->Branch("gen_angle2D", &gen_angle2D_);

        genCandidateNtuple->Branch("gen_pTD1", &gen_pTD1_);
        genCandidateNtuple->Branch("gen_etaD1", &gen_etaD1_);
        genCandidateNtuple->Branch("gen_phiD1", &gen_phiD1_);
        genCandidateNtuple->Branch("gen_massD1", &gen_massD1_);
        genCandidateNtuple->Branch("gen_yD1", &gen_yD1_);
        genCandidateNtuple->Branch("gen_chargeD1", &gen_chargeD1_);
        genCandidateNtuple->Branch("gen_pdgIdD1", &gen_pdgIdD1_);

        genCandidateNtuple->Branch("gen_pTD2", &gen_pTD2_);
        genCandidateNtuple->Branch("gen_etaD2", &gen_etaD2_);
        genCandidateNtuple->Branch("gen_phiD2", &gen_phiD2_);
        genCandidateNtuple->Branch("gen_massD2", &gen_massD2_);
        genCandidateNtuple->Branch("gen_yD2", &gen_yD2_);
        genCandidateNtuple->Branch("gen_chargeD2", &gen_chargeD2_);
        genCandidateNtuple->Branch("gen_pdgIdD2", &gen_pdgIdD2_);
    }
};
void VertexCompositeNtupleProducer2::initHistogram() {
  for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++) {
    for(unsigned int iy=0;iy<yBins_.size()-1;iy++) {
      hMassVsMVA[iy][ipt] = fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d",iy,ipt),";mva;mass(GeV)",100,-1.,1.,massHistBins_,massHistPeak_-massHistWidth_,massHistPeak_+massHistWidth_);
  //  h3DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;3D DCA;",100,-1.,1.,1000,0,10);
  //  h2DDCAVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;2D DCA;",100,-1.,1.,1000,0,10);

      if(saveAllHistogram_) {
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

        if(threeProngDecay_) {
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
};

//define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeNtupleProducer2);