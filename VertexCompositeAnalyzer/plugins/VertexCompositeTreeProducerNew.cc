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
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416
#define MAXCAN 50000

using namespace std;

class VertexCompositeTreeProducerNew : public edm::EDAnalyzer {
public:
  explicit VertexCompositeTreeProducerNew(const edm::ParameterSet&);
  ~VertexCompositeTreeProducerNew();

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

  // ----------member data ---------------------------
    
    edm::Service<TFileService> fs;

    TTree* VertexCompositeNtuple;
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
    TH2F*  hzDCASignificancedaughter1VsMVA[6][10];
    TH2F*  hxyDCASignificancedaughter1VsMVA[6][10];
    TH2F*  hNHitD1VsMVA[6][10];
    TH2F*  hpTD1VsMVA[6][10];
    TH2F*  hpTerrD1VsMVA[6][10];
    TH2F*  hEtaD1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsMVA[6][10];
    TH2F*  hdedxHarmonic2D1VsP[6][10];
    TH2F*  hzDCASignificancedaughter2VsMVA[6][10];
    TH2F*  hxyDCASignificancedaughter2VsMVA[6][10];
    TH2F*  hNHitD2VsMVA[6][10];
    TH2F*  hpTD2VsMVA[6][10];
    TH2F*  hpTerrD2VsMVA[6][10];
    TH2F*  hEtaD2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsMVA[6][10];
    TH2F*  hdedxHarmonic2D2VsP[6][10];
    TH2F*  hzDCASignificancedaughter3VsMVA[6][10];
    TH2F*  hxyDCASignificancedaughter3VsMVA[6][10];
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
    bool doGenDoubleDecay_;
    bool hasSwap_;
    bool decayInGen_;
    bool twoLayerDecay_;
    bool doubleCand_;
    bool threeProngDecay_;
    bool doMuon_;
    bool doMuonFull_;
    bool debug_;
    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    int PID_dau3_;

    int PID_dau1_grand1_;
    int PID_dau1_grand2_;
    int PID_dau2_grand1_;
    int PID_dau2_grand2_;
    
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
    float HFsumET;
    float bestvx;
    float bestvy;
    float bestvz;
    int candSize;
    float ephfpAngle[3];
    float ephfmAngle[3];
    float ephfpQ[3];
    float ephfmQ[3];
    float ephfpSumW;
    float ephfmSumW;
    
    //Composite candidate info
    float mva[MAXCAN];
    float pt[MAXCAN];
    float eta[MAXCAN];
    float phi[MAXCAN];
    float flavor[MAXCAN];
    float flavor1[MAXCAN]; // For double D decay 
    float flavor2[MAXCAN]; // For double D decay
    float y[MAXCAN];
    float mass[MAXCAN];
    float VtxProb[MAXCAN];
    float dlos[MAXCAN];
    float dl[MAXCAN];
    float dlerror[MAXCAN];
    float agl[MAXCAN];
    float dca3D[MAXCAN];
    float dcaErr3D[MAXCAN];
    float vtxChi2[MAXCAN];
    float ndf[MAXCAN];
    float agl_abs[MAXCAN];
    float agl2D[MAXCAN];
    float agl2D_abs[MAXCAN];
    float dlos2D[MAXCAN];
    float dl2D[MAXCAN];
    bool isSwap[MAXCAN];
    bool matchGEN[MAXCAN];
    int idmom_reco[MAXCAN];

    bool matchGEN1[MAXCAN]; // For Double Decay
    bool matchGEN2[MAXCAN]; // For Double Decay
    bool isSwap1[MAXCAN]; // For Double Decay
    bool isSwap2[MAXCAN]; // For Double Decay
    
    //dau candidate info
    float grand_mass[MAXCAN];
    float grand_VtxProb[MAXCAN];
    float grand_dlos[MAXCAN];
    float grand_dl[MAXCAN];
    float grand_dlerror[MAXCAN];
    float grand_agl[MAXCAN];
    float grand_vtxChi2[MAXCAN];
    float grand_ndf[MAXCAN];
    float grand_agl_abs[MAXCAN];
    float grand_agl2D[MAXCAN];
    float grand_agl2D_abs[MAXCAN];
    float grand_dlos2D[MAXCAN];
    float mva1[MAXCAN];

    float grand_mass2[MAXCAN];
    float grand_VtxProb2[MAXCAN];
    float grand_dlos2[MAXCAN];
    float grand_dl2[MAXCAN];
    float grand_dlerror2[MAXCAN];
    float grand_agl2[MAXCAN];
    float grand_vtxChi22[MAXCAN];
    float grand_ndf2[MAXCAN];
    float grand_agl_abs2[MAXCAN];
    float grand_agl2D2[MAXCAN];
    float grand_agl2D_abs2[MAXCAN];
    float grand_dlos2D2[MAXCAN];
    float mva2[MAXCAN];

    //dau info
    float dzos1[MAXCAN];
    float dzos2[MAXCAN];
    float dzos3[MAXCAN];
    float dxyos1[MAXCAN];
    float dxyos2[MAXCAN];
    float dxyos3[MAXCAN];
    float nhit1[MAXCAN];
    float nhit2[MAXCAN];
    float nhit3[MAXCAN];
    bool trkquality1[MAXCAN];
    bool trkquality2[MAXCAN];
    bool trkquality3[MAXCAN];
    float pt1[MAXCAN];
    float pt2[MAXCAN];
    float pt3[MAXCAN];
    float ptErr1[MAXCAN];
    float ptErr2[MAXCAN];
    float ptErr3[MAXCAN];
    float p1[MAXCAN];
    float p2[MAXCAN];
    float p3[MAXCAN];
    float eta1[MAXCAN];
    float eta2[MAXCAN];
    float eta3[MAXCAN];
    float phi1[MAXCAN];
    float phi2[MAXCAN];
    float phi3[MAXCAN];
    int charge1[MAXCAN];
    int charge2[MAXCAN];
    int charge3[MAXCAN];
    int pid1[MAXCAN];
    int pid2[MAXCAN];
    int pid3[MAXCAN];
    float tof1[MAXCAN];
    float tof2[MAXCAN];
    float tof3[MAXCAN];
    float H2dedx1[MAXCAN];
    float H2dedx2[MAXCAN];
    float H2dedx3[MAXCAN];
    float T4dedx1[MAXCAN];
    float T4dedx2[MAXCAN];
    float T4dedx3[MAXCAN];
    float trkChi1[MAXCAN];
    float trkChi2[MAXCAN];
    float trkChi3[MAXCAN];
   
    //grand-dau info
    float grand_dzos1[MAXCAN];
    float grand_dzos2[MAXCAN];
    float grand_dxyos1[MAXCAN];
    float grand_dxyos2[MAXCAN];
    float grand_nhit1[MAXCAN];
    float grand_nhit2[MAXCAN];
    bool grand_trkquality1[MAXCAN];
    bool grand_trkquality2[MAXCAN];
    float grand_pt1[MAXCAN];
    float grand_pt2[MAXCAN];
    float grand_ptErr1[MAXCAN];
    float grand_ptErr2[MAXCAN];
    float grand_p1[MAXCAN];
    float grand_p2[MAXCAN];
    float grand_eta1[MAXCAN];
    float grand_eta2[MAXCAN];
    int grand_charge1[MAXCAN];
    int grand_charge2[MAXCAN];
    float grand_H2dedx1[MAXCAN];
    float grand_H2dedx2[MAXCAN];
    float grand_T4dedx1[MAXCAN];
    float grand_T4dedx2[MAXCAN];
    float grand_trkChi1[MAXCAN];
    float grand_trkChi2[MAXCAN];

    float grand_dzos21[MAXCAN];
    float grand_dzos22[MAXCAN];
    float grand_dxyos21[MAXCAN];
    float grand_dxyos22[MAXCAN];
    float grand_nhit21[MAXCAN];
    float grand_nhit22[MAXCAN];
    bool grand_trkquality21[MAXCAN];
    bool grand_trkquality22[MAXCAN];
    float grand_pt21[MAXCAN];
    float grand_pt22[MAXCAN];
    float grand_ptErr21[MAXCAN];
    float grand_ptErr22[MAXCAN];
    float grand_p21[MAXCAN];
    float grand_p22[MAXCAN];
    float grand_eta21[MAXCAN];
    float grand_eta22[MAXCAN];
    int grand_charge21[MAXCAN];
    int grand_charge22[MAXCAN];
    float grand_H2dedx21[MAXCAN];
    float grand_H2dedx22[MAXCAN];
    float grand_T4dedx21[MAXCAN];
    float grand_T4dedx22[MAXCAN];
    float grand_trkChi21[MAXCAN];
    float grand_trkChi22[MAXCAN];
    
    //dau muon info
    bool  onestmuon1[MAXCAN];
    bool  onestmuon2[MAXCAN];
    bool  pfmuon1[MAXCAN];
    bool  pfmuon2[MAXCAN];
    bool  glbmuon1[MAXCAN];
    bool  glbmuon2[MAXCAN];
    bool  trkmuon1[MAXCAN];
    bool  trkmuon2[MAXCAN];
    bool  calomuon1[MAXCAN];
    bool  calomuon2[MAXCAN];
    bool  softmuon1[MAXCAN];
    bool  softmuon2[MAXCAN];
    float nmatchedst1[MAXCAN];
    float nmatchedch1[MAXCAN];
    float ntrackerlayer1[MAXCAN];
    float npixellayer1[MAXCAN];
    float matchedenergy1[MAXCAN];
    float nmatchedst2[MAXCAN];
    float nmatchedch2[MAXCAN];
    float ntrackerlayer2[MAXCAN];
    float npixellayer2[MAXCAN];
    float matchedenergy2[MAXCAN];
    float dx1_seg_[MAXCAN];
    float dy1_seg_[MAXCAN];
    float dxSig1_seg_[MAXCAN];
    float dySig1_seg_[MAXCAN];
    float ddxdz1_seg_[MAXCAN];
    float ddydz1_seg_[MAXCAN];
    float ddxdzSig1_seg_[MAXCAN];
    float ddydzSig1_seg_[MAXCAN];
    float dx2_seg_[MAXCAN];
    float dy2_seg_[MAXCAN];
    float dxSig2_seg_[MAXCAN];
    float dySig2_seg_[MAXCAN];
    float ddxdz2_seg_[MAXCAN];
    float ddydz2_seg_[MAXCAN];
    float ddxdzSig2_seg_[MAXCAN];
    float ddydzSig2_seg_[MAXCAN];

    // gen info    
    int candSize_gen;
    int idself[MAXCAN];
    float mass_gen[MAXCAN];
    float pt_gen[MAXCAN];
    float eta_gen[MAXCAN];
    float phi_gen[MAXCAN];
    int status_gen[MAXCAN];
    int idmom[MAXCAN];
    float ptmom[MAXCAN];
    float ymom[MAXCAN];
    float etamom[MAXCAN];
    float phimom[MAXCAN];
    int statusmom[MAXCAN];
    float y_gen[MAXCAN];
    int iddau1[MAXCAN];
    int iddau2[MAXCAN];
    int iddau3[MAXCAN];

    int idself1[MAXCAN];
    float mass_gen1[MAXCAN];
    float pt_gen1[MAXCAN];
    float eta_gen1[MAXCAN];
    float phi_gen1[MAXCAN];
    int status_gen1[MAXCAN];

    int idself2[MAXCAN];
    float mass_gen2[MAXCAN];
    float pt_gen2[MAXCAN];
    float eta_gen2[MAXCAN];
    float phi_gen2[MAXCAN];
    int status_gen2[MAXCAN];

    //vector for gen match
    vector< vector<double> > *pVect;

    vector< vector<double> > *pVectg1;
    vector< vector<double> > *pVectg2;
    vector<double> *Dvector1;
    vector<double> *Dvector2;
    vector<double> *Dvector3;
    
    vector<double> *D1gvector1;
    vector<double> *D1gvector2;
    vector<double> *D2gvector1;
    vector<double> *D2gvector2;
    vector<int> *pVectIDmom;
    
    bool useAnyMVA_;
    bool useDCA_;
    bool isSkimMVA_;
    bool isCentrality_;
    bool isEventPlane_;

    edm::Handle<int> cbin_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token2_;

    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

    edm::EDGetTokenT<int> tok_centBinLabel_;
    edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

    edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

    // for DCA
    edm::EDGetTokenT<std::vector<float > > tok_DCAVal_;
    edm::EDGetTokenT<std::vector<float > > tok_DCAErr_;
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

VertexCompositeTreeProducerNew::VertexCompositeTreeProducerNew(const edm::ParameterSet& iConfig)
{
    //options
    doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
    doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
    twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
    threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
    doubleCand_ = iConfig.getUntrackedParameter<bool>("doubleCand");
    doGenDoubleDecay_ = iConfig.getUntrackedParameter<bool>("doGenDoubleDecay");
    doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
    doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
    hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
    decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
    doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
    doMuonFull_ = iConfig.getUntrackedParameter<bool>("doMuonFull");
    debug_ = iConfig.getUntrackedParameter<bool>("debug");
    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    if(threeProngDecay_) PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");
    if(doGenDoubleDecay_){
      PID_dau1_grand1_ = iConfig.getUntrackedParameter<int>("PID_dau1_grand1") ;
      PID_dau1_grand2_ = iConfig.getUntrackedParameter<int>("PID_dau1_grand2") ;
      PID_dau2_grand1_ = iConfig.getUntrackedParameter<int>("PID_dau2_grand1") ;
      PID_dau2_grand2_ = iConfig.getUntrackedParameter<int>("PID_dau2_grand2") ;
    }
    
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
    MVAValues_Token2_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection2"));
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
    if(useAnyMVA_ && iConfig.exists("MVACollection2"))
      MVAValues_Token2_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection2"));
    if(iConfig.exists("DCAValCollection") && iConfig.exists("DCAErrCollection")) {
      useDCA_ = true;
      tok_DCAVal_ = consumes<std::vector<float > >(iConfig.getParameter<edm::InputTag>("DCAValCollection"));
      tok_DCAErr_ = consumes<std::vector<float > >(iConfig.getParameter<edm::InputTag>("DCAErrCollection"));
    }
}


VertexCompositeTreeProducerNew::~VertexCompositeTreeProducerNew()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeTreeProducerNew::analyze(const edm::Event& iEvent, const edm::EventSetup&
iSetup)
{
    using std::vector;
    using namespace edm;
    using namespace reco;

    if(doGenNtuple_) fillGEN(iEvent,iSetup);
    if(doRecoNtuple_) fillRECO(iEvent,iSetup);

    if(saveTree_) VertexCompositeNtuple->Fill();
}

void
VertexCompositeTreeProducerNew::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    const reco::VertexCompositeCandidateCollection * v0candidates_ = v0candidates.product();
    
    edm::Handle<MVACollection> mvavalues;
    edm::Handle<MVACollection> mvavalues2;
    if(useAnyMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
      if( doubleCand_ ){
        iEvent.getByToken(MVAValues_Token2_,mvavalues2);
        assert( (*mvavalues2).size() == v0candidates->size() );
      }
    }
    edm::Handle<std::vector<float > > dcaValues;
    edm::Handle<std::vector<float > > dcaErrors;
    if(useDCA_)
    {
      iEvent.getByToken(tok_DCAVal_, dcaValues);
      iEvent.getByToken(tok_DCAErr_ , dcaErrors);
      assert( (*dcaValues).size() == v0candidates->size() );
      assert( (*dcaErrors).size() == v0candidates->size() );
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

      HFsumET = cent->EtHFtowerSum();
      Npixel = cent->multiplicityPixel();
//      int ntrk = cent->Ntracks();
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

    //Gen info for matching
    std::vector<reco::GenParticleRef> genRefs;
    if(doGenMatching_)
    {
        pVect = new vector< vector<double>>;
        pVectIDmom = new vector<int>;
        if( doGenDoubleDecay_ ){
          pVectg1 = new vector< vector<double>>;
          pVectg2 = new vector< vector<double>>;
        }
        if(!genpars.isValid()) cout<<"Gen matching cannot be done without Gen collection!!"<<endl; return;
        int count = 0;
        double cache_pt_sameMotherCheck = 0.00;
        double cache_eta_sameMotherCheck = 0.00;
        double cache_phi_sameMotherCheck = 0.00;
        bool cache_sameMotherCheck = false;
        for(unsigned it=0; it<genpars->size(); ++it){
          const reco::GenParticle & trk = (*genpars)[it];
          int id = trk.pdgId();
          int idmom_tmp = -77;
          const reco::Candidate * Dd1 = trk.daughter(0);
          const reco::Candidate * Dd2 = trk.daughter(1);
          if( debug_ && ( Dd1 != nullptr && Dd2 != nullptr)){
          	if( abs(Dd1->pdgId()) == 421 && abs(Dd2->pdgId()) == 421) {
              // std::cout << " Gen dau PDG ID : " << Dd1->pdgId() << ", " << Dd2->pdgId()<< std::endl;
              // std::cout << " Gen trk PDG ID : " << id << ", status " << trk.status() << std::endl;
              if(trk.numberOfMothers()!=0)
              {
                const reco::Candidate * mom = trk.mother();
                idmom_tmp = mom->pdgId();
                if( debug_ ) std::cout << " Gen Mother (nMom " << trk.numberOfMothers() << ") trk PDG ID : " << idmom_tmp << std::endl;
          		  idmom_tmp = -77;
              }
            }
          } 
          //if( fabs(id)!=PID_) continue; //check is target
          if( !doGenDoubleDecay_ &&fabs(id)!=PID_) continue; //check is target
          if( Dd1 == nullptr || Dd2 == nullptr ) continue; //check is target
          if( doGenDoubleDecay_ && !(abs(Dd1->pdgId()) == 421 && abs(Dd2->pdgId()) == 421) ) continue; //check is target
          if(decayInGen_ && trk.numberOfDaughters()!=2 && !threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
          if(decayInGen_ && trk.numberOfDaughters()!=3 && threeProngDecay_) continue; //check 2-pron decay if target decays in Gen
          if(trk.numberOfMothers()!=0)
          {
              const reco::Candidate * mom = trk.mother();
              idmom_tmp = mom->pdgId();
          }
          
          const reco::Candidate * Dd3 = 0;            
          const reco::Candidate * Dd1g1 = 0;
          const reco::Candidate * Dd1g2 = 0;
          const reco::Candidate * Dd2g1 = 0;
          const reco::Candidate * Dd2g2 = 0;


          if(!threeProngDecay_ && !(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau2_) && !(fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd1->pdgId())==PID_dau2_)) continue; //check daughter id                


          if( abs(id) == 4){
            if( cache_sameMotherCheck == false ){
              cache_pt_sameMotherCheck = Dd1->pt();
              cache_eta_sameMotherCheck = Dd1->eta();
              cache_phi_sameMotherCheck = Dd1->phi();
            }
            if( cache_sameMotherCheck == true ){
              if ( cache_pt_sameMotherCheck == Dd1->pt()  && cache_eta_sameMotherCheck == Dd1->eta()  && cache_phi_sameMotherCheck == Dd1->phi() ){
                cache_sameMotherCheck = false;
                continue;
              }
            }
          }

          if(threeProngDecay_)
          {
            Dd3 = trk.daughter(2);
            if(!(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau2_ && fabs(Dd3->pdgId())==PID_dau3_)
            && !(fabs(Dd1->pdgId())==PID_dau1_ && fabs(Dd2->pdgId())==PID_dau3_ && fabs(Dd3->pdgId())==PID_dau2_)
            && !(fabs(Dd1->pdgId())==PID_dau2_ && fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd3->pdgId())==PID_dau3_)
            && !(fabs(Dd1->pdgId())==PID_dau2_ && fabs(Dd2->pdgId())==PID_dau3_ && fabs(Dd3->pdgId())==PID_dau1_)
            && !(fabs(Dd1->pdgId())==PID_dau3_ && fabs(Dd2->pdgId())==PID_dau1_ && fabs(Dd3->pdgId())==PID_dau2_)
            && !(fabs(Dd1->pdgId())==PID_dau3_ && fabs(Dd2->pdgId())==PID_dau2_ && fabs(Dd3->pdgId())==PID_dau1_) ) continue;
          }

          if( doGenDoubleDecay_){ 
            Dd1g1 = Dd1->daughter(0);
            Dd1g2 = Dd1->daughter(1);
            Dd2g1 = Dd2->daughter(0);
            Dd2g2 = Dd2->daughter(1);
            if( ! (
              ((fabs(Dd1g1->pdgId())==PID_dau1_grand1_) && (fabs(Dd1g2->pdgId())==PID_dau1_grand2_)) ||
              ((fabs(Dd1g1->pdgId())==PID_dau1_grand2_) && (fabs(Dd1g2->pdgId())==PID_dau1_grand1_)) ||
              ((fabs(Dd2g1->pdgId())==PID_dau2_grand1_) && (fabs(Dd2g2->pdgId())==PID_dau1_grand2_)) ||
              ((fabs(Dd2g1->pdgId())==PID_dau2_grand2_) && (fabs(Dd2g2->pdgId())==PID_dau1_grand1_))
            ) ) continue;
            if( debug_ ) std::cout << " PDG ID pass " << std::endl;
          }

            Dvector1 = new vector<double>;
            Dvector2 = new vector<double>;


            if( doGenDoubleDecay_){ 
              D1gvector1 = new vector<double>;
              D1gvector2 = new vector<double>;
              D2gvector1 = new vector<double>;
              D2gvector2 = new vector<double>;

              D1gvector1->push_back(Dd1g1->pt());
              D1gvector1->push_back(Dd1g1->eta());
              D1gvector1->push_back(Dd1g1->phi());
              D1gvector1->push_back(Dd1g1->charge());
              D1gvector1->push_back(Dd1g1->mass());

              D1gvector2->push_back(Dd1g2->pt());
              D1gvector2->push_back(Dd1g2->eta());
              D1gvector2->push_back(Dd1g2->phi());
              D1gvector2->push_back(Dd1g2->charge());
              D1gvector2->push_back(Dd1g2->mass());

              D2gvector1->push_back(Dd2g1->pt());
              D2gvector1->push_back(Dd2g1->eta());
              D2gvector1->push_back(Dd2g1->phi());
              D2gvector1->push_back(Dd2g1->charge());
              D2gvector1->push_back(Dd2g1->mass());

              D2gvector2->push_back(Dd2g2->pt());
              D2gvector2->push_back(Dd2g2->eta());
              D2gvector2->push_back(Dd2g2->phi());
              D2gvector2->push_back(Dd2g2->charge());
              D2gvector2->push_back(Dd2g2->mass());

              pVectg1->push_back(*D1gvector1);
              pVectg1->push_back(*D2gvector1);
              pVectg2->push_back(*D1gvector2);
              pVectg2->push_back(*D2gvector2);

              delete D1gvector1;
              delete D1gvector2;
              delete D2gvector1;
              delete D2gvector2;
            }
            if( debug_ ){            
            std::cout << "Dau pT :  "<< Dd1->pt() << ", " << Dd2->pt() << std::endl;
            std::cout << "Dau eta :  "<< Dd1->eta() << ", " << Dd2->eta() << std::endl;
            std::cout << "Dau phi :  "<< Dd1->phi() << ", " << Dd2->phi() << std::endl;
            std::cout << "Dau chg :  "<< Dd1->charge() << ", " << Dd2->charge() << std::endl;
            std::cout << "Dau mass :  "<< Dd1->mass() << ", " << Dd2->mass() << std::endl;
            }
            Dvector1->push_back(Dd1->pt());
            Dvector1->push_back(Dd1->eta());
            Dvector1->push_back(Dd1->phi());
            Dvector1->push_back(Dd1->charge());
            Dvector1->push_back(Dd1->mass());
            
            Dvector2->push_back(Dd2->pt());
            Dvector2->push_back(Dd2->eta());
            Dvector2->push_back(Dd2->phi());
            Dvector2->push_back(Dd2->charge());
            Dvector2->push_back(Dd2->mass());
            
            pVect->push_back(*Dvector1);
            pVect->push_back(*Dvector2);
            
            pVectIDmom->push_back(idmom_tmp);
            
            delete Dvector1;
            delete Dvector2;
count++;

            if(threeProngDecay_)
            {
              Dvector3 = new vector<double>;

              Dvector3->push_back(Dd3->pt());
              Dvector3->push_back(Dd3->eta());
              Dvector3->push_back(Dd3->phi());
              Dvector3->push_back(Dd3->charge());
              Dvector3->push_back(Dd3->mass());

              pVect->push_back(*Dvector3);
              delete Dvector3;
            }
            genRefs.push_back(reco::GenParticleRef(genpars, it));
        }

      if(debug_ && count > 0 ) std::cout << "Filled gen : " << count << " vs. " << genRefs.size() << std::endl;
    }


    //RECO Candidate info
    candSize = v0candidates_->size();
    for(unsigned it=0; it<v0candidates_->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();

        eta[it] = trk.eta();
        y[it] = trk.rapidity();
        pt[it] = trk.pt();
        phi[it] = trk.phi();
        flavor[it] = trk.pdgId()/abs(trk.pdgId());

        mva[it] = 0.0;
        if(useAnyMVA_) mva[it] = (*mvavalues)[it];

        dca3D[it] = -1.0;
        dcaErr3D[it] = -1.0;
        if(useDCA_) {
          dca3D[it] = dcaValues->at(it);
          dcaErr3D[it] = dcaErrors->at(it);
        }

        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass[it] = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        if(doubleCand_ ){

          flavor1[it] = d1->pdgId()/abs(d1->pdgId());
          flavor2[it] = d2->pdgId()/abs(d2->pdgId());
if( debug_ && d1->pt()== d2->pt()  ) std::cout << "Two daughter is same" << std::endl; 
        }
        const reco::Candidate * d3 = 0;        
        if(threeProngDecay_) d3 = trk.daughter(2);

        //Gen match
        if(doGenMatching_)
        {
            if( !doGenDoubleDecay_ ){
              matchGEN[it] = false;
              int nGenDau = (int)pVect->size();
              isSwap[it] = false;
              idmom_reco[it] = -77;
            
              for(int i=0;i<nGenDau;i++)
              {
                  vector<double> Dvector1_ = (*pVect)[i]; //get GEN daughter vector
                  if(d1->charge()!=Dvector1_.at(3)) continue; //check match charge
                  double deltaR = sqrt(pow(d1->eta()-Dvector1_.at(1),2)+pow(d1->phi()-Dvector1_.at(2),2));

                  if(deltaR > deltaR_) continue; //check deltaR matching
                  if(fabs((d1->pt()-Dvector1_.at(0))/d1->pt()) > 0.5) continue; //check deltaPt matching
                  double d1massGEN = Dvector1_.at(4);
                  double d1mass = d1->mass();
                  double d2massGEN=0, d2mass=0;
                  double d3massGEN=0, d3mass=0;

                  if(nGenDau==2)
                  {
                    if(i%2==0)
                    {
                      vector<double> Dvector2 = (*pVect)[i+1]; //get GEN daughter vector for track2
                      if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                      double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));

                      if(deltaR > deltaR_) continue; //check deltaR matching
                      if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                      d2massGEN = Dvector2.at(4);
                      d2mass = d2->mass();

                      matchGEN[it] = true; //matched gen
                    }

                    if(i%2==1)
                    {
                      vector<double> Dvector2 = (*pVect)[i-1]; //get GEN daughter vector for track2
                      if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                      double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));

                      if(deltaR > deltaR_) continue; //check deltaR matching
                      if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                      d2massGEN = Dvector2.at(4);
                      d2mass = d2->mass();

                      matchGEN[it] = true; //matched gen
                    }

                    //check swap
                    if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap[it] = true;

                    //check prompt & record mom id
                    idmom_reco[it] = pVectIDmom->at(i/2);
                  }

                  if(nGenDau==3)
                  {
                    if(i%3==0)
                    {
                      vector<double> Dvector2 = (*pVect)[i+1]; //get GEN daughter vector for track2
                      vector<double> Dvector3 = (*pVect)[i+2]; //get GEN daughter vector for track3

                      if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3)) 
                      && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                      double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                      double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                      double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                      double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                      if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                      double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                      double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                      double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                      double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                      if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                      d2massGEN = Dvector2.at(4);
                      d2mass = d2->mass();
                      d3massGEN = Dvector3.at(4);
                      d3mass = d3->mass();

                      matchGEN[it] = true; //matched gen
                    }

                    if(i%3==1)
                    {
                      vector<double> Dvector2 = (*pVect)[i-1]; //get GEN daughter vector for track2
                      vector<double> Dvector3 = (*pVect)[i+1]; //get GEN daughter vector for track3

                      if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3))
                      && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                      double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                      double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                      double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                      double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                      if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                      double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                      double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                      double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                      double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                      if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                      d2massGEN = Dvector2.at(4);
                      d2mass = d2->mass();
                      d3massGEN = Dvector3.at(4);
                      d3mass = d3->mass();

                      matchGEN[it] = true; //matched gen
                    }

                    if(i%3==2)
                    { 
                      vector<double> Dvector2 = (*pVect)[i-2]; //get GEN daughter vector for track2
                      vector<double> Dvector3 = (*pVect)[i-1]; //get GEN daughter vector for track3

                      if(!(d2->charge()==Dvector2.at(3) && d3->charge()==Dvector3.at(3))
                      && !(d3->charge()==Dvector2.at(3) && d2->charge()==Dvector3.at(3))) continue; //check match charge

                      double deltaR22 = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                      double deltaR33 = sqrt(pow(d3->eta()-Dvector3.at(1),2)+pow(d3->phi()-Dvector3.at(2),2));
                      double deltaR23 = sqrt(pow(d2->eta()-Dvector3.at(1),2)+pow(d2->phi()-Dvector3.at(2),2));
                      double deltaR32 = sqrt(pow(d3->eta()-Dvector2.at(1),2)+pow(d3->phi()-Dvector2.at(2),2));

                      if(!(deltaR22 < deltaR_ && deltaR33 < deltaR_) && !(deltaR23 < deltaR_ && deltaR32 < deltaR_) ) continue;

                      double deltaPt22 = fabs((d2->pt()-Dvector2.at(0))/d2->pt());
                      double deltaPt33 = fabs((d3->pt()-Dvector3.at(0))/d3->pt());
                      double deltaPt23 = fabs((d2->pt()-Dvector3.at(0))/d2->pt());
                      double deltaPt32 = fabs((d3->pt()-Dvector2.at(0))/d3->pt());

                      if( !(deltaPt22 < 0.5 && deltaPt33 < 0.5) && !(deltaPt23 < 0.5 && deltaPt32 < 0.5) ) continue; //check deltaPt matching

                      d2massGEN = Dvector2.at(4);
                      d2mass = d2->mass();
                      d3massGEN = Dvector3.at(4);
                      d3mass = d3->mass();

                      matchGEN[it] = true; //matched gen
                    }

                    //check swap
                    if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01 || abs(d3massGEN - d3mass)>0.01) isSwap[it] = true;

                    //check prompt & record mom id
                    idmom_reco[it] = pVectIDmom->at(i/3);
                  }

              }
            }
            if( doGenDoubleDecay_ ){
              matchGEN1[it] = false;
              matchGEN2[it] = false;
              isSwap1[it] = false;
              isSwap2[it] = false;
              idmom_reco[it] = -77;
              const auto nGen = genRefs.size();
              int matchGenIdx = -1;
              std::pair<unsigned int, unsigned int> idxMatch1 = {99999, 999};
              std::pair<unsigned int, unsigned int> idxMatch2 = {99999, 999};
              for( unsigned int igen=0; igen<nGen; igen++){
                auto const& theGen = genRefs.at(igen);
                // Only works for 2 body two layer decay
                reco::Candidate const* recoDaus[2] = {nullptr, nullptr};
                reco::Candidate const* recoGDaus1[2] = {nullptr, nullptr};
                reco::Candidate const* recoGDaus2[2] = {nullptr, nullptr};

                reco::Candidate const* genDaus[2] = {nullptr, nullptr};
                reco::Candidate const* genGDaus1[2] = {nullptr, nullptr};
                reco::Candidate const* genGDaus2[2] = {nullptr, nullptr};

                const auto nGenDau = theGen->numberOfDaughters();
                const auto nGenGDau1 = theGen->daughter(0)->numberOfDaughters();
                const auto nGenGDau2 = theGen->daughter(1)->numberOfDaughters();
                std::vector<unsigned int> perm = {0, 1};
                std::vector<unsigned int> perm1 = {0, 1};
                std::vector<unsigned int> perm2 = {0, 1};
                do{
                  matchGEN[it] = false;
                  for( unsigned int iDau=0; iDau<nGenDau; ++iDau){
                    genDaus[iDau] = theGen->daughter(perm.at(iDau));
                    recoDaus[iDau] = trk.daughter(iDau);
                  }
                  do{
                    matchGEN1[it] = false;
                    for( unsigned int iGDau1=0; iGDau1<nGenGDau1; ++iGDau1){
                      genGDaus1[iGDau1] = genDaus[0]->daughter(perm1.at(iGDau1));
                      recoGDaus1[iGDau1] = recoDaus[0]->daughter(iGDau1);
                    }
                    for( unsigned int iGDau1=0; iGDau1<nGenGDau1; ++iGDau1){
                      const double dR = reco::deltaR(genGDaus1[iGDau1]->eta(), genGDaus1[iGdau1]->phi(),
                          recoGDaus1[iGdau1]->eta(), regoGDaus1[iGdau1]->phi());
                      const double dPt = abs(genGDaus1[iGdau1]->pt()-regoGDaus1[iGdau1]->pt())/regoGDaus1[iGdau1]->pt();
                      const bool unMatchCharge = genGDaus1[iGdau1]->charge() != regoGDaus1[iGdau1]->charge();
                      const bool unMatchDR = dR > deltaR_;
                      const bool unMatchDPt = dPt > 0.5;
                      matchGEN1[it] = !(unMatchCharge || unMatchDR || unMatchDPt);
                      if(matchGEN1[it]) idxMatch1 = {perm.at(0), iGDau1};
                    }
                    if(matchGEN1[it]) break;
                  } while( std::next_permutation(perm1.begin(), perm1.end()));
                  do{
                    matchGEN2[it] = false;
                    for( unsigned int iGDau2=0; iGDau2<nGenGDau2; ++iGDau2){
                      genGDaus2[iGDau2] = genDaus[1]->daughter(perm2.at(iGDau2));
                      recoGDaus2[iGDau2] = recoDaus[1]->daughter(iGDau2);
                    }
                    for( unsigned int iGDau2=0; iGDau2<nGenGDau2; ++iGDau2){
                      const double dR = reco::deltaR(genGDaus2[iGdau2]->eta(), genGDaus2[iGdau2]->phi(),
                          recoGDaus2[iGdau2]->eta(), regoGDaus2[iGdau2]->phi());
                      const double dPt = abs(genGDaus2[iGdau2]->pt()-regoGDaus2[iGdau2]->pt())/regoGDaus2[iGdau2]->pt();
                      const bool unMatchCharge = genGDaus2[iGdau2]->charge() != regoGDaus2[iGdau2]->charge();
                      const bool unMatchDR = dR > deltaR_;
                      const bool unMatchDPt = dPt > 0.5;
                      matchGEN2[it] = !(unMatchCharge || unMatchDR || unMatchDPt);
                      if(matchGEN2[it]) idxMatch2 = {perm.at(1), iGDau2};
                    }
                    if(matchGEN2[it]) break;
                  } while( std::next_permutation(perm2.begin(), perm2.end()));
                  if( matchGEN1[it] &&matchGEN2[it]){
                    matchGEN[it] = (idxMatch1.first != idxMatch2.first && idxMatch1.second != idxMatch2.second) ? true : false;
                    if(matchGEN[it]) matchGenIdx = igen;
                  }
                  std::sort(perm1.begin(), perm1.end());
                  std::sort(perm2.begin(), perm2.end());
                } while( std::next_permutation(perm.begin(), perm.end()));
                if( matchGEN[it] ) break;

              } // END for nGen
              if( matchGEN1[it] && matchGEN){
                auto const& theGen = genRefs.at(matchGenIdx);
                auto idgenDau1 = theGen->daughter(idxMatch1.first)->pdgId();
                auto idgenDau2 = theGen->daughter(idxMatch2.first)->pdgId();
                auto idrecoDau1 = trk.daughter(idxMatch1.first)->pdgId();
                auto idrecoDau2 = trk.daughter(idxMatch2.first)->pdgId();
                isSwap1[it]  = (idgenDau1 == idrecoDau1);
                isSwap2[it]  = (idgenDau2 == idrecoDau2);
              }
            } // END if doGenDoubleDecay_

//             if( doGenDoubleDecay_ ){
//               matchGEN[it] = false;
//               int nGenDau = (int)pVect->size();
// if( debug_ )	std::cout << " Checking Match... (nDau) " << nGenDau << std::endl;
//               isSwap[it] = false;
//               idmom_reco[it] = -77;
            
//               for(int i=0;i<nGenDau;i++)
//               {

//               	  matchGEN1[it] = false;
//               	  matchGEN2[it] = false;
			
//                   vector<double> Dvector1_ = (*pVect)[i]; //get GEN daughter vector, a D meson 
//                   //if(d1->charge()!=Dvector1_.at(3)) continue; //check match charge
//                   //double deltaR = sqrt(pow(d1->eta()-Dvector1_.at(1),2)+pow(d1->phi()-Dvector1_.at(2),2));

//                   //if(deltaR > deltaR_) continue; //check deltaR matching
//                   //if(fabs((d1->pt()-Dvector1_.at(0))/d1->pt()) > 0.5) continue; //check deltaPt matching
//                   double d1massGEN = Dvector1_.at(4);
//                   double d1mass = d1->mass();
//                   double d2massGEN=0, d2mass=0;

//                   double d1gmassGEN2 = 0;
//                   double d1gmass2 = 0;
//                   double d1gmassGEN1 = 0;
//                   double d1gmass1 = 0;
//                   double d2gmassGEN2 = 0;
//                   double d2gmass2 = 0;
//                   double d2gmassGEN1 = 0;
//                   double d2gmass1 = 0;
//                   bool matchChargeGD1 = false;
//                   bool matchChargeGD2 = false;
//                   bool matchDRGD1 = false;
//                   bool matchDRGD2 = false;
//                   bool matchDPTGD1 = false;
//                   bool matchDPTGD2 = false;

//                   for( int ii=0;ii<2;ii++){
//                     matchChargeGD1=false; matchDRGD1=false;matchDPTGD1=false;
//                     matchChargeGD2=false; matchDRGD2=false;matchDPTGD2=false;
//                     const reco::Candidate * gd11 = d1->daughter(ii);
//                     const reco::Candidate * gd12 = d1->daughter(1-ii);
//                     vector<double>& D1gvector1_ = (*pVectg1)[i]; //get GEN grand daughter vector, a D meson's track 1
//                     // double deltaRg1 = sqrt(pow(gd11->eta()-D1gvector1_.at(1),2)+pow(gd11->phi()-D1gvector1_.at(2),2));
//                     matchChargeGD1 = gd11->charge()==D1gvector1_.at(3);
//                     matchDRGD1 = sqrt(pow(gd11->eta()-D1gvector1_.at(1),2)+pow(gd11->phi()-D1gvector1_.at(2),2)) <= deltaR_;
//                     matchDPTGD1 = fabs((gd11->pt()-D1gvector1_.at(0))/gd11->pt()) <= 0.5;

	
// if( debug_ )	std::cout   << " D1 g1 : charge ("<< gd11->charge() << ", "<< D1gvector1_.at(3) << "), dR ("<<sqrt(pow(gd11->eta()-D1gvector1_.at(1),2)+pow(gd11->phi()-D1gvector1_.at(2),2)) << ", " << deltaR_ <<"), dPt ("<< fabs((gd11->pt()-D1gvector1_.at(0))/gd11->pt()) <<","<< "0.5" <<") : " << matchChargeGD1 << ", "<< matchDRGD1 <<", "<< matchDPTGD1 << std::endl;

//                     vector<double> & D1gvector2_ = (*pVectg2)[i]; //get GEN grand daughter vector, a D meson's track 1
//                     // double deltaRg2 = sqrt(pow(gd12->eta()-D1gvector2_.at(1),2)+pow(gd12->phi()-D1gvector2_.at(2),2));
//                     matchChargeGD2 = gd12->charge()==D1gvector2_.at(3);
//                     matchDRGD2 = sqrt(pow(gd12->eta()-D1gvector2_.at(1),2)+pow(gd12->phi()-D1gvector2_.at(2),2)) <= deltaR_;
//                     matchDPTGD2 = fabs((gd12->pt()-D1gvector2_.at(0))/gd12->pt()) <= 0.5;
// if( debug_ )	std::cout   << " D1 g2 : charge ("<< gd12->charge() << ", "<< D1gvector2_.at(3) << "), dR ("<<sqrt(pow(gd12->eta()-D1gvector2_.at(1),2)+pow(gd12->phi()-D1gvector2_.at(2),2)) << ", " << deltaR_ <<"), dPt ("<< fabs((gd12->pt()-D1gvector2_.at(0))/gd12->pt()) <<","<< "0.5" <<") : " << matchChargeGD2 << ", "<< matchDRGD2 <<", "<< matchDPTGD2 << std::endl;

//                     matchGEN1[it] = matchGEN1[it] || (matchChargeGD1 && matchChargeGD2 && matchDRGD1 && matchDRGD2 && matchDPTGD1 && matchDPTGD2 );

// if( debug_ ) std::cout << "match : " <<  matchGEN1[it]  << std::endl;
//                     if( matchGEN1[it] ) {
//                       d1gmassGEN1 = D1gvector1_.at(4);
//                       d1gmass1 = gd11->mass();
//                       d1gmassGEN2 = D1gvector2_.at(4);
//                       d1gmass2 = gd12->mass();
//                     }
//                   }

//                   if((nGenDau%2) == 0)
//                   {
//                     int i2 = (i%2==0) ? i+1 : i-1;
//                     vector<double> Dvector2 = (*pVect)[i2]; //get GEN daughter vector for D meson 2
//                     //if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
//                     //double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));

//                     //if(deltaR > deltaR_) continue; //check deltaR matching
//                     //if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
//                     d2massGEN = Dvector2.at(4);
//                     d2mass = d2->mass();

//                     // deltaRg1 =999.9; deltaRg2 = 999.9;
//                   for( int ii=0;ii<2;ii++){
//                     matchChargeGD1=false; matchDRGD1=false;matchDPTGD1=false;
//                     matchChargeGD2=false; matchDRGD2=false;matchDPTGD2=false;
//                     const reco::Candidate * gd21 = d2->daughter(ii);
//                     const reco::Candidate * gd22 = d2->daughter(1-ii);
//                     vector<double> D2gvector1_ = (*pVectg1)[i2]; //get GEN grand daughter vector, a D meson's track 1
//                     // deltaRg1 = sqrt(pow(gd21->eta()-D2gvector1_.at(1),2)+pow(gd21->phi()-D2gvector1_.at(2),2));
//                     matchChargeGD1 = gd21->charge()==D2gvector1_.at(3);
//                     matchDRGD1 = sqrt(pow(gd21->eta()-D2gvector1_.at(1),2)+pow(gd21->phi()-D2gvector1_.at(2),2)) <= deltaR_;
//                     matchDPTGD1 = fabs((gd21->pt()-D2gvector1_.at(0))/gd21->pt()) <= 0.5;

// if( debug_ )	std::cout << " D2 g1 : charge ("<< gd21->charge() << ", "<< D2gvector1_.at(3) << "), dR ("<<sqrt(pow(gd21->eta()-D2gvector1_.at(1),2)+pow(gd21->phi()-D2gvector1_.at(2),2)) << ", " << deltaR_ <<"), dPt ("<< fabs((gd21->pt()-D2gvector1_.at(0))/gd21->pt()) <<","<< "0.5" <<") : " << matchChargeGD1 << ", "<< matchDRGD1 <<", "<< matchDPTGD1 << std::endl;

//                     vector<double> D2gvector2_ = (*pVectg2)[i2]; //get GEN grand daughter vector, a D meson's track 1
//                     // deltaRg2 = sqrt(pow(gd22->eta()-D2gvector2_.at(1),2)+pow(gd22->phi()-D2gvector2_.at(2),2));
//                     matchChargeGD2 = gd22->charge()==D2gvector2_.at(3);
//                     matchDRGD2 = sqrt(pow(gd22->eta()-D2gvector2_.at(1),2)+pow(gd22->phi()-D2gvector2_.at(2),2)) <= deltaR_;
//                     matchDPTGD2 = fabs((gd22->pt()-D2gvector2_.at(0))/gd22->pt()) <= 0.5;

// if( debug_ )	std::cout << " D2 g2 : charge ("<< gd22->charge() << ", "<< D2gvector2_.at(3) << "), dR ("<<sqrt(pow(gd22->eta()-D2gvector2_.at(1),2)+pow(gd22->phi()-D2gvector2_.at(2),2)) << ", " << deltaR_ <<"), dPt ("<< fabs((gd22->pt()-D2gvector2_.at(0))/gd22->pt()) <<","<< "0.5" <<") : " << matchChargeGD2 << ", "<< matchDRGD2 <<", "<< matchDPTGD2 << std::endl;

//                     matchGEN2[it] = matchGEN2[it] || (matchChargeGD1 && matchChargeGD2 && matchDRGD1 && matchDRGD2 && matchDPTGD1 && matchDPTGD2 );
//                     if( matchGEN2[it] ) {
//                       d2gmassGEN1 = D2gvector1_.at(4);
//                       d2gmass1 = gd21->mass();
//                       d2gmassGEN2 = D2gvector2_.at(4);
//                       d2gmass2 = gd22->mass();
//                     }
//                   }

//                     matchGEN[it] = matchGEN1[it] && matchGEN2[it]; //matched gen

//                     //check swap
//                     if(abs(d1gmassGEN1 - d1gmass1)>0.01 || abs(d1gmassGEN2 - d1gmass2)>0.01) isSwap1[it] = true;
//                     if(abs(d2gmassGEN1 - d2gmass1)>0.01 || abs(d2gmassGEN2 - d2gmass2)>0.01) isSwap2[it] = true;
//                     if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap[it] = true;

//                     //check prompt & record mom id
//                     idmom_reco[it] = pVectIDmom->at(i/2);
//                   }
//               }  // end GEN dau loop
//             } // end double decay 
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
        pt1[it] = d1->pt();
        pt2[it] = d2->pt();
        
        //momentum
        p1[it] = d1->p();
        p2[it] = d2->p();
        
        //eta
        eta1[it] = d1->eta();
        eta2[it] = d2->eta();
        
        //phi
        phi1[it] = d1->phi();
        phi2[it] = d2->phi();
        
        //charge
        charge1[it] = d1->charge();
        charge2[it] = d2->charge();
        
        double pxd3 = -999.9;
        double pyd3 = -999.9;
        double pzd3 = -999.9;
        if(threeProngDecay_ && d3)
        {
          pxd3 = d3->px();
          pyd3 = d3->py();
          pzd3 = d3->pz();
          pt3[it] = d3->pt();
          p3[it] = d3->p();
          eta3[it] = d3->eta();
          phi3[it] = d3->phi();
          charge3[it] = d3->charge();
        }
        TVector3 dauvec3(pxd3,pyd3,pzd3);

        pid1[it] = -99999;
        pid2[it] = -99999;
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
                if(deltaR < deltaR_ && fabs((trk.pt()-pt1[it])/pt1[it]) < 0.5 && trk.charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id;
                } 

                // matching daughter 2
                deltaR = trkvect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((trk.pt()-pt2[it])/pt2[it]) < 0.5 && trk.charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id;
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
                if(deltaR < deltaR_ && fabs((Dd1->pt()-pt1[it])/pt1[it]) < 0.5 && Dd1->charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id1;
                }
                deltaR = d2vect.DeltaR(dauvec1);
                if(deltaR < deltaR_ && fabs((Dd2->pt()-pt1[it])/pt1[it]) < 0.5 && Dd2->charge()==charge1[it] && pid1[it]==-99999)
                {
                  pid1[it] = id1;
                }

                deltaR = d1vect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((Dd1->pt()-pt2[it])/pt2[it]) < 0.5 && Dd1->charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id2;
                }
                deltaR = d2vect.DeltaR(dauvec2);
                if(deltaR < deltaR_ && fabs((Dd2->pt()-pt2[it])/pt2[it]) < 0.5 && Dd2->charge()==charge2[it] && pid2[it]==-99999)
                {
                  pid2[it] = id2;
                }
              }

              if(pid1[it]!=-99999 && pid2[it]!=-99999) break;
          }
        }

        //vtxChi2
        vtxChi2[it] = trk.vertexChi2();
        ndf[it] = trk.vertexNdof();
        VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
        
        agl[it] = cos(secvec.Angle(ptosvec));
        agl_abs[it] = secvec.Angle(ptosvec);
        
        agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
        
        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl[it] = ROOT::Math::Mag(distanceVector);
        dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl[it];
        
        dlos[it] = dl[it]/dlerror[it];

        // correct way for both DCA and its Error
        // std::cout << "By cur3DIP                " << dca3D[it] << " +/- " << dcaErr3D[it] <<"\n";
        // incorrect way for DCA error
        // std::cout << "By decay length and alpha " << std::sin(agl_abs[it])*dl[it]<< " +/- " << dlerror[it]* std::sin(agl_abs[it]) <<"\n";
        // std::cout << "\n";
        
        //Decay length 2D
        SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
        SVector6 v2(trk.vertexCovariance(0,0), trk.vertexCovariance(0,1),trk.vertexCovariance(1,1),0,0,0);
        
        SMatrixSym3D sv1(v1);
        SMatrixSym3D sv2(v2);
        
        SMatrixSym3D totalCov2D = sv1 + sv2;
        SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
        
        dl2D[it] = ROOT::Math::Mag(distanceVector2D);
        double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];
        
        dlos2D[it] = dl2D[it]/dl2Derror;

        //trk info
        auto dau1 = d1->get<reco::TrackRef>();
        if(!twoLayerDecay_)
        {
            //trk quality
            trkquality1[it] = dau1->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            H2dedx1[it] = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                H2dedx1[it] = dEdxTrack[dau1].dEdx();
            }
            
            T4dedx1[it] = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                T4dedx1[it] = dEdxTrack[dau1].dEdx();
            }
            
            //track Chi2
            trkChi1[it] = dau1->normalizedChi2();
            
            //track pT error
            ptErr1[it] = dau1->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            
            //trkNHits
            nhit1[it] = dau1->numberOfValidHits();
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzbest1 = dau1->dz(bestvtx);
            double dxybest1 = dau1->dxy(bestvtx);
            double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
            double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
            
            dzos1[it] = dzbest1/dzerror1;
            dxyos1[it] = dxybest1/dxyerror1;
        }
        if( !doubleCand_){
        auto dau2 = d2->get<reco::TrackRef>();
        
        //trk quality
        trkquality2[it] = dau2->quality(reco::TrackBase::highPurity);
        
        //trk dEdx
        H2dedx2[it] = -999.9;
        
        if(dEdxHandle1.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
            H2dedx2[it] = dEdxTrack[dau2].dEdx();
        }
        
        T4dedx2[it] = -999.9;
        
        if(dEdxHandle2.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
            T4dedx2[it] = dEdxTrack[dau2].dEdx();
        }
        
        //track Chi2
        trkChi2[it] = dau2->normalizedChi2();
        
        //track pT error
        ptErr2[it] = dau2->ptError();
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        
        //trkNHits
        nhit2[it] = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2[it] = dzbest2/dzerror2;
        dxyos2[it] = dxybest2/dxyerror2;
        
        if(doMuon_)
        {
          edm::Handle<reco::MuonCollection> theMuonHandle;
          iEvent.getByToken(tok_muon_, theMuonHandle);
            
          nmatchedch1[it] = -1;
          nmatchedst1[it] = -1;
          matchedenergy1[it] = -1;
          nmatchedch2[it] = -1;
          nmatchedst2[it] = -1;
          matchedenergy2[it] = -1;
            
          double x_exp = -999.;
          double y_exp = -999.;
          double xerr_exp = -999.;
          double yerr_exp = -999.;
          double dxdz_exp = -999.;
          double dydz_exp = -999.;
          double dxdzerr_exp = -999.;
          double dydzerr_exp = -999.;
            
          double x_seg = -999.;
          double y_seg = -999.;
          double xerr_seg = -999.;
          double yerr_seg = -999.;
          double dxdz_seg = -999.;
          double dydz_seg = -999.;
          double dxdzerr_seg = -999.;
          double dydzerr_seg = -999.;
            
          double dx_seg = 999.;
          double dy_seg = 999.;
          double dxerr_seg = 999.;
          double dyerr_seg = 999.;
          double dxSig_seg = 999.;
          double dySig_seg = 999.;
          double ddxdz_seg = 999.;
          double ddydz_seg = 999.;
          double ddxdzerr_seg = 999.;
          double ddydzerr_seg = 999.;
          double ddxdzSig_seg = 999.;
          double ddydzSig_seg = 999.;
            
          onestmuon1[it] = false;
          pfmuon1[it] = false;
          glbmuon1[it] = false;
          trkmuon1[it] = false;
          calomuon1[it] = false; 
          softmuon1[it] = false;
          onestmuon2[it] = false;
          pfmuon2[it] = false;
          glbmuon2[it] = false;
          trkmuon2[it] = false;
          calomuon2[it] = false;
          softmuon2[it] = false;

          const int muId1 = muAssocToTrack( dau1, theMuonHandle );
          const int muId2 = muAssocToTrack( dau2, theMuonHandle );

          if( muId1 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId1];

            onestmuon1[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
            pfmuon1[it] =  cand.isPFMuon();
            glbmuon1[it] =  cand.isGlobalMuon();
            trkmuon1[it] =  cand.isTrackerMuon();
            calomuon1[it] =  cand.isCaloMuon();

            if( 
                glbmuon1[it] && trkmuon1[it] &&
                cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && 
                cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 && 
                fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
                fabs(cand.innerTrack()->dz(vtx.position())) < 20.
              ) softmuon1[it] = true;
          }

          if( muId2 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId2];

            onestmuon2[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
            pfmuon2[it] =  cand.isPFMuon();
            glbmuon2[it] =  cand.isGlobalMuon();
            trkmuon2[it] =  cand.isTrackerMuon();
            calomuon2[it] =  cand.isCaloMuon();

            if(
                glbmuon2[it] && trkmuon2[it] &&
                cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
                fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
                fabs(cand.innerTrack()->dz(vtx.position())) < 20.
              ) softmuon2[it] = true;
          }

          if(doMuonFull_)
          {

          if( muId1 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId1];

            nmatchedch1[it] = cand.numberOfMatches();
            nmatchedst1[it] = cand.numberOfMatchedStations();
                   
            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy1[it] = muenergy.hadMax;
                    
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
                      
              dx1_seg_[it]=dx_seg;
              dy1_seg_[it]=dy_seg;
              dxSig1_seg_[it]=dxSig_seg;
              dySig1_seg_[it]=dySig_seg;
              ddxdz1_seg_[it]=ddxdz_seg;
              ddydz1_seg_[it]=ddydz_seg;
              ddxdzSig1_seg_[it]=ddxdzSig_seg;
              ddydzSig1_seg_[it]=ddydzSig_seg;
            }
          } 

          if( muId2 != -1 )
          {
            const reco::Muon& cand = (*theMuonHandle)[muId2];

            nmatchedch2[it] = cand.numberOfMatches();
            nmatchedst2[it] = cand.numberOfMatchedStations();
                    
            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy2[it] = muenergy.hadMax;
                    
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
                        
              dx2_seg_[it]=dx_seg;
              dy2_seg_[it]=dy_seg;
              dxSig2_seg_[it]=dxSig_seg;
              dySig2_seg_[it]=dySig_seg;
              ddxdz2_seg_[it]=ddxdz_seg;
              ddydz2_seg_[it]=ddydz_seg;
              ddxdzSig2_seg_[it]=ddxdzSig_seg;
              ddydzSig2_seg_[it]=ddydzSig_seg;
            }
          }
          } // doMuonFull
        }
        }
        
        if(twoLayerDecay_)
        {
            grand_mass[it] = d1->mass();
            mva1[it] = (*mvavalues)[it];
            
            const reco::Candidate * gd1 = d1->daughter(0);
            const reco::Candidate * gd2 = d1->daughter(1);
            
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
            
            grand_trkquality1[it] = gdau1->quality(reco::TrackBase::highPurity);
            grand_trkquality2[it] = gdau2->quality(reco::TrackBase::highPurity);
            
            //trk dEdx
            grand_H2dedx1[it] = -999.9;
            grand_H2dedx2[it] = -999.9;
            
            if(dEdxHandle1.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                grand_H2dedx1[it] = dEdxTrack[gdau1].dEdx();
                grand_H2dedx2[it] = dEdxTrack[gdau2].dEdx();
            }
            
            grand_T4dedx1[it] = -999.9;
            grand_T4dedx2[it] = -999.9;
            
            if(dEdxHandle2.isValid()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                grand_T4dedx1[it] = dEdxTrack[gdau1].dEdx();
                grand_T4dedx2[it] = dEdxTrack[gdau2].dEdx();
            }
            
            //track pt
            grand_pt1[it] = gd1->pt();
            grand_pt2[it] = gd2->pt();
            
            //track momentum
            grand_p1[it] = gd1->p();
            grand_p2[it] = gd2->p();
            
            //track eta
            grand_eta1[it] = gd1->eta();
            grand_eta2[it] = gd2->eta();
            
            //track charge
            grand_charge1[it] = gd1->charge();
            grand_charge2[it] = gd2->charge();
            
            //track Chi2
            grand_trkChi1[it] = gdau1->normalizedChi2();
            grand_trkChi2[it] = gdau2->normalizedChi2();
            
            //track pT error
            grand_ptErr1[it] = gdau1->ptError();
            grand_ptErr2[it] = gdau2->ptError();
            
            //vertexCovariance 00-xError 11-y 22-z
            secvz = d1->vz(); secvx = d1->vx(); secvy = d1->vy();
            
            //trkNHits
            grand_nhit1[it] = gdau1->numberOfValidHits();
            grand_nhit2[it] = gdau2->numberOfValidHits();
            
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double gdzbest1 = gdau1->dz(bestvtx);
            double gdxybest1 = gdau1->dxy(bestvtx);
            double gdzerror1 = sqrt(gdau1->dzError()*gdau1->dzError()+bestvzError*bestvzError);
            double gdxyerror1 = sqrt(gdau1->d0Error()*gdau1->d0Error()+bestvxError*bestvyError);
            
            grand_dzos1[it] = gdzbest1/gdzerror1;
            grand_dxyos1[it] = gdxybest1/gdxyerror1;
            
            double gdzbest2 = gdau2->dz(bestvtx);
            double gdxybest2 = gdau2->dxy(bestvtx);
            double gdzerror2 = sqrt(gdau2->dzError()*gdau2->dzError()+bestvzError*bestvzError);
            double gdxyerror2 = sqrt(gdau2->d0Error()*gdau2->d0Error()+bestvxError*bestvyError);
            
            grand_dzos2[it] = gdzbest2/gdzerror2;
            grand_dxyos2[it] = gdxybest2/gdxyerror2;
            
            //vtxChi2
            grand_vtxChi2[it] = d1->vertexChi2();
            grand_ndf[it] = d1->vertexNdof();
            grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it],grand_ndf[it]);
            
            //PAngle
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(d1->px(),d1->py(),d1->pz());
            
            TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
            TVector3 secvec2D(d1->px(),d1->py(),0);
            
            grand_agl[it] = cos(secvec.Angle(ptosvec));
            grand_agl_abs[it] = secvec.Angle(ptosvec);
            
            grand_agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
            grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
            
            //Decay length 3D
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            typedef ROOT::Math::SVector<double, 6> SVector6;
            
            SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            grand_dl[it] = ROOT::Math::Mag(distanceVector);
            grand_dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl[it];
            
            grand_dlos[it] = grand_dl[it]/grand_dlerror[it];
            
            //Decay length 2D
            SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
            SVector6 v2(d1->vertexCovariance(0,0), d1->vertexCovariance(0,1),d1->vertexCovariance(1,1),0,0,0);
            
            SMatrixSym3D sv1(v1);
            SMatrixSym3D sv2(v2);
            
            SMatrixSym3D totalCov2D = sv1 + sv2;
            SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
            
            double gdl2D = ROOT::Math::Mag(distanceVector2D);
            double gdl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
            
            grand_dlos2D[it] = gdl2D/gdl2Derror;
            if( doubleCand_ )
            {
                grand_mass2[it] = d2->mass();
                mva2[it] = (*mvavalues2)[it];

                const reco::Candidate * gd21 = d2->daughter(0);
                const reco::Candidate * gd22 = d2->daughter(1);

                double gpxd21 = gd21->px();
                double gpyd21 = gd21->py();
                double gpzd21 = gd21->pz();
                double gpxd22 = gd22->px();
                double gpyd22 = gd22->py();
                double gpzd22 = gd22->pz();

                TVector3 gdauvec21(gpxd21,gpyd21,gpzd21);
                TVector3 gdauvec22(gpxd22,gpyd22,gpzd22);

                auto gdau21 = gd21->get<reco::TrackRef>();
                auto gdau22 = gd22->get<reco::TrackRef>();

                //trk quality

                grand_trkquality21[it] = gdau21->quality(reco::TrackBase::highPurity);
                grand_trkquality22[it] = gdau22->quality(reco::TrackBase::highPurity);

                //trk dEdx
                grand_H2dedx21[it] = -999.9;
                grand_H2dedx22[it] = -999.9;

                if(dEdxHandle1.isValid()){
                    const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                    grand_H2dedx21[it] = dEdxTrack[gdau21].dEdx();
                    grand_H2dedx22[it] = dEdxTrack[gdau22].dEdx();
                }

                grand_T4dedx21[it] = -999.9;
                grand_T4dedx22[it] = -999.9;

                if(dEdxHandle2.isValid()){
                    const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                    grand_T4dedx21[it] = dEdxTrack[gdau21].dEdx();
                    grand_T4dedx22[it] = dEdxTrack[gdau22].dEdx();
                }

                //track pt
                grand_pt21[it] = gd21->pt();
                grand_pt22[it] = gd22->pt();

                //track momentum
                grand_p21[it] = gd21->p();
                grand_p22[it] = gd22->p();

                //track eta
                grand_eta21[it] = gd21->eta();
                grand_eta22[it] = gd22->eta();

                //track charge
                grand_charge21[it] = gd21->charge();
                grand_charge22[it] = gd22->charge();

                //track Chi2
                grand_trkChi21[it] = gdau21->normalizedChi2();
                grand_trkChi22[it] = gdau22->normalizedChi2();

                //track pT error
                grand_ptErr21[it] = gdau21->ptError();
                grand_ptErr22[it] = gdau22->ptError();

                //vertexCovariance 00-xError 11-y 22-z
                secvz = d2->vz(); secvx = d2->vx(); secvy = d2->vy();

                //trkNHits
                grand_nhit21[it] = gdau21->numberOfValidHits();
                grand_nhit22[it] = gdau22->numberOfValidHits();

                //DCA
                // math::XYZPoint bestvtx2(bestvx,bestvy,bestvz);

                double gdzbest21 = gdau21->dz(bestvtx);
                double gdxybest21 = gdau21->dxy(bestvtx);
                double gdzerror21 = sqrt(gdau21->dzError()*gdau21->dzError()+bestvzError*bestvzError);
                double gdxyerror21 = sqrt(gdau21->d0Error()*gdau21->d0Error()+bestvxError*bestvyError);

                grand_dzos21[it] = gdzbest21/gdzerror21;
                grand_dxyos21[it] = gdxybest21/gdxyerror21;

                double gdzbest22 = gdau22->dz(bestvtx);
                double gdxybest22 = gdau22->dxy(bestvtx);
                double gdzerror22 = sqrt(gdau22->dzError()*gdau22->dzError()+bestvzError*bestvzError);
                double gdxyerror22 = sqrt(gdau22->d0Error()*gdau22->d0Error()+bestvxError*bestvyError);

                grand_dzos22[it] = gdzbest22/gdzerror22;
                grand_dxyos22[it] = gdxybest22/gdxyerror22;

                //vtxChi2
                grand_vtxChi22[it] = d2->vertexChi2();
                grand_ndf2[it] = d2->vertexNdof();
                grand_VtxProb2[it] = TMath::Prob(grand_vtxChi22[it],grand_ndf2[it]);

                //PAngle
                TVector3 ptosvec2(secvx-bestvx,secvy-bestvy,secvz-bestvz);
                TVector3 secvec2(d2->px(),d2->py(),d2->pz());

                TVector3 ptosvec2D2(secvx-bestvx,secvy-bestvy,0);
                TVector3 secvec2D2(d2->px(),d2->py(),0);

                grand_agl2[it] = cos(secvec.Angle(ptosvec2));
                grand_agl_abs2[it] = secvec.Angle(ptosvec2);

                grand_agl2D2[it] = cos(secvec2D2.Angle(ptosvec2D2));
                grand_agl2D_abs2[it] = secvec2D2.Angle(ptosvec2D2);

                //Decay length 3D
                typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                typedef ROOT::Math::SVector<double, 3> SVector3;
                typedef ROOT::Math::SVector<double, 6> SVector6;

                SMatrixSym3D totalCov2 = vtx.covariance() + d2->vertexCovariance();
                SVector3 distanceVector2(secvx-bestvx,secvy-bestvy,secvz-bestvz);

                grand_dl2[it] = ROOT::Math::Mag(distanceVector2);
                grand_dlerror2[it] = sqrt(ROOT::Math::Similarity(totalCov2, distanceVector2))/grand_dl2[it];

                grand_dlos2[it] = grand_dl2[it]/grand_dlerror2[it];

                //Decay length 2D
                SVector6 v21(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
                SVector6 v22(d1->vertexCovariance(0,0), d1->vertexCovariance(0,1),d1->vertexCovariance(1,1),0,0,0);

                SMatrixSym3D sv21(v1);
                SMatrixSym3D sv22(v2);

                SMatrixSym3D totalCov2D2 = sv21 + sv22;
                SVector3 distanceVector2D2(secvx-bestvx,secvy-bestvy,0);

                double gdl2D2 = ROOT::Math::Mag(distanceVector2D2);
                double gdl2Derror2 = sqrt(ROOT::Math::Similarity(totalCov2D2, distanceVector2D2))/gdl2D2;

                grand_dlos2D2[it] = gdl2D2/gdl2Derror2;
            }
        }

        if(saveHistogram_)
        {
          for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
            for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
            {
              if(pt[it]<pTBins_[ipt+1] && pt[it]>pTBins_[ipt] && y[it]<yBins_[iy+1] && y[it]>yBins_[iy])
              {
                hMassVsMVA[iy][ipt]->Fill(mva[it],mass[it]);
//                h3DDCAVsMVA[iy][ipt]->Fill(mva[it],dl[it]*sin(agl_abs[it]));
//                h2DDCAVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]*sin(agl2D_abs[it]));

                if(saveAllHistogram_)
                {
                hpTVsMVA[iy][ipt]->Fill(mva[it],pt[it]);
                hetaVsMVA[iy][ipt]->Fill(mva[it],eta[it]);
                hyVsMVA[iy][ipt]->Fill(mva[it],y[it]);
                hVtxProbVsMVA[iy][ipt]->Fill(mva[it],VtxProb[it]);
                h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl[it]);
                h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl_abs[it]);
                h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl2D[it]);
                h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it],agl2D_abs[it]);
                h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it],dlos[it]);
                h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it],dl[it]);
                h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it],dlos2D[it]);
                h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]);
                hzDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it],dzos1[it]);
                hxyDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it],dxyos1[it]);
                hNHitD1VsMVA[iy][ipt]->Fill(mva[it],nhit1[it]);
                hpTD1VsMVA[iy][ipt]->Fill(mva[it],pt1[it]);
                hpTerrD1VsMVA[iy][ipt]->Fill(mva[it],ptErr1[it]/pt1[it]);
                hEtaD1VsMVA[iy][ipt]->Fill(mva[it],eta1[it]);
                hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva[it],H2dedx1[it]);
                hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1[it],H2dedx1[it]);
                hzDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it],dzos2[it]);
                hxyDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it],dxyos2[it]);
                hNHitD2VsMVA[iy][ipt]->Fill(mva[it],nhit2[it]);
                hpTD2VsMVA[iy][ipt]->Fill(mva[it],pt2[it]);
                hpTerrD2VsMVA[iy][ipt]->Fill(mva[it],ptErr2[it]/pt2[it]);
                hEtaD2VsMVA[iy][ipt]->Fill(mva[it],eta2[it]);
                hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva[it],H2dedx2[it]);
                hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2[it],H2dedx2[it]);
                if(threeProngDecay_)
                {
                  hzDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it],dzos3[it]);
                  hxyDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it],dxyos3[it]);
                  hNHitD3VsMVA[iy][ipt]->Fill(mva[it],nhit3[it]);
                  hpTD3VsMVA[iy][ipt]->Fill(mva[it],pt3[it]);
                  hpTerrD3VsMVA[iy][ipt]->Fill(mva[it],ptErr3[it]/pt3[it]);
                  hEtaD3VsMVA[iy][ipt]->Fill(mva[it],eta3[it]);
                  hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva[it],H2dedx3[it]);
                  hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1[it],H2dedx3[it]);
                }

                }
              }
            }
        }

    }
}

void
VertexCompositeTreeProducerNew::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);

    candSize_gen = 0;
    for(unsigned it=0; it<genpars->size(); ++it){

        const reco::GenParticle & trk = (*genpars)[it];

        const reco::Candidate * Dd1 = trk.daughter(0);
        const reco::Candidate * Dd2 = trk.daughter(1);

        int id = trk.pdgId();
        //if( fabs(id)!=PID_) continue; //check is target
        if( !doGenDoubleDecay_ &&fabs(id)!=PID_) continue; //check is target
        if( Dd1 == nullptr || Dd2 == nullptr ) continue; //check is target
        if( doGenDoubleDecay_ && !(abs(Dd1->pdgId()) == 421 && abs(Dd2->pdgId()) == 421) ) continue; //check is target
if( debug_ )std::cout << "pass id, decay dau : " << trk.numberOfDaughters() << std::endl;

        if(decayInGen_ && (trk.numberOfDaughters()!=2 && trk.numberOfDaughters()!=3)) continue; //check 2-pron decay if target decays in Gen

if( debug_ ) std::cout << "pass decay" << std::endl;

        candSize_gen+=1;

        pt_gen[candSize_gen-1] = trk.pt();
        mass_gen[candSize_gen-1] = trk.mass();
        eta_gen[candSize_gen-1] = trk.eta();
        phi_gen[candSize_gen-1] = trk.phi();
        status_gen[candSize_gen-1] = trk.status();
        idself[candSize_gen-1] = trk.pdgId();
        idmom[candSize_gen-1] = -77;
        y_gen[candSize_gen-1] = trk.rapidity();
        ptmom[candSize_gen-1] = -999.0;
        etamom[candSize_gen-1] = -999.0;
        phimom[candSize_gen-1] = -999.0;
        ymom[candSize_gen-1] = -999.0;
        statusmom[candSize_gen-1] = -999;

        if(trk.numberOfMothers()!=0)
        {
            const reco::Candidate * mom = trk.mother();
            idmom[candSize_gen-1] = mom->pdgId();
            ptmom[candSize_gen-1] = mom->pt();
            etamom[candSize_gen-1] = mom->eta();
            phimom[candSize_gen-1] = mom->phi();
            ymom[candSize_gen-1] = mom->rapidity();
            statusmom[candSize_gen-1] = mom->status();
        }

        if(!decayInGen_) continue;

        const reco::Candidate * Dd3 = trk.daughter(2);

        iddau1[candSize_gen-1] = fabs(Dd1->pdgId());
        iddau2[candSize_gen-1] = fabs(Dd2->pdgId());
        if(Dd3) iddau3[candSize_gen-1] = fabs(Dd3->pdgId());
        if( doGenDoubleDecay_){
          pt_gen1[candSize_gen-1] = Dd1->pt();
          mass_gen1[candSize_gen-1] = Dd1->mass();
          eta_gen1[candSize_gen-1] = Dd1->eta();
          phi_gen1[candSize_gen-1] = Dd1->phi();
          status_gen1[candSize_gen-1] = Dd1->status();
          idself1[candSize_gen-1] = Dd1->pdgId();

          pt_gen2[candSize_gen-1] = Dd2->pt();
          mass_gen2[candSize_gen-1] = Dd2->mass();
          eta_gen2[candSize_gen-1] = Dd2->eta();
          phi_gen2[candSize_gen-1] = Dd2->phi();
          status_gen2[candSize_gen-1] = Dd2->status();
          idself2[candSize_gen-1] = Dd2->pdgId();
        }
    }
}

// ------------ method called once each job just before starting event
//loop  ------------
void
VertexCompositeTreeProducerNew::beginJob()
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
VertexCompositeTreeProducerNew::initHistogram()
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
   hzDCASignificancedaughter1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificancedaughter1VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificancedaughter1;",100,-1.,1.,100,-10,10);
   hxyDCASignificancedaughter1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificancedaughter1VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificancedaughter1;",100,-1.,1.,100,-10,10);
   hNHitD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD1VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD1;",100,-1.,1.,100,0,100);
   hpTD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTD1;",100,-1.,1.,100,0,10);
   hpTerrD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD1VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD1;",100,-1.,1.,50,0,0.5);
   hEtaD1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD1VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD1;",100,-1.,1.,40,-4,4);
   hdedxHarmonic2D1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D1;",100,-1.,1.,100,0,10);
   hdedxHarmonic2D1VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D1",100,0,10,100,0,10);
   hzDCASignificancedaughter2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificancedaughter2VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificancedaughter2;",100,-1.,1.,100,-10,10);
   hxyDCASignificancedaughter2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificancedaughter2VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificancedaughter2;",100,-1.,1.,100,-10,10);
   hNHitD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hNHitD2VsMVA_y%d_pt%d",iy,ipt),";mva;NHitD2;",100,-1.,1.,100,0,100);
   hpTD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTD2;",100,-1.,1.,100,0,10);
   hpTerrD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTerrD2VsMVA_y%d_pt%d",iy,ipt),";mva;pTerrD2;",100,-1.,1.,50,0,0.5);
   hEtaD2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hEtaD2VsMVA_y%d_pt%d",iy,ipt),";mva;EtaD2;",100,-1.,1.,40,-4,4);
   hdedxHarmonic2D2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsMVA_y%d_pt%d",iy,ipt),";mva;dedxHarmonic2D2;",100,-1.,1.,100,0,10);
   hdedxHarmonic2D2VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsP_y%d_pt%d",iy,ipt),";p (GeV);dedxHarmonic2D2",100,0,10,100,0,10);

   if(threeProngDecay_)
   {
     hzDCASignificancedaughter3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hzDCASignificancedaughter3VsMVA_y%d_pt%d",iy,ipt),";mva;zDCASignificancedaughter3;",100,-1.,1.,100,-10,10);
     hxyDCASignificancedaughter3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificancedaughter3VsMVA_y%d_pt%d",iy,ipt),";mva;xyDCASignificancedaughter3;",100,-1.,1.,100,-10,10);
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
VertexCompositeTreeProducerNew::initTree()
{ 
    VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");
    
    if(doRecoNtuple_) 
    { 
  
    // Event info
    VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    VertexCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
    VertexCompositeNtuple->Branch("HFsumET",&HFsumET,"HFsumET/F");
    VertexCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    VertexCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    VertexCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    VertexCompositeNtuple->Branch("candSize",&candSize,"candSize/I");
    if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality,"centrality/I");
    if(isEventPlane_) 
    {
      VertexCompositeNtuple->Branch("ephfpAngle",&ephfpAngle,"ephfpAngle[3]/F");
      VertexCompositeNtuple->Branch("ephfmAngle",&ephfmAngle,"ephfmAngle[3]/F");
      VertexCompositeNtuple->Branch("ephfpQ",&ephfpQ,"ephfpQ[3]/F");
      VertexCompositeNtuple->Branch("ephfmQ",&ephfmQ,"ephfmQ[3]/F");
      VertexCompositeNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
      VertexCompositeNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
    }

    // particle info
    VertexCompositeNtuple->Branch("pT",&pt,"pT[candSize]/F");
    VertexCompositeNtuple->Branch("y",&y,"y[candSize]/F");
    VertexCompositeNtuple->Branch("eta",&eta,"eta[candSize]/F");
    VertexCompositeNtuple->Branch("phi",&phi,"phi[candSize]/F");
    VertexCompositeNtuple->Branch("mass",&mass,"mass[candSize]/F");
    if(useAnyMVA_){
      VertexCompositeNtuple->Branch("mva",&mva,"mva[candSize]/F");
      if( doubleCand_ ){
        VertexCompositeNtuple->Branch("mvaDaughter1",&mva1,"mvaDaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("mvaDaughter2",&mva2,"mvaDaughter2[candSize]/F");

      }
    }
    if(useDCA_) {
      VertexCompositeNtuple->Branch("dca3D", &dca3D, "dca3D[candSize]/F");
      VertexCompositeNtuple->Branch("dcaErr3D", &dcaErr3D, "dcaErr3D[candSize]/F");
    }

    if(!isSkimMVA_)  
    {
        //Composite candidate info RECO
        VertexCompositeNtuple->Branch("flavor",&flavor,"flavor[candSize]/F");
        VertexCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb[candSize]/F");
//        VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2[candSize]/F");
//        VertexCompositeNtuple->Branch("VtxNDF",&ndf,"VtxNDF[candSize]/F");
        VertexCompositeNtuple->Branch("3DCosPointingAngle",&agl,"3DCosPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("2DCosPointingAngle",&agl2D,"2DCosPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength[candSize]/F");
        VertexCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance[candSize]/F");
        VertexCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength[candSize]/F");
    
        if(doGenMatching_)
        {
            VertexCompositeNtuple->Branch("isSwap",&isSwap,"isSwap[candSize]/O");
            VertexCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco[candSize]/I");
            VertexCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN[candSize]/O");
            if( doGenDoubleDecay_){
              VertexCompositeNtuple->Branch("isSwap1",&isSwap1,"isSwap1[candSize]/O");
              VertexCompositeNtuple->Branch("matchGEN1",&matchGEN1,"matchGEN1[candSize]/O");
              VertexCompositeNtuple->Branch("isSwap2",&isSwap2,"isSwap2[candSize]/O");
              VertexCompositeNtuple->Branch("matchGEN2",&matchGEN2,"matchGEN2[candSize]/O");
            }
        }
        
        if(doGenMatchingTOF_)
        {
          VertexCompositeNtuple->Branch("PIDD1",&pid1,"PIDD1[candSize]/I");
          VertexCompositeNtuple->Branch("PIDD2",&pid1,"PIDD2[candSize]/I");
          VertexCompositeNtuple->Branch("TOFD1",&tof1,"TOFD1[candSize]/F");
          VertexCompositeNtuple->Branch("TOFD2",&tof1,"TOFD2[candSize]/F");
        }

        //daughter & grand daughter info
        if(twoLayerDecay_)
        {
            VertexCompositeNtuple->Branch("flavordaughter1",&flavor1,"flavordaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("massdaughter1",&grand_mass,"massdaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1[candSize]/F");
            VertexCompositeNtuple->Branch("PhiD1",&phi1,"PhiD1[candSize]/F");
            VertexCompositeNtuple->Branch("VtxProbdaughter1",&grand_VtxProb,"VtxProbdaughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxNDFdaughter1",&grand_ndf,"VtxNDFdaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DCosPointingAngle1daughter1",&grand_agl,"3DCosPointingAngledaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DPointingAngledaughter1",&grand_agl_abs,"3DPointingAngledaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("2DCosPointingAngledaughter1",&grand_agl2D,"2DCosPointingAngledaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("2DPointingAngledaughter1",&grand_agl2D_abs,"2DPointingAngledaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthSignificancedaughter1",&grand_dlos,"3DDecayLengthSignificancedaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthdaughter1",&grand_dl,"3DDecayLengthdaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthErrordaughter1",&grand_dlerror,"3DDecayLengthErrordaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("2DDecayLengthSignificancedaughter1",&grand_dlos2D,"2DDecayLengthSignificancedaughter1[candSize]/F");
            
            VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter11",&grand_dzos1,"zDCASignificanceGranddaughter11[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter12",&grand_dzos2,"zDCASignificanceGranddaughter12[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter11",&grand_dxyos1,"xyDCASignificanceGranddaughter11[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter12",&grand_dxyos2,"xyDCASignificanceGranddaughter12[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD11",&grand_nhit1,"NHitGrandD11[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD12",&grand_nhit2,"NHitGrandD12[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityGranddaughter11",&grand_trkquality1,"HighPurityGranddaughter11[candSize]/O");
            VertexCompositeNtuple->Branch("HighPurityGranddaughter12",&grand_trkquality2,"HighPurityGranddaughter12[candSize]/O");
            VertexCompositeNtuple->Branch("pTGrandD11",&grand_pt1,"pTGrandD11[candSize]/F");
            VertexCompositeNtuple->Branch("pTGrandD12",&grand_pt2,"pTGrandD12[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD11",&grand_ptErr1,"pTerrGrandD11[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD12",&grand_ptErr2,"pTerrGrandD12[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD11",&grand_p1,"pGrandD11[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD12",&grand_p2,"pGrandD12[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD11",&grand_eta1,"EtaGrandD11[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD12",&grand_eta2,"EtaGrandD12[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeGrandD11",&grand_charge1,"chargeGrandD11[candSize]/I");
//            VertexCompositeNtuple->Branch("chargeGrandD12",&grand_charge2,"chargeGrandD12[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD11",&grand_H2dedx1,"dedxHarmonic2GrandD11[candSize]/F");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD12",&grand_H2dedx2,"dedxHarmonic2GrandD12[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter1",&grand_T4dedx1,"dedxTruncated40Granddaughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter2",&grand_T4dedx2,"dedxTruncated40Granddaughter2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter1",&grand_trkChi1,"normalizedChi2Granddaughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter2",&grand_trkChi2,"normalizedChi2Granddaughter2[candSize]/F");
            if( doubleCand_ ){
            VertexCompositeNtuple->Branch("flavordaughter2",&flavor2,"flavordaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("massdaughter2",&grand_mass2,"massdaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2[candSize]/F");
            VertexCompositeNtuple->Branch("PhiD2",&phi2,"PhiD2[candSize]/F");
            VertexCompositeNtuple->Branch("VtxProbdaughter2",&grand_VtxProb2,"VtxProbdaughter2[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("VtxNDFdaughter1",&grand_ndf,"VtxNDFdaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("3DCosPointingAngle1daughter2",&grand_agl2,"3DCosPointingAngledaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("3DPointingAngledaughter2",&grand_agl_abs2,"3DPointingAngledaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("2DCosPointingAngledaughter2",&grand_agl2D2,"2DCosPointingAngledaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("2DPointingAngledaughter2",&grand_agl2D_abs2,"2DPointingAngledaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthSignificancedaughter2",&grand_dlos2,"3DDecayLengthSignificancedaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthdaughter2",&grand_dl2,"3DDecayLengthdaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("3DDecayLengthErrordaughter2",&grand_dlerror2,"3DDecayLengthErrordaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("2DDecayLengthSignificancedaughter2",&grand_dlos2D2,"2DDecayLengthSignificancedaughter2[candSize]/F");

            VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter21",&grand_dzos21,"zDCASignificanceGranddaughter21[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter22",&grand_dzos22,"zDCASignificanceGranddaughter22[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter21",&grand_dxyos21,"xyDCASignificanceGranddaughter21[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter22",&grand_dxyos22,"xyDCASignificanceGranddaughter22[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD21",&grand_nhit21,"NHitGrandD21[candSize]/F");
            VertexCompositeNtuple->Branch("NHitGrandD22",&grand_nhit22,"NHitGrandD22[candSize]/F");
            VertexCompositeNtuple->Branch("HighPurityGranddaughter21",&grand_trkquality21,"HighPurityGranddaughter21[candSize]/O");
            VertexCompositeNtuple->Branch("HighPurityGranddaughter22",&grand_trkquality22,"HighPurityGranddaughter22[candSize]/O");
            VertexCompositeNtuple->Branch("pTGrandD21",&grand_pt21,"pTGrandD21[candSize]/F");
            VertexCompositeNtuple->Branch("pTGrandD22",&grand_pt22,"pTGrandD22[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD21",&grand_ptErr21,"pTerrGrandD21[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrGrandD22",&grand_ptErr22,"pTerrGrandD22[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD21",&grand_p21,"pGrandD1[candSize]/F");
//            VertexCompositeNtuple->Branch("pGrandD22",&grand_p22,"pGrandD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD21",&grand_eta21,"EtaGrandD21[candSize]/F");
            VertexCompositeNtuple->Branch("EtaGrandD22",&grand_eta22,"EtaGrandD22[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeGrandD21",&grand_charge21,"chargeGrandD1[candSize]/I");
//            VertexCompositeNtuple->Branch("chargeGrandD22",&grand_charge22,"chargeGrandD2[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD21",&grand_H2dedx21,"dedxHarmonic2GrandD21[candSize]/F");
            VertexCompositeNtuple->Branch("dedxHarmonic2GrandD22",&grand_H2dedx22,"dedxHarmonic2GrandD22[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter1",&grand_T4dedx21,"dedxTruncated40Granddaughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter2",&grand_T4dedx22,"dedxTruncated40Granddaughter2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter1",&grand_trkChi21,"normalizedChi2Granddaughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter2",&grand_trkChi22,"normalizedChi2Granddaughter2[candSize]/F");
            }
        }
        else
        {
            VertexCompositeNtuple->Branch("zDCASignificancedaughter1",&dzos1,"zDCASignificancedaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificancedaughter1",&dxyos1,"xyDCASignificancedaughter1[candSize]/F");
            VertexCompositeNtuple->Branch("NHitD1",&nhit1,"NHitD1[candSize]/F");
            VertexCompositeNtuple->Branch("HighPuritydaughter1",&trkquality1,"HighPuritydaughter1[candSize]/O");
            VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrD1",&ptErr1,"pTerrD1[candSize]/F");
//            VertexCompositeNtuple->Branch("pD1",&p1,"pD1[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1[candSize]/F");
            VertexCompositeNtuple->Branch("PhiD1",&phi1,"PhiD1[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeD1",&charge1,"chargeD1[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D1",&H2dedx1,"dedxHarmonic2D1[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40daughter1",&T4dedx1,"dedxTruncated40daughter1[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2daughter1",&trkChi1,"normalizedChi2daughter1[candSize]/F");
            VertexCompositeNtuple->Branch("zDCASignificancedaughter2",&dzos2,"zDCASignificancedaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("xyDCASignificancedaughter2",&dxyos2,"xyDCASignificancedaughter2[candSize]/F");
            VertexCompositeNtuple->Branch("NHitD2",&nhit2,"NHitD2[candSize]/F");
            VertexCompositeNtuple->Branch("HighPuritydaughter2",&trkquality2,"HighPuritydaughter2[candSize]/O");
            VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2[candSize]/F");
            VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2[candSize]/F");
//            VertexCompositeNtuple->Branch("pD2",&p2,"pD2[candSize]/F");
            VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2[candSize]/F");
            VertexCompositeNtuple->Branch("PhiD2",&phi2,"PhiD2[candSize]/F");
//            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2[candSize]/I");
            VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2[candSize]/F");
//            VertexCompositeNtuple->Branch("dedxTruncated40daughter2",&T4dedx2,"dedxTruncated40daughter2[candSize]/F");
//            VertexCompositeNtuple->Branch("normalizedChi2daughter2",&trkChi2,"normalizedChi2daughter2[candSize]/F");
            if(threeProngDecay_)
            {
              VertexCompositeNtuple->Branch("zDCASignificancedaughter3",&dzos3,"zDCASignificancedaughter3[candSize]/F");
              VertexCompositeNtuple->Branch("xyDCASignificancedaughter3",&dxyos3,"xyDCASignificancedaughter3[candSize]/F");
              VertexCompositeNtuple->Branch("NHitD3",&nhit3,"NHitD3[candSize]/F");
              VertexCompositeNtuple->Branch("HighPuritydaughter3",&trkquality3,"HighPuritydaughter3[candSize]/O");
              VertexCompositeNtuple->Branch("pTD3",&pt1,"pTD3[candSize]/F");
              VertexCompositeNtuple->Branch("pTerrD3",&ptErr3,"pTerrD3[candSize]/F");
              VertexCompositeNtuple->Branch("EtaD3",&eta1,"EtaD3[candSize]/F");
              VertexCompositeNtuple->Branch("dedxHarmonic2D3",&H2dedx1,"dedxHarmonic2D3[candSize]/F");
            }
        }
        
        if(doMuon_)
        {
            VertexCompositeNtuple->Branch("OneStMuon1",&onestmuon1,"OneStMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("OneStMuon2",&onestmuon2,"OneStMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("PFMuon1",&pfmuon1,"PFMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("PFMuon2",&pfmuon2,"PFMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("GlbMuon1",&glbmuon1,"GlbMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("GlbMuon2",&glbmuon2,"GlbMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("trkMuon1",&trkmuon1,"trkMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("trkMuon2",&trkmuon2,"trkMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("caloMuon1",&calomuon1,"caloMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("caloMuon2",&calomuon2,"caloMuon2[candSize]/O");
            VertexCompositeNtuple->Branch("SoftMuon1",&softmuon1,"SoftMuon1[candSize]/O");
            VertexCompositeNtuple->Branch("SoftMuon2",&softmuon2,"SoftMuon2[candSize]/O");

            if(doMuonFull_)
            {
              VertexCompositeNtuple->Branch("nMatchedChamberD1",&nmatchedch1,"nMatchedChamberD1[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedStationD1",&nmatchedst1,"nMatchedStationD1[candSize]/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD1",&matchedenergy1,"EnergyDepositionD1[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedChamberD2",&nmatchedch2,"nMatchedChamberD2[candSize]/F");
              VertexCompositeNtuple->Branch("nMatchedStationD2",&nmatchedst2,"nMatchedStationD2[candSize]/F");
              VertexCompositeNtuple->Branch("EnergyDepositionD2",&matchedenergy2,"EnergyDepositionD2[candSize]/F");
              VertexCompositeNtuple->Branch("dx1_seg",        &dx1_seg_, "dx1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dy1_seg",        &dy1_seg_, "dy1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dxSig1_seg",     &dxSig1_seg_, "dxSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dySig1_seg",     &dySig1_seg_, "dySig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdz1_seg",     &ddxdz1_seg_, "ddxdz1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydz1_seg",     &ddydz1_seg_, "ddydz1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdzSig1_seg",  &ddxdzSig1_seg_, "ddxdzSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydzSig1_seg",  &ddydzSig1_seg_, "ddydzSig1_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dx2_seg",        &dx2_seg_, "dx2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dy2_seg",        &dy2_seg_, "dy2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dxSig2_seg",     &dxSig2_seg_, "dxSig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("dySig2_seg",     &dySig2_seg_, "dySig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdz2_seg",     &ddxdz2_seg_, "ddxdz2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydz2_seg",     &ddydz2_seg_, "ddydz2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddxdzSig2_seg",  &ddxdzSig2_seg_, "ddxdzSig2_seg[candSize]/F");
              VertexCompositeNtuple->Branch("ddydzSig2_seg",  &ddydzSig2_seg_, "ddydzSig2_seg[candSize]/F");
           }
        }
    }

    } // doRecoNtuple_

    if(doGenNtuple_)
    {
        VertexCompositeNtuple->Branch("candSize_gen",&candSize_gen,"candSize_gen/I");
        VertexCompositeNtuple->Branch("id_gen",&idself,"id_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("mass_gen",&mass_gen,"mass_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("pT_gen",&pt_gen,"pT_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("eta_gen",&eta_gen,"eta_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("phi_gen",&phi_gen,"phi_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("y_gen",&y_gen,"y_gen[candSize_gen]/F");
        VertexCompositeNtuple->Branch("status_gen",&status_gen,"status_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherID_gen",&idmom,"MotherID_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherPt_gen",&ptmom,"MotherPt_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherEta_gen",&etamom,"MotherEta_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherPhi_gen",&phimom,"MotherPhi_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherY_gen",&ymom,"MotherY_gen[candSize_gen]/I");
        VertexCompositeNtuple->Branch("MotherStatus_gen",&statusmom,"MotherStatus_gen[candSize_gen]/I");
        if(doGenDoubleDecay_){
          VertexCompositeNtuple->Branch("id_gen1",&idself1,"id_gen1[candSize_gen]/I");
          VertexCompositeNtuple->Branch("mass_gen1",&mass_gen1,"mass_gen1[candSize_gen]/I");
          VertexCompositeNtuple->Branch("pt_gen1",&pt_gen1,"pt_gen1[candSize_gen]/I");
          VertexCompositeNtuple->Branch("eta_gen1",&eta_gen1,"eta_gen1[candSize_gen]/I");
          VertexCompositeNtuple->Branch("phi_gen1",&phi_gen1,"phi_gen1[candSize_gen]/I");
          VertexCompositeNtuple->Branch("status_gen1",&status_gen1,"status_gen1[candSize_gen]/I");

          VertexCompositeNtuple->Branch("id_gen2",&idself2,"id_gen2[candSize_gen]/I");
          VertexCompositeNtuple->Branch("mass_gen2",&mass_gen2,"mass_gen2[candSize_gen]/I");
          VertexCompositeNtuple->Branch("pt_gen2",&pt_gen2,"pt_gen2[candSize_gen]/I");
          VertexCompositeNtuple->Branch("eta_gen2",&eta_gen2,"eta_gen2[candSize_gen]/I");
          VertexCompositeNtuple->Branch("phi_gen2",&phi_gen2,"phi_gen2[candSize_gen]/I");
          VertexCompositeNtuple->Branch("status_gen2",&status_gen2,"status_gen2[candSize_gen]/I");

        }

        if(decayInGen_)
        {

            VertexCompositeNtuple->Branch("DauID1_gen",&iddau1,"DauID1_gen[candSize_gen]/I");
            VertexCompositeNtuple->Branch("DauID2_gen",&iddau2,"DauID2_gen[candSize_gen]/I");
            VertexCompositeNtuple->Branch("DauID3_gen",&iddau3,"DauID3_gen[candSize_gen]/I");
        }
    }
}

int VertexCompositeTreeProducerNew::
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
VertexCompositeTreeProducerNew::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducerNew);
