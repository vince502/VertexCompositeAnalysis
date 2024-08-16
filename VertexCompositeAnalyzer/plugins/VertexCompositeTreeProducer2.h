
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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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

//#define DEBUG true


//
// class decleration
//

#define PI 3.1416
#define MAXCAN 50000

using namespace std;

class VertexCompositeTreeProducer2 : public edm::one::EDAnalyzer<> {
public:
  explicit VertexCompositeTreeProducer2(const edm::ParameterSet&);
  ~VertexCompositeTreeProducer2();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();
  void genDecayLength(const uint&, const reco::GenParticle&);

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
    float y[MAXCAN];
    float mass[MAXCAN];
    float VtxProb[MAXCAN];
    float dlos[MAXCAN];
    float dl[MAXCAN];
    float dlerror[MAXCAN];
    float agl[MAXCAN];
    float vtxChi2[MAXCAN];
    float ndf[MAXCAN];
    float agl_abs[MAXCAN];
    float agl2D[MAXCAN];
    float agl2D_abs[MAXCAN];
    float dlos2D[MAXCAN];
    float dl2D[MAXCAN];
    bool isSwap[MAXCAN];
    bool matchGEN[MAXCAN];
    int pionFlavor[MAXCAN];
    int idmom_reco[MAXCAN];
    float gen_agl_abs[MAXCAN];
    float gen_agl2D_abs[MAXCAN];
    float gen_dl[MAXCAN];
    float gen_dl2D[MAXCAN];
    
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
    float pt_gen[MAXCAN];
    float eta_gen[MAXCAN];
    int status_gen[MAXCAN];
    int idmom[MAXCAN];
    float y_gen[MAXCAN];
    int iddau1[MAXCAN];
    int iddau2[MAXCAN];
    int iddau3[MAXCAN];

    //vector for gen match
    vector< vector<double> > *pVect;
    vector< vector<double> > *gpVect;
    vector<double> *Dvector1;
    vector<double> *GDvector1;
    vector<double> *Dvector2;
    vector<double> *GDvector2;
    vector<double> *Dvector3;
    vector<int> *pVectIDmom;
    
    bool useAnyMVA_;
    bool isSkimMVA_;
    bool isCentrality_;
    bool isEventPlane_;

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

    edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;
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

