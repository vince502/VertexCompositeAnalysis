// system include files
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// #include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
// #include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
// #include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

// #define DEBUG2

//
// class decleration
//

#define PI 3.1416
#define MAXCAN 2000

using namespace std;

class VertexCompositeTreeProducerNew : public edm::EDAnalyzer {
public:
  explicit VertexCompositeTreeProducerNew(const edm::ParameterSet &);
  ~VertexCompositeTreeProducerNew();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void fillRECO(const edm::Event &, const edm::EventSetup &);
  virtual void fillGEN(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  virtual void initHistogram();
  virtual void initTree();

  int muAssocToTrack(const reco::TrackRef &trackref, const edm::Handle<reco::MuonCollection> &muonh) const;

  reco::GenParticleRef findMother(const reco::GenParticleRef &);
  void genDecayLength(const reco::GenParticle &gCand, float &gen_decayLength2D_, float &gen_decayLength3D_,
                      float &gen_angle2D_, float &gen_angle3D_);
  void getAncestorId(const reco::GenParticle &gCand, int &gen_ancestorId_, int &gen_ancestorFlavor_);

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree *VertexCompositeNtuple;
  TH2F *hMassVsMVA[6][10];
  TH2F *hpTVsMVA[6][10];
  TH2F *hetaVsMVA[6][10];
  TH2F *hyVsMVA[6][10];
  TH2F *hVtxProbVsMVA[6][10];
  TH2F *h3DCosPointingAngleVsMVA[6][10];
  TH2F *h3DPointingAngleVsMVA[6][10];
  TH2F *h2DCosPointingAngleVsMVA[6][10];
  TH2F *h2DPointingAngleVsMVA[6][10];
  TH2F *h3DDecayLengthSignificanceVsMVA[6][10];
  TH2F *h3DDecayLengthVsMVA[6][10];
  TH2F *h2DDecayLengthSignificanceVsMVA[6][10];
  TH2F *h2DDecayLengthVsMVA[6][10];
  TH2F *h3DDCAVsMVA[6][10];
  TH2F *h2DDCAVsMVA[6][10];
  TH2F *hzDCASignificancedaughter1VsMVA[6][10];
  TH2F *hxyDCASignificancedaughter1VsMVA[6][10];
  TH2F *hNHitD1VsMVA[6][10];
  TH2F *hpTD1VsMVA[6][10];
  TH2F *hpTerrD1VsMVA[6][10];
  TH2F *hEtaD1VsMVA[6][10];
  TH2F *hdedxHarmonic2D1VsMVA[6][10];
  TH2F *hdedxHarmonic2D1VsP[6][10];
  TH2F *hzDCASignificancedaughter2VsMVA[6][10];
  TH2F *hxyDCASignificancedaughter2VsMVA[6][10];
  TH2F *hNHitD2VsMVA[6][10];
  TH2F *hpTD2VsMVA[6][10];
  TH2F *hpTerrD2VsMVA[6][10];
  TH2F *hEtaD2VsMVA[6][10];
  TH2F *hdedxHarmonic2D2VsMVA[6][10];
  TH2F *hdedxHarmonic2D2VsP[6][10];
  TH2F *hzDCASignificancedaughter3VsMVA[6][10];
  TH2F *hxyDCASignificancedaughter3VsMVA[6][10];
  TH2F *hNHitD3VsMVA[6][10];
  TH2F *hpTD3VsMVA[6][10];
  TH2F *hpTerrD3VsMVA[6][10];
  TH2F *hEtaD3VsMVA[6][10];
  TH2F *hdedxHarmonic2D3VsMVA[6][10];
  TH2F *hdedxHarmonic2D3VsP[6][10];

  bool saveTree_;
  bool saveHistogram_;
  bool saveAllHistogram_;
  double massHistPeak_;
  double massHistWidth_;
  int massHistBins_;

  // options
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

  // cut variables
  double multMax_;
  double multMin_;
  double deltaR_; // deltaR for Gen matching

  vector<double> pTBins_;
  vector<double> yBins_;

  // tree branches
  // event info
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

  // Composite candidate info
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
  int idmom_reco1[MAXCAN];
  int idmom_reco2[MAXCAN];
  int idBAnc_reco1[MAXCAN];
  int idBAnc_reco2[MAXCAN];

  bool matchGEN1[MAXCAN];       // For Double Decay
  bool matchGEN2[MAXCAN];       // For Double Decay
  unsigned matchToGen1[MAXCAN]; // For Double Decay
  unsigned matchToGen2[MAXCAN]; // For Double Decay
  bool isSwap1[MAXCAN];         // For Double Decay
  bool isSwap2[MAXCAN];         // For Double Decay

  // dau candidate info
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

  // dau info
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

  // grand-dau info
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
  float grand_phi1[MAXCAN];
  float grand_phi2[MAXCAN];
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
  float grand_phi21[MAXCAN];
  float grand_phi22[MAXCAN];
  int grand_charge21[MAXCAN];
  int grand_charge22[MAXCAN];
  float grand_H2dedx21[MAXCAN];
  float grand_H2dedx22[MAXCAN];
  float grand_T4dedx21[MAXCAN];
  float grand_T4dedx22[MAXCAN];
  float grand_trkChi21[MAXCAN];
  float grand_trkChi22[MAXCAN];

  // dau muon info
  bool onestmuon1[MAXCAN];
  bool onestmuon2[MAXCAN];
  bool pfmuon1[MAXCAN];
  bool pfmuon2[MAXCAN];
  bool glbmuon1[MAXCAN];
  bool glbmuon2[MAXCAN];
  bool trkmuon1[MAXCAN];
  bool trkmuon2[MAXCAN];
  bool calomuon1[MAXCAN];
  bool calomuon2[MAXCAN];
  bool softmuon1[MAXCAN];
  bool softmuon2[MAXCAN];
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

  // gen information for # of daughters == 2
  int gen_D1ancestorFlavor_[MAXCAN];
  int gen_D1ancestorId_[MAXCAN];
  float gen_D1PVx_[MAXCAN];
  float gen_D1PVy_[MAXCAN];
  float gen_D1PVz_[MAXCAN];

  float gen_D1pT_[MAXCAN];
  float gen_D1eta_[MAXCAN];
  float gen_D1phi_[MAXCAN];
  float gen_D1mass_[MAXCAN];
  float gen_D1y_[MAXCAN];
  int gen_D1charge_[MAXCAN];
  int gen_D1pdgId_[MAXCAN];

  float gen_D1decayLength3D_[MAXCAN];
  float gen_D1decayLength2D_[MAXCAN];
  float gen_D1angle3D_[MAXCAN];
  float gen_D1angle2D_[MAXCAN];

  float gen_D1pTD1_[MAXCAN];
  float gen_D1etaD1_[MAXCAN];
  float gen_D1phiD1_[MAXCAN];
  float gen_D1massD1_[MAXCAN];
  float gen_D1yD1_[MAXCAN];
  float gen_D1chargeD1_[MAXCAN];
  float gen_D1pdgIdD1_[MAXCAN];

  float gen_D1pTD2_[MAXCAN];
  float gen_D1etaD2_[MAXCAN];
  float gen_D1phiD2_[MAXCAN];
  float gen_D1massD2_[MAXCAN];
  float gen_D1yD2_[MAXCAN];
  float gen_D1chargeD2_[MAXCAN];
  float gen_D1pdgIdD2_[MAXCAN];

  int gen_D2ancestorFlavor_[MAXCAN];
  int gen_D2ancestorId_[MAXCAN];
  float gen_D2PVx_[MAXCAN];
  float gen_D2PVy_[MAXCAN];
  float gen_D2PVz_[MAXCAN];

  float gen_D2pT_[MAXCAN];
  float gen_D2eta_[MAXCAN];
  float gen_D2phi_[MAXCAN];
  float gen_D2mass_[MAXCAN];
  float gen_D2y_[MAXCAN];
  int gen_D2charge_[MAXCAN];
  int gen_D2pdgId_[MAXCAN];

  float gen_D2decayLength3D_[MAXCAN];
  float gen_D2decayLength2D_[MAXCAN];
  float gen_D2angle3D_[MAXCAN];
  float gen_D2angle2D_[MAXCAN];

  float gen_D2pTD1_[MAXCAN];
  float gen_D2etaD1_[MAXCAN];
  float gen_D2phiD1_[MAXCAN];
  float gen_D2massD1_[MAXCAN];
  float gen_D2yD1_[MAXCAN];
  int gen_D2chargeD1_[MAXCAN];
  int gen_D2pdgIdD1_[MAXCAN];

  float gen_D2pTD2_[MAXCAN];
  float gen_D2etaD2_[MAXCAN];
  float gen_D2phiD2_[MAXCAN];
  float gen_D2massD2_[MAXCAN];
  float gen_D2yD2_[MAXCAN];
  int gen_D2chargeD2_[MAXCAN];
  int gen_D2pdgIdD2_[MAXCAN];

  int idself[MAXCAN];
  float mass_gen[MAXCAN];
  float pt_gen[MAXCAN];
  float eta_gen[MAXCAN];
  float phi_gen[MAXCAN];
  float dl2D_gen[MAXCAN];
  float dl3D_gen[MAXCAN];
  float angle2D_gen[MAXCAN];
  float angle3D_gen[MAXCAN];
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

  // vector for gen match
  vector<vector<double>> *pVect;

  vector<vector<double>> *pVectg1;
  vector<vector<double>> *pVectg2;
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

  reco::Particle::Point genVertex_;

  // tokens
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
  edm::EDGetTokenT<MVACollection> MVAValues_Token_;
  edm::EDGetTokenT<MVACollection> MVAValues_Token2_;

  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token1_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token2_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
  edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

  edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

  // for DCA
  edm::EDGetTokenT<std::vector<float>> tok_DCAVal_;
  edm::EDGetTokenT<std::vector<float>> tok_DCAErr_;
};

//
// constants, enums and typedefs
//
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// static data member definitions
//

//
// constructors and destructor
//

VertexCompositeTreeProducerNew::VertexCompositeTreeProducerNew(const edm::ParameterSet &iConfig) {
  // options
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
  if (threeProngDecay_)
    PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");
  if (doGenDoubleDecay_) {
    PID_dau1_grand1_ = iConfig.getUntrackedParameter<int>("PID_dau1_grand1");
    PID_dau1_grand2_ = iConfig.getUntrackedParameter<int>("PID_dau1_grand2");
    PID_dau2_grand1_ = iConfig.getUntrackedParameter<int>("PID_dau2_grand1");
    PID_dau2_grand2_ = iConfig.getUntrackedParameter<int>("PID_dau2_grand2");
  }

  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");
  saveAllHistogram_ = iConfig.getUntrackedParameter<bool>("saveAllHistogram");
  massHistPeak_ = iConfig.getUntrackedParameter<double>("massHistPeak");
  massHistWidth_ = iConfig.getUntrackedParameter<double>("massHistWidth");
  massHistBins_ = iConfig.getUntrackedParameter<int>("massHistBins");

  useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
  isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA");

  // cut variables
  multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
  multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
  deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);

  pTBins_ = iConfig.getUntrackedParameter<std::vector<double>>("pTBins");
  yBins_ = iConfig.getUntrackedParameter<std::vector<double>>("yBins");

  // input tokens
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
  recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(
      iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
  MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  MVAValues_Token2_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection2"));
  tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
  Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
  Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxTruncated40"));
  tok_genParticle_ = consumes<reco::GenParticleCollection>(
      edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

  isCentrality_ = false;
  if (iConfig.exists("isCentrality"))
    isCentrality_ = iConfig.getParameter<bool>("isCentrality");
  if (isCentrality_) {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  isEventPlane_ = false;
  if (iConfig.exists("isEventPlane"))
    isEventPlane_ = iConfig.getParameter<bool>("isEventPlane");
  if (isEventPlane_) {
    tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));
  }

  if (useAnyMVA_ && iConfig.exists("MVACollection"))
    MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  if (useAnyMVA_ && iConfig.exists("MVACollection2"))
    MVAValues_Token2_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection2"));
  if (iConfig.exists("DCAValCollection") && iConfig.exists("DCAErrCollection")) {
    useDCA_ = true;
    tok_DCAVal_ = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("DCAValCollection"));
    tok_DCAErr_ = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("DCAErrCollection"));
  }
}

VertexCompositeTreeProducerNew::~VertexCompositeTreeProducerNew() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------
void VertexCompositeTreeProducerNew::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;

  if (doGenNtuple_)
    fillGEN(iEvent, iSetup);
  if (doRecoNtuple_)
    fillRECO(iEvent, iSetup);

  if (saveTree_)
    VertexCompositeNtuple->Fill();
}

void VertexCompositeTreeProducerNew::fillRECO(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tok_generalTrk_, tracks);

  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
  iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_, v0candidates);
  const reco::VertexCompositeCandidateCollection *v0candidates_ = v0candidates.product();

  edm::Handle<MVACollection> mvavalues;
  edm::Handle<MVACollection> mvavalues2;
  if (useAnyMVA_) {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    assert((*mvavalues).size() == v0candidates->size());
    if (doubleCand_) {
      iEvent.getByToken(MVAValues_Token2_, mvavalues2);
      assert((*mvavalues2).size() == v0candidates->size());
    }
  }
  edm::Handle<std::vector<float>> dcaValues;
  edm::Handle<std::vector<float>> dcaErrors;
  if (useDCA_) {
    iEvent.getByToken(tok_DCAVal_, dcaValues);
    iEvent.getByToken(tok_DCAErr_, dcaErrors);
    assert((*dcaValues).size() == v0candidates->size());
    assert((*dcaErrors).size() == v0candidates->size());
  }

  edm::Handle<reco::GenParticleCollection> genpars;
  if (doGenMatching_ || doGenMatchingTOF_)
    iEvent.getByToken(tok_genParticle_, genpars);

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle1;
  iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle2;
  iEvent.getByToken(Dedx_Token2_, dEdxHandle2);

  centrality = -1;
  if (isCentrality_) {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);

    iEvent.getByToken(tok_centBinLabel_, cbin_);
    centrality = *cbin_;

    HFsumET = cent->EtHFtowerSum();
    Npixel = cent->multiplicityPixel();
    //      int ntrk = cent->Ntracks();
  }

  if (isEventPlane_) {
    edm::Handle<reco::EvtPlaneCollection> eventplanes;
    iEvent.getByToken(tok_eventplaneSrc_, eventplanes);

    const reco::EvtPlane &ephfp1 = (*eventplanes)[0];
    const reco::EvtPlane &ephfm1 = (*eventplanes)[1];
    const reco::EvtPlane &ephfp2 = (*eventplanes)[6];
    const reco::EvtPlane &ephfm2 = (*eventplanes)[7];
    const reco::EvtPlane &ephfp3 = (*eventplanes)[13];
    const reco::EvtPlane &ephfm3 = (*eventplanes)[14];

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

  // best vertex
  bestvz = -999.9;
  bestvx = -999.9;
  bestvy = -999.9;
  double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
  const reco::Vertex &vtx = (*vertices)[0];
  bestvz = vtx.z();
  bestvx = vtx.x();
  bestvy = vtx.y();
  bestvzError = vtx.zError();
  bestvxError = vtx.xError();
  bestvyError = vtx.yError();

  // Ntrkoffline
  Ntrkoffline = 0;
  if (multMax_ != -1 && multMin_ != -1) {
    for (unsigned it = 0; it < tracks->size(); ++it) {

      const reco::Track &trk = (*tracks)[it];

      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzvtx = trk.dz(bestvtx);
      double dxyvtx = trk.dxy(bestvtx);
      double dzerror = sqrt(trk.dzError() * trk.dzError() + bestvzError * bestvzError);
      double dxyerror = sqrt(trk.d0Error() * trk.d0Error() + bestvxError * bestvyError);

      if (!trk.quality(reco::TrackBase::highPurity))
        continue;
      if (fabs(trk.ptError()) / trk.pt() > 0.10)
        continue;
      if (fabs(dzvtx / dzerror) > 3)
        continue;
      if (fabs(dxyvtx / dxyerror) > 3)
        continue;

      double eta = trk.eta();
      double pt = trk.pt();

      if (fabs(eta) > 2.4)
        continue;
      if (pt <= 0.4)
        continue;
      Ntrkoffline++;
    }
  }

  // Gen info for matching
  std::vector<reco::GenParticleRef> genRefs;
  if (doGenMatching_) {
    if (!genpars.isValid()) {
      cout << "Gen matching cannot be done without Gen collection!!" << endl;
      return;
    }
    for (unsigned int it = 0; it < genpars->size(); ++it) {
      const reco::GenParticle &trk = (*genpars)[it];
      int id = trk.pdgId();
      if (fabs(id) != PID_)
        continue; // check is target
      if (decayInGen_ && trk.numberOfDaughters() != 2 && !threeProngDecay_)
        continue; // check 2-pron decay if target decays in Gen
      if (decayInGen_ && trk.numberOfDaughters() != 3 && threeProngDecay_)
        continue; // check 2-pron decay if target decays in Gen

      // wrong when considering two layer decay
      int nDau = threeProngDecay_ ? 3 : 2;
      std::vector<unsigned int> idxs;
      std::vector<unsigned int> permutations(nDau);
      std::iota(permutations.begin(), permutations.end(), 0);
      std::sort(permutations.begin(), permutations.end());
      if (!threeProngDecay_) {
        do {
          auto Dd1 = trk.daughter(permutations.at(0));
          auto Dd2 = trk.daughter(permutations.at(1));
          if (abs(Dd1->pdgId()) == PID_dau1_ && abs(Dd2->pdgId()) == PID_dau2_) {
            idxs = permutations;
            break;
          }
        } while (std::next_permutation(permutations.begin(), permutations.end()));
      } else {
        do {
          auto Dd1 = trk.daughter(permutations.at(0));
          auto Dd2 = trk.daughter(permutations.at(1));
          auto Dd3 = trk.daughter(permutations.at(2));

          if (abs(Dd1->pdgId()) == PID_dau1_ && abs(Dd2->pdgId()) == PID_dau2_ && abs(Dd3->pdgId() == PID_dau3_)) {
            idxs = permutations;
            break;
          }
        } while (std::next_permutation(permutations.begin(), permutations.end()));
      }
      if (decayInGen_ && idxs.empty())
        continue;
      genRefs.push_back(reco::GenParticleRef(genpars, it));
    }
    // if (genRefs.size()>1) std::cout << "More than one target of generated
    // particles\n";
  }

  // RECO Candidate info
  candSize = v0candidates_->size();
  if (debug_)
    std::cout << "Reco cand Size : " << candSize << std::endl;
  for (unsigned it = 0; it < v0candidates_->size(); ++it) {

    const reco::VertexCompositeCandidate &trk = (*v0candidates_)[it];

    double secvz = -999.9, secvx = -999.9, secvy = -999.9;
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();

    eta[it] = trk.eta();
    y[it] = trk.rapidity();
    pt[it] = trk.pt();
    phi[it] = trk.phi();
    flavor[it] = trk.pdgId() / abs(trk.pdgId());

    mva[it] = 0.0;
    if (useAnyMVA_)
      mva[it] = (*mvavalues)[it];

    dca3D[it] = -1.0;
    dcaErr3D[it] = -1.0;
    if (useDCA_) {
      dca3D[it] = dcaValues->at(it);
      dcaErr3D[it] = dcaErrors->at(it);
    }

    double px = trk.px();
    double py = trk.py();
    double pz = trk.pz();
    mass[it] = trk.mass();

    const reco::Candidate *d1 = trk.daughter(0);
    const reco::Candidate *d2 = trk.daughter(1);
    if (doubleCand_) {

      flavor1[it] = d1->pdgId() / abs(d1->pdgId());
      flavor2[it] = d2->pdgId() / abs(d2->pdgId());
      if (debug_ && d1->pt() == d2->pt())
        std::cout << "Two daughter is same" << std::endl;
    }
    const reco::Candidate *d3 = 0;
    if (threeProngDecay_)
      d3 = trk.daughter(2);

    // Gen match
    if (doGenMatching_) {
      const auto nGen = genRefs.size();
      if (!doGenDoubleDecay_) {
        matchGEN[it] = false;
        isSwap[it] = false;
        idmom_reco1[it] = -77;
        idBAnc_reco1[it] = -77;

        for (unsigned int igen = 0; igen < nGen; igen++) {
          const auto &genRef = genRefs.at(igen);

          reco::Candidate const *recoDaus[3] = {nullptr, nullptr, nullptr};
          reco::Candidate const *genDaus[3] = {nullptr, nullptr, nullptr};

          const auto nGenDau = genRef->numberOfDaughters();
          std::vector<unsigned int> permutations(nGenDau);
          std::iota(permutations.begin(), permutations.end(), 0);
          std::sort(permutations.begin(), permutations.end());

          do {
            matchGEN[it] = false;
            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              genDaus[iDau] = genRef->daughter(permutations.at(iDau));
              recoDaus[iDau] = trk.daughter(iDau);
            }

            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              const double dR = reco::deltaR(genDaus[iDau]->eta(), genDaus[iDau]->phi(), recoDaus[iDau]->eta(),
                                             recoDaus[iDau]->phi());
              const double dPt = abs(genDaus[iDau]->pt() - recoDaus[iDau]->pt()) / recoDaus[iDau]->pt();
              const bool unMatchCharge = genDaus[iDau]->charge() != recoDaus[iDau]->charge();
              const bool unMatchDR = dR > deltaR_;
              const bool unMatchDPt = dPt > 0.5;
              matchGEN[it] = matchGEN[it] || unMatchCharge || unMatchDR || unMatchDPt;
              if (!matchGEN[it]) {
                isSwap[it] = (recoDaus[iDau]->pdgId() != genRef->pdgId());
                break;
              };
            }
            matchGEN[it] = !matchGEN[it];
          } while (std::next_permutation(permutations.begin(), permutations.end()));

          if (matchGEN[it]) {
            auto mom_ref = findMother(genRef);
            if (mom_ref.isNonnull())
              idmom_reco1[it] = mom_ref->pdgId();
            int __count_anc__ = 0;
            auto __ref_anc__ = mom_ref;
            while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
              __ref_anc__ = findMother(__ref_anc__);
              if (__ref_anc__.isNonnull()) {
                if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                  idBAnc_reco1[it] = __ref_anc__->pdgId();
                  break;
                }
              }
              __count_anc__++;
            }

            gen_D1pT_[it] = genRef->pt();
            gen_D1eta_[it] = genRef->eta();
            gen_D1phi_[it] = genRef->phi();
            gen_D1mass_[it] = genRef->mass();
            gen_D1y_[it] = genRef->rapidity();
            gen_D1charge_[it] = genRef->charge();
            gen_D1pdgId_[it] = genRef->pdgId();

            // all done in genDecayLength
            // gen_decayLength3D_;
            // gen_decayLength2D_;
            // gen_angle3D_;
            // gen_angle2D_;
            genDecayLength(*genRef, gen_D1decayLength2D_[it], gen_D1decayLength3D_[it], gen_D1angle2D_[it],
                           gen_D1angle3D_[it]);
            getAncestorId(*genRef, gen_D1ancestorId_[it], gen_D1ancestorFlavor_[it]);

            gen_D1pTD1_[it] = genDaus[0]->pt();
            gen_D1etaD1_[it] = genDaus[0]->eta();
            gen_D1phiD1_[it] = genDaus[0]->phi();
            gen_D1massD1_[it] = genDaus[0]->mass();
            gen_D1yD1_[it] = genDaus[0]->rapidity();
            gen_D1chargeD1_[it] = genDaus[0]->charge();
            gen_D1pdgIdD1_[it] = genDaus[0]->pdgId();

            gen_D1pTD2_[it] = genDaus[1]->pt();
            gen_D1etaD2_[it] = genDaus[1]->eta();
            gen_D1phiD2_[it] = genDaus[1]->phi();
            gen_D1massD2_[it] = genDaus[1]->mass();
            gen_D1yD2_[it] = genDaus[1]->rapidity();
            gen_D1chargeD2_[it] = genDaus[1]->charge();
            gen_D1pdgIdD2_[it] = genDaus[1]->pdgId();
            break;
          }
        }
        if (!matchGEN[it]) {
          gen_D1pT_[it] = -99;
          gen_D1eta_[it] = -99;
          gen_D1phi_[it] = -99;
          gen_D1mass_[it] = -99;
          gen_D1y_[it] = -99;
          gen_D1decayLength3D_[it] = -99;
          gen_D1decayLength2D_[it] = -99;
          gen_D1angle3D_[it] = -99;
          gen_D1angle2D_[it] = -99;
          gen_D1pTD1_[it] = -99;
          gen_D1etaD1_[it] = -99;
          gen_D1phiD1_[it] = -99;
          gen_D1massD1_[it] = -99;
          gen_D1yD1_[it] = -99;
          gen_D1pTD2_[it] = -99;
          gen_D1etaD2_[it] = -99;
          gen_D1phiD2_[it] = -99;
          gen_D1massD2_[it] = -99;
          gen_D1yD2_[it] = -99;
        }
      }
      if (doGenDoubleDecay_) {
        matchGEN1[it] = false;
        matchGEN2[it] = false;
        isSwap1[it] = false;
        isSwap2[it] = false;
        idmom_reco1[it] = -77;
        idmom_reco2[it] = -77;
        idBAnc_reco1[it] = -77;
        idBAnc_reco2[it] = -77;
        matchToGen1[it] = MAXCAN + 1;
        matchToGen2[it] = MAXCAN + 1;
        if (debug_)
          std::cout << "nGen : " << nGen << std::endl;
        matchGEN[it] = false;
        for (unsigned int igen = 0; igen < nGen; igen++) {
          auto const &theGen = genRefs.at(igen);
          // Only works for 2 body two layer decay
          reco::Candidate const *recoDaus1[2] = {nullptr, nullptr};
          reco::Candidate const *recoDaus2[2] = {nullptr, nullptr};

          reco::Candidate const *genDaus[2] = {nullptr, nullptr};

          const auto nGenDau = theGen->numberOfDaughters();
          if (debug_)
            std::cout << "nGenDau: " << nGenDau << std::endl;
          std::vector<unsigned int> perm = {0, 1};
          do {
            bool _matchGEN_ = false;
            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              genDaus[iDau] = theGen->daughter(perm.at(iDau));
              recoDaus1[iDau] = trk.daughter(0)->daughter(iDau);
            }
            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              const double dR = reco::deltaR(genDaus[iDau]->eta(), genDaus[iDau]->phi(), recoDaus1[iDau]->eta(),
                                             recoDaus1[iDau]->phi());
              const double dPt = abs(genDaus[iDau]->pt() - recoDaus1[iDau]->pt()) / recoDaus1[iDau]->pt();
              const bool unMatchCharge = genDaus[iDau]->charge() != recoDaus1[iDau]->charge();
              const bool unMatchDR = dR > deltaR_;
              const bool unMatchDPt = dPt > 0.5;
              _matchGEN_ = !(unMatchCharge || unMatchDR || unMatchDPt);
#ifdef DEBUG2
              cout << Form("gen %d etas: (%.4f, %.4f), phis: (%.4f, %.4f), "
                           "match : %d",
                           igen, recoDaus1[iDau]->eta(), genDaus[iDau]->eta(), recoDaus1[iDau]->phi(),
                           recoDaus1[iDau]->phi(), _matchGEN_)
                   << endl;
#endif
              if (_matchGEN_) {
                // isSwap1[it] =
                // ((recoDaus1[iDau]->pdgId()/abs(recoDaus1[iDau]->pdgId())) !=
                // (theGen->pdgId()/abs(theGen->pdgId())));
                isSwap1[it] = ((trk.daughter(0)->pdgId()) != (theGen->pdgId()));
                if (matchToGen1[it] < MAXCAN + 1)
                  std::cout << "Double matching! occurred with igen " << matchToGen1[it] << " and " << igen
                            << std::endl;
                matchToGen1[it] = igen;
                break;
              }
#ifdef DEBUG2
              cout << Form("it %d match : %d, swap: %d", it, _matchGEN_, isSwap1[it]) << nGen << endl;
#endif
            }
            if (matchGEN1[it])
              break;
            matchGEN1[it] = _matchGEN_;
            if (matchGEN1[it]) {
              auto mom_ref = findMother(theGen);
              if (mom_ref.isNonnull())
                idmom_reco1[it] = mom_ref->pdgId();
              int __count_anc__ = 0;
              auto __ref_anc__ = mom_ref;
              while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
                __ref_anc__ = findMother(__ref_anc__);
                if (__ref_anc__.isNonnull()) {
                  if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                    idBAnc_reco1[it] = __ref_anc__->pdgId();
                    break;
                  }
                }
                __count_anc__++;
              }

              gen_D1pT_[it] = theGen->pt();
              gen_D1eta_[it] = theGen->eta();
              gen_D1phi_[it] = theGen->phi();
              gen_D1mass_[it] = theGen->mass();
              gen_D1y_[it] = theGen->rapidity();
              gen_D1charge_[it] = theGen->charge();
              gen_D1pdgId_[it] = theGen->pdgId();

              // all done in genDecayLength
              // gen_decayLength3D_;
              // gen_decayLength2D_;
              // gen_angle3D_;
              // gen_angle2D_;
              genDecayLength(*theGen, gen_D1decayLength2D_[it], gen_D1decayLength3D_[it], gen_D1angle2D_[it],
                             gen_D1angle3D_[it]);
              getAncestorId(*theGen, gen_D1ancestorId_[it], gen_D1ancestorFlavor_[it]);

              gen_D1pTD1_[it] = genDaus[0]->pt();
              gen_D1etaD1_[it] = genDaus[0]->eta();
              gen_D1phiD1_[it] = genDaus[0]->phi();
              gen_D1massD1_[it] = genDaus[0]->mass();
              gen_D1yD1_[it] = genDaus[0]->rapidity();
              gen_D1chargeD1_[it] = genDaus[0]->charge();
              gen_D1pdgIdD1_[it] = genDaus[0]->pdgId();

              gen_D1pTD2_[it] = genDaus[1]->pt();
              gen_D1etaD2_[it] = genDaus[1]->eta();
              gen_D1phiD2_[it] = genDaus[1]->phi();
              gen_D1massD2_[it] = genDaus[1]->mass();
              gen_D1yD2_[it] = genDaus[1]->rapidity();
              gen_D1chargeD2_[it] = genDaus[1]->charge();
              gen_D1pdgIdD2_[it] = genDaus[1]->pdgId();
              break;
            }
          } while (std::next_permutation(perm.begin(), perm.end()));
          std::sort(perm.begin(), perm.end());
          do {
            bool _matchGEN_ = false;
            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              genDaus[iDau] = theGen->daughter(perm.at(iDau));
              recoDaus2[iDau] = trk.daughter(1)->daughter(iDau);
            }
            for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
              const double dR = reco::deltaR(genDaus[iDau]->eta(), genDaus[iDau]->phi(), recoDaus2[iDau]->eta(),
                                             recoDaus2[iDau]->phi());
              const double dPt = abs(genDaus[iDau]->pt() - recoDaus2[iDau]->pt()) / recoDaus2[iDau]->pt();
              const bool unMatchCharge = genDaus[iDau]->charge() != recoDaus2[iDau]->charge();
              const bool unMatchDR = dR > deltaR_;
              const bool unMatchDPt = dPt > 0.5;
              _matchGEN_ = !(unMatchCharge || unMatchDR || unMatchDPt);
              if (_matchGEN_) {
                isSwap2[it] = ((trk.daughter(1)->pdgId()) != (theGen->pdgId()));
                // isSwap2[it] =
                // ((recoDaus2[iDau]->pdgId()/abs(recoDaus2[iDau]->pdgId())) !=
                // (theGen->pdgId()/abs(theGen->pdgId())));
                if (matchToGen2[it] < MAXCAN + 1)
                  std::cout << "Double matching! occurred with igen " << matchToGen2[it] << " and " << igen
                            << std::endl;
                matchToGen2[it] = igen;
                break;
              }
            }
            if (matchGEN2[it])
              break;
            matchGEN2[it] = _matchGEN_;
            if (matchGEN2[it]) {
              auto mom_ref = findMother(theGen);
              if (mom_ref.isNonnull())
                idmom_reco2[it] = mom_ref->pdgId();
              int __count_anc__ = 0;
              auto __ref_anc__ = mom_ref;
              while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
                __ref_anc__ = findMother(__ref_anc__);
                if (__ref_anc__.isNonnull()) {
                  if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                    idBAnc_reco2[it] = __ref_anc__->pdgId();
                    break;
                  }
                }
                __count_anc__++;
              }

              gen_D2pT_[it] = theGen->pt();
              gen_D2eta_[it] = theGen->eta();
              gen_D2phi_[it] = theGen->phi();
              gen_D2mass_[it] = theGen->mass();
              gen_D2y_[it] = theGen->rapidity();
              gen_D2charge_[it] = theGen->charge();
              gen_D2pdgId_[it] = theGen->pdgId();

              // all done in genDecayLength
              // gen_decayLength3D_;
              // gen_decayLength2D_;
              // gen_angle3D_;
              // gen_angle2D_;
              genDecayLength(*theGen, gen_D2decayLength2D_[it], gen_D2decayLength3D_[it], gen_D2angle2D_[it],
                             gen_D2angle3D_[it]);
              getAncestorId(*theGen, gen_D2ancestorId_[it], gen_D2ancestorFlavor_[it]);

              gen_D2pTD1_[it] = genDaus[0]->pt();
              gen_D2etaD1_[it] = genDaus[0]->eta();
              gen_D2phiD1_[it] = genDaus[0]->phi();
              gen_D2massD1_[it] = genDaus[0]->mass();
              gen_D2yD1_[it] = genDaus[0]->rapidity();
              gen_D2chargeD1_[it] = genDaus[0]->charge();
              gen_D2pdgIdD1_[it] = genDaus[0]->pdgId();

              gen_D2pTD2_[it] = genDaus[1]->pt();
              gen_D2etaD2_[it] = genDaus[1]->eta();
              gen_D2phiD2_[it] = genDaus[1]->phi();
              gen_D2massD2_[it] = genDaus[1]->mass();
              gen_D2yD2_[it] = genDaus[1]->rapidity();
              gen_D2chargeD2_[it] = genDaus[1]->charge();
              gen_D2pdgIdD2_[it] = genDaus[1]->pdgId();
              break;
            }
          } while (std::next_permutation(perm.begin(), perm.end()));

        } // END for nGen
        if (!matchGEN1[it]) {
          gen_D1pT_[it] = -99;
          gen_D1eta_[it] = -99;
          gen_D1phi_[it] = -99;
          gen_D1mass_[it] = -99;
          gen_D1y_[it] = -99;
          gen_D1decayLength3D_[it] = -99;
          gen_D1decayLength2D_[it] = -99;
          gen_D1angle3D_[it] = -99;
          gen_D1angle2D_[it] = -99;
          gen_D1pTD1_[it] = -99;
          gen_D1etaD1_[it] = -99;
          gen_D1phiD1_[it] = -99;
          gen_D1massD1_[it] = -99;
          gen_D1yD1_[it] = -99;
          gen_D1pTD2_[it] = -99;
          gen_D1etaD2_[it] = -99;
          gen_D1phiD2_[it] = -99;
          gen_D1massD2_[it] = -99;
          gen_D1yD2_[it] = -99;
        }
        if (!matchGEN2[it]) {
          gen_D2pT_[it] = -99;
          gen_D2eta_[it] = -99;
          gen_D2phi_[it] = -99;
          gen_D2mass_[it] = -99;
          gen_D2y_[it] = -99;
          gen_D2decayLength3D_[it] = -99;
          gen_D2decayLength2D_[it] = -99;
          gen_D2angle3D_[it] = -99;
          gen_D2angle2D_[it] = -99;
          gen_D2pTD1_[it] = -99;
          gen_D2etaD1_[it] = -99;
          gen_D2phiD1_[it] = -99;
          gen_D2massD1_[it] = -99;
          gen_D2yD1_[it] = -99;
          gen_D2pTD2_[it] = -99;
          gen_D2etaD2_[it] = -99;
          gen_D2phiD2_[it] = -99;
          gen_D2massD2_[it] = -99;
          gen_D2yD2_[it] = -99;
        }
      } // END if doGenDoubleDecay_
      matchGEN[it] = (matchGEN1[it] && matchGEN2[it]);
    }

    double pxd1 = d1->px();
    double pyd1 = d1->py();
    double pzd1 = d1->pz();
    double pxd2 = d2->px();
    double pyd2 = d2->py();
    double pzd2 = d2->pz();

    TVector3 dauvec1(pxd1, pyd1, pzd1);
    TVector3 dauvec2(pxd2, pyd2, pzd2);

    // pt
    pt1[it] = d1->pt();
    pt2[it] = d2->pt();

    // momentum
    p1[it] = d1->p();
    p2[it] = d2->p();

    // eta
    eta1[it] = d1->eta();
    eta2[it] = d2->eta();

    // phi
    phi1[it] = d1->phi();
    phi2[it] = d2->phi();

    // charge
    charge1[it] = d1->charge();
    charge2[it] = d2->charge();

    double pxd3 = -999.9;
    double pyd3 = -999.9;
    double pzd3 = -999.9;
    if (threeProngDecay_ && d3) {
      pxd3 = d3->px();
      pyd3 = d3->py();
      pzd3 = d3->pz();
      pt3[it] = d3->pt();
      p3[it] = d3->p();
      eta3[it] = d3->eta();
      phi3[it] = d3->phi();
      charge3[it] = d3->charge();
    }
    TVector3 dauvec3(pxd3, pyd3, pzd3);

    pid1[it] = -99999;
    pid2[it] = -99999;
    if (doGenMatchingTOF_) {
      for (unsigned it = 0; it < genpars->size(); ++it) {

        const reco::GenParticle &trk = (*genpars)[it];

        if (trk.pt() < 0.001)
          continue;

        int id = trk.pdgId();
        TVector3 trkvect(trk.px(), trk.py(), trk.pz());

        if (fabs(id) != PID_ && trk.charge()) {
          // matching daughter 1
          double deltaR = trkvect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt1[it]) / pt1[it]) < 0.5 && trk.charge() == charge1[it] &&
              pid1[it] == -99999) {
            pid1[it] = id;
          }

          // matching daughter 2
          deltaR = trkvect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt2[it]) / pt2[it]) < 0.5 && trk.charge() == charge2[it] &&
              pid2[it] == -99999) {
            pid2[it] = id;
          }
        }

        if (fabs(id) == PID_ && trk.numberOfDaughters() == 2) {
          const reco::Candidate *Dd1 = trk.daughter(0);
          const reco::Candidate *Dd2 = trk.daughter(1);
          TVector3 d1vect(Dd1->px(), Dd1->py(), Dd1->pz());
          TVector3 d2vect(Dd2->px(), Dd2->py(), Dd2->pz());
          int id1 = Dd1->pdgId();
          int id2 = Dd2->pdgId();

          double deltaR = d1vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt1[it]) / pt1[it]) < 0.5 && Dd1->charge() == charge1[it] &&
              pid1[it] == -99999) {
            pid1[it] = id1;
          }
          deltaR = d2vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt1[it]) / pt1[it]) < 0.5 && Dd2->charge() == charge1[it] &&
              pid1[it] == -99999) {
            pid1[it] = id1;
          }

          deltaR = d1vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt2[it]) / pt2[it]) < 0.5 && Dd1->charge() == charge2[it] &&
              pid2[it] == -99999) {
            pid2[it] = id2;
          }
          deltaR = d2vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt2[it]) / pt2[it]) < 0.5 && Dd2->charge() == charge2[it] &&
              pid2[it] == -99999) {
            pid2[it] = id2;
          }
        }

        if (pid1[it] != -99999 && pid2[it] != -99999)
          break;
      }
    }

    // vtxChi2
    vtxChi2[it] = trk.vertexChi2();
    ndf[it] = trk.vertexNdof();
    VtxProb[it] = TMath::Prob(vtxChi2[it], ndf[it]);

    // PAngle
    TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
    TVector3 secvec(px, py, pz);

    TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
    TVector3 secvec2D(px, py, 0);

    agl[it] = cos(secvec.Angle(ptosvec));
    agl_abs[it] = secvec.Angle(ptosvec);

    agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

    // Decay length 3D
    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;
    typedef ROOT::Math::SVector<double, 6> SVector6;

    SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
    SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

    dl[it] = ROOT::Math::Mag(distanceVector);
    dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl[it];

    dlos[it] = dl[it] / dlerror[it];

    // correct way for both DCA and its Error
    // std::cout << "By cur3DIP                " << dca3D[it] << " +/- " <<
    // dcaErr3D[it] <<"\n"; incorrect way for DCA error std::cout << "By decay
    // length and alpha " << std::sin(agl_abs[it])*dl[it]<< " +/- " <<
    // dlerror[it]* std::sin(agl_abs[it]) <<"\n"; std::cout << "\n";

    // Decay length 2D
    SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
    SVector6 v2(trk.vertexCovariance(0, 0), trk.vertexCovariance(0, 1), trk.vertexCovariance(1, 1), 0, 0, 0);

    SMatrixSym3D sv1(v1);
    SMatrixSym3D sv2(v2);

    SMatrixSym3D totalCov2D = sv1 + sv2;
    SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

    dl2D[it] = ROOT::Math::Mag(distanceVector2D);
    double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / dl2D[it];

    dlos2D[it] = dl2D[it] / dl2Derror;

    // trk info
    auto dau1 = d1->get<reco::TrackRef>();
    if (!twoLayerDecay_) {
      // trk quality
      trkquality1[it] = dau1->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx1[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx1[it] = dEdxTrack[dau1].dEdx();
      }

      T4dedx1[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx1[it] = dEdxTrack[dau1].dEdx();
      }

      // track Chi2
      trkChi1[it] = dau1->normalizedChi2();

      // track pT error
      ptErr1[it] = dau1->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit1[it] = dau1->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest1 = dau1->dz(bestvtx);
      double dxybest1 = dau1->dxy(bestvtx);
      double dzerror1 = sqrt(dau1->dzError() * dau1->dzError() + bestvzError * bestvzError);
      double dxyerror1 = sqrt(dau1->d0Error() * dau1->d0Error() + bestvxError * bestvyError);

      dzos1[it] = dzbest1 / dzerror1;
      dxyos1[it] = dxybest1 / dxyerror1;
    }
    if (!doubleCand_) {
      auto dau2 = d2->get<reco::TrackRef>();

      // trk quality
      trkquality2[it] = dau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx2[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx2[it] = dEdxTrack[dau2].dEdx();
      }

      T4dedx2[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx2[it] = dEdxTrack[dau2].dEdx();
      }

      // track Chi2
      trkChi2[it] = dau2->normalizedChi2();

      // track pT error
      ptErr2[it] = dau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit2[it] = dau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest2 = dau2->dz(bestvtx);
      double dxybest2 = dau2->dxy(bestvtx);
      double dzerror2 = sqrt(dau2->dzError() * dau2->dzError() + bestvzError * bestvzError);
      double dxyerror2 = sqrt(dau2->d0Error() * dau2->d0Error() + bestvxError * bestvyError);

      dzos2[it] = dzbest2 / dzerror2;
      dxyos2[it] = dxybest2 / dxyerror2;

      if (doMuon_) {
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

        const int muId1 = muAssocToTrack(dau1, theMuonHandle);
        const int muId2 = muAssocToTrack(dau2, theMuonHandle);

        if (muId1 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId1];

          onestmuon1[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
          pfmuon1[it] = cand.isPFMuon();
          glbmuon1[it] = cand.isGlobalMuon();
          trkmuon1[it] = cand.isTrackerMuon();
          calomuon1[it] = cand.isCaloMuon();

          if (glbmuon1[it] && trkmuon1[it] && cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
              cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
              fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 && fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
            softmuon1[it] = true;
        }

        if (muId2 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId2];

          onestmuon2[it] = muon::isGoodMuon(cand, muon::selectionTypeFromString("TMOneStationTight"));
          pfmuon2[it] = cand.isPFMuon();
          glbmuon2[it] = cand.isGlobalMuon();
          trkmuon2[it] = cand.isTrackerMuon();
          calomuon2[it] = cand.isCaloMuon();

          if (glbmuon2[it] && trkmuon2[it] && cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
              cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
              fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 && fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
            softmuon2[it] = true;
        }

        if (doMuonFull_) {

          if (muId1 != -1) {
            const reco::Muon &cand = (*theMuonHandle)[muId1];

            nmatchedch1[it] = cand.numberOfMatches();
            nmatchedst1[it] = cand.numberOfMatchedStations();

            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy1[it] = muenergy.hadMax;

            const std::vector<reco::MuonChamberMatch> &muchmatches = cand.matches();

            for (unsigned int ich = 0; ich < muchmatches.size(); ich++) {
              x_exp = muchmatches[ich].x;
              y_exp = muchmatches[ich].y;
              xerr_exp = muchmatches[ich].xErr;
              yerr_exp = muchmatches[ich].yErr;
              dxdz_exp = muchmatches[ich].dXdZ;
              dydz_exp = muchmatches[ich].dYdZ;
              dxdzerr_exp = muchmatches[ich].dXdZErr;
              dydzerr_exp = muchmatches[ich].dYdZErr;

              std::vector<reco::MuonSegmentMatch> musegmatches = muchmatches[ich].segmentMatches;

              if (!musegmatches.size())
                continue;
              for (unsigned int jseg = 0; jseg < musegmatches.size(); jseg++) {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;

                if (sqrt((x_seg - x_exp) * (x_seg - x_exp) + (y_seg - y_exp) * (y_seg - y_exp)) <
                    sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                  dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }

              dx1_seg_[it] = dx_seg;
              dy1_seg_[it] = dy_seg;
              dxSig1_seg_[it] = dxSig_seg;
              dySig1_seg_[it] = dySig_seg;
              ddxdz1_seg_[it] = ddxdz_seg;
              ddydz1_seg_[it] = ddydz_seg;
              ddxdzSig1_seg_[it] = ddxdzSig_seg;
              ddydzSig1_seg_[it] = ddydzSig_seg;
            }
          }

          if (muId2 != -1) {
            const reco::Muon &cand = (*theMuonHandle)[muId2];

            nmatchedch2[it] = cand.numberOfMatches();
            nmatchedst2[it] = cand.numberOfMatchedStations();

            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy2[it] = muenergy.hadMax;

            const std::vector<reco::MuonChamberMatch> &muchmatches = cand.matches();
            for (unsigned int ich = 0; ich < muchmatches.size(); ich++)
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

              if (!musegmatches.size())
                continue;
              for (unsigned int jseg = 0; jseg < musegmatches.size(); jseg++) {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;

                if (sqrt((x_seg - x_exp) * (x_seg - x_exp) + (y_seg - y_exp) * (y_seg - y_exp)) <
                    sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                  dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }

              dx2_seg_[it] = dx_seg;
              dy2_seg_[it] = dy_seg;
              dxSig2_seg_[it] = dxSig_seg;
              dySig2_seg_[it] = dySig_seg;
              ddxdz2_seg_[it] = ddxdz_seg;
              ddydz2_seg_[it] = ddydz_seg;
              ddxdzSig2_seg_[it] = ddxdzSig_seg;
              ddydzSig2_seg_[it] = ddydzSig_seg;
            }
          }
        } // doMuonFull
      }
    }

    if (twoLayerDecay_) {
      grand_mass[it] = d1->mass();
      mva1[it] = (*mvavalues)[it];

      const reco::Candidate *gd1 = d1->daughter(0);
      const reco::Candidate *gd2 = d1->daughter(1);

      double gpxd1 = gd1->px();
      double gpyd1 = gd1->py();
      double gpzd1 = gd1->pz();
      double gpxd2 = gd2->px();
      double gpyd2 = gd2->py();
      double gpzd2 = gd2->pz();

      TVector3 gdauvec1(gpxd1, gpyd1, gpzd1);
      TVector3 gdauvec2(gpxd2, gpyd2, gpzd2);

      auto gdau1 = gd1->get<reco::TrackRef>();
      auto gdau2 = gd2->get<reco::TrackRef>();

      // trk quality

      grand_trkquality1[it] = gdau1->quality(reco::TrackBase::highPurity);
      grand_trkquality2[it] = gdau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      grand_H2dedx1[it] = -999.9;
      grand_H2dedx2[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        grand_H2dedx1[it] = dEdxTrack[gdau1].dEdx();
        grand_H2dedx2[it] = dEdxTrack[gdau2].dEdx();
      }

      grand_T4dedx1[it] = -999.9;
      grand_T4dedx2[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        grand_T4dedx1[it] = dEdxTrack[gdau1].dEdx();
        grand_T4dedx2[it] = dEdxTrack[gdau2].dEdx();
      }

      // track pt
      grand_pt1[it] = gd1->pt();
      grand_pt2[it] = gd2->pt();

      // track momentum
      grand_p1[it] = gd1->p();
      grand_p2[it] = gd2->p();

      // track eta
      grand_eta1[it] = gd1->eta();
      grand_eta2[it] = gd2->eta();

      grand_phi1[it] = gd1->phi();
      grand_phi2[it] = gd2->phi();

      // track charge
      grand_charge1[it] = gd1->charge();
      grand_charge2[it] = gd2->charge();

      // track Chi2
      grand_trkChi1[it] = gdau1->normalizedChi2();
      grand_trkChi2[it] = gdau2->normalizedChi2();

      // track pT error
      grand_ptErr1[it] = gdau1->ptError();
      grand_ptErr2[it] = gdau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = d1->vz();
      secvx = d1->vx();
      secvy = d1->vy();

      // trkNHits
      grand_nhit1[it] = gdau1->numberOfValidHits();
      grand_nhit2[it] = gdau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double gdzbest1 = gdau1->dz(bestvtx);
      double gdxybest1 = gdau1->dxy(bestvtx);
      double gdzerror1 = sqrt(gdau1->dzError() * gdau1->dzError() + bestvzError * bestvzError);
      double gdxyerror1 = sqrt(gdau1->d0Error() * gdau1->d0Error() + bestvxError * bestvyError);

      grand_dzos1[it] = gdzbest1 / gdzerror1;
      grand_dxyos1[it] = gdxybest1 / gdxyerror1;

      double gdzbest2 = gdau2->dz(bestvtx);
      double gdxybest2 = gdau2->dxy(bestvtx);
      double gdzerror2 = sqrt(gdau2->dzError() * gdau2->dzError() + bestvzError * bestvzError);
      double gdxyerror2 = sqrt(gdau2->d0Error() * gdau2->d0Error() + bestvxError * bestvyError);

      grand_dzos2[it] = gdzbest2 / gdzerror2;
      grand_dxyos2[it] = gdxybest2 / gdxyerror2;

      // vtxChi2
      grand_vtxChi2[it] = d1->vertexChi2();
      grand_ndf[it] = d1->vertexNdof();
      grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it], grand_ndf[it]);

      // PAngle
      TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
      TVector3 secvec(d1->px(), d1->py(), d1->pz());

      TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
      TVector3 secvec2D(d1->px(), d1->py(), 0);

      grand_agl[it] = cos(secvec.Angle(ptosvec));
      grand_agl_abs[it] = secvec.Angle(ptosvec);

      grand_agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

      // Decay length 3D
      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;
      typedef ROOT::Math::SVector<double, 6> SVector6;

      SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
      SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

      grand_dl[it] = ROOT::Math::Mag(distanceVector);
      grand_dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / grand_dl[it];

      grand_dlos[it] = grand_dl[it] / grand_dlerror[it];

      // Decay length 2D
      SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
      SVector6 v2(d1->vertexCovariance(0, 0), d1->vertexCovariance(0, 1), d1->vertexCovariance(1, 1), 0, 0, 0);

      SMatrixSym3D sv1(v1);
      SMatrixSym3D sv2(v2);

      SMatrixSym3D totalCov2D = sv1 + sv2;
      SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

      double gdl2D = ROOT::Math::Mag(distanceVector2D);
      double gdl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / gdl2D;

      grand_dlos2D[it] = gdl2D / gdl2Derror;
      if (doubleCand_) {
        grand_mass2[it] = d2->mass();
        mva2[it] = (*mvavalues2)[it];

        const reco::Candidate *gd21 = d2->daughter(0);
        const reco::Candidate *gd22 = d2->daughter(1);

        double gpxd21 = gd21->px();
        double gpyd21 = gd21->py();
        double gpzd21 = gd21->pz();
        double gpxd22 = gd22->px();
        double gpyd22 = gd22->py();
        double gpzd22 = gd22->pz();

        TVector3 gdauvec21(gpxd21, gpyd21, gpzd21);
        TVector3 gdauvec22(gpxd22, gpyd22, gpzd22);

        auto gdau21 = gd21->get<reco::TrackRef>();
        auto gdau22 = gd22->get<reco::TrackRef>();

        // trk quality

        grand_trkquality21[it] = gdau21->quality(reco::TrackBase::highPurity);
        grand_trkquality22[it] = gdau22->quality(reco::TrackBase::highPurity);

        // trk dEdx
        grand_H2dedx21[it] = -999.9;
        grand_H2dedx22[it] = -999.9;

        if (dEdxHandle1.isValid()) {
          const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
          grand_H2dedx21[it] = dEdxTrack[gdau21].dEdx();
          grand_H2dedx22[it] = dEdxTrack[gdau22].dEdx();
        }

        grand_T4dedx21[it] = -999.9;
        grand_T4dedx22[it] = -999.9;

        if (dEdxHandle2.isValid()) {
          const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
          grand_T4dedx21[it] = dEdxTrack[gdau21].dEdx();
          grand_T4dedx22[it] = dEdxTrack[gdau22].dEdx();
        }

        // track pt
        grand_pt21[it] = gd21->pt();
        grand_pt22[it] = gd22->pt();

        // track momentum
        grand_p21[it] = gd21->p();
        grand_p22[it] = gd22->p();

        // track eta
        grand_eta21[it] = gd21->eta();
        grand_eta22[it] = gd22->eta();

        grand_phi21[it] = gd21->phi();
        grand_phi22[it] = gd22->phi();

        // track charge
        grand_charge21[it] = gd21->charge();
        grand_charge22[it] = gd22->charge();

        // track Chi2
        grand_trkChi21[it] = gdau21->normalizedChi2();
        grand_trkChi22[it] = gdau22->normalizedChi2();

        // track pT error
        grand_ptErr21[it] = gdau21->ptError();
        grand_ptErr22[it] = gdau22->ptError();

        // vertexCovariance 00-xError 11-y 22-z
        secvz = d2->vz();
        secvx = d2->vx();
        secvy = d2->vy();

        // trkNHits
        grand_nhit21[it] = gdau21->numberOfValidHits();
        grand_nhit22[it] = gdau22->numberOfValidHits();

        // DCA
        //  math::XYZPoint bestvtx2(bestvx,bestvy,bestvz);

        double gdzbest21 = gdau21->dz(bestvtx);
        double gdxybest21 = gdau21->dxy(bestvtx);
        double gdzerror21 = sqrt(gdau21->dzError() * gdau21->dzError() + bestvzError * bestvzError);
        double gdxyerror21 = sqrt(gdau21->d0Error() * gdau21->d0Error() + bestvxError * bestvyError);

        grand_dzos21[it] = gdzbest21 / gdzerror21;
        grand_dxyos21[it] = gdxybest21 / gdxyerror21;

        double gdzbest22 = gdau22->dz(bestvtx);
        double gdxybest22 = gdau22->dxy(bestvtx);
        double gdzerror22 = sqrt(gdau22->dzError() * gdau22->dzError() + bestvzError * bestvzError);
        double gdxyerror22 = sqrt(gdau22->d0Error() * gdau22->d0Error() + bestvxError * bestvyError);

        grand_dzos22[it] = gdzbest22 / gdzerror22;
        grand_dxyos22[it] = gdxybest22 / gdxyerror22;

        // vtxChi2
        grand_vtxChi22[it] = d2->vertexChi2();
        grand_ndf2[it] = d2->vertexNdof();
        grand_VtxProb2[it] = TMath::Prob(grand_vtxChi22[it], grand_ndf2[it]);

        // PAngle
        TVector3 ptosvec2(secvx - bestvx, secvy - bestvy, secvz - bestvz);
        TVector3 secvec2(d2->px(), d2->py(), d2->pz());

        TVector3 ptosvec2D2(secvx - bestvx, secvy - bestvy, 0);
        TVector3 secvec2D2(d2->px(), d2->py(), 0);

        grand_agl2[it] = cos(secvec2.Angle(ptosvec2));
        grand_agl_abs2[it] = secvec2.Angle(ptosvec2);

        grand_agl2D2[it] = cos(secvec2D2.Angle(ptosvec2D2));
        grand_agl2D_abs2[it] = secvec2D2.Angle(ptosvec2D2);

        // Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;

        SMatrixSym3D totalCov2 = vtx.covariance() + d2->vertexCovariance();
        SVector3 distanceVector2(secvx - bestvx, secvy - bestvy, secvz - bestvz);

        grand_dl2[it] = ROOT::Math::Mag(distanceVector2);
        grand_dlerror2[it] = sqrt(ROOT::Math::Similarity(totalCov2, distanceVector2)) / grand_dl2[it];

        grand_dlos2[it] = grand_dl2[it] / grand_dlerror2[it];

        // Decay length 2D
        SVector6 v21(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
        SVector6 v22(d2->vertexCovariance(0, 0), d2->vertexCovariance(0, 1), d2->vertexCovariance(1, 1), 0, 0, 0);

        SMatrixSym3D sv21(v1);
        SMatrixSym3D sv22(v2);

        SMatrixSym3D totalCov2D2 = sv21 + sv22;
        SVector3 distanceVector2D2(secvx - bestvx, secvy - bestvy, 0);

        double gdl2D2 = ROOT::Math::Mag(distanceVector2D2);
        double gdl2Derror2 = sqrt(ROOT::Math::Similarity(totalCov2D2, distanceVector2D2)) / gdl2D2;

        grand_dlos2D2[it] = gdl2D2 / gdl2Derror2;
      }
    }

    if (saveHistogram_) {
      for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++)
        for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
          if (pt[it] < pTBins_[ipt + 1] && pt[it] > pTBins_[ipt] && y[it] < yBins_[iy + 1] && y[it] > yBins_[iy]) {
            hMassVsMVA[iy][ipt]->Fill(mva[it], mass[it]);
            //                h3DDCAVsMVA[iy][ipt]->Fill(mva[it],dl[it]*sin(agl_abs[it]));
            //                h2DDCAVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]*sin(agl2D_abs[it]));

            if (saveAllHistogram_) {
              hpTVsMVA[iy][ipt]->Fill(mva[it], pt[it]);
              hetaVsMVA[iy][ipt]->Fill(mva[it], eta[it]);
              hyVsMVA[iy][ipt]->Fill(mva[it], y[it]);
              hVtxProbVsMVA[iy][ipt]->Fill(mva[it], VtxProb[it]);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl[it]);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl_abs[it]);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D[it]);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D_abs[it]);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos[it]);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl[it]);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos2D[it]);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl2D[it]);
              hzDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it], dzos1[it]);
              hxyDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it], dxyos1[it]);
              hNHitD1VsMVA[iy][ipt]->Fill(mva[it], nhit1[it]);
              hpTD1VsMVA[iy][ipt]->Fill(mva[it], pt1[it]);
              hpTerrD1VsMVA[iy][ipt]->Fill(mva[it], ptErr1[it] / pt1[it]);
              hEtaD1VsMVA[iy][ipt]->Fill(mva[it], eta1[it]);
              hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva[it], H2dedx1[it]);
              hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1[it], H2dedx1[it]);
              hzDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it], dzos2[it]);
              hxyDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it], dxyos2[it]);
              hNHitD2VsMVA[iy][ipt]->Fill(mva[it], nhit2[it]);
              hpTD2VsMVA[iy][ipt]->Fill(mva[it], pt2[it]);
              hpTerrD2VsMVA[iy][ipt]->Fill(mva[it], ptErr2[it] / pt2[it]);
              hEtaD2VsMVA[iy][ipt]->Fill(mva[it], eta2[it]);
              hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva[it], H2dedx2[it]);
              hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2[it], H2dedx2[it]);
              if (threeProngDecay_) {
                hzDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it], dzos3[it]);
                hxyDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it], dxyos3[it]);
                hNHitD3VsMVA[iy][ipt]->Fill(mva[it], nhit3[it]);
                hpTD3VsMVA[iy][ipt]->Fill(mva[it], pt3[it]);
                hpTerrD3VsMVA[iy][ipt]->Fill(mva[it], ptErr3[it] / pt3[it]);
                hEtaD3VsMVA[iy][ipt]->Fill(mva[it], eta3[it]);
                hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva[it], H2dedx3[it]);
                hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1[it], H2dedx3[it]);
              }
            }
          }
        }
    }
  }
}

void VertexCompositeTreeProducerNew::fillGEN(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);

  candSize_gen = 0;
  for (unsigned it = 0; it < genpars->size(); ++it) {

    const reco::GenParticle &trk = (*genpars)[it];

    const reco::Candidate *Dd1 = trk.daughter(0);
    const reco::Candidate *Dd2 = trk.daughter(1);

    int id = trk.pdgId();
    // if( fabs(id)!=PID_) continue; //check is target
    if (fabs(id) != PID_)
      continue; // check is target
    if (Dd1 == nullptr || Dd2 == nullptr)
      continue; // check is target
    // if( doGenDoubleDecay_ && !(abs(Dd1->pdgId()) == 421 && abs(Dd2->pdgId())
    // == 421) ) continue; //check is target if(!(trk.statusFlags().isLastCopy()
    // && (trk.pdgId()==21 || std::abs(trk.pdgId())<=6))) continue;
    if (debug_)
      std::cout << "pass id, decay dau : " << trk.numberOfDaughters() << std::endl;

    if (decayInGen_ && (trk.numberOfDaughters() != 2 && trk.numberOfDaughters() != 3))
      continue; // check 2-pron decay if target decays in Gen

    if (debug_)
      std::cout << "pass decay" << std::endl;

    candSize_gen += 1;

    pt_gen[candSize_gen - 1] = trk.pt();
    mass_gen[candSize_gen - 1] = trk.mass();
    eta_gen[candSize_gen - 1] = trk.eta();
    phi_gen[candSize_gen - 1] = trk.phi();
    status_gen[candSize_gen - 1] = trk.status();
    idself[candSize_gen - 1] = trk.pdgId();
    idmom[candSize_gen - 1] = -77;
    y_gen[candSize_gen - 1] = trk.rapidity();
    ptmom[candSize_gen - 1] = -999.0;
    etamom[candSize_gen - 1] = -999.0;
    phimom[candSize_gen - 1] = -999.0;
    ymom[candSize_gen - 1] = -999.0;
    statusmom[candSize_gen - 1] = -999;
    genDecayLength(trk, dl2D_gen[candSize_gen - 1], dl3D_gen[candSize_gen - 1], angle2D_gen[candSize_gen - 1],
                   angle3D_gen[candSize_gen - 1]);

    if (trk.numberOfMothers() != 0) {
      const reco::Candidate *mom = trk.mother();
      idmom[candSize_gen - 1] = mom->pdgId();
      ptmom[candSize_gen - 1] = mom->pt();
      etamom[candSize_gen - 1] = mom->eta();
      phimom[candSize_gen - 1] = mom->phi();
      ymom[candSize_gen - 1] = mom->rapidity();
      statusmom[candSize_gen - 1] = mom->status();
    }

    if (!decayInGen_)
      continue;

    const reco::Candidate *Dd3 = trk.daughter(2);

    iddau1[candSize_gen - 1] = fabs(Dd1->pdgId());
    iddau2[candSize_gen - 1] = fabs(Dd2->pdgId());
    if (Dd3)
      iddau3[candSize_gen - 1] = fabs(Dd3->pdgId());
    pt_gen1[candSize_gen - 1] = Dd1->pt();
    mass_gen1[candSize_gen - 1] = Dd1->mass();
    eta_gen1[candSize_gen - 1] = Dd1->eta();
    phi_gen1[candSize_gen - 1] = Dd1->phi();
    status_gen1[candSize_gen - 1] = Dd1->status();
    idself1[candSize_gen - 1] = Dd1->pdgId();

    pt_gen2[candSize_gen - 1] = Dd2->pt();
    mass_gen2[candSize_gen - 1] = Dd2->mass();
    eta_gen2[candSize_gen - 1] = Dd2->eta();
    phi_gen2[candSize_gen - 1] = Dd2->phi();
    status_gen2[candSize_gen - 1] = Dd2->status();
    idself2[candSize_gen - 1] = Dd2->pdgId();
  }
}

// ------------ method called once each job just before starting event
// loop  ------------
void VertexCompositeTreeProducerNew::beginJob() {
  TH1D::SetDefaultSumw2();

  if (!doRecoNtuple_ && !doGenNtuple_) {
    cout << "No output for either RECO or GEN!! Fix config!!" << endl;
    return;
  }

  if (twoLayerDecay_ && doMuon_) {
    cout << "Muons cannot be coming from two layer decay!! Fix config!!" << endl;
    return;
  }

  if (saveHistogram_)
    initHistogram();
  if (saveTree_)
    initTree();
}

void VertexCompositeTreeProducerNew::initHistogram() {
  for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++) {
    for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
      hMassVsMVA[iy][ipt] =
          fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d", iy, ipt), ";mva;mass(GeV)", 100, -1., 1., massHistBins_,
                         massHistPeak_ - massHistWidth_, massHistPeak_ + massHistWidth_);
      //   h3DDCAVsMVA[iy][ipt] =
      //   fs->make<TH2F>(Form("h3DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;3D
      //   DCA;",100,-1.,1.,1000,0,10); h2DDCAVsMVA[iy][ipt] =
      //   fs->make<TH2F>(Form("h2DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;2D
      //   DCA;",100,-1.,1.,1000,0,10);

      if (saveAllHistogram_) {
        hpTVsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTVsMVA_y%d_pt%d", iy, ipt), ";mva;pT;", 100, -1, 1, 100, 0, 10);
        hetaVsMVA[iy][ipt] = fs->make<TH2F>(Form("hetaVsMVA_y%d_pt%d", iy, ipt), ";mva;eta;", 100, -1., 1., 40, -4, 4);
        hyVsMVA[iy][ipt] = fs->make<TH2F>(Form("hyVsMVA_y%d_pt%d", iy, ipt), ";mva;y;", 100, -1., 1., 40, -4, 4);
        hVtxProbVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hVtxProbVsMVA_y%d_pt%d", iy, ipt), ";mva;VtxProb;", 100, -1., 1., 100, 0, 1);
        h3DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DCosPointingAngleVsMVA_y%d_pt%d", iy, ipt),
                                                           ";mva;3DCosPointingAngle;", 100, -1., 1., 100, -1, 1);
        h3DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DPointingAngleVsMVA_y%d_pt%d", iy, ipt),
                                                        ";mva;3DPointingAngle;", 100, -1., 1., 50, -3.14, 3.14);
        h2DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DCosPointingAngleVsMVA_y%d_pt%d", iy, ipt),
                                                           ";mva;2DCosPointingAngle;", 100, -1., 1., 100, -1, 1);
        h2DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DPointingAngleVsMVA_y%d_pt%d", iy, ipt),
                                                        ";mva;2DPointingAngle;", 100, -1., 1., 50, -3.14, 3.14);
        h3DDecayLengthSignificanceVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("h3DDecayLengthSignificanceVsMVA_y%d_pt%d", iy, ipt), ";mva;3DDecayLengthSignificance;",
                           100, -1., 1., 300, 0, 30);
        h2DDecayLengthSignificanceVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("h2DDecayLengthSignificanceVsMVA_y%d_pt%d", iy, ipt), ";mva;2DDecayLengthSignificance;",
                           100, -1., 1., 300, 0, 30);
        h3DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthVsMVA_y%d_pt%d", iy, ipt),
                                                      ";mva;3DDecayLength;", 100, -1., 1., 300, 0, 30);
        h2DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthVsMVA_y%d_pt%d", iy, ipt),
                                                      ";mva;2DDecayLength;", 100, -1., 1., 300, 0, 30);
        hzDCASignificancedaughter1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hzDCASignificancedaughter1VsMVA_y%d_pt%d", iy, ipt), ";mva;zDCASignificancedaughter1;",
                           100, -1., 1., 100, -10, 10);
        hxyDCASignificancedaughter1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hxyDCASignificancedaughter1VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;xyDCASignificancedaughter1;", 100, -1., 1., 100, -10, 10);
        hNHitD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hNHitD1VsMVA_y%d_pt%d", iy, ipt), ";mva;NHitD1;", 100, -1., 1., 100, 0, 100);
        hpTD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTD1VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD1;", 100, -1., 1., 100, 0, 10);
        hpTerrD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTerrD1VsMVA_y%d_pt%d", iy, ipt), ";mva;pTerrD1;", 100, -1., 1., 50, 0, 0.5);
        hEtaD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hEtaD1VsMVA_y%d_pt%d", iy, ipt), ";mva;EtaD1;", 100, -1., 1., 40, -4, 4);
        hdedxHarmonic2D1VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsMVA_y%d_pt%d", iy, ipt),
                                                        ";mva;dedxHarmonic2D1;", 100, -1., 1., 100, 0, 10);
        hdedxHarmonic2D1VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D1VsP_y%d_pt%d", iy, ipt),
                                                      ";p (GeV);dedxHarmonic2D1", 100, 0, 10, 100, 0, 10);
        hzDCASignificancedaughter2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hzDCASignificancedaughter2VsMVA_y%d_pt%d", iy, ipt), ";mva;zDCASignificancedaughter2;",
                           100, -1., 1., 100, -10, 10);
        hxyDCASignificancedaughter2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hxyDCASignificancedaughter2VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;xyDCASignificancedaughter2;", 100, -1., 1., 100, -10, 10);
        hNHitD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hNHitD2VsMVA_y%d_pt%d", iy, ipt), ";mva;NHitD2;", 100, -1., 1., 100, 0, 100);
        hpTD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTD2VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD2;", 100, -1., 1., 100, 0, 10);
        hpTerrD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTerrD2VsMVA_y%d_pt%d", iy, ipt), ";mva;pTerrD2;", 100, -1., 1., 50, 0, 0.5);
        hEtaD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hEtaD2VsMVA_y%d_pt%d", iy, ipt), ";mva;EtaD2;", 100, -1., 1., 40, -4, 4);
        hdedxHarmonic2D2VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsMVA_y%d_pt%d", iy, ipt),
                                                        ";mva;dedxHarmonic2D2;", 100, -1., 1., 100, 0, 10);
        hdedxHarmonic2D2VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D2VsP_y%d_pt%d", iy, ipt),
                                                      ";p (GeV);dedxHarmonic2D2", 100, 0, 10, 100, 0, 10);

        if (threeProngDecay_) {
          hzDCASignificancedaughter3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hzDCASignificancedaughter3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;zDCASignificancedaughter3;", 100, -1., 1., 100, -10, 10);
          hxyDCASignificancedaughter3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hxyDCASignificancedaughter3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;xyDCASignificancedaughter3;", 100, -1., 1., 100, -10, 10);
          hNHitD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hNHitD3VsMVA_y%d_pt%d", iy, ipt), ";mva;NHitD3;", 100, -1., 1., 100, 0, 100);
          hpTD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hpTD3VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD3;", 100, -1., 1., 100, 0, 10);
          hpTerrD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hpTerrD3VsMVA_y%d_pt%d", iy, ipt), ";mva;pTerrD3;", 100, -1., 1., 50, 0, 0.5);
          hEtaD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hEtaD3VsMVA_y%d_pt%d", iy, ipt), ";mva;EtaD3;", 100, -1., 1., 40, -4, 4);
          hdedxHarmonic2D3VsMVA[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsMVA_y%d_pt%d", iy, ipt),
                                                          ";mva;dedxHarmonic2D3;", 100, -1., 1., 100, 0, 10);
          hdedxHarmonic2D3VsP[iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D3VsP_y%d_pt%d", iy, ipt),
                                                        ";p (GeV);dedxHarmonic2D3", 100, 0, 10, 100, 0, 10);
        }
      }
    }
  }
}

void VertexCompositeTreeProducerNew::initTree() {
  VertexCompositeNtuple = fs->make<TTree>("VertexCompositeNtuple", "VertexCompositeNtuple");

  if (doRecoNtuple_) {

    // Event info
    VertexCompositeNtuple->Branch("Ntrkoffline", &Ntrkoffline, "Ntrkoffline/I");
    VertexCompositeNtuple->Branch("Npixel", &Npixel, "Npixel/I");
    VertexCompositeNtuple->Branch("HFsumET", &HFsumET, "HFsumET/F");
    VertexCompositeNtuple->Branch("bestvtxX", &bestvx, "bestvtxX/F");
    VertexCompositeNtuple->Branch("bestvtxY", &bestvy, "bestvtxY/F");
    VertexCompositeNtuple->Branch("bestvtxZ", &bestvz, "bestvtxZ/F");
    VertexCompositeNtuple->Branch("candSize", &candSize, "candSize/I");
    if (isCentrality_)
      VertexCompositeNtuple->Branch("centrality", &centrality, "centrality/I");
    if (isEventPlane_) {
      VertexCompositeNtuple->Branch("ephfpAngle", &ephfpAngle, "ephfpAngle[3]/F");
      VertexCompositeNtuple->Branch("ephfmAngle", &ephfmAngle, "ephfmAngle[3]/F");
      VertexCompositeNtuple->Branch("ephfpQ", &ephfpQ, "ephfpQ[3]/F");
      VertexCompositeNtuple->Branch("ephfmQ", &ephfmQ, "ephfmQ[3]/F");
      VertexCompositeNtuple->Branch("ephfpSumW", &ephfpSumW, "ephfpSumW/F");
      VertexCompositeNtuple->Branch("ephfmSumW", &ephfmSumW, "ephfmSumW/F");
    }

    // particle info
    VertexCompositeNtuple->Branch("pT", &pt, "pT[candSize]/F");
    VertexCompositeNtuple->Branch("y", &y, "y[candSize]/F");
    VertexCompositeNtuple->Branch("eta", &eta, "eta[candSize]/F");
    VertexCompositeNtuple->Branch("phi", &phi, "phi[candSize]/F");
    VertexCompositeNtuple->Branch("mass", &mass, "mass[candSize]/F");
    if (useAnyMVA_) {
      VertexCompositeNtuple->Branch("mva", &mva, "mva[candSize]/F");
      if (doubleCand_) {
        VertexCompositeNtuple->Branch("mvaDaughter1", &mva1, "mvaDaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("mvaDaughter2", &mva2, "mvaDaughter2[candSize]/F");
      }
    }
    if (useDCA_) {
      VertexCompositeNtuple->Branch("dca3D", &dca3D, "dca3D[candSize]/F");
      VertexCompositeNtuple->Branch("dcaErr3D", &dcaErr3D, "dcaErr3D[candSize]/F");
    }

    if (!isSkimMVA_) {
      // Composite candidate info RECO
      VertexCompositeNtuple->Branch("flavor", &flavor, "flavor[candSize]/F");
      VertexCompositeNtuple->Branch("VtxProb", &VtxProb, "VtxProb[candSize]/F");
      //        VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2[candSize]/F");
      //        VertexCompositeNtuple->Branch("VtxNDF",&ndf,"VtxNDF[candSize]/F");
      VertexCompositeNtuple->Branch("3DCosPointingAngle", &agl, "3DCosPointingAngle[candSize]/F");
      VertexCompositeNtuple->Branch("3DPointingAngle", &agl_abs, "3DPointingAngle[candSize]/F");
      VertexCompositeNtuple->Branch("2DCosPointingAngle", &agl2D, "2DCosPointingAngle[candSize]/F");
      VertexCompositeNtuple->Branch("2DPointingAngle", &agl2D_abs, "2DPointingAngle[candSize]/F");
      VertexCompositeNtuple->Branch("3DDecayLengthSignificance", &dlos, "3DDecayLengthSignificance[candSize]/F");
      VertexCompositeNtuple->Branch("3DDecayLength", &dl, "3DDecayLength[candSize]/F");
      VertexCompositeNtuple->Branch("2DDecayLengthSignificance", &dlos2D, "2DDecayLengthSignificance[candSize]/F");
      VertexCompositeNtuple->Branch("2DDecayLength", &dl2D, "2DDecayLength[candSize]/F");

      if (doGenMatching_) {
        VertexCompositeNtuple->Branch("isSwap", &isSwap, "isSwap[candSize]/O");
        VertexCompositeNtuple->Branch("idmom_reco", &idmom_reco1, "idmom_reco[candSize]/I");
        VertexCompositeNtuple->Branch("matchGEN", &matchGEN, "matchGEN[candSize]/O");
        VertexCompositeNtuple->Branch("idmom_reco1", &idmom_reco1, "idmom_reco1[candSize]/I");
        VertexCompositeNtuple->Branch("idBAnc_reco1", &idBAnc_reco1, "idBAnc_reco1[candSize]/I");
        VertexCompositeNtuple->Branch("matchToGen1", &matchToGen1, "matchToGen1[candSize]/I");
        if (doGenDoubleDecay_) {
          VertexCompositeNtuple->Branch("isSwap1", &isSwap1, "isSwap1[candSize]/O");
          VertexCompositeNtuple->Branch("isSwap2", &isSwap2, "isSwap2[candSize]/O");
          VertexCompositeNtuple->Branch("matchGEN1", &matchGEN1, "matchGEN1[candSize]/O");
          VertexCompositeNtuple->Branch("matchGEN2", &matchGEN2, "matchGEN2[candSize]/O");
          VertexCompositeNtuple->Branch("idmom_reco2", &idmom_reco2, "idmom_reco2[candSize]/I");
          VertexCompositeNtuple->Branch("idBAnc_reco2", &idBAnc_reco2, "idBAnc_reco2[candSize]/I");
          VertexCompositeNtuple->Branch("matchToGen2", &matchToGen2, "matchToGen2[candSize]/I");

          VertexCompositeNtuple->Branch("gen_D1ancestorFlavor", &gen_D1ancestorFlavor_,
                                        "gen_D1ancestorFlavor[candSize]I");
          VertexCompositeNtuple->Branch("gen_D1ancestorId", &gen_D1ancestorId_, "gen_D1ancestorId[candSize]I");
          VertexCompositeNtuple->Branch("gen_D1PVx", &gen_D1PVx_, "gen_D1PVx[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1PVy", &gen_D1PVy_, "gen_D1PVy[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1PVz", &gen_D1PVz_, "gen_D1PVz[candSize]F");

          VertexCompositeNtuple->Branch("gen_D1pT", &gen_D1pT_, "gen_D1pT[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1eta", &gen_D1eta_, "gen_D1eta[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1phi", &gen_D1phi_, "gen_D1phi[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1mass", &gen_D1mass_, "gen_D1mass[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1y", &gen_D1y_, "gen_D1y[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1charge", &gen_D1charge_, "gen_D1charge[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1pdgId", &gen_D1pdgId_, "gen_D1pdgId[candSize]I");

          VertexCompositeNtuple->Branch("gen_D1decayLength3D", &gen_D1decayLength3D_, "gen_D1decayLength3D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1decayLength2D", &gen_D1decayLength2D_, "gen_D1decayLength2D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1angle3D", &gen_D1angle3D_, "gen_D1angle3D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1angle2D", &gen_D1angle2D_, "gen_D1angle2D[candSize]F");

          VertexCompositeNtuple->Branch("gen_D1pTD1", &gen_D1pTD1_, "gen_D1pTD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1etaD1", &gen_D1etaD1_, "gen_D1etaD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1phiD1", &gen_D1phiD1_, "gen_D1phiD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1massD1", &gen_D1massD1_, "gen_D1massD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1yD1", &gen_D1yD1_, "gen_D1yD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1chargeD1", &gen_D1chargeD1_, "gen_D1chargeD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1pdgIdD1", &gen_D1pdgIdD1_, "gen_D1pdgIdD1[candSize]I");

          VertexCompositeNtuple->Branch("gen_D1pTD2", &gen_D1pTD2_, "gen_D1pTD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1etaD2", &gen_D1etaD2_, "gen_D1etaD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1phiD2", &gen_D1phiD2_, "gen_D1phiD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1massD2", &gen_D1massD2_, "gen_D1massD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1yD2", &gen_D1yD2_, "gen_D1yD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1chargeD2", &gen_D1chargeD2_, "gen_D1chargeD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D1pdgIdD2", &gen_D1pdgIdD2_, "gen_D1pdgIdD2[candSize]I");

          VertexCompositeNtuple->Branch("gen_D2ancestorFlavor", &gen_D2ancestorFlavor_,
                                        "gen_D2ancestorFlavor[candSize]I");
          VertexCompositeNtuple->Branch("gen_D2ancestorId", &gen_D2ancestorId_, "gen_D2ancestorId[candSize]I");
          VertexCompositeNtuple->Branch("gen_D2PVx", &gen_D2PVx_, "gen_D2PVx[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2PVy", &gen_D2PVy_, "gen_D2PVy[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2PVz", &gen_D2PVz_, "gen_D2PVz[candSize]F");

          VertexCompositeNtuple->Branch("gen_D2pT", &gen_D2pT_, "gen_D2pT[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2eta", &gen_D2eta_, "gen_D2eta[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2phi", &gen_D2phi_, "gen_D2phi[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2mass", &gen_D2mass_, "gen_D2mass[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2y", &gen_D2y_, "gen_D2y[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2charge", &gen_D2charge_, "gen_D2charge[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2pdgId", &gen_D2pdgId_, "gen_D2pdgId[candSize]I");

          VertexCompositeNtuple->Branch("gen_D2decayLength3D", &gen_D2decayLength3D_, "gen_D2decayLength3D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2decayLength2D", &gen_D2decayLength2D_, "gen_D2decayLength2D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2angle3D", &gen_D2angle3D_, "gen_D2angle3D[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2angle2D", &gen_D2angle2D_, "gen_D2angle2D[candSize]F");

          VertexCompositeNtuple->Branch("gen_D2pTD1", &gen_D2pTD1_, "gen_D2pTD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2etaD1", &gen_D2etaD1_, "gen_D2etaD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2phiD1", &gen_D2phiD1_, "gen_D2phiD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2massD1", &gen_D2massD1_, "gen_D2massD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2yD1", &gen_D2yD1_, "gen_D2yD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2chargeD1", &gen_D2chargeD1_, "gen_D2chargeD1[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2pdgIdD1", &gen_D2pdgIdD1_, "gen_D2pdgIdD1[candSize]I");

          VertexCompositeNtuple->Branch("gen_D2pTD2", &gen_D2pTD2_, "gen_D2pTD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2etaD2", &gen_D2etaD2_, "gen_D2etaD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2phiD2", &gen_D2phiD2_, "gen_D2phiD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2massD2", &gen_D2massD2_, "gen_D2massD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2yD2", &gen_D2yD2_, "gen_D2yD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2chargeD2", &gen_D2chargeD2_, "gen_D2chargeD2[candSize]F");
          VertexCompositeNtuple->Branch("gen_D2pdgIdD2", &gen_D2pdgIdD2_, "gen_D2pdgIdD2[candSize]I");
        }
      }

      if (doGenMatchingTOF_) {
        VertexCompositeNtuple->Branch("PIDD1", &pid1, "PIDD1[candSize]/I");
        VertexCompositeNtuple->Branch("PIDD2", &pid1, "PIDD2[candSize]/I");
        VertexCompositeNtuple->Branch("TOFD1", &tof1, "TOFD1[candSize]/F");
        VertexCompositeNtuple->Branch("TOFD2", &tof1, "TOFD2[candSize]/F");
      }

      // daughter & grand daughter info
      if (twoLayerDecay_) {
        VertexCompositeNtuple->Branch("flavordaughter1", &flavor1, "flavordaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("massdaughter1", &grand_mass, "massdaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("pTD1", &pt1, "pTD1[candSize]/F");
        VertexCompositeNtuple->Branch("EtaD1", &eta1, "EtaD1[candSize]/F");
        VertexCompositeNtuple->Branch("PhiD1", &phi1, "PhiD1[candSize]/F");
        VertexCompositeNtuple->Branch("VtxProbdaughter1", &grand_VtxProb, "VtxProbdaughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("VtxNDFdaughter1",&grand_ndf,"VtxNDFdaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("3DCosPointingAngle1daughter1", &grand_agl,
                                      "3DCosPointingAngledaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("3DPointingAngledaughter1", &grand_agl_abs,
                                      "3DPointingAngledaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("2DCosPointingAngledaughter1", &grand_agl2D,
                                      "2DCosPointingAngledaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("2DPointingAngledaughter1", &grand_agl2D_abs,
                                      "2DPointingAngledaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLengthSignificancedaughter1", &grand_dlos,
                                      "3DDecayLengthSignificancedaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLengthdaughter1", &grand_dl, "3DDecayLengthdaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("3DDecayLengthErrordaughter1", &grand_dlerror,
                                      "3DDecayLengthErrordaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("2DDecayLengthSignificancedaughter1", &grand_dlos2D,
                                      "2DDecayLengthSignificancedaughter1[candSize]/F");

        VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter11", &grand_dzos1,
                                      "zDCASignificanceGranddaughter11[candSize]/F");
        VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter12", &grand_dzos2,
                                      "zDCASignificanceGranddaughter12[candSize]/F");
        VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter11", &grand_dxyos1,
                                      "xyDCASignificanceGranddaughter11[candSize]/F");
        VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter12", &grand_dxyos2,
                                      "xyDCASignificanceGranddaughter12[candSize]/F");
        VertexCompositeNtuple->Branch("NHitGrandD11", &grand_nhit1, "NHitGrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("NHitGrandD12", &grand_nhit2, "NHitGrandD12[candSize]/F");
        VertexCompositeNtuple->Branch("HighPurityGranddaughter11", &grand_trkquality1,
                                      "HighPurityGranddaughter11[candSize]/O");
        VertexCompositeNtuple->Branch("HighPurityGranddaughter12", &grand_trkquality2,
                                      "HighPurityGranddaughter12[candSize]/O");
        VertexCompositeNtuple->Branch("pTGrandD11", &grand_pt1, "pTGrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("pTGrandD12", &grand_pt2, "pTGrandD12[candSize]/F");
        VertexCompositeNtuple->Branch("pTerrGrandD11", &grand_ptErr1, "pTerrGrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("pTerrGrandD12", &grand_ptErr2, "pTerrGrandD12[candSize]/F");
        //            VertexCompositeNtuple->Branch("pGrandD11",&grand_p1,"pGrandD11[candSize]/F");
        //            VertexCompositeNtuple->Branch("pGrandD12",&grand_p2,"pGrandD12[candSize]/F");
        VertexCompositeNtuple->Branch("EtaGrandD11", &grand_eta1, "EtaGrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("EtaGrandD12", &grand_eta2, "EtaGrandD12[candSize]/F");
        VertexCompositeNtuple->Branch("PhiGrandD11", &grand_phi1, "PhiGrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("PhiGrandD12", &grand_phi2, "PhiGrandD12[candSize]/F");
        //            VertexCompositeNtuple->Branch("chargeGrandD11",&grand_charge1,"chargeGrandD11[candSize]/I");
        //            VertexCompositeNtuple->Branch("chargeGrandD12",&grand_charge2,"chargeGrandD12[candSize]/I");
        VertexCompositeNtuple->Branch("dedxHarmonic2GrandD11", &grand_H2dedx1, "dedxHarmonic2GrandD11[candSize]/F");
        VertexCompositeNtuple->Branch("dedxHarmonic2GrandD12", &grand_H2dedx2, "dedxHarmonic2GrandD12[candSize]/F");
        //            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter1",&grand_T4dedx1,"dedxTruncated40Granddaughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter2",&grand_T4dedx2,"dedxTruncated40Granddaughter2[candSize]/F");
        //            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter1",&grand_trkChi1,"normalizedChi2Granddaughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter2",&grand_trkChi2,"normalizedChi2Granddaughter2[candSize]/F");
        if (doubleCand_) {
          VertexCompositeNtuple->Branch("flavordaughter2", &flavor2, "flavordaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("massdaughter2", &grand_mass2, "massdaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("pTD2", &pt2, "pTD2[candSize]/F");
          VertexCompositeNtuple->Branch("EtaD2", &eta2, "EtaD2[candSize]/F");
          VertexCompositeNtuple->Branch("PhiD2", &phi2, "PhiD2[candSize]/F");
          VertexCompositeNtuple->Branch("VtxProbdaughter2", &grand_VtxProb2, "VtxProbdaughter2[candSize]/F");
          //            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
          //            VertexCompositeNtuple->Branch("VtxNDFdaughter1",&grand_ndf,"VtxNDFdaughter1[candSize]/F");
          VertexCompositeNtuple->Branch("3DCosPointingAngle1daughter2", &grand_agl2,
                                        "3DCosPointingAngledaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("3DPointingAngledaughter2", &grand_agl_abs2,
                                        "3DPointingAngledaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("2DCosPointingAngledaughter2", &grand_agl2D2,
                                        "2DCosPointingAngledaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("2DPointingAngledaughter2", &grand_agl2D_abs2,
                                        "2DPointingAngledaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("3DDecayLengthSignificancedaughter2", &grand_dlos2,
                                        "3DDecayLengthSignificancedaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("3DDecayLengthdaughter2", &grand_dl2, "3DDecayLengthdaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("3DDecayLengthErrordaughter2", &grand_dlerror2,
                                        "3DDecayLengthErrordaughter2[candSize]/F");
          VertexCompositeNtuple->Branch("2DDecayLengthSignificancedaughter2", &grand_dlos2D2,
                                        "2DDecayLengthSignificancedaughter2[candSize]/F");

          VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter21", &grand_dzos21,
                                        "zDCASignificanceGranddaughter21[candSize]/F");
          VertexCompositeNtuple->Branch("zDCASignificanceGranddaughter22", &grand_dzos22,
                                        "zDCASignificanceGranddaughter22[candSize]/F");
          VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter21", &grand_dxyos21,
                                        "xyDCASignificanceGranddaughter21[candSize]/F");
          VertexCompositeNtuple->Branch("xyDCASignificanceGranddaughter22", &grand_dxyos22,
                                        "xyDCASignificanceGranddaughter22[candSize]/F");
          VertexCompositeNtuple->Branch("NHitGrandD21", &grand_nhit21, "NHitGrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("NHitGrandD22", &grand_nhit22, "NHitGrandD22[candSize]/F");
          VertexCompositeNtuple->Branch("HighPurityGranddaughter21", &grand_trkquality21,
                                        "HighPurityGranddaughter21[candSize]/O");
          VertexCompositeNtuple->Branch("HighPurityGranddaughter22", &grand_trkquality22,
                                        "HighPurityGranddaughter22[candSize]/O");
          VertexCompositeNtuple->Branch("pTGrandD21", &grand_pt21, "pTGrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("pTGrandD22", &grand_pt22, "pTGrandD22[candSize]/F");
          VertexCompositeNtuple->Branch("pTerrGrandD21", &grand_ptErr21, "pTerrGrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("pTerrGrandD22", &grand_ptErr22, "pTerrGrandD22[candSize]/F");
          //            VertexCompositeNtuple->Branch("pGrandD21",&grand_p21,"pGrandD1[candSize]/F");
          //            VertexCompositeNtuple->Branch("pGrandD22",&grand_p22,"pGrandD2[candSize]/F");
          VertexCompositeNtuple->Branch("EtaGrandD21", &grand_eta21, "EtaGrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("EtaGrandD22", &grand_eta22, "EtaGrandD22[candSize]/F");
          VertexCompositeNtuple->Branch("PhiGrandD21", &grand_phi21, "PhiGrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("PhiGrandD22", &grand_phi22, "PhiGrandD22[candSize]/F");
          //            VertexCompositeNtuple->Branch("chargeGrandD21",&grand_charge21,"chargeGrandD1[candSize]/I");
          //            VertexCompositeNtuple->Branch("chargeGrandD22",&grand_charge22,"chargeGrandD2[candSize]/I");
          VertexCompositeNtuple->Branch("dedxHarmonic2GrandD21", &grand_H2dedx21, "dedxHarmonic2GrandD21[candSize]/F");
          VertexCompositeNtuple->Branch("dedxHarmonic2GrandD22", &grand_H2dedx22, "dedxHarmonic2GrandD22[candSize]/F");
          //            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter1",&grand_T4dedx21,"dedxTruncated40Granddaughter1[candSize]/F");
          //            VertexCompositeNtuple->Branch("dedxTruncated40Granddaughter2",&grand_T4dedx22,"dedxTruncated40Granddaughter2[candSize]/F");
          //            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter1",&grand_trkChi21,"normalizedChi2Granddaughter1[candSize]/F");
          //            VertexCompositeNtuple->Branch("normalizedChi2Granddaughter2",&grand_trkChi22,"normalizedChi2Granddaughter2[candSize]/F");
        }
      } else {
        VertexCompositeNtuple->Branch("zDCASignificancedaughter1", &dzos1, "zDCASignificancedaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("xyDCASignificancedaughter1", &dxyos1, "xyDCASignificancedaughter1[candSize]/F");
        VertexCompositeNtuple->Branch("NHitD1", &nhit1, "NHitD1[candSize]/F");
        VertexCompositeNtuple->Branch("HighPuritydaughter1", &trkquality1, "HighPuritydaughter1[candSize]/O");
        VertexCompositeNtuple->Branch("pTD1", &pt1, "pTD1[candSize]/F");
        VertexCompositeNtuple->Branch("pTerrD1", &ptErr1, "pTerrD1[candSize]/F");
        //            VertexCompositeNtuple->Branch("pD1",&p1,"pD1[candSize]/F");
        VertexCompositeNtuple->Branch("EtaD1", &eta1, "EtaD1[candSize]/F");
        VertexCompositeNtuple->Branch("PhiD1", &phi1, "PhiD1[candSize]/F");
        //            VertexCompositeNtuple->Branch("chargeD1",&charge1,"chargeD1[candSize]/I");
        VertexCompositeNtuple->Branch("dedxHarmonic2D1", &H2dedx1, "dedxHarmonic2D1[candSize]/F");
        //            VertexCompositeNtuple->Branch("dedxTruncated40daughter1",&T4dedx1,"dedxTruncated40daughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("normalizedChi2daughter1",&trkChi1,"normalizedChi2daughter1[candSize]/F");
        VertexCompositeNtuple->Branch("zDCASignificancedaughter2", &dzos2, "zDCASignificancedaughter2[candSize]/F");
        VertexCompositeNtuple->Branch("xyDCASignificancedaughter2", &dxyos2, "xyDCASignificancedaughter2[candSize]/F");
        VertexCompositeNtuple->Branch("NHitD2", &nhit2, "NHitD2[candSize]/F");
        VertexCompositeNtuple->Branch("HighPuritydaughter2", &trkquality2, "HighPuritydaughter2[candSize]/O");
        VertexCompositeNtuple->Branch("pTD2", &pt2, "pTD2[candSize]/F");
        VertexCompositeNtuple->Branch("pTerrD2", &ptErr2, "pTerrD2[candSize]/F");
        //            VertexCompositeNtuple->Branch("pD2",&p2,"pD2[candSize]/F");
        VertexCompositeNtuple->Branch("EtaD2", &eta2, "EtaD2[candSize]/F");
        VertexCompositeNtuple->Branch("PhiD2", &phi2, "PhiD2[candSize]/F");
        //            VertexCompositeNtuple->Branch("chargeD2",&charge2,"chargeD2[candSize]/I");
        VertexCompositeNtuple->Branch("dedxHarmonic2D2", &H2dedx2, "dedxHarmonic2D2[candSize]/F");
        //            VertexCompositeNtuple->Branch("dedxTruncated40daughter2",&T4dedx2,"dedxTruncated40daughter2[candSize]/F");
        //            VertexCompositeNtuple->Branch("normalizedChi2daughter2",&trkChi2,"normalizedChi2daughter2[candSize]/F");
        if (threeProngDecay_) {
          VertexCompositeNtuple->Branch("zDCASignificancedaughter3", &dzos3, "zDCASignificancedaughter3[candSize]/F");
          VertexCompositeNtuple->Branch("xyDCASignificancedaughter3", &dxyos3,
                                        "xyDCASignificancedaughter3[candSize]/F");
          VertexCompositeNtuple->Branch("NHitD3", &nhit3, "NHitD3[candSize]/F");
          VertexCompositeNtuple->Branch("HighPuritydaughter3", &trkquality3, "HighPuritydaughter3[candSize]/O");
          VertexCompositeNtuple->Branch("pTD3", &pt1, "pTD3[candSize]/F");
          VertexCompositeNtuple->Branch("pTerrD3", &ptErr3, "pTerrD3[candSize]/F");
          VertexCompositeNtuple->Branch("EtaD3", &eta1, "EtaD3[candSize]/F");
          VertexCompositeNtuple->Branch("dedxHarmonic2D3", &H2dedx1, "dedxHarmonic2D3[candSize]/F");
        }
      }

      if (doMuon_) {
        VertexCompositeNtuple->Branch("OneStMuon1", &onestmuon1, "OneStMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("OneStMuon2", &onestmuon2, "OneStMuon2[candSize]/O");
        VertexCompositeNtuple->Branch("PFMuon1", &pfmuon1, "PFMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("PFMuon2", &pfmuon2, "PFMuon2[candSize]/O");
        VertexCompositeNtuple->Branch("GlbMuon1", &glbmuon1, "GlbMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("GlbMuon2", &glbmuon2, "GlbMuon2[candSize]/O");
        VertexCompositeNtuple->Branch("trkMuon1", &trkmuon1, "trkMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("trkMuon2", &trkmuon2, "trkMuon2[candSize]/O");
        VertexCompositeNtuple->Branch("caloMuon1", &calomuon1, "caloMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("caloMuon2", &calomuon2, "caloMuon2[candSize]/O");
        VertexCompositeNtuple->Branch("SoftMuon1", &softmuon1, "SoftMuon1[candSize]/O");
        VertexCompositeNtuple->Branch("SoftMuon2", &softmuon2, "SoftMuon2[candSize]/O");

        if (doMuonFull_) {
          VertexCompositeNtuple->Branch("nMatchedChamberD1", &nmatchedch1, "nMatchedChamberD1[candSize]/F");
          VertexCompositeNtuple->Branch("nMatchedStationD1", &nmatchedst1, "nMatchedStationD1[candSize]/F");
          VertexCompositeNtuple->Branch("EnergyDepositionD1", &matchedenergy1, "EnergyDepositionD1[candSize]/F");
          VertexCompositeNtuple->Branch("nMatchedChamberD2", &nmatchedch2, "nMatchedChamberD2[candSize]/F");
          VertexCompositeNtuple->Branch("nMatchedStationD2", &nmatchedst2, "nMatchedStationD2[candSize]/F");
          VertexCompositeNtuple->Branch("EnergyDepositionD2", &matchedenergy2, "EnergyDepositionD2[candSize]/F");
          VertexCompositeNtuple->Branch("dx1_seg", &dx1_seg_, "dx1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dy1_seg", &dy1_seg_, "dy1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dxSig1_seg", &dxSig1_seg_, "dxSig1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dySig1_seg", &dySig1_seg_, "dySig1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddxdz1_seg", &ddxdz1_seg_, "ddxdz1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddydz1_seg", &ddydz1_seg_, "ddydz1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddxdzSig1_seg", &ddxdzSig1_seg_, "ddxdzSig1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddydzSig1_seg", &ddydzSig1_seg_, "ddydzSig1_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dx2_seg", &dx2_seg_, "dx2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dy2_seg", &dy2_seg_, "dy2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dxSig2_seg", &dxSig2_seg_, "dxSig2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("dySig2_seg", &dySig2_seg_, "dySig2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddxdz2_seg", &ddxdz2_seg_, "ddxdz2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddydz2_seg", &ddydz2_seg_, "ddydz2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddxdzSig2_seg", &ddxdzSig2_seg_, "ddxdzSig2_seg[candSize]/F");
          VertexCompositeNtuple->Branch("ddydzSig2_seg", &ddydzSig2_seg_, "ddydzSig2_seg[candSize]/F");
        }
      }
    }

  } // doRecoNtuple_

  if (doGenNtuple_) {
    VertexCompositeNtuple->Branch("candSize_gen", &candSize_gen, "candSize_gen/I");
    VertexCompositeNtuple->Branch("id_gen", &idself, "id_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("mass_gen", &mass_gen, "mass_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("pT_gen", &pt_gen, "pT_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("eta_gen", &eta_gen, "eta_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("phi_gen", &phi_gen, "phi_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("y_gen", &y_gen, "y_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("status_gen", &status_gen, "status_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherID_gen", &idmom, "MotherID_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherPt_gen", &ptmom, "MotherPt_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherEta_gen", &etamom, "MotherEta_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherPhi_gen", &phimom, "MotherPhi_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherY_gen", &ymom, "MotherY_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("MotherStatus_gen", &statusmom, "MotherStatus_gen[candSize_gen]/I");
    VertexCompositeNtuple->Branch("dl2D_gen", &dl2D_gen, "dl2D_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("dl3D_gen", &dl3D_gen, "dl3D_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("angle2D_gen", &angle2D_gen, "angle2D_gen[candSize_gen]/F");
    VertexCompositeNtuple->Branch("angle3D_gen", &angle3D_gen, "angle3D_gen[candSize_gen]/F");
    if (doGenDoubleDecay_) {
      VertexCompositeNtuple->Branch("id_gen1", &idself1, "id_gen1[candSize_gen]/I");
      VertexCompositeNtuple->Branch("mass_gen1", &mass_gen1, "mass_gen1[candSize_gen]/F");
      VertexCompositeNtuple->Branch("pt_gen1", &pt_gen1, "pt_gen1[candSize_gen]/F");
      VertexCompositeNtuple->Branch("eta_gen1", &eta_gen1, "eta_gen1[candSize_gen]/F");
      VertexCompositeNtuple->Branch("phi_gen1", &phi_gen1, "phi_gen1[candSize_gen]/F");
      VertexCompositeNtuple->Branch("status_gen1", &status_gen1, "status_gen1[candSize_gen]/F");

      VertexCompositeNtuple->Branch("id_gen2", &idself2, "id_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("mass_gen2", &mass_gen2, "mass_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("pt_gen2", &pt_gen2, "pt_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("eta_gen2", &eta_gen2, "eta_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("phi_gen2", &phi_gen2, "phi_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("status_gen2", &status_gen2, "status_gen2[candSize_gen]/I");
    }

    if (decayInGen_) {

      VertexCompositeNtuple->Branch("DauID1_gen", &iddau1, "DauID1_gen[candSize_gen]/I");
      VertexCompositeNtuple->Branch("DauID2_gen", &iddau2, "DauID2_gen[candSize_gen]/I");
      VertexCompositeNtuple->Branch("DauID3_gen", &iddau3, "DauID3_gen[candSize_gen]/I");
    }
  }
}

int VertexCompositeTreeProducerNew::muAssocToTrack(const reco::TrackRef &trackref,
                                                   const edm::Handle<reco::MuonCollection> &muonh) const {
  auto muon = std::find_if(muonh->cbegin(), muonh->cend(),
                           [&](const reco::Muon &m) { return (m.track().isNonnull() && m.track() == trackref); });
  return (muon != muonh->cend() ? std::distance(muonh->cbegin(), muon) : -1);
}

// ------------ method called once each job just after ending the event
// loop  ------------
void VertexCompositeTreeProducerNew::endJob() {}

reco::GenParticleRef VertexCompositeTreeProducerNew::findMother(const reco::GenParticleRef &genParRef) {
  if (genParRef.isNull())
    return genParRef;
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId();
  const int pdg_OLD = pdg;
  while (pdg == pdg_OLD && genMomRef->numberOfMothers() > 0) {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  if (pdg == pdg_OLD)
    genMomRef = reco::GenParticleRef();
  return genMomRef;
}

void VertexCompositeTreeProducerNew::genDecayLength(const reco::GenParticle &gCand, float &gen_decayLength2D_,
                                                    float &gen_decayLength3D_, float &gen_angle2D_,
                                                    float &gen_angle3D_) {
  gen_decayLength2D_ = -99.;
  gen_decayLength3D_ = -99.;
  gen_angle2D_ = -99;
  gen_angle3D_ = -99;

  if (gCand.numberOfDaughters() == 0 || !gCand.daughter(0))
    return;
  const auto &dauVtx = gCand.daughter(0)->vertex();
  TVector3 ptosvec(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), dauVtx.Z() - genVertex_.z());
  TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
  gen_angle3D_ = secvec.Angle(ptosvec);
  gen_decayLength3D_ = ptosvec.Mag();
  TVector3 ptosvec2D(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), 0.0);
  TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
  gen_angle2D_ = secvec2D.Angle(ptosvec2D);
  gen_decayLength2D_ = ptosvec2D.Mag();
}

void VertexCompositeTreeProducerNew::getAncestorId(const reco::GenParticle &gCand, int &gen_ancestorId_,
                                                   int &gen_ancestorFlavor_) {
  gen_ancestorId_ = 0;
  gen_ancestorFlavor_ = 0;
  for (auto mothers = gCand.motherRefVector(); !mothers.empty();) {
    auto mom = mothers.at(0);
    mothers = mom->motherRefVector();
    gen_ancestorId_ = mom->pdgId();
    const auto idstr = std::to_string(std::abs(gen_ancestorId_));
    gen_ancestorFlavor_ = std::stoi(std::string{idstr.begin(), idstr.begin() + 1});
    if (idstr[0] == '5') {
      break;
    }
    if (std::abs(gen_ancestorId_) <= 40)
      break;
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducerNew);
