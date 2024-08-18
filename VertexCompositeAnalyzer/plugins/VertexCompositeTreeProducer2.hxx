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
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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

// #define DEBUG

#define PI 3.1415926
#define MAXCAN 2000

using namespace std;

class VertexCompositeTreeProducer2 : public edm::EDAnalyzer {
public:
  explicit VertexCompositeTreeProducer2(const edm::ParameterSet &);
  ~VertexCompositeTreeProducer2();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void fillRECO(const edm::Event &, const edm::EventSetup &);
  virtual void fillGEN(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  virtual void initHistogram();
  virtual void initTree();

  bool matchHadron(const reco::Candidate *_dmeson_, const reco::GenParticle &_gen_) const;
  bool checkSwap(const reco::Candidate *_dmeson_, const reco::GenParticle &_gen_) const;
  bool matchTrackdR(const reco::Candidate *_recoTrk_, const reco::Candidate *_genTrk_, bool chkchrg) const;

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
  float y1[MAXCAN];
  float y2[MAXCAN];
  float y3[MAXCAN];
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

VertexCompositeTreeProducer2::VertexCompositeTreeProducer2(const edm::ParameterSet &iConfig) {
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
};

VertexCompositeTreeProducer2::~VertexCompositeTreeProducer2() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

};

void VertexCompositeTreeProducer2::endJob() {

};

void VertexCompositeTreeProducer2::initHistogram() {
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
};

void VertexCompositeTreeProducer2::initTree() {
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
        VertexCompositeNtuple->Branch("isSwap1", &isSwap1, "isSwap1[candSize]/O");
        VertexCompositeNtuple->Branch("matchGEN1", &matchGEN1, "matchGEN1[candSize]/O");
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
        if (doGenDoubleDecay_) {
          VertexCompositeNtuple->Branch("isSwap2", &isSwap2, "isSwap2[candSize]/O");
          VertexCompositeNtuple->Branch("matchGEN2", &matchGEN2, "matchGEN2[candSize]/O");
          VertexCompositeNtuple->Branch("idmom_reco2", &idmom_reco2, "idmom_reco2[candSize]/I");
          VertexCompositeNtuple->Branch("idBAnc_reco2", &idBAnc_reco2, "idBAnc_reco2[candSize]/I");
          VertexCompositeNtuple->Branch("matchToGen2", &matchToGen2, "matchToGen2[candSize]/I");

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
        VertexCompositeNtuple->Branch("YD1", &y1, "YD1[candSize]/F");
        VertexCompositeNtuple->Branch("PhiD1", &phi1, "PhiD1[candSize]/F");
        VertexCompositeNtuple->Branch("VtxProbdaughter1", &grand_VtxProb, "VtxProbdaughter1[candSize]/F");
        //            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
        VertexCompositeNtuple->Branch("VtxNDFdaughter1", &grand_ndf, "VtxNDFdaughter1[candSize]/F");
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
          VertexCompositeNtuple->Branch("YD2", &y2, "YD2[candSize]/F");
          VertexCompositeNtuple->Branch("PhiD2", &phi2, "PhiD2[candSize]/F");
          VertexCompositeNtuple->Branch("VtxProbdaughter2", &grand_VtxProb2, "VtxProbdaughter2[candSize]/F");
          //            VertexCompositeNtuple->Branch("VtxChi2daughter1",&grand_vtxChi2,"VtxChi2daughter1[candSize]/F");
          VertexCompositeNtuple->Branch("VtxNDFdaughter2", &grand_ndf2, "VtxNDFdaughter2[candSize]/F");
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
        VertexCompositeNtuple->Branch("YD1", &y1, "Y D1[candSize]/F");
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
        VertexCompositeNtuple->Branch("YD2", &y2, "Y D2[candSize]/F");
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
      VertexCompositeNtuple->Branch("status_gen1", &status_gen1, "status_gen1[candSize_gen]/I");

      VertexCompositeNtuple->Branch("id_gen2", &idself2, "id_gen2[candSize_gen]/I");
      VertexCompositeNtuple->Branch("mass_gen2", &mass_gen2, "mass_gen2[candSize_gen]/F");
      VertexCompositeNtuple->Branch("pt_gen2", &pt_gen2, "pt_gen2[candSize_gen]/F");
      VertexCompositeNtuple->Branch("eta_gen2", &eta_gen2, "eta_gen2[candSize_gen]/F");
      VertexCompositeNtuple->Branch("phi_gen2", &phi_gen2, "phi_gen2[candSize_gen]/F");
      VertexCompositeNtuple->Branch("status_gen2", &status_gen2, "status_gen2[candSize_gen]/I");
    }

    if (decayInGen_) {
      VertexCompositeNtuple->Branch("DauID1_gen", &iddau1, "DauID1_gen[candSize_gen]/I");
      VertexCompositeNtuple->Branch("DauID2_gen", &iddau2, "DauID2_gen[candSize_gen]/I");
      VertexCompositeNtuple->Branch("DauID3_gen", &iddau3, "DauID3_gen[candSize_gen]/I");
    }
  }
};

bool VertexCompositeTreeProducer2::matchHadron(const reco::Candidate *_dmeson_, const reco::GenParticle &_gen_) const {
  reco::Candidate const *reco_trk1 = _dmeson_->daughter(0);
  reco::Candidate const *reco_trk2 = _dmeson_->daughter(1);

  reco::Candidate const *gen_trk1 = _gen_.daughter(0);
  reco::Candidate const *gen_trk2 = _gen_.daughter(1);

  bool match = false;
#ifdef DEBUG
  // cout << "Match 1 on 1 " ;
#endif
  if (matchTrackdR(reco_trk1, gen_trk1, true)) {
#ifdef DEBUG
    // cout << "Match 2 on 2 " ;
#endif
    if (matchTrackdR(reco_trk2, gen_trk2, true)) {
      match = true;
#ifdef DEBUG
      // cout << endl ;
#endif
      return match;
    }
  }
#ifdef DEBUG
  // cout << "Match 2 on 1 " ;
#endif
  if (matchTrackdR(reco_trk2, gen_trk1, true)) {
#ifdef DEBUG
    // cout << "Match 1 on 2 " ;
#endif
    if (matchTrackdR(reco_trk1, gen_trk2, true)) {
      match = true;
#ifdef DEBUG
      // cout << endl ;
#endif
      return match;
    }
  }
#ifdef DEBUG
  // cout << endl ;
#endif
  return match;
};

bool VertexCompositeTreeProducer2::checkSwap(const reco::Candidate *_dmeson_, const reco::GenParticle &_gen_) const {
  return _dmeson_->pdgId() != _gen_.pdgId();
};

bool VertexCompositeTreeProducer2::matchTrackdR(const reco::Candidate *_recoTrk_, const reco::Candidate *_genTrk_,
                                                bool chkchrg) const {
  bool pass = false;
  // deltaR_
#ifdef DEBUG
  // cout << "Charge : " << _recoTrk_->charge() << ", " <<  _genTrk_->charge() ;
#endif
  if (chkchrg && (_recoTrk_->charge() != _genTrk_->charge())) {
#ifdef DEBUG
    // cout << endl ;
#endif
    return false;
  }
  const double dR = reco::deltaR(*_recoTrk_, *_genTrk_);
#ifdef DEBUG
  // cout << ", mathTrk : (" << _recoTrk_->eta() << ", " << _genTrk_->eta() <<
  // "), (" << _recoTrk_->phi() << ", " << _genTrk_->phi() << ") -> dR = "<< dR
  // ;
#endif
  if (dR < deltaR_)
    pass = true;
#ifdef DEBUG
    // cout << Form(", dR, deltaR_ pass : %.3f %.3f %d", dR, deltaR_, pass) <<
    // endl;
#endif
  return pass;
};

// define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducer2);
