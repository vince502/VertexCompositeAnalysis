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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

//
// class decleration
//

#define PI 3.1416

using namespace std;

class VertexCompositeNtupleProducer : public edm::EDAnalyzer {
public:
  explicit VertexCompositeNtupleProducer(const edm::ParameterSet &);
  ~VertexCompositeNtupleProducer();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void fillRECO(const edm::Event &, const edm::EventSetup &);
  virtual void fillGEN(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  virtual void initHistogram();
  virtual void initTree();

  int muAssocToTrack(const reco::TrackRef &trackref,
                     const edm::Handle<reco::MuonCollection> &muonh) const;

  reco::GenParticleRef findMother(const reco::GenParticleRef &);
  void genDecayLength(const reco::GenParticle &);
  void getAncestorId(const reco::GenParticle &);

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree *VertexCompositeNtuple;
  TTree *genCandidateNtuple;
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
  TH2F *hzDCASignificanceDaugther1VsMVA[6][10];
  TH2F *hxyDCASignificanceDaugther1VsMVA[6][10];
  TH2F *hNHitD1VsMVA[6][10];
  TH2F *hpTD1VsMVA[6][10];
  TH2F *hpTerrD1VsMVA[6][10];
  TH2F *hEtaD1VsMVA[6][10];
  TH2F *hdedxHarmonic2D1VsMVA[6][10];
  TH2F *hdedxHarmonic2D1VsP[6][10];
  TH2F *hzDCASignificanceDaugther2VsMVA[6][10];
  TH2F *hxyDCASignificanceDaugther2VsMVA[6][10];
  TH2F *hNHitD2VsMVA[6][10];
  TH2F *hpTD2VsMVA[6][10];
  TH2F *hpTerrD2VsMVA[6][10];
  TH2F *hEtaD2VsMVA[6][10];
  TH2F *hdedxHarmonic2D2VsMVA[6][10];
  TH2F *hdedxHarmonic2D2VsP[6][10];
  TH2F *hzDCASignificanceDaugther3VsMVA[6][10];
  TH2F *hxyDCASignificanceDaugther3VsMVA[6][10];
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
  float HFsumETPlus;
  float HFsumETMinus;
  float ZDCPlus;
  float ZDCMinus;
  float bestvx;
  float bestvy;
  float bestvz;

  // Composite candidate info
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

  // dau candidate info
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

  // dau info
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

  // grand-dau info
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

  // dau muon info
  bool onestmuon1;
  bool onestmuon2;
  bool pfmuon1;
  bool pfmuon2;
  bool glbmuon1;
  bool glbmuon2;
  bool trkmuon1;
  bool trkmuon2;
  bool calomuon1;
  bool calomuon2;
  bool softmuon1;
  bool softmuon2;
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
  int gen_ancestorFlavor_;
  int gen_ancestorId_;
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

  // tokens
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection>
      recoVertexCompositeCandidateCollection_Token_;
  edm::EDGetTokenT<MVACollection> MVAValues_Token_;

  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token1_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token2_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
  edm::EDGetTokenT<reco::MuonCollection> tok_muon_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

  // for DCA
  edm::EDGetTokenT<std::vector<float>> tok_DCAVal_;
  edm::EDGetTokenT<std::vector<float>> tok_DCAErr_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>
    SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// static data member definitions
//

//
// constructors and destructor
//

VertexCompositeNtupleProducer::VertexCompositeNtupleProducer(
    const edm::ParameterSet &iConfig) {
  // options
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
  if (threeProngDecay_)
    PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");

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
  tok_offlinePV_ = consumes<reco::VertexCollection>(
      iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  tok_generalTrk_ = consumes<reco::TrackCollection>(
      iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
  recoVertexCompositeCandidateCollection_Token_ =
      consumes<reco::VertexCompositeCandidateCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>(
              "VertexCompositeCollection"));
  MVAValues_Token_ = consumes<MVACollection>(
      iConfig.getParameter<edm::InputTag>("MVACollection"));
  tok_muon_ = consumes<reco::MuonCollection>(
      iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
  Dedx_Token1_ =
      consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
  Dedx_Token2_ =
      consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxTruncated40"));
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(
      iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

  isCentrality_ = false;
  if (iConfig.exists("isCentrality"))
    isCentrality_ = iConfig.getParameter<bool>("isCentrality");
  if (isCentrality_) {
    tok_centBinLabel_ = consumes<int>(
        iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(
        iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  if (useAnyMVA_ && iConfig.exists("MVACollection"))
    MVAValues_Token_ = consumes<MVACollection>(
        iConfig.getParameter<edm::InputTag>("MVACollection"));
  if (iConfig.exists("DCAValCollection") &&
      iConfig.exists("DCAErrCollection")) {
    useDCA_ = true;
    tok_DCAVal_ = consumes<std::vector<float>>(
        iConfig.getParameter<edm::InputTag>("DCAValCollection"));
    tok_DCAErr_ = consumes<std::vector<float>>(
        iConfig.getParameter<edm::InputTag>("DCAErrCollection"));
  }
}

VertexCompositeNtupleProducer::~VertexCompositeNtupleProducer() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------
void VertexCompositeNtupleProducer::analyze(const edm::Event &iEvent,
                                            const edm::EventSetup &iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;

  if (doGenNtuple_)
    fillGEN(iEvent, iSetup);
  if (doRecoNtuple_)
    fillRECO(iEvent, iSetup);
}

void VertexCompositeNtupleProducer::fillRECO(const edm::Event &iEvent,
                                             const edm::EventSetup &iSetup) {
  // get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tok_generalTrk_, tracks);

  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
  iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,
                    v0candidates);
  const reco::VertexCompositeCandidateCollection *v0candidates_ =
      v0candidates.product();

  edm::Handle<MVACollection> mvavalues;
  if (useAnyMVA_) {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    assert((*mvavalues).size() == v0candidates->size());
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

    HFsumETPlus = cent->EtHFtowerSumPlus();
    HFsumETMinus = cent->EtHFtowerSumMinus();
    ZDCPlus = cent->zdcSumPlus();
    ZDCMinus = cent->zdcSumMinus();
    Npixel = cent->multiplicityPixel();
    //      int ntrk = cent->Ntracks();
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
      double dzerror =
          sqrt(trk.dzError() * trk.dzError() + bestvzError * bestvzError);
      double dxyerror =
          sqrt(trk.d0Error() * trk.d0Error() + bestvxError * bestvyError);

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

          if (abs(Dd1->pdgId()) == PID_dau1_ &&
              abs(Dd2->pdgId()) == PID_dau2_) {
            idxs = permutations;
            break;
          }
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));
      } else {
        do {
          auto Dd1 = trk.daughter(permutations.at(0));
          auto Dd2 = trk.daughter(permutations.at(1));
          auto Dd3 = trk.daughter(permutations.at(2));

          if (abs(Dd1->pdgId()) == PID_dau1_ &&
              abs(Dd2->pdgId()) == PID_dau2_ &&
              abs(Dd3->pdgId() == PID_dau3_)) {
            idxs = permutations;
            break;
          }
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));
      }

      if (decayInGen_ && idxs.empty())
        continue;
      genRefs.push_back(reco::GenParticleRef(genpars, it));
    }
    // if (genRefs.size()>1) std::cout << "More than one target of generated
    // particles\n";
  }

  // RECO Candidate info
  for (unsigned int it = 0; it < v0candidates_->size(); ++it) {

    const reco::VertexCompositeCandidate &trk = (*v0candidates_)[it];

    double secvz = -999.9, secvx = -999.9, secvy = -999.9;
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();

    eta = trk.eta();
    phi = trk.phi();
    y = trk.rapidity();
    pt = trk.pt();
    flavor = trk.pdgId() / abs(trk.pdgId());

    mva = 0.0;
    if (useAnyMVA_)
      mva = (*mvavalues)[it];

    dca3D = -1.0;
    dcaErr3D = -1.0;
    if (useDCA_) {
      dca3D = dcaValues->at(it);
      dcaErr3D = dcaErrors->at(it);
    }

    double px = trk.px();
    double py = trk.py();
    double pz = trk.pz();
    mass = trk.mass();

    const reco::Candidate *d1 = trk.daughter(0);
    const reco::Candidate *d2 = trk.daughter(1);
    const reco::Candidate *d3 = 0;
    if (threeProngDecay_)
      d3 = trk.daughter(2);

    // Gen match
    if (doGenMatching_) {
      matchGEN = false;
      isSwap = false;
      idmom_reco = -77;

      // something wrong here because I cannot gaunrantee there is only one
      // generated particle is matched here.
      const auto nGen = genRefs.size();
      for (unsigned int igen = 0; igen < nGen; igen++) {
        const auto &genRef = genRefs.at(igen);

        reco::Candidate const *recoDaus[3] = {nullptr, nullptr, nullptr};
        reco::Candidate const *genDaus[3] = {nullptr, nullptr, nullptr};

        const auto nGenDau = genRef->numberOfDaughters();
        std::vector<unsigned int> permutations(nGenDau);
        std::iota(permutations.begin(), permutations.end(), 0);
        std::sort(permutations.begin(), permutations.end());

        do {
          matchGEN = false;
          for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
            genDaus[iDau] = genRef->daughter(permutations.at(iDau));
            recoDaus[iDau] = trk.daughter(iDau);
          }

          for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
            const double dR =
                reco::deltaR(genDaus[iDau]->eta(), genDaus[iDau]->phi(),
                             recoDaus[iDau]->eta(), recoDaus[iDau]->phi());
            const double dPt = abs(genDaus[iDau]->pt() - recoDaus[iDau]->pt()) /
                               recoDaus[iDau]->pt();
            const bool unMatchCharge =
                genDaus[iDau]->charge() != recoDaus[iDau]->charge();
            const bool unMatchDR = dR > deltaR_;
            const bool unMatchDPt = dPt > 0.5;
            matchGEN = matchGEN || unMatchCharge || unMatchDR || unMatchDPt;
          }
          matchGEN = !matchGEN;
          if (matchGEN)
            break;
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));

        for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
          const double diffMass =
              abs(genDaus[iDau]->mass() - recoDaus[iDau]->mass());
          isSwap = isSwap || diffMass > 0.01;
        }
        if (matchGEN) {
          auto mom_ref = findMother(genRef);
          if (mom_ref.isNonnull())
            idmom_reco = mom_ref->pdgId();

          gen_pT_ = genRef->pt();
          gen_eta_ = genRef->eta();
          gen_phi_ = genRef->phi();
          gen_mass_ = genRef->mass();
          gen_y_ = genRef->rapidity();
          gen_charge_ = genRef->charge();
          gen_pdgId_ = genRef->pdgId();

          // all done in genDecayLength
          // gen_decayLength3D_;
          // gen_decayLength2D_;
          // gen_angle3D_;
          // gen_angle2D_;
          genDecayLength(*genRef);
          getAncestorId(*genRef);

          gen_pTD1_ = genDaus[0]->pt();
          gen_etaD1_ = genDaus[0]->eta();
          gen_phiD1_ = genDaus[0]->phi();
          gen_massD1_ = genDaus[0]->mass();
          gen_yD1_ = genDaus[0]->rapidity();
          gen_chargeD1_ = genDaus[0]->charge();
          gen_pdgIdD1_ = genDaus[0]->pdgId();

          gen_pTD2_ = genDaus[1]->pt();
          gen_etaD2_ = genDaus[1]->eta();
          gen_phiD2_ = genDaus[1]->phi();
          gen_massD2_ = genDaus[1]->mass();
          gen_yD2_ = genDaus[1]->rapidity();
          gen_chargeD2_ = genDaus[1]->charge();
          gen_pdgIdD2_ = genDaus[1]->pdgId();
        } else {
          gen_pT_ = -99;
          gen_eta_ = -99;
          gen_phi_ = -99;
          gen_mass_ = -99;
          gen_y_ = -99;

          gen_decayLength3D_ = -99;
          gen_decayLength2D_ = -99;
          gen_angle3D_ = -99;
          gen_angle2D_ = -99;

          gen_pTD1_ = -99;
          gen_etaD1_ = -99;
          gen_phiD1_ = -99;
          gen_massD1_ = -99;
          gen_yD1_ = -99;

          gen_pTD2_ = -99;
          gen_etaD2_ = -99;
          gen_phiD2_ = -99;
          gen_massD2_ = -99;
          gen_yD2_ = -99;
        }
        if (matchGEN)
          break;
      }
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
    pt1 = d1->pt();
    pt2 = d2->pt();

    // momentum
    p1 = d1->p();
    p2 = d2->p();

    // eta
    eta1 = d1->eta();
    eta2 = d2->eta();

    // phi
    phi1 = d1->phi();
    phi2 = d2->phi();

    // charge
    charge1 = d1->charge();
    charge2 = d2->charge();

    double pxd3 = -999.9;
    double pyd3 = -999.9;
    double pzd3 = -999.9;
    if (threeProngDecay_ && d3) {
      pxd3 = d3->px();
      pyd3 = d3->py();
      pzd3 = d3->pz();
      pt3 = d3->pt();
      p3 = d3->p();
      eta3 = d3->eta();
      phi3 = d3->phi();
      charge3 = d3->charge();
    }
    TVector3 dauvec3(pxd3, pyd3, pzd3);

    // useless
    pid1 = -99999;
    pid2 = -99999;
    if (doGenMatchingTOF_) {
      for (unsigned it = 0; it < genpars->size(); ++it) {

        const reco::GenParticle &trk = (*genpars)[it];

        if (trk.pt() < 0.001)
          continue;

        int id = trk.pdgId();
        TVector3 trkvect(trk.px(), trk.py(), trk.pz());

        // cout<<trk.pdgId()<<" "<<trk.numberOfDaughters()<<"
        // "<<trk.charge()<<endl;
        if (fabs(id) != PID_ && trk.charge()) {
          // matching daughter 1
          double deltaR = trkvect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt1) / pt1) < 0.5 &&
              trk.charge() == charge1 && pid1 == -99999) {
            pid1 = id;
            //      tof1 = ;
          }

          // matching daughter 2
          deltaR = trkvect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt2) / pt2) < 0.5 &&
              trk.charge() == charge2 && pid2 == -99999) {
            pid2 = id;
            //      tof2 = ;
          }
        }

        if (fabs(id) == PID_ && trk.numberOfDaughters() == 2) {
          // cout<<trk.pdgId()<<" "<<trk.numberOfDaughters()<<"
          // "<<trk.charge()<<endl;
          const reco::Candidate *Dd1 = trk.daughter(0);
          const reco::Candidate *Dd2 = trk.daughter(1);
          TVector3 d1vect(Dd1->px(), Dd1->py(), Dd1->pz());
          TVector3 d2vect(Dd2->px(), Dd2->py(), Dd2->pz());
          int id1 = Dd1->pdgId();
          int id2 = Dd2->pdgId();

          double deltaR = d1vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt1) / pt1) < 0.5 &&
              Dd1->charge() == charge1 && pid1 == -99999) {
            pid1 = id1;
          }
          deltaR = d2vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt1) / pt1) < 0.5 &&
              Dd2->charge() == charge1 && pid1 == -99999) {
            pid1 = id1;
          }

          deltaR = d1vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt2) / pt2) < 0.5 &&
              Dd1->charge() == charge2 && pid2 == -99999) {
            pid2 = id2;
          }
          deltaR = d2vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt2) / pt2) < 0.5 &&
              Dd2->charge() == charge2 && pid2 == -99999) {
            pid2 = id2;
          }
        }

        if (pid1 != -99999 && pid2 != -99999)
          break;
      }
    }

    // vtxChi2
    vtxChi2 = trk.vertexChi2();
    ndf = trk.vertexNdof();
    VtxProb = TMath::Prob(vtxChi2, ndf);

    // PAngle
    TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
    TVector3 secvec(px, py, pz);

    TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
    TVector3 secvec2D(px, py, 0);

    agl = cos(secvec.Angle(ptosvec));
    agl_abs = secvec.Angle(ptosvec);

    agl2D = cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs = secvec2D.Angle(ptosvec2D);

    // Decay length 3D

    SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
    SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

    dl = ROOT::Math::Mag(distanceVector);
    dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl;

    dlos = dl / dlerror;

    // Decay length 2D
    SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1),
                vtx.covariance(1, 1), 0, 0, 0);
    SVector6 v2(trk.vertexCovariance(0, 0), trk.vertexCovariance(0, 1),
                trk.vertexCovariance(1, 1), 0, 0, 0);

    SMatrixSym3D sv1(v1);
    SMatrixSym3D sv2(v2);

    SMatrixSym3D totalCov2D = sv1 + sv2;
    SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

    dl2D = ROOT::Math::Mag(distanceVector2D);
    double dl2Derror =
        sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / dl2D;

    dlos2D = dl2D / dl2Derror;

    // trk info
    auto dau1 = d1->get<reco::TrackRef>();
    if (!twoLayerDecay_) {
      // trk quality
      trkquality1 = dau1->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx1 = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx1 = dEdxTrack[dau1].dEdx();
      }

      T4dedx1 = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx1 = dEdxTrack[dau1].dEdx();
      }

      // track Chi2
      trkChi1 = dau1->normalizedChi2();

      // track pT error
      ptErr1 = dau1->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit1 = dau1->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest1 = dau1->dz(bestvtx);
      double dxybest1 = dau1->dxy(bestvtx);
      double dzerror1 =
          sqrt(dau1->dzError() * dau1->dzError() + bestvzError * bestvzError);
      double dxyerror1 =
          sqrt(dau1->d0Error() * dau1->d0Error() + bestvxError * bestvyError);

      dzos1 = dzbest1 / dzerror1;
      dxyos1 = dxybest1 / dxyerror1;
    }

    auto dau2 = d2->get<reco::TrackRef>();

    // trk quality
    trkquality2 = dau2->quality(reco::TrackBase::highPurity);

    // trk dEdx
    H2dedx2 = -999.9;

    if (dEdxHandle1.isValid()) {
      const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
      H2dedx2 = dEdxTrack[dau2].dEdx();
    }

    T4dedx2 = -999.9;

    if (dEdxHandle2.isValid()) {
      const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
      T4dedx2 = dEdxTrack[dau2].dEdx();
    }

    // track Chi2
    trkChi2 = dau2->normalizedChi2();

    // track pT error
    ptErr2 = dau2->ptError();

    // vertexCovariance 00-xError 11-y 22-z
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();

    // trkNHits
    nhit2 = dau2->numberOfValidHits();

    // DCA
    math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

    double dzbest2 = dau2->dz(bestvtx);
    double dxybest2 = dau2->dxy(bestvtx);
    double dzerror2 =
        sqrt(dau2->dzError() * dau2->dzError() + bestvzError * bestvzError);
    double dxyerror2 =
        sqrt(dau2->d0Error() * dau2->d0Error() + bestvxError * bestvyError);

    dzos2 = dzbest2 / dzerror2;
    dxyos2 = dxybest2 / dxyerror2;

    if (doMuon_) {
      edm::Handle<reco::MuonCollection> theMuonHandle;
      iEvent.getByToken(tok_muon_, theMuonHandle);

      nmatchedch1 = -1;
      nmatchedst1 = -1;
      matchedenergy1 = -1;
      nmatchedch2 = -1;
      nmatchedst2 = -1;
      matchedenergy2 = -1;

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

      onestmuon1 = false;
      pfmuon1 = false;
      glbmuon1 = false;
      trkmuon1 = false;
      calomuon1 = false;
      softmuon1 = false;
      onestmuon2 = false;
      pfmuon2 = false;
      glbmuon2 = false;
      trkmuon2 = false;
      calomuon2 = false;
      softmuon2 = false;

      const int muId1 = muAssocToTrack(dau1, theMuonHandle);
      const int muId2 = muAssocToTrack(dau2, theMuonHandle);

      if (muId1 != -1) {
        const reco::Muon &cand = (*theMuonHandle)[muId1];

        onestmuon1 = muon::isGoodMuon(
            cand, muon::selectionTypeFromString("TMOneStationTight"));
        pfmuon1 = cand.isPFMuon();
        glbmuon1 = cand.isGlobalMuon();
        trkmuon1 = cand.isTrackerMuon();
        calomuon1 = cand.isCaloMuon();

        if (
            // glbmuon1 &&
            trkmuon1 &&
            cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                5 &&
            cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
            fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
            fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
          softmuon1 = true;
      }

      if (muId2 != -1) {
        const reco::Muon &cand = (*theMuonHandle)[muId2];

        onestmuon2 = muon::isGoodMuon(
            cand, muon::selectionTypeFromString("TMOneStationTight"));
        pfmuon2 = cand.isPFMuon();
        glbmuon2 = cand.isGlobalMuon();
        trkmuon2 = cand.isTrackerMuon();
        calomuon2 = cand.isCaloMuon();

        if (
            // glbmuon2 &&
            trkmuon2 &&
            cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                5 &&
            cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
            fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
            fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
          softmuon2 = true;
      }

      if (doMuonFull_) {

        if (muId1 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId1];

          nmatchedch1 = cand.numberOfMatches();
          nmatchedst1 = cand.numberOfMatchedStations();

          reco::MuonEnergy muenergy = cand.calEnergy();
          matchedenergy1 = muenergy.hadMax;

          const std::vector<reco::MuonChamberMatch> &muchmatches =
              cand.matches();

          for (unsigned int ich = 0; ich < muchmatches.size(); ich++) {
            x_exp = muchmatches[ich].x;
            y_exp = muchmatches[ich].y;
            xerr_exp = muchmatches[ich].xErr;
            yerr_exp = muchmatches[ich].yErr;
            dxdz_exp = muchmatches[ich].dXdZ;
            dydz_exp = muchmatches[ich].dYdZ;
            dxdzerr_exp = muchmatches[ich].dXdZErr;
            dydzerr_exp = muchmatches[ich].dYdZErr;

            std::vector<reco::MuonSegmentMatch> musegmatches =
                muchmatches[ich].segmentMatches;

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

              if (sqrt((x_seg - x_exp) * (x_seg - x_exp) +
                       (y_seg - y_exp) * (y_seg - y_exp)) <
                  sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                dx_seg = x_seg - x_exp;
                dy_seg = y_seg - y_exp;
                dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                dxSig_seg = dx_seg / dxerr_seg;
                dySig_seg = dy_seg / dyerr_seg;
                ddxdz_seg = dxdz_seg - dxdz_exp;
                ddydz_seg = dydz_seg - dydz_exp;
                ddxdzerr_seg =
                    sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                ddydzerr_seg =
                    sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                ddydzSig_seg = ddydz_seg / ddydzerr_seg;
              }
            }

            dx1_seg_ = dx_seg;
            dy1_seg_ = dy_seg;
            dxSig1_seg_ = dxSig_seg;
            dySig1_seg_ = dySig_seg;
            ddxdz1_seg_ = ddxdz_seg;
            ddydz1_seg_ = ddydz_seg;
            ddxdzSig1_seg_ = ddxdzSig_seg;
            ddydzSig1_seg_ = ddydzSig_seg;
          }
        }

        if (muId2 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId2];

          nmatchedch2 = cand.numberOfMatches();
          nmatchedst2 = cand.numberOfMatchedStations();

          reco::MuonEnergy muenergy = cand.calEnergy();
          matchedenergy2 = muenergy.hadMax;

          const std::vector<reco::MuonChamberMatch> &muchmatches =
              cand.matches();
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

            std::vector<reco::MuonSegmentMatch> musegmatches =
                muchmatches[ich].segmentMatches;

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

              if (sqrt((x_seg - x_exp) * (x_seg - x_exp) +
                       (y_seg - y_exp) * (y_seg - y_exp)) <
                  sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                dx_seg = x_seg - x_exp;
                dy_seg = y_seg - y_exp;
                dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                dxSig_seg = dx_seg / dxerr_seg;
                dySig_seg = dy_seg / dyerr_seg;
                ddxdz_seg = dxdz_seg - dxdz_exp;
                ddydz_seg = dydz_seg - dydz_exp;
                ddxdzerr_seg =
                    sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                ddydzerr_seg =
                    sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                ddydzSig_seg = ddydz_seg / ddydzerr_seg;
              }
            }

            dx2_seg_ = dx_seg;
            dy2_seg_ = dy_seg;
            dxSig2_seg_ = dxSig_seg;
            dySig2_seg_ = dySig_seg;
            ddxdz2_seg_ = ddxdz_seg;
            ddydz2_seg_ = ddydz_seg;
            ddxdzSig2_seg_ = ddxdzSig_seg;
            ddydzSig2_seg_ = ddydzSig_seg;
          }
        }
      } // doMuonFull
    }

    if (twoLayerDecay_) {
      grand_mass = d1->mass();

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

      grand_trkquality1 = gdau1->quality(reco::TrackBase::highPurity);
      grand_trkquality2 = gdau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      grand_H2dedx1 = -999.9;
      grand_H2dedx2 = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        grand_H2dedx1 = dEdxTrack[gdau1].dEdx();
        grand_H2dedx2 = dEdxTrack[gdau2].dEdx();
      }

      grand_T4dedx1 = -999.9;
      grand_T4dedx2 = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        grand_T4dedx1 = dEdxTrack[gdau1].dEdx();
        grand_T4dedx2 = dEdxTrack[gdau2].dEdx();
      }

      // track pt
      grand_pt1 = gd1->pt();
      grand_pt2 = gd2->pt();

      // track momentum
      grand_p1 = gd1->p();
      grand_p2 = gd2->p();

      // track eta
      grand_eta1 = gd1->eta();
      grand_eta2 = gd2->eta();

      // track charge
      grand_charge1 = gd1->charge();
      grand_charge2 = gd2->charge();

      // track Chi2
      grand_trkChi1 = gdau1->normalizedChi2();
      grand_trkChi2 = gdau2->normalizedChi2();

      // track pT error
      grand_ptErr1 = gdau1->ptError();
      grand_ptErr2 = gdau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = d1->vz();
      secvx = d1->vx();
      secvy = d1->vy();

      // trkNHits
      grand_nhit1 = gdau1->numberOfValidHits();
      grand_nhit2 = gdau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double gdzbest1 = gdau1->dz(bestvtx);
      double gdxybest1 = gdau1->dxy(bestvtx);
      double gdzerror1 =
          sqrt(gdau1->dzError() * gdau1->dzError() + bestvzError * bestvzError);
      double gdxyerror1 =
          sqrt(gdau1->d0Error() * gdau1->d0Error() + bestvxError * bestvyError);

      grand_dzos1 = gdzbest1 / gdzerror1;
      grand_dxyos1 = gdxybest1 / gdxyerror1;

      double gdzbest2 = gdau2->dz(bestvtx);
      double gdxybest2 = gdau2->dxy(bestvtx);
      double gdzerror2 =
          sqrt(gdau2->dzError() * gdau2->dzError() + bestvzError * bestvzError);
      double gdxyerror2 =
          sqrt(gdau2->d0Error() * gdau2->d0Error() + bestvxError * bestvyError);

      grand_dzos2 = gdzbest2 / gdzerror2;
      grand_dxyos2 = gdxybest2 / gdxyerror2;

      // vtxChi2
      grand_vtxChi2 = d1->vertexChi2();
      grand_ndf = d1->vertexNdof();
      grand_VtxProb = TMath::Prob(grand_vtxChi2, grand_ndf);

      // PAngle
      TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
      TVector3 secvec(d1->px(), d1->py(), d1->pz());

      TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
      TVector3 secvec2D(d1->px(), d1->py(), 0);

      grand_agl = cos(secvec.Angle(ptosvec));
      grand_agl_abs = secvec.Angle(ptosvec);

      grand_agl2D = cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs = secvec2D.Angle(ptosvec2D);

      // Decay length 3D
      SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
      SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

      grand_dl = ROOT::Math::Mag(distanceVector);
      grand_dlerror =
          sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / grand_dl;

      grand_dlos = grand_dl / grand_dlerror;

      // Decay length 2D
      SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1),
                  vtx.covariance(1, 1), 0, 0, 0);
      SVector6 v2(d1->vertexCovariance(0, 0), d1->vertexCovariance(0, 1),
                  d1->vertexCovariance(1, 1), 0, 0, 0);

      SMatrixSym3D sv1(v1);
      SMatrixSym3D sv2(v2);

      SMatrixSym3D totalCov2D = sv1 + sv2;
      SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

      double gdl2D = ROOT::Math::Mag(distanceVector2D);
      double gdl2Derror =
          sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / gdl2D;

      grand_dlos2D = gdl2D / gdl2Derror;
    }

    if (saveTree_)
      VertexCompositeNtuple->Fill();
    if (saveHistogram_) {
      for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++)
        for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
          if (pt < pTBins_[ipt + 1] && pt > pTBins_[ipt] &&
              y < yBins_[iy + 1] && y > yBins_[iy]) {
            hMassVsMVA[iy][ipt]->Fill(mva, mass);
            //                h3DDCAVsMVA[iy][ipt]->Fill(mva,dl*sin(agl_abs));
            //                h2DDCAVsMVA[iy][ipt]->Fill(mva,dl2D*sin(agl2D_abs));

            if (saveAllHistogram_) {
              hpTVsMVA[iy][ipt]->Fill(mva, pt);
              hetaVsMVA[iy][ipt]->Fill(mva, eta);
              hyVsMVA[iy][ipt]->Fill(mva, y);
              hVtxProbVsMVA[iy][ipt]->Fill(mva, VtxProb);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl_abs);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D_abs);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos2D);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl2D);
              hzDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva, dzos1);
              hxyDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva, dxyos1);
              hNHitD1VsMVA[iy][ipt]->Fill(mva, nhit1);
              hpTD1VsMVA[iy][ipt]->Fill(mva, pt1);
              hpTerrD1VsMVA[iy][ipt]->Fill(mva, ptErr1 / pt1);
              hEtaD1VsMVA[iy][ipt]->Fill(mva, eta1);
              hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva, H2dedx1);
              hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1, H2dedx1);
              hzDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva, dzos2);
              hxyDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva, dxyos2);
              hNHitD2VsMVA[iy][ipt]->Fill(mva, nhit2);
              hpTD2VsMVA[iy][ipt]->Fill(mva, pt2);
              hpTerrD2VsMVA[iy][ipt]->Fill(mva, ptErr2 / pt2);
              hEtaD2VsMVA[iy][ipt]->Fill(mva, eta2);
              hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva, H2dedx2);
              hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2, H2dedx2);
              if (threeProngDecay_) {
                hzDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva, dzos3);
                hxyDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva, dxyos3);
                hNHitD3VsMVA[iy][ipt]->Fill(mva, nhit3);
                hpTD3VsMVA[iy][ipt]->Fill(mva, pt3);
                hpTerrD3VsMVA[iy][ipt]->Fill(mva, ptErr3 / pt3);
                hEtaD3VsMVA[iy][ipt]->Fill(mva, eta3);
                hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva, H2dedx3);
                hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1, H2dedx3);
              }
            }
          }
        }
    }
  }
}

void VertexCompositeNtupleProducer::fillGEN(const edm::Event &iEvent,
                                            const edm::EventSetup &iSetup) {

  genVertex_ = reco::Vertex();

  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);
  // generated primary vertex information
  for (const auto &p : *genpars) {
    if (p.statusFlags().isLastCopy() &&
        (p.pdgId() == 21 || std::abs(p.pdgId()) <= 6)) {
      genVertex_ = p.vertex();
      break;
    }
  }
  gen_PVx_ = genVertex_.x();
  gen_PVy_ = genVertex_.y();
  gen_PVz_ = genVertex_.z();

  for (unsigned it = 0; it < genpars->size(); ++it) {

    const reco::GenParticle &trk = (*genpars)[it];

    int id = trk.pdgId();
    if (fabs(id) != PID_)
      continue; // check is target
    if (decayInGen_ && trk.numberOfDaughters() != 2)
      continue; // check 2-pron decay if target decays in Gen

    pt_gen = trk.pt();
    eta_gen = trk.eta();
    phi_gen = trk.phi();
    status_gen = trk.status();
    idmom = -77;
    y_gen = trk.rapidity();
    ptmom = -999.0;
    etamom = -999.0;
    phimom = -999.0;
    ymom = -999.0;
    statusmom = -999;

    if (trk.numberOfMothers() != 0) {
      const reco::Candidate *mom = trk.mother();
      idmom = mom->pdgId();
      ptmom = mom->pt();
      etamom = mom->eta();
      phimom = mom->phi();
      ymom = mom->rapidity();
      statusmom = mom->status();
    }

    if (!decayInGen_)
      continue;

    const reco::Candidate *Dd1 = trk.daughter(0);
    const reco::Candidate *Dd2 = trk.daughter(1);
    const reco::Candidate *Dd3 = trk.daughter(2);

    iddau1 = fabs(Dd1->pdgId());
    iddau2 = fabs(Dd2->pdgId());
    if (Dd3)
      iddau3 = fabs(Dd3->pdgId());

    getAncestorId(trk);
    genDecayLength(trk);

    gen_pT_ = trk.pt();
    gen_eta_ = trk.eta();
    gen_phi_ = trk.phi();
    gen_mass_ = trk.mass();
    gen_y_ = trk.rapidity();
    gen_charge_ = trk.charge();
    gen_pdgId_ = trk.pdgId();

    gen_pTD1_ = Dd1->pt();
    gen_etaD1_ = Dd1->eta();
    gen_phiD1_ = Dd1->phi();
    gen_massD1_ = Dd1->mass();
    gen_yD1_ = Dd1->rapidity();
    gen_chargeD1_ = Dd1->charge();
    gen_pdgIdD1_ = Dd1->pdgId();

    gen_pTD2_ = Dd2->pt();
    gen_etaD2_ = Dd2->eta();
    gen_phiD2_ = Dd2->phi();
    gen_massD2_ = Dd2->mass();
    gen_yD2_ = Dd2->rapidity();
    gen_chargeD2_ = Dd2->charge();
    gen_pdgIdD2_ = Dd2->pdgId();

    genCandidateNtuple->Fill();
  }
}

// ------------ method called once each job just before starting event
// loop  ------------
void VertexCompositeNtupleProducer::beginJob() {
  TH1D::SetDefaultSumw2();

  if (!doRecoNtuple_ && !doGenNtuple_) {
    cout << "No output for either RECO or GEN!! Fix config!!" << endl;
    return;
  }

  if (twoLayerDecay_ && doMuon_) {
    cout << "Muons cannot be coming from two layer decay!! Fix config!!"
         << endl;
    return;
  }

  if (saveHistogram_)
    initHistogram();
  if (saveTree_)
    initTree();
}

void VertexCompositeNtupleProducer::initHistogram() {
  for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++) {
    for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
      hMassVsMVA[iy][ipt] = fs->make<TH2F>(
          Form("hMassVsMVA_y%d_pt%d", iy, ipt), ";mva;mass(GeV)", 100, -1., 1.,
          massHistBins_, massHistPeak_ - massHistWidth_,
          massHistPeak_ + massHistWidth_);
      //   h3DDCAVsMVA[iy][ipt] =
      //   fs->make<TH2F>(Form("h3DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;3D
      //   DCA;",100,-1.,1.,1000,0,10); h2DDCAVsMVA[iy][ipt] =
      //   fs->make<TH2F>(Form("h2DDCAVsMVA_y%d_pt%d",iy,ipt),";mva;2D
      //   DCA;",100,-1.,1.,1000,0,10);

      if (saveAllHistogram_) {
        hpTVsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTVsMVA_y%d_pt%d", iy, ipt),
                                           ";mva;pT;", 100, -1, 1, 100, 0, 10);
        hetaVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hetaVsMVA_y%d_pt%d", iy, ipt), ";mva;eta;",
                           100, -1., 1., 40, -4, 4);
        hyVsMVA[iy][ipt] = fs->make<TH2F>(Form("hyVsMVA_y%d_pt%d", iy, ipt),
                                          ";mva;y;", 100, -1., 1., 40, -4, 4);
        hVtxProbVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hVtxProbVsMVA_y%d_pt%d", iy, ipt),
                           ";mva;VtxProb;", 100, -1., 1., 100, 0, 1);
        h3DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h3DCosPointingAngleVsMVA_y%d_pt%d", iy, ipt),
            ";mva;3DCosPointingAngle;", 100, -1., 1., 100, -1, 1);
        h3DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h3DPointingAngleVsMVA_y%d_pt%d", iy, ipt),
            ";mva;3DPointingAngle;", 100, -1., 1., 50, -3.14, 3.14);
        h2DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h2DCosPointingAngleVsMVA_y%d_pt%d", iy, ipt),
            ";mva;2DCosPointingAngle;", 100, -1., 1., 100, -1, 1);
        h2DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h2DPointingAngleVsMVA_y%d_pt%d", iy, ipt),
            ";mva;2DPointingAngle;", 100, -1., 1., 50, -3.14, 3.14);
        h3DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h3DDecayLengthSignificanceVsMVA_y%d_pt%d", iy, ipt),
            ";mva;3DDecayLengthSignificance;", 100, -1., 1., 300, 0, 30);
        h2DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(
            Form("h2DDecayLengthSignificanceVsMVA_y%d_pt%d", iy, ipt),
            ";mva;2DDecayLengthSignificance;", 100, -1., 1., 300, 0, 30);
        h3DDecayLengthVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("h3DDecayLengthVsMVA_y%d_pt%d", iy, ipt),
                           ";mva;3DDecayLength;", 100, -1., 1., 300, 0, 30);
        h2DDecayLengthVsMVA[iy][ipt] =
            fs->make<TH2F>(Form("h2DDecayLengthVsMVA_y%d_pt%d", iy, ipt),
                           ";mva;2DDecayLength;", 100, -1., 1., 300, 0, 30);
        hzDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(
            Form("hzDCASignificanceDaugther1VsMVA_y%d_pt%d", iy, ipt),
            ";mva;zDCASignificanceDaugther1;", 100, -1., 1., 100, -10, 10);
        hxyDCASignificanceDaugther1VsMVA[iy][ipt] = fs->make<TH2F>(
            Form("hxyDCASignificanceDaugther1VsMVA_y%d_pt%d", iy, ipt),
            ";mva;xyDCASignificanceDaugther1;", 100, -1., 1., 100, -10, 10);
        hNHitD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hNHitD1VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;NHitD1;", 100, -1., 1., 100, 0, 100);
        hpTD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTD1VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD1;",
                           100, -1., 1., 100, 0, 10);
        hpTerrD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTerrD1VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;pTerrD1;", 100, -1., 1., 50, 0, 0.5);
        hEtaD1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hEtaD1VsMVA_y%d_pt%d", iy, ipt), ";mva;EtaD1;",
                           100, -1., 1., 40, -4, 4);
        hdedxHarmonic2D1VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hdedxHarmonic2D1VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;dedxHarmonic2D1;", 100, -1., 1., 100, 0, 10);
        hdedxHarmonic2D1VsP[iy][ipt] =
            fs->make<TH2F>(Form("hdedxHarmonic2D1VsP_y%d_pt%d", iy, ipt),
                           ";p (GeV);dedxHarmonic2D1", 100, 0, 10, 100, 0, 10);
        hzDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(
            Form("hzDCASignificanceDaugther2VsMVA_y%d_pt%d", iy, ipt),
            ";mva;zDCASignificanceDaugther2;", 100, -1., 1., 100, -10, 10);
        hxyDCASignificanceDaugther2VsMVA[iy][ipt] = fs->make<TH2F>(
            Form("hxyDCASignificanceDaugther2VsMVA_y%d_pt%d", iy, ipt),
            ";mva;xyDCASignificanceDaugther2;", 100, -1., 1., 100, -10, 10);
        hNHitD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hNHitD2VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;NHitD2;", 100, -1., 1., 100, 0, 100);
        hpTD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTD2VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD2;",
                           100, -1., 1., 100, 0, 10);
        hpTerrD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hpTerrD2VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;pTerrD2;", 100, -1., 1., 50, 0, 0.5);
        hEtaD2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hEtaD2VsMVA_y%d_pt%d", iy, ipt), ";mva;EtaD2;",
                           100, -1., 1., 40, -4, 4);
        hdedxHarmonic2D2VsMVA[iy][ipt] =
            fs->make<TH2F>(Form("hdedxHarmonic2D2VsMVA_y%d_pt%d", iy, ipt),
                           ";mva;dedxHarmonic2D2;", 100, -1., 1., 100, 0, 10);
        hdedxHarmonic2D2VsP[iy][ipt] =
            fs->make<TH2F>(Form("hdedxHarmonic2D2VsP_y%d_pt%d", iy, ipt),
                           ";p (GeV);dedxHarmonic2D2", 100, 0, 10, 100, 0, 10);

        if (threeProngDecay_) {
          hzDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(
              Form("hzDCASignificanceDaugther3VsMVA_y%d_pt%d", iy, ipt),
              ";mva;zDCASignificanceDaugther3;", 100, -1., 1., 100, -10, 10);
          hxyDCASignificanceDaugther3VsMVA[iy][ipt] = fs->make<TH2F>(
              Form("hxyDCASignificanceDaugther3VsMVA_y%d_pt%d", iy, ipt),
              ";mva;xyDCASignificanceDaugther3;", 100, -1., 1., 100, -10, 10);
          hNHitD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hNHitD3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;NHitD3;", 100, -1., 1., 100, 0, 100);
          hpTD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hpTD3VsMVA_y%d_pt%d", iy, ipt), ";mva;pTD3;",
                             100, -1., 1., 100, 0, 10);
          hpTerrD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hpTerrD3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;pTerrD3;", 100, -1., 1., 50, 0, 0.5);
          hEtaD3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hEtaD3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;EtaD3;", 100, -1., 1., 40, -4, 4);
          hdedxHarmonic2D3VsMVA[iy][ipt] =
              fs->make<TH2F>(Form("hdedxHarmonic2D3VsMVA_y%d_pt%d", iy, ipt),
                             ";mva;dedxHarmonic2D3;", 100, -1., 1., 100, 0, 10);
          hdedxHarmonic2D3VsP[iy][ipt] = fs->make<TH2F>(
              Form("hdedxHarmonic2D3VsP_y%d_pt%d", iy, ipt),
              ";p (GeV);dedxHarmonic2D3", 100, 0, 10, 100, 0, 10);
        }
      }
    }
  }
}

void VertexCompositeNtupleProducer::initTree() {
  VertexCompositeNtuple =
      fs->make<TTree>("VertexCompositeNtuple", "VertexCompositeNtuple");
  genCandidateNtuple =
      fs->make<TTree>("genCandidateNtuple", "genCandidateNtuple");

  VertexCompositeNtuple->Branch("pT", &pt, "pT/F");
  VertexCompositeNtuple->Branch("y", &y, "y/F");
  VertexCompositeNtuple->Branch("eta", &eta, "eta/F");
  VertexCompositeNtuple->Branch("phi", &phi, "phi/F");
  VertexCompositeNtuple->Branch("mass", &mass, "mass/F");

  if (useAnyMVA_)
    VertexCompositeNtuple->Branch("mva", &mva, "mva/F");

  if (useDCA_) {
    VertexCompositeNtuple->Branch("dca3D", &dca3D, "dca3D/F");
    VertexCompositeNtuple->Branch("dcaErr3D", &dcaErr3D, "dcaErr3D/F");
  }

  if (isCentrality_)
    VertexCompositeNtuple->Branch("centrality", &centrality, "centrality/I");

  if (!isSkimMVA_) {
    // Event info
    VertexCompositeNtuple->Branch("Ntrkoffline", &Ntrkoffline, "Ntrkoffline/I");
    VertexCompositeNtuple->Branch("Npixel", &Npixel, "Npixel/I");
    VertexCompositeNtuple->Branch("HFsumETPlus", &HFsumETPlus, "HFsumETPlus/F");
    VertexCompositeNtuple->Branch("HFsumETMinus", &HFsumETMinus,
                                  "HFsumETMinus/F");
    VertexCompositeNtuple->Branch("ZDCPlus", &ZDCPlus, "ZDCPlus/F");
    VertexCompositeNtuple->Branch("ZDCMinus", &ZDCMinus, "ZDCMinus/F");
    VertexCompositeNtuple->Branch("bestvtxX", &bestvx, "bestvtxX/F");
    VertexCompositeNtuple->Branch("bestvtxY", &bestvy, "bestvtxY/F");
    VertexCompositeNtuple->Branch("bestvtxZ", &bestvz, "bestvtxZ/F");

    // Composite candidate info RECO
    VertexCompositeNtuple->Branch("flavor", &flavor, "flavor/F");
    VertexCompositeNtuple->Branch("VtxProb", &VtxProb, "VtxProb/F");
    //        VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2/F");
    //        VertexCompositeNtuple->Branch("VtxNDF",&ndf,"VtxNDF/F");
    VertexCompositeNtuple->Branch("3DCosPointingAngle", &agl,
                                  "3DCosPointingAngle/F");
    VertexCompositeNtuple->Branch("3DPointingAngle", &agl_abs,
                                  "3DPointingAngle/F");
    VertexCompositeNtuple->Branch("2DCosPointingAngle", &agl2D,
                                  "2DCosPointingAngle/F");
    VertexCompositeNtuple->Branch("2DPointingAngle", &agl2D_abs,
                                  "2DPointingAngle/F");
    VertexCompositeNtuple->Branch("3DDecayLengthSignificance", &dlos,
                                  "3DDecayLengthSignificance/F");
    VertexCompositeNtuple->Branch("3DDecayLength", &dl, "3DDecayLength/F");
    //        VertexCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError/F");
    VertexCompositeNtuple->Branch("2DDecayLengthSignificance", &dlos2D,
                                  "2DDecayLengthSignificance/F");
    VertexCompositeNtuple->Branch("2DDecayLength", &dl2D, "2DDecayLength/F");

    if (doGenMatching_) {
      VertexCompositeNtuple->Branch("isSwap", &isSwap, "isSwap/O");
      VertexCompositeNtuple->Branch("idmom_reco", &idmom_reco, "idmom_reco/I");
      VertexCompositeNtuple->Branch("matchGEN", &matchGEN, "matchGEN/O");

      if (!twoLayerDecay_ && !threeProngDecay_) {
        VertexCompositeNtuple->Branch("gen_ancestorFlavor",
                                      &gen_ancestorFlavor_);
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

    if (doGenMatchingTOF_) {
      VertexCompositeNtuple->Branch("PIDD1", &pid1, "PIDD1/I");
      VertexCompositeNtuple->Branch("PIDD2", &pid1, "PIDD2/I");
      VertexCompositeNtuple->Branch("TOFD1", &tof1, "TOFD1/F");
      VertexCompositeNtuple->Branch("TOFD2", &tof1, "TOFD2/F");
    }

    // daugther & grand daugther info
    if (twoLayerDecay_) {
      VertexCompositeNtuple->Branch("massDaugther1", &grand_mass,
                                    "massDaugther1/F");
      VertexCompositeNtuple->Branch("pTD1", &pt1, "pTD1/F");
      VertexCompositeNtuple->Branch("EtaD1", &eta1, "EtaD1/F");
      VertexCompositeNtuple->Branch("PhiD1", &phi1, "PhiD1/F");
      VertexCompositeNtuple->Branch("VtxProbDaugther1", &grand_VtxProb,
                                    "VtxProbDaugther1/F");
      VertexCompositeNtuple->Branch("VtxChi2Daugther1", &grand_vtxChi2,
                                    "VtxChi2Daugther1/F");
      VertexCompositeNtuple->Branch("VtxNDFDaugther1", &grand_ndf,
                                    "VtxNDFDaugther1/F");
      VertexCompositeNtuple->Branch("3DCosPointingAngleDaugther1", &grand_agl,
                                    "3DCosPointingAngleDaugther1/F");
      VertexCompositeNtuple->Branch("3DPointingAngleDaugther1", &grand_agl_abs,
                                    "3DPointingAngleDaugther1/F");
      VertexCompositeNtuple->Branch("2DCosPointingAngleDaugther1", &grand_agl2D,
                                    "2DCosPointingAngleDaugther1/F");
      VertexCompositeNtuple->Branch("2DPointingAngleDaugther1",
                                    &grand_agl2D_abs,
                                    "2DPointingAngleDaugther1/F");
      VertexCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",
                                    &grand_dlos,
                                    "3DDecayLengthSignificanceDaugther1/F");
      VertexCompositeNtuple->Branch("3DDecayLengthDaugther1", &grand_dl,
                                    "3DDecayLengthDaugther1/F");
      VertexCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",
                                    &grand_dlerror,
                                    "3DDecayLengthErrorDaugther1/F");
      VertexCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",
                                    &grand_dlos2D,
                                    "2DDecayLengthSignificanceDaugther1/F");
      VertexCompositeNtuple->Branch("zDCASignificanceDaugther2", &dzos2,
                                    "zDCASignificanceDaugther2/F");
      VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2", &dxyos2,
                                    "xyDCASignificanceDaugther2/F");
      VertexCompositeNtuple->Branch("NHitD2", &nhit2, "NHitD2/F");
      VertexCompositeNtuple->Branch("HighPurityDaugther2", &trkquality2,
                                    "HighPurityDaugther2/O");
      VertexCompositeNtuple->Branch("pTD2", &pt2, "pTD2/F");
      VertexCompositeNtuple->Branch("pTerrD2", &ptErr2, "pTerrD2/F");
      VertexCompositeNtuple->Branch("pD2", &p2, "pD2/F");
      VertexCompositeNtuple->Branch("EtaD2", &eta2, "EtaD2/F");
      VertexCompositeNtuple->Branch("PhiD2", &phi2, "PhiD2/F");
      VertexCompositeNtuple->Branch("chargeD2", &charge2, "chargeD2/I");
      VertexCompositeNtuple->Branch("dedxHarmonic2D2", &H2dedx2,
                                    "dedxHarmonic2D2/F");
      VertexCompositeNtuple->Branch("dedxTruncated40Daugther2", &T4dedx2,
                                    "dedxTruncated40Daugther2/F");
      VertexCompositeNtuple->Branch("normalizedChi2Daugther2", &trkChi2,
                                    "normalizedChi2Daugther2/F");
      VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther1",
                                    &grand_dzos1,
                                    "zDCASignificanceGrandDaugther1/F");
      VertexCompositeNtuple->Branch("zDCASignificanceGrandDaugther2",
                                    &grand_dzos2,
                                    "zDCASignificanceGrandDaugther2/F");
      VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther1",
                                    &grand_dxyos1,
                                    "xyDCASignificanceGrandDaugther1/F");
      VertexCompositeNtuple->Branch("xyDCASignificanceGrandDaugther2",
                                    &grand_dxyos2,
                                    "xyDCASignificanceGrandDaugther2/F");
      VertexCompositeNtuple->Branch("NHitGrandD1", &grand_nhit1,
                                    "NHitGrandD1/F");
      VertexCompositeNtuple->Branch("NHitGrandD2", &grand_nhit2,
                                    "NHitGrandD2/F");
      VertexCompositeNtuple->Branch("HighPurityGrandDaugther1",
                                    &grand_trkquality1,
                                    "HighPurityGrandDaugther1/O");
      VertexCompositeNtuple->Branch("HighPurityGrandDaugther2",
                                    &grand_trkquality2,
                                    "HighPurityGrandDaugther2/O");
      VertexCompositeNtuple->Branch("pTGrandD1", &grand_pt1, "pTGrandD1/F");
      VertexCompositeNtuple->Branch("pTGrandD2", &grand_pt2, "pTGrandD2/F");
      VertexCompositeNtuple->Branch("pTerrGrandD1", &grand_ptErr1,
                                    "pTerrGrandD1/F");
      VertexCompositeNtuple->Branch("pTerrGrandD2", &grand_ptErr2,
                                    "pTerrGrandD2/F");
      VertexCompositeNtuple->Branch("pGrandD1", &grand_p1, "pGrandD1/F");
      VertexCompositeNtuple->Branch("pGrandD2", &grand_p2, "pGrandD2/F");
      VertexCompositeNtuple->Branch("EtaGrandD1", &grand_eta1, "EtaGrandD1/F");
      VertexCompositeNtuple->Branch("EtaGrandD2", &grand_eta2, "EtaGrandD2/F");
      VertexCompositeNtuple->Branch("chargeGrandD1", &grand_charge1,
                                    "chargeGrandD1/I");
      VertexCompositeNtuple->Branch("chargeGrandD2", &grand_charge2,
                                    "chargeGrandD2/I");
      VertexCompositeNtuple->Branch("dedxHarmonic2GrandD1", &grand_H2dedx1,
                                    "dedxHarmonic2GrandD1/F");
      VertexCompositeNtuple->Branch("dedxHarmonic2GrandD2", &grand_H2dedx2,
                                    "dedxHarmonic2GrandD2/F");
      VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther1",
                                    &grand_T4dedx1,
                                    "dedxTruncated40GrandDaugther1/F");
      VertexCompositeNtuple->Branch("dedxTruncated40GrandDaugther2",
                                    &grand_T4dedx2,
                                    "dedxTruncated40GrandDaugther2/F");
      VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther1",
                                    &grand_trkChi1,
                                    "normalizedChi2GrandDaugther1/F");
      VertexCompositeNtuple->Branch("normalizedChi2GrandDaugther2",
                                    &grand_trkChi2,
                                    "normalizedChi2GrandDaugther2/F");
    } else {
      VertexCompositeNtuple->Branch("zDCASignificanceDaugther1", &dzos1,
                                    "zDCASignificanceDaugther1/F");
      VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1", &dxyos1,
                                    "xyDCASignificanceDaugther1/F");
      VertexCompositeNtuple->Branch("NHitD1", &nhit1, "NHitD1/F");
      VertexCompositeNtuple->Branch("HighPurityDaugther1", &trkquality1,
                                    "HighPurityDaugther1/O");
      VertexCompositeNtuple->Branch("pTD1", &pt1, "pTD1/F");
      VertexCompositeNtuple->Branch("pTerrD1", &ptErr1, "pTerrD1/F");
      VertexCompositeNtuple->Branch("pD1", &p1, "pD1/F");
      VertexCompositeNtuple->Branch("EtaD1", &eta1, "EtaD1/F");
      VertexCompositeNtuple->Branch("PhiD1", &eta1, "PhiD1/F");
      VertexCompositeNtuple->Branch("chargeD1", &charge1, "chargeD1/I");
      VertexCompositeNtuple->Branch("dedxHarmonic2D1", &H2dedx1,
                                    "dedxHarmonic2D1/F");
      VertexCompositeNtuple->Branch("dedxTruncated40Daugther1", &T4dedx1,
                                    "dedxTruncated40Daugther1/F");
      VertexCompositeNtuple->Branch("normalizedChi2Daugther1", &trkChi1,
                                    "normalizedChi2Daugther1/F");
      VertexCompositeNtuple->Branch("zDCASignificanceDaugther2", &dzos2,
                                    "zDCASignificanceDaugther2/F");
      VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2", &dxyos2,
                                    "xyDCASignificanceDaugther2/F");
      VertexCompositeNtuple->Branch("NHitD2", &nhit2, "NHitD2/F");
      VertexCompositeNtuple->Branch("HighPurityDaugther2", &trkquality2,
                                    "HighPurityDaugther2/O");
      VertexCompositeNtuple->Branch("pTD2", &pt2, "pTD2/F");
      VertexCompositeNtuple->Branch("pTerrD2", &ptErr2, "pTerrD2/F");
      VertexCompositeNtuple->Branch("pD2", &p2, "pD2/F");
      VertexCompositeNtuple->Branch("EtaD2", &eta2, "EtaD2/F");
      VertexCompositeNtuple->Branch("PhiD2", &eta2, "PhiD2/F");
      VertexCompositeNtuple->Branch("chargeD2", &charge2, "chargeD2/I");
      VertexCompositeNtuple->Branch("dedxHarmonic2D2", &H2dedx2,
                                    "dedxHarmonic2D2/F");
      VertexCompositeNtuple->Branch("dedxTruncated40Daugther2", &T4dedx2,
                                    "dedxTruncated40Daugther2/F");
      VertexCompositeNtuple->Branch("normalizedChi2Daugther2", &trkChi2,
                                    "normalizedChi2Daugther2/F");
      if (threeProngDecay_) {
        VertexCompositeNtuple->Branch("zDCASignificanceDaugther3", &dzos3, "zDCASignificanceDaugther3/F");
        VertexCompositeNtuple->Branch("xyDCASignificanceDaugther3", &dxyos3, "xyDCASignificanceDaugther3/F");
        VertexCompositeNtuple->Branch("NHitD3", &nhit3, "NHitD3/F");
        VertexCompositeNtuple->Branch("HighPurityDaugther3", &trkquality3, "HighPurityDaugther3/O");
        VertexCompositeNtuple->Branch("pTD3", &pt1, "pTD3/F");
        VertexCompositeNtuple->Branch("pTerrD3", &ptErr3, "pTerrD3/F");
        VertexCompositeNtuple->Branch("EtaD3", &eta1, "EtaD3/F");
        VertexCompositeNtuple->Branch("dedxHarmonic2D3", &H2dedx1, "dedxHarmonic2D3/F");
      }
    }

    if (doMuon_) {
      VertexCompositeNtuple->Branch("OneStMuon1", &onestmuon1, "OneStMuon1/O");
      VertexCompositeNtuple->Branch("OneStMuon2", &onestmuon2, "OneStMuon2/O");
      VertexCompositeNtuple->Branch("PFMuon1", &pfmuon1, "PFMuon1/O");
      VertexCompositeNtuple->Branch("PFMuon2", &pfmuon2, "PFMuon2/O");
      VertexCompositeNtuple->Branch("GlbMuon1", &glbmuon1, "GlbMuon1/O");
      VertexCompositeNtuple->Branch("GlbMuon2", &glbmuon2, "GlbMuon2/O");
      VertexCompositeNtuple->Branch("trkMuon1", &trkmuon1, "trkMuon1/O");
      VertexCompositeNtuple->Branch("trkMuon2", &trkmuon2, "trkMuon2/O");
      VertexCompositeNtuple->Branch("caloMuon1", &calomuon1, "caloMuon1/O");
      VertexCompositeNtuple->Branch("caloMuon2", &calomuon2, "caloMuon2/O");
      VertexCompositeNtuple->Branch("SoftMuon1", &softmuon1, "SoftMuon1/O");
      VertexCompositeNtuple->Branch("SoftMuon2", &softmuon2, "SoftMuon2/O");
      if (doMuonFull_) {
        VertexCompositeNtuple->Branch("nMatchedChamberD1", &nmatchedch1, "nMatchedChamberD1/F");
        VertexCompositeNtuple->Branch("nMatchedStationD1", &nmatchedst1, "nMatchedStationD1/F");
        VertexCompositeNtuple->Branch("EnergyDepositionD1", &matchedenergy1, "EnergyDepositionD1/F");
        VertexCompositeNtuple->Branch("nMatchedChamberD2", &nmatchedch2, "nMatchedChamberD2/F");
        VertexCompositeNtuple->Branch("nMatchedStationD2", &nmatchedst2, "nMatchedStationD2/F");
        VertexCompositeNtuple->Branch("EnergyDepositionD2", &matchedenergy2, "EnergyDepositionD2/F");
        VertexCompositeNtuple->Branch("dx1_seg", &dx1_seg_, "dx1_seg/F");
        VertexCompositeNtuple->Branch("dy1_seg", &dy1_seg_, "dy1_seg/F");
        VertexCompositeNtuple->Branch("dxSig1_seg", &dxSig1_seg_, "dxSig1_seg/F");
        VertexCompositeNtuple->Branch("dySig1_seg", &dySig1_seg_, "dySig1_seg/F");
        VertexCompositeNtuple->Branch("ddxdz1_seg", &ddxdz1_seg_, "ddxdz1_seg/F");
        VertexCompositeNtuple->Branch("ddydz1_seg", &ddydz1_seg_, "ddydz1_seg/F");
        VertexCompositeNtuple->Branch("ddxdzSig1_seg", &ddxdzSig1_seg_, "ddxdzSig1_seg/F");
        VertexCompositeNtuple->Branch("ddydzSig1_seg", &ddydzSig1_seg_, "ddydzSig1_seg/F");
        VertexCompositeNtuple->Branch("dx2_seg", &dx2_seg_, "dx2_seg/F");
        VertexCompositeNtuple->Branch("dy2_seg", &dy2_seg_, "dy2_seg/F");
        VertexCompositeNtuple->Branch("dxSig2_seg", &dxSig2_seg_, "dxSig2_seg/F");
        VertexCompositeNtuple->Branch("dySig2_seg", &dySig2_seg_, "dySig2_seg/F");
        VertexCompositeNtuple->Branch("ddxdz2_seg", &ddxdz2_seg_, "ddxdz2_seg/F");
        VertexCompositeNtuple->Branch("ddydz2_seg", &ddydz2_seg_, "ddydz2_seg/F");
        VertexCompositeNtuple->Branch("ddxdzSig2_seg", &ddxdzSig2_seg_, "ddxdzSig2_seg/F");
        VertexCompositeNtuple->Branch("ddydzSig2_seg", &ddydzSig2_seg_, "ddydzSig2_seg/F");
      }
    }
  }

  if (doGenNtuple_) {
    genCandidateNtuple->Branch("pT_gen", &pt_gen, "pT_gen/F");
    genCandidateNtuple->Branch("eta_gen", &eta_gen, "eta_gen/F");
    genCandidateNtuple->Branch("y_gen", &y_gen, "y_gen/F");
    genCandidateNtuple->Branch("status_gen", &status_gen, "status_gen/I");
    genCandidateNtuple->Branch("MotherID_gen", &idmom, "MotherID_gen/I");
    genCandidateNtuple->Branch("MotherPt_gen", &ptmom, "MotherPt_gen/F");
    genCandidateNtuple->Branch("MotherEta_gen", &etamom, "MotherEta_gen/F");
    genCandidateNtuple->Branch("MotherPhi_gen", &phimom, "MotherPhi_gen/F");
    genCandidateNtuple->Branch("MotherY_gen", &ymom, "MotherY_gen/F");
    genCandidateNtuple->Branch("MotherStatus_gen", &statusmom,
                               "MotherStatus_gen/I");

    if (decayInGen_) {
      genCandidateNtuple->Branch("DauID1_gen", &iddau1, "DauID1_gen/I");
      genCandidateNtuple->Branch("DauID2_gen", &iddau2, "DauID2_gen/I");
      genCandidateNtuple->Branch("DauID3_gen", &iddau3, "DauID3_gen/I");
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
}

int VertexCompositeNtupleProducer::muAssocToTrack(
    const reco::TrackRef &trackref,
    const edm::Handle<reco::MuonCollection> &muonh) const {
  auto muon =
      std::find_if(muonh->cbegin(), muonh->cend(), [&](const reco::Muon &m) {
        return (m.track().isNonnull() && m.track() == trackref);
      });
  return (muon != muonh->cend() ? std::distance(muonh->cbegin(), muon) : -1);
}

// ------------ method called once each job just after ending the event
// loop  ------------
void VertexCompositeNtupleProducer::endJob() {}

reco::GenParticleRef VertexCompositeNtupleProducer::findMother(
    const reco::GenParticleRef &genParRef) {
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

void VertexCompositeNtupleProducer::genDecayLength(
    const reco::GenParticle &gCand) {
  gen_decayLength2D_ = -99.;
  gen_decayLength3D_ = -99.;
  gen_angle2D_ = -99;
  gen_angle3D_ = -99;

  if (gCand.numberOfDaughters() == 0 || !gCand.daughter(0))
    return;
  const auto &dauVtx = gCand.daughter(0)->vertex();
  TVector3 ptosvec(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(),
                   dauVtx.Z() - genVertex_.z());
  TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
  gen_angle3D_ = secvec.Angle(ptosvec);
  gen_decayLength3D_ = ptosvec.Mag();
  TVector3 ptosvec2D(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(),
                     0.0);
  TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
  gen_angle2D_ = secvec2D.Angle(ptosvec2D);
  gen_decayLength2D_ = ptosvec2D.Mag();
}

void VertexCompositeNtupleProducer::getAncestorId(
    const reco::GenParticle &gCand) {
  gen_ancestorId_ = 0;
  gen_ancestorFlavor_ = 0;
  for (auto mothers = gCand.motherRefVector(); !mothers.empty();) {
    auto mom = mothers.at(0);
    mothers = mom->motherRefVector();
    gen_ancestorId_ = mom->pdgId();
    const auto idstr = std::to_string(std::abs(gen_ancestorId_));
    gen_ancestorFlavor_ =
        std::stoi(std::string{idstr.begin(), idstr.begin() + 1});
    if (idstr[0] == '5') {
      break;
    }
    if (std::abs(gen_ancestorId_) <= 40)
      break;
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeNtupleProducer);
