
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
  virtual void clearVectors();
  void genDecayLength(const uint&, const reco::GenParticle&);

  bool matchHadron(const reco::Candidate* _dmeson_, const reco::GenParticle& _gen_, bool isMatchD0) const;
  bool matchHadron(const reco::Candidate* _dmeson_, const reco::Candidate& _gen_, bool isMatchD0) const;
  bool matchTrackdR(const reco::Candidate* _recoTrk_, const reco::Candidate* _genTrk_, bool chkchrg) const;
  bool checkSwap(const reco::Candidate* _dmeson_, const reco::GenParticle& _gen_) const;
  bool checkSwap(const reco::Candidate* _dmeson_, const reco::Candidate& _gen_) const;

  int muAssocToTrack( const reco::TrackRef& trackref, const edm::Handle<reco::MuonCollection>& muonh) const;

  reco::GenParticleRef findMother(const reco::GenParticleRef&);
  void genDecayLength(const reco::Candidate& gCand, float& gen_decayLength2D_, float& gen_decayLength3D_, float& gen_angle2D_, float& gen_angle3D_);
  void getAncestorId(const reco::Candidate& gCand, int& gen_ancestorId_, int& gen_ancestorFlavor_ );

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
    std::vector<float> mva;
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> flavor;
    std::vector<float> y;
    std::vector<float> mass;
    std::vector<float> VtxProb;
    std::vector<float> dlos;
    std::vector<float> dl;
    std::vector<float> dlerror;
    std::vector<float> agl;
    std::vector<float> vtxChi2;
    std::vector<float> ndf;
    std::vector<float> agl_abs;
    std::vector<float> agl2D;
    std::vector<float> agl2D_abs;
    std::vector<float> dlos2D;
    std::vector<float> dl2D;
    std::vector<bool> isSwap;
    std::vector<bool> matchGEN;
    std::vector<int> idBAnc_reco;
    std::vector<int> pionFlavor;
    std::vector<int> idmom_reco;
    std::vector<float> gen_agl_abs;
    std::vector<float> gen_agl2D_abs;
    std::vector<float> gen_dl;
    std::vector<float> gen_dl2D;

    //dau candidate info
    std::vector<float> grand_mass;
    std::vector<float> grand_VtxProb;
    std::vector<float> grand_dlos;
    std::vector<float> grand_dl;
    std::vector<float> grand_dlerror;
    std::vector<float> grand_agl;
    std::vector<float> grand_vtxChi2;
    std::vector<float> grand_ndf;
    std::vector<float> grand_agl_abs;
    std::vector<float> grand_agl2D;
    std::vector<float> grand_agl2D_abs;
    std::vector<float> grand_dlos2D;

    //dau info
    std::vector<float> dzos1;
    std::vector<float> dzos2;
    std::vector<float> dzos3;
    std::vector<float> dxyos1;
    std::vector<float> dxyos2;
    std::vector<float> dxyos3;
    std::vector<float> nhit1;
    std::vector<float> nhit2;
    std::vector<float> nhit3;
    std::vector<bool> trkquality1;
    std::vector<bool> trkquality2;
    std::vector<bool> trkquality3;
    std::vector<float> pt1;
    std::vector<float> pt2;
    std::vector<float> pt3;
    std::vector<float> ptErr1;
    std::vector<float> ptErr2;
    std::vector<float> ptErr3;
    std::vector<float> p1;
    std::vector<float> p2;
    std::vector<float> p3;
    std::vector<float> eta1;
    std::vector<float> eta2;
    std::vector<float> eta3;
    std::vector<float> phi1;
    std::vector<float> phi2;
    std::vector<float> phi3;
    std::vector<int> charge1;
    std::vector<int> charge2;
    std::vector<int> charge3;
    std::vector<int> pid1;
    std::vector<int> pid2;
    std::vector<int> pid3;
    std::vector<float> tof1;
    std::vector<float> tof2;
    std::vector<float> tof3;
    std::vector<float> H2dedx1;
    std::vector<float> H2dedx2;
    std::vector<float> H2dedx3;
    std::vector<float> T4dedx1;
    std::vector<float> T4dedx2;
    std::vector<float> T4dedx3;
    std::vector<float> trkChi1;
    std::vector<float> trkChi2;
    std::vector<float> trkChi3;

    //grand-dau info
    std::vector<float> grand_dzos1;
    std::vector<float> grand_dzos2;
    std::vector<float> grand_dxyos1;
    std::vector<float> grand_dxyos2;
    std::vector<float> grand_nhit1;
    std::vector<float> grand_nhit2;
    std::vector<bool> grand_trkquality1;
    std::vector<bool> grand_trkquality2;
    std::vector<float> grand_pt1;
    std::vector<float> grand_pt2;
    std::vector<float> grand_ptErr1;
    std::vector<float> grand_ptErr2;
    std::vector<float> grand_p1;
    std::vector<float> grand_p2;
    std::vector<float> grand_eta1;
    std::vector<float> grand_eta2;
    std::vector<int> grand_charge1;
    std::vector<int> grand_charge2;
    std::vector<float> grand_H2dedx1;
    std::vector<float> grand_H2dedx2;
    std::vector<float> grand_T4dedx1;
    std::vector<float> grand_T4dedx2;
    std::vector<float> grand_trkChi1;
    std::vector<float> grand_trkChi2;

    //dau muon info
    std::vector<bool> onestmuon1;
    std::vector<bool> onestmuon2;
    std::vector<bool> pfmuon1;
    std::vector<bool> pfmuon2;
    std::vector<bool> glbmuon1;
    std::vector<bool> glbmuon2;
    std::vector<bool> trkmuon1;
    std::vector<bool> trkmuon2;
    std::vector<bool> calomuon1;
    std::vector<bool> calomuon2;
    std::vector<bool> softmuon1;
    std::vector<bool> softmuon2;
    std::vector<float> nmatchedst1;
    std::vector<float> nmatchedch1;
    std::vector<float> ntrackerlayer1;
    std::vector<float> npixellayer1;
    std::vector<float> matchedenergy1;
    std::vector<float> nmatchedst2;
    std::vector<float> nmatchedch2;
    std::vector<float> ntrackerlayer2;
    std::vector<float> npixellayer2;
    std::vector<float> matchedenergy2;
    std::vector<float> dx1_seg_;
    std::vector<float> dy1_seg_;
    std::vector<float> dxSig1_seg_;
    std::vector<float> dySig1_seg_;
    std::vector<float> ddxdz1_seg_;
    std::vector<float> ddydz1_seg_;
    std::vector<float> ddxdzSig1_seg_;
    std::vector<float> ddydzSig1_seg_;
    std::vector<float> dx2_seg_;
    std::vector<float> dy2_seg_;
    std::vector<float> dxSig2_seg_;
    std::vector<float> dySig2_seg_;
    std::vector<float> ddxdz2_seg_;
    std::vector<float> ddydz2_seg_;
    std::vector<float> ddxdzSig2_seg_;
    std::vector<float> ddydzSig2_seg_;

    // gen info    
    int candSize_gen;
    std::vector<float> mass_gen;
    std::vector<float> pt_gen;
    std::vector<float> eta_gen;
    std::vector<float> phi_gen;
    std::vector<int> status_gen;
    std::vector<int> idmom;
    std::vector<float> y_gen;
    std::vector<int> iddau1;
    std::vector<int> iddau2;
    std::vector<int> iddau3;

    std::vector<float> matchGen_D0pT_;
    std::vector<float> matchGen_D0eta_;
    std::vector<float> matchGen_D0phi_;
    std::vector<float> matchGen_D0mass_;
    std::vector<float> matchGen_D0y_;
    std::vector<int> matchGen_D0charge_;
    std::vector<int> matchGen_D0pdgId_;

    std::vector<float> matchGen_D0Dau1_pT_;
    std::vector<float> matchGen_D0Dau1_eta_;
    std::vector<float> matchGen_D0Dau1_phi_;
    std::vector<float> matchGen_D0Dau1_mass_;
    std::vector<float> matchGen_D0Dau1_y_;
    std::vector<int> matchGen_D0Dau1_charge_;
    std::vector<int> matchGen_D0Dau1_pdgId_;

    std::vector<float> matchGen_D0Dau2_pT_;
    std::vector<float> matchGen_D0Dau2_eta_;
    std::vector<float> matchGen_D0Dau2_phi_;
    std::vector<float> matchGen_D0Dau2_mass_;
    std::vector<float> matchGen_D0Dau2_y_;
    std::vector<int> matchGen_D0Dau2_charge_;
    std::vector<int> matchGen_D0Dau2_pdgId_;

    std::vector<float> matchGen_D1pT_;
    std::vector<float> matchGen_D1eta_;
    std::vector<float> matchGen_D1phi_;
    std::vector<float> matchGen_D1mass_;
    std::vector<float> matchGen_D1y_;
    std::vector<float> matchGen_D1decayLength2D_;
    std::vector<float> matchGen_D1decayLength3D_;
    std::vector<float> matchGen_D1angle2D_;
    std::vector<float> matchGen_D1angle3D_;
    std::vector<int> matchGen_D1ancestorId_;
    std::vector<int> matchGen_D1ancestorFlavor_;
    std::vector<int> matchGen_D1charge_;
    std::vector<int> matchGen_D1pdgId_;

    std::vector<float> gen_D0pT_;
    std::vector<float> gen_D0eta_;
    std::vector<float> gen_D0phi_;
    std::vector<float> gen_D0mass_;
    std::vector<float> gen_D0y_;
    std::vector<int> gen_D0charge_;
    std::vector<int> gen_D0pdgId_;

    std::vector<float> gen_D0Dau1_pT_;
    std::vector<float> gen_D0Dau1_eta_;
    std::vector<float> gen_D0Dau1_phi_;
    std::vector<float> gen_D0Dau1_mass_;
    std::vector<float> gen_D0Dau1_y_;
    std::vector<int> gen_D0Dau1_charge_;
    std::vector<int> gen_D0Dau1_pdgId_;

    std::vector<float> gen_D0Dau2_pT_;
    std::vector<float> gen_D0Dau2_eta_;
    std::vector<float> gen_D0Dau2_phi_;
    std::vector<float> gen_D0Dau2_mass_;
    std::vector<float> gen_D0Dau2_y_;
    std::vector<int> gen_D0Dau2_charge_;
    std::vector<int> gen_D0Dau2_pdgId_;

    std::vector<float> gen_D1pT_;
    std::vector<float> gen_D1eta_;
    std::vector<float> gen_D1phi_;
    std::vector<float> gen_D1mass_;
    std::vector<float> gen_D1y_;
    std::vector<int> gen_D1charge_;
    std::vector<int> gen_D1pdgId_;

    //vector for gen match
    // vector< vector<double> > *pVect;
    // vector< vector<double> > *gpVect;
    // vector<double> *Dvector1;
    // vector<double> *GDvector1;
    // vector<double> *Dvector2;
    // vector<double> *GDvector2;
    // vector<double> *Dvector3;
    // vector<int> *pVectIDmom;
    
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



bool VertexCompositeTreeProducer2::matchHadron(const reco::Candidate* _dmeson_, const reco::GenParticle& _gen_, bool isMatchD0) const {
  bool match = false;
  if(isMatchD0){
    reco::Candidate const* reco_trk1 = _dmeson_->daughter(0);
    reco::Candidate const* reco_trk2 = _dmeson_->daughter(1);

    reco::Candidate const* gen_trk1 = _gen_.daughter(0);
    reco::Candidate const* gen_trk2 = _gen_.daughter(1);

    if( matchTrackdR(reco_trk1, gen_trk1, true)){
        if( matchTrackdR(reco_trk2, gen_trk2,true)) {
            match = true;
            return match;
        }
    }    
    if( matchTrackdR(reco_trk2, gen_trk1, true)){
        if( matchTrackdR(reco_trk1, gen_trk2,true)) {
            match = true;
            return match;
        }
    }    
  }
  if(!isMatchD0){
    if(matchTrackdR(_dmeson_, &_gen_,true)) match = true;
  }
  return match;
};
bool VertexCompositeTreeProducer2::matchHadron(const reco::Candidate* _dmeson_, const reco::Candidate& _gen_, bool isMatchD0) const {
  bool match = false;
  if(isMatchD0){
    reco::Candidate const* reco_trk1 = _dmeson_->daughter(0);
    reco::Candidate const* reco_trk2 = _dmeson_->daughter(1);

    reco::Candidate const* gen_trk1 = _gen_.daughter(0);
    reco::Candidate const* gen_trk2 = _gen_.daughter(1);

    if( matchTrackdR(reco_trk1, gen_trk1, true)){
        if( matchTrackdR(reco_trk2, gen_trk2,true)) {
            match = true;
            return match;
        }
    }    
    if( matchTrackdR(reco_trk2, gen_trk1, true)){
        if( matchTrackdR(reco_trk1, gen_trk2,true)) {
            match = true;
            return match;
        }
    }    
  }
  if(!isMatchD0){
    if(matchTrackdR(_dmeson_, &_gen_,true)) match = true;
  }
  return match;
};

bool VertexCompositeTreeProducer2::checkSwap(const reco::Candidate* _dmeson_, const reco::GenParticle& _gen_) const {
    cout <<"_dmeson_ pdg ID : "<<_dmeson_->pdgId() << endl;
    cout <<"_gen_ pdgId : " << _gen_.pdgId() << endl;
    return _dmeson_->pdgId() != _gen_.pdgId();
};
bool VertexCompositeTreeProducer2::checkSwap(const reco::Candidate* _dmeson_, const reco::Candidate& _gen_) const {
    cout <<"_dmeson_ pdg ID : "<<_dmeson_->pdgId() << endl;
    cout <<"_gen_ pdgId : " << _gen_.pdgId() << endl;
    return _dmeson_->pdgId() != _gen_.pdgId();
};

bool VertexCompositeTreeProducer2::matchTrackdR(const reco::Candidate* _recoTrk_, const reco::Candidate* _genTrk_, bool chkchrg=true) const {
    bool pass= false;
    // _deltaR_
    if(chkchrg && (_recoTrk_->charge() != _genTrk_->charge())) return false;
    const double dR = reco::deltaR(*_recoTrk_, *_genTrk_);
    if(dR < deltaR_) pass = true;
    return pass;
};


reco::GenParticleRef VertexCompositeTreeProducer2::findMother(const reco::GenParticleRef& genParRef)
{
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

void VertexCompositeTreeProducer2::genDecayLength(const reco::Candidate& gCand, vector<float> &gen_decayLength2D_, vector<float> &gen_decayLength3D_, vector<float> &gen_angle2D_, vector<float> &gen_angle3D_){
  _gen_decayLength2D_ = -99.;
  _gen_decayLength3D_ = -99.;
  _gen_angle2D_ = -99;
  _gen_angle3D_ = -99;

  if(gCand.numberOfDaughters()==0 || !gCand.daughter(0)){
    gen_decayLength2D_.push_back(_gen_decayLength2D_);
    gen_decayLength3D_.push_back(_gen_decayLength3D_);
    gen_angle2D_.push_back(_gen_angle2D_);
    gen_angle3D_.push_back(_gen_angle3D_);
    return;
  }
  const auto& dauVtx = gCand.daughter(0)->vertex();
  const auto& genVertex_ = gCand.vertex();
  TVector3 ptosvec(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), dauVtx.Z() - genVertex_.z());
  TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
  gen_angle3D_.push_back(secvec.Angle(ptosvec));
  gen_decayLength3D_.push_back(ptosvec.Mag());
  TVector3 ptosvec2D(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(), 0.0);
  TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
  gen_angle2D_.push_back(secvec2D.Angle(ptosvec2D));
  gen_decayLength2D_.push_back(ptosvec2D.Mag());
};

void VertexCompositeTreeProducer2::getAncestorId(const reco::Candidate& gCand, vector<int> &gen_ancestorId_, vector<int> &gen_ancestorFlavor_){
  _gen_ancestorId_ = 0;
  _gen_ancestorFlavor_ = 0;
//  reco::GenParticle gCand1(gCand.charge(),gCand.p4(),gCand.vertex(),421,2,true);
  //for (auto mothers = gCand.motherRefVector();
  //    !mothers.empty(); ) {
  //  auto mom = mothers.at(0);
  //  mothers = mom->motherRefVector();
  //  gen_ancestorId_ = mom->pdgId();
  //  cout << "gen_ancestorId_ : " << gen_ancestorId_ << endl;
  //  const auto idstr = std::to_string(std::abs(gen_ancestorId_));
  //  gen_ancestorFlavor_ = std::stoi(std::string{idstr.begin(), idstr.begin()+1});
  //  cout << "gen_ancestorFlavor_ : " << gen_ancestorFlavor_ << endl;
  //  if (idstr[0] == '5') {
  //    break;
  //  }
  //  if (std::abs(gen_ancestorId_) <= 40) break;
  //}
  if((mom==nullptr)){
    gen_ancestorId_.push_back(_gen_ancestorId_);
    gen_ancestorFlavor_.push_back(_gen_ancestorFlavor_);
    return;
  }
  for (auto mom = gCand.mother(); !(mom==nullptr);){
          gen_ancestorId_.push_back(mom->pdgId());
    const auto idstr = std::to_string(std::abs(gen_ancestorId_));
    gen_ancestorFlavor_.push_back(std::stoi(std::string{idstr.begin(), idstr.begin()+1}));
    if (idstr[0] == '5') {
      break;
    }
    if (std::abs(gen_ancestorId_) <= 40) break;
          mom = mom->mother();
  }

          
};
void VertexCompositeTreeProducer2::clearVectors() {
  // Clear event info vectors
  mva.clear();
  pt.clear();
  eta.clear();
  phi.clear();
  flavor.clear();
  y.clear();
  mass.clear();
  VtxProb.clear();
  dlos.clear();
  dl.clear();
  dlerror.clear();
  agl.clear();
  vtxChi2.clear();
  ndf.clear();
  agl_abs.clear();
  agl2D.clear();
  agl2D_abs.clear();
  dlos2D.clear();
  dl2D.clear();
  isSwap.clear();
  matchGEN.clear();
  idBAnc_reco.clear();
  pionFlavor.clear();
  idmom_reco.clear();
  gen_agl_abs.clear();
  gen_agl2D_abs.clear();
  gen_dl.clear();
  gen_dl2D.clear();

  // Clear dau candidate info vectors
  grand_mass.clear();
  grand_VtxProb.clear();
  grand_dlos.clear();
  grand_dl.clear();
  grand_dlerror.clear();
  grand_agl.clear();
  grand_vtxChi2.clear();
  grand_ndf.clear();
  grand_agl_abs.clear();
  grand_agl2D.clear();
  grand_agl2D_abs.clear();
  grand_dlos2D.clear();

  // Clear dau info vectors
  dzos1.clear();
  dzos2.clear();
  dzos3.clear();
  dxyos1.clear();
  dxyos2.clear();
  dxyos3.clear();
  nhit1.clear();
  nhit2.clear();
  nhit3.clear();
  trkquality1.clear();
  trkquality2.clear();
  trkquality3.clear();
  pt1.clear();
  pt2.clear();
  pt3.clear();
  ptErr1.clear();
  ptErr2.clear();
  ptErr3.clear();
  p1.clear();
  p2.clear();
  p3.clear();
  eta1.clear();
  eta2.clear();
  eta3.clear();
  phi1.clear();
  phi2.clear();
  phi3.clear();
  charge1.clear();
  charge2.clear();
  charge3.clear();
  pid1.clear();
  pid2.clear();
  pid3.clear();
  tof1.clear();
  tof2.clear();
  tof3.clear();
  H2dedx1.clear();
  H2dedx2.clear();
  H2dedx3.clear();
  T4dedx1.clear();
  T4dedx2.clear();
  T4dedx3.clear();
  trkChi1.clear();
  trkChi2.clear();
  trkChi3.clear();

  // Clear grand-dau info vectors
  grand_dzos1.clear();
  grand_dzos2.clear();
  grand_dxyos1.clear();
  grand_dxyos2.clear();
  grand_nhit1.clear();
  grand_nhit2.clear();
  grand_trkquality1.clear();
  grand_trkquality2.clear();
  grand_pt1.clear();
  grand_pt2.clear();
  grand_ptErr1.clear();
  grand_ptErr2.clear();
  grand_p1.clear();
  grand_p2.clear();
  grand_eta1.clear();
  grand_eta2.clear();
  grand_charge1.clear();
  grand_charge2.clear();
  grand_H2dedx1.clear();
  grand_H2dedx2.clear();
  grand_T4dedx1.clear();
  grand_T4dedx2.clear();
  grand_trkChi1.clear();
  grand_trkChi2.clear();

  // Clear dau muon info vectors
  onestmuon1.clear();
  onestmuon2.clear();
  pfmuon1.clear();
  pfmuon2.clear();
  glbmuon1.clear();
  glbmuon2.clear();
  trkmuon1.clear();
  trkmuon2.clear();
  calomuon1.clear();
  calomuon2.clear();
  softmuon1.clear();
  softmuon2.clear();
  nmatchedst1.clear();
  nmatchedch1.clear();
  ntrackerlayer1.clear();
  npixellayer1.clear();
  matchedenergy1.clear();
  nmatchedst2.clear();
  nmatchedch2.clear();
  ntrackerlayer2.clear();
  npixellayer2.clear();
  matchedenergy2.clear();
  dx1_seg_.clear();
  dy1_seg_.clear();
  dxSig1_seg_.clear();
  dySig1_seg_.clear();
  ddxdz1_seg_.clear();
  ddydz1_seg_.clear();
  ddxdzSig1_seg_.clear();
  ddydzSig1_seg_.clear();
  dx2_seg_.clear();
  dy2_seg_.clear();
  dxSig2_seg_.clear();
  dySig2_seg_.clear();
  ddxdz2_seg_.clear();
  ddydz2_seg_.clear();
  ddxdzSig2_seg_.clear();
  ddydzSig2_seg_.clear();

  // Clear gen info vectors
  mass_gen.clear();
  pt_gen.clear();
  eta_gen.clear();
  phi_gen.clear();
  status_gen.clear();
  idmom.clear();
  y_gen.clear();
  iddau1.clear();
  iddau2.clear();
  iddau3.clear();

  matchGen_D0pT_.clear();
  matchGen_D0eta_.clear();
  matchGen_D0phi_.clear();
  matchGen_D0mass_.clear();
  matchGen_D0y_.clear();
  matchGen_D0charge_.clear();
  matchGen_D0pdgId_.clear();

  matchGen_D0Dau1_pT_.clear();
  matchGen_D0Dau1_eta_.clear();
  matchGen_D0Dau1_phi_.clear();
  matchGen_D0Dau1_mass_.clear();
  matchGen_D0Dau1_y_.clear();
  matchGen_D0Dau1_charge_.clear();
  matchGen_D0Dau1_pdgId_.clear();

  matchGen_D0Dau2_pT_.clear();
  matchGen_D0Dau2_eta_.clear();
  matchGen_D0Dau2_phi_.clear();
  matchGen_D0Dau2_mass_.clear();
  matchGen_D0Dau2_y_.clear();
  matchGen_D0Dau2_charge_.clear();
  matchGen_D0Dau2_pdgId_.clear();

  matchGen_D1pT_.clear();
  matchGen_D1eta_.clear();
  matchGen_D1phi_.clear();
  matchGen_D1mass_.clear();
  matchGen_D1y_.clear();
  matchGen_D1decayLength2D_.clear();
  matchGen_D1decayLength3D_.clear();
  matchGen_D1angle2D_.clear();
  matchGen_D1angle3D_.clear();
  matchGen_D1ancestorId_.clear();
  matchGen_D1ancestorFlavor_.clear();
  matchGen_D1charge_.clear();
  matchGen_D1pdgId_.clear();

  gen_D0pT_.clear();
  gen_D0eta_.clear();
  gen_D0phi_.clear();
  gen_D0mass_.clear();
  gen_D0y_.clear();
  gen_D0charge_.clear();
  gen_D0pdgId_.clear();

  gen_D0Dau1_pT_.clear();
  gen_D0Dau1_eta_.clear();
  gen_D0Dau1_phi_.clear();
  gen_D0Dau1_mass_.clear();
  gen_D0Dau1_y_.clear();
  gen_D0Dau1_charge_.clear();
  gen_D0Dau1_pdgId_.clear();

  gen_D0Dau2_pT_.clear();
  gen_D0Dau2_eta_.clear();
  gen_D0Dau2_phi_.clear();
  gen_D0Dau2_mass_.clear();
  gen_D0Dau2_y_.clear();
  gen_D0Dau2_charge_.clear();
  gen_D0Dau2_pdgId_.clear();

  gen_D1pT_.clear();
  gen_D1eta_.clear();
  gen_D1phi_.clear();
  gen_D1mass_.clear();
  gen_D1y_.clear();
  gen_D1charge_.clear();
  gen_D1pdgId_.clear();
};
