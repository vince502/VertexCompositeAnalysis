// -*- C++ -*-
//
// Package:    VertexCompositeAnalyzer
// Class:      ChiCNtupleProducer
//
// Description: Simple ntuple producer for chi_c -> J/psi + gamma candidates
//              Reads userFloats/userInts from OniaPhotonProducer / PoniaOniaPhotonProducer output
//              Includes MC generator-level chi_c information
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TTree.h>
#include <TVector3.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class ChiCNtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCNtupleProducer(const edm::ParameterSet&);
  ~ChiCNtupleProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void clearVectors();

  void fillGenInfo(const edm::Event&);
  void fillConversions(const edm::Event&);
  void fillPATConversions(const edm::Event&);
  void fillRecoJpsi(const edm::Event&);

  // Tokens
  edm::EDGetTokenT<pat::CompositeCandidateCollection> candToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<int> centBinToken_;
  edm::EDGetTokenT<reco::Centrality> centSrcToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsiToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> patConversionToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> refitToken_;

  bool isCentrality_;
  bool isMC_;
  bool storeConversions_;
  bool storePATConversions_;
  bool storeRecoJpsi_;
  bool storeRefits_;
  std::string treeName_;
  std::vector<int> genChiPDGs_;
  std::vector<int> genOniaPDGs_;

  // Tree
  TTree* tree_;

  // Event info
  unsigned int runNb_;
  unsigned long long eventNb_;
  unsigned int lsNb_;
  short centrality_;
  int Ntrkoffline_;

  // Candidate vectors
  std::vector<float> chi_mass_;
  std::vector<float> chi_pt_;
  std::vector<float> chi_eta_;
  std::vector<float> chi_phi_;
  std::vector<float> chi_y_;

  // Matching indices
  std::vector<int> oniaIdx_;
  std::vector<int> oniaPdgId_;
  std::vector<float> nominalMass_;
  std::vector<int> convIdx_;

  // Matching quantities
  std::vector<float> dz_;
  std::vector<float> deltaM_;
  std::vector<float> correctedMass_;

  // Onia info
  std::vector<float> onia_mass_;
  std::vector<float> onia_pt_;
  std::vector<float> onia_eta_;
  std::vector<float> onia_phi_;
  std::vector<float> onia_y_;
  std::vector<float> onia_vProb_;

  // Conversion photon info
  std::vector<float> conv_pt_;
  std::vector<float> conv_eta_;
  std::vector<float> conv_phi_;
  std::vector<float> conv_mass_;

  // Conversion vertex
  std::vector<float> convVtx_x_;
  std::vector<float> convVtx_y_;
  std::vector<float> convVtx_z_;
  std::vector<float> convVtx_rho_;
  std::vector<int> conv_flags_;

  // Conversion electron 1 (higher pT)
  std::vector<float> ele1_pt_;
  std::vector<float> ele1_eta_;
  std::vector<float> ele1_phi_;
  std::vector<int> ele1_charge_;
  std::vector<float> ele1_chi2_;
  std::vector<int> ele1_nHits_;
  std::vector<float> ele1_dxy_;
  std::vector<float> ele1_dz_;

  // Conversion electron 2 (lower pT)
  std::vector<float> ele2_pt_;
  std::vector<float> ele2_eta_;
  std::vector<float> ele2_phi_;
  std::vector<int> ele2_charge_;
  std::vector<float> ele2_chi2_;
  std::vector<int> ele2_nHits_;
  std::vector<float> ele2_dxy_;
  std::vector<float> ele2_dz_;

  // Refit info (if available)
  std::vector<float> refit_mass_;
  std::vector<float> refit_vProb_;
  std::vector<float> refit_ctauPV_;
  std::vector<float> refit_ctauErrPV_;
  std::vector<float> refit_cosAlpha_;

  // ============ Generator-level info (MC only) ============
  int Gen_chi_size_;
  std::vector<float> Gen_chi_pt_;
  std::vector<float> Gen_chi_eta_;
  std::vector<float> Gen_chi_phi_;
  std::vector<float> Gen_chi_y_;
  std::vector<float> Gen_chi_mass_;
  std::vector<int> Gen_chi_pdgId_;
  std::vector<int> Gen_chi_status_;
  std::vector<float> Gen_chi_ctau_;
  std::vector<float> Gen_chi_ctau3D_;
  std::vector<float> Gen_chi_vx_;
  std::vector<float> Gen_chi_vy_;
  std::vector<float> Gen_chi_vz_;
  std::vector<int> Gen_chi_jpsiIdx_;
  std::vector<int> Gen_chi_gammaIdx_;
  std::vector<int> Gen_chi_whichRec_;

  int Gen_jpsi_size_;
  std::vector<float> Gen_jpsi_pt_;
  std::vector<float> Gen_jpsi_eta_;
  std::vector<float> Gen_jpsi_phi_;
  std::vector<float> Gen_jpsi_y_;
  std::vector<float> Gen_jpsi_mass_;
  std::vector<int> Gen_jpsi_momIdx_;
  std::vector<int> Gen_jpsi_pdgId_;

  int Gen_mu_size_;
  std::vector<float> Gen_mu_pt_;
  std::vector<float> Gen_mu_eta_;
  std::vector<float> Gen_mu_phi_;
  std::vector<int> Gen_mu_charge_;
  std::vector<int> Gen_mu_jpsiIdx_;

  int Gen_gamma_size_;
  std::vector<float> Gen_gamma_pt_;
  std::vector<float> Gen_gamma_eta_;
  std::vector<float> Gen_gamma_phi_;
  std::vector<int> Gen_gamma_momIdx_;

  // ============ Reco Conversions ============
  int Reco_conv_size_;
  std::vector<float> Reco_conv_pt_;
  std::vector<float> Reco_conv_eta_;
  std::vector<float> Reco_conv_phi_;
  std::vector<float> Reco_conv_mass_;
  std::vector<float> Reco_conv_vtxX_;
  std::vector<float> Reco_conv_vtxY_;
  std::vector<float> Reco_conv_vtxZ_;
  std::vector<float> Reco_conv_vtxRho_;
  std::vector<float> Reco_conv_vtxChi2_;
  std::vector<int> Reco_conv_nTracks_;
  std::vector<int> Reco_conv_quality_;
  std::vector<float> Reco_conv_pairCotThetaSep_;
  std::vector<float> Reco_conv_pairInvMass_;
  std::vector<float> Reco_conv_distOfMinApproach_;
  std::vector<float> Reco_conv_dPhiTracksAtVtx_;
  std::vector<float> Reco_conv_trk1_pt_;
  std::vector<float> Reco_conv_trk1_eta_;
  std::vector<float> Reco_conv_trk1_phi_;
  std::vector<int> Reco_conv_trk1_charge_;
  std::vector<float> Reco_conv_trk1_dxy_;
  std::vector<float> Reco_conv_trk1_dz_;
  std::vector<int> Reco_conv_trk1_nHits_;
  std::vector<float> Reco_conv_trk2_pt_;
  std::vector<float> Reco_conv_trk2_eta_;
  std::vector<float> Reco_conv_trk2_phi_;
  std::vector<int> Reco_conv_trk2_charge_;
  std::vector<float> Reco_conv_trk2_dxy_;
  std::vector<float> Reco_conv_trk2_dz_;
  std::vector<int> Reco_conv_trk2_nHits_;

  // ============ Reco J/psi ============
  int Reco_jpsi_size_;
  std::vector<float> Reco_jpsi_mass_;
  std::vector<float> Reco_jpsi_pt_;
  std::vector<float> Reco_jpsi_eta_;
  std::vector<float> Reco_jpsi_phi_;
  std::vector<float> Reco_jpsi_y_;
  std::vector<float> Reco_jpsi_vProb_;
  std::vector<float> Reco_jpsi_ctauPV_;
  std::vector<float> Reco_jpsi_ctauErrPV_;
  std::vector<float> Reco_mupl_pt_;
  std::vector<float> Reco_mupl_eta_;
  std::vector<float> Reco_mupl_phi_;
  std::vector<float> Reco_mumi_pt_;
  std::vector<float> Reco_mumi_eta_;
  std::vector<float> Reco_mumi_phi_;
};

ChiCNtupleProducer::ChiCNtupleProducer(const edm::ParameterSet& iConfig)
    : candToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("candidateSrc"))),
      pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      centBinToken_(),
      centSrcToken_(),
      genParticleToken_(),
      conversionToken_(),
      jpsiToken_(),
      patConversionToken_(),
      refitToken_(),
      isCentrality_(iConfig.getUntrackedParameter<bool>("isCentrality", false)),
      isMC_(iConfig.getUntrackedParameter<bool>("isMC", false)),
      storeConversions_(iConfig.getUntrackedParameter<bool>("storeConversions", false)),
      storePATConversions_(iConfig.getUntrackedParameter<bool>("storePATConversions", false)),
      storeRecoJpsi_(iConfig.getUntrackedParameter<bool>("storeRecoJpsi", false)),
      storeRefits_(false),
      treeName_(iConfig.getUntrackedParameter<std::string>("treeName", "chiTree")),
      genChiPDGs_(iConfig.getParameter<std::vector<int>>("genChiPDGs")),
      genOniaPDGs_(iConfig.getParameter<std::vector<int>>("genOniaPDGs")),
      tree_(nullptr),
      runNb_(0),
      eventNb_(0),
      lsNb_(0),
      centrality_(-1),
      Ntrkoffline_(-1),
      Gen_chi_size_(0),
      Gen_jpsi_size_(0),
      Gen_mu_size_(0),
      Gen_gamma_size_(0),
      Reco_conv_size_(0),
      Reco_jpsi_size_(0) {
  usesResource("TFileService");

  const auto refitSrc = iConfig.getParameter<edm::InputTag>("refitSrc");
  storeRefits_ = !refitSrc.label().empty();
  if (storeRefits_) {
    refitToken_ = consumes<pat::CompositeCandidateCollection>(refitSrc);
  }
  if (storePATConversions_) {
    patConversionToken_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("patConversions"));
  }
  if (storePATConversions_ && storeConversions_) {
    throw cms::Exception("Configuration") << "Choose PAT or reco conversions explicitly, not both";
  }

  if (isCentrality_) {
    centBinToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    centSrcToken_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  if (isMC_) {
    genParticleToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  }

  if (storeConversions_) {
    conversionToken_ = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
  }

  if (storeRecoJpsi_) {
    jpsiToken_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("jpsiSrc"));
  }
}

void ChiCNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("candidateSrc", edm::InputTag("chiCandidatesMuMu"));
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.addUntracked<bool>("isCentrality", false);
  desc.add<edm::InputTag>("centralityBinLabel", edm::InputTag("centralityBin", "HFtowers"));
  desc.add<edm::InputTag>("centralitySrc", edm::InputTag("hiCentrality"));
  desc.addUntracked<bool>("isMC", false);
  desc.add<edm::InputTag>("genParticles", edm::InputTag("prunedGenParticles"));
  desc.addUntracked<bool>("storeConversions", false);
  desc.addUntracked<bool>("storePATConversions", false);
  desc.add<edm::InputTag>("patConversions", edm::InputTag("oniaPhotonCandidates", "conversions"));
  desc.add<edm::InputTag>("refitSrc", edm::InputTag(""));
  desc.add<edm::InputTag>("conversions", edm::InputTag("reducedEgamma", "reducedConversions"));
  desc.addUntracked<bool>("storeRecoJpsi", false);
  desc.add<edm::InputTag>("jpsiSrc", edm::InputTag("onia2MuMuPAT"));
  desc.addUntracked<std::string>("treeName", "chiTree");
  desc.add<std::vector<int>>("genChiPDGs", {10441, 20443, 445, 10551, 20553, 555});
  desc.add<std::vector<int>>("genOniaPDGs", {443, 553});
  descriptions.add("chiCNtupleProducer", desc);
}

void ChiCNtupleProducer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>(treeName_.c_str(), "Chi_c candidates");

  tree_->Branch("runNb", &runNb_);
  tree_->Branch("eventNb", &eventNb_);
  tree_->Branch("lsNb", &lsNb_);
  tree_->Branch("centrality", &centrality_);
  tree_->Branch("Ntrkoffline", &Ntrkoffline_);

  tree_->Branch("chi_mass", &chi_mass_);
  tree_->Branch("chi_pt", &chi_pt_);
  tree_->Branch("chi_eta", &chi_eta_);
  tree_->Branch("chi_phi", &chi_phi_);
  tree_->Branch("chi_y", &chi_y_);

  tree_->Branch("oniaIdx", &oniaIdx_);
  tree_->Branch("oniaPdgId", &oniaPdgId_);
  tree_->Branch("nominalMass", &nominalMass_);
  tree_->Branch("convIdx", &convIdx_);

  tree_->Branch("dz", &dz_);
  tree_->Branch("deltaM", &deltaM_);
  tree_->Branch("correctedMass", &correctedMass_);

  tree_->Branch("onia_mass", &onia_mass_);
  tree_->Branch("onia_pt", &onia_pt_);
  tree_->Branch("onia_eta", &onia_eta_);
  tree_->Branch("onia_phi", &onia_phi_);
  tree_->Branch("onia_y", &onia_y_);
  tree_->Branch("onia_vProb", &onia_vProb_);

  tree_->Branch("conv_pt", &conv_pt_);
  tree_->Branch("conv_eta", &conv_eta_);
  tree_->Branch("conv_phi", &conv_phi_);
  tree_->Branch("conv_mass", &conv_mass_);

  tree_->Branch("convVtx_x", &convVtx_x_);
  tree_->Branch("convVtx_y", &convVtx_y_);
  tree_->Branch("convVtx_z", &convVtx_z_);
  tree_->Branch("convVtx_rho", &convVtx_rho_);
  tree_->Branch("conv_flags", &conv_flags_);

  tree_->Branch("ele1_pt", &ele1_pt_);
  tree_->Branch("ele1_eta", &ele1_eta_);
  tree_->Branch("ele1_phi", &ele1_phi_);
  tree_->Branch("ele1_charge", &ele1_charge_);
  tree_->Branch("ele1_chi2", &ele1_chi2_);
  tree_->Branch("ele1_nHits", &ele1_nHits_);
  tree_->Branch("ele1_dxy", &ele1_dxy_);
  tree_->Branch("ele1_dz", &ele1_dz_);

  tree_->Branch("ele2_pt", &ele2_pt_);
  tree_->Branch("ele2_eta", &ele2_eta_);
  tree_->Branch("ele2_phi", &ele2_phi_);
  tree_->Branch("ele2_charge", &ele2_charge_);
  tree_->Branch("ele2_chi2", &ele2_chi2_);
  tree_->Branch("ele2_nHits", &ele2_nHits_);
  tree_->Branch("ele2_dxy", &ele2_dxy_);
  tree_->Branch("ele2_dz", &ele2_dz_);

  tree_->Branch("refit_mass", &refit_mass_);
  tree_->Branch("refit_vProb", &refit_vProb_);
  tree_->Branch("refit_ctauPV", &refit_ctauPV_);
  tree_->Branch("refit_ctauErrPV", &refit_ctauErrPV_);
  tree_->Branch("refit_cosAlpha", &refit_cosAlpha_);

  if (isMC_) {
    tree_->Branch("Gen_chi_size", &Gen_chi_size_);
    tree_->Branch("Gen_chi_pt", &Gen_chi_pt_);
    tree_->Branch("Gen_chi_eta", &Gen_chi_eta_);
    tree_->Branch("Gen_chi_phi", &Gen_chi_phi_);
    tree_->Branch("Gen_chi_y", &Gen_chi_y_);
    tree_->Branch("Gen_chi_mass", &Gen_chi_mass_);
    tree_->Branch("Gen_chi_pdgId", &Gen_chi_pdgId_);
    tree_->Branch("Gen_chi_status", &Gen_chi_status_);
    tree_->Branch("Gen_chi_ctau", &Gen_chi_ctau_);
    tree_->Branch("Gen_chi_ctau3D", &Gen_chi_ctau3D_);
    tree_->Branch("Gen_chi_vx", &Gen_chi_vx_);
    tree_->Branch("Gen_chi_vy", &Gen_chi_vy_);
    tree_->Branch("Gen_chi_vz", &Gen_chi_vz_);
    tree_->Branch("Gen_chi_jpsiIdx", &Gen_chi_jpsiIdx_);
    tree_->Branch("Gen_chi_gammaIdx", &Gen_chi_gammaIdx_);
    tree_->Branch("Gen_chi_whichRec", &Gen_chi_whichRec_);

    tree_->Branch("Gen_jpsi_size", &Gen_jpsi_size_);
    tree_->Branch("Gen_jpsi_pt", &Gen_jpsi_pt_);
    tree_->Branch("Gen_jpsi_eta", &Gen_jpsi_eta_);
    tree_->Branch("Gen_jpsi_phi", &Gen_jpsi_phi_);
    tree_->Branch("Gen_jpsi_y", &Gen_jpsi_y_);
    tree_->Branch("Gen_jpsi_mass", &Gen_jpsi_mass_);
    tree_->Branch("Gen_jpsi_momIdx", &Gen_jpsi_momIdx_);
    tree_->Branch("Gen_jpsi_pdgId", &Gen_jpsi_pdgId_);

    tree_->Branch("Gen_mu_size", &Gen_mu_size_);
    tree_->Branch("Gen_mu_pt", &Gen_mu_pt_);
    tree_->Branch("Gen_mu_eta", &Gen_mu_eta_);
    tree_->Branch("Gen_mu_phi", &Gen_mu_phi_);
    tree_->Branch("Gen_mu_charge", &Gen_mu_charge_);
    tree_->Branch("Gen_mu_jpsiIdx", &Gen_mu_jpsiIdx_);

    tree_->Branch("Gen_gamma_size", &Gen_gamma_size_);
    tree_->Branch("Gen_gamma_pt", &Gen_gamma_pt_);
    tree_->Branch("Gen_gamma_eta", &Gen_gamma_eta_);
    tree_->Branch("Gen_gamma_phi", &Gen_gamma_phi_);
    tree_->Branch("Gen_gamma_momIdx", &Gen_gamma_momIdx_);
  }

  if (storeConversions_ || storePATConversions_) {
    tree_->Branch("Reco_conv_size", &Reco_conv_size_);
    tree_->Branch("Reco_conv_pt", &Reco_conv_pt_);
    tree_->Branch("Reco_conv_eta", &Reco_conv_eta_);
    tree_->Branch("Reco_conv_phi", &Reco_conv_phi_);
    tree_->Branch("Reco_conv_mass", &Reco_conv_mass_);
    tree_->Branch("Reco_conv_vtxX", &Reco_conv_vtxX_);
    tree_->Branch("Reco_conv_vtxY", &Reco_conv_vtxY_);
    tree_->Branch("Reco_conv_vtxZ", &Reco_conv_vtxZ_);
    tree_->Branch("Reco_conv_vtxRho", &Reco_conv_vtxRho_);
    tree_->Branch("Reco_conv_vtxChi2", &Reco_conv_vtxChi2_);
    tree_->Branch("Reco_conv_nTracks", &Reco_conv_nTracks_);
    tree_->Branch("Reco_conv_quality", &Reco_conv_quality_);
    tree_->Branch("Reco_conv_pairCotThetaSep", &Reco_conv_pairCotThetaSep_);
    tree_->Branch("Reco_conv_pairInvMass", &Reco_conv_pairInvMass_);
    tree_->Branch("Reco_conv_distOfMinApproach", &Reco_conv_distOfMinApproach_);
    tree_->Branch("Reco_conv_dPhiTracksAtVtx", &Reco_conv_dPhiTracksAtVtx_);
    tree_->Branch("Reco_conv_trk1_pt", &Reco_conv_trk1_pt_);
    tree_->Branch("Reco_conv_trk1_eta", &Reco_conv_trk1_eta_);
    tree_->Branch("Reco_conv_trk1_phi", &Reco_conv_trk1_phi_);
    tree_->Branch("Reco_conv_trk1_charge", &Reco_conv_trk1_charge_);
    tree_->Branch("Reco_conv_trk1_dxy", &Reco_conv_trk1_dxy_);
    tree_->Branch("Reco_conv_trk1_dz", &Reco_conv_trk1_dz_);
    tree_->Branch("Reco_conv_trk1_nHits", &Reco_conv_trk1_nHits_);
    tree_->Branch("Reco_conv_trk2_pt", &Reco_conv_trk2_pt_);
    tree_->Branch("Reco_conv_trk2_eta", &Reco_conv_trk2_eta_);
    tree_->Branch("Reco_conv_trk2_phi", &Reco_conv_trk2_phi_);
    tree_->Branch("Reco_conv_trk2_charge", &Reco_conv_trk2_charge_);
    tree_->Branch("Reco_conv_trk2_dxy", &Reco_conv_trk2_dxy_);
    tree_->Branch("Reco_conv_trk2_dz", &Reco_conv_trk2_dz_);
    tree_->Branch("Reco_conv_trk2_nHits", &Reco_conv_trk2_nHits_);
  }

  if (storeRecoJpsi_) {
    tree_->Branch("Reco_jpsi_size", &Reco_jpsi_size_);
    tree_->Branch("Reco_jpsi_mass", &Reco_jpsi_mass_);
    tree_->Branch("Reco_jpsi_pt", &Reco_jpsi_pt_);
    tree_->Branch("Reco_jpsi_eta", &Reco_jpsi_eta_);
    tree_->Branch("Reco_jpsi_phi", &Reco_jpsi_phi_);
    tree_->Branch("Reco_jpsi_y", &Reco_jpsi_y_);
    tree_->Branch("Reco_jpsi_vProb", &Reco_jpsi_vProb_);
    tree_->Branch("Reco_jpsi_ctauPV", &Reco_jpsi_ctauPV_);
    tree_->Branch("Reco_jpsi_ctauErrPV", &Reco_jpsi_ctauErrPV_);
    tree_->Branch("Reco_mupl_pt", &Reco_mupl_pt_);
    tree_->Branch("Reco_mupl_eta", &Reco_mupl_eta_);
    tree_->Branch("Reco_mupl_phi", &Reco_mupl_phi_);
    tree_->Branch("Reco_mumi_pt", &Reco_mumi_pt_);
    tree_->Branch("Reco_mumi_eta", &Reco_mumi_eta_);
    tree_->Branch("Reco_mumi_phi", &Reco_mumi_phi_);
  }
}

void ChiCNtupleProducer::clearVectors() {
  chi_mass_.clear();
  chi_pt_.clear();
  chi_eta_.clear();
  chi_phi_.clear();
  chi_y_.clear();
  oniaIdx_.clear();
  oniaPdgId_.clear();
  nominalMass_.clear();
  convIdx_.clear();
  dz_.clear();
  deltaM_.clear();
  correctedMass_.clear();

  onia_mass_.clear();
  onia_pt_.clear();
  onia_eta_.clear();
  onia_phi_.clear();
  onia_y_.clear();
  onia_vProb_.clear();

  conv_pt_.clear();
  conv_eta_.clear();
  conv_phi_.clear();
  conv_mass_.clear();
  convVtx_x_.clear();
  convVtx_y_.clear();
  convVtx_z_.clear();
  convVtx_rho_.clear();
  conv_flags_.clear();

  ele1_pt_.clear();
  ele1_eta_.clear();
  ele1_phi_.clear();
  ele1_charge_.clear();
  ele1_chi2_.clear();
  ele1_nHits_.clear();
  ele1_dxy_.clear();
  ele1_dz_.clear();

  ele2_pt_.clear();
  ele2_eta_.clear();
  ele2_phi_.clear();
  ele2_charge_.clear();
  ele2_chi2_.clear();
  ele2_nHits_.clear();
  ele2_dxy_.clear();
  ele2_dz_.clear();

  refit_mass_.clear();
  refit_vProb_.clear();
  refit_ctauPV_.clear();
  refit_ctauErrPV_.clear();
  refit_cosAlpha_.clear();

  Gen_chi_size_ = 0;
  Gen_chi_pt_.clear();
  Gen_chi_eta_.clear();
  Gen_chi_phi_.clear();
  Gen_chi_y_.clear();
  Gen_chi_mass_.clear();
  Gen_chi_pdgId_.clear();
  Gen_chi_status_.clear();
  Gen_chi_ctau_.clear();
  Gen_chi_ctau3D_.clear();
  Gen_chi_vx_.clear();
  Gen_chi_vy_.clear();
  Gen_chi_vz_.clear();
  Gen_chi_jpsiIdx_.clear();
  Gen_chi_gammaIdx_.clear();
  Gen_chi_whichRec_.clear();

  Gen_jpsi_size_ = 0;
  Gen_jpsi_pt_.clear();
  Gen_jpsi_eta_.clear();
  Gen_jpsi_phi_.clear();
  Gen_jpsi_y_.clear();
  Gen_jpsi_mass_.clear();
  Gen_jpsi_momIdx_.clear();
  Gen_jpsi_pdgId_.clear();

  Gen_mu_size_ = 0;
  Gen_mu_pt_.clear();
  Gen_mu_eta_.clear();
  Gen_mu_phi_.clear();
  Gen_mu_charge_.clear();
  Gen_mu_jpsiIdx_.clear();

  Gen_gamma_size_ = 0;
  Gen_gamma_pt_.clear();
  Gen_gamma_eta_.clear();
  Gen_gamma_phi_.clear();
  Gen_gamma_momIdx_.clear();

  Reco_conv_size_ = 0;
  Reco_conv_pt_.clear();
  Reco_conv_eta_.clear();
  Reco_conv_phi_.clear();
  Reco_conv_mass_.clear();
  Reco_conv_vtxX_.clear();
  Reco_conv_vtxY_.clear();
  Reco_conv_vtxZ_.clear();
  Reco_conv_vtxRho_.clear();
  Reco_conv_vtxChi2_.clear();
  Reco_conv_nTracks_.clear();
  Reco_conv_quality_.clear();
  Reco_conv_pairCotThetaSep_.clear();
  Reco_conv_pairInvMass_.clear();
  Reco_conv_distOfMinApproach_.clear();
  Reco_conv_dPhiTracksAtVtx_.clear();
  Reco_conv_trk1_pt_.clear();
  Reco_conv_trk1_eta_.clear();
  Reco_conv_trk1_phi_.clear();
  Reco_conv_trk1_charge_.clear();
  Reco_conv_trk1_dxy_.clear();
  Reco_conv_trk1_dz_.clear();
  Reco_conv_trk1_nHits_.clear();
  Reco_conv_trk2_pt_.clear();
  Reco_conv_trk2_eta_.clear();
  Reco_conv_trk2_phi_.clear();
  Reco_conv_trk2_charge_.clear();
  Reco_conv_trk2_dxy_.clear();
  Reco_conv_trk2_dz_.clear();
  Reco_conv_trk2_nHits_.clear();

  Reco_jpsi_size_ = 0;
  Reco_jpsi_mass_.clear();
  Reco_jpsi_pt_.clear();
  Reco_jpsi_eta_.clear();
  Reco_jpsi_phi_.clear();
  Reco_jpsi_y_.clear();
  Reco_jpsi_vProb_.clear();
  Reco_jpsi_ctauPV_.clear();
  Reco_jpsi_ctauErrPV_.clear();
  Reco_mupl_pt_.clear();
  Reco_mupl_eta_.clear();
  Reco_mupl_phi_.clear();
  Reco_mumi_pt_.clear();
  Reco_mumi_eta_.clear();
  Reco_mumi_phi_.clear();
}

void ChiCNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  clearVectors();

  runNb_ = iEvent.id().run();
  eventNb_ = iEvent.id().event();
  lsNb_ = iEvent.luminosityBlock();

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  centrality_ = -1;
  Ntrkoffline_ = -1;
  if (isCentrality_) {
    edm::Handle<int> centBin;
    iEvent.getByToken(centBinToken_, centBin);
    if (centBin.isValid())
      centrality_ = *centBin;

    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(centSrcToken_, cent);
    if (cent.isValid())
      Ntrkoffline_ = cent->Ntracks();
  }

  edm::Handle<pat::CompositeCandidateCollection> candidates;
  iEvent.getByToken(candToken_, candidates);

  std::map<int, const pat::CompositeCandidate*> refits;
  if (storeRefits_) {
    edm::Handle<pat::CompositeCandidateCollection> refitHandle;
    iEvent.getByToken(refitToken_, refitHandle);
    if (refitHandle.isValid()) {
      for (const auto& fit : *refitHandle) {
        const int idx = fit.hasUserInt("Index") ? fit.userInt("Index")
                      : (fit.hasUserInt("rawIndex") ? fit.userInt("rawIndex") : -1);
        if (idx >= 0)
          refits[idx] = &fit;
      }
    }
  }

  size_t candidateIndex = 0;
  if (candidates.isValid()) {
    for (const auto& inputCand : *candidates) {
      pat::CompositeCandidate cand(inputCand);
      const auto* onia = dynamic_cast<const pat::CompositeCandidate*>(cand.daughter("dimuon"));
      const auto* photon = dynamic_cast<const pat::CompositeCandidate*>(cand.daughter("photon"));
      if (onia && photon) {
        const float nomMass = cand.hasUserFloat("nominalMass") ? cand.userFloat("nominalMass") : 3.0969f;
        if (!cand.hasUserFloat("oniaMass"))
          cand.addUserFloat("oniaMass", onia->mass());
        if (!cand.hasUserFloat("oniaPt"))
          cand.addUserFloat("oniaPt", onia->pt());
        if (!cand.hasUserFloat("oniaEta"))
          cand.addUserFloat("oniaEta", onia->eta());
        if (!cand.hasUserFloat("oniaPhi"))
          cand.addUserFloat("oniaPhi", onia->phi());
        if (!cand.hasUserFloat("oniaY"))
          cand.addUserFloat("oniaY", onia->rapidity());
        if (!cand.hasUserFloat("oniaVProb") && onia->hasUserFloat("vProb"))
          cand.addUserFloat("oniaVProb", onia->userFloat("vProb"));
        if (!cand.hasUserFloat("deltaM"))
          cand.addUserFloat("deltaM", cand.mass() - onia->mass());
        if (!cand.hasUserFloat("correctedMass"))
          cand.addUserFloat("correctedMass", cand.mass() - onia->mass() + nomMass);
        if (!cand.hasUserFloat("convPt"))
          cand.addUserFloat("convPt", photon->pt());
        if (!cand.hasUserFloat("convEta"))
          cand.addUserFloat("convEta", photon->eta());
        if (!cand.hasUserFloat("convPhi"))
          cand.addUserFloat("convPhi", photon->phi());
        if (!cand.hasUserFloat("convMass"))
          cand.addUserFloat("convMass", photon->mass());
        if (!cand.hasUserFloat("convVtxX"))
          cand.addUserFloat("convVtxX", photon->vx());
        if (!cand.hasUserFloat("convVtxY"))
          cand.addUserFloat("convVtxY", photon->vy());
        if (!cand.hasUserFloat("convVtxZ"))
          cand.addUserFloat("convVtxZ", photon->vz());
        if (!cand.hasUserFloat("convVtxRho"))
          cand.addUserFloat("convVtxRho", photon->vertex().rho());
        if (!cand.hasUserInt("convFlags"))
          cand.addUserInt("convFlags", photon->hasUserInt("flags") ? photon->userInt("flags") : -1);

        const auto* track0 = photon->hasUserData("track0") ? photon->userData<reco::Track>("track0") : nullptr;
        const auto* track1 = photon->hasUserData("track1") ? photon->userData<reco::Track>("track1") : nullptr;
        if (track0 && track1 && !cand.hasUserFloat("ele1_pt")) {
          if (track0->pt() < track1->pt())
            std::swap(track0, track1);
          reco::Candidate::Point oniaVtx = onia->vertex();
          if (onia->hasUserData("commonVertex")) {
            const reco::Vertex* vtxPtr = onia->userData<reco::Vertex>("commonVertex");
            if (vtxPtr != nullptr)
              oniaVtx = vtxPtr->position();
          }
          unsigned leg = 0;
          for (const auto* track : {track0, track1}) {
            const std::string prefix = "ele" + std::to_string(++leg) + "_";
            cand.addUserFloat(prefix + "pt", track->pt());
            cand.addUserFloat(prefix + "eta", track->eta());
            cand.addUserFloat(prefix + "phi", track->phi());
            cand.addUserInt(prefix + "charge", track->charge());
            cand.addUserFloat(prefix + "chi2", track->normalizedChi2());
            cand.addUserInt(prefix + "nHits", track->numberOfValidHits());
            cand.addUserFloat(prefix + "dxy", track->dxy(oniaVtx));
            cand.addUserFloat(prefix + "dz", track->dz(oniaVtx));
          }
        }
      }

      chi_mass_.push_back(cand.mass());
      chi_pt_.push_back(cand.pt());
      chi_eta_.push_back(cand.eta());
      chi_phi_.push_back(cand.phi());
      chi_y_.push_back(cand.rapidity());

      const int resolvedOniaIdx = cand.hasUserInt("oniaIdx") ? cand.userInt("oniaIdx")
                                : (cand.hasUserInt("oniaIndex") ? cand.userInt("oniaIndex") : -1);
      const int resolvedConvIdx = cand.hasUserInt("convIdx") ? cand.userInt("convIdx")
                                : (cand.hasUserInt("photonSourceIndex") ? cand.userInt("photonSourceIndex") : -1);
      oniaIdx_.push_back(resolvedOniaIdx);
      convIdx_.push_back(resolvedConvIdx);
      oniaPdgId_.push_back(cand.hasUserInt("oniaPdgId") ? cand.userInt("oniaPdgId") : 0);
      nominalMass_.push_back(cand.hasUserFloat("nominalMass") ? cand.userFloat("nominalMass") : -999.f);

      dz_.push_back(cand.hasUserFloat("dz") ? cand.userFloat("dz") : -999.f);
      deltaM_.push_back(cand.hasUserFloat("deltaM") ? cand.userFloat("deltaM") : -999.f);
      correctedMass_.push_back(cand.hasUserFloat("correctedMass") ? cand.userFloat("correctedMass") : -999.f);

      onia_mass_.push_back(cand.hasUserFloat("oniaMass") ? cand.userFloat("oniaMass") : -999.f);
      onia_pt_.push_back(cand.hasUserFloat("oniaPt") ? cand.userFloat("oniaPt") : -999.f);
      onia_eta_.push_back(cand.hasUserFloat("oniaEta") ? cand.userFloat("oniaEta") : -999.f);
      onia_phi_.push_back(cand.hasUserFloat("oniaPhi") ? cand.userFloat("oniaPhi") : -999.f);
      onia_y_.push_back(cand.hasUserFloat("oniaY") ? cand.userFloat("oniaY") : -999.f);
      onia_vProb_.push_back(cand.hasUserFloat("oniaVProb") ? cand.userFloat("oniaVProb") : -999.f);

      conv_pt_.push_back(cand.hasUserFloat("convPt") ? cand.userFloat("convPt") : -999.f);
      conv_eta_.push_back(cand.hasUserFloat("convEta") ? cand.userFloat("convEta") : -999.f);
      conv_phi_.push_back(cand.hasUserFloat("convPhi") ? cand.userFloat("convPhi") : -999.f);
      conv_mass_.push_back(cand.hasUserFloat("convMass") ? cand.userFloat("convMass") : -999.f);

      convVtx_x_.push_back(cand.hasUserFloat("convVtxX") ? cand.userFloat("convVtxX") : -999.f);
      convVtx_y_.push_back(cand.hasUserFloat("convVtxY") ? cand.userFloat("convVtxY") : -999.f);
      convVtx_z_.push_back(cand.hasUserFloat("convVtxZ") ? cand.userFloat("convVtxZ") : -999.f);
      convVtx_rho_.push_back(cand.hasUserFloat("convVtxRho") ? cand.userFloat("convVtxRho") : -999.f);
      conv_flags_.push_back(cand.hasUserInt("convFlags") ? cand.userInt("convFlags") : -1);

      ele1_pt_.push_back(cand.hasUserFloat("ele1_pt") ? cand.userFloat("ele1_pt") : -999.f);
      ele1_eta_.push_back(cand.hasUserFloat("ele1_eta") ? cand.userFloat("ele1_eta") : -999.f);
      ele1_phi_.push_back(cand.hasUserFloat("ele1_phi") ? cand.userFloat("ele1_phi") : -999.f);
      ele1_charge_.push_back(cand.hasUserInt("ele1_charge") ? cand.userInt("ele1_charge") : -999);
      ele1_chi2_.push_back(cand.hasUserFloat("ele1_chi2") ? cand.userFloat("ele1_chi2") : -999.f);
      ele1_nHits_.push_back(cand.hasUserInt("ele1_nHits") ? cand.userInt("ele1_nHits") : -999);
      ele1_dxy_.push_back(cand.hasUserFloat("ele1_dxy") ? cand.userFloat("ele1_dxy") : -999.f);
      ele1_dz_.push_back(cand.hasUserFloat("ele1_dz") ? cand.userFloat("ele1_dz") : -999.f);

      ele2_pt_.push_back(cand.hasUserFloat("ele2_pt") ? cand.userFloat("ele2_pt") : -999.f);
      ele2_eta_.push_back(cand.hasUserFloat("ele2_eta") ? cand.userFloat("ele2_eta") : -999.f);
      ele2_phi_.push_back(cand.hasUserFloat("ele2_phi") ? cand.userFloat("ele2_phi") : -999.f);
      ele2_charge_.push_back(cand.hasUserInt("ele2_charge") ? cand.userInt("ele2_charge") : -999);
      ele2_chi2_.push_back(cand.hasUserFloat("ele2_chi2") ? cand.userFloat("ele2_chi2") : -999.f);
      ele2_nHits_.push_back(cand.hasUserInt("ele2_nHits") ? cand.userInt("ele2_nHits") : -999);
      ele2_dxy_.push_back(cand.hasUserFloat("ele2_dxy") ? cand.userFloat("ele2_dxy") : -999.f);
      ele2_dz_.push_back(cand.hasUserFloat("ele2_dz") ? cand.userFloat("ele2_dz") : -999.f);

      const auto found = refits.find(static_cast<int>(candidateIndex++));
      const auto* fit = (found == refits.end()) ? nullptr : found->second;
      refit_mass_.push_back(fit ? static_cast<float>(fit->mass()) : -999.f);
      refit_vProb_.push_back((fit && fit->hasUserFloat("vProb")) ? fit->userFloat("vProb") : -999.f);
      refit_ctauPV_.push_back((fit && fit->hasUserFloat("ctauPV")) ? fit->userFloat("ctauPV") : -999.f);
      refit_ctauErrPV_.push_back((fit && fit->hasUserFloat("ctauErrPV")) ? fit->userFloat("ctauErrPV") : -999.f);
      refit_cosAlpha_.push_back((fit && fit->hasUserFloat("cosAlpha")) ? fit->userFloat("cosAlpha") : -999.f);
    }
  }

  if (isMC_) {
    fillGenInfo(iEvent);
  }

  if (storeConversions_) {
    fillConversions(iEvent);
  }

  if (storePATConversions_) {
    fillPATConversions(iEvent);
  }

  if (storeRecoJpsi_) {
    fillRecoJpsi(iEvent);
  }

  tree_->Fill();
}

void ChiCNtupleProducer::fillGenInfo(const edm::Event& iEvent) {
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  if (!genParticles.isValid()) {
    return;
  }

  const auto& chiPDGs = genChiPDGs_;

  for (size_t i = 0; i < genParticles->size(); ++i) {
    const auto& gen = (*genParticles)[i];
    const int absPdgId = std::abs(gen.pdgId());

    const bool isChi = std::find(chiPDGs.begin(), chiPDGs.end(), absPdgId) != chiPDGs.end();
    if (!isChi)
      continue;

    if (gen.status() != 2)
      continue;
    if (gen.numberOfDaughters() < 2)
      continue;

    const reco::GenParticle* jpsiDau = nullptr;
    const reco::GenParticle* gammaDau = nullptr;

    for (size_t d = 0; d < gen.numberOfDaughters(); ++d) {
      const reco::Candidate* dau = gen.daughter(d);
      if (!dau)
        continue;

      const int dauPdgId = std::abs(dau->pdgId());
      if (std::find(genOniaPDGs_.begin(), genOniaPDGs_.end(), dauPdgId) != genOniaPDGs_.end()) {
        jpsiDau = dynamic_cast<const reco::GenParticle*>(dau);
      } else if (dauPdgId == 22) {
        gammaDau = dynamic_cast<const reco::GenParticle*>(dau);
      }
    }

    if (!jpsiDau || !gammaDau)
      continue;

    const int chiIdx = Gen_chi_size_;

    Gen_chi_pt_.push_back(gen.pt());
    Gen_chi_eta_.push_back(gen.eta());
    Gen_chi_phi_.push_back(gen.phi());
    Gen_chi_y_.push_back(gen.rapidity());
    Gen_chi_mass_.push_back(gen.mass());
    Gen_chi_pdgId_.push_back(gen.pdgId());
    Gen_chi_status_.push_back(gen.status());
    Gen_chi_vx_.push_back(gen.vx());
    Gen_chi_vy_.push_back(gen.vy());
    Gen_chi_vz_.push_back(gen.vz());

    float ctau = 0.f;
    float ctau3D = 0.f;
    if (gen.numberOfMothers() > 0) {
      const reco::Candidate* mom = gen.mother();
      if (mom) {
        TVector3 prodVtx(mom->vx(), mom->vy(), mom->vz());
        TVector3 decayVtx(gen.vx(), gen.vy(), gen.vz());
        TVector3 momentum(gen.px(), gen.py(), gen.pz());
        TVector3 diff = decayVtx - prodVtx;

        if (momentum.Perp() > 0.)
          ctau = static_cast<float>(diff.Perp() * gen.mass() / momentum.Perp());
        if (momentum.Mag() > 0.)
          ctau3D = static_cast<float>(diff.Mag() * gen.mass() / momentum.Mag());
      }
    }
    Gen_chi_ctau_.push_back(ctau * 10.0f);
    Gen_chi_ctau3D_.push_back(ctau3D * 10.0f);

    const int jpsiIdx = Gen_jpsi_size_;
    Gen_chi_jpsiIdx_.push_back(jpsiIdx);

    Gen_jpsi_pt_.push_back(jpsiDau->pt());
    Gen_jpsi_eta_.push_back(jpsiDau->eta());
    Gen_jpsi_phi_.push_back(jpsiDau->phi());
    Gen_jpsi_y_.push_back(jpsiDau->rapidity());
    Gen_jpsi_mass_.push_back(jpsiDau->mass());
    Gen_jpsi_momIdx_.push_back(chiIdx);
    Gen_jpsi_pdgId_.push_back(jpsiDau->pdgId());
    Gen_jpsi_size_++;

    const int gammaIdx = Gen_gamma_size_;
    Gen_chi_gammaIdx_.push_back(gammaIdx);

    Gen_gamma_pt_.push_back(gammaDau->pt());
    Gen_gamma_eta_.push_back(gammaDau->eta());
    Gen_gamma_phi_.push_back(gammaDau->phi());
    Gen_gamma_momIdx_.push_back(chiIdx);
    Gen_gamma_size_++;

    for (size_t dd = 0; dd < jpsiDau->numberOfDaughters(); ++dd) {
      const reco::Candidate* mu = jpsiDau->daughter(dd);
      if (!mu)
        continue;
      if (std::abs(mu->pdgId()) != 13)
        continue;
      if (mu->status() != 1)
        continue;

      Gen_mu_pt_.push_back(mu->pt());
      Gen_mu_eta_.push_back(mu->eta());
      Gen_mu_phi_.push_back(mu->phi());
      Gen_mu_charge_.push_back(mu->pdgId() == 13 ? -1 : 1);
      Gen_mu_jpsiIdx_.push_back(jpsiIdx);
      Gen_mu_size_++;
    }

    Gen_chi_whichRec_.push_back(-1);
    Gen_chi_size_++;
  }

  for (int igen = 0; igen < Gen_chi_size_; ++igen) {
    float minDR = 0.15f;
    int bestReco = -1;

    for (size_t irec = 0; irec < chi_mass_.size(); ++irec) {
      if (oniaPdgId_[irec] != 0 &&
          std::abs(oniaPdgId_[irec]) != std::abs(Gen_jpsi_pdgId_[Gen_chi_jpsiIdx_[igen]]))
        continue;
      const float dEta = Gen_chi_eta_[igen] - chi_eta_[irec];
      float dPhi = Gen_chi_phi_[igen] - chi_phi_[irec];
      while (dPhi > M_PI)
        dPhi -= 2 * M_PI;
      while (dPhi < -M_PI)
        dPhi += 2 * M_PI;
      const float dR = std::sqrt(dEta * dEta + dPhi * dPhi);

      const float dMass = std::abs(Gen_chi_mass_[igen] - chi_mass_[irec]);
      if (dR < minDR && dMass < 0.5f) {
        minDR = dR;
        bestReco = static_cast<int>(irec);
      }
    }
    Gen_chi_whichRec_[igen] = bestReco;
  }
}

void ChiCNtupleProducer::fillConversions(const edm::Event& iEvent) {
  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionToken_, conversions);

  if (!conversions.isValid()) {
    return;
  }

  for (const auto& conv : *conversions) {
    const auto& p4 = conv.refittedPair4Momentum();
    Reco_conv_pt_.push_back(p4.pt());
    Reco_conv_eta_.push_back(p4.eta());
    Reco_conv_phi_.push_back(p4.phi());
    Reco_conv_mass_.push_back(p4.mass());

    const auto& vtx = conv.conversionVertex();
    Reco_conv_vtxX_.push_back(vtx.x());
    Reco_conv_vtxY_.push_back(vtx.y());
    Reco_conv_vtxZ_.push_back(vtx.z());
    Reco_conv_vtxRho_.push_back(std::sqrt(vtx.x() * vtx.x() + vtx.y() * vtx.y()));
    Reco_conv_vtxChi2_.push_back(vtx.chi2());

    Reco_conv_nTracks_.push_back(conv.nTracks());

    int quality = 0;
    if (conv.isConverted())
      quality |= 1;
    if (conv.quality(reco::Conversion::generalTracksOnly))
      quality |= 2;
    if (conv.quality(reco::Conversion::arbitratedEcalSeeded))
      quality |= 4;
    if (conv.quality(reco::Conversion::arbitratedMerged))
      quality |= 8;
    if (conv.quality(reco::Conversion::highPurity))
      quality |= 16;
    Reco_conv_quality_.push_back(quality);

    Reco_conv_pairCotThetaSep_.push_back(conv.pairCotThetaSeparation());
    Reco_conv_pairInvMass_.push_back(conv.pairInvariantMass());
    Reco_conv_distOfMinApproach_.push_back(conv.distOfMinimumApproach());
    Reco_conv_dPhiTracksAtVtx_.push_back(conv.dPhiTracksAtVtx());

    const auto& tracks = conv.tracks();
    bool trk1Valid = false;
    bool trk2Valid = false;

    if (tracks.size() >= 1 && tracks[0].isAvailable()) {
      try {
        const auto& trk1 = tracks[0];
        Reco_conv_trk1_pt_.push_back(trk1->pt());
        Reco_conv_trk1_eta_.push_back(trk1->eta());
        Reco_conv_trk1_phi_.push_back(trk1->phi());
        Reco_conv_trk1_charge_.push_back(trk1->charge());
        Reco_conv_trk1_dxy_.push_back(trk1->dxy());
        Reco_conv_trk1_dz_.push_back(trk1->dz());
        Reco_conv_trk1_nHits_.push_back(trk1->numberOfValidHits());
        trk1Valid = true;
      } catch (...) {
      }
    }
    if (!trk1Valid) {
      const auto& pin = conv.tracksPin();
      if (pin.size() >= 1) {
        Reco_conv_trk1_pt_.push_back(pin[0].Rho());
        Reco_conv_trk1_eta_.push_back(pin[0].Eta());
        Reco_conv_trk1_phi_.push_back(pin[0].Phi());
      } else {
        Reco_conv_trk1_pt_.push_back(-999.f);
        Reco_conv_trk1_eta_.push_back(-999.f);
        Reco_conv_trk1_phi_.push_back(-999.f);
      }
      Reco_conv_trk1_charge_.push_back(0);
      Reco_conv_trk1_dxy_.push_back(-999.f);
      Reco_conv_trk1_dz_.push_back(-999.f);
      Reco_conv_trk1_nHits_.push_back(0);
    }

    if (tracks.size() >= 2 && tracks[1].isAvailable()) {
      try {
        const auto& trk2 = tracks[1];
        Reco_conv_trk2_pt_.push_back(trk2->pt());
        Reco_conv_trk2_eta_.push_back(trk2->eta());
        Reco_conv_trk2_phi_.push_back(trk2->phi());
        Reco_conv_trk2_charge_.push_back(trk2->charge());
        Reco_conv_trk2_dxy_.push_back(trk2->dxy());
        Reco_conv_trk2_dz_.push_back(trk2->dz());
        Reco_conv_trk2_nHits_.push_back(trk2->numberOfValidHits());
        trk2Valid = true;
      } catch (...) {
      }
    }
    if (!trk2Valid) {
      const auto& pin = conv.tracksPin();
      if (pin.size() >= 2) {
        Reco_conv_trk2_pt_.push_back(pin[1].Rho());
        Reco_conv_trk2_eta_.push_back(pin[1].Eta());
        Reco_conv_trk2_phi_.push_back(pin[1].Phi());
      } else {
        Reco_conv_trk2_pt_.push_back(-999.f);
        Reco_conv_trk2_eta_.push_back(-999.f);
        Reco_conv_trk2_phi_.push_back(-999.f);
      }
      Reco_conv_trk2_charge_.push_back(0);
      Reco_conv_trk2_dxy_.push_back(-999.f);
      Reco_conv_trk2_dz_.push_back(-999.f);
      Reco_conv_trk2_nHits_.push_back(0);
    }

    ++Reco_conv_size_;
  }
}

void ChiCNtupleProducer::fillPATConversions(const edm::Event& event) {
  const auto patConversions = event.getHandle(patConversionToken_);
  if (!patConversions.isValid()) {
    return;
  }
  for (const auto& photon : *patConversions) {
    Reco_conv_pt_.push_back(photon.pt());
    Reco_conv_eta_.push_back(photon.eta());
    Reco_conv_phi_.push_back(photon.phi());
    Reco_conv_mass_.push_back(photon.mass());
    Reco_conv_vtxX_.push_back(photon.vx());
    Reco_conv_vtxY_.push_back(photon.vy());
    Reco_conv_vtxZ_.push_back(photon.vz());
    Reco_conv_vtxRho_.push_back(photon.vertex().rho());
    Reco_conv_vtxChi2_.push_back(-999.f);
    const auto* track0 = photon.hasUserData("track0") ? photon.userData<reco::Track>("track0") : nullptr;
    const auto* track1 = photon.hasUserData("track1") ? photon.userData<reco::Track>("track1") : nullptr;
    Reco_conv_nTracks_.push_back((track0 != nullptr) + (track1 != nullptr));
    Reco_conv_quality_.push_back(photon.hasUserInt("flags") ? photon.userInt("flags") : -1);
    Reco_conv_pairCotThetaSep_.push_back(-999.f);
    Reco_conv_pairInvMass_.push_back(photon.mass());
    Reco_conv_distOfMinApproach_.push_back(-999.f);
    Reco_conv_dPhiTracksAtVtx_.push_back(-999.f);
    Reco_conv_trk1_pt_.push_back(track0 ? static_cast<float>(track0->pt()) : -999.f);
    Reco_conv_trk1_eta_.push_back(track0 ? static_cast<float>(track0->eta()) : -999.f);
    Reco_conv_trk1_phi_.push_back(track0 ? static_cast<float>(track0->phi()) : -999.f);
    Reco_conv_trk1_charge_.push_back(track0 ? track0->charge() : 0);
    Reco_conv_trk1_dxy_.push_back(track0 ? static_cast<float>(track0->dxy()) : -999.f);
    Reco_conv_trk1_dz_.push_back(track0 ? static_cast<float>(track0->dz()) : -999.f);
    Reco_conv_trk1_nHits_.push_back(track0 ? track0->numberOfValidHits() : 0);
    Reco_conv_trk2_pt_.push_back(track1 ? static_cast<float>(track1->pt()) : -999.f);
    Reco_conv_trk2_eta_.push_back(track1 ? static_cast<float>(track1->eta()) : -999.f);
    Reco_conv_trk2_phi_.push_back(track1 ? static_cast<float>(track1->phi()) : -999.f);
    Reco_conv_trk2_charge_.push_back(track1 ? track1->charge() : 0);
    Reco_conv_trk2_dxy_.push_back(track1 ? static_cast<float>(track1->dxy()) : -999.f);
    Reco_conv_trk2_dz_.push_back(track1 ? static_cast<float>(track1->dz()) : -999.f);
    Reco_conv_trk2_nHits_.push_back(track1 ? track1->numberOfValidHits() : 0);
    ++Reco_conv_size_;
  }
}

void ChiCNtupleProducer::fillRecoJpsi(const edm::Event& iEvent) {
  edm::Handle<pat::CompositeCandidateCollection> jpsis;
  iEvent.getByToken(jpsiToken_, jpsis);
  if (!jpsis.isValid()) {
    return;
  }

  for (const auto& jpsi : *jpsis) {
    Reco_jpsi_mass_.push_back(jpsi.mass());
    Reco_jpsi_pt_.push_back(jpsi.pt());
    Reco_jpsi_eta_.push_back(jpsi.eta());
    Reco_jpsi_phi_.push_back(jpsi.phi());
    Reco_jpsi_y_.push_back(jpsi.rapidity());
    Reco_jpsi_vProb_.push_back(jpsi.hasUserFloat("vProb") ? jpsi.userFloat("vProb") : -999.f);
    Reco_jpsi_ctauPV_.push_back(jpsi.hasUserFloat("ppdlPV") ? jpsi.userFloat("ppdlPV")
                                : (jpsi.hasUserFloat("ctauPV") ? jpsi.userFloat("ctauPV") : -999.f));
    Reco_jpsi_ctauErrPV_.push_back(jpsi.hasUserFloat("ppdlErrPV") ? jpsi.userFloat("ppdlErrPV")
                                   : (jpsi.hasUserFloat("ctauErrPV") ? jpsi.userFloat("ctauErrPV") : -999.f));

    const reco::Candidate* dau1 = jpsi.daughter("muon1");
    const reco::Candidate* dau2 = jpsi.daughter("muon2");
    if ((!dau1 || !dau2) && jpsi.numberOfDaughters() >= 2) {
      dau1 = jpsi.daughter(0);
      dau2 = jpsi.daughter(1);
    }
    const reco::Candidate* mupl =
        (dau1 && dau1->charge() > 0) ? dau1 : ((dau2 && dau2->charge() > 0) ? dau2 : nullptr);
    const reco::Candidate* mumi =
        (dau1 && dau1->charge() < 0) ? dau1 : ((dau2 && dau2->charge() < 0) ? dau2 : nullptr);

    Reco_mupl_pt_.push_back(mupl ? static_cast<float>(mupl->pt()) : -999.f);
    Reco_mupl_eta_.push_back(mupl ? static_cast<float>(mupl->eta()) : -999.f);
    Reco_mupl_phi_.push_back(mupl ? static_cast<float>(mupl->phi()) : -999.f);
    Reco_mumi_pt_.push_back(mumi ? static_cast<float>(mumi->pt()) : -999.f);
    Reco_mumi_eta_.push_back(mumi ? static_cast<float>(mumi->eta()) : -999.f);
    Reco_mumi_phi_.push_back(mumi ? static_cast<float>(mumi->phi()) : -999.f);
    ++Reco_jpsi_size_;
  }
}

DEFINE_FWK_MODULE(ChiCNtupleProducer);
