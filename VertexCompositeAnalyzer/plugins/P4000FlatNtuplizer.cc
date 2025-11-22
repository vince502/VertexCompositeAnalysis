// -*- C++ -*-
//
// Package:    VertexCompositeAnalyzer
// Class:      P4000FlatNtuplizer
//

#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/P4000FlatNtuplizer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

P4000FlatNtuplizer::P4000FlatNtuplizer(const edm::ParameterSet& cfg)
  : p4000Token_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("p4000Collection"))),
    pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices")))
{
  usesResource("TFileService");
  fileService_ = edm::Service<TFileService>();
}

P4000FlatNtuplizer::~P4000FlatNtuplizer() = default;

void P4000FlatNtuplizer::beginJob() {
  tree_ = fileService_->make<TTree>("P4000Tree", "P4000 -> J/psi + phi");

  // Event info
  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");

  // P(4000) candidate
  tree_->Branch("p4000Mass", &p4000Mass_, "p4000Mass/F");
  tree_->Branch("p4000Pt", &p4000Pt_, "p4000Pt/F");
  tree_->Branch("p4000Eta", &p4000Eta_, "p4000Eta/F");
  tree_->Branch("p4000Phi", &p4000Phi_, "p4000Phi/F");
  tree_->Branch("p4000Rapidity", &p4000Rapidity_, "p4000Rapidity/F");
  tree_->Branch("p4000Vx", &p4000Vx_, "p4000Vx/F");
  tree_->Branch("p4000Vy", &p4000Vy_, "p4000Vy/F");
  tree_->Branch("p4000Vz", &p4000Vz_, "p4000Vz/F");
  tree_->Branch("p4000VtxChi2", &p4000VtxChi2_, "p4000VtxChi2/F");
  tree_->Branch("p4000VtxNdof", &p4000VtxNdof_, "p4000VtxNdof/F");
  tree_->Branch("p4000VtxProb", &p4000VtxProb_, "p4000VtxProb/F");

  // J/ψ
  tree_->Branch("jpsiMass", &jpsiMass_, "jpsiMass/F");
  tree_->Branch("jpsiPt", &jpsiPt_, "jpsiPt/F");
  tree_->Branch("jpsiEta", &jpsiEta_, "jpsiEta/F");
  tree_->Branch("jpsiPhi", &jpsiPhi_, "jpsiPhi/F");
  tree_->Branch("jpsiVx", &jpsiVx_, "jpsiVx/F");
  tree_->Branch("jpsiVy", &jpsiVy_, "jpsiVy/F");
  tree_->Branch("jpsiVz", &jpsiVz_, "jpsiVz/F");

  // φ
  tree_->Branch("phiMass", &phiMass_, "phiMass/F");
  tree_->Branch("phiPt", &phiPt_, "phiPt/F");
  tree_->Branch("phiEta", &phiEta_, "phiEta/F");
  tree_->Branch("phiPhi", &phiPhi_, "phiPhi/F");
  tree_->Branch("phiVx", &phiVx_, "phiVx/F");
  tree_->Branch("phiVy", &phiVy_, "phiVy/F");
  tree_->Branch("phiVz", &phiVz_, "phiVz/F");

  // Muons (J/ψ daughters)
  tree_->Branch("mu1Pt", &muPt_[0], "mu1Pt/F");
  tree_->Branch("mu2Pt", &muPt_[1], "mu2Pt/F");
  tree_->Branch("mu1Eta", &muEta_[0], "mu1Eta/F");
  tree_->Branch("mu2Eta", &muEta_[1], "mu2Eta/F");
  tree_->Branch("mu1Phi", &muPhi_[0], "mu1Phi/F");
  tree_->Branch("mu2Phi", &muPhi_[1], "mu2Phi/F");
  tree_->Branch("mu1Charge", &muCharge_[0], "mu1Charge/I");
  tree_->Branch("mu2Charge", &muCharge_[1], "mu2Charge/I");
  tree_->Branch("mu1Dxy", &muDxy_[0], "mu1Dxy/F");
  tree_->Branch("mu2Dxy", &muDxy_[1], "mu2Dxy/F");
  tree_->Branch("mu1Dz", &muDz_[0], "mu1Dz/F");
  tree_->Branch("mu2Dz", &muDz_[1], "mu2Dz/F");

  // Kaons (φ daughters)
  tree_->Branch("kaon1Pt", &kaonPt_[0], "kaon1Pt/F");
  tree_->Branch("kaon2Pt", &kaonPt_[1], "kaon2Pt/F");
  tree_->Branch("kaon1Eta", &kaonEta_[0], "kaon1Eta/F");
  tree_->Branch("kaon2Eta", &kaonEta_[1], "kaon2Eta/F");
  tree_->Branch("kaon1Phi", &kaonPhi_[0], "kaon1Phi/F");
  tree_->Branch("kaon2Phi", &kaonPhi_[1], "kaon2Phi/F");
  tree_->Branch("kaon1Charge", &kaonCharge_[0], "kaon1Charge/I");
  tree_->Branch("kaon2Charge", &kaonCharge_[1], "kaon2Charge/I");
  tree_->Branch("kaon1Dxy", &kaonDxy_[0], "kaon1Dxy/F");
  tree_->Branch("kaon2Dxy", &kaonDxy_[1], "kaon2Dxy/F");
  tree_->Branch("kaon1Dz", &kaonDz_[0], "kaon1Dz/F");
  tree_->Branch("kaon2Dz", &kaonDz_[1], "kaon2Dz/F");

  // Primary vertex
  tree_->Branch("nPV", &nPV_, "nPV/I");
  tree_->Branch("pvX", &pvX_, "pvX/F");
  tree_->Branch("pvY", &pvY_, "pvY/F");
  tree_->Branch("pvZ", &pvZ_, "pvZ/F");
  tree_->Branch("pvNTracks", &pvNTracks_, "pvNTracks/I");
}

void P4000FlatNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  edm::Handle<pat::CompositeCandidateCollection> p4000Handle;
  event.getByToken(p4000Token_, p4000Handle);

  if (!p4000Handle.isValid() || p4000Handle->empty()) return;

  edm::Handle<reco::VertexCollection> pvHandle;
  event.getByToken(pvToken_, pvHandle);

  const reco::Vertex* primaryVertex = nullptr;
  nPV_ = 0;
  pvNTracks_ = 0;
  pvX_ = pvY_ = pvZ_ = 0.f;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    nPV_ = static_cast<int>(pvHandle->size());
    primaryVertex = &pvHandle->front();
    pvX_ = static_cast<float>(primaryVertex->x());
    pvY_ = static_cast<float>(primaryVertex->y());
    pvZ_ = static_cast<float>(primaryVertex->z());
    pvNTracks_ = static_cast<int>(primaryVertex->tracksSize());
  }

  for (const auto& p4000 : *p4000Handle) {
    resetBranches();
    fillCandidate(p4000, primaryVertex);
    tree_->Fill();
  }
}

void P4000FlatNtuplizer::resetBranches() {
  p4000Mass_ = -1.f;
  p4000Pt_ = -1.f;
  p4000Eta_ = -99.f;
  p4000Phi_ = -99.f;
  p4000Rapidity_ = -99.f;
  p4000Vx_ = p4000Vy_ = p4000Vz_ = 0.f;
  p4000VtxChi2_ = -1.f;
  p4000VtxNdof_ = -1.f;
  p4000VtxProb_ = -1.f;

  jpsiMass_ = -1.f;
  jpsiPt_ = -1.f;
  jpsiEta_ = -99.f;
  jpsiPhi_ = -99.f;
  jpsiVx_ = jpsiVy_ = jpsiVz_ = 0.f;

  phiMass_ = -1.f;
  phiPt_ = -1.f;
  phiEta_ = -99.f;
  phiPhi_ = -99.f;
  phiVx_ = phiVy_ = phiVz_ = 0.f;

  muPt_.fill(-1.f);
  muEta_.fill(-99.f);
  muPhi_.fill(-99.f);
  muCharge_.fill(0);
  muDxy_.fill(0.f);
  muDz_.fill(0.f);

  kaonPt_.fill(-1.f);
  kaonEta_.fill(-99.f);
  kaonPhi_.fill(-99.f);
  kaonCharge_.fill(0);
  kaonDxy_.fill(0.f);
  kaonDz_.fill(0.f);
}

void P4000FlatNtuplizer::endJob() {}

void P4000FlatNtuplizer::fillCandidate(const pat::CompositeCandidate& cand, const reco::Vertex* primaryVertex) {
  // P(4000) info
  p4000Mass_ = static_cast<float>(cand.mass());
  p4000Pt_ = static_cast<float>(cand.pt());
  p4000Eta_ = static_cast<float>(cand.eta());
  p4000Phi_ = static_cast<float>(cand.phi());
  p4000Rapidity_ = static_cast<float>(cand.rapidity());
  p4000Vx_ = static_cast<float>(cand.vx());
  p4000Vy_ = static_cast<float>(cand.vy());
  p4000Vz_ = static_cast<float>(cand.vz());

  if (cand.hasUserFloat("VtxChi2")) {
    p4000VtxChi2_ = cand.userFloat("VtxChi2");
  }
  if (cand.hasUserFloat("VtxNdof")) {
    p4000VtxNdof_ = cand.userFloat("VtxNdof");
  }
  if (cand.hasUserFloat("VtxProb")) {
    p4000VtxProb_ = cand.userFloat("VtxProb");
  }

  // J/ψ info
  const reco::Candidate* jpsi = cand.daughter("Jpsi");
  if (jpsi) {
    jpsiMass_ = static_cast<float>(jpsi->mass());
    jpsiPt_ = static_cast<float>(jpsi->pt());
    jpsiEta_ = static_cast<float>(jpsi->eta());
    jpsiPhi_ = static_cast<float>(jpsi->phi());
    jpsiVx_ = static_cast<float>(jpsi->vx());
    jpsiVy_ = static_cast<float>(jpsi->vy());
    jpsiVz_ = static_cast<float>(jpsi->vz());

    // J/ψ daughters (muons)
    if (jpsi->numberOfDaughters() >= 2) {
      for (size_t i = 0; i < 2 && i < jpsi->numberOfDaughters(); ++i) {
        const auto* mu = jpsi->daughter(i);
        muPt_[i] = static_cast<float>(mu->pt());
        muEta_[i] = static_cast<float>(mu->eta());
        muPhi_[i] = static_cast<float>(mu->phi());
        muCharge_[i] = mu->charge();

        const auto* muTrack = dynamic_cast<const reco::RecoChargedCandidate*>(mu);
        if (muTrack && muTrack->track().isNonnull() && primaryVertex) {
          muDxy_[i] = static_cast<float>(muTrack->track()->dxy(primaryVertex->position()));
          muDz_[i] = static_cast<float>(muTrack->track()->dz(primaryVertex->position()));
        }
      }
    }
  }

  // φ info
  const reco::Candidate* phi = cand.daughter("Phi");
  if (phi) {
    phiMass_ = static_cast<float>(phi->mass());
    phiPt_ = static_cast<float>(phi->pt());
    phiEta_ = static_cast<float>(phi->eta());
    phiPhi_ = static_cast<float>(phi->phi());
    phiVx_ = static_cast<float>(phi->vx());
    phiVy_ = static_cast<float>(phi->vy());
    phiVz_ = static_cast<float>(phi->vz());

    // φ daughters (kaons)
    if (phi->numberOfDaughters() >= 2) {
      for (size_t i = 0; i < 2 && i < phi->numberOfDaughters(); ++i) {
        const auto* kaon = phi->daughter(i);
        kaonPt_[i] = static_cast<float>(kaon->pt());
        kaonEta_[i] = static_cast<float>(kaon->eta());
        kaonPhi_[i] = static_cast<float>(kaon->phi());
        kaonCharge_[i] = kaon->charge();

        const auto* kaonTrack = dynamic_cast<const reco::RecoChargedCandidate*>(kaon);
        if (kaonTrack && kaonTrack->track().isNonnull() && primaryVertex) {
          kaonDxy_[i] = static_cast<float>(kaonTrack->track()->dxy(primaryVertex->position()));
          kaonDz_[i] = static_cast<float>(kaonTrack->track()->dz(primaryVertex->position()));
        }
      }
    }
  }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(P4000FlatNtuplizer);
