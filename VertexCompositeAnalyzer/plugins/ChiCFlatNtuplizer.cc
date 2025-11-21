#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCFlatNtuplizer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TTree.h"

#include <Eigen/Dense>

namespace {
struct EventShapeResult {
  double sphericity{0.0};
  std::array<double, 3> eigenvalues{{0.0, 0.0, 0.0}};
};

EventShapeResult computeEventShape(const std::array<const reco::Candidate*, 4>& daughters) {
  EventShapeResult result;
  Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
  double sumP2 = 0.0;

  for (const auto* dau : daughters) {
    if (!dau)
      continue;
    Eigen::Vector3d p(dau->px(), dau->py(), dau->pz());
    tensor += p * p.transpose();
    sumP2 += p.squaredNorm();
  }

  if (sumP2 <= 0.0)
    return result;

  tensor /= sumP2;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(tensor);
  if (solver.info() != Eigen::Success)
    return result;

  result.eigenvalues = {{solver.eigenvalues()(0), solver.eigenvalues()(1), solver.eigenvalues()(2)}};
  std::sort(result.eigenvalues.begin(), result.eigenvalues.end(), std::greater<double>());
  result.sphericity = 1.5 * (result.eigenvalues[1] + result.eigenvalues[2]);
  return result;
}

double computeAcoplanarity(const reco::Candidate* lead, const reco::Candidate* sublead) {
  if (!lead || !sublead)
    return 0.0;
  const double delta = reco::deltaPhi(lead->phi(), sublead->phi());
  return 1.0 - std::abs(delta) / M_PI;
}

// Compute q2 vector magnitude for 4 pions: q2 = (1/N) * sum(exp(2i*phi))
// Returns magnitude of the q2 vector
double computeQ2Magnitude(const std::array<const reco::Candidate*, 4>& daughters) {
  double q2x = 0.0;
  double q2y = 0.0;
  int nValid = 0;
  
  for (const auto* dau : daughters) {
    if (!dau)
      continue;
    const double phi = dau->phi();
    q2x += std::cos(2.0 * phi);
    q2y += std::sin(2.0 * phi);
    ++nValid;
  }
  
  if (nValid == 0)
    return 0.0;
  
  q2x /= nValid;
  q2y /= nValid;
  
  return std::sqrt(q2x * q2x + q2y * q2y);
}

// Compute q2 phase
double computeQ2Phase(const std::array<const reco::Candidate*, 4>& daughters) {
  double q2x = 0.0;
  double q2y = 0.0;
  int nValid = 0;
  
  for (const auto* dau : daughters) {
    if (!dau)
      continue;
    const double phi = dau->phi();
    q2x += std::cos(2.0 * phi);
    q2y += std::sin(2.0 * phi);
    ++nValid;
  }
  
  if (nValid == 0)
    return 0.0;
  
  q2x /= nValid;
  q2y /= nValid;
  
  return std::atan2(q2y, q2x);
}
}  // namespace

ChiCFlatNtuplizer::ChiCFlatNtuplizer(const edm::ParameterSet& cfg)
    : useGenMatching_(false), tree_(nullptr) {
  usesResource("TFileService");

  treeName_ = cfg.getUntrackedParameter<std::string>("treeName", "ChiCFlatNtuple");
  const auto pvTag = cfg.getParameter<edm::InputTag>("primaryVertices");
  pvToken_ = consumes<VertexCollection>(pvTag);

  // Optional: generator particle matching
  if (cfg.existsAs<edm::InputTag>("genParticles")) {
    const auto genTag = cfg.getParameter<edm::InputTag>("genParticles");
    genToken_ = consumes<GenParticleCollection>(genTag);
    useGenMatching_ = true;
  }

  const auto& sourcePsets = cfg.getParameter<std::vector<edm::ParameterSet>>("sources");
  if (sourcePsets.empty()) {
    throw cms::Exception("Configuration") << "ChiCFlatNtuplizer requires at least one source.";
  }

  sources_.reserve(sourcePsets.size());
  for (size_t idx = 0; idx < sourcePsets.size(); ++idx) {
    const auto& ps = sourcePsets[idx];
    SourceConfig sc;
    sc.tag = ps.getParameter<edm::InputTag>("collection");
    sc.name = ps.getParameter<std::string>("name");
    sc.pdgId = ps.getParameter<int>("pdgId");
    sc.index = static_cast<unsigned int>(idx);
    sc.token = consumes<pat::CompositeCandidateCollection>(sc.tag);
    sources_.push_back(sc);
  }
}

void ChiCFlatNtuplizer::beginJob() {
  if (!fileService_.isAvailable()) {
    throw cms::Exception("Configuration") << "TFileService is required for ChiCFlatNtuplizer.";
  }

  tree_ = fileService_->make<TTree>(treeName_.c_str(), treeName_.c_str());

  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");
  tree_->Branch("source_index", &sourceIndex_, "source_index/I");
  tree_->Branch("source_pdgId", &sourcePdgId_, "source_pdgId/I");
  tree_->Branch("source_label", &sourceLabel_);

  tree_->Branch("chic_mass", &candMass_, "chic_mass/F");
  tree_->Branch("chic_pt", &candPt_, "chic_pt/F");
  tree_->Branch("chic_eta", &candEta_, "chic_eta/F");
  tree_->Branch("chic_phi", &candPhi_, "chic_phi/F");
  tree_->Branch("chic_y", &candRapidity_, "chic_y/F");
  tree_->Branch("chic_charge", &candCharge_, "chic_charge/I");
  tree_->Branch("chic_vx", &candVx_, "chic_vx/F");
  tree_->Branch("chic_vy", &candVy_, "chic_vy/F");
  tree_->Branch("chic_vz", &candVz_, "chic_vz/F");
  tree_->Branch("chic_d3d", &candD3D_, "chic_d3d/F");

  tree_->Branch("acoplanarity", &candAcoplanarity_, "acoplanarity/F");
  tree_->Branch("sphericity", &candSphericity_, "sphericity/F");
  tree_->Branch("pca_lambda1", &candLambda1_, "pca_lambda1/F");
  tree_->Branch("pca_lambda2", &candLambda2_, "pca_lambda2/F");
  tree_->Branch("pca_lambda3", &candLambda3_, "pca_lambda3/F");

  tree_->Branch("pi1_pt", &dauPt_[0], "pi1_pt/F");
  tree_->Branch("pi2_pt", &dauPt_[1], "pi2_pt/F");
  tree_->Branch("pi3_pt", &dauPt_[2], "pi3_pt/F");
  tree_->Branch("pi4_pt", &dauPt_[3], "pi4_pt/F");

  tree_->Branch("pi1_eta", &dauEta_[0], "pi1_eta/F");
  tree_->Branch("pi2_eta", &dauEta_[1], "pi2_eta/F");
  tree_->Branch("pi3_eta", &dauEta_[2], "pi3_eta/F");
  tree_->Branch("pi4_eta", &dauEta_[3], "pi4_eta/F");

  tree_->Branch("pi1_phi", &dauPhi_[0], "pi1_phi/F");
  tree_->Branch("pi2_phi", &dauPhi_[1], "pi2_phi/F");
  tree_->Branch("pi3_phi", &dauPhi_[2], "pi3_phi/F");
  tree_->Branch("pi4_phi", &dauPhi_[3], "pi4_phi/F");

  tree_->Branch("pi1_charge", &dauCharge_[0], "pi1_charge/I");
  tree_->Branch("pi2_charge", &dauCharge_[1], "pi2_charge/I");
  tree_->Branch("pi3_charge", &dauCharge_[2], "pi3_charge/I");
  tree_->Branch("pi4_charge", &dauCharge_[3], "pi4_charge/I");

  tree_->Branch("pi1_dxy", &dauDxy_[0], "pi1_dxy/F");
  tree_->Branch("pi2_dxy", &dauDxy_[1], "pi2_dxy/F");
  tree_->Branch("pi3_dxy", &dauDxy_[2], "pi3_dxy/F");
  tree_->Branch("pi4_dxy", &dauDxy_[3], "pi4_dxy/F");

  tree_->Branch("pi1_dz", &dauDz_[0], "pi1_dz/F");
  tree_->Branch("pi2_dz", &dauDz_[1], "pi2_dz/F");
  tree_->Branch("pi3_dz", &dauDz_[2], "pi3_dz/F");
  tree_->Branch("pi4_dz", &dauDz_[3], "pi4_dz/F");

  tree_->Branch("pi1_d3d", &dauD3d_[0], "pi1_d3d/F");
  tree_->Branch("pi2_d3d", &dauD3d_[1], "pi2_d3d/F");
  tree_->Branch("pi3_d3d", &dauD3d_[2], "pi3_d3d/F");
  tree_->Branch("pi4_d3d", &dauD3d_[3], "pi4_d3d/F");

  tree_->Branch("pi1_mass", &dauMass_[0], "pi1_mass/F");
  tree_->Branch("pi2_mass", &dauMass_[1], "pi2_mass/F");
  tree_->Branch("pi3_mass", &dauMass_[2], "pi3_mass/F");
  tree_->Branch("pi4_mass", &dauMass_[3], "pi4_mass/F");

  tree_->Branch("pi1_dedx", &dauDeDx_[0], "pi1_dedx/F");
  tree_->Branch("pi2_dedx", &dauDeDx_[1], "pi2_dedx/F");
  tree_->Branch("pi3_dedx", &dauDeDx_[2], "pi3_dedx/F");
  tree_->Branch("pi4_dedx", &dauDeDx_[3], "pi4_dedx/F");

  tree_->Branch("dca_12", &pairDca_[0], "dca_12/F");
  tree_->Branch("dca_13", &pairDca_[1], "dca_13/F");
  tree_->Branch("dca_14", &pairDca_[2], "dca_14/F");
  tree_->Branch("dca_23", &pairDca_[3], "dca_23/F");
  tree_->Branch("dca_24", &pairDca_[4], "dca_24/F");
  tree_->Branch("dca_34", &pairDca_[5], "dca_34/F");

  // Generator-level and q-vector variables
  tree_->Branch("genMatch", &genMatch_, "genMatch/I");
  tree_->Branch("genMass", &genMass_, "genMass/F");
  tree_->Branch("genPt", &genPt_, "genPt/F");
  tree_->Branch("genEta", &genEta_, "genEta/F");
  tree_->Branch("genPhi", &genPhi_, "genPhi/F");
  tree_->Branch("genY", &genY_, "genY/F");
  tree_->Branch("q2Magnitude", &q2Magnitude_, "q2Magnitude/F");
  tree_->Branch("q2Phase", &q2Phase_, "q2Phase/F");
}

void ChiCFlatNtuplizer::resetBranches() {
  sourceIndex_ = -1;
  sourcePdgId_ = 0;
  sourceLabel_.clear();

  candCharge_ = 0;
  candMass_ = -1.f;
  candPt_ = -1.f;
  candEta_ = -99.f;
  candPhi_ = -99.f;
  candRapidity_ = -99.f;
  candVx_ = candVy_ = candVz_ = 0.f;
  candD3D_ = 0.f;

  candAcoplanarity_ = 0.f;
  candSphericity_ = 0.f;
  candLambda1_ = candLambda2_ = candLambda3_ = 0.f;

  dauPt_.fill(-1.f);
  dauEta_.fill(-99.f);
  dauPhi_.fill(-99.f);
  dauCharge_.fill(0);
  dauDxy_.fill(0.f);
  dauDz_.fill(0.f);
  dauD3d_.fill(0.f);
  dauMass_.fill(-1.f);
  dauDeDx_.fill(-1.f);

  pairDca_.fill(0.f);

  genMatch_ = 0;
  genMass_ = -1.f;
  genPt_ = -1.f;
  genEta_ = -99.f;
  genPhi_ = -99.f;
  genY_ = -99.f;
  q2Magnitude_ = 0.f;
  q2Phase_ = 0.f;
}

void ChiCFlatNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  const reco::Vertex* primaryVertex = nullptr;
  edm::Handle<VertexCollection> pvHandle;
  if (event.getByToken(pvToken_, pvHandle) && pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = &pvHandle->front();
  }

  for (const auto& src : sources_) {
    edm::Handle<pat::CompositeCandidateCollection> handle;
    event.getByToken(src.token, handle);
    if (!handle.isValid())
      continue;

    for (const auto& cand : *handle) {
      resetBranches();
      fillCandidate(src, cand, primaryVertex);
      tree_->Fill();
    }
  }
}

void ChiCFlatNtuplizer::fillCandidate(const SourceConfig& src,
                                       const pat::CompositeCandidate& cand,
                                       const reco::Vertex* primaryVertex,
                                       const edm::Handle<GenParticleCollection>* genParticles) {
  sourceIndex_ = static_cast<int>(src.index);
  sourcePdgId_ = src.pdgId;
  sourceLabel_ = src.name;

  candCharge_ = cand.charge();
  candMass_ = static_cast<float>(cand.mass());
  candPt_ = static_cast<float>(cand.pt());
  candEta_ = static_cast<float>(cand.eta());
  candPhi_ = static_cast<float>(cand.phi());
  candRapidity_ = static_cast<float>(cand.rapidity());

  const auto& vtx = cand.vertex();
  candVx_ = static_cast<float>(vtx.x());
  candVy_ = static_cast<float>(vtx.y());
  candVz_ = static_cast<float>(vtx.z());

  if (primaryVertex) {
    const double dx = vtx.x() - primaryVertex->x();
    const double dy = vtx.y() - primaryVertex->y();
    const double dz = vtx.z() - primaryVertex->z();
    candD3D_ = static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
  } else {
    candD3D_ = -1.f;
  }

  struct DaughterInfo {
    const reco::Candidate* cand{nullptr};
    const reco::Track* track{nullptr};
    double pt{0.0};
  };

  std::vector<DaughterInfo> daughters;
  daughters.reserve(4);
  const unsigned int totalDau = cand.numberOfDaughters();
  for (unsigned int i = 0; i < totalDau; ++i) {
    const auto* dau = cand.daughter(i);
    if (!dau)
      continue;

    const reco::Track* trackPtr = nullptr;
    if (const auto* recoDau = dynamic_cast<const reco::RecoChargedCandidate*>(dau)) {
      const auto trackRef = recoDau->track();
      if (!trackRef.isNull())
        trackPtr = trackRef.get();
    }
    if (!trackPtr)
      trackPtr = dau->bestTrack();

    daughters.push_back({dau, trackPtr, dau->pt()});
  }

  std::sort(daughters.begin(), daughters.end(), [](const DaughterInfo& lhs, const DaughterInfo& rhs) {
    return lhs.pt > rhs.pt;
  });

  if (daughters.size() > 4)
    daughters.resize(4);

  std::array<const reco::Candidate*, 4> orderedDaughters{{nullptr, nullptr, nullptr, nullptr}};
  for (std::size_t i = 0; i < daughters.size(); ++i) {
    const auto& info = daughters[i];
    orderedDaughters[i] = info.cand;

    dauPt_[i] = static_cast<float>(info.cand->pt());
    dauEta_[i] = static_cast<float>(info.cand->eta());
    dauPhi_[i] = static_cast<float>(info.cand->phi());
    dauCharge_[i] = info.cand->charge();
    dauMass_[i] = static_cast<float>(info.cand->mass());

    if (info.track && primaryVertex) {
      const auto& pvPos = primaryVertex->position();
      const double dxy = info.track->dxy(pvPos);
      const double dz = info.track->dz(pvPos);
      dauDxy_[i] = static_cast<float>(dxy);
      dauDz_[i] = static_cast<float>(dz);
      dauD3d_[i] = static_cast<float>(std::sqrt(dxy * dxy + dz * dz));
      
      // dE/dx not available (IsolatedTrack collection removed)
      dauDeDx_[i] = -1.f;
    } else {
      dauDxy_[i] = 0.f;
      dauDz_[i] = 0.f;
      dauD3d_[i] = 0.f;
      dauDeDx_[i] = -1.f;
    }
  }

  const auto computePairDistance = [](const reco::Candidate* first, const reco::Candidate* second) {
    if (!first || !second)
      return 0.f;
    const auto& v1 = first->vertex();
    const auto& v2 = second->vertex();
    const double dx = v1.x() - v2.x();
    const double dy = v1.y() - v2.y();
    const double dz = v1.z() - v2.z();
    return static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
  };

  pairDca_[0] = computePairDistance(orderedDaughters[0], orderedDaughters[1]);
  pairDca_[1] = computePairDistance(orderedDaughters[0], orderedDaughters[2]);
  pairDca_[2] = computePairDistance(orderedDaughters[0], orderedDaughters[3]);
  pairDca_[3] = computePairDistance(orderedDaughters[1], orderedDaughters[2]);
  pairDca_[4] = computePairDistance(orderedDaughters[1], orderedDaughters[3]);
  pairDca_[5] = computePairDistance(orderedDaughters[2], orderedDaughters[3]);

  if (cand.hasUserFloat("acoplanarity")) {
    candAcoplanarity_ = cand.userFloat("acoplanarity");
  } else {
    candAcoplanarity_ = static_cast<float>(computeAcoplanarity(orderedDaughters[0], orderedDaughters[1]));
  }

  if (cand.hasUserFloat("sphericity") && cand.hasUserFloat("pca_lambda1")) {
    candSphericity_ = cand.userFloat("sphericity");
    candLambda1_ = cand.userFloat("pca_lambda1");
    candLambda2_ = cand.userFloat("pca_lambda2");
    candLambda3_ = cand.userFloat("pca_lambda3");
  } else {
    const auto es = computeEventShape(orderedDaughters);
    candSphericity_ = static_cast<float>(es.sphericity);
    candLambda1_ = static_cast<float>(es.eigenvalues[0]);
    candLambda2_ = static_cast<float>(es.eigenvalues[1]);
    candLambda3_ = static_cast<float>(es.eigenvalues[2]);
  }

  // Compute q2 vector for 4 pions
  q2Magnitude_ = static_cast<float>(computeQ2Magnitude(orderedDaughters));
  q2Phase_ = static_cast<float>(computeQ2Phase(orderedDaughters));

  // Generator matching
  if (genParticles && genParticles->isValid()) {
    const reco::GenParticle* genMatch = findGenMatch(cand, *genParticles);
    if (genMatch) {
      genMatch_ = 1;
      genMass_ = static_cast<float>(genMatch->mass());
      genPt_ = static_cast<float>(genMatch->pt());
      genEta_ = static_cast<float>(genMatch->eta());
      genPhi_ = static_cast<float>(genMatch->phi());
      genY_ = static_cast<float>(genMatch->rapidity());
    }
  }
}

const reco::GenParticle* ChiCFlatNtuplizer::findGenMatch(
    const pat::CompositeCandidate& cand,
    const edm::Handle<GenParticleCollection>& genParticles) const {
  const double recoPt = cand.pt();
  const double recoEta = cand.eta();
  const double recoPhi = cand.phi();
  const double recoMass = cand.mass();

  const double maxDeltaR = 0.3;
  const double maxDeltaPt = 0.3;
  const double maxDeltaMass = 0.2;

  const reco::GenParticle* bestMatch = nullptr;
  double bestDeltaR = maxDeltaR;

  for (const auto& genPart : *genParticles) {
    // Look for ChiC (PDG ID 445) or similar
    if (std::abs(genPart.pdgId()) != 445)
      continue;

    // Check if it decays to 4 pions
    std::vector<const reco::GenParticle*> pions;
    collectStablePions(genPart, pions);
    
    int nPlus = 0, nMinus = 0;
    bool valid = true;
    for (const auto* pion : pions) {
      if (std::abs(pion->pdgId()) != 211) {
        valid = false;
        break;
      }
      if (pion->pdgId() > 0) ++nPlus;
      else ++nMinus;
    }
    
    if (!valid || pions.size() != 4 || nPlus != 2 || nMinus != 2)
      continue;

    const double deltaR = reco::deltaR(recoEta, recoPhi, genPart.eta(), genPart.phi());
    const double deltaPt = std::abs(recoPt - genPart.pt()) / recoPt;
    const double deltaMass = std::abs(recoMass - genPart.mass());

    if (deltaR < maxDeltaR && deltaPt < maxDeltaPt && deltaMass < maxDeltaMass) {
      if (deltaR < bestDeltaR) {
        bestDeltaR = deltaR;
        bestMatch = &genPart;
      }
    }
  }

  return bestMatch;
}

void ChiCFlatNtuplizer::collectStablePions(
    const reco::GenParticle& particle,
    std::vector<const reco::GenParticle*>& pions) const {
  const auto nDau = particle.numberOfDaughters();
  if (nDau == 0)
    return;

  for (size_t i = 0; i < nDau; ++i) {
    const auto* dauCandidate = particle.daughter(i);
    if (!dauCandidate)
      continue;

    const auto* dau = dynamic_cast<const reco::GenParticle*>(dauCandidate);
    if (!dau)
      continue;

    if (dau->status() == 1) {
      if (std::abs(dau->pdgId()) == 211) {
        pions.push_back(dau);
      }
    } else {
      collectStablePions(*dau, pions);
    }
  }
}


DEFINE_FWK_MODULE(ChiCFlatNtuplizer);
