#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCFlatNtuplizer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

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
}  // namespace

ChiCFlatNtuplizer::ChiCFlatNtuplizer(const edm::ParameterSet& cfg)
    : tree_(nullptr) {
  usesResource("TFileService");

  treeName_ = cfg.getUntrackedParameter<std::string>("treeName", "ChiCFlatNtuple");
  const auto pvTag = cfg.getParameter<edm::InputTag>("primaryVertices");
  pvToken_ = consumes<VertexCollection>(pvTag);

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

  tree_->Branch("dca_12", &pairDca_[0], "dca_12/F");
  tree_->Branch("dca_13", &pairDca_[1], "dca_13/F");
  tree_->Branch("dca_14", &pairDca_[2], "dca_14/F");
  tree_->Branch("dca_23", &pairDca_[3], "dca_23/F");
  tree_->Branch("dca_24", &pairDca_[4], "dca_24/F");
  tree_->Branch("dca_34", &pairDca_[5], "dca_34/F");
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

  pairDca_.fill(0.f);
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
                                       const reco::Vertex* primaryVertex) {
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

    if (info.track && primaryVertex) {
      const auto& pvPos = primaryVertex->position();
      const double dxy = info.track->dxy(pvPos);
      const double dz = info.track->dz(pvPos);
      dauDxy_[i] = static_cast<float>(dxy);
      dauDz_[i] = static_cast<float>(dz);
      dauD3d_[i] = static_cast<float>(std::sqrt(dxy * dxy + dz * dz));
    } else {
      dauDxy_[i] = 0.f;
      dauDz_[i] = 0.f;
      dauD3d_[i] = 0.f;
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
}

DEFINE_FWK_MODULE(ChiCFlatNtuplizer);
