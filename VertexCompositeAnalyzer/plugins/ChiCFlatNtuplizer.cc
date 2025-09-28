#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCFlatNtuplizer.h"

#include <algorithm>
#include <array>
#include <cmath>

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

  tree_->Branch("pi_pt", dauPt_.data(), "pi_pt[4]/F");
  tree_->Branch("pi_eta", dauEta_.data(), "pi_eta[4]/F");
  tree_->Branch("pi_phi", dauPhi_.data(), "pi_phi[4]/F");
  tree_->Branch("pi_charge", dauCharge_.data(), "pi_charge[4]/I");
  tree_->Branch("pi_dxy", dauDxy_.data(), "pi_dxy[4]/F");
  tree_->Branch("pi_dz", dauDz_.data(), "pi_dz[4]/F");
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

  std::array<const reco::Candidate*, 4> daughters{{nullptr, nullptr, nullptr, nullptr}};
  const unsigned int nDau = std::min<unsigned int>(4, cand.numberOfDaughters());
  for (unsigned int i = 0; i < nDau; ++i) {
    const auto* dau = cand.daughter(i);
    daughters[i] = dau;
    dauPt_[i] = static_cast<float>(dau->pt());
    dauEta_[i] = static_cast<float>(dau->eta());
    dauPhi_[i] = static_cast<float>(dau->phi());
    dauCharge_[i] = dau->charge();

    const reco::Track* trackPtr = nullptr;
    if (const auto* recoDau = dynamic_cast<const reco::RecoChargedCandidate*>(dau)) {
      const auto trackRef = recoDau->track();
      if (!trackRef.isNull())
        trackPtr = trackRef.get();
    }
    if (!trackPtr)
      trackPtr = dau->bestTrack();

    if (trackPtr && primaryVertex) {
      dauDxy_[i] = static_cast<float>(trackPtr->dxy(primaryVertex->position()));
      dauDz_[i] = static_cast<float>(trackPtr->dz(primaryVertex->position()));
    } else {
      dauDxy_[i] = 0.f;
      dauDz_[i] = 0.f;
    }
  }

  if (cand.hasUserFloat("acoplanarity")) {
    candAcoplanarity_ = cand.userFloat("acoplanarity");
  } else {
    candAcoplanarity_ = static_cast<float>(computeAcoplanarity(daughters[0], daughters[1]));
  }

  if (cand.hasUserFloat("sphericity") && cand.hasUserFloat("pca_lambda1")) {
    candSphericity_ = cand.userFloat("sphericity");
    candLambda1_ = cand.userFloat("pca_lambda1");
    candLambda2_ = cand.userFloat("pca_lambda2");
    candLambda3_ = cand.userFloat("pca_lambda3");
  } else {
    const auto es = computeEventShape(daughters);
    candSphericity_ = static_cast<float>(es.sphericity);
    candLambda1_ = static_cast<float>(es.eigenvalues[0]);
    candLambda2_ = static_cast<float>(es.eigenvalues[1]);
    candLambda3_ = static_cast<float>(es.eigenvalues[2]);
  }
}

DEFINE_FWK_MODULE(ChiCFlatNtuplizer);
