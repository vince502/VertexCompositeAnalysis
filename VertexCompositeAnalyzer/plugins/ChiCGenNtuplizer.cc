#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCGenNtuplizer.h"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {
constexpr float kInvalidFloat = -1.0e3f;
}

ChiCGenNtuplizer::ChiCGenNtuplizer(const edm::ParameterSet& cfg)
    : tree_(nullptr),
      run_(0),
      lumi_(0),
      event_(0),
      chic_found_(0),
      chic_pt_(kInvalidFloat),
      chic_eta_(kInvalidFloat),
      chic_phi_(kInvalidFloat),
      chic_y_(kInvalidFloat),
      chic_mass_(kInvalidFloat),
      chic_d3d_(kInvalidFloat) {
  usesResource("TFileService");

  const auto genTag = cfg.getParameter<edm::InputTag>("genParticles");
  genToken_ = consumes<GenParticleCollection>(genTag);

  if (cfg.existsAs<edm::InputTag>("primaryVertices")) {
    const auto pvTag = cfg.getParameter<edm::InputTag>("primaryVertices");
    pvToken_ = consumes<VertexCollection>(pvTag);
  }

  treeName_ = cfg.getUntrackedParameter<std::string>("treeName", "ChiCGenNtuple");

  resetBranches();
}

void ChiCGenNtuplizer::beginJob() {
  if (!fileService_.isAvailable()) {
    throw cms::Exception("Configuration")
        << "TFileService is not available. Please add it to the process.";
  }

  tree_ = fileService_->make<TTree>(treeName_.c_str(), treeName_.c_str());

  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");

  tree_->Branch("chic_found", &chic_found_, "chic_found/I");
  tree_->Branch("chic_pt", &chic_pt_, "chic_pt/F");
  tree_->Branch("chic_eta", &chic_eta_, "chic_eta/F");
  tree_->Branch("chic_phi", &chic_phi_, "chic_phi/F");
  tree_->Branch("chic_y", &chic_y_, "chic_y/F");
  tree_->Branch("chic_mass", &chic_mass_, "chic_mass/F");
  tree_->Branch("chic_d3d", &chic_d3d_, "chic_d3d/F");

  tree_->Branch("pi1_pt", &pion_pt_[0], "pi1_pt/F");
  tree_->Branch("pi1_eta", &pion_eta_[0], "pi1_eta/F");
  tree_->Branch("pi1_phi", &pion_phi_[0], "pi1_phi/F");
  tree_->Branch("pi1_d3d", &pion_d3d_[0], "pi1_d3d/F");

  tree_->Branch("pi2_pt", &pion_pt_[1], "pi2_pt/F");
  tree_->Branch("pi2_eta", &pion_eta_[1], "pi2_eta/F");
  tree_->Branch("pi2_phi", &pion_phi_[1], "pi2_phi/F");
  tree_->Branch("pi2_d3d", &pion_d3d_[1], "pi2_d3d/F");

  tree_->Branch("pi3_pt", &pion_pt_[2], "pi3_pt/F");
  tree_->Branch("pi3_eta", &pion_eta_[2], "pi3_eta/F");
  tree_->Branch("pi3_phi", &pion_phi_[2], "pi3_phi/F");
  tree_->Branch("pi3_d3d", &pion_d3d_[2], "pi3_d3d/F");

  tree_->Branch("pi4_pt", &pion_pt_[3], "pi4_pt/F");
  tree_->Branch("pi4_eta", &pion_eta_[3], "pi4_eta/F");
  tree_->Branch("pi4_phi", &pion_phi_[3], "pi4_phi/F");
  tree_->Branch("pi4_d3d", &pion_d3d_[3], "pi4_d3d/F");

  tree_->Branch("dca_12", &dca_pairs_[0], "dca_12/F");
  tree_->Branch("dca_13", &dca_pairs_[1], "dca_13/F");
  tree_->Branch("dca_14", &dca_pairs_[2], "dca_14/F");
  tree_->Branch("dca_23", &dca_pairs_[3], "dca_23/F");
  tree_->Branch("dca_24", &dca_pairs_[4], "dca_24/F");
  tree_->Branch("dca_34", &dca_pairs_[5], "dca_34/F");
}

void ChiCGenNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  resetBranches();

  run_ = event.id().run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  edm::Handle<GenParticleCollection> genHandle;
  event.getByToken(genToken_, genHandle);
  if (!genHandle.isValid()) {
    tree_->Fill();
    std::cout << "GenHandle wrong!" << std::endl;
    return;
  }

  XYZPoint pvPoint(0., 0., 0.);
  if (!pvToken_.isUninitialized()) {
    edm::Handle<VertexCollection> pvHandle;
    event.getByToken(pvToken_, pvHandle);
    if (pvHandle.isValid() && !pvHandle->empty()) {
      pvPoint = pvHandle->front().position();
    }
  }

  const reco::GenParticle* bestChi = nullptr;
  std::vector<const reco::GenParticle*> bestPions;

  for (const auto& particle : *genHandle) {
    if (particle.pdgId() != 445)
      continue;
    std::cout << "Found Chic! 445" << std::endl;

    std::vector<const reco::GenParticle*> pionBuffer;
    collectStablePions(particle, pionBuffer);
    if (pionBuffer.size() != 4)
      continue;

    int nPlus = 0;
    int nMinus = 0;
    bool valid = true;
    for (const auto* pion : pionBuffer) {
      const int id = pion->pdgId();
      if (std::abs(id) != 211) {
        valid = false;
        break;
      }
      if (id > 0) {
        ++nPlus;
      } else {
        ++nMinus;
      }
    }
    if (!valid || nPlus != 2 || nMinus != 2)
      continue;

    if (!bestChi || particle.pt() > bestChi->pt()) {
      bestChi = &particle;
      bestPions = pionBuffer;
    }
  }

  if (bestChi) {
    chic_found_ = 1;
    chic_pt_ = static_cast<float>(bestChi->pt());
    chic_eta_ = static_cast<float>(bestChi->eta());
    chic_phi_ = static_cast<float>(bestChi->phi());
    chic_y_ = static_cast<float>(bestChi->rapidity());
    chic_mass_ = static_cast<float>(bestChi->mass());
    chic_d3d_ = distance3D(bestChi->vertex(), pvPoint);

    std::sort(bestPions.begin(), bestPions.end(),
              [](const reco::GenParticle* lhs, const reco::GenParticle* rhs) {
                return lhs->pt() > rhs->pt();
              });

    for (size_t idx = 0; idx < bestPions.size(); ++idx) {
      const auto* pion = bestPions[idx];
      pion_pt_[idx] = static_cast<float>(pion->pt());
      pion_eta_[idx] = static_cast<float>(pion->eta());
      pion_phi_[idx] = static_cast<float>(pion->phi());
      pion_d3d_[idx] = distance3D(pion->vertex(), pvPoint);
    }

    // Fill pairwise DCA values for ordered combinations (1<->2, 1<->3, ...)
    const std::array<std::pair<size_t, size_t>, 6> combos = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    for (size_t combIdx = 0; combIdx < combos.size(); ++combIdx) {
      const auto [i, j] = combos[combIdx];
      dca_pairs_[combIdx] = pairDca(*bestPions[i], *bestPions[j]);
    }
  }

  tree_->Fill();
}

void ChiCGenNtuplizer::resetBranches() {
  chic_found_ = 0;
  chic_pt_ = kInvalidFloat;
  chic_eta_ = kInvalidFloat;
  chic_phi_ = kInvalidFloat;
  chic_y_ = kInvalidFloat;
  chic_mass_ = kInvalidFloat;
  chic_d3d_ = kInvalidFloat;

  pion_pt_.fill(kInvalidFloat);
  pion_eta_.fill(kInvalidFloat);
  pion_phi_.fill(kInvalidFloat);
  pion_d3d_.fill(kInvalidFloat);

  dca_pairs_.fill(kInvalidFloat);
}

void ChiCGenNtuplizer::collectStablePions(const reco::GenParticle& particle,
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

float ChiCGenNtuplizer::distance3D(const XYZPoint& from, const XYZPoint& to) {
  const double dx = from.x() - to.x();
  const double dy = from.y() - to.y();
  const double dz = from.z() - to.z();
  return static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
}

float ChiCGenNtuplizer::pairDca(const reco::GenParticle& first, const reco::GenParticle& second) {
  const XYZPoint v1 = first.vertex();
  const XYZPoint v2 = second.vertex();

  XYZVector u1 = first.momentum();
  XYZVector u2 = second.momentum();

  if (u1.Mag2() > 0.) {
    u1 = u1.Unit();
  }
  if (u2.Mag2() > 0.) {
    u2 = u2.Unit();
  }

  const XYZVector w0(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());

  const double a = u1.Dot(u1);
  const double b = u1.Dot(u2);
  const double c = u2.Dot(u2);
  const double d = u1.Dot(w0);
  const double e = u2.Dot(w0);
  const double denom = a * c - b * b;

  double sc = 0.0;
  double tc = 0.0;

  if (std::abs(denom) > 1.0e-12) {
    sc = (b * e - c * d) / denom;
    tc = (a * e - b * d) / denom;
  } else {
    // Parallel case; project difference onto one of the directions
    tc = (c > 1.0e-12) ? e / c : 0.0;
  }

  const XYZVector dP = w0 + (sc * u1) - (tc * u2);
  return static_cast<float>(std::sqrt(dP.Mag2()));
}

DEFINE_FWK_MODULE(ChiCGenNtuplizer);
