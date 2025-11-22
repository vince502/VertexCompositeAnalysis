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
  tree_->Branch("chic_decayLength3D", &candDecayLength3D_, "chic_decayLength3D/F");
  tree_->Branch("chic_decayLength2D", &candDecayLength2D_, "chic_decayLength2D/F");
  tree_->Branch("chic_pointingAngle3D", &candPointingAngle3D_, "chic_pointingAngle3D/F");
  tree_->Branch("chic_pointingAngle2D", &candPointingAngle2D_, "chic_pointingAngle2D/F");
  tree_->Branch("chic_cosPointingAngle3D", &candCosPointingAngle3D_, "chic_cosPointingAngle3D/F");
  tree_->Branch("chic_cosPointingAngle2D", &candCosPointingAngle2D_, "chic_cosPointingAngle2D/F");
  tree_->Branch("chic_dxy", &candDxy_, "chic_dxy/F");
  tree_->Branch("chic_dz", &candDz_, "chic_dz/F");
  tree_->Branch("chic_dxySig", &candDxySig_, "chic_dxySig/F");
  tree_->Branch("chic_dzSig", &candDzSig_, "chic_dzSig/F");

  tree_->Branch("acoplanarity", &candAcoplanarity_, "acoplanarity/F");
  tree_->Branch("sphericity", &candSphericity_, "sphericity/F");
  tree_->Branch("pca_lambda1", &candLambda1_, "pca_lambda1/F");
  tree_->Branch("pca_lambda2", &candLambda2_, "pca_lambda2/F");
  tree_->Branch("pca_lambda3", &candLambda3_, "pca_lambda3/F");

  // Vertex fit quality variables
  tree_->Branch("vtxChi2", &vtxChi2_, "vtxChi2/F");
  tree_->Branch("vtxNdof", &vtxNdof_, "vtxNdof/F");
  tree_->Branch("vtxProb", &vtxProb_, "vtxProb/F");
  tree_->Branch("vtxNormalizedChi2", &vtxNormalizedChi2_, "vtxNormalizedChi2/F");
  
  // Decay length and pointing from vertex fit
  tree_->Branch("decayLength2D", &decayLength2D_, "decayLength2D/F");
  tree_->Branch("decayLength3D", &decayLength3D_, "decayLength3D/F");
  tree_->Branch("decayLengthSig2D", &decayLengthSig2D_, "decayLengthSig2D/F");
  tree_->Branch("decayLengthSig3D", &decayLengthSig3D_, "decayLengthSig3D/F");
  tree_->Branch("alpha2D", &alpha2D_, "alpha2D/F");
  tree_->Branch("alpha3D", &alpha3D_, "alpha3D/F");
  tree_->Branch("cosAlpha2D", &cosAlpha2D_, "cosAlpha2D/F");
  tree_->Branch("cosAlpha3D", &cosAlpha3D_, "cosAlpha3D/F");
  
  // Track DCA variables
  tree_->Branch("trackDCA_min", &trackDCA_min_, "trackDCA_min/F");
  tree_->Branch("trackDCA_max", &trackDCA_max_, "trackDCA_max/F");
  tree_->Branch("trackDCA_avg", &trackDCA_avg_, "trackDCA_avg/F");

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

  tree_->Branch("pi1_dxySig", &dauDxySig_[0], "pi1_dxySig/F");
  tree_->Branch("pi2_dxySig", &dauDxySig_[1], "pi2_dxySig/F");
  tree_->Branch("pi3_dxySig", &dauDxySig_[2], "pi3_dxySig/F");
  tree_->Branch("pi4_dxySig", &dauDxySig_[3], "pi4_dxySig/F");

  tree_->Branch("pi1_dzSig", &dauDzSig_[0], "pi1_dzSig/F");
  tree_->Branch("pi2_dzSig", &dauDzSig_[1], "pi2_dzSig/F");
  tree_->Branch("pi3_dzSig", &dauDzSig_[2], "pi3_dzSig/F");
  tree_->Branch("pi4_dzSig", &dauDzSig_[3], "pi4_dzSig/F");

  tree_->Branch("pi1_dxSig", &dauDxSig_[0], "pi1_dxSig/F");
  tree_->Branch("pi2_dxSig", &dauDxSig_[1], "pi2_dxSig/F");
  tree_->Branch("pi3_dxSig", &dauDxSig_[2], "pi3_dxSig/F");
  tree_->Branch("pi4_dxSig", &dauDxSig_[3], "pi4_dxSig/F");

  tree_->Branch("pi1_dySig", &dauDySig_[0], "pi1_dySig/F");
  tree_->Branch("pi2_dySig", &dauDySig_[1], "pi2_dySig/F");
  tree_->Branch("pi3_dySig", &dauDySig_[2], "pi3_dySig/F");
  tree_->Branch("pi4_dySig", &dauDySig_[3], "pi4_dySig/F");

  tree_->Branch("pi1_chi2", &dauChi2_[0], "pi1_chi2/F");
  tree_->Branch("pi2_chi2", &dauChi2_[1], "pi2_chi2/F");
  tree_->Branch("pi3_chi2", &dauChi2_[2], "pi3_chi2/F");
  tree_->Branch("pi4_chi2", &dauChi2_[3], "pi4_chi2/F");

  tree_->Branch("pi1_nhits", &dauNhits_[0], "pi1_nhits/I");
  tree_->Branch("pi2_nhits", &dauNhits_[1], "pi2_nhits/I");
  tree_->Branch("pi3_nhits", &dauNhits_[2], "pi3_nhits/I");
  tree_->Branch("pi4_nhits", &dauNhits_[3], "pi4_nhits/I");

  tree_->Branch("pi1_npixHits", &dauNpixHits_[0], "pi1_npixHits/I");
  tree_->Branch("pi2_npixHits", &dauNpixHits_[1], "pi2_npixHits/I");
  tree_->Branch("pi3_npixHits", &dauNpixHits_[2], "pi3_npixHits/I");
  tree_->Branch("pi4_npixHits", &dauNpixHits_[3], "pi4_npixHits/I");

  tree_->Branch("pi1_ptErr", &dauPtErr_[0], "pi1_ptErr/F");
  tree_->Branch("pi2_ptErr", &dauPtErr_[1], "pi2_ptErr/F");
  tree_->Branch("pi3_ptErr", &dauPtErr_[2], "pi3_ptErr/F");
  tree_->Branch("pi4_ptErr", &dauPtErr_[3], "pi4_ptErr/F");

  tree_->Branch("pi1_etaErr", &dauEtaErr_[0], "pi1_etaErr/F");
  tree_->Branch("pi2_etaErr", &dauEtaErr_[1], "pi2_etaErr/F");
  tree_->Branch("pi3_etaErr", &dauEtaErr_[2], "pi3_etaErr/F");
  tree_->Branch("pi4_etaErr", &dauEtaErr_[3], "pi4_etaErr/F");

  tree_->Branch("pi1_phiErr", &dauPhiErr_[0], "pi1_phiErr/F");
  tree_->Branch("pi2_phiErr", &dauPhiErr_[1], "pi2_phiErr/F");
  tree_->Branch("pi3_phiErr", &dauPhiErr_[2], "pi3_phiErr/F");
  tree_->Branch("pi4_phiErr", &dauPhiErr_[3], "pi4_phiErr/F");

  tree_->Branch("dca_12", &pairDca_[0], "dca_12/F");
  tree_->Branch("dca_13", &pairDca_[1], "dca_13/F");
  tree_->Branch("dca_14", &pairDca_[2], "dca_14/F");
  tree_->Branch("dca_23", &pairDca_[3], "dca_23/F");
  tree_->Branch("dca_24", &pairDca_[4], "dca_24/F");
  tree_->Branch("dca_34", &pairDca_[5], "dca_34/F");

  tree_->Branch("pairMass_12", &pairMass_[0], "pairMass_12/F");
  tree_->Branch("pairMass_13", &pairMass_[1], "pairMass_13/F");
  tree_->Branch("pairMass_14", &pairMass_[2], "pairMass_14/F");
  tree_->Branch("pairMass_23", &pairMass_[3], "pairMass_23/F");
  tree_->Branch("pairMass_24", &pairMass_[4], "pairMass_24/F");
  tree_->Branch("pairMass_34", &pairMass_[5], "pairMass_34/F");

  tree_->Branch("pairPt_12", &pairPt_[0], "pairPt_12/F");
  tree_->Branch("pairPt_13", &pairPt_[1], "pairPt_13/F");
  tree_->Branch("pairPt_14", &pairPt_[2], "pairPt_14/F");
  tree_->Branch("pairPt_23", &pairPt_[3], "pairPt_23/F");
  tree_->Branch("pairPt_24", &pairPt_[4], "pairPt_24/F");
  tree_->Branch("pairPt_34", &pairPt_[5], "pairPt_34/F");

  tree_->Branch("pairQ2_12", &pairQ2_[0], "pairQ2_12/F");
  tree_->Branch("pairQ2_13", &pairQ2_[1], "pairQ2_13/F");
  tree_->Branch("pairQ2_14", &pairQ2_[2], "pairQ2_14/F");
  tree_->Branch("pairQ2_23", &pairQ2_[3], "pairQ2_23/F");
  tree_->Branch("pairQ2_24", &pairQ2_[4], "pairQ2_24/F");
  tree_->Branch("pairQ2_34", &pairQ2_[5], "pairQ2_34/F");

  tree_->Branch("chic_decay", &chicDecay_, "chic_decay/I");
  tree_->Branch("nPV", &nPV_, "nPV/I");
  tree_->Branch("pvX", &pvX_, "pvX/F");
  tree_->Branch("pvY", &pvY_, "pvY/F");
  tree_->Branch("pvZ", &pvZ_, "pvZ/F");
  tree_->Branch("pvNdof", &pvNdof_, "pvNdof/F");
  tree_->Branch("pvChi2", &pvChi2_, "pvChi2/F");
  tree_->Branch("pvNTracks", &pvNTracks_, "pvNTracks/I");
  tree_->Branch("pvXError", &pvXError_, "pvXError/F");
  tree_->Branch("pvYError", &pvYError_, "pvYError/F");
  tree_->Branch("pvZError", &pvZError_, "pvZError/F");

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
  candD3D_ = -1.f;
  candDecayLength3D_ = -1.f;
  candDecayLength2D_ = -1.f;
  candPointingAngle3D_ = -99.f;
  candPointingAngle2D_ = -99.f;
  candCosPointingAngle3D_ = -99.f;
  candCosPointingAngle2D_ = -99.f;
  candDxy_ = 0.f;
  candDz_ = 0.f;
  candDxySig_ = 0.f;
  candDzSig_ = 0.f;

  candAcoplanarity_ = 0.f;
  candSphericity_ = 0.f;
  candLambda1_ = candLambda2_ = candLambda3_ = 0.f;

  // Vertex fit quality variables
  vtxChi2_ = -1.f;
  vtxNdof_ = -1.f;
  vtxProb_ = -1.f;
  vtxNormalizedChi2_ = -1.f;
  
  // Decay length and pointing from vertex fit
  decayLength2D_ = -1.f;
  decayLength3D_ = -1.f;
  decayLengthSig2D_ = -1.f;
  decayLengthSig3D_ = -1.f;
  alpha2D_ = -99.f;
  alpha3D_ = -99.f;
  cosAlpha2D_ = -99.f;
  cosAlpha3D_ = -99.f;
  
  // Track DCA variables
  trackDCA_min_ = -1.f;
  trackDCA_max_ = -1.f;
  trackDCA_avg_ = -1.f;
  
  // PV error variables
  pvXError_ = -1.f;
  pvYError_ = -1.f;
  pvZError_ = -1.f;

  dauPt_.fill(-1.f);
  dauEta_.fill(-99.f);
  dauPhi_.fill(-99.f);
  dauCharge_.fill(0);
  dauDxy_.fill(0.f);
  dauDz_.fill(0.f);
  dauD3d_.fill(0.f);
  dauMass_.fill(-1.f);
  dauDeDx_.fill(-1.f);
  dauDxySig_.fill(0.f);
  dauDzSig_.fill(0.f);
  dauDxSig_.fill(0.f);
  dauDySig_.fill(0.f);
  dauChi2_.fill(-1.f);
  dauNhits_.fill(0);
  dauNpixHits_.fill(0);
  dauPtErr_.fill(-1.f);
  dauEtaErr_.fill(-1.f);
  dauPhiErr_.fill(-1.f);

  pairDca_.fill(0.f);
  pairMass_.fill(-1.f);
  pairPt_.fill(-1.f);
  pairEta_.fill(-99.f);
  pairPhi_.fill(-99.f);
  pairQ2_.fill(0.f);
  
  chicDecay_ = -1;

  nPV_ = 0;
  pvNTracks_ = 0;
  pvX_ = pvY_ = pvZ_ = 0.f;
  pvNdof_ = -1.f;
  pvChi2_ = -1.f;

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
  nPV_ = 0;  // Initialize to 0
  pvNTracks_ = 0;  // Initialize to 0
  
  if (event.getByToken(pvToken_, pvHandle) && pvHandle.isValid()) {
    nPV_ = static_cast<int>(pvHandle->size());
    if (!pvHandle->empty()) {
      primaryVertex = &pvHandle->front();
      pvX_ = static_cast<float>(primaryVertex->x());
      pvY_ = static_cast<float>(primaryVertex->y());
      pvZ_ = static_cast<float>(primaryVertex->z());
      pvNdof_ = static_cast<float>(primaryVertex->ndof());
      pvChi2_ = static_cast<float>(primaryVertex->chi2());
      pvNTracks_ = static_cast<int>(primaryVertex->tracksSize());  // Track multiplicity
    }
  }

  edm::Handle<GenParticleCollection> genHandle;
  const edm::Handle<GenParticleCollection>* genParticles = nullptr;
  if (useGenMatching_ && event.getByToken(genToken_, genHandle) && genHandle.isValid()) {
    genParticles = &genHandle;
  }

  for (const auto& src : sources_) {
    edm::Handle<pat::CompositeCandidateCollection> handle;
    event.getByToken(src.token, handle);
    if (!handle.isValid() || handle->empty()) {
      // Empty collections are expected - skip silently
      continue;
    }

    for (const auto& cand : *handle) {
      // Mass cut: only store candidates in range 2.5 - 4.2 GeV
      const double candMass = cand.mass();
      if (candMass < 2.5 || candMass > 4.2)
        continue;
      
      resetBranches();
      // Restore event-level info that was reset (nPV and PV info are event-level, not per-candidate)
      nPV_ = static_cast<int>(pvHandle.isValid() ? pvHandle->size() : 0);
      if (primaryVertex) {
        pvX_ = static_cast<float>(primaryVertex->x());
        pvY_ = static_cast<float>(primaryVertex->y());
        pvZ_ = static_cast<float>(primaryVertex->z());
        pvNdof_ = static_cast<float>(primaryVertex->ndof());
        pvChi2_ = static_cast<float>(primaryVertex->chi2());
        pvNTracks_ = static_cast<int>(primaryVertex->tracksSize());  // Track multiplicity
      } else {
        pvNTracks_ = 0;
      }
      fillCandidate(src, cand, primaryVertex, genParticles);
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
  
  // Read vertex fit quality variables from userFloats
  if (cand.hasUserFloat("VtxChi2")) {
    vtxChi2_ = cand.userFloat("VtxChi2");
  }
  if (cand.hasUserFloat("VtxNdof")) {
    vtxNdof_ = cand.userFloat("VtxNdof");
  }
  if (cand.hasUserFloat("VtxProb")) {
    vtxProb_ = cand.userFloat("VtxProb");
  }
  if (cand.hasUserFloat("VtxNormalizedChi2")) {
    vtxNormalizedChi2_ = cand.userFloat("VtxNormalizedChi2");
  }
  
  // Read decay length and pointing from vertex fit
  if (cand.hasUserFloat("decaylength2D")) {
    decayLength2D_ = cand.userFloat("decaylength2D");
  }
  if (cand.hasUserFloat("decaylength3D")) {
    decayLength3D_ = cand.userFloat("decaylength3D");
  }
  if (cand.hasUserFloat("decaylengthsignif2D")) {
    decayLengthSig2D_ = cand.userFloat("decaylengthsignif2D");
  }
  if (cand.hasUserFloat("decaylengthsignif3D")) {
    decayLengthSig3D_ = cand.userFloat("decaylengthsignif3D");
  }
  if (cand.hasUserFloat("alpha2D")) {
    alpha2D_ = cand.userFloat("alpha2D");
  }
  if (cand.hasUserFloat("alpha3D")) {
    alpha3D_ = cand.userFloat("alpha3D");
  }
  if (cand.hasUserFloat("cosAlpha2D")) {
    cosAlpha2D_ = cand.userFloat("cosAlpha2D");
  }
  if (cand.hasUserFloat("cosAlpha3D")) {
    cosAlpha3D_ = cand.userFloat("cosAlpha3D");
  }
  
  // Read track DCA variables
  if (cand.hasUserFloat("trackDCA_min")) {
    trackDCA_min_ = cand.userFloat("trackDCA_min");
  }
  if (cand.hasUserFloat("trackDCA_max")) {
    trackDCA_max_ = cand.userFloat("trackDCA_max");
  }
  if (cand.hasUserFloat("trackDCA_avg")) {
    trackDCA_avg_ = cand.userFloat("trackDCA_avg");
  }
  
  // Read PV error variables
  if (cand.hasUserFloat("pvXError")) {
    pvXError_ = cand.userFloat("pvXError");
  }
  if (cand.hasUserFloat("pvYError")) {
    pvYError_ = cand.userFloat("pvYError");
  }
  if (cand.hasUserFloat("pvZError")) {
    pvZError_ = cand.userFloat("pvZError");
  }
  
  // Set decay mode identifier based on source name
  if (src.name.find("4Pi") != std::string::npos || src.name.find("ChiC_4Pi") != std::string::npos) {
    chicDecay_ = 0;  // 4pi
  } else if (src.name.find("DiKaon") != std::string::npos || src.name.find("ChiC_DiKaon") != std::string::npos) {
    chicDecay_ = 1;  // di-kaon
  } else if (src.name.find("KsKs") != std::string::npos || src.name.find("ChiC_KsKs") != std::string::npos) {
    chicDecay_ = 2;  // di-Ks
  } else if (src.name.find("PP") != std::string::npos || src.name.find("EtaC_PP") != std::string::npos) {
    chicDecay_ = 3;  // p+p-
  } else {
    chicDecay_ = -1;  // unknown
  }

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
    const auto& pvPos = primaryVertex->position();
    const double dx = vtx.x() - pvPos.x();
    const double dy = vtx.y() - pvPos.y();
    const double dz = vtx.z() - pvPos.z();
    
    candD3D_ = static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
    candDecayLength3D_ = candD3D_;
    candDecayLength2D_ = static_cast<float>(std::sqrt(dx * dx + dy * dy));
    
    // Pointing angle calculation
    const math::XYZVector flightDir(dx, dy, dz);
    const math::XYZVector candMom(cand.px(), cand.py(), cand.pz());
    const double flightDirMag = std::sqrt(flightDir.mag2());
    const double candMomMag = std::sqrt(candMom.mag2());
    
    if (flightDirMag > 0 && candMomMag > 0) {
      const double cosAngle3D = flightDir.Dot(candMom) / (flightDirMag * candMomMag);
      candCosPointingAngle3D_ = static_cast<float>(cosAngle3D);
      candPointingAngle3D_ = static_cast<float>(std::acos(std::max(-1.0, std::min(1.0, cosAngle3D))));
      
      const math::XYZVector flightDir2D(dx, dy, 0);
      const math::XYZVector candMom2D(cand.px(), cand.py(), 0);
      const double flightDir2DMag = std::sqrt(flightDir2D.mag2());
      const double candMom2DMag = std::sqrt(candMom2D.mag2());
      if (flightDir2DMag > 0 && candMom2DMag > 0) {
        const double cosAngle2D = flightDir2D.Dot(candMom2D) / (flightDir2DMag * candMom2DMag);
        candCosPointingAngle2D_ = static_cast<float>(cosAngle2D);
        candPointingAngle2D_ = static_cast<float>(std::acos(std::max(-1.0, std::min(1.0, cosAngle2D))));
      }
    }
    
    // Impact parameters (using best track if available)
    const reco::Track* candTrack = cand.bestTrack();
    if (candTrack) {
      candDxy_ = static_cast<float>(candTrack->dxy(pvPos));
      candDz_ = static_cast<float>(candTrack->dz(pvPos));
      const double dxyErr = candTrack->dxyError();
      const double dzErr = candTrack->dzError();
      if (dxyErr > 0) candDxySig_ = static_cast<float>(std::abs(candDxy_) / dxyErr);
      if (dzErr > 0) candDzSig_ = static_cast<float>(std::abs(candDz_) / dzErr);
    }
  } else {
    candD3D_ = -1.f;
    candDecayLength3D_ = -1.f;
    candDecayLength2D_ = -1.f;
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

  // Ensure we don't exceed array bounds
  if (daughters.size() > 4) {
    daughters.resize(4);
  }

  std::array<const reco::Candidate*, 4> orderedDaughters{{nullptr, nullptr, nullptr, nullptr}};
  const std::size_t nDaughters = std::min(daughters.size(), static_cast<std::size_t>(4));
  for (std::size_t i = 0; i < nDaughters; ++i) {
    const auto& info = daughters[i];
    orderedDaughters[i] = info.cand;

    dauPt_[i] = static_cast<float>(info.cand->pt());
    dauEta_[i] = static_cast<float>(info.cand->eta());
    dauPhi_[i] = static_cast<float>(info.cand->phi());
    dauCharge_[i] = info.cand->charge();
    dauMass_[i] = static_cast<float>(info.cand->mass());

    if (info.track) {
      dauChi2_[i] = static_cast<float>(info.track->normalizedChi2());
      dauNhits_[i] = static_cast<int>(info.track->numberOfValidHits());
      dauNpixHits_[i] = static_cast<int>(info.track->hitPattern().numberOfValidPixelHits());
      dauPtErr_[i] = static_cast<float>(info.track->ptError());
      dauEtaErr_[i] = static_cast<float>(info.track->etaError());
      dauPhiErr_[i] = static_cast<float>(info.track->phiError());
      
      if (primaryVertex && info.track) {
        const auto& pvPos = primaryVertex->position();
        const double dxy = info.track->dxy(pvPos);
        const double dz = info.track->dz(pvPos);
        dauDxy_[i] = static_cast<float>(dxy);
        dauDz_[i] = static_cast<float>(dz);
        dauD3d_[i] = static_cast<float>(std::sqrt(dxy * dxy + dz * dz));
        
        const double dxyErr = info.track->dxyError();
        const double dzErr = info.track->dzError();
        if (dxyErr > 0) {
          dauDxySig_[i] = static_cast<float>(std::abs(dxy) / dxyErr);
        }
        if (dzErr > 0) {
          dauDzSig_[i] = static_cast<float>(std::abs(dz) / dzErr);
        }
        
        // Compute dx and dy from dxy and phi angle (safer than referencePoint)
        // dxy is signed impact parameter, phi is track azimuthal angle
        const double phi = info.track->phi();
        const double dx = dxy * std::sin(phi);
        const double dy = -dxy * std::cos(phi);
        
        // Use dxy error as approximation for dx and dy errors
        // More accurate would require full covariance matrix diagonalization
        if (dxyErr > 0) {
          // Approximate: dxErr ≈ dxyErr * |sin(phi)|, dyErr ≈ dxyErr * |cos(phi)|
          const double dxErr = dxyErr * std::abs(std::sin(phi));
          const double dyErr = dxyErr * std::abs(std::cos(phi));
          if (dxErr > 1e-10) dauDxSig_[i] = static_cast<float>(std::abs(dx) / dxErr);
          if (dyErr > 1e-10) dauDySig_[i] = static_cast<float>(std::abs(dy) / dyErr);
        }
      }
      
      // Try to read userFloat values from producer (more accurate, computed with vertex fitting)
      if (i < 4) {  // Safety check
        // Build keys efficiently
        const char digit = char('1' + i);
        std::string longKey = "dauLongImpactSig_";
        longKey += digit;
        std::string transKey = "dauTransImpactSig_";
        transKey += digit;
        
        if (cand.hasUserFloat(longKey)) {
          dauDzSig_[i] = cand.userFloat(longKey);
        }
        if (cand.hasUserFloat(transKey)) {
          dauDxySig_[i] = cand.userFloat(transKey);
        }
      }
      
      // dE/dx not directly available from PAT candidates without ValueMap
      // Set to -1 (can be filled later if dE/dx ValueMap is provided)
      dauDeDx_[i] = -1.f;
    } else {
      dauChi2_[i] = -1.f;
      dauNhits_[i] = 0;
      dauNpixHits_[i] = 0;
      dauPtErr_[i] = -1.f;
      dauEtaErr_[i] = -1.f;
      dauPhiErr_[i] = -1.f;
      dauDeDx_[i] = -1.f;
    }
  }

  // Compute pair DCA and invariant masses
  const auto computePairInfo = [](const reco::Candidate* first, const reco::Candidate* second) -> std::pair<float, float> {
    if (!first || !second)
      return std::make_pair(0.f, -1.f);
    
    const auto& v1 = first->vertex();
    const auto& v2 = second->vertex();
    const double dx = v1.x() - v2.x();
    const double dy = v1.y() - v2.y();
    const double dz = v1.z() - v2.z();
    const float dca = static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
    
    // Invariant mass of the pair
    const double px = first->px() + second->px();
    const double py = first->py() + second->py();
    const double pz = first->pz() + second->pz();
    const double e = first->energy() + second->energy();
    const float mass = static_cast<float>(std::sqrt(e * e - px * px - py * py - pz * pz));
    
    return std::make_pair(dca, mass);
  };
  
  const std::array<std::pair<const reco::Candidate*, const reco::Candidate*>, 6> pairs = {{
    {orderedDaughters[0], orderedDaughters[1]},
    {orderedDaughters[0], orderedDaughters[2]},
    {orderedDaughters[0], orderedDaughters[3]},
    {orderedDaughters[1], orderedDaughters[2]},
    {orderedDaughters[1], orderedDaughters[3]},
    {orderedDaughters[2], orderedDaughters[3]}
  }};
  
  // Compute q2 for each pair: q2 = (1/N) * sum(exp(2i*phi)) for the pair
  const auto computePairQ2 = [](const reco::Candidate* first, const reco::Candidate* second) -> float {
    if (!first || !second)
      return 0.f;
    const double phi1 = first->phi();
    const double phi2 = second->phi();
    const double q2x = 0.5 * (std::cos(2.0 * phi1) + std::cos(2.0 * phi2));
    const double q2y = 0.5 * (std::sin(2.0 * phi1) + std::sin(2.0 * phi2));
    return static_cast<float>(std::sqrt(q2x * q2x + q2y * q2y));
  };
  
  for (size_t i = 0; i < 6; ++i) {
    const auto& pair = pairs[i];
    const auto info = computePairInfo(pair.first, pair.second);
    pairDca_[i] = info.first;
    pairMass_[i] = info.second;
    pairQ2_[i] = computePairQ2(pair.first, pair.second);
    
    if (pair.first && pair.second) {
      const double px = pair.first->px() + pair.second->px();
      const double py = pair.first->py() + pair.second->py();
      pairPt_[i] = static_cast<float>(std::sqrt(px * px + py * py));
      const double pz = pair.first->pz() + pair.second->pz();
      const double p = std::sqrt(px * px + py * py + pz * pz);
      if (p > 0) {
        pairEta_[i] = static_cast<float>(0.5 * std::log((p + pz) / (p - pz)));
        pairPhi_[i] = static_cast<float>(std::atan2(py, px));
      }
    }
  }

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
  if (genParticles && genParticles->isValid() && useGenMatching_) {
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
