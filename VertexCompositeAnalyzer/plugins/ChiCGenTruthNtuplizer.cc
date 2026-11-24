#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"

namespace {
  constexpr float kMissing = std::numeric_limits<float>::quiet_NaN();

  int sourceInt(const pat::CompositeCandidate& candidate, const std::string& name, int fallback = -1) {
    return candidate.hasUserInt(name) ? candidate.userInt(name) : fallback;
  }

  float sourceFloat(const pat::CompositeCandidate& candidate, const std::string& name, float fallback = kMissing) {
    return candidate.hasUserFloat(name) ? candidate.userFloat(name) : fallback;
  }

  bool isChiC(int pdgId) {
    const int id = std::abs(pdgId);
    return id == 10441 || id == 20443 || id == 445;
  }
}  // namespace

// MC-only companion tree. It keeps the generator decay graph and records
// every matching decision, so unmatched candidates are diagnosable rather
// than silently labelled as wrong objects. Association matching is preferred
// for packed candidates. Photon matching is restricted to the radiative
// sibling of the matched J/psi in the same chi_c decay graph, preventing a
// nearby photon from another chi_c or from J/psi final-state radiation from
// being labelled as the signal photon.
class ChiCGenTruthNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCGenTruthNtuplizer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  using GenAssociation = edm::Association<reco::GenParticleCollection>;
  using CandidateCollection = pat::CompositeCandidateCollection;
  struct Match {
    int index = -1;
    float deltaR = kMissing;
    float ptRatio = kMissing;
  };
  struct JpsiMatch {
    Match jpsi;
    std::array<Match, 2> muons;
  };
  // associationStatus: 0=no association, 1=association in focused ChiC graph,
  // 2=association to another pruned-gen particle (normally an unrelated track).
  struct AssociationMatch {
    const reco::GenParticle* particle = nullptr;
    int sourceIndex = -1;
    int outputIndex = -1;
    int status = 0;
  };

  Match matchPdg(const reco::Candidate&, int absPdgId, float maxDeltaR, bool requirePt) const;
  Match matchRadiativePhoton(const reco::Candidate&, int chiIndex, float maxDeltaR, bool requirePt) const;
  Match matchMuonFromJpsi(const reco::Candidate&, int jpsiIndex) const;
  JpsiMatch matchJpsi(const pat::CompositeCandidate&) const;
  bool hasAncestor(int childIndex, int ancestorIndex) const;
  int chiCAncestor(int) const;
  int commonChiCAncestor(int, int) const;
  bool isRadiativeChiPhoton(int, int) const;
  AssociationMatch associatedPackedGen(int,
                                       const edm::Handle<pat::PackedCandidateCollection>&,
                                       const edm::Handle<GenAssociation>&) const;
  AssociationMatch associatedLegGen(const pat::CompositeCandidate&,
                                    const edm::Handle<pat::PackedCandidateCollection>&,
                                    const edm::Handle<pat::PackedCandidateCollection>&,
                                    const edm::Handle<GenAssociation>&,
                                    const edm::Handle<GenAssociation>&) const;
  void clear();

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<CandidateCollection> oniaToken_;
  edm::EDGetTokenT<CandidateCollection> conventionalToken_;
  edm::EDGetTokenT<CandidateCollection> oneProngToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::EDGetTokenT<GenAssociation> packedAssociationToken_;
  edm::EDGetTokenT<GenAssociation> lostAssociationToken_;
  bool doOneProng_;
  TTree* tree_ = nullptr;
  const reco::GenParticleCollection* gen_ = nullptr;
  std::unordered_map<const reco::Candidate*, int> genIndex_;
  std::vector<const reco::Candidate*> selectedGen_;

  std::vector<unsigned int> eventRun_, eventLumi_, eventNumber_;
  std::vector<int> genPdgId_, genStatus_, genCharge_, genIsLastCopy_, genIsPromptFinalState_, genFromHardProcess_;
  std::vector<float> genPx_, genPy_, genPz_, genE_, genMass_, genVX_, genVY_, genVZ_;
  std::vector<unsigned int> genMotherBegin_, genMotherSize_, genDaughterBegin_, genDaughterSize_;
  std::vector<int> genMotherIndex_, genDaughterIndex_;
  std::vector<int> truthMuonJpsiIndex_, truthMuonSlot_, truthMuonGenIndex_, truthMuonGenPdgId_, truthMuonMatchStatus_;
  std::vector<float> truthMuonDeltaR_, truthMuonPtRatio_;
  std::vector<int> truthOniaIndex_, truthOniaGenIndex_, truthOniaChiGenIndex_, truthOniaMatchStatus_;
  std::vector<float> truthOniaDeltaR_, truthOniaPtRatio_;
  std::vector<int> truthRecoChiIndex_, truthRecoChiKind_, truthJpsiGenIndex_, truthPhotonGenIndex_, truthChiGenIndex_,
      truthObservedLegGenIndex_, truthObservedLegGenPdgId_, truthObservedLegMatchMethod_, truthMatchStatus_;
  std::vector<int> truthMatchIsBest_, truthMatchMultiplicity_;
  // photonMatchMethod: 0=none, 1=packed association, 2=restricted P4 match,
  // 3=restricted one-prong direction match.
  std::vector<int> truthPhotonMatchMethod_, truthPhotonAssociationStatus_, truthPhotonAssociatedGenSourceIndex_,
      truthPhotonAssociatedGenPdgId_;
  std::vector<int> truthObservedLegSourceCollection_, truthObservedLegSourceKey_, truthObservedLegRecoPdgId_,
      truthObservedLegCharge_, truthObservedLegAssociationStatus_, truthObservedLegAssociatedGenSourceIndex_,
      truthObservedLegAssociatedGenPdgId_, truthMissingPartnerStatus_;
  std::vector<float> truthJpsiDeltaR_, truthPhotonDeltaR_, truthPhotonPtRatio_, truthObservedLegDeltaR_,
      truthObservedLegPtRatio_, truthObservedLegEnergyFraction_, truthObservedLegAssociatedGenPx_,
      truthObservedLegAssociatedGenPy_, truthObservedLegAssociatedGenPz_, truthObservedLegAssociatedGenE_,
      truthObservedLegAssociatedGenVX_, truthObservedLegAssociatedGenVY_, truthObservedLegAssociatedGenVZ_;
};

ChiCGenTruthNtuplizer::ChiCGenTruthNtuplizer(const edm::ParameterSet& cfg)
    : genToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
      oniaToken_(consumes<CandidateCollection>(cfg.getParameter<edm::InputTag>("oniaStates"))),
      conventionalToken_(consumes<CandidateCollection>(cfg.getParameter<edm::InputTag>("conventionalCandidates"))),
      oneProngToken_(mayConsume<CandidateCollection>(cfg.getParameter<edm::InputTag>("oneProngCandidates"))),
      packedToken_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostToken_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      packedAssociationToken_(
          mayConsume<GenAssociation>(cfg.getParameter<edm::InputTag>("packedToGenAssociation"))),
      lostAssociationToken_(mayConsume<GenAssociation>(cfg.getParameter<edm::InputTag>("lostToGenAssociation"))),
      doOneProng_(cfg.getParameter<bool>("doOneProng")) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> file;
  tree_ = file->make<TTree>("ChiCGenTruth", "MiniAODSIM generator decay graph and reco-truth matching");
#define BRANCH(name) tree_->Branch(#name, &name)
  BRANCH(eventRun_); BRANCH(eventLumi_); BRANCH(eventNumber_);
  BRANCH(genPdgId_); BRANCH(genStatus_); BRANCH(genCharge_); BRANCH(genIsLastCopy_); BRANCH(genIsPromptFinalState_); BRANCH(genFromHardProcess_); BRANCH(genPx_); BRANCH(genPy_); BRANCH(genPz_); BRANCH(genE_); BRANCH(genMass_); BRANCH(genVX_); BRANCH(genVY_); BRANCH(genVZ_); BRANCH(genMotherBegin_); BRANCH(genMotherSize_); BRANCH(genMotherIndex_); BRANCH(genDaughterBegin_); BRANCH(genDaughterSize_); BRANCH(genDaughterIndex_);
  BRANCH(truthMuonJpsiIndex_); BRANCH(truthMuonSlot_); BRANCH(truthMuonGenIndex_); BRANCH(truthMuonGenPdgId_); BRANCH(truthMuonMatchStatus_); BRANCH(truthMuonDeltaR_); BRANCH(truthMuonPtRatio_);
  BRANCH(truthOniaIndex_); BRANCH(truthOniaGenIndex_); BRANCH(truthOniaChiGenIndex_); BRANCH(truthOniaMatchStatus_); BRANCH(truthOniaDeltaR_); BRANCH(truthOniaPtRatio_);
  BRANCH(truthRecoChiIndex_); BRANCH(truthRecoChiKind_); BRANCH(truthJpsiGenIndex_); BRANCH(truthPhotonGenIndex_); BRANCH(truthChiGenIndex_); BRANCH(truthMatchIsBest_); BRANCH(truthMatchMultiplicity_); BRANCH(truthPhotonMatchMethod_); BRANCH(truthPhotonAssociationStatus_); BRANCH(truthPhotonAssociatedGenSourceIndex_); BRANCH(truthPhotonAssociatedGenPdgId_); BRANCH(truthObservedLegGenIndex_); BRANCH(truthObservedLegGenPdgId_); BRANCH(truthObservedLegMatchMethod_); BRANCH(truthMatchStatus_); BRANCH(truthObservedLegSourceCollection_); BRANCH(truthObservedLegSourceKey_); BRANCH(truthObservedLegRecoPdgId_); BRANCH(truthObservedLegCharge_); BRANCH(truthObservedLegAssociationStatus_); BRANCH(truthObservedLegAssociatedGenSourceIndex_); BRANCH(truthObservedLegAssociatedGenPdgId_); BRANCH(truthMissingPartnerStatus_); BRANCH(truthJpsiDeltaR_); BRANCH(truthPhotonDeltaR_); BRANCH(truthPhotonPtRatio_); BRANCH(truthObservedLegDeltaR_); BRANCH(truthObservedLegPtRatio_); BRANCH(truthObservedLegEnergyFraction_); BRANCH(truthObservedLegAssociatedGenPx_); BRANCH(truthObservedLegAssociatedGenPy_); BRANCH(truthObservedLegAssociatedGenPz_); BRANCH(truthObservedLegAssociatedGenE_); BRANCH(truthObservedLegAssociatedGenVX_); BRANCH(truthObservedLegAssociatedGenVY_); BRANCH(truthObservedLegAssociatedGenVZ_);
#undef BRANCH
}

void ChiCGenTruthNtuplizer::clear() {
#define CLEAR(name) name.clear()
  CLEAR(eventRun_); CLEAR(eventLumi_); CLEAR(eventNumber_); CLEAR(genPdgId_); CLEAR(genStatus_); CLEAR(genCharge_); CLEAR(genIsLastCopy_); CLEAR(genIsPromptFinalState_); CLEAR(genFromHardProcess_); CLEAR(genPx_); CLEAR(genPy_); CLEAR(genPz_); CLEAR(genE_); CLEAR(genMass_); CLEAR(genVX_); CLEAR(genVY_); CLEAR(genVZ_); CLEAR(genMotherBegin_); CLEAR(genMotherSize_); CLEAR(genMotherIndex_); CLEAR(genDaughterBegin_); CLEAR(genDaughterSize_); CLEAR(genDaughterIndex_); CLEAR(truthMuonJpsiIndex_); CLEAR(truthMuonSlot_); CLEAR(truthMuonGenIndex_); CLEAR(truthMuonGenPdgId_); CLEAR(truthMuonMatchStatus_); CLEAR(truthMuonDeltaR_); CLEAR(truthMuonPtRatio_); CLEAR(truthOniaIndex_); CLEAR(truthOniaGenIndex_); CLEAR(truthOniaChiGenIndex_); CLEAR(truthOniaMatchStatus_); CLEAR(truthOniaDeltaR_); CLEAR(truthOniaPtRatio_); CLEAR(truthRecoChiIndex_); CLEAR(truthRecoChiKind_); CLEAR(truthJpsiGenIndex_); CLEAR(truthPhotonGenIndex_); CLEAR(truthChiGenIndex_); CLEAR(truthMatchIsBest_); CLEAR(truthMatchMultiplicity_); CLEAR(truthPhotonMatchMethod_); CLEAR(truthPhotonAssociationStatus_); CLEAR(truthPhotonAssociatedGenSourceIndex_); CLEAR(truthPhotonAssociatedGenPdgId_); CLEAR(truthObservedLegGenIndex_); CLEAR(truthObservedLegGenPdgId_); CLEAR(truthObservedLegMatchMethod_); CLEAR(truthMatchStatus_); CLEAR(truthObservedLegSourceCollection_); CLEAR(truthObservedLegSourceKey_); CLEAR(truthObservedLegRecoPdgId_); CLEAR(truthObservedLegCharge_); CLEAR(truthObservedLegAssociationStatus_); CLEAR(truthObservedLegAssociatedGenSourceIndex_); CLEAR(truthObservedLegAssociatedGenPdgId_); CLEAR(truthMissingPartnerStatus_); CLEAR(truthJpsiDeltaR_); CLEAR(truthPhotonDeltaR_); CLEAR(truthPhotonPtRatio_); CLEAR(truthObservedLegDeltaR_); CLEAR(truthObservedLegPtRatio_); CLEAR(truthObservedLegEnergyFraction_); CLEAR(truthObservedLegAssociatedGenPx_); CLEAR(truthObservedLegAssociatedGenPy_); CLEAR(truthObservedLegAssociatedGenPz_); CLEAR(truthObservedLegAssociatedGenE_); CLEAR(truthObservedLegAssociatedGenVX_); CLEAR(truthObservedLegAssociatedGenVY_); CLEAR(truthObservedLegAssociatedGenVZ_);
#undef CLEAR
}

ChiCGenTruthNtuplizer::Match ChiCGenTruthNtuplizer::matchPdg(
    const reco::Candidate& reco, int absPdgId, float maxDeltaR, bool requirePt) const {
  Match best;
  if (gen_ == nullptr)
    return best;
  for (size_t i = 0; i < gen_->size(); ++i) {
    const auto& gen = gen_->at(i);
    if (std::abs(gen.pdgId()) != absPdgId || gen.pt() <= 0.f)
      continue;
    const float dr = reco::deltaR(reco, gen);
    const float ratio = reco.pt() / gen.pt();
    if (dr > maxDeltaR || (requirePt && (ratio < 0.2f || ratio > 5.f)) || (best.index >= 0 && dr >= best.deltaR))
      continue;
    const auto output = genIndex_.find(&gen);
    if (output == genIndex_.end())
      continue;
    best = {output->second, dr, ratio};
  }
  return best;
}

bool ChiCGenTruthNtuplizer::hasAncestor(int childIndex, int ancestorIndex) const {
  if (childIndex < 0 || ancestorIndex < 0 ||
      static_cast<size_t>(childIndex) >= selectedGen_.size() ||
      static_cast<size_t>(ancestorIndex) >= selectedGen_.size())
    return false;
  std::vector<const reco::Candidate*> pending{selectedGen_.at(childIndex)};
  std::unordered_set<const reco::Candidate*> visited;
  while (!pending.empty()) {
    const auto* current = pending.back();
    pending.pop_back();
    if (current == nullptr || !visited.insert(current).second)
      continue;
    const auto found = genIndex_.find(current);
    if (found != genIndex_.end() && found->second == ancestorIndex)
      return true;
    for (size_t mother = 0; mother < current->numberOfMothers(); ++mother)
      pending.push_back(current->mother(mother));
  }
  return false;
}

ChiCGenTruthNtuplizer::Match ChiCGenTruthNtuplizer::matchMuonFromJpsi(
    const reco::Candidate& reco, int jpsiIndex) const {
  Match best;
  for (const auto* particle : selectedGen_) {
    const auto& gen = *static_cast<const reco::GenParticle*>(particle);
    const int muonIndex = genIndex_.at(particle);
    if (std::abs(gen.pdgId()) != 13 || gen.charge() != reco.charge() || gen.pt() <= 0.f ||
        !hasAncestor(muonIndex, jpsiIndex))
      continue;
    const float deltaR = reco::deltaR(reco, gen);
    const float ptRatio = reco.pt() / gen.pt();
    if (deltaR > .03f || ptRatio < .2f || ptRatio > 5.f || (best.index >= 0 && deltaR >= best.deltaR))
      continue;
    best = {muonIndex, deltaR, ptRatio};
  }
  return best;
}

ChiCGenTruthNtuplizer::JpsiMatch ChiCGenTruthNtuplizer::matchJpsi(const pat::CompositeCandidate& reco) const {
  JpsiMatch best;
  const auto* muon0 = reco.daughter("muon1");
  const auto* muon1 = reco.daughter("muon2");
  if (muon0 == nullptr || muon1 == nullptr)
    return best;
  float bestScore = std::numeric_limits<float>::infinity();
  for (const auto* particle : selectedGen_) {
    const auto& gen = *static_cast<const reco::GenParticle*>(particle);
    const int jpsiIndex = genIndex_.at(particle);
    if (std::abs(gen.pdgId()) != 443 || gen.pt() <= 0.f)
      continue;
    const float deltaR = reco::deltaR(reco, gen);
    const float ptRatio = reco.pt() / gen.pt();
    if (deltaR > .05f || ptRatio < .2f || ptRatio > 5.f)
      continue;
    const Match muonMatch0 = matchMuonFromJpsi(*muon0, jpsiIndex);
    const Match muonMatch1 = matchMuonFromJpsi(*muon1, jpsiIndex);
    if (muonMatch0.index < 0 || muonMatch1.index < 0 || muonMatch0.index == muonMatch1.index)
      continue;
    const float score = deltaR + muonMatch0.deltaR + muonMatch1.deltaR;
    if (score >= bestScore)
      continue;
    bestScore = score;
    best.jpsi = {jpsiIndex, deltaR, ptRatio};
    best.muons = {muonMatch0, muonMatch1};
  }
  return best;
}

int ChiCGenTruthNtuplizer::chiCAncestor(int index) const {
  if (index < 0 || static_cast<size_t>(index) >= selectedGen_.size())
    return -1;
  std::vector<const reco::Candidate*> pending{selectedGen_.at(index)};
  std::unordered_set<const reco::Candidate*> visited;
  while (!pending.empty()) {
    const auto* current = pending.back();
    pending.pop_back();
    if (current == nullptr || !visited.insert(current).second)
      continue;
    const auto found = genIndex_.find(current);
    if (found != genIndex_.end() && isChiC(current->pdgId()))
      return found->second;
    for (size_t mother = 0; mother < current->numberOfMothers(); ++mother)
      pending.push_back(current->mother(mother));
  }
  return -1;
}

int ChiCGenTruthNtuplizer::commonChiCAncestor(int first, int second) const {
  const int firstAncestor = chiCAncestor(first);
  const int secondAncestor = chiCAncestor(second);
  return firstAncestor >= 0 && firstAncestor == secondAncestor ? firstAncestor : -1;
}

bool ChiCGenTruthNtuplizer::isRadiativeChiPhoton(int photon, int chi) const {
  if (photon < 0 || chi < 0 || static_cast<size_t>(photon) >= selectedGen_.size())
    return false;
  struct Node {
    const reco::Candidate* particle;
    bool passedJpsi;
  };
  std::vector<Node> pending{{selectedGen_.at(photon), false}};
  std::unordered_set<const reco::Candidate*> visitedWithoutJpsi;
  std::unordered_set<const reco::Candidate*> visitedWithJpsi;
  while (!pending.empty()) {
    const auto node = pending.back();
    pending.pop_back();
    if (node.particle == nullptr)
      continue;
    auto& visited = node.passedJpsi ? visitedWithJpsi : visitedWithoutJpsi;
    if (!visited.insert(node.particle).second)
      continue;
    const auto found = genIndex_.find(node.particle);
    if (found != genIndex_.end() && found->second == chi)
      return !node.passedJpsi;
    const bool passedJpsi = node.passedJpsi || std::abs(node.particle->pdgId()) == 443;
    for (size_t mother = 0; mother < node.particle->numberOfMothers(); ++mother)
      pending.push_back({node.particle->mother(mother), passedJpsi});
  }
  return false;
}

ChiCGenTruthNtuplizer::Match ChiCGenTruthNtuplizer::matchRadiativePhoton(
    const reco::Candidate& reco, int chiIndex, float maxDeltaR, bool requirePt) const {
  Match best;
  if (chiIndex < 0)
    return best;
  for (const auto* particle : selectedGen_) {
    const auto& gen = *static_cast<const reco::GenParticle*>(particle);
    const int photonIndex = genIndex_.at(particle);
    if (gen.pdgId() != 22 || gen.pt() <= 0.f || chiCAncestor(photonIndex) != chiIndex ||
        !isRadiativeChiPhoton(photonIndex, chiIndex))
      continue;
    const float deltaR = reco::deltaR(reco, gen);
    const float ptRatio = reco.pt() / gen.pt();
    if (deltaR > maxDeltaR || (requirePt && (ptRatio < 0.2f || ptRatio > 5.f)) ||
        (best.index >= 0 && deltaR >= best.deltaR))
      continue;
    best = {photonIndex, deltaR, ptRatio};
  }
  return best;
}

ChiCGenTruthNtuplizer::AssociationMatch ChiCGenTruthNtuplizer::associatedPackedGen(
    int key,
    const edm::Handle<pat::PackedCandidateCollection>& packed,
    const edm::Handle<GenAssociation>& packedAssociation) const {
  AssociationMatch match;
  if (!packed.isValid() || !packedAssociation.isValid() || key < 0 ||
      static_cast<size_t>(key) >= packed->size() || !packedAssociation->contains(packed.id()))
    return match;
  const reco::GenParticleRef result =
      (*packedAssociation)[edm::Ref<pat::PackedCandidateCollection>(packed, key)];
  if (result.isNull())
    return match;
  match.particle = result.get();
  match.sourceIndex = result.key();
  const auto found = genIndex_.find(result.get());
  if (found == genIndex_.end()) {
    match.status = 2;
    return match;
  }
  match.outputIndex = found->second;
  match.status = 1;
  return match;
}

ChiCGenTruthNtuplizer::AssociationMatch ChiCGenTruthNtuplizer::associatedLegGen(
    const pat::CompositeCandidate& leg,
    const edm::Handle<pat::PackedCandidateCollection>& packed,
    const edm::Handle<pat::PackedCandidateCollection>& lost,
    const edm::Handle<GenAssociation>& packedAssociation,
    const edm::Handle<GenAssociation>& lostAssociation) const {
  AssociationMatch match;
  const int source = sourceInt(leg, "sourceCollection"), key = sourceInt(leg, "sourceKey");
  reco::GenParticleRef result;
  if (source == 0 && packed.isValid() && packedAssociation.isValid() && key >= 0 &&
      static_cast<size_t>(key) < packed->size() && packedAssociation->contains(packed.id())) {
    result = (*packedAssociation)[edm::Ref<pat::PackedCandidateCollection>(packed, key)];
  }
  if (source == 1 && lost.isValid() && lostAssociation.isValid() && key >= 0 &&
      static_cast<size_t>(key) < lost->size() && lostAssociation->contains(lost.id())) {
    result = (*lostAssociation)[edm::Ref<pat::PackedCandidateCollection>(lost, key)];
  }
  if (result.isNull())
    return match;
  match.particle = result.get();
  match.sourceIndex = result.key();
  const auto found = genIndex_.find(result.get());
  if (found == genIndex_.end()) {
    match.status = 2;
    return match;
  }
  match.outputIndex = found->second;
  match.status = 1;
  return match;
}

void ChiCGenTruthNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  clear();
  const auto& gen = event.get(genToken_);
  const auto& onias = event.get(oniaToken_);
  const auto& conventional = event.get(conventionalToken_);
  const auto oneProngHandle = doOneProng_ ? event.getHandle(oneProngToken_) : edm::Handle<CandidateCollection>{};
  const auto packed = event.getHandle(packedToken_);
  const auto lost = event.getHandle(lostToken_);
  const auto packedAssociation = event.getHandle(packedAssociationToken_);
  const auto lostAssociation = event.getHandle(lostAssociationToken_);

  gen_ = &gen;
  genIndex_.clear();
  selectedGen_.clear();
  auto fromChiC = [](const reco::Candidate* particle) {
    std::vector<const reco::Candidate*> pending{particle};
    std::unordered_set<const reco::Candidate*> visited;
    while (!pending.empty()) {
      const auto* current = pending.back();
      pending.pop_back();
      if (current == nullptr || !visited.insert(current).second)
        continue;
      if (isChiC(current->pdgId()))
        return true;
      for (size_t mother = 0; mother < current->numberOfMothers(); ++mother)
        pending.push_back(current->mother(mother));
    }
    return false;
  };
  for (const auto& p : gen) {
    if (fromChiC(&p)) {
      genIndex_[&p] = static_cast<int>(selectedGen_.size());
      selectedGen_.push_back(&p);
    }
  }
  for (const auto* ptr : selectedGen_) {
    const auto& p = *static_cast<const reco::GenParticle*>(ptr);
    genPdgId_.push_back(p.pdgId());
    genStatus_.push_back(p.status());
    genCharge_.push_back(p.charge());
    genIsLastCopy_.push_back(p.statusFlags().isLastCopy());
    genIsPromptFinalState_.push_back(p.statusFlags().isPrompt());
    genFromHardProcess_.push_back(p.statusFlags().fromHardProcess());
    genPx_.push_back(p.px());
    genPy_.push_back(p.py());
    genPz_.push_back(p.pz());
    genE_.push_back(p.energy());
    genMass_.push_back(p.mass());
    genVX_.push_back(p.vx());
    genVY_.push_back(p.vy());
    genVZ_.push_back(p.vz());
    genMotherBegin_.push_back(genMotherIndex_.size());
    for (size_t j = 0; j < p.numberOfMothers(); ++j) {
      auto it = genIndex_.find(p.mother(j));
      if (it != genIndex_.end())
        genMotherIndex_.push_back(it->second);
    }
    genMotherSize_.push_back(genMotherIndex_.size() - genMotherBegin_.back());
    genDaughterBegin_.push_back(genDaughterIndex_.size());
    for (size_t j = 0; j < p.numberOfDaughters(); ++j) {
      auto it = genIndex_.find(p.daughter(j));
      if (it != genIndex_.end())
        genDaughterIndex_.push_back(it->second);
    }
    genDaughterSize_.push_back(genDaughterIndex_.size() - genDaughterBegin_.back());
  }

  std::vector<JpsiMatch> oniaMatches;
  oniaMatches.reserve(onias.size());
  for (size_t j = 0; j < onias.size(); ++j) {
    const JpsiMatch match = matchJpsi(onias[j]);
    oniaMatches.push_back(match);
    const int chiIndex = chiCAncestor(match.jpsi.index);
    const int oniaStatus = (match.jpsi.index >= 0 ? 1 : 0) |
                           (match.muons[0].index >= 0 && match.muons[1].index >= 0 ? 2 : 0) |
                           (chiIndex >= 0 ? 4 : 0);
    truthOniaIndex_.push_back(static_cast<int>(j));
    truthOniaGenIndex_.push_back(match.jpsi.index);
    truthOniaChiGenIndex_.push_back(chiIndex);
    truthOniaMatchStatus_.push_back(oniaStatus);
    truthOniaDeltaR_.push_back(match.jpsi.deltaR);
    truthOniaPtRatio_.push_back(match.jpsi.ptRatio);
    for (unsigned int slot = 0; slot != 2; ++slot) {
      const Match& muonMatch = match.muons[slot];
      truthMuonJpsiIndex_.push_back(static_cast<int>(j));
      truthMuonSlot_.push_back(slot);
      truthMuonGenIndex_.push_back(muonMatch.index);
      truthMuonGenPdgId_.push_back(muonMatch.index >= 0 ? selectedGen_[muonMatch.index]->pdgId() : 0);
      truthMuonMatchStatus_.push_back(muonMatch.index >= 0 ? 1 : 0);
      truthMuonDeltaR_.push_back(muonMatch.deltaR);
      truthMuonPtRatio_.push_back(muonMatch.ptRatio);
    }
  }

  auto appendCandidate = [&](const pat::CompositeCandidate& chi, int kind, int recoIndex) {
    const auto* jpsi = dynamic_cast<const pat::CompositeCandidate*>(chi.daughter(kind == 0 ? "dimuon" : "jpsi"));
    const auto* leg = kind == 1 ? dynamic_cast<const pat::CompositeCandidate*>(chi.daughter("leg")) : nullptr;
    const int candidateOniaIndex =
        sourceInt(chi, "oniaIndex", sourceInt(chi, "oniaIdx", sourceInt(chi, "jpsiIndex")));
    JpsiMatch jpsiMatch;
    if (candidateOniaIndex >= 0 && static_cast<size_t>(candidateOniaIndex) < oniaMatches.size())
      jpsiMatch = oniaMatches[candidateOniaIndex];
    else if (jpsi != nullptr)
      jpsiMatch = matchJpsi(*jpsi);
    const Match& jm = jpsiMatch.jpsi;
    Match gm;
    AssociationMatch photonAssociation;
    int photonMethod = 0;
    const int jpsiChi = chiCAncestor(jm.index);
    if (kind == 0) {
      const auto* photon = chi.daughter("photon");
      const int packedPhotonIndex = sourceInt(chi, "packedPhotonIndex");
      if (photon && packedPhotonIndex >= 0) {
        photonAssociation = associatedPackedGen(packedPhotonIndex, packed, packedAssociation);
        const int associatedIndex = photonAssociation.outputIndex;
        if (associatedIndex >= 0 && jpsiChi >= 0 && photonAssociation.particle->pdgId() == 22 &&
            chiCAncestor(associatedIndex) == jpsiChi && isRadiativeChiPhoton(associatedIndex, jpsiChi)) {
          gm = {associatedIndex,
                float(reco::deltaR(*photon, *photonAssociation.particle)),
                photonAssociation.particle->pt() > 0.f
                    ? float(photon->pt() / photonAssociation.particle->pt())
                    : kMissing};
          photonMethod = 1;
        }
      }
      if (photon && gm.index < 0 && jpsiChi >= 0) {
        gm = matchRadiativePhoton(*photon, jpsiChi, .05f, true);
        if (gm.index >= 0)
          photonMethod = 2;
      }
    } else {
      const float x = sourceFloat(chi, "freePhotonDirectionX"),
                  y = sourceFloat(chi, "freePhotonDirectionY"),
                  z = sourceFloat(chi, "freePhotonDirectionZ");
      const float pt = std::hypot(x, y), eta = pt > 0.f ? std::asinh(z / pt) : kMissing, phi = std::atan2(y, x);
      if (std::isfinite(eta)) {
        for (const auto* ptr : selectedGen_) {
          const auto& g = *static_cast<const reco::GenParticle*>(ptr);
          const int photonIndex = genIndex_.at(ptr);
          if (g.pdgId() != 22 || jpsiChi < 0 || chiCAncestor(photonIndex) != jpsiChi ||
              !isRadiativeChiPhoton(photonIndex, jpsiChi))
            continue;
          const float dr = reco::deltaR(eta, phi, float(g.eta()), float(g.phi()));
          if (gm.index < 0 || (std::isfinite(dr) && (!std::isfinite(gm.deltaR) || dr < gm.deltaR)))
            gm = {photonIndex, dr, kMissing};
        }
      }
      if (gm.index >= 0)
        photonMethod = 3;
    }

    int legGen = -1, legMethod = 0;
    Match lm;
    AssociationMatch association;
    if (leg) {
      association = associatedLegGen(*leg, packed, lost, packedAssociation, lostAssociation);
      if (association.particle != nullptr) {
        legMethod = 1;
        legGen = association.outputIndex;
        lm.deltaR = reco::deltaR(*leg, *association.particle);
        lm.ptRatio = association.particle->pt() > 0 ? leg->pt() / association.particle->pt() : kMissing;
      } else {
        lm = matchPdg(*leg, 11, .05f, true);
        legGen = lm.index;
        legMethod = legGen >= 0 ? 2 : 0;
      }
    }

    const int chiGen = commonChiCAncestor(jm.index, gm.index);
    const int status = (jm.index >= 0 ? 1 : 0) | (gm.index >= 0 ? 2 : 0) | (chiGen >= 0 ? 4 : 0) |
                       (kind == 1 && std::isfinite(gm.deltaR) && gm.deltaR < .05f ? 8 : 0) |
                       (legMethod == 1 ? 16 : 0) | (legMethod == 2 ? 32 : 0);
    truthRecoChiIndex_.push_back(recoIndex);
    truthRecoChiKind_.push_back(kind);
    truthJpsiGenIndex_.push_back(jm.index);
    truthPhotonGenIndex_.push_back(gm.index);
    truthChiGenIndex_.push_back(chiGen);
    truthPhotonMatchMethod_.push_back(photonMethod);
    truthPhotonAssociationStatus_.push_back(photonAssociation.status);
    truthPhotonAssociatedGenSourceIndex_.push_back(photonAssociation.sourceIndex);
    truthPhotonAssociatedGenPdgId_.push_back(photonAssociation.particle ? photonAssociation.particle->pdgId() : 0);
    truthObservedLegGenIndex_.push_back(legGen);
    truthObservedLegGenPdgId_.push_back(legGen >= 0 ? selectedGen_[legGen]->pdgId() : 0);
    truthObservedLegMatchMethod_.push_back(legMethod);
    truthMatchStatus_.push_back(status);
    truthObservedLegSourceCollection_.push_back(leg ? sourceInt(*leg, "sourceCollection") : -1);
    truthObservedLegSourceKey_.push_back(leg ? sourceInt(*leg, "sourceKey") : -1);
    truthObservedLegRecoPdgId_.push_back(leg ? leg->pdgId() : 0);
    truthObservedLegCharge_.push_back(leg ? leg->charge() : 0);
    truthObservedLegAssociationStatus_.push_back(association.status);
    truthObservedLegAssociatedGenSourceIndex_.push_back(association.sourceIndex);
    truthObservedLegAssociatedGenPdgId_.push_back(association.particle ? association.particle->pdgId() : 0);
    truthMissingPartnerStatus_.push_back(kind == 1 ? 1 : 0);
    truthJpsiDeltaR_.push_back(jm.deltaR);
    truthPhotonDeltaR_.push_back(gm.deltaR);
    truthPhotonPtRatio_.push_back(gm.ptRatio);
    truthObservedLegDeltaR_.push_back(lm.deltaR);
    truthObservedLegPtRatio_.push_back(lm.ptRatio);
    truthObservedLegEnergyFraction_.push_back(
        leg && gm.index >= 0 && selectedGen_[gm.index]->energy() > 0
            ? leg->energy() / selectedGen_[gm.index]->energy()
            : kMissing);
    truthObservedLegAssociatedGenPx_.push_back(association.particle ? association.particle->px() : kMissing);
    truthObservedLegAssociatedGenPy_.push_back(association.particle ? association.particle->py() : kMissing);
    truthObservedLegAssociatedGenPz_.push_back(association.particle ? association.particle->pz() : kMissing);
    truthObservedLegAssociatedGenE_.push_back(association.particle ? association.particle->energy() : kMissing);
    truthObservedLegAssociatedGenVX_.push_back(association.particle ? association.particle->vx() : kMissing);
    truthObservedLegAssociatedGenVY_.push_back(association.particle ? association.particle->vy() : kMissing);
    truthObservedLegAssociatedGenVZ_.push_back(association.particle ? association.particle->vz() : kMissing);
  };

  int recoIndex = 0;
  for (const auto& candidate : conventional)
    appendCandidate(candidate, 0, recoIndex++);
  if (doOneProng_ && oneProngHandle.isValid()) {
    for (const auto& candidate : *oneProngHandle)
      appendCandidate(candidate, 1, recoIndex++);
  }

  truthMatchIsBest_.assign(truthMatchStatus_.size(), 0);
  truthMatchMultiplicity_.assign(truthMatchStatus_.size(), 0);
  std::map<std::tuple<int, int, int>, std::vector<size_t>> matchGroups;
  for (size_t index = 0; index < truthMatchStatus_.size(); ++index) {
    if ((truthMatchStatus_[index] & 0x7) != 0x7)
      continue;
    matchGroups[{truthRecoChiKind_[index], truthJpsiGenIndex_[index], truthPhotonGenIndex_[index]}].push_back(index);
  }
  for (const auto& group : matchGroups) {
    const auto& indices = group.second;
    size_t bestIndex = indices.front();
    auto score = [&](size_t index) {
      const int methodRank = truthPhotonMatchMethod_[index] == 1 ? 0 : 1;
      const float deltaR = std::isfinite(truthPhotonDeltaR_[index]) ? truthPhotonDeltaR_[index]
                                                                    : std::numeric_limits<float>::infinity();
      const float ptResidual = std::isfinite(truthPhotonPtRatio_[index])
                                   ? std::abs(truthPhotonPtRatio_[index] - 1.f)
                                   : std::numeric_limits<float>::infinity();
      return std::make_tuple(methodRank, deltaR, ptResidual, truthRecoChiIndex_[index]);
    };
    for (const size_t index : indices) {
      truthMatchMultiplicity_[index] = static_cast<int>(indices.size());
      if (score(index) < score(bestIndex))
        bestIndex = index;
    }
    truthMatchIsBest_[bestIndex] = 1;
  }

  eventRun_.push_back(event.id().run());
  eventLumi_.push_back(event.luminosityBlock());
  eventNumber_.push_back(event.id().event());
  tree_->Fill();
}

void ChiCGenTruthNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription d;
  d.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
  d.add<edm::InputTag>("oniaStates", edm::InputTag("chiCSharedOniaStates"));
  d.add<edm::InputTag>("conventionalCandidates", edm::InputTag("chiCConventionalPhotonStates"));
  d.add<edm::InputTag>("oneProngCandidates", edm::InputTag("chiCSharedOniaOneProngCandidates"));
  d.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  d.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  d.add<edm::InputTag>("packedToGenAssociation", edm::InputTag("packedPFCandidateToGenAssociation"));
  d.add<edm::InputTag>("lostToGenAssociation", edm::InputTag("lostTracksToGenAssociation"));
  d.add<bool>("doOneProng", true);
  descriptions.addWithDefaultLabel(d);
}

DEFINE_FWK_MODULE(ChiCGenTruthNtuplizer);
