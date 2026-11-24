#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/interface/ChiCEventActivity.h"
#include "TTree.h"

namespace {
  constexpr float kMissing = std::numeric_limits<float>::quiet_NaN();

  std::vector<unsigned int> selectedIndices(const std::vector<std::string>& names,
                                            const std::vector<std::string>& patterns) {
    std::set<unsigned int> selected;
    for (const auto& pattern : patterns) {
      for (const auto match : edm::regexMatch(names, pattern))
        selected.insert(std::distance(names.begin(), match));
    }
    return std::vector<unsigned int>(selected.begin(), selected.end());
  }

  std::string parameterSetID(const edm::TriggerResults& results) {
    std::string value;
    results.parameterSetID().toString(value);
    return value;
  }
}  // namespace

// Compact all-events companion to ChiCCovarianceNtuplizer. Trigger names are
// written once per distinct HLT/RECO menu; each event stores only byte-sized
// decision vectors and a menu index.
class ChiCEventInfoNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCEventInfoNtuplizer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  struct MenuCacheEntry {
    int index = -1;
    std::vector<unsigned int> hltIndices;
    std::vector<unsigned int> filterIndices;
  };

  MenuCacheEntry const& menu(const edm::Event&,
                             const edm::Handle<edm::TriggerResults>&,
                             const edm::Handle<edm::TriggerResults>&);
  void clear();

  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandidatesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<HFRecHitCollection> hfRecHitsToken_;
  edm::EDGetTokenT<reco::Centrality> centralityToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupToken_;
  edm::EDGetTokenT<edm::TriggerResults> hltResultsToken_;
  edm::EDGetTokenT<edm::TriggerResults> eventFilterResultsToken_;
  chic::TrackActivityConfig trackActivityConfig_;
  std::vector<std::string> hltPathPatterns_;
  std::vector<std::string> eventFilterPatterns_;
  bool isMC_;

  TTree* eventTree_ = nullptr;
  TTree* menuTree_ = nullptr;
  std::map<std::string, MenuCacheEntry> menuCache_;

  unsigned int run_ = 0, lumi_ = 0;
  unsigned long long event_ = 0;
  bool isRealData_ = false;
  unsigned int nPV_ = 0;
  int trackStatus_ = 0, beamSpotStatus_ = 0, rhoStatus_ = 0, hfStatus_ = 0;
  int centralityStatus_ = 0, pileupStatus_ = 0, hltStatus_ = 0, eventFilterStatus_ = 0;
  std::vector<float> pvX_, pvY_, pvZ_, pvNdof_;
  std::vector<unsigned char> pvIsFake_;
  std::vector<unsigned int> pvNTrack_;
  float beamSpotX_ = kMissing, beamSpotY_ = kMissing, beamSpotZ_ = kMissing, beamSpotSigmaZ_ = kMissing;
  float beamSpotWidthX_ = kMissing, beamSpotWidthY_ = kMissing, beamSpotDxdz_ = kMissing, beamSpotDydz_ = kMissing;
  float rho_ = kMissing;
  float hfRecHitEnergyPlus_ = kMissing, hfRecHitEnergyMinus_ = kMissing;
  unsigned int hfRecHitCountPlus_ = 0, hfRecHitCountMinus_ = 0;
  float centralityHFHitSumPlus_ = kMissing, centralityHFHitSumMinus_ = kMissing;
  float centralityHFTowerSumPlus_ = kMissing, centralityHFTowerSumMinus_ = kMissing;
  float centralityNTrack_ = kMissing;
  int pileupObserved_ = -1;
  float pileupTrue_ = kMissing;
  int menuIndex_ = -1;
  std::vector<unsigned char> hltAccept_, eventFilterAccept_;

  int menuTreeIndex_ = -1;
  std::string menuHLTParameterSet_, menuEventFilterParameterSet_;
  std::vector<std::string> menuHLTNames_, menuEventFilterNames_;
};

ChiCEventInfoNtuplizer::ChiCEventInfoNtuplizer(const edm::ParameterSet& cfg)
    : primaryVerticesToken_(mayConsume<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      packedCandidatesToken_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostTracksToken_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      beamSpotToken_(mayConsume<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
      rhoToken_(mayConsume<double>(cfg.getParameter<edm::InputTag>("rho"))),
      hfRecHitsToken_(mayConsume<HFRecHitCollection>(cfg.getParameter<edm::InputTag>("hfRecHits"))),
      centralityToken_(mayConsume<reco::Centrality>(cfg.getParameter<edm::InputTag>("centrality"))),
      pileupToken_(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("pileupInfo"))),
      hltResultsToken_(mayConsume<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("hltResults"))),
      eventFilterResultsToken_(mayConsume<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("eventFilterResults"))),
      trackActivityConfig_{cfg.getParameter<double>("trackMinPt"),
                           cfg.getParameter<double>("trackMaxAbsEta"),
                           cfg.getParameter<double>("trackMaxRelPtError"),
                           cfg.getParameter<double>("trackMaxDzSignificance"),
                           cfg.getParameter<double>("trackMaxDxySignificance")},
      hltPathPatterns_(cfg.getParameter<std::vector<std::string>>("hltPathPatterns")),
      eventFilterPatterns_(cfg.getParameter<std::vector<std::string>>("eventFilterPatterns")),
      isMC_(cfg.getParameter<bool>("isMC")) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> file;
  eventTree_ = file->make<TTree>("EventInfo", "All-input-event information for chi_c analysis constraints");
#define EVENT_BRANCH(name) eventTree_->Branch(#name, &name)
  EVENT_BRANCH(run_); EVENT_BRANCH(lumi_); EVENT_BRANCH(event_); EVENT_BRANCH(isRealData_);
  EVENT_BRANCH(nPV_); EVENT_BRANCH(pvX_); EVENT_BRANCH(pvY_); EVENT_BRANCH(pvZ_); EVENT_BRANCH(pvNdof_); EVENT_BRANCH(pvIsFake_); EVENT_BRANCH(pvNTrack_);
  EVENT_BRANCH(trackStatus_); EVENT_BRANCH(beamSpotStatus_); EVENT_BRANCH(rhoStatus_); EVENT_BRANCH(hfStatus_); EVENT_BRANCH(centralityStatus_); EVENT_BRANCH(pileupStatus_); EVENT_BRANCH(hltStatus_); EVENT_BRANCH(eventFilterStatus_);
  EVENT_BRANCH(beamSpotX_); EVENT_BRANCH(beamSpotY_); EVENT_BRANCH(beamSpotZ_); EVENT_BRANCH(beamSpotSigmaZ_); EVENT_BRANCH(beamSpotWidthX_); EVENT_BRANCH(beamSpotWidthY_); EVENT_BRANCH(beamSpotDxdz_); EVENT_BRANCH(beamSpotDydz_);
  EVENT_BRANCH(rho_); EVENT_BRANCH(hfRecHitEnergyPlus_); EVENT_BRANCH(hfRecHitEnergyMinus_); EVENT_BRANCH(hfRecHitCountPlus_); EVENT_BRANCH(hfRecHitCountMinus_);
  EVENT_BRANCH(centralityHFHitSumPlus_); EVENT_BRANCH(centralityHFHitSumMinus_); EVENT_BRANCH(centralityHFTowerSumPlus_); EVENT_BRANCH(centralityHFTowerSumMinus_); EVENT_BRANCH(centralityNTrack_);
  EVENT_BRANCH(pileupObserved_); EVENT_BRANCH(pileupTrue_); EVENT_BRANCH(menuIndex_); EVENT_BRANCH(hltAccept_); EVENT_BRANCH(eventFilterAccept_);
#undef EVENT_BRANCH

  menuTree_ = file->make<TTree>("TriggerMenus", "Names corresponding to EventInfo trigger/filter decision vectors");
#define MENU_BRANCH(name) menuTree_->Branch(#name, &name)
  MENU_BRANCH(menuTreeIndex_); MENU_BRANCH(menuHLTParameterSet_); MENU_BRANCH(menuEventFilterParameterSet_);
  MENU_BRANCH(menuHLTNames_); MENU_BRANCH(menuEventFilterNames_);
#undef MENU_BRANCH
}

void ChiCEventInfoNtuplizer::clear() {
  nPV_ = 0;
  trackStatus_ = beamSpotStatus_ = rhoStatus_ = hfStatus_ = 0;
  centralityStatus_ = pileupStatus_ = hltStatus_ = eventFilterStatus_ = 0;
  pvX_.clear(); pvY_.clear(); pvZ_.clear(); pvNdof_.clear(); pvIsFake_.clear(); pvNTrack_.clear();
  beamSpotX_ = beamSpotY_ = beamSpotZ_ = beamSpotSigmaZ_ = kMissing;
  beamSpotWidthX_ = beamSpotWidthY_ = beamSpotDxdz_ = beamSpotDydz_ = kMissing;
  rho_ = hfRecHitEnergyPlus_ = hfRecHitEnergyMinus_ = kMissing;
  hfRecHitCountPlus_ = hfRecHitCountMinus_ = 0;
  centralityHFHitSumPlus_ = centralityHFHitSumMinus_ = kMissing;
  centralityHFTowerSumPlus_ = centralityHFTowerSumMinus_ = centralityNTrack_ = kMissing;
  pileupObserved_ = -1; pileupTrue_ = kMissing;
  menuIndex_ = -1; hltAccept_.clear(); eventFilterAccept_.clear();
}

ChiCEventInfoNtuplizer::MenuCacheEntry const& ChiCEventInfoNtuplizer::menu(
    const edm::Event& event,
    const edm::Handle<edm::TriggerResults>& hlt,
    const edm::Handle<edm::TriggerResults>& filters) {
  const std::string hltID = hlt.isValid() ? parameterSetID(*hlt) : "missing";
  const std::string filterID = filters.isValid() ? parameterSetID(*filters) : "missing";
  const std::string key = hltID + ":" + filterID;
  auto found = menuCache_.find(key);
  if (found != menuCache_.end())
    return found->second;

  MenuCacheEntry entry;
  entry.index = static_cast<int>(menuCache_.size());
  menuTreeIndex_ = entry.index;
  menuHLTParameterSet_ = hltID;
  menuEventFilterParameterSet_ = filterID;
  menuHLTNames_.clear();
  menuEventFilterNames_.clear();
  if (hlt.isValid()) {
    const auto& names = event.triggerNames(*hlt).triggerNames();
    entry.hltIndices = selectedIndices(names, hltPathPatterns_);
    for (const auto index : entry.hltIndices)
      menuHLTNames_.push_back(names[index]);
  }
  if (filters.isValid()) {
    const auto& names = event.triggerNames(*filters).triggerNames();
    entry.filterIndices = selectedIndices(names, eventFilterPatterns_);
    for (const auto index : entry.filterIndices)
      menuEventFilterNames_.push_back(names[index]);
  }
  menuTree_->Fill();
  return menuCache_.emplace(key, std::move(entry)).first->second;
}

void ChiCEventInfoNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  clear();
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();
  isRealData_ = event.isRealData();

  const auto vertices = event.getHandle(primaryVerticesToken_);
  const auto packedCandidates = event.getHandle(packedCandidatesToken_);
  const auto lostTracks = event.getHandle(lostTracksToken_);
  if (vertices.isValid()) {
    nPV_ = vertices->size();
    for (const auto& vertex : *vertices) {
      pvX_.push_back(vertex.x());
      pvY_.push_back(vertex.y());
      pvZ_.push_back(vertex.z());
      pvNdof_.push_back(vertex.ndof());
      pvIsFake_.push_back(vertex.isFake());
    }
  }
  if (vertices.isValid() && packedCandidates.isValid() && lostTracks.isValid()) {
    trackStatus_ = 1;
    pvNTrack_ = chic::countAssociatedTracks(*vertices, *packedCandidates, *lostTracks, trackActivityConfig_);
  } else if (vertices.isValid()) {
    pvNTrack_.assign(vertices->size(), 0);
  }

  const auto beamSpot = event.getHandle(beamSpotToken_);
  if (beamSpot.isValid()) {
    beamSpotStatus_ = 1;
    beamSpotX_ = beamSpot->x0();
    beamSpotY_ = beamSpot->y0();
    beamSpotZ_ = beamSpot->z0();
    beamSpotSigmaZ_ = beamSpot->sigmaZ();
    beamSpotWidthX_ = beamSpot->BeamWidthX();
    beamSpotWidthY_ = beamSpot->BeamWidthY();
    beamSpotDxdz_ = beamSpot->dxdz();
    beamSpotDydz_ = beamSpot->dydz();
  }
  const auto rho = event.getHandle(rhoToken_);
  if (rho.isValid()) {
    rhoStatus_ = 1;
    rho_ = static_cast<float>(*rho);
  }

  const auto hfRecHits = event.getHandle(hfRecHitsToken_);
  if (hfRecHits.isValid()) {
    hfStatus_ = 1;
    hfRecHitEnergyPlus_ = hfRecHitEnergyMinus_ = 0.f;
    for (const auto& hit : *hfRecHits) {
      const HcalDetId id(hit.id());
      if (id.zside() > 0) {
        hfRecHitEnergyPlus_ += hit.energy();
        ++hfRecHitCountPlus_;
      } else {
        hfRecHitEnergyMinus_ += hit.energy();
        ++hfRecHitCountMinus_;
      }
    }
  }
  const auto centrality = event.getHandle(centralityToken_);
  if (centrality.isValid()) {
    centralityStatus_ = 1;
    centralityHFHitSumPlus_ = centrality->EtHFhitSumPlus();
    centralityHFHitSumMinus_ = centrality->EtHFhitSumMinus();
    centralityHFTowerSumPlus_ = centrality->EtHFtowerSumPlus();
    centralityHFTowerSumMinus_ = centrality->EtHFtowerSumMinus();
    centralityNTrack_ = centrality->Ntracks();
  }
  if (isMC_) {
    const auto pileup = event.getHandle(pileupToken_);
    if (pileup.isValid()) {
      for (const auto& bunchCrossing : *pileup) {
        if (bunchCrossing.getBunchCrossing() != 0)
          continue;
        pileupStatus_ = 1;
        pileupObserved_ = bunchCrossing.getPU_NumInteractions();
        pileupTrue_ = bunchCrossing.getTrueNumInteractions();
        break;
      }
    }
  }

  const auto hlt = event.getHandle(hltResultsToken_);
  const auto filters = event.getHandle(eventFilterResultsToken_);
  hltStatus_ = hlt.isValid();
  eventFilterStatus_ = filters.isValid();
  const auto& cachedMenu = menu(event, hlt, filters);
  menuIndex_ = cachedMenu.index;
  if (hlt.isValid())
    for (const auto index : cachedMenu.hltIndices)
      hltAccept_.push_back(hlt->accept(index));
  if (filters.isValid())
    for (const auto index : cachedMenu.filterIndices)
      eventFilterAccept_.push_back(filters->accept(index));

  eventTree_->Fill();
}

void ChiCEventInfoNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  description.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  description.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  description.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
  description.add<edm::InputTag>("hfRecHits", edm::InputTag("slimmedHcalRecHits", "reducedHcalRecHits"));
  description.add<edm::InputTag>("centrality", edm::InputTag("hiCentrality"));
  description.add<edm::InputTag>("pileupInfo", edm::InputTag("slimmedAddPileupInfo"));
  description.add<edm::InputTag>("hltResults", edm::InputTag("TriggerResults", "", "HLT"));
  description.add<edm::InputTag>("eventFilterResults", edm::InputTag("TriggerResults", "", "RECO"));
  description.add<double>("trackMinPt", 0.4);
  description.add<double>("trackMaxAbsEta", 2.4);
  description.add<double>("trackMaxRelPtError", 0.1);
  description.add<double>("trackMaxDzSignificance", 3.0);
  description.add<double>("trackMaxDxySignificance", 3.0);
  description.add<std::vector<std::string>>(
      "hltPathPatterns", {"HLT_PPRefL*Mu*", "Dataset_PPRefDoubleMuon*", "Dataset_PPRefSingleMuon*"});
  description.add<std::vector<std::string>>("eventFilterPatterns", {"Flag_*"});
  description.add<bool>("isMC", false);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCEventInfoNtuplizer);
