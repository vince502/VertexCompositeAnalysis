#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/interface/ChiCEventActivity.h"
#include "TTree.h"

// Aligned companion to the established covariance and EventInfo trees.
// Stores cut inputs rather than applying offline ID, acceptance, trigger,
// vertex-probability, conversion-radius or mass-difference requirements.
class OniaPhotonOfflineNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit OniaPhotonOfflineNtuplizer(const edm::ParameterSet& cfg)
      : onias_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("onias"))),
        chis_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("chis"))),
        muons_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        vertices_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
        packed_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
        lost_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
        activity_{cfg.getParameter<double>("trackMinPt"),
                  cfg.getParameter<double>("trackMaxAbsEta"),
                  cfg.getParameter<double>("trackMaxRelPtError"),
                  cfg.getParameter<double>("trackMaxDzSignificance"),
                  cfg.getParameter<double>("trackMaxDxySignificance")},
        triggers_(mayConsume<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("triggerObjects"))),
        results_(mayConsume<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"))) {
    usesResource(TFileService::kSharedResource);
    edm::Service<TFileService> file;
    tree_ = file->make<TTree>("Offline", "Post-production cut inputs, aligned with ChiCCovariance and EventInfo");
    tree_->Branch("run", &run_);
    tree_->Branch("lumi", &lumi_);
    tree_->Branch("event", &event_);
    tree_->Branch("triggerStatus", &triggerStatus_);
    tree_->Branch("nSourceMuons", &nSourceMuons_);
    for (const std::string name : {
             "oniaCharge",        "oniaInputIndex",       "oniaPVIndex",       "oniaPVNTrack",
             "oniaMuonBegin",     "oniaMuonSize",         "muonOniaIndex",     "muonSourceIndex",
             "muonSourceMatches", "muonIsGlobal",         "muonIsTracker",     "muonIsPF",
             "muonTMOneStationTight", "muonSoftId",       "muonTightId",       "muonLooseId",
             "muonHighPurity",    "muonTrackerLayers",    "muonPixelLayers",   "muonTrackerHits",
             "muonPixelHits",     "muonGlobalHits",       "muonStations",      "photonSourceIndex",
             "photonFlags",       "triggerPathBegin",     "triggerPathSize",   "triggerFilterBegin",
             "triggerFilterSize", "triggerIdBegin",       "triggerIdSize",     "triggerIds",
             "triggerPathLastFilter", "triggerPathL3Filter"})
      tree_->Branch(name.c_str(), &integers_[name]);
    for (const std::string name : {
             "oniaVertexProbability", "oniaMassError", "oniaDCA",         "muonPtError",
             "muonDxyPV",             "muonDzPV",      "muonDxyError",    "muonDzError",
             "muonTrackNormChi2",     "muonGlobalNormChi2", "photonInputPx", "photonInputPy",
             "photonInputPz",         "photonInputE",  "photonFitChi2",   "photonFitNdof",
             "photonOniaDz",          "triggerPt",     "triggerEta",      "triggerPhi",
             "triggerMass"})
      tree_->Branch(name.c_str(), &floats_[name]);
    tree_->Branch("triggerPaths", &triggerPaths_);
    tree_->Branch("triggerFilters", &triggerFilters_);
    tree_->Branch("triggerCollections", &triggerCollections_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription d;
    d.add<edm::InputTag>("onias", edm::InputTag("chiCSharedOniaStates"));
    d.add<edm::InputTag>("chis", edm::InputTag("chiCConventionalPhotonStates"));
    d.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
    d.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    d.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
    d.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
    d.add<double>("trackMinPt", 0.4);
    d.add<double>("trackMaxAbsEta", 2.4);
    d.add<double>("trackMaxRelPtError", 0.1);
    d.add<double>("trackMaxDzSignificance", 3.0);
    d.add<double>("trackMaxDxySignificance", 3.0);
    d.add<edm::InputTag>("triggerObjects", edm::InputTag("slimmedPatTrigger"));
    d.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
    descriptions.addWithDefaultLabel(d);
  }

  void analyze(const edm::Event& event, const edm::EventSetup&) override {
    for (auto& field : integers_)
      field.second.clear();
    for (auto& field : floats_)
      field.second.clear();
    triggerPaths_.clear();
    triggerFilters_.clear();
    triggerCollections_.clear();
    run_ = event.id().run();
    lumi_ = event.id().luminosityBlock();
    event_ = event.id().event();

    const auto& onias = event.get(onias_);
    const auto& chis = event.get(chis_);
    const auto& muons = event.get(muons_);
    nSourceMuons_ = muons.size();
    const auto& vertices = event.get(vertices_);

    const auto packedHandle = event.getHandle(packed_);
    const auto lostHandle = event.getHandle(lost_);
    const bool hasTrackActivity = packedHandle.isValid() && lostHandle.isValid();
    const std::vector<unsigned int> ntracks =
        hasTrackActivity ? chic::countAssociatedTracks(vertices, *packedHandle, *lostHandle, activity_)
                         : std::vector<unsigned int>{};

    for (size_t index = 0; index < onias.size(); ++index) {
      const auto& onia = onias[index];
      const int pv = onia.hasUserInt("pvIndex") ? onia.userInt("pvIndex") : -1;
      const auto* vertex = (pv >= 0 && static_cast<size_t>(pv) < vertices.size()) ? &vertices[pv] : nullptr;
      const int inputIdx = onia.hasUserInt("oniaInputIndex")
                               ? onia.userInt("oniaInputIndex")
                               : (onia.hasUserInt("oniaStateIndex") ? onia.userInt("oniaStateIndex")
                                                                    : static_cast<int>(index));
      integers_["oniaCharge"].push_back(onia.charge());
      integers_["oniaInputIndex"].push_back(inputIdx);
      integers_["oniaPVIndex"].push_back(pv);
      integers_["oniaPVNTrack"].push_back(
          (vertex && hasTrackActivity && static_cast<size_t>(pv) < ntracks.size()) ? static_cast<int>(ntracks[pv]) : -1);
      integers_["oniaMuonBegin"].push_back(integers_["muonOniaIndex"].size());
      floats_["oniaVertexProbability"].push_back(value(onia, "vProb"));
      floats_["oniaMassError"].push_back(value(onia, "MassErr"));
      floats_["oniaDCA"].push_back(value(onia, "DCA"));
      int size = 0;
      for (const char* label : {"muon1", "muon2"}) {
        const auto* muon = dynamic_cast<const pat::Muon*>(onia.daughter(label));
        if (!muon)
          continue;
        ++size;
        const int source =
            muon->hasUserInt("sourceSlimmedMuonIndex") ? muon->userInt("sourceSlimmedMuonIndex") : -1;
        const bool sourceValid = source >= 0 && static_cast<size_t>(source) < muons.size();
        integers_["muonSourceIndex"].push_back(sourceValid ? source : -1);
        integers_["muonSourceMatches"].push_back(sourceValid ? 1 : 0);
        integers_["muonOniaIndex"].push_back(static_cast<int>(index));
        integers_["muonIsGlobal"].push_back(muon->isGlobalMuon());
        integers_["muonIsTracker"].push_back(muon->isTrackerMuon());
        integers_["muonIsPF"].push_back(muon->isPFMuon());
        integers_["muonTMOneStationTight"].push_back(muon->muonID("TMOneStationTight"));
        integers_["muonLooseId"].push_back(muon::isLooseMuon(*muon));
        integers_["muonSoftId"].push_back(vertex ? muon::isSoftMuon(*muon, *vertex) : -1);
        integers_["muonTightId"].push_back(vertex ? muon::isTightMuon(*muon, *vertex) : -1);
        integers_["muonStations"].push_back(muon->numberOfMatchedStations());
        const auto track = muon->innerTrack();
        integers_["muonHighPurity"].push_back(track.isNonnull() ? track->quality(reco::TrackBase::highPurity) : -1);
        integers_["muonTrackerLayers"].push_back(
            track.isNonnull() ? track->hitPattern().trackerLayersWithMeasurement() : -1);
        integers_["muonPixelLayers"].push_back(
            track.isNonnull() ? track->hitPattern().pixelLayersWithMeasurement() : -1);
        integers_["muonTrackerHits"].push_back(track.isNonnull() ? track->numberOfValidHits() : -1);
        integers_["muonPixelHits"].push_back(track.isNonnull() ? track->hitPattern().numberOfValidPixelHits() : -1);
        floats_["muonPtError"].push_back(track.isNonnull() ? track->ptError() : missing());
        floats_["muonDxyPV"].push_back(track.isNonnull() && vertex ? track->dxy(vertex->position()) : missing());
        floats_["muonDzPV"].push_back(track.isNonnull() && vertex ? track->dz(vertex->position()) : missing());
        floats_["muonDxyError"].push_back(track.isNonnull() ? track->dxyError() : missing());
        floats_["muonDzError"].push_back(track.isNonnull() ? track->dzError() : missing());
        floats_["muonTrackNormChi2"].push_back(track.isNonnull() ? track->normalizedChi2() : missing());
        const auto global = muon->globalTrack();
        integers_["muonGlobalHits"].push_back(global.isNonnull() ? global->hitPattern().numberOfValidMuonHits() : -1);
        floats_["muonGlobalNormChi2"].push_back(global.isNonnull() ? global->normalizedChi2() : missing());
      }
      integers_["oniaMuonSize"].push_back(size);
    }

    for (const auto& chi : chis) {
      const auto* photon = dynamic_cast<const pat::CompositeCandidate*>(chi.daughter("photon"));
      const int phoSrcIdx = chi.hasUserInt("photonSourceIndex")
                                ? chi.userInt("photonSourceIndex")
                                : (chi.hasUserInt("convIdx") ? chi.userInt("convIdx") : -1);
      integers_["photonSourceIndex"].push_back(phoSrcIdx);
      integers_["photonFlags"].push_back(photon && photon->hasUserInt("flags") ? photon->userInt("flags") : -1);
      floats_["photonInputPx"].push_back(photon ? photon->px() : missing());
      floats_["photonInputPy"].push_back(photon ? photon->py() : missing());
      floats_["photonInputPz"].push_back(photon ? photon->pz() : missing());
      floats_["photonInputE"].push_back(photon ? photon->energy() : missing());
      floats_["photonFitChi2"].push_back(value(chi, "conventionalPhotonFitChi2"));
      floats_["photonFitNdof"].push_back(value(chi, "conventionalPhotonFitNdof"));
      floats_["photonOniaDz"].push_back(value(chi, "dz"));
    }

    const auto results = event.getHandle(results_);
    const auto triggers = event.getHandle(triggers_);
    triggerStatus_ = (results.isValid() ? 1 : 0) | (triggers.isValid() ? 2 : 0);
    if (results.isValid() && triggers.isValid()) {
      const auto& names = event.triggerNames(*results);
      for (auto object : *triggers) {
        object.unpackPathNames(names);
        object.unpackFilterLabels(event, *results);
        floats_["triggerPt"].push_back(object.pt());
        floats_["triggerEta"].push_back(object.eta());
        floats_["triggerPhi"].push_back(object.phi());
        floats_["triggerMass"].push_back(object.mass());
        triggerCollections_.push_back(object.collection());
        integers_["triggerPathBegin"].push_back(triggerPaths_.size());
        const auto paths = object.pathNames(false, false);
        integers_["triggerPathSize"].push_back(paths.size());
        for (const auto& path : paths) {
          triggerPaths_.push_back(path);
          integers_["triggerPathLastFilter"].push_back(object.hasPathName(path, true, false));
          integers_["triggerPathL3Filter"].push_back(object.hasPathName(path, false, true));
        }
        integers_["triggerFilterBegin"].push_back(triggerFilters_.size());
        integers_["triggerFilterSize"].push_back(object.filterLabels().size());
        for (const auto& filter : object.filterLabels())
          triggerFilters_.push_back(filter);
        integers_["triggerIdBegin"].push_back(integers_["triggerIds"].size());
        integers_["triggerIdSize"].push_back(object.filterIds().size());
        for (int id : object.filterIds())
          integers_["triggerIds"].push_back(id);
      }
    }
    tree_->Fill();
  }

private:
  static float missing() { return std::numeric_limits<float>::quiet_NaN(); }
  static float value(const pat::CompositeCandidate& c, const char* key) {
    return c.hasUserFloat(key) ? c.userFloat(key) : missing();
  }

  edm::EDGetTokenT<pat::CompositeCandidateCollection> onias_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> chis_;
  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packed_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lost_;
  chic::TrackActivityConfig activity_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggers_;
  edm::EDGetTokenT<edm::TriggerResults> results_;
  unsigned int run_ = 0, lumi_ = 0;
  unsigned int nSourceMuons_ = 0;
  unsigned long long event_ = 0;
  int triggerStatus_ = 0;
  std::map<std::string, std::vector<int>> integers_;
  std::map<std::string, std::vector<float>> floats_;
  std::vector<std::string> triggerPaths_, triggerFilters_, triggerCollections_;
  TTree* tree_ = nullptr;
};

DEFINE_FWK_MODULE(OniaPhotonOfflineNtuplizer);
