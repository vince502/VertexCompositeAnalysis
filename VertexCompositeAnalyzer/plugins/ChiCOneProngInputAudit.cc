#include <cstdint>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class ChiCOneProngInputAudit : public edm::one::EDAnalyzer<> {
public:
  explicit ChiCOneProngInputAudit(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

private:
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  std::uint64_t events_ = 0;
  std::uint64_t chargedPacked_ = 0;
  std::uint64_t packedWithTrackDetails_ = 0;
  std::uint64_t chargedLost_ = 0;
  std::uint64_t lostWithTrackDetails_ = 0;
  double lowestTrackPt_ = -1.;
};

ChiCOneProngInputAudit::ChiCOneProngInputAudit(const edm::ParameterSet& cfg)
    : muonsToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
      packedToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      beamSpotToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))) {}

void ChiCOneProngInputAudit::analyze(const edm::Event& event, const edm::EventSetup&) {
  ++events_;
  const auto& packed = event.get(packedToken_);
  const auto& lost = event.get(lostToken_);
  // Force retrieval of the remaining mandatory MiniAOD inputs on every event.
  (void)event.get(muonsToken_);
  (void)event.get(verticesToken_);
  (void)event.get(beamSpotToken_);
  const auto count = [this](const pat::PackedCandidateCollection& collection,
                            std::uint64_t& charged,
                            std::uint64_t& detailed) {
    for (const auto& candidate : collection) {
      if (candidate.charge() == 0)
        continue;
      ++charged;
      if (!candidate.hasTrackDetails())
        continue;
      ++detailed;
      if (lowestTrackPt_ < 0. || candidate.pt() < lowestTrackPt_)
        lowestTrackPt_ = candidate.pt();
    }
  };
  count(packed, chargedPacked_, packedWithTrackDetails_);
  count(lost, chargedLost_, lostWithTrackDetails_);
}

void ChiCOneProngInputAudit::endJob() {
  edm::LogVerbatim("ChiCOneProngInputAudit")
      << "MiniAOD input audit: events=" << events_ << " chargedPacked=" << chargedPacked_
      << " packedWithTrackDetails=" << packedWithTrackDetails_ << " chargedLost=" << chargedLost_
      << " lostWithTrackDetails=" << lostWithTrackDetails_ << " lowestTrackPt=" << lowestTrackPt_
      << ". The presence of every configured muon, vertex, and beam-spot product was required event-by-event.";
}

void ChiCOneProngInputAudit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  description.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  description.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCOneProngInputAudit);
