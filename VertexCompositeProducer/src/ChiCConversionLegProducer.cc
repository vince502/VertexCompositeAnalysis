#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace {
  constexpr double kElectronMass = 0.00051099895;

  template <typename T, typename = void>
  struct HasPseudoPosDefTrack : std::false_type {};
  template <typename T>
  struct HasPseudoPosDefTrack<T, std::void_t<decltype(std::declval<const T&>().pseudoPosDefTrack())>> : std::true_type {};

  template <typename T>
  reco::Track extractPackedPseudoTrack(const T& cand) {
    if constexpr (HasPseudoPosDefTrack<T>::value) {
      return cand.pseudoPosDefTrack();
    } else {
      return cand.pseudoTrack();
    }
  }

  reco::Candidate::LorentzVector electronP4(const pat::PackedCandidate& candidate) {
    const double momentum2 =
        candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz();
    return {candidate.px(), candidate.py(), candidate.pz(), std::sqrt(momentum2 + kElectronMass * kElectronMass)};
  }

  void saveTrackCovariance(const reco::Track& track, pat::CompositeCandidate& leg) {
    // The MiniAOD packed covariance is a covariance of the five perigee
    // parameters. Persist the complete lower triangle; neither this producer
    // nor the ntuplizer is allowed to replace it with diagonal errors.
    for (unsigned int row = 0; row != 5; ++row) {
      for (unsigned int column = 0; column <= row; ++column) {
        leg.addUserFloat("trackCov" + std::to_string(row) + std::to_string(column), track.covariance(row, column));
      }
    }
  }
}  // namespace

class ChiCConversionLegProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCConversionLegProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  void append(const edm::Handle<pat::PackedCandidateCollection>&,
              unsigned int,
              const reco::Vertex::Point&,
              pat::CompositeCandidateCollection&) const;

  edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  double minPt_;
  double maxAbsEta_;
  bool requireHighPurity_;
};

ChiCConversionLegProducer::ChiCConversionLegProducer(const edm::ParameterSet& cfg)
    : packedToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      minPt_(cfg.getParameter<double>("minPt")),
      maxAbsEta_(cfg.getParameter<double>("maxAbsEta")),
      requireHighPurity_(cfg.getParameter<bool>("requireHighPurity")) {
  produces<pat::CompositeCandidateCollection>();
}

void ChiCConversionLegProducer::append(const edm::Handle<pat::PackedCandidateCollection>& inputHandle,
                                       unsigned int sourceCollection,
                                       const reco::Vertex::Point& pv,
                                       pat::CompositeCandidateCollection& output) const {
  if (!inputHandle.isValid())
    return;
  const auto& input = *inputHandle;
  for (size_t index = 0; index < input.size(); ++index) {
    const auto& track = input[index];
    if (track.charge() == 0 || !track.hasTrackDetails() || track.pt() < minPt_ || std::abs(track.eta()) > maxAbsEta_ ||
        (requireHighPurity_ && !track.trackHighPurity()))
      continue;
    const auto pseudoTrack = extractPackedPseudoTrack(track);
    pat::CompositeCandidate leg;
    leg.setP4(electronP4(track));
    leg.setCharge(track.charge());
    leg.setVertex(track.vertex());
    leg.addDaughter(track, "source");
    leg.addUserInt("sourceCollection", sourceCollection);
    leg.addUserInt("sourceKey", static_cast<int>(index));
    leg.addUserInt("sourcePdgId", track.pdgId());
    leg.addUserInt("trackHighPurity", track.trackHighPurity());
    leg.addUserInt("nTrackerHits", track.numberOfHits());
    leg.addUserInt("nPixelHits", track.numberOfPixelHits());
    leg.addUserInt("nTrackerLayers", track.trackerLayersWithMeasurement());
    leg.addUserInt("nPixelLayers", track.pixelLayersWithMeasurement());
    leg.addUserInt("lostInnerHits", static_cast<int>(track.lostInnerHits()));
    leg.addUserInt("firstHit", static_cast<int>(track.firstHit()));
    // PackedCandidate keeps the packed track parameters and covariance, but
    // not the original reco::Track algorithm label. pseudoTrack() therefore
    // carries the undefined enum value (zero); do not present it as physics.
    leg.addUserInt("trackAlgo", -1);
    leg.addUserInt("trackOriginalAlgo", -1);
    leg.addUserInt("trackAlgoAvailable", 0);
    leg.addUserFloat("trackNormalizedChi2", pseudoTrack.normalizedChi2());
    leg.addUserFloat("trackDxy", track.dxy(pv));
    leg.addUserFloat("trackDz", track.dz(pv));
    leg.addUserFloat("trackDxyError", track.dxyError());
    leg.addUserFloat("trackDzError", track.dzError());
    for (unsigned int parameter = 0; parameter != 5; ++parameter) {
      leg.addUserFloat("trackParameter" + std::to_string(parameter), pseudoTrack.parameter(parameter));
    }
    saveTrackCovariance(pseudoTrack, leg);
    output.push_back(leg);
  }
}

void ChiCConversionLegProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  edm::Handle<pat::PackedCandidateCollection> packed;
  edm::Handle<pat::PackedCandidateCollection> lost;
  event.getByToken(packedToken_, packed);
  event.getByToken(lostToken_, lost);
  const auto& vertices = event.get(verticesToken_);
  const reco::Vertex::Point pv = vertices.empty() ? reco::Vertex::Point(0., 0., 0.) : vertices.front().position();
  auto output = std::make_unique<pat::CompositeCandidateCollection>();
  append(packed, 0, pv, *output);
  append(lost, 1, pv, *output);
  event.put(std::move(output));
}

void ChiCConversionLegProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  description.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<double>("minPt", 0.2);
  description.add<double>("maxAbsEta", 2.5);
  description.add<bool>("requireHighPurity", false);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCConversionLegProducer);
