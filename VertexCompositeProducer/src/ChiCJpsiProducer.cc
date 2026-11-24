#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TMath.h"

class ChiCJpsiProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCJpsiProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
  double minMuonPt_;
  double maxMuonAbsEta_;
  double minMass_;
  double maxMass_;
  double minVertexProbability_;
};

ChiCJpsiProducer::ChiCJpsiProducer(const edm::ParameterSet& cfg)
    : muonsToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      beamSpotToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      minMuonPt_(cfg.getParameter<double>("minMuonPt")),
      maxMuonAbsEta_(cfg.getParameter<double>("maxMuonAbsEta")),
      minMass_(cfg.getParameter<double>("minMass")),
      maxMass_(cfg.getParameter<double>("maxMass")),
      minVertexProbability_(cfg.getParameter<double>("minVertexProbability")) {
  produces<pat::CompositeCandidateCollection>();
}

void ChiCJpsiProducer::produce(edm::StreamID,
                               edm::Event& event,
                               const edm::EventSetup& setup) const {
  edm::Handle<pat::MuonCollection> muonsHandle;
  event.getByToken(muonsToken_, muonsHandle);
  const auto& muons = *muonsHandle;
  const auto& vertices = event.get(verticesToken_);
  const auto& beamSpot = event.get(beamSpotToken_);
  const auto& transientTrackBuilder = setup.getData(transientTrackBuilderToken_);
  auto output = std::make_unique<pat::CompositeCandidateCollection>();

  for (size_t first = 0; first < muons.size(); ++first) {
    const auto& muon1 = muons[first];
    if (muon1.pt() < minMuonPt_ || std::abs(muon1.eta()) > maxMuonAbsEta_ || muon1.bestTrack() == nullptr)
      continue;
    for (size_t second = first + 1; second < muons.size(); ++second) {
      const auto& muon2 = muons[second];
      if (muon2.pt() < minMuonPt_ || std::abs(muon2.eta()) > maxMuonAbsEta_ || muon2.bestTrack() == nullptr ||
          muon1.charge() * muon2.charge() >= 0)
        continue;

      const auto p4 = muon1.p4() + muon2.p4();
      if (p4.mass() < minMass_ || p4.mass() > maxMass_)
        continue;
      std::vector<reco::TransientTrack> tracks;
      tracks.emplace_back(transientTrackBuilder.build(muon1.bestTrack()));
      tracks.emplace_back(transientTrackBuilder.build(muon2.bestTrack()));
      KalmanVertexFitter fitter(true);
      const TransientVertex transientVertex = fitter.vertex(tracks);
      if (!transientVertex.isValid())
        continue;
      const reco::Vertex vertex(transientVertex);
      const double probability = TMath::Prob(vertex.chi2(), vertex.ndof());
      if (probability < minVertexProbability_)
        continue;

      size_t pvIndex = 0;
      double bestDz = std::numeric_limits<double>::infinity();
      for (size_t index = 0; index < vertices.size(); ++index) {
        const double dz = std::abs(vertex.z() - vertices[index].z());
        if (dz < bestDz) {
          bestDz = dz;
          pvIndex = index;
        }
      }
      const auto origin = vertices.empty() ? beamSpot.position() : vertices[pvIndex].position();

      pat::CompositeCandidate candidate;
      candidate.setP4(p4);
      candidate.setCharge(0);
      candidate.setVertex(vertex.position());
      candidate.addDaughter(muon1, "muon1");
      candidate.addDaughter(muon2, "muon2");
      candidate.addUserInt("muon1Index", static_cast<int>(first));
      candidate.addUserInt("muon2Index", static_cast<int>(second));
      candidate.addUserFloat("unconstrainedMass", p4.mass());
      candidate.addUserFloat("vertexChi2", vertex.chi2());
      candidate.addUserFloat("vertexNdof", vertex.ndof());
      candidate.addUserFloat("vertexProbability", probability);
      candidate.addUserFloat("vertexXError", vertex.xError());
      candidate.addUserFloat("vertexYError", vertex.yError());
      candidate.addUserFloat("vertexZError", vertex.zError());
      candidate.addUserFloat("vertexCovXY", vertex.covariance(0, 1));
      candidate.addUserFloat("vertexCovXZ", vertex.covariance(0, 2));
      candidate.addUserFloat("vertexCovYZ", vertex.covariance(1, 2));
      candidate.addUserFloat("pvDz", vertex.z() - origin.z());
      candidate.addUserInt("primaryVertexIndex", static_cast<int>(pvIndex));
      candidate.addUserInt("jpsiStatus", vertices.empty() ? 0x2 : 0x1);
      output->push_back(candidate);
    }
  }
  event.put(std::move(output));
}

void ChiCJpsiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  description.add<double>("minMuonPt", 2.0);
  description.add<double>("maxMuonAbsEta", 2.4);
  description.add<double>("minMass", 2.7);
  description.add<double>("maxMass", 3.5);
  description.add<double>("minVertexProbability", 0.005);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCJpsiProducer);
