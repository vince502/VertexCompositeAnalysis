#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class VertexToPatCompositeCandidateConverter : public edm::stream::EDProducer<> {
public:
  explicit VertexToPatCompositeCandidateConverter(const edm::ParameterSet&);
  ~VertexToPatCompositeCandidateConverter() override = default;

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> vertexToken_;
};

VertexToPatCompositeCandidateConverter::VertexToPatCompositeCandidateConverter(const edm::ParameterSet& iConfig) {
  vertexToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("src"));
  produces<pat::CompositeCandidateCollection>();
}

void VertexToPatCompositeCandidateConverter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto output = std::make_unique<pat::CompositeCandidateCollection>();

  edm::Handle<reco::VertexCompositeCandidateCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  if (!vertices.isValid()) {
    // edm::LogWarning("VertexToPatCompositeCandidateConverter") << "No reco::VertexCompositeCandidate collection found!";
    return;
  }

  for (const auto& vtxCand : *vertices) {
    pat::CompositeCandidate patCand(vtxCand); // Copy base Candidate members

    // Add extra VertexCompositeCandidate-specific members as userData
    patCand.addUserFloat("vertexChi2", vtxCand.vertexChi2());
    patCand.addUserFloat("vertexNdof", vtxCand.vertexNdof());
    patCand.addUserFloat("vProb", vtxCand.vertexNormalizedChi2());
    patCand.addUserFloat("x", vtxCand.vertex().x());
    patCand.addUserFloat("y", vtxCand.vertex().y());
    patCand.addUserFloat("z", vtxCand.vertex().z());
    patCand.addUserFloat("xError", std::sqrt(vtxCand.vertexCovariance(0, 0)));
    patCand.addUserFloat("yError", std::sqrt(vtxCand.vertexCovariance(1, 1)));
    patCand.addUserFloat("zError", std::sqrt(vtxCand.vertexCovariance(2, 2)));

    output->emplace_back(std::move(patCand));
  }

  iEvent.put(std::move(output));
}

// Define as a plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexToPatCompositeCandidateConverter);