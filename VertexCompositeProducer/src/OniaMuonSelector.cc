#include <memory>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// Preserve the original slimmedMuons index before PAT copying/track embedding.
// The two production requirements are explicit; IDs and acceptance stay offline.
class OniaMuonSelector : public edm::global::EDProducer<> {
public:
  explicit OniaMuonSelector(const edm::ParameterSet& cfg)
      : source_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("src"))),
        minPt_(cfg.getParameter<double>("minPt")) {
    produces<pat::MuonCollection>();
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& input = event.get(source_);
    auto output = std::make_unique<pat::MuonCollection>();
    for (size_t index = 0; index < input.size(); ++index) {
      const auto& muon = input[index];
      if (muon.pt() <= minPt_ || muon.innerTrack().isNull())
        continue;
      pat::Muon selected(muon);
      selected.addUserInt("sourceSlimmedMuonIndex", static_cast<int>(index));
      output->push_back(selected);
    }
    event.put(std::move(output));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription d;
    d.add<edm::InputTag>("src", edm::InputTag("slimmedMuons"));
    d.add<double>("minPt", 2.0);
    descriptions.addWithDefaultLabel(d);
  }

private:
  edm::EDGetTokenT<pat::MuonCollection> source_;
  double minPt_;
};

DEFINE_FWK_MODULE(OniaMuonSelector);
