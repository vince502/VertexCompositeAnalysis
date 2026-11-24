#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

class DiMuonFilter : public edm::stream::EDProducer<> {
public:
  explicit DiMuonFilter(const edm::ParameterSet&);
  ~DiMuonFilter() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  unsigned int isTriggerMatched(const pat::CompositeCandidate* diMuon_cand) const;

  edm::EDGetTokenT<pat::CompositeCandidateCollection> theOnias_;
  StringCutObjectSelector<reco::Candidate, true> SingleMuonSelection_;
  StringCutObjectSelector<reco::Candidate, true> DiMuonSelection_;
  bool do_trigger_match_;
  std::vector<std::string> HLTFilters_;
};

DiMuonFilter::DiMuonFilter(const edm::ParameterSet& iConfig)
    : theOnias_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaTag"))),
      SingleMuonSelection_(iConfig.getParameter<std::string>("singlemuonSelection")),
      DiMuonSelection_(iConfig.getParameter<std::string>("dimuonSelection")),
      do_trigger_match_(iConfig.getParameter<bool>("do_trigger_match")),
      HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters")) {
  produces<pat::CompositeCandidateCollection>();
}

void DiMuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("OniaTag", edm::InputTag("chiCSharedOnia2MuMu"));
  desc.add<std::string>("singlemuonSelection", "");
  desc.add<std::string>("dimuonSelection", "");
  desc.add<bool>("do_trigger_match", false);
  desc.add<std::vector<std::string>>("HLTFilters", {});
  descriptions.addWithDefaultLabel(desc);
}

unsigned int DiMuonFilter::isTriggerMatched(const pat::CompositeCandidate* diMuon_cand) const {
  if (diMuon_cand == nullptr)
    return 0;
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  if (muon1 == nullptr || muon2 == nullptr)
    return 0;
  unsigned int matched = 0;
  for (unsigned int iTr = 0; iTr < HLTFilters_.size(); ++iTr) {
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) {
      matched += (1u << iTr);
    }
  }
  return matched;
}

void DiMuonFilter::produce(edm::Event& iEvent, const edm::EventSetup&) {
  auto mumuOutput = std::make_unique<pat::CompositeCandidateCollection>();
  edm::Handle<pat::CompositeCandidateCollection> onias;
  iEvent.getByToken(theOnias_, onias);
  if (onias.isValid() && !onias->empty()) {
    for (size_t ii = 0, nn = onias->size(); ii < nn; ++ii) {
      const pat::CompositeCandidate& onia = onias->at(ii);
      const reco::Candidate* mu1 = onia.daughter("muon1");
      const reco::Candidate* mu2 = onia.daughter("muon2");
      if (mu1 == nullptr || mu2 == nullptr)
        continue;
      if (!DiMuonSelection_(onia) || !SingleMuonSelection_(*mu1) || !SingleMuonSelection_(*mu2))
        continue;
      const unsigned int matched = do_trigger_match_ ? isTriggerMatched(&onia) : 0u;
      if (do_trigger_match_ && matched == 0u)
        continue;
      pat::CompositeCandidate outCand(onia);
      if (do_trigger_match_ && !outCand.hasUserInt("isTriggerMatched")) {
        outCand.addUserInt("isTriggerMatched", static_cast<int>(matched));
      }
      mumuOutput->push_back(outCand);
    }
  }
  iEvent.put(std::move(mumuOutput));
}

DEFINE_FWK_MODULE(DiMuonFilter);
