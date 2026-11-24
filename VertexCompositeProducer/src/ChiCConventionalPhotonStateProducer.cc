#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace {
  constexpr double kElectronMass = 0.00051099895;

  void saveVertexCovariance(const GlobalError& covariance, pat::CompositeCandidate& candidate) {
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        candidate.addUserFloat("conventionalConversionCov" + std::to_string(row) + std::to_string(column),
                               covariance.matrix()(row, column));
  }
}  // namespace

// Rebuilds the two-prong conversion state from the persistent MiniAOD tracks.
// It intentionally does not claim a hit-level/Gsf refit: the fitted covariance
// is the kinematic fit covariance of the two available reco::Track states.
class ChiCConventionalPhotonStateProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCConventionalPhotonStateProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<pat::CompositeCandidateCollection> candidatesToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
};

ChiCConventionalPhotonStateProducer::ChiCConventionalPhotonStateProducer(const edm::ParameterSet& cfg)
    : candidatesToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("candidates"))),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))) {
  produces<pat::CompositeCandidateCollection>();
}

void ChiCConventionalPhotonStateProducer::produce(edm::StreamID,
                                                  edm::Event& event,
                                                  const edm::EventSetup& setup) const {
  const auto& input = event.get(candidatesToken_);
  const auto& builder = setup.getData(transientTrackBuilderToken_);
  auto output = std::make_unique<pat::CompositeCandidateCollection>();
  for (const auto& original : input) {
    const auto* photon = dynamic_cast<const pat::CompositeCandidate*>(original.daughter("photon"));
    const auto* jpsi = original.daughter("dimuon");
    if (photon == nullptr || jpsi == nullptr || !photon->hasUserData("track0") || !photon->hasUserData("track1"))
      continue;
    const reco::Track* track0 = photon->userData<reco::Track>("track0");
    const reco::Track* track1 = photon->userData<reco::Track>("track1");
    if (track0 == nullptr || track1 == nullptr)
      continue;
    const float massUncertainty = static_cast<float>(kElectronMass * 1.e-6);
    KinematicParticleFactoryFromTransientTrack factory;
    std::vector<RefCountedKinematicParticle> particles;
    particles.emplace_back(factory.particle(builder.build(*track0), kElectronMass, 0.f, 0.f, massUncertainty));
    particles.emplace_back(factory.particle(builder.build(*track1), kElectronMass, 0.f, 0.f, massUncertainty));
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree tree = fitter.fit(particles);
    if (tree->isEmpty() || !tree->isValid())
      continue;
    tree->movePointerToTheTop();
    const auto fittedPhoton = tree->currentParticle();
    const auto fittedVertex = tree->currentDecayVertex();
    if (!fittedPhoton || !fittedVertex || !fittedPhoton->currentState().isValid() || !fittedVertex->vertexIsValid())
      continue;
    const auto& state = fittedPhoton->currentState();
    const auto& parameters = state.kinematicParameters();
    const auto& covariance = state.kinematicParametersError().matrix();
    const auto vertex = fittedVertex->position();
    pat::CompositeCandidate candidate(original);
    const reco::Candidate::LorentzVector photonP4(parameters.momentum().x(),
                                                  parameters.momentum().y(),
                                                  parameters.momentum().z(),
                                                  parameters.energy());
    candidate.setP4(jpsi->p4() + photonP4);
    candidate.addUserFloat("conventionalPhotonPx", photonP4.px());
    candidate.addUserFloat("conventionalPhotonPy", photonP4.py());
    candidate.addUserFloat("conventionalPhotonPz", photonP4.pz());
    candidate.addUserFloat("conventionalPhotonE", photonP4.energy());
    candidate.addUserFloat("conventionalPhotonMass", parameters.mass());
    candidate.addUserFloat("conventionalConversionX", vertex.x());
    candidate.addUserFloat("conventionalConversionY", vertex.y());
    candidate.addUserFloat("conventionalConversionZ", vertex.z());
    candidate.addUserFloat("conventionalPhotonFitChi2", fittedVertex->chiSquared());
    candidate.addUserFloat("conventionalPhotonFitNdof", fittedVertex->degreesOfFreedom());
    saveVertexCovariance(fittedVertex->error(), candidate);
    for (unsigned int row = 0; row != 7; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        candidate.addUserFloat("conventionalPhotonStateCov" + std::to_string(row) + std::to_string(column),
                               covariance(row, column));
    candidate.addUserInt("conventionalPhotonStateStatus", 0x1);
    output->push_back(candidate);
  }
  event.put(std::move(output));
}

void ChiCConventionalPhotonStateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("candidates", edm::InputTag("chiCConventionalCandidates"));
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCConventionalPhotonStateProducer);
