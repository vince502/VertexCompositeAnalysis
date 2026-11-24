#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/SingleLegConversionGeometry.h"

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

  reco::Candidate::LorentzVector visibleP4(const GlobalVector& direction, double momentum) {
    const GlobalVector unit = direction.mag2() > 0. ? direction.unit() : GlobalVector(0., 0., 0.);
    return {momentum * unit.x(),
            momentum * unit.y(),
            momentum * unit.z(),
            std::sqrt(momentum * momentum + kElectronMass * kElectronMass)};
  }
}  // namespace

class ChiCSingleLegGeometryProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCSingleLegGeometryProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const pat::PackedCandidate* sourceTrack(const pat::CompositeCandidate&,
                                          const pat::PackedCandidateCollection&,
                                          const pat::PackedCandidateCollection&) const;

  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsisToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> legsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
  double maxDeltaRToJpsi_;
  vca::SingleLegConversionGeometry geometry_;
};

ChiCSingleLegGeometryProducer::ChiCSingleLegGeometryProducer(const edm::ParameterSet& cfg)
    : jpsisToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("jpsis"))),
      legsToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("legs"))),
      packedToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxDeltaRToJpsi_(cfg.getParameter<double>("maxDeltaRToJpsi")),
      geometry_(cfg.getParameter<double>("minConversionRadius"),
                cfg.getParameter<double>("maxConversionRadius"),
                cfg.getParameter<double>("maxAbsConversionZ"),
                cfg.getParameter<double>("maxBackwardPath"),
                cfg.getParameter<unsigned int>("nGeometrySteps"),
                cfg.getParameter<double>("minDirectionCosine"),
                cfg.getParameter<unsigned int>("maxRoots")) {
  produces<pat::CompositeCandidateCollection>();
}

const pat::PackedCandidate* ChiCSingleLegGeometryProducer::sourceTrack(
    const pat::CompositeCandidate& leg,
    const pat::PackedCandidateCollection& packed,
    const pat::PackedCandidateCollection& lost) const {
  if (!leg.hasUserInt("sourceCollection") || !leg.hasUserInt("sourceKey"))
    return nullptr;
  const int source = leg.userInt("sourceCollection");
  const int key = leg.userInt("sourceKey");
  if (key < 0)
    return nullptr;
  if (source == 0 && static_cast<size_t>(key) < packed.size())
    return &packed[key];
  if (source == 1 && static_cast<size_t>(key) < lost.size())
    return &lost[key];
  return nullptr;
}

void ChiCSingleLegGeometryProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const {
  edm::Handle<pat::CompositeCandidateCollection> jpsisHandle;
  edm::Handle<pat::CompositeCandidateCollection> legsHandle;
  edm::Handle<pat::PackedCandidateCollection> packedHandle;
  edm::Handle<pat::PackedCandidateCollection> lostHandle;
  event.getByToken(jpsisToken_, jpsisHandle);
  event.getByToken(legsToken_, legsHandle);
  event.getByToken(packedToken_, packedHandle);
  event.getByToken(lostToken_, lostHandle);
  const auto& jpsis = *jpsisHandle;
  const auto& legs = *legsHandle;
  const auto& packed = *packedHandle;
  const auto& lost = *lostHandle;
  const auto& builder = setup.getData(transientTrackBuilderToken_);
  auto output = std::make_unique<pat::CompositeCandidateCollection>();

  for (size_t jpsiIndex = 0; jpsiIndex < jpsis.size(); ++jpsiIndex) {
    const auto& jpsi = jpsis[jpsiIndex];
    const GlobalPoint origin(jpsi.vx(), jpsi.vy(), jpsi.vz());
    for (size_t legIndex = 0; legIndex < legs.size(); ++legIndex) {
      const auto& leg = legs[legIndex];
      if (reco::deltaR(jpsi, leg) > maxDeltaRToJpsi_)
        continue;
      const auto* source = sourceTrack(leg, packed, lost);
      if (source == nullptr || !source->hasTrackDetails())
        continue;
      const reco::TransientTrack transientTrack = builder.build(extractPackedPseudoTrack(*source));
      // PackedCandidate pseudo-tracks do not carry a TrackExtra in MiniAOD.
      // Accessing innerPosition()/innerMomentum() would therefore dereference
      // an invalid TrackExtra reference for some valid charged candidates.
      // The packed candidate vertex and momentum are persistent MiniAOD data
      // and provide the required reference point and tangent for this seed.
      const auto& position = source->vertex();
      const auto& momentum = source->momentum();
      const vca::SingleLegGeometryInput input{origin,
                                              GlobalPoint(position.x(), position.y(), position.z()),
                                              GlobalVector(momentum.x(), momentum.y(), momentum.z())};
      const auto hypotheses = geometry_.solve(input);
      for (size_t root = 0; root < hypotheses.size(); ++root) {
        const auto& hypothesis = hypotheses[root];
        // The seed solver finds all geometrically allowed roots. At each root
        // use the MiniAOD pseudo-track in the magnetic field to obtain the
        // local helix tangent and its full Cartesian covariance. This is not
        // a hit refit, which MiniAOD cannot support.
        const auto closest = transientTrack.trajectoryStateClosestToPoint(hypothesis.conversionPoint);
        if (!closest.isValid() || !closest.hasError())
          continue;
        const auto conversionPoint = closest.position();
        const auto electronMomentum = closest.momentum();
        const GlobalVector flightVec(conversionPoint.x() - origin.x(),
                                     conversionPoint.y() - origin.y(),
                                     conversionPoint.z() - origin.z());
        if (flightVec.mag2() == 0.)
          continue;
        const GlobalVector photonDirection = flightVec.unit();
        const auto cartesianCovariance = closest.theState().cartesianError().matrix();
        const double tangentCosine =
            electronMomentum.mag() > 0. ? electronMomentum.unit().dot(photonDirection) : -2.;
        pat::CompositeCandidate geometry;
        geometry.setP4(visibleP4(photonDirection, leg.p()));
        geometry.setCharge(0);
        geometry.setVertex(reco::Candidate::Point(conversionPoint.x(), conversionPoint.y(), conversionPoint.z()));
        geometry.addDaughter(jpsi, "jpsi");
        geometry.addDaughter(leg, "leg");
        geometry.addUserInt("jpsiIndex", static_cast<int>(jpsiIndex));
        geometry.addUserInt("legIndex", static_cast<int>(legIndex));
        geometry.addUserInt("geometryRoot", static_cast<int>(root));
        geometry.addUserInt("geometryStatus", static_cast<int>(hypothesis.status));
        geometry.addUserFloat("conversionRadius", std::hypot(conversionPoint.x(), conversionPoint.y()));
        geometry.addUserFloat("conversionZ", conversionPoint.z());
        geometry.addUserFloat("conversionX", conversionPoint.x());
        geometry.addUserFloat("conversionY", conversionPoint.y());
        geometry.addUserFloat("geometryPathLength", hypothesis.pathLength);
        geometry.addUserFloat("directionCosine", tangentCosine);
        geometry.addUserFloat("angularMismatch",
                              tangentCosine >= -1. && tangentCosine <= 1.
                                  ? std::acos(tangentCosine)
                                  : std::numeric_limits<float>::quiet_NaN());
        geometry.addUserFloat("photonDirectionX", photonDirection.x());
        geometry.addUserFloat("photonDirectionY", photonDirection.y());
        geometry.addUserFloat("photonDirectionZ", photonDirection.z());
        geometry.addUserFloat("electronAtConversionPx", electronMomentum.x());
        geometry.addUserFloat("electronAtConversionPy", electronMomentum.y());
        geometry.addUserFloat("electronAtConversionPz", electronMomentum.z());
        for (unsigned int row = 0; row != 3; ++row)
          for (unsigned int column = 0; column <= row; ++column)
            geometry.addUserFloat("conversionCov" + std::to_string(row) + std::to_string(column),
                                  cartesianCovariance(row, column));
        geometry.addUserInt("conversionCovarianceStatus", 0x1);
        output->push_back(geometry);
      }
    }
  }
  event.put(std::move(output));
}

void ChiCSingleLegGeometryProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("jpsis", edm::InputTag("chiCJpsis"));
  description.add<edm::InputTag>("legs", edm::InputTag("chiCConversionLegs"));
  description.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  description.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  description.add<double>("maxDeltaRToJpsi", 1.5);
  description.add<double>("minConversionRadius", 1.0);
  description.add<double>("maxConversionRadius", 110.0);
  description.add<double>("maxAbsConversionZ", 280.0);
  description.add<double>("maxBackwardPath", 100.0);
  description.add<unsigned int>("nGeometrySteps", 64);
  description.add<double>("minDirectionCosine", 0.90);
  description.add<unsigned int>("maxRoots", 3);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCSingleLegGeometryProducer);
