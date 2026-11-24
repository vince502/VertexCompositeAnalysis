#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace {
  constexpr double kMuonMass = 0.1056583745;

  using Matrix3 = std::array<std::array<double, 3>, 3>;
  using Vector3 = std::array<double, 3>;

  Matrix3 vertexCovariance(const reco::Vertex& vertex) {
    Matrix3 result{};
    for (unsigned int i = 0; i != 3; ++i)
      for (unsigned int j = 0; j != 3; ++j)
        result[i][j] = vertex.covariance(i, j);
    return result;
  }

  bool invertCovariance3(const Matrix3& input, Matrix3& inverse) {
    // Covariances must be symmetric positive definite. A determinant-only
    // inverse can accept indefinite matrices and is unstable near singular
    // ones, so use a checked Cholesky decomposition instead.
    Matrix3 lower{};
    double largestDiagonal = 0.;
    for (unsigned int i = 0; i != 3; ++i) {
      if (!std::isfinite(input[i][i]) || input[i][i] <= 0.)
        return false;
      largestDiagonal = std::max(largestDiagonal, input[i][i]);
    }
    const double minimumPivot = largestDiagonal * 1.e-14;
    for (unsigned int row = 0; row != 3; ++row) {
      for (unsigned int column = 0; column <= row; ++column) {
        const double symmetricInput = 0.5 * (input[row][column] + input[column][row]);
        if (!std::isfinite(symmetricInput))
          return false;
        double value = symmetricInput;
        for (unsigned int inner = 0; inner < column; ++inner)
          value -= lower[row][inner] * lower[column][inner];
        if (row == column) {
          if (!std::isfinite(value) || value <= minimumPivot)
            return false;
          lower[row][column] = std::sqrt(value);
        } else {
          lower[row][column] = value / lower[column][column];
        }
      }
    }

    inverse = Matrix3{};
    for (unsigned int rhs = 0; rhs != 3; ++rhs) {
      Vector3 forward{};
      for (unsigned int row = 0; row != 3; ++row) {
        double value = row == rhs ? 1. : 0.;
        for (unsigned int column = 0; column < row; ++column)
          value -= lower[row][column] * forward[column];
        forward[row] = value / lower[row][row];
      }
      for (int row = 2; row >= 0; --row) {
        double value = forward[row];
        for (unsigned int column = row + 1; column != 3; ++column)
          value -= lower[column][row] * inverse[column][rhs];
        inverse[row][rhs] = value / lower[row][row];
      }
    }
    for (unsigned int row = 0; row != 3; ++row) {
      for (unsigned int column = 0; column < row; ++column) {
        const double symmetric = 0.5 * (inverse[row][column] + inverse[column][row]);
        inverse[row][column] = symmetric;
        inverse[column][row] = symmetric;
      }
      for (unsigned int column = 0; column != 3; ++column)
        if (!std::isfinite(inverse[row][column]))
          return false;
    }
    return true;
  }

  Matrix3 add(const Matrix3& left, const Matrix3& right) {
    Matrix3 result{};
    for (unsigned int i = 0; i != 3; ++i)
      for (unsigned int j = 0; j != 3; ++j)
        result[i][j] = left[i][j] + right[i][j];
    return result;
  }

  Vector3 multiply(const Matrix3& matrix, const Vector3& vector) {
    Vector3 result{};
    for (unsigned int i = 0; i != 3; ++i)
      for (unsigned int j = 0; j != 3; ++j)
        result[i] += matrix[i][j] * vector[j];
    return result;
  }

  double dot(const Vector3& left, const Vector3& right) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
  }

  void saveCovariance(const std::string& prefix, const Matrix3& covariance, pat::CompositeCandidate& candidate) {
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        candidate.addUserFloat(prefix + "Cov" + std::to_string(row) + std::to_string(column), covariance[row][column]);
  }

  void saveVertex(const std::string& prefix,
                  const GlobalPoint& point,
                  const Matrix3& covariance,
                  pat::CompositeCandidate& candidate) {
    candidate.addUserFloat(prefix + "X", point.x());
    candidate.addUserFloat(prefix + "Y", point.y());
    candidate.addUserFloat(prefix + "Z", point.z());
    saveCovariance(prefix, covariance, candidate);
  }
}  // namespace

// Materialises the fit state that legacy Onia2MuMu candidates keep only as
// opaque user data. The output has the same muon daughters and selection as
// the conventional branch, while carrying an explicit 7x7 (x,y,z,px,py,pz,m)
// covariance and a separately labelled prompt common-origin hypothesis.
class ChiCOniaStateProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCOniaStateProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsisToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
  bool useEmbeddedPVIndex_;
  bool useOutputIndices_;
  std::vector<int> oniaPDGs_;
  std::vector<double> massRanges_;
  std::vector<double> nominalMasses_;
  bool oppositeSign_;
};

ChiCOniaStateProducer::ChiCOniaStateProducer(const edm::ParameterSet& cfg)
    : jpsisToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("jpsis"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      useEmbeddedPVIndex_(cfg.getParameter<bool>("useEmbeddedPVIndex")),
      useOutputIndices_(cfg.getParameter<bool>("useOutputIndices")),
      oniaPDGs_(cfg.getParameter<std::vector<int>>("oniaPDGs")),
      massRanges_(cfg.getParameter<std::vector<double>>("massRanges")),
      nominalMasses_(cfg.getParameter<std::vector<double>>("nominalMasses")),
      oppositeSign_(cfg.getParameter<bool>("oppositeSign")) {
  if (!oniaPDGs_.empty()) {
    if (massRanges_.size() != 2 * oniaPDGs_.size() || nominalMasses_.size() != oniaPDGs_.size()) {
      throw cms::Exception("Configuration")
          << "One dimuon window and nominal mass are required per Onia hypothesis";
    }
  }
  produces<pat::CompositeCandidateCollection>();
}

void ChiCOniaStateProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const {
  const auto& input = event.get(jpsisToken_);
  const auto& vertices = event.get(verticesToken_);
  const auto& builder = setup.getData(transientTrackBuilderToken_);
  auto output = std::make_unique<pat::CompositeCandidateCollection>();

  for (size_t index = 0; index < input.size(); ++index) {
    const auto& original = input[index];

    int hypothesis = -1;
    if (!oniaPDGs_.empty()) {
      if (oppositeSign_ && original.charge() != 0)
        continue;
      for (size_t h = 0; h < oniaPDGs_.size(); ++h) {
        if (original.mass() > massRanges_[2 * h] && original.mass() < massRanges_[2 * h + 1]) {
          if (hypothesis != -1)
            throw cms::Exception("Configuration") << "Overlapping Onia mass hypotheses";
          hypothesis = static_cast<int>(h);
        }
      }
      if (hypothesis == -1)
        continue;
    }

    const auto* muon1 = dynamic_cast<const pat::Muon*>(original.daughter("muon1"));
    const auto* muon2 = dynamic_cast<const pat::Muon*>(original.daughter("muon2"));
    if (muon1 == nullptr || muon2 == nullptr || muon1->innerTrack().isNull() || muon2->innerTrack().isNull())
      continue;

    std::vector<reco::TransientTrack> tracks;
    tracks.emplace_back(builder.build(muon1->innerTrack()));
    tracks.emplace_back(builder.build(muon2->innerTrack()));
    const float massUncertainty = static_cast<float>(kMuonMass * 1.e-6);
    KinematicParticleFactoryFromTransientTrack factory;
    std::vector<RefCountedKinematicParticle> particles;
    particles.emplace_back(factory.particle(tracks[0], kMuonMass, 0.f, 0.f, massUncertainty));
    particles.emplace_back(factory.particle(tracks[1], kMuonMass, 0.f, 0.f, massUncertainty));
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree tree = fitter.fit(particles);
    if (tree->isEmpty() || !tree->isValid())
      continue;
    tree->movePointerToTheTop();
    const auto fitted = tree->currentParticle();
    const auto decayVertex = tree->currentDecayVertex();
    if (!fitted || !decayVertex || !fitted->currentState().isValid() || !decayVertex->vertexIsValid())
      continue;
    const auto& state = fitted->currentState();
    const auto& parameters = state.kinematicParameters();
    const auto& error = state.kinematicParametersError().matrix();
    const GlobalPoint psiPoint = parameters.position();
    const GlobalVector psiMomentum = parameters.momentum();
    Matrix3 psiCov{};
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        psiCov[row][column] = error(row, column);

    int pvIndex = -1;
    reco::Vertex selectedPV;
    bool selectedPVValid = false;
    const reco::Vertex* embeddedPV =
        original.hasUserData("PVwithmuons") ? original.userData<reco::Vertex>("PVwithmuons") : nullptr;
    if (embeddedPV != nullptr) {
      selectedPV = *embeddedPV;
      selectedPVValid = true;
    }

    if (useEmbeddedPVIndex_ && selectedPVValid) {
      // PVwithmuons is stored by value. Recover its collection index from
      // that value itself; never attach the index of a different, merely
      // J/psi-compatible vertex to the stored coordinates.
      constexpr double kVertexMatchTolerance2 = 1.e-20;      // (1e-10 cm)^2
      constexpr double kCovarianceMatchTolerance2 = 1.e-28;  // (1e-14 cm^2)^2
      double bestCovarianceDistance2 = std::numeric_limits<double>::infinity();
      for (size_t pv = 0; pv < vertices.size(); ++pv) {
        const double dx = selectedPV.x() - vertices[pv].x();
        const double dy = selectedPV.y() - vertices[pv].y();
        const double dz = selectedPV.z() - vertices[pv].z();
        const double distance2 = dx * dx + dy * dy + dz * dz;
        if (distance2 > kVertexMatchTolerance2)
          continue;
        double covarianceDistance2 = 0.;
        for (unsigned int row = 0; row != 3; ++row)
          for (unsigned int column = 0; column <= row; ++column) {
            const double difference = selectedPV.covariance(row, column) - vertices[pv].covariance(row, column);
            covarianceDistance2 += difference * difference;
          }
        if (covarianceDistance2 < bestCovarianceDistance2) {
          bestCovarianceDistance2 = covarianceDistance2;
          pvIndex = static_cast<int>(pv);
        }
      }
      if (pvIndex >= 0 && bestCovarianceDistance2 <= kCovarianceMatchTolerance2)
        selectedPV = vertices[pvIndex];  // make value/index identity exact
      else
        pvIndex = -1;
    } else {
      double bestChi2 = std::numeric_limits<double>::infinity();
      for (size_t pv = 0; pv < vertices.size(); ++pv) {
        const auto sum = add(psiCov, vertexCovariance(vertices[pv]));
        Matrix3 inverse{};
        if (!invertCovariance3(sum, inverse))
          continue;
        const Vector3 displacement{{psiPoint.x() - vertices[pv].x(),
                                    psiPoint.y() - vertices[pv].y(),
                                    psiPoint.z() - vertices[pv].z()}};
        const double chi2 = dot(displacement, multiply(inverse, displacement));
        if (std::isfinite(chi2) && chi2 >= 0. && chi2 < bestChi2) {
          bestChi2 = chi2;
          pvIndex = static_cast<int>(pv);
          if (!selectedPVValid)
            selectedPV = vertices[pv];
        }
      }
      if (!selectedPVValid && pvIndex >= 0)
        selectedPVValid = true;
    }

    pat::CompositeCandidate candidate(original);
    candidate.setP4(reco::Candidate::LorentzVector(
        psiMomentum.x(), psiMomentum.y(), psiMomentum.z(), parameters.energy()));
    candidate.setVertex(reco::Candidate::Point(psiPoint.x(), psiPoint.y(), psiPoint.z()));
    candidate.addUserInt("oniaStateIndex", static_cast<int>(useOutputIndices_ ? output->size() : index));
    candidate.addUserInt("oniaInputIndex", static_cast<int>(index));
    if (hypothesis >= 0) {
      candidate.addUserInt("oniaPdgId", oniaPDGs_[hypothesis]);
      candidate.addUserFloat("nominalMass", nominalMasses_[hypothesis]);
    }
    candidate.addUserInt("jpsiFitStatus", 0x1 | (selectedPVValid ? 0x2 : 0));
    candidate.addUserInt("pvIndex", pvIndex);
    candidate.addUserInt("pvJpsiCrossCovarianceStatus", 0);  // MiniAOD legacy PV may include the muons; no cross block is available.
    candidate.addUserFloat("rawJpsiPx", original.px());
    candidate.addUserFloat("rawJpsiPy", original.py());
    candidate.addUserFloat("rawJpsiPz", original.pz());
    candidate.addUserFloat("rawJpsiE", original.energy());
    candidate.addUserFloat("jpsiFitChi2", decayVertex->chiSquared());
    candidate.addUserFloat("jpsiFitNdof", decayVertex->degreesOfFreedom());
    candidate.addUserFloat("jpsiFitMass", parameters.mass());
    candidate.addUserFloat("jpsiFitPx", psiMomentum.x());
    candidate.addUserFloat("jpsiFitPy", psiMomentum.y());
    candidate.addUserFloat("jpsiFitPz", psiMomentum.z());
    candidate.addUserFloat("jpsiFitE", parameters.energy());
    for (unsigned int row = 0; row != 7; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        candidate.addUserFloat("jpsiStateCov" + std::to_string(row) + std::to_string(column), error(row, column));
    saveVertex("jpsiVertex", psiPoint, psiCov, candidate);

    if (selectedPVValid) {
      const Matrix3 pvCov = vertexCovariance(selectedPV);
      const Matrix3 sum = add(psiCov, pvCov);  // Explicit block-diagonal C_{PV,psi} approximation.
      Matrix3 sumInverse{};
      const Vector3 displacement{{psiPoint.x() - selectedPV.x(),
                                  psiPoint.y() - selectedPV.y(),
                                  psiPoint.z() - selectedPV.z()}};
      saveVertex("primaryVertex", GlobalPoint(selectedPV.x(), selectedPV.y(), selectedPV.z()), pvCov, candidate);
      Matrix3 psiInverse{};
      Matrix3 pvInverse{};
      if (invertCovariance3(sum, sumInverse) && invertCovariance3(psiCov, psiInverse) &&
          invertCovariance3(pvCov, pvInverse)) {
        const Matrix3 information = add(psiInverse, pvInverse);
        Matrix3 commonCov{};
        if (invertCovariance3(information, commonCov)) {
          const Vector3 rhsPsi = multiply(psiInverse, Vector3{{psiPoint.x(), psiPoint.y(), psiPoint.z()}});
          const Vector3 rhsPV = multiply(pvInverse, Vector3{{selectedPV.x(), selectedPV.y(), selectedPV.z()}});
          const Vector3 common =
              multiply(commonCov, Vector3{{rhsPsi[0] + rhsPV[0], rhsPsi[1] + rhsPV[1], rhsPsi[2] + rhsPV[2]}});
          const double compatibilityChi2 = dot(displacement, multiply(sumInverse, displacement));
          if (std::isfinite(common[0]) && std::isfinite(common[1]) && std::isfinite(common[2]) &&
              std::isfinite(compatibilityChi2) && compatibilityChi2 >= -1.e-9) {
            saveVertex("promptOrigin", GlobalPoint(common[0], common[1], common[2]), commonCov, candidate);
            candidate.addUserFloat("promptCompatibilityChi2", std::max(0., compatibilityChi2));
            candidate.addUserInt("promptOriginStatus", 0x1 | 0x4);  // valid, block-diagonal PV--J/psi covariance
          } else {
            candidate.addUserInt("promptOriginStatus", 0x8);
          }
        } else {
          candidate.addUserInt("promptOriginStatus", 0x8);
        }
      } else {
        candidate.addUserInt("promptOriginStatus", 0x8);
      }
    } else {
      candidate.addUserInt("promptOriginStatus", 0x8);
    }
    output->push_back(candidate);
  }
  event.put(std::move(output));
}

void ChiCOniaStateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("jpsis", edm::InputTag("chiCSharedOniaJpsis"));
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<bool>("useEmbeddedPVIndex", true);
  description.add<bool>("useOutputIndices", false);
  description.add<std::vector<int>>("oniaPDGs", {});
  description.add<std::vector<double>>("massRanges", {});
  description.add<std::vector<double>>("nominalMasses", {});
  description.add<bool>("oppositeSign", false);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCOniaStateProducer);
