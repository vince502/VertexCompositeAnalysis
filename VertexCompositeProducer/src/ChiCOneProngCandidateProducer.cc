#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace {
  constexpr std::array<double, 3> kChiCMasses{{3.41471, 3.51067, 3.55617}};

  using Matrix3 = std::array<std::array<double, 3>, 3>;
  using Matrix4 = std::array<std::array<double, 4>, 4>;
  using Vector3 = std::array<double, 3>;

  double userFloatOr(const pat::CompositeCandidate& candidate, const std::string& name, double fallback) {
    return candidate.hasUserFloat(name) ? candidate.userFloat(name) : fallback;
  }

  int userIntOr(const pat::CompositeCandidate& candidate, const std::string& name, int fallback) {
    return candidate.hasUserInt(name) ? candidate.userInt(name) : fallback;
  }

  Matrix3 covariance3(const pat::CompositeCandidate& candidate, const std::string& prefix) {
    Matrix3 result{};
    for (unsigned int row = 0; row != 3; ++row) {
      for (unsigned int column = 0; column <= row; ++column) {
        const double value = userFloatOr(candidate, prefix + "Cov" + std::to_string(row) + std::to_string(column), 0.);
        result[row][column] = value;
        result[column][row] = value;
      }
    }
    return result;
  }

  Matrix4 jpsiMomentumMassCovariance(const pat::CompositeCandidate& candidate) {
    Matrix4 result{};
    for (unsigned int row = 0; row != 4; ++row) {
      for (unsigned int column = 0; column <= row; ++column) {
        const double value = userFloatOr(
            candidate, "jpsiStateCov" + std::to_string(row + 3) + std::to_string(column + 3), 0.);
        result[row][column] = value;
        result[column][row] = value;
      }
    }
    return result;
  }

  Matrix3 add(const Matrix3& left, const Matrix3& right) {
    Matrix3 result{};
    for (unsigned int i = 0; i != 3; ++i)
      for (unsigned int j = 0; j != 3; ++j)
        result[i][j] = left[i][j] + right[i][j];
    return result;
  }

  double dot(const Vector3& left, const Vector3& right) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
  }

  Matrix3 directionCovariance(const Vector3& direction, double length, const Matrix3& pointCovariance) {
    Matrix3 jacobian{};
    if (!(length > 0.) || !std::isfinite(length))
      return jacobian;
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        jacobian[row][column] = ((row == column ? 1. : 0.) - direction[row] * direction[column]) / length;
    Matrix3 intermediate{};
    Matrix3 result{};
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        for (unsigned int k = 0; k != 3; ++k)
          intermediate[row][column] += jacobian[row][k] * pointCovariance[k][column];
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        for (unsigned int k = 0; k != 3; ++k)
          result[row][column] += intermediate[row][k] * jacobian[column][k];
    return result;
  }

  bool photonEnergy(double parentMass, const std::array<double, 4>& psi, const Vector3& direction, double& energy) {
    const double psiEnergy = std::sqrt(psi[0] * psi[0] + psi[1] * psi[1] + psi[2] * psi[2] + psi[3] * psi[3]);
    const double psiMass2 = psiEnergy * psiEnergy - psi[0] * psi[0] - psi[1] * psi[1] - psi[2] * psi[2];
    const double denominator = 2. * (psiEnergy - psi[0] * direction[0] - psi[1] * direction[1] - psi[2] * direction[2]);
    if (!(denominator > 0.) || !std::isfinite(denominator))
      return false;
    energy = (parentMass * parentMass - psiMass2) / denominator;
    return std::isfinite(energy) && energy > 0.;
  }

  double photonEnergyVariance(double parentMass,
                              const std::array<double, 4>& psi,
                              const Matrix4& psiCovariance,
                              const Vector3& direction,
                              const Matrix3& directionCovariance) {
    std::array<double, 4> momentumGradient{};
    Vector3 directionGradient{};
    for (unsigned int component = 0; component != 4; ++component) {
      const double step = std::max(1.e-6, std::abs(psi[component]) * 1.e-5);
      auto plus = psi;
      auto minus = psi;
      plus[component] += step;
      minus[component] -= step;
      double ePlus = 0., eMinus = 0.;
      if (photonEnergy(parentMass, plus, direction, ePlus) && photonEnergy(parentMass, minus, direction, eMinus))
        momentumGradient[component] = (ePlus - eMinus) / (2. * step);
    }
    for (unsigned int component = 0; component != 3; ++component) {
      const double step = 1.e-5;
      auto plus = direction;
      auto minus = direction;
      plus[component] += step;
      minus[component] -= step;
      double ePlus = 0., eMinus = 0.;
      if (photonEnergy(parentMass, psi, plus, ePlus) && photonEnergy(parentMass, psi, minus, eMinus))
        directionGradient[component] = (ePlus - eMinus) / (2. * step);
    }
    double variance = 0.;
    for (unsigned int row = 0; row != 4; ++row)
      for (unsigned int column = 0; column != 4; ++column)
        variance += momentumGradient[row] * psiCovariance[row][column] * momentumGradient[column];
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        variance += directionGradient[row] * directionCovariance[row][column] * directionGradient[column];
    return variance > 0. && std::isfinite(variance) ? variance : std::numeric_limits<double>::quiet_NaN();
  }

  void saveDirection(const std::string& prefix,
                     const Vector3& direction,
                     const Matrix3& covariance,
                     pat::CompositeCandidate& candidate) {
    candidate.addUserFloat(prefix + "X", direction[0]);
    candidate.addUserFloat(prefix + "Y", direction[1]);
    candidate.addUserFloat(prefix + "Z", direction[2]);
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        candidate.addUserFloat(prefix + "Cov" + std::to_string(row) + std::to_string(column), covariance[row][column]);
  }

  void saveMassHypotheses(const std::string& prefix,
                          const std::array<double, 4>& psi,
                          const Matrix4& psiCovariance,
                          const Vector3& direction,
                          const Matrix3& directionCovariance,
                          double observedElectronEnergy,
                          pat::CompositeCandidate& candidate) {
    for (unsigned int state = 0; state != kChiCMasses.size(); ++state) {
      double energy = std::numeric_limits<double>::quiet_NaN();
      const bool valid = photonEnergy(kChiCMasses[state], psi, direction, energy);
      const std::string label = prefix + "ChiC" + std::to_string(state);
      const double var =
          valid ? photonEnergyVariance(kChiCMasses[state], psi, psiCovariance, direction, directionCovariance)
                : std::numeric_limits<double>::quiet_NaN();
      candidate.addUserFloat(label + "PhotonEnergy", energy);
      candidate.addUserFloat(label + "PhotonEnergyError",
                             (valid && std::isfinite(var) && var >= 0.) ? std::sqrt(var)
                                                                        : std::numeric_limits<double>::quiet_NaN());
      candidate.addUserFloat(label + "MissingElectronEnergyProxy",
                             valid ? energy - observedElectronEnergy : std::numeric_limits<double>::quiet_NaN());
      candidate.addUserFloat(label + "ObservedEnergyFractionProxy",
                             (valid && energy > 0.) ? observedElectronEnergy / energy
                                                    : std::numeric_limits<double>::quiet_NaN());
      candidate.addUserInt(label + "MassHypothesisStatus", valid ? 0x1 : 0x2);
      if (valid) {
        candidate.addUserFloat(label + "PhotonPx", energy * direction[0]);
        candidate.addUserFloat(label + "PhotonPy", energy * direction[1]);
        candidate.addUserFloat(label + "PhotonPz", energy * direction[2]);
        candidate.addUserFloat(label + "ChiCPx", psi[0] + energy * direction[0]);
        candidate.addUserFloat(label + "ChiCPy", psi[1] + energy * direction[1]);
        candidate.addUserFloat(label + "ChiCPz", psi[2] + energy * direction[2]);
        candidate.addUserFloat(
            label + "ChiCE",
            std::sqrt(psi[0] * psi[0] + psi[1] * psi[1] + psi[2] * psi[2] + psi[3] * psi[3]) + energy);
      }
    }
  }
}  // namespace

class ChiCOneProngCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit ChiCOneProngCandidateProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::EDGetTokenT<pat::CompositeCandidateCollection> geometriesToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsisToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> legsToken_;
};

ChiCOneProngCandidateProducer::ChiCOneProngCandidateProducer(const edm::ParameterSet& cfg)
    : geometriesToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("geometries"))),
      jpsisToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("jpsis"))),
      legsToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("legs"))) {
  produces<pat::CompositeCandidateCollection>();
  produces<std::vector<int>>("eventSummary");
}

void ChiCOneProngCandidateProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& geometries = event.get(geometriesToken_);
  const auto& jpsis = event.get(jpsisToken_);
  const auto& legs = event.get(legsToken_);
  auto output = std::make_unique<pat::CompositeCandidateCollection>();
  for (const auto& geometry : geometries) {
    const auto* jpsi = dynamic_cast<const pat::CompositeCandidate*>(geometry.daughter("jpsi"));
    const auto* leg = dynamic_cast<const pat::CompositeCandidate*>(geometry.daughter("leg"));
    if (jpsi == nullptr || leg == nullptr)
      continue;
    const Vector3 conversion{{geometry.vx(), geometry.vy(), geometry.vz()}};
    const Vector3 freeOrigin{{jpsi->vx(), jpsi->vy(), jpsi->vz()}};
    const Vector3 promptOrigin{{userFloatOr(*jpsi, "promptOriginX", std::numeric_limits<double>::quiet_NaN()),
                                userFloatOr(*jpsi, "promptOriginY", std::numeric_limits<double>::quiet_NaN()),
                                userFloatOr(*jpsi, "promptOriginZ", std::numeric_limits<double>::quiet_NaN())}};
    const auto conversionCovariance = covariance3(geometry, "conversion");
    const auto freeCovariance = covariance3(*jpsi, "jpsiVertex");
    const auto promptCovariance = covariance3(*jpsi, "promptOrigin");
    const auto makeDirection = [](const Vector3& point, const Vector3& origin, Vector3& direction, double& length) {
      direction = {{point[0] - origin[0], point[1] - origin[1], point[2] - origin[2]}};
      const double len2 = dot(direction, direction);
      if (len2 > 0. && std::isfinite(len2)) {
        length = std::sqrt(len2);
        for (double& component : direction)
          component /= length;
        return true;
      }
      length = 0.;
      return false;
    };
    Vector3 freeDirection{}, promptDirection{};
    double freeLength = 0., promptLength = 0.;
    const bool validFree = makeDirection(conversion, freeOrigin, freeDirection, freeLength);
    const bool validPrompt = (userIntOr(*jpsi, "promptOriginStatus", 0) & 0x1) &&
                             makeDirection(conversion, promptOrigin, promptDirection, promptLength);
    if (!validFree)
      continue;
    const auto freeDirectionCovariance =
        directionCovariance(freeDirection, freeLength, add(conversionCovariance, freeCovariance));
    const auto promptDirectionCovariance =
        validPrompt ? directionCovariance(promptDirection, promptLength, add(conversionCovariance, promptCovariance))
                    : Matrix3{};

    pat::CompositeCandidate candidate(geometry);
    const auto rawP4 = jpsi->p4() + leg->p4();
    const auto visibleP4 = jpsi->p4() + geometry.p4();
    candidate.setP4(visibleP4);
    candidate.addUserFloat("rawThreeBodyMass", rawP4.mass());
    candidate.addUserFloat("rawPartialMassDifference", rawP4.mass() - jpsi->mass());
    candidate.addUserFloat("geometryVisibleMass", visibleP4.mass());
    candidate.addUserFloat("geometryVisibleMassDifference", visibleP4.mass() - jpsi->mass());
    candidate.addUserFloat("observedElectronEnergy", leg->energy());
    candidate.addUserFloat("visiblePhotonEnergy", geometry.energy());
    candidate.addUserInt("missingEnergyUnresolved", 1);
    candidate.addUserInt("missingElectronMeasured", 0);
    candidate.addUserInt("photonDirectionCrossCovarianceStatus", 0);  // C_{conversion,origin} unavailable in MiniAOD.
    saveDirection("freePhotonDirection", freeDirection, freeDirectionCovariance, candidate);
    if (validPrompt)
      saveDirection("promptPhotonDirection", promptDirection, promptDirectionCovariance, candidate);
    const std::array<double, 4> psi{{jpsi->px(), jpsi->py(), jpsi->pz(), jpsi->mass()}};
    const auto psiCovariance = jpsiMomentumMassCovariance(*jpsi);
    saveMassHypotheses("free", psi, psiCovariance, freeDirection, freeDirectionCovariance, leg->energy(), candidate);
    if (validPrompt)
      saveMassHypotheses(
          "prompt", psi, psiCovariance, promptDirection, promptDirectionCovariance, leg->energy(), candidate);
    const double directionDot = validPrompt ? std::max(-1., std::min(1., dot(freeDirection, promptDirection)))
                                            : std::numeric_limits<double>::quiet_NaN();
    candidate.addUserFloat("originDirectionDelta",
                           validPrompt ? std::acos(directionDot) : std::numeric_limits<double>::quiet_NaN());
    candidate.addUserFloat(
        "promptCompatibilityChi2",
        userFloatOr(*jpsi, "promptCompatibilityChi2", std::numeric_limits<double>::quiet_NaN()));
    candidate.addUserInt("freeHypothesisStatus", 0x1 | 0x4);  // valid, Gaussian covariance propagation (not a hit refit)
    candidate.addUserInt("promptHypothesisStatus", validPrompt ? (0x1 | 0x4 | 0x8) : 0x2);  // 0x8: block-diagonal PV/Jpsi covariance
    candidate.addUserInt("oneProngStatus",
                         userIntOr(geometry, "geometryStatus", 0) | 0x100 | (validPrompt ? 0x200 : 0));
    output->push_back(candidate);
  }
  const int nGeometry = static_cast<int>(output->size());
  event.put(std::move(output));
  event.put(std::make_unique<std::vector<int>>(std::initializer_list<int>{
                static_cast<int>(jpsis.size()), static_cast<int>(legs.size()), nGeometry}),
            "eventSummary");
}

void ChiCOneProngCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("geometries", edm::InputTag("chiCSingleLegGeometries"));
  description.add<edm::InputTag>("jpsis", edm::InputTag("chiCOniaStates"));
  description.add<edm::InputTag>("legs", edm::InputTag("chiCConversionLegs"));
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCOneProngCandidateProducer);
