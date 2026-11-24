#ifndef VertexCompositeAnalysis_VertexCompositeProducer_SingleLegConversionGeometry_h
#define VertexCompositeAnalysis_VertexCompositeProducer_SingleLegConversionGeometry_h

#include <algorithm>
#include <cmath>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

namespace vca {

  struct SingleLegGeometryInput {
    GlobalPoint photonOrigin;
    GlobalPoint firstTrackPoint;
    GlobalVector trackMomentum;
  };

  struct SingleLegGeometryHypothesis {
    GlobalPoint conversionPoint;
    GlobalVector photonDirection;
    double pathLength = 0.;
    double directionCosine = -2.;
    unsigned int status = 0;
  };

  /**
   * A deterministic MiniAOD seed for a one-leg conversion point.
   *
   * MiniAOD has no hit-level trajectory, so this class scans backward from the
   * first available pseudo-track point along the measured tangent. It returns
   * every local maximum of the photon--electron directional compatibility.
   * The status word marks this as a straight-line seed; it is deliberately kept
   * separate from the EDProducer so it can be replaced by a material-aware
   * transient-track propagator after two-prong closure validation.
   */
  class SingleLegConversionGeometry {
  public:
    SingleLegConversionGeometry(double minRadius,
                                double maxRadius,
                                double maxAbsZ,
                                double maxBackwardPath,
                                unsigned int nSteps,
                                double minDirectionCosine,
                                unsigned int maxRoots)
        : minRadius_(minRadius),
          maxRadius_(maxRadius),
          maxAbsZ_(maxAbsZ),
          maxBackwardPath_(maxBackwardPath),
          nSteps_(std::max(2u, nSteps)),
          minDirectionCosine_(minDirectionCosine),
          maxRoots_(std::max(1u, maxRoots)) {}

    std::vector<SingleLegGeometryHypothesis> solve(const SingleLegGeometryInput& input) const {
      std::vector<SingleLegGeometryHypothesis> accepted;
      const double momMag = input.trackMomentum.mag();
      if (!(momMag > 0.) || !std::isfinite(momMag))
        return accepted;
      const GlobalVector tangent = input.trackMomentum.unit();
      std::vector<SingleLegGeometryHypothesis> scan;
      scan.reserve(nSteps_ + 1);
      for (unsigned int step = 0; step <= nSteps_; ++step) {
        const double path = maxBackwardPath_ * static_cast<double>(step) / nSteps_;
        const GlobalPoint point(input.firstTrackPoint.x() - path * tangent.x(),
                                input.firstTrackPoint.y() - path * tangent.y(),
                                input.firstTrackPoint.z() - path * tangent.z());
        const double radius = std::hypot(point.x(), point.y());
        if (radius < minRadius_ || radius > maxRadius_ || std::abs(point.z()) > maxAbsZ_)
          continue;
        const GlobalVector flight = point - input.photonOrigin;
        if (!(flight.mag2() > 0.) || !std::isfinite(flight.mag2()) || flight.dot(tangent) <= 0.)
          continue;
        const GlobalVector unitFlight = flight.unit();
        const double cosine = unitFlight.dot(tangent);
        if (cosine >= minDirectionCosine_)
          scan.push_back({point, unitFlight, path, cosine, 0x7u});
      }
      for (size_t index = 0; index < scan.size(); ++index) {
        const bool left = index == 0 || scan[index].directionCosine >= scan[index - 1].directionCosine;
        const bool right = index + 1 == scan.size() || scan[index].directionCosine > scan[index + 1].directionCosine;
        if (left && right)
          accepted.push_back(scan[index]);
      }
      std::sort(accepted.begin(), accepted.end(), [](const auto& a, const auto& b) {
        return a.directionCosine > b.directionCosine;
      });
      if (accepted.size() > maxRoots_)
        accepted.resize(maxRoots_);
      return accepted;
    }

  private:
    double minRadius_;
    double maxRadius_;
    double maxAbsZ_;
    double maxBackwardPath_;
    unsigned int nSteps_;
    double minDirectionCosine_;
    unsigned int maxRoots_;
  };

}  // namespace vca

#endif
