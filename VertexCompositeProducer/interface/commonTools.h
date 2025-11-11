// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      commonTools
// 
/**\class commonTools commonTools.h VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h

 Description: Common utility functions for VertexCompositeProducer

 Implementation:
     Collection of helper functions and tools for vertex and track operations
*/
//
// Original Author:  Your Name
//
//

#ifndef VertexCompositeAnalysis__COMMON_TOOLS_H
#define VertexCompositeAnalysis__COMMON_TOOLS_H

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <utility>
#include <functional>
#include <cmath>

namespace VertexCompositeProducerCommonTools {

/**
 * @brief Criteria for selecting the best vertex
 */
enum class VertexSelectionCriteria {
  CLOSEST_DCA_3D,      // Select vertex closest to track in 3D DCA
  CLOSEST_DXY,         // Select vertex with smallest dxy to track
  CLOSEST_DZ,          // Select vertex with smallest dz to track
  WEIGHTED_SCORE,      // Weighted combination of DCA and vertex quality score
  BEST_VERTEX_SCORE    // Select based purely on vertex quality score (default)
};

/**
 * @brief Finds the best primary vertex from a collection
 * 
 * @param vertices The vertex collection to search
 * @param beamSpot The beam spot (used as fallback if no good PV found)
 * @param minTracks Minimum number of tracks required for valid PV (default: 5)
 * @return std::pair<math::XYZPoint, unsigned int> Best vertex position and its index
 *         Returns beamSpot position and index=0 if no valid PV found
 */
inline std::pair<math::XYZPoint, unsigned int> 
getBestVertex(const reco::VertexCollection& vertices, 
              const reco::BeamSpot& beamSpot,
              unsigned int minTracks = 5) 
{
  // Lambda to check if a vertex is valid PV
  auto isValidPV = [&minTracks](const reco::Vertex& vtx) -> bool {
    return !vtx.isFake() && vtx.tracksSize() >= minTracks;
  };

  // Lambda to score vertex quality (higher is better)
  auto vertexScore = [](const reco::Vertex& vtx) -> double {
    // Score based on number of tracks and chi2/ndof
    double nTracks = static_cast<double>(vtx.tracksSize());
    double chi2ndof = (vtx.ndof() > 0) ? vtx.chi2() / vtx.ndof() : 999.0;
    return nTracks / (1.0 + chi2ndof);
  };

  // Lambda to get vertex position as XYZPoint
  auto getPosition = [](const reco::Vertex& vtx) -> math::XYZPoint {
    return math::XYZPoint(vtx.x(), vtx.y(), vtx.z());
  };

  // Find best vertex
  if (vertices.empty()) {
    return std::make_pair(
      math::XYZPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()),
      0u
    );
  }

  // Check first vertex (usually the best)
  if (isValidPV(vertices[0])) {
    return std::make_pair(getPosition(vertices[0]), 0u);
  }

  // If first vertex not valid, search for best one
  unsigned int bestIdx = 0;
  double bestScore = -1.0;
  bool foundValid = false;

  for (unsigned int i = 0; i < vertices.size(); ++i) {
    if (isValidPV(vertices[i])) {
      double score = vertexScore(vertices[i]);
      if (score > bestScore) {
        bestScore = score;
        bestIdx = i;
        foundValid = true;
      }
    }
  }

  if (foundValid) {
    return std::make_pair(getPosition(vertices[bestIdx]), bestIdx);
  }

  // No valid PV found, use beamspot
  return std::make_pair(
    math::XYZPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()),
    0u
  );
}

/**
 * @brief Finds the best primary vertex from a collection using specified criteria
 * 
 * @param vertices The vertex collection to search
 * @param beamSpot The beam spot (used as fallback if no good PV found)
 * @param track Reference track for DCA-based criteria (can be nullptr for BEST_VERTEX_SCORE)
 * @param criteria Selection criteria to use
 * @param minTracks Minimum number of tracks required for valid PV (default: 5)
 * @param dcaWeight Weight for DCA in WEIGHTED_SCORE mode (default: 0.5)
 * @return std::pair<math::XYZPoint, unsigned int> Best vertex position and its index
 */
inline std::pair<math::XYZPoint, unsigned int> 
getBestVertex(const reco::VertexCollection& vertices, 
              const reco::BeamSpot& beamSpot,
              const reco::Track* track,
              VertexSelectionCriteria criteria = VertexSelectionCriteria::BEST_VERTEX_SCORE,
              unsigned int minTracks = 5,
              double dcaWeight = 0.5) 
{
  // Lambda to check if a vertex is valid PV
  auto isValidPV = [&minTracks](const reco::Vertex& vtx) -> bool {
    return !vtx.isFake() && vtx.tracksSize() >= minTracks;
  };

  // Lambda to score vertex quality (higher is better)
  auto vertexScore = [](const reco::Vertex& vtx) -> double {
    double nTracks = static_cast<double>(vtx.tracksSize());
    double chi2ndof = (vtx.ndof() > 0) ? vtx.chi2() / vtx.ndof() : 999.0;
    return nTracks / (1.0 + chi2ndof);
  };

  // Lambda to get vertex position as XYZPoint
  auto getPosition = [](const reco::Vertex& vtx) -> math::XYZPoint {
    return math::XYZPoint(vtx.x(), vtx.y(), vtx.z());
  };

  // Lambda to calculate DCA 3D
  auto calcDCA3D = [](const reco::Track* trk, const math::XYZPoint& vtxPos) -> double {
    if (!trk) return 999999.0;
    double dx = trk->vx() - vtxPos.x();
    double dy = trk->vy() - vtxPos.y();
    double dz = trk->vz() - vtxPos.z();
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  };

  // Lambda to calculate dxy
  auto calcDxy = [](const reco::Track* trk, const math::XYZPoint& vtxPos) -> double {
    if (!trk) return 999999.0;
    return std::abs(trk->dxy(vtxPos));
  };

  // Lambda to calculate dz
  auto calcDz = [](const reco::Track* trk, const math::XYZPoint& vtxPos) -> double {
    if (!trk) return 999999.0;
    return std::abs(trk->dz(vtxPos));
  };

  // Fallback to beamspot if no vertices
  if (vertices.empty()) {
    return std::make_pair(
      math::XYZPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()),
      0u
    );
  }

  // Select scoring function based on criteria
  unsigned int bestIdx = 0;
  double bestMetric = 999999.0;
  bool foundValid = false;
  bool minimizeMetric = true; // Most metrics should be minimized

  for (unsigned int i = 0; i < vertices.size(); ++i) {
    if (!isValidPV(vertices[i])) continue;

    math::XYZPoint vtxPos = getPosition(vertices[i]);
    double metric = 0.0;

    switch (criteria) {
      case VertexSelectionCriteria::CLOSEST_DCA_3D:
        metric = calcDCA3D(track, vtxPos);
        minimizeMetric = true;
        break;

      case VertexSelectionCriteria::CLOSEST_DXY:
        metric = calcDxy(track, vtxPos);
        minimizeMetric = true;
        break;

      case VertexSelectionCriteria::CLOSEST_DZ:
        metric = calcDz(track, vtxPos);
        minimizeMetric = true;
        break;

      case VertexSelectionCriteria::WEIGHTED_SCORE:
        {
          double dca3d = calcDCA3D(track, vtxPos);
          double vtxQuality = vertexScore(vertices[i]);
          // Normalize: lower DCA is better, higher vertex score is better
          // Combine: weighted sum where lower is better
          metric = dcaWeight * dca3d - (1.0 - dcaWeight) * vtxQuality;
          minimizeMetric = true;
        }
        break;

      case VertexSelectionCriteria::BEST_VERTEX_SCORE:
      default:
        metric = vertexScore(vertices[i]);
        minimizeMetric = false; // Higher score is better
        break;
    }

    // Update best vertex based on metric
    bool isBetter = minimizeMetric ? (metric < bestMetric) : (metric > bestMetric);
    
    if (!foundValid || isBetter) {
      bestMetric = metric;
      bestIdx = i;
      foundValid = true;
    }
  }

  if (foundValid) {
    return std::make_pair(getPosition(vertices[bestIdx]), bestIdx);
  }

  // No valid PV found, use beamspot
  return std::make_pair(
    math::XYZPoint(beamSpot.position().x(), beamSpot.position().y(), beamSpot.position().z()),
    0u
  );
}

/**
 * @brief Gets vertex position and errors
 * 
 * @param vtx The vertex
 * @param beamSpot The beam spot (used if vertex is invalid)
 * @param useBeamSpot Whether the beamspot should be used
 * @return std::tuple<math::XYZPoint, double, double, double> Position, xError, yError, zError
 */
inline std::tuple<math::XYZPoint, double, double, double>
getVertexPositionAndErrors(const reco::Vertex& vtx, 
                          const reco::BeamSpot& beamSpot,
                          bool useBeamSpot = false)
{
  if (useBeamSpot || vtx.isFake()) {
    return std::make_tuple(
      math::XYZPoint(beamSpot.position().x(), beamSpot.position().y(), 0.0),
      beamSpot.BeamWidthX(),
      beamSpot.BeamWidthY(),
      0.0
    );
  }

  return std::make_tuple(
    math::XYZPoint(vtx.x(), vtx.y(), vtx.z()),
    vtx.xError(),
    vtx.yError(),
    vtx.zError()
  );
}

/**
 * @brief Checks if vertex passes quality criteria
 * 
 * @param vtx The vertex to check
 * @param minTracks Minimum number of tracks
 * @param maxChi2NDF Maximum chi2/ndof
 * @return bool True if vertex passes criteria
 */
inline bool isGoodVertex(const reco::Vertex& vtx,
                        unsigned int minTracks = 5,
                        double maxChi2NDF = 10.0)
{
  if (vtx.isFake()) return false;
  if (vtx.tracksSize() < minTracks) return false;
  if (vtx.ndof() <= 0) return false;
  
  double chi2ndf = vtx.chi2() / vtx.ndof();
  if (chi2ndf > maxChi2NDF) return false;
  
  return true;
}

} // namespace VertexCompositeProducerCommonTools

#endif

