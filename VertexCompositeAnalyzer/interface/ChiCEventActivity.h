#ifndef VertexCompositeAnalysis_VertexCompositeAnalyzer_ChiCEventActivity_h
#define VertexCompositeAnalysis_VertexCompositeAnalyzer_ChiCEventActivity_h

#include <cmath>
#include <vector>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

namespace chic {
  struct TrackActivityConfig {
    double minPt;
    double maxAbsEta;
    double maxRelPtError;
    double maxDzSignificance;
    double maxDxySignificance;
  };

  inline std::vector<unsigned int> countAssociatedTracks(
      const reco::VertexCollection& vertices,
      const pat::PackedCandidateCollection& packedCandidates,
      const pat::PackedCandidateCollection& lostTracks,
      const TrackActivityConfig& config) {
    std::vector<unsigned int> result(vertices.size(), 0);
    const auto countCollection = [&](const pat::PackedCandidateCollection& candidates) {
      for (const auto& candidate : candidates) {
        if (candidate.charge() == 0 || !candidate.hasTrackDetails())
          continue;
        const reco::Track track = candidate.pseudoTrack();
        if (!track.quality(reco::TrackBase::highPurity) || track.pt() <= 0. || track.pt() <= config.minPt ||
            std::abs(track.eta()) >= config.maxAbsEta ||
            std::abs(track.ptError()) / track.pt() >= config.maxRelPtError)
          continue;
        for (size_t vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex) {
          const auto& vertex = vertices[vertexIndex];
          const double dz = track.dz(vertex.position());
          const double dzError2 = track.dzError() * track.dzError() + vertex.zError() * vertex.zError();
          if (!(dzError2 > 0.) || !std::isfinite(dz) ||
              dz * dz >= config.maxDzSignificance * config.maxDzSignificance * dzError2)
            continue;
          const double dxy = track.dxy(vertex.position());
          const double dxyError2 = track.dxyError() * track.dxyError() + vertex.xError() * vertex.yError();
          if (!(dxyError2 > 0.) || !std::isfinite(dxy) ||
              dxy * dxy >= config.maxDxySignificance * config.maxDxySignificance * dxyError2)
            continue;
          ++result[vertexIndex];
        }
      }
    };
    countCollection(packedCandidates);
    countCollection(lostTracks);
    return result;
  }
}  // namespace chic

#endif
