#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/PoniaOniaPhotonProducer.h"

#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

PoniaOniaPhotonProducer::PoniaOniaPhotonProducer(const edm::ParameterSet& ps)
    : dimuon_Label_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("dimuons"))),
      photon_Label_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("conversions"))),
      pi0OnlineSwitch_(ps.getParameter<bool>("pi0OnlineSwitch")),
      deltaMass_(ps.getParameter<std::vector<double>>("deltaMass")),
      dzMax_(ps.getParameter<double>("dzmax")),
      triggerMatch_(ps.getParameter<bool>("triggerMatch")),
      minPhotonPt_(ps.getParameter<double>("minPhotonPt")),
      applyPairCuts_(ps.getParameter<bool>("applyPairCuts")),
      requireOniaVertex_(ps.getParameter<bool>("requireOniaVertex")),
      candidates_(0),
      delta_mass_fail_(0),
      dz_cut_fail_(0),
      pizero_fail_(0),
      missing_vertex_fail_(0) {
  if (applyPairCuts_ && deltaMass_.size() < 2) {
    throw cms::Exception("Configuration")
        << "PoniaOniaPhotonProducer: deltaMass must contain at least 2 elements [min, max] when applyPairCuts is true.";
  }
  produces<pat::CompositeCandidateCollection>();
}

void PoniaOniaPhotonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("dimuons", edm::InputTag("chiCSharedOniaStates"));
  desc.add<edm::InputTag>("conversions", edm::InputTag("oniaPhotonCandidates", "conversions"));
  desc.add<bool>("pi0OnlineSwitch", false);
  desc.add<std::vector<double>>("deltaMass", {0.0, 1.5});
  desc.add<double>("dzmax", 1.0);
  desc.add<bool>("triggerMatch", false);
  desc.add<double>("minPhotonPt", 0.0);
  desc.add<bool>("applyPairCuts", true);
  desc.add<bool>("requireOniaVertex", true);
  descriptions.addWithDefaultLabel(desc);
}

void PoniaOniaPhotonProducer::produce(edm::Event& event, const edm::EventSetup&) {
  auto chiCandColl = std::make_unique<pat::CompositeCandidateCollection>();

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByToken(dimuon_Label_, dimuons);

  edm::Handle<pat::CompositeCandidateCollection> conversions;
  event.getByToken(photon_Label_, conversions);

  if (!dimuons.isValid() || !conversions.isValid()) {
    event.put(std::move(chiCandColl));
    return;
  }

  int oniaIndex = 0;
  for (auto dimuonCand = dimuons->begin(); dimuonCand != dimuons->end(); ++dimuonCand, ++oniaIndex) {
    if (triggerMatch_) {
      if (!dimuonCand->hasUserInt("isTriggerMatched") || !dimuonCand->userInt("isTriggerMatched")) {
        continue;
      }
    }

    reco::Candidate::Point oniaVtxPos;
    bool hasValidOniaVertex = false;
    if (dimuonCand->hasUserData("commonVertex")) {
      const reco::Vertex* ipv = dimuonCand->userData<reco::Vertex>("commonVertex");
      if (ipv != nullptr) {
        oniaVtxPos = ipv->position();
        hasValidOniaVertex = true;
      }
    }
    if (!hasValidOniaVertex && dimuonCand->hasUserInt("jpsiFitStatus") &&
        (dimuonCand->userInt("jpsiFitStatus") & 0x1)) {
      oniaVtxPos = dimuonCand->vertex();
      hasValidOniaVertex = true;
    }
    if (!hasValidOniaVertex) {
      ++missing_vertex_fail_;
      if (requireOniaVertex_) {
        continue;
      }
    }

    const int resolvedOniaIdx =
        dimuonCand->hasUserInt("oniaStateIndex") ? dimuonCand->userInt("oniaStateIndex") : oniaIndex;

    int convIdx = 0;
    for (auto conv = conversions->begin(); conv != conversions->end(); ++conv, ++convIdx) {
      if (minPhotonPt_ > 0. && conv->pt() <= minPhotonPt_) {
        continue;
      }

      pat::CompositeCandidate chiCand =
          makeChiCandidate(*dimuonCand, *conv, hasValidOniaVertex ? &oniaVtxPos : nullptr);
      chiCand.addUserInt("oniaIdx", resolvedOniaIdx);
      chiCand.addUserInt("oniaIndex", resolvedOniaIdx);
      chiCand.addUserInt("convIdx", convIdx);
      chiCand.addUserInt("photonSourceIndex", convIdx);
      if (dimuonCand->hasUserFloat("nominalMass")) {
        chiCand.addUserFloat("nominalMass", dimuonCand->userFloat("nominalMass"));
      }
      if (dimuonCand->hasUserInt("oniaPdgId")) {
        chiCand.addUserInt("oniaPdgId", dimuonCand->userInt("oniaPdgId"));
      }

      if (applyPairCuts_ && !cutDeltaMass(chiCand, *dimuonCand)) {
        ++delta_mass_fail_;
        continue;
      }

      if (hasValidOniaVertex) {
        const float rawDz = Getdz(*conv, oniaVtxPos);
        const float dz = std::isfinite(rawDz) ? std::abs(rawDz) : std::numeric_limits<float>::quiet_NaN();
        chiCand.addUserFloat("dz", dz);
        if (applyPairCuts_ && !cutdz(dz)) {
          ++dz_cut_fail_;
          continue;
        }
      } else {
        chiCand.addUserFloat("dz", std::numeric_limits<float>::quiet_NaN());
      }

      const int flags = conv->hasUserInt("flags") ? (conv->userInt("flags") % 32) : 0;
      const bool pi0_fail = (flags & 8) != 0;
      if (pi0OnlineSwitch_ && pi0_fail) {
        ++pizero_fail_;
        continue;
      }

      chiCandColl->push_back(chiCand);
      ++candidates_;
    }
  }
  event.put(std::move(chiCandColl));
}

float PoniaOniaPhotonProducer::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point& p) const {
  const reco::Candidate::LorentzVector& mom = c.p4();
  const double rho = mom.Rho();
  if (rho <= 0.) {
    return std::numeric_limits<float>::quiet_NaN();
  }
  const reco::Candidate::Point& vtx = c.vertex();
  const double dz =
      (vtx.Z() - p.Z()) - ((vtx.X() - p.X()) * mom.X() + (vtx.Y() - p.Y()) * mom.Y()) / rho * mom.Z() / rho;
  return static_cast<float>(dz);
}

void PoniaOniaPhotonProducer::endStream() {
  edm::LogVerbatim("PoniaOniaPhotonProducer")
      << "###########################\n"
      << "Chi Candidate producer report:\n"
      << "###########################\n"
      << "Delta mass fail: " << delta_mass_fail_ << "\n"
      << "Dz fail:         " << dz_cut_fail_ << "\n"
      << "Pi0 fail:        " << pizero_fail_ << "\n"
      << "Missing vertex:  " << missing_vertex_fail_ << "\n"
      << "###########################\n"
      << "Found " << candidates_ << " Chi candidates.\n"
      << "###########################";
}

pat::CompositeCandidate PoniaOniaPhotonProducer::makeChiCandidate(
    const pat::CompositeCandidate& dimuon,
    const pat::CompositeCandidate& photon,
    const reco::Candidate::Point* oniaVertexPos) const {
  pat::CompositeCandidate chiCand;
  chiCand.addDaughter(dimuon, "dimuon");
  chiCand.addDaughter(photon, "photon");
  if (oniaVertexPos != nullptr) {
    chiCand.setVertex(*oniaVertexPos);
  }
  const reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
  chiCand.setP4(vChic);
  return chiCand;
}

bool PoniaOniaPhotonProducer::cutDeltaMass(const pat::CompositeCandidate& chiCand,
                                           const pat::CompositeCandidate& dimuonCand) const {
  const double deltam = chiCand.p4().M() - dimuonCand.p4().M();
  return (deltam > deltaMass_[0] && deltam < deltaMass_[1]);
}

DEFINE_FWK_MODULE(PoniaOniaPhotonProducer);
