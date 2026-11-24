#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/interface/ChiCEventActivity.h"
#include "TTree.h"

namespace {
  constexpr float kElectronMass = 0.00051099895f;
  constexpr float kMissing = std::numeric_limits<float>::quiet_NaN();

  float value(const pat::CompositeCandidate& candidate, const std::string& name, float fallback = kMissing) {
    return candidate.hasUserFloat(name) ? candidate.userFloat(name) : fallback;
  }

  int intValue(const pat::CompositeCandidate& candidate, const std::string& name, int fallback = -1) {
    return candidate.hasUserInt(name) ? candidate.userInt(name) : fallback;
  }

  float electronEnergy(float px, float py, float pz) {
    return std::sqrt(px * px + py * py + pz * pz + kElectronMass * kElectronMass);
  }

  // The free and prompt photon-direction covariances are stored separately.
  // Their cross block is unavailable in MiniAOD and is explicitly flagged in
  // photonDirectionCrossCovStatus_. This returns the corresponding
  // block-diagonal uncertainty of acos(n_free . n_prompt).
  float originAngleError(const pat::CompositeCandidate& candidate) {
    const float angle = value(candidate, "originDirectionDelta");
    if (!std::isfinite(angle) || std::abs(std::sin(angle)) < 1.e-6f)
      return kMissing;
    float freeDirection[3]{};
    float promptDirection[3]{};
    float freeCovariance[3][3]{};
    float promptCovariance[3][3]{};
    for (unsigned int row = 0; row != 3; ++row) {
      freeDirection[row] = value(candidate, "freePhotonDirection" + std::string(1, "XYZ"[row]));
      promptDirection[row] = value(candidate, "promptPhotonDirection" + std::string(1, "XYZ"[row]));
      for (unsigned int column = 0; column != 3; ++column) {
        const std::string freeName =
            "freePhotonDirectionCov" + std::to_string(std::max(row, column)) + std::to_string(std::min(row, column));
        const std::string promptName =
            "promptPhotonDirectionCov" + std::to_string(std::max(row, column)) + std::to_string(std::min(row, column));
        freeCovariance[row][column] = value(candidate, freeName);
        promptCovariance[row][column] = value(candidate, promptName);
      }
    }
    const float cosine = std::cos(angle);
    float freeGradient[3]{};
    float promptGradient[3]{};
    for (unsigned int i = 0; i != 3; ++i) {
      freeGradient[i] = -(promptDirection[i] - cosine * freeDirection[i]) / std::sin(angle);
      promptGradient[i] = -(freeDirection[i] - cosine * promptDirection[i]) / std::sin(angle);
    }
    float variance = 0.f;
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column != 3; ++column)
        variance += freeGradient[row] * freeCovariance[row][column] * freeGradient[column] +
                    promptGradient[row] * promptCovariance[row][column] * promptGradient[column];
    return variance > 0.f && std::isfinite(variance) ? std::sqrt(variance) : kMissing;
  }

  template <typename Covariance>
  void appendLowerTriangle(const Covariance& covariance, unsigned int dimension, std::vector<float>& values) {
    for (unsigned int row = 0; row != dimension; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        values.push_back(covariance(row, column));
  }
}  // namespace

// A deliberately plain ROOT schema: every physics payload is a std::vector
// and every flattened relation has begin/size branches. It needs no custom
// dictionaries or analysis-side CMSSW classes and can be read directly with
// uproot/awkward arrays.
class ChiCCovarianceNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCCovarianceNtuplizer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  void clear();
  void appendTrackCovariance(const reco::Track&,
                             std::vector<unsigned int>& covBegin,
                             std::vector<unsigned int>& covSize,
                             std::vector<float>& covariances);
  void appendMuon(const pat::Muon&, unsigned int jpsiIndex);
  int oniaMuonIndex(const pat::Muon&) const;
  int appendOniaLinks(int oniaIndex);
  void appendElectronFromTrack(const reco::Track&, unsigned int photonIndex, int source, int sourceKey);
  void appendElectronFromOneProng(const pat::CompositeCandidate&, unsigned int photonIndex);
  int oniaIndex(const reco::Candidate*);
  void appendConventional(const pat::CompositeCandidate&);
  void appendOneProng(const pat::CompositeCandidate&);
  void appendTruth(const pat::CompositeCandidate&, int);

  struct TruthMatch {
    int jpsi = -1;
    int photon = -1;
    int chi = -1;
    float jpsiDeltaR = kMissing;
    float jpsiPtRatio = kMissing;
    float photonFreeDeltaR = kMissing;
    float photonPromptDeltaR = kMissing;
    float photonPtRatio = kMissing;
    int directionKind = -1;  // 0=measured conventional P4, 1=free, 2=prompt
    bool directionCompatible = false;
  };
  TruthMatch matchTruth(const pat::CompositeCandidate&, int) const;
  int chiAncestor(int) const;
  int commonChiAncestor(int, int) const;
  bool isRadiativeChiPhoton(int, int) const;

  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniaStatesToken_;
  edm::EDGetTokenT<pat::MuonCollection> oniaMuonsToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> conventionalToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oneProngToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandidatesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  chic::TrackActivityConfig trackActivityConfig_;
  bool doOneProng_;
  bool isMC_;
  const reco::GenParticleCollection* genParticles_ = nullptr;
  const pat::MuonCollection* oniaMuons_ = nullptr;
  std::unordered_map<const reco::Candidate*, int> genIndex_;
  TTree* tree_ = nullptr;

  std::vector<unsigned int> eventRun_, eventLumi_, eventNumber_, eventNPV_, eventNJpsi_, eventNMuon_, eventNChi_,
      eventNPhoton_, eventNElectron_;
  std::vector<float> pvX_, pvY_, pvZ_;
  std::vector<unsigned int> pvCovBegin_, pvCovSize_, pvNTrack_;
  std::vector<float> pvCov_;
  std::vector<float> jpsiPx_, jpsiPy_, jpsiPz_, jpsiE_, jpsiRawPx_, jpsiRawPy_, jpsiRawPz_, jpsiRawE_, jpsiMass_,
      jpsiFitChi2_, jpsiFitNdof_, jpsiPromptChi2_;
  std::vector<int> jpsiFitStatus_, jpsiPVIndex_, jpsiPVCrossCovStatus_, jpsiPromptStatus_, jpsiOniaIndex_,
      jpsiOniaMuPlusIndex_, jpsiOniaMuMinusIndex_;
  std::vector<unsigned int> jpsiMuonBegin_, jpsiMuonSize_, jpsiSVCovBegin_, jpsiSVCovSize_, jpsiPVCovBegin_,
      jpsiPVCovSize_, jpsiPVDeltaCovBegin_, jpsiPVDeltaCovSize_, jpsiPromptCovBegin_, jpsiPromptCovSize_,
      jpsiStateCovBegin_, jpsiStateCovSize_;
  std::vector<float> jpsiSVX_, jpsiSVY_, jpsiSVZ_, jpsiPVX_, jpsiPVY_, jpsiPVZ_, jpsiPVDeltaX_, jpsiPVDeltaY_,
      jpsiPVDeltaZ_, jpsiPromptX_, jpsiPromptY_, jpsiPromptZ_, jpsiSVCov_, jpsiPVCov_, jpsiPVDeltaCov_, jpsiPromptCov_,
      jpsiStateCov_;
  std::vector<float> muonPx_, muonPy_, muonPz_, muonE_, muonVX_, muonVY_, muonVZ_, muonDxy_, muonDz_;
  std::vector<int> muonCharge_, muonJpsiIndex_, muonTrackAvailable_;
  std::vector<unsigned int> muonTrackCovBegin_, muonTrackCovSize_;
  std::vector<float> muonTrackCov_;
  std::vector<int> chiKind_, chiJpsiIndex_, chiPVIndex_, chiPVNTrack_, chiOniaIndex_, chiOniaMuPlusIndex_,
      chiOniaMuMinusIndex_, chiFreeVertexStatus_, chiPromptVertexStatus_;
  std::vector<unsigned int> chiPhotonBegin_, chiPhotonSize_, chiFreeCovBegin_, chiFreeCovSize_, chiPromptCovBegin_,
      chiPromptCovSize_;
  std::vector<float> chiPx_, chiPy_, chiPz_, chiE_, chiRawPx_, chiRawPy_, chiRawPz_, chiRawE_, chiFreeX_, chiFreeY_,
      chiFreeZ_, chiPromptX_, chiPromptY_, chiPromptZ_, chiFreeCov_, chiPromptCov_;
  std::vector<int> chiGenMatchStatus_, chiGenIndex_, chiGenPdgId_, chiGenJpsiIndex_, chiGenPhotonIndex_,
      chiGenPhotonDirectionKind_;
  std::vector<float> chiGenPx_, chiGenPy_, chiGenPz_, chiGenE_, chiGenMass_, chiGenVX_, chiGenVY_, chiGenVZ_;
  std::vector<float> chiGenJpsiPx_, chiGenJpsiPy_, chiGenJpsiPz_, chiGenJpsiE_, chiGenJpsiDeltaR_, chiGenJpsiPtRatio_;
  std::vector<float> chiGenPhotonPx_, chiGenPhotonPy_, chiGenPhotonPz_, chiGenPhotonE_, chiGenPhotonVX_,
      chiGenPhotonVY_, chiGenPhotonVZ_, chiGenPhotonFreeDeltaR_, chiGenPhotonPromptDeltaR_, chiGenPhotonPtRatio_;
  std::vector<int> photonKind_, photonExpectedElectronSize_, photonMissingElectronMeasured_, photonGeometryRoot_,
      photonGeometryStatus_, photonSourceCollection_, photonSourceKey_;
  std::vector<unsigned int> photonElectronBegin_, photonElectronSize_, photonConversionCovBegin_,
      photonConversionCovSize_, photonStateCovBegin_, photonStateCovSize_, photonFreeDirectionCovBegin_,
      photonFreeDirectionCovSize_, photonPromptDirectionCovBegin_, photonPromptDirectionCovSize_;
  std::vector<float> photonPx_, photonPy_, photonPz_, photonE_, photonVisiblePx_, photonVisiblePy_, photonVisiblePz_,
      photonVisibleE_, photonConversionX_, photonConversionY_, photonConversionZ_, photonConversionCov_,
      photonStateCov_, photonFreeDirectionX_, photonFreeDirectionY_, photonFreeDirectionZ_, photonFreeDirectionCov_,
      photonPromptDirectionX_, photonPromptDirectionY_, photonPromptDirectionZ_, photonPromptDirectionCov_,
      photonObservedElectronE_, photonRawThreeBodyMass_, photonRawDeltaM_, photonVisibleMass_, photonVisibleDeltaM_,
      photonPromptCompatibilityChi2_, photonOriginDirectionDelta_, photonOriginDirectionDeltaError_,
      photonOriginAngularSignificance_, photonOriginDeltaChiC0E_, photonOriginDeltaChiC1E_, photonOriginDeltaChiC2E_;
  std::vector<int> photonFreeStatus_, photonPromptStatus_, photonDirectionCrossCovStatus_;
  std::vector<float> photonFreeChiC0E_, photonFreeChiC1E_, photonFreeChiC2E_, photonFreeChiC0Err_, photonFreeChiC1Err_,
      photonFreeChiC2Err_, photonPromptChiC0E_, photonPromptChiC1E_, photonPromptChiC2E_, photonPromptChiC0Err_,
      photonPromptChiC1Err_, photonPromptChiC2Err_;
  std::vector<int> photonGenMatchStatus_, photonGenIndex_;
  std::vector<float> photonGenPx_, photonGenPy_, photonGenPz_, photonGenE_;
  std::vector<float> photonFreeChiC0EResidualToGen_, photonFreeChiC1EResidualToGen_, photonFreeChiC2EResidualToGen_,
      photonPromptChiC0EResidualToGen_, photonPromptChiC1EResidualToGen_, photonPromptChiC2EResidualToGen_;
  std::vector<int> photonFreeChiC0Status_, photonFreeChiC1Status_, photonFreeChiC2Status_, photonPromptChiC0Status_,
      photonPromptChiC1Status_, photonPromptChiC2Status_;
  std::vector<float> electronPx_, electronPy_, electronPz_, electronE_, electronVX_, electronVY_, electronVZ_,
      electronDxy_, electronDz_, electronNormChi2_;
  // electronSource_: 0=packed PF pseudo-track, 1=lost-track pseudo-track,
  // 2=conventional conversion reco::Track.
  std::vector<int> electronCharge_, electronPhotonIndex_, electronSource_, electronSourceKey_, electronHighPurity_,
      electronTrackerHits_, electronPixelHits_, electronTrackerLayers_, electronPixelLayers_, electronLostInnerHits_,
      electronFirstHit_, electronTrackAlgo_, electronOriginalTrackAlgo_, electronTrackAlgoAvailable_;
  std::vector<unsigned int> electronTrackCovBegin_, electronTrackCovSize_;
  std::vector<float> electronTrackCov_;
};

ChiCCovarianceNtuplizer::ChiCCovarianceNtuplizer(const edm::ParameterSet& cfg)
    : oniaStatesToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("oniaStates"))),
      oniaMuonsToken_(mayConsume<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("oniaMuons"))),
      conventionalToken_(
          consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("conventionalCandidates"))),
      oneProngToken_(
          mayConsume<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("oneProngCandidates"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))),
      packedCandidatesToken_(
          mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
      lostTracksToken_(mayConsume<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      genParticlesToken_(mayConsume<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
      trackActivityConfig_{cfg.getParameter<double>("trackMinPt"),
                           cfg.getParameter<double>("trackMaxAbsEta"),
                           cfg.getParameter<double>("trackMaxRelPtError"),
                           cfg.getParameter<double>("trackMaxDzSignificance"),
                           cfg.getParameter<double>("trackMaxDxySignificance")},
      doOneProng_(cfg.getParameter<bool>("doOneProng")),
      isMC_(cfg.getParameter<bool>("isMC")) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> file;
  tree_ = file->make<TTree>("ChiCCovariance", "Nested conventional and one-prong chi_c covariance ntuple");
#define BRANCH(name) tree_->Branch(#name, &name)
  BRANCH(eventRun_); BRANCH(eventLumi_); BRANCH(eventNumber_); BRANCH(eventNPV_); BRANCH(eventNJpsi_); BRANCH(eventNMuon_); BRANCH(eventNChi_); BRANCH(eventNPhoton_); BRANCH(eventNElectron_);
  BRANCH(pvX_); BRANCH(pvY_); BRANCH(pvZ_); BRANCH(pvCovBegin_); BRANCH(pvCovSize_); BRANCH(pvCov_); BRANCH(pvNTrack_);
  BRANCH(jpsiPx_); BRANCH(jpsiPy_); BRANCH(jpsiPz_); BRANCH(jpsiE_); BRANCH(jpsiRawPx_); BRANCH(jpsiRawPy_); BRANCH(jpsiRawPz_); BRANCH(jpsiRawE_); BRANCH(jpsiMass_); BRANCH(jpsiFitChi2_); BRANCH(jpsiFitNdof_); BRANCH(jpsiPromptChi2_); BRANCH(jpsiFitStatus_); BRANCH(jpsiPVIndex_); BRANCH(jpsiPVCrossCovStatus_); BRANCH(jpsiPromptStatus_); BRANCH(jpsiOniaIndex_); BRANCH(jpsiOniaMuPlusIndex_); BRANCH(jpsiOniaMuMinusIndex_); BRANCH(jpsiMuonBegin_); BRANCH(jpsiMuonSize_); BRANCH(jpsiSVCovBegin_); BRANCH(jpsiSVCovSize_); BRANCH(jpsiPVCovBegin_); BRANCH(jpsiPVCovSize_); BRANCH(jpsiPVDeltaCovBegin_); BRANCH(jpsiPVDeltaCovSize_); BRANCH(jpsiPromptCovBegin_); BRANCH(jpsiPromptCovSize_); BRANCH(jpsiStateCovBegin_); BRANCH(jpsiStateCovSize_); BRANCH(jpsiSVX_); BRANCH(jpsiSVY_); BRANCH(jpsiSVZ_); BRANCH(jpsiPVX_); BRANCH(jpsiPVY_); BRANCH(jpsiPVZ_); BRANCH(jpsiPVDeltaX_); BRANCH(jpsiPVDeltaY_); BRANCH(jpsiPVDeltaZ_); BRANCH(jpsiPromptX_); BRANCH(jpsiPromptY_); BRANCH(jpsiPromptZ_); BRANCH(jpsiSVCov_); BRANCH(jpsiPVCov_); BRANCH(jpsiPVDeltaCov_); BRANCH(jpsiPromptCov_); BRANCH(jpsiStateCov_);
  BRANCH(muonPx_); BRANCH(muonPy_); BRANCH(muonPz_); BRANCH(muonE_); BRANCH(muonVX_); BRANCH(muonVY_); BRANCH(muonVZ_); BRANCH(muonDxy_); BRANCH(muonDz_); BRANCH(muonCharge_); BRANCH(muonJpsiIndex_); BRANCH(muonTrackAvailable_); BRANCH(muonTrackCovBegin_); BRANCH(muonTrackCovSize_); BRANCH(muonTrackCov_);
  BRANCH(chiKind_); BRANCH(chiJpsiIndex_); BRANCH(chiPVIndex_); BRANCH(chiPVNTrack_); BRANCH(chiOniaIndex_); BRANCH(chiOniaMuPlusIndex_); BRANCH(chiOniaMuMinusIndex_); BRANCH(chiFreeVertexStatus_); BRANCH(chiPromptVertexStatus_); BRANCH(chiPhotonBegin_); BRANCH(chiPhotonSize_); BRANCH(chiFreeCovBegin_); BRANCH(chiFreeCovSize_); BRANCH(chiPromptCovBegin_); BRANCH(chiPromptCovSize_); BRANCH(chiPx_); BRANCH(chiPy_); BRANCH(chiPz_); BRANCH(chiE_); BRANCH(chiRawPx_); BRANCH(chiRawPy_); BRANCH(chiRawPz_); BRANCH(chiRawE_); BRANCH(chiFreeX_); BRANCH(chiFreeY_); BRANCH(chiFreeZ_); BRANCH(chiPromptX_); BRANCH(chiPromptY_); BRANCH(chiPromptZ_); BRANCH(chiFreeCov_); BRANCH(chiPromptCov_);
  BRANCH(chiGenMatchStatus_); BRANCH(chiGenIndex_); BRANCH(chiGenPdgId_); BRANCH(chiGenJpsiIndex_); BRANCH(chiGenPhotonIndex_); BRANCH(chiGenPhotonDirectionKind_); BRANCH(chiGenPx_); BRANCH(chiGenPy_); BRANCH(chiGenPz_); BRANCH(chiGenE_); BRANCH(chiGenMass_); BRANCH(chiGenVX_); BRANCH(chiGenVY_); BRANCH(chiGenVZ_); BRANCH(chiGenJpsiPx_); BRANCH(chiGenJpsiPy_); BRANCH(chiGenJpsiPz_); BRANCH(chiGenJpsiE_); BRANCH(chiGenJpsiDeltaR_); BRANCH(chiGenJpsiPtRatio_); BRANCH(chiGenPhotonPx_); BRANCH(chiGenPhotonPy_); BRANCH(chiGenPhotonPz_); BRANCH(chiGenPhotonE_); BRANCH(chiGenPhotonVX_); BRANCH(chiGenPhotonVY_); BRANCH(chiGenPhotonVZ_); BRANCH(chiGenPhotonFreeDeltaR_); BRANCH(chiGenPhotonPromptDeltaR_); BRANCH(chiGenPhotonPtRatio_);
  BRANCH(photonKind_); BRANCH(photonExpectedElectronSize_); BRANCH(photonMissingElectronMeasured_); BRANCH(photonGeometryRoot_); BRANCH(photonGeometryStatus_); BRANCH(photonSourceCollection_); BRANCH(photonSourceKey_); BRANCH(photonElectronBegin_); BRANCH(photonElectronSize_); BRANCH(photonConversionCovBegin_); BRANCH(photonConversionCovSize_); BRANCH(photonStateCovBegin_); BRANCH(photonStateCovSize_); BRANCH(photonFreeDirectionCovBegin_); BRANCH(photonFreeDirectionCovSize_); BRANCH(photonPromptDirectionCovBegin_); BRANCH(photonPromptDirectionCovSize_); BRANCH(photonPx_); BRANCH(photonPy_); BRANCH(photonPz_); BRANCH(photonE_); BRANCH(photonVisiblePx_); BRANCH(photonVisiblePy_); BRANCH(photonVisiblePz_); BRANCH(photonVisibleE_); BRANCH(photonConversionX_); BRANCH(photonConversionY_); BRANCH(photonConversionZ_); BRANCH(photonConversionCov_); BRANCH(photonStateCov_); BRANCH(photonFreeDirectionX_); BRANCH(photonFreeDirectionY_); BRANCH(photonFreeDirectionZ_); BRANCH(photonFreeDirectionCov_); BRANCH(photonPromptDirectionX_); BRANCH(photonPromptDirectionY_); BRANCH(photonPromptDirectionZ_); BRANCH(photonPromptDirectionCov_); BRANCH(photonObservedElectronE_); BRANCH(photonRawThreeBodyMass_); BRANCH(photonRawDeltaM_); BRANCH(photonVisibleMass_); BRANCH(photonVisibleDeltaM_); BRANCH(photonPromptCompatibilityChi2_); BRANCH(photonOriginDirectionDelta_); BRANCH(photonOriginDirectionDeltaError_); BRANCH(photonOriginAngularSignificance_); BRANCH(photonOriginDeltaChiC0E_); BRANCH(photonOriginDeltaChiC1E_); BRANCH(photonOriginDeltaChiC2E_); BRANCH(photonFreeStatus_); BRANCH(photonPromptStatus_); BRANCH(photonDirectionCrossCovStatus_); BRANCH(photonFreeChiC0E_); BRANCH(photonFreeChiC1E_); BRANCH(photonFreeChiC2E_); BRANCH(photonFreeChiC0Err_); BRANCH(photonFreeChiC1Err_); BRANCH(photonFreeChiC2Err_); BRANCH(photonPromptChiC0E_); BRANCH(photonPromptChiC1E_); BRANCH(photonPromptChiC2E_); BRANCH(photonPromptChiC0Err_); BRANCH(photonPromptChiC1Err_); BRANCH(photonPromptChiC2Err_); BRANCH(photonFreeChiC0Status_); BRANCH(photonFreeChiC1Status_); BRANCH(photonFreeChiC2Status_); BRANCH(photonPromptChiC0Status_); BRANCH(photonPromptChiC1Status_); BRANCH(photonPromptChiC2Status_);
  BRANCH(photonGenMatchStatus_); BRANCH(photonGenIndex_); BRANCH(photonGenPx_); BRANCH(photonGenPy_); BRANCH(photonGenPz_); BRANCH(photonGenE_); BRANCH(photonFreeChiC0EResidualToGen_); BRANCH(photonFreeChiC1EResidualToGen_); BRANCH(photonFreeChiC2EResidualToGen_); BRANCH(photonPromptChiC0EResidualToGen_); BRANCH(photonPromptChiC1EResidualToGen_); BRANCH(photonPromptChiC2EResidualToGen_);
  BRANCH(electronPx_); BRANCH(electronPy_); BRANCH(electronPz_); BRANCH(electronE_); BRANCH(electronVX_); BRANCH(electronVY_); BRANCH(electronVZ_); BRANCH(electronDxy_); BRANCH(electronDz_); BRANCH(electronNormChi2_); BRANCH(electronCharge_); BRANCH(electronPhotonIndex_); BRANCH(electronSource_); BRANCH(electronSourceKey_); BRANCH(electronHighPurity_); BRANCH(electronTrackerHits_); BRANCH(electronPixelHits_); BRANCH(electronTrackerLayers_); BRANCH(electronPixelLayers_); BRANCH(electronLostInnerHits_); BRANCH(electronFirstHit_); BRANCH(electronTrackAlgo_); BRANCH(electronOriginalTrackAlgo_); BRANCH(electronTrackAlgoAvailable_); BRANCH(electronTrackCovBegin_); BRANCH(electronTrackCovSize_); BRANCH(electronTrackCov_);
#undef BRANCH
}

void ChiCCovarianceNtuplizer::clear() {
#define CLEAR(name) name.clear()
  CLEAR(eventRun_); CLEAR(eventLumi_); CLEAR(eventNumber_); CLEAR(eventNPV_); CLEAR(eventNJpsi_); CLEAR(eventNMuon_); CLEAR(eventNChi_); CLEAR(eventNPhoton_); CLEAR(eventNElectron_);
  CLEAR(pvX_); CLEAR(pvY_); CLEAR(pvZ_); CLEAR(pvCovBegin_); CLEAR(pvCovSize_); CLEAR(pvCov_); CLEAR(pvNTrack_);
  CLEAR(jpsiPx_); CLEAR(jpsiPy_); CLEAR(jpsiPz_); CLEAR(jpsiE_); CLEAR(jpsiRawPx_); CLEAR(jpsiRawPy_); CLEAR(jpsiRawPz_); CLEAR(jpsiRawE_); CLEAR(jpsiMass_); CLEAR(jpsiFitChi2_); CLEAR(jpsiFitNdof_); CLEAR(jpsiPromptChi2_); CLEAR(jpsiFitStatus_); CLEAR(jpsiPVIndex_); CLEAR(jpsiPVCrossCovStatus_); CLEAR(jpsiPromptStatus_); CLEAR(jpsiOniaIndex_); CLEAR(jpsiOniaMuPlusIndex_); CLEAR(jpsiOniaMuMinusIndex_); CLEAR(jpsiMuonBegin_); CLEAR(jpsiMuonSize_); CLEAR(jpsiSVCovBegin_); CLEAR(jpsiSVCovSize_); CLEAR(jpsiPVCovBegin_); CLEAR(jpsiPVCovSize_); CLEAR(jpsiPVDeltaCovBegin_); CLEAR(jpsiPVDeltaCovSize_); CLEAR(jpsiPromptCovBegin_); CLEAR(jpsiPromptCovSize_); CLEAR(jpsiStateCovBegin_); CLEAR(jpsiStateCovSize_); CLEAR(jpsiSVX_); CLEAR(jpsiSVY_); CLEAR(jpsiSVZ_); CLEAR(jpsiPVX_); CLEAR(jpsiPVY_); CLEAR(jpsiPVZ_); CLEAR(jpsiPVDeltaX_); CLEAR(jpsiPVDeltaY_); CLEAR(jpsiPVDeltaZ_); CLEAR(jpsiPromptX_); CLEAR(jpsiPromptY_); CLEAR(jpsiPromptZ_); CLEAR(jpsiSVCov_); CLEAR(jpsiPVCov_); CLEAR(jpsiPVDeltaCov_); CLEAR(jpsiPromptCov_); CLEAR(jpsiStateCov_);
  CLEAR(muonPx_); CLEAR(muonPy_); CLEAR(muonPz_); CLEAR(muonE_); CLEAR(muonVX_); CLEAR(muonVY_); CLEAR(muonVZ_); CLEAR(muonDxy_); CLEAR(muonDz_); CLEAR(muonCharge_); CLEAR(muonJpsiIndex_); CLEAR(muonTrackAvailable_); CLEAR(muonTrackCovBegin_); CLEAR(muonTrackCovSize_); CLEAR(muonTrackCov_);
  CLEAR(chiKind_); CLEAR(chiJpsiIndex_); CLEAR(chiPVIndex_); CLEAR(chiPVNTrack_); CLEAR(chiOniaIndex_); CLEAR(chiOniaMuPlusIndex_); CLEAR(chiOniaMuMinusIndex_); CLEAR(chiFreeVertexStatus_); CLEAR(chiPromptVertexStatus_); CLEAR(chiPhotonBegin_); CLEAR(chiPhotonSize_); CLEAR(chiFreeCovBegin_); CLEAR(chiFreeCovSize_); CLEAR(chiPromptCovBegin_); CLEAR(chiPromptCovSize_); CLEAR(chiPx_); CLEAR(chiPy_); CLEAR(chiPz_); CLEAR(chiE_); CLEAR(chiRawPx_); CLEAR(chiRawPy_); CLEAR(chiRawPz_); CLEAR(chiRawE_); CLEAR(chiFreeX_); CLEAR(chiFreeY_); CLEAR(chiFreeZ_); CLEAR(chiPromptX_); CLEAR(chiPromptY_); CLEAR(chiPromptZ_); CLEAR(chiFreeCov_); CLEAR(chiPromptCov_);
  CLEAR(chiGenMatchStatus_); CLEAR(chiGenIndex_); CLEAR(chiGenPdgId_); CLEAR(chiGenJpsiIndex_); CLEAR(chiGenPhotonIndex_); CLEAR(chiGenPhotonDirectionKind_); CLEAR(chiGenPx_); CLEAR(chiGenPy_); CLEAR(chiGenPz_); CLEAR(chiGenE_); CLEAR(chiGenMass_); CLEAR(chiGenVX_); CLEAR(chiGenVY_); CLEAR(chiGenVZ_); CLEAR(chiGenJpsiPx_); CLEAR(chiGenJpsiPy_); CLEAR(chiGenJpsiPz_); CLEAR(chiGenJpsiE_); CLEAR(chiGenJpsiDeltaR_); CLEAR(chiGenJpsiPtRatio_); CLEAR(chiGenPhotonPx_); CLEAR(chiGenPhotonPy_); CLEAR(chiGenPhotonPz_); CLEAR(chiGenPhotonE_); CLEAR(chiGenPhotonVX_); CLEAR(chiGenPhotonVY_); CLEAR(chiGenPhotonVZ_); CLEAR(chiGenPhotonFreeDeltaR_); CLEAR(chiGenPhotonPromptDeltaR_); CLEAR(chiGenPhotonPtRatio_);
  CLEAR(photonKind_); CLEAR(photonExpectedElectronSize_); CLEAR(photonMissingElectronMeasured_); CLEAR(photonGeometryRoot_); CLEAR(photonGeometryStatus_); CLEAR(photonSourceCollection_); CLEAR(photonSourceKey_); CLEAR(photonElectronBegin_); CLEAR(photonElectronSize_); CLEAR(photonConversionCovBegin_); CLEAR(photonConversionCovSize_); CLEAR(photonStateCovBegin_); CLEAR(photonStateCovSize_); CLEAR(photonFreeDirectionCovBegin_); CLEAR(photonFreeDirectionCovSize_); CLEAR(photonPromptDirectionCovBegin_); CLEAR(photonPromptDirectionCovSize_); CLEAR(photonPx_); CLEAR(photonPy_); CLEAR(photonPz_); CLEAR(photonE_); CLEAR(photonVisiblePx_); CLEAR(photonVisiblePy_); CLEAR(photonVisiblePz_); CLEAR(photonVisibleE_); CLEAR(photonConversionX_); CLEAR(photonConversionY_); CLEAR(photonConversionZ_); CLEAR(photonConversionCov_); CLEAR(photonStateCov_); CLEAR(photonFreeDirectionX_); CLEAR(photonFreeDirectionY_); CLEAR(photonFreeDirectionZ_); CLEAR(photonFreeDirectionCov_); CLEAR(photonPromptDirectionX_); CLEAR(photonPromptDirectionY_); CLEAR(photonPromptDirectionZ_); CLEAR(photonPromptDirectionCov_); CLEAR(photonObservedElectronE_); CLEAR(photonRawThreeBodyMass_); CLEAR(photonRawDeltaM_); CLEAR(photonVisibleMass_); CLEAR(photonVisibleDeltaM_); CLEAR(photonPromptCompatibilityChi2_); CLEAR(photonOriginDirectionDelta_); CLEAR(photonOriginDirectionDeltaError_); CLEAR(photonOriginAngularSignificance_); CLEAR(photonOriginDeltaChiC0E_); CLEAR(photonOriginDeltaChiC1E_); CLEAR(photonOriginDeltaChiC2E_); CLEAR(photonFreeStatus_); CLEAR(photonPromptStatus_); CLEAR(photonDirectionCrossCovStatus_); CLEAR(photonFreeChiC0E_); CLEAR(photonFreeChiC1E_); CLEAR(photonFreeChiC2E_); CLEAR(photonFreeChiC0Err_); CLEAR(photonFreeChiC1Err_); CLEAR(photonFreeChiC2Err_); CLEAR(photonPromptChiC0E_); CLEAR(photonPromptChiC1E_); CLEAR(photonPromptChiC2E_); CLEAR(photonPromptChiC0Err_); CLEAR(photonPromptChiC1Err_); CLEAR(photonPromptChiC2Err_); CLEAR(photonFreeChiC0Status_); CLEAR(photonFreeChiC1Status_); CLEAR(photonFreeChiC2Status_); CLEAR(photonPromptChiC0Status_); CLEAR(photonPromptChiC1Status_); CLEAR(photonPromptChiC2Status_);
  CLEAR(photonGenMatchStatus_); CLEAR(photonGenIndex_); CLEAR(photonGenPx_); CLEAR(photonGenPy_); CLEAR(photonGenPz_); CLEAR(photonGenE_); CLEAR(photonFreeChiC0EResidualToGen_); CLEAR(photonFreeChiC1EResidualToGen_); CLEAR(photonFreeChiC2EResidualToGen_); CLEAR(photonPromptChiC0EResidualToGen_); CLEAR(photonPromptChiC1EResidualToGen_); CLEAR(photonPromptChiC2EResidualToGen_);
  CLEAR(electronPx_); CLEAR(electronPy_); CLEAR(electronPz_); CLEAR(electronE_); CLEAR(electronVX_); CLEAR(electronVY_); CLEAR(electronVZ_); CLEAR(electronDxy_); CLEAR(electronDz_); CLEAR(electronNormChi2_); CLEAR(electronCharge_); CLEAR(electronPhotonIndex_); CLEAR(electronSource_); CLEAR(electronSourceKey_); CLEAR(electronHighPurity_); CLEAR(electronTrackerHits_); CLEAR(electronPixelHits_); CLEAR(electronTrackerLayers_); CLEAR(electronPixelLayers_); CLEAR(electronLostInnerHits_); CLEAR(electronFirstHit_); CLEAR(electronTrackAlgo_); CLEAR(electronOriginalTrackAlgo_); CLEAR(electronTrackAlgoAvailable_); CLEAR(electronTrackCovBegin_); CLEAR(electronTrackCovSize_); CLEAR(electronTrackCov_);
#undef CLEAR
}

void ChiCCovarianceNtuplizer::appendTrackCovariance(const reco::Track& track,
                                                    std::vector<unsigned int>& covBegin,
                                                    std::vector<unsigned int>& covSize,
                                                    std::vector<float>& covariances) {
  covBegin.push_back(covariances.size());
  appendLowerTriangle(track.covariance(), 5, covariances);
  covSize.push_back(15);
}

void ChiCCovarianceNtuplizer::appendMuon(const pat::Muon& muon, unsigned int jpsiIndex) {
  muonPx_.push_back(muon.px());
  muonPy_.push_back(muon.py());
  muonPz_.push_back(muon.pz());
  muonE_.push_back(muon.energy());
  muonVX_.push_back(muon.vx());
  muonVY_.push_back(muon.vy());
  muonVZ_.push_back(muon.vz());
  muonCharge_.push_back(muon.charge());
  muonJpsiIndex_.push_back(jpsiIndex);
  if (muon.innerTrack().isNonnull()) {
    const auto& track = *muon.innerTrack();
    muonTrackAvailable_.push_back(1);
    muonDxy_.push_back(track.dxy());
    muonDz_.push_back(track.dz());
    appendTrackCovariance(track, muonTrackCovBegin_, muonTrackCovSize_, muonTrackCov_);
  } else {
    muonTrackAvailable_.push_back(0);
    muonDxy_.push_back(kMissing);
    muonDz_.push_back(kMissing);
    muonTrackCovBegin_.push_back(muonTrackCov_.size());
    muonTrackCovSize_.push_back(0);
  }
}

int ChiCCovarianceNtuplizer::oniaMuonIndex(const pat::Muon& target) const {
  if (oniaMuons_ == nullptr)
    return -1;

  // HiOnia's Reco_QQ_mupl_idx/Reco_QQ_mumi_idx are produced through
  // mapMuonMomToIndex_, whose key is the pT reduced to [0.01, 0.1] and
  // truncated at 1.e-7. The map overwrites equal keys, so the last matching
  // muon wins. Reproduce that convention here: these branches are links to
  // the HiOnia tree, not independent physical-object matches.
  auto oniaMomentumKey = [](float pt) -> long int {
    float reduced = std::abs(pt);
    if (!std::isfinite(reduced))
      return -1L;
    while (reduced > 0.1f)
      reduced /= 10.f;
    return static_cast<long int>(10000000.f * reduced);
  };
  const long int targetKey = oniaMomentumKey(target.pt());
  if (targetKey < 0)
    return -1;
  int matchingIndex = -1;
  for (size_t index = 0; index < oniaMuons_->size(); ++index) {
    const auto& source = oniaMuons_->at(index);
    if (oniaMomentumKey(source.pt()) == targetKey)
      matchingIndex = static_cast<int>(index);
  }
  return matchingIndex;
}

int ChiCCovarianceNtuplizer::appendOniaLinks(int oniaIdxVal) {
  chiOniaIndex_.push_back(oniaIdxVal);
  const auto found = std::find(jpsiOniaIndex_.begin(), jpsiOniaIndex_.end(), oniaIdxVal);
  if (found != jpsiOniaIndex_.end()) {
    const size_t stateIndex = std::distance(jpsiOniaIndex_.begin(), found);
    chiOniaMuPlusIndex_.push_back(jpsiOniaMuPlusIndex_[stateIndex]);
    chiOniaMuMinusIndex_.push_back(jpsiOniaMuMinusIndex_[stateIndex]);
    return static_cast<int>(stateIndex);
  } else {
    chiOniaMuPlusIndex_.push_back(-1);
    chiOniaMuMinusIndex_.push_back(-1);
    return -1;
  }
}

void ChiCCovarianceNtuplizer::appendElectronFromTrack(const reco::Track& track,
                                                      unsigned int photonIndex,
                                                      int source,
                                                      int sourceKey) {
  electronPx_.push_back(track.px());
  electronPy_.push_back(track.py());
  electronPz_.push_back(track.pz());
  electronE_.push_back(electronEnergy(track.px(), track.py(), track.pz()));
  electronVX_.push_back(track.vx());
  electronVY_.push_back(track.vy());
  electronVZ_.push_back(track.vz());
  electronDxy_.push_back(track.dxy());
  electronDz_.push_back(track.dz());
  electronNormChi2_.push_back(track.normalizedChi2());
  electronCharge_.push_back(track.charge());
  electronPhotonIndex_.push_back(photonIndex);
  electronSource_.push_back(source);
  electronSourceKey_.push_back(sourceKey);
  electronHighPurity_.push_back(track.quality(reco::TrackBase::highPurity));
  electronTrackerHits_.push_back(track.numberOfValidHits());
  electronPixelHits_.push_back(track.hitPattern().numberOfValidPixelHits());
  electronTrackerLayers_.push_back(track.hitPattern().trackerLayersWithMeasurement());
  electronPixelLayers_.push_back(track.hitPattern().pixelLayersWithMeasurement());
  electronLostInnerHits_.push_back(-1);
  electronFirstHit_.push_back(-1);
  electronTrackAlgo_.push_back(track.algo());
  electronOriginalTrackAlgo_.push_back(track.originalAlgo());
  electronTrackAlgoAvailable_.push_back(1);
  appendTrackCovariance(track, electronTrackCovBegin_, electronTrackCovSize_, electronTrackCov_);
}

void ChiCCovarianceNtuplizer::appendElectronFromOneProng(const pat::CompositeCandidate& leg, unsigned int photonIndex) {
  electronPx_.push_back(leg.px());
  electronPy_.push_back(leg.py());
  electronPz_.push_back(leg.pz());
  electronE_.push_back(leg.energy());
  electronVX_.push_back(leg.vx());
  electronVY_.push_back(leg.vy());
  electronVZ_.push_back(leg.vz());
  electronDxy_.push_back(value(leg, "trackDxy"));
  electronDz_.push_back(value(leg, "trackDz"));
  electronNormChi2_.push_back(value(leg, "trackNormalizedChi2"));
  electronCharge_.push_back(leg.charge());
  electronPhotonIndex_.push_back(photonIndex);
  electronSource_.push_back(intValue(leg, "sourceCollection"));
  electronSourceKey_.push_back(intValue(leg, "sourceKey"));
  electronHighPurity_.push_back(intValue(leg, "trackHighPurity"));
  electronTrackerHits_.push_back(intValue(leg, "nTrackerHits"));
  electronPixelHits_.push_back(intValue(leg, "nPixelHits"));
  electronTrackerLayers_.push_back(intValue(leg, "nTrackerLayers"));
  electronPixelLayers_.push_back(intValue(leg, "nPixelLayers"));
  electronLostInnerHits_.push_back(intValue(leg, "lostInnerHits"));
  electronFirstHit_.push_back(intValue(leg, "firstHit"));
  electronTrackAlgo_.push_back(intValue(leg, "trackAlgo"));
  electronOriginalTrackAlgo_.push_back(intValue(leg, "trackOriginalAlgo"));
  electronTrackAlgoAvailable_.push_back(intValue(leg, "trackAlgoAvailable", 0));
  electronTrackCovBegin_.push_back(electronTrackCov_.size());
  for (unsigned int row = 0; row != 5; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      electronTrackCov_.push_back(value(leg, "trackCov" + std::to_string(row) + std::to_string(column)));
  electronTrackCovSize_.push_back(15);
}

int ChiCCovarianceNtuplizer::oniaIndex(const reco::Candidate* candidate) {
  const auto* onia = dynamic_cast<const pat::CompositeCandidate*>(candidate);
  return onia == nullptr ? -1 : intValue(*onia, "oniaStateIndex");
}

int ChiCCovarianceNtuplizer::chiAncestor(int index) const {
  if (genParticles_ == nullptr || index < 0 || static_cast<size_t>(index) >= genParticles_->size())
    return -1;
  std::vector<const reco::Candidate*> pending{&genParticles_->at(index)};
  std::unordered_set<const reco::Candidate*> visited;
  while (!pending.empty()) {
    const reco::Candidate* current = pending.back();
    pending.pop_back();
    if (current == nullptr || !visited.insert(current).second)
      continue;
    const int absId = std::abs(current->pdgId());
    if (absId == 10441 || absId == 20443 || absId == 445) {
      const auto found = genIndex_.find(current);
      if (found != genIndex_.end())
        return found->second;
    }
    for (size_t mother = 0; mother < current->numberOfMothers(); ++mother)
      pending.push_back(current->mother(mother));
  }
  return -1;
}

int ChiCCovarianceNtuplizer::commonChiAncestor(int jpsi, int photon) const {
  const int first = chiAncestor(jpsi);
  const int second = chiAncestor(photon);
  return first >= 0 && first == second ? first : -1;
}

bool ChiCCovarianceNtuplizer::isRadiativeChiPhoton(int photon, int chi) const {
  if (genParticles_ == nullptr || photon < 0 || chi < 0 ||
      static_cast<size_t>(photon) >= genParticles_->size() ||
      static_cast<size_t>(chi) >= genParticles_->size())
    return false;
  struct Node {
    const reco::Candidate* particle;
    bool passedJpsi;
  };
  std::vector<Node> pending{{&genParticles_->at(photon), false}};
  std::unordered_set<const reco::Candidate*> visitedWithoutJpsi;
  std::unordered_set<const reco::Candidate*> visitedWithJpsi;
  while (!pending.empty()) {
    const auto node = pending.back();
    pending.pop_back();
    if (node.particle == nullptr)
      continue;
    auto& visited = node.passedJpsi ? visitedWithJpsi : visitedWithoutJpsi;
    if (!visited.insert(node.particle).second)
      continue;
    const auto found = genIndex_.find(node.particle);
    if (found != genIndex_.end() && found->second == chi)
      return !node.passedJpsi;
    const bool passedJpsi = node.passedJpsi || std::abs(node.particle->pdgId()) == 443;
    for (size_t mother = 0; mother < node.particle->numberOfMothers(); ++mother)
      pending.push_back({node.particle->mother(mother), passedJpsi});
  }
  return false;
}

ChiCCovarianceNtuplizer::TruthMatch ChiCCovarianceNtuplizer::matchTruth(
    const pat::CompositeCandidate& candidate, int kind) const {
  TruthMatch best;
  if (!isMC_ || genParticles_ == nullptr)
    return best;
  const auto* jpsi = candidate.daughter(kind == 0 ? "dimuon" : "jpsi");
  const auto* measuredPhoton = kind == 0 ? candidate.daughter("photon") : nullptr;
  if (jpsi == nullptr || (kind == 0 && measuredPhoton == nullptr))
    return best;

  struct JpsiMatch {
    int index;
    float deltaR;
    float ptRatio;
  };
  struct PhotonMatch {
    int index;
    float freeDeltaR;
    float promptDeltaR;
    float ptRatio;
    int directionKind;
    float matchDeltaR;
  };
  std::vector<JpsiMatch> jpsiMatches;
  std::vector<PhotonMatch> photonMatches;

  float freeEta = kMissing, freePhi = kMissing, promptEta = kMissing, promptPhi = kMissing;
  if (kind == 1) {
    const float freeX = value(candidate, "freePhotonDirectionX");
    const float freeY = value(candidate, "freePhotonDirectionY");
    const float freeZ = value(candidate, "freePhotonDirectionZ");
    const float freePt = std::hypot(freeX, freeY);
    if (freePt > 0.f) {
      freeEta = std::asinh(freeZ / freePt);
      freePhi = std::atan2(freeY, freeX);
    }
    const float promptX = value(candidate, "promptPhotonDirectionX");
    const float promptY = value(candidate, "promptPhotonDirectionY");
    const float promptZ = value(candidate, "promptPhotonDirectionZ");
    const float promptPt = std::hypot(promptX, promptY);
    if (promptPt > 0.f) {
      promptEta = std::asinh(promptZ / promptPt);
      promptPhi = std::atan2(promptY, promptX);
    }
  }

  for (size_t index = 0; index < genParticles_->size(); ++index) {
    const auto& gen = genParticles_->at(index);
    if (chiAncestor(index) < 0)
      continue;
    if (std::abs(gen.pdgId()) == 443 && gen.pt() > 0.f) {
      const float deltaR = reco::deltaR(*jpsi, gen);
      const float ptRatio = jpsi->pt() / gen.pt();
      if (deltaR < .05f && ptRatio > .2f && ptRatio < 5.f)
        jpsiMatches.push_back({static_cast<int>(index), deltaR, ptRatio});
    }
    if (gen.pdgId() != 22)
      continue;
    const int photonChi = chiAncestor(index);
    if (photonChi < 0 || !isRadiativeChiPhoton(index, photonChi))
      continue;
    if (kind == 0) {
      const float deltaR = reco::deltaR(*measuredPhoton, gen);
      const float ptRatio = gen.pt() > 0.f ? measuredPhoton->pt() / gen.pt() : kMissing;
      if (deltaR < .05f && ptRatio > .2f && ptRatio < 5.f)
        photonMatches.push_back({static_cast<int>(index), deltaR, kMissing, ptRatio, 0, deltaR});
    } else {
      const float freeDeltaR =
          std::isfinite(freeEta) ? reco::deltaR(freeEta, freePhi, float(gen.eta()), float(gen.phi())) : kMissing;
      const float promptDeltaR =
          std::isfinite(promptEta) ? reco::deltaR(promptEta, promptPhi, float(gen.eta()), float(gen.phi())) : kMissing;
      const bool useFree = std::isfinite(freeDeltaR) && (!std::isfinite(promptDeltaR) || freeDeltaR <= promptDeltaR);
      const float matchDeltaR = useFree ? freeDeltaR : promptDeltaR;
      if (std::isfinite(matchDeltaR))
        photonMatches.push_back(
            {static_cast<int>(index), freeDeltaR, promptDeltaR, kMissing, useFree ? 1 : 2, matchDeltaR});
    }
  }

  float bestScore = std::numeric_limits<float>::max();
  for (const auto& jpsiMatch : jpsiMatches) {
    for (const auto& photonMatch : photonMatches) {
      const int chi = commonChiAncestor(jpsiMatch.index, photonMatch.index);
      if (chi < 0)
        continue;
      const float score =
          jpsiMatch.deltaR + (kind == 0 ? photonMatch.matchDeltaR : 1.e-3f * std::min(photonMatch.matchDeltaR, 10.f));
      if (score >= bestScore)
        continue;
      bestScore = score;
      best.jpsi = jpsiMatch.index;
      best.photon = photonMatch.index;
      best.chi = chi;
      best.jpsiDeltaR = jpsiMatch.deltaR;
      best.jpsiPtRatio = jpsiMatch.ptRatio;
      best.photonFreeDeltaR = photonMatch.freeDeltaR;
      best.photonPromptDeltaR = photonMatch.promptDeltaR;
      best.photonPtRatio = photonMatch.ptRatio;
      best.directionKind = photonMatch.directionKind;
      best.directionCompatible = photonMatch.matchDeltaR < .05f;
    }
  }

  if (best.jpsi < 0 && !jpsiMatches.empty()) {
    const auto closest =
        std::min_element(jpsiMatches.begin(), jpsiMatches.end(), [](const auto& first, const auto& second) {
          return first.deltaR < second.deltaR;
        });
    best.jpsi = closest->index;
    best.jpsiDeltaR = closest->deltaR;
    best.jpsiPtRatio = closest->ptRatio;
    best.chi = chiAncestor(closest->index);
  }
  return best;
}

void ChiCCovarianceNtuplizer::appendTruth(const pat::CompositeCandidate& candidate, int kind) {
  const TruthMatch match = matchTruth(candidate, kind);
  const bool hasChi = match.chi >= 0;
  const bool hasJpsi = match.jpsi >= 0;
  const bool hasPhoton = match.photon >= 0;
  const auto* chi = hasChi ? &genParticles_->at(match.chi) : nullptr;
  const auto* jpsi = hasJpsi ? &genParticles_->at(match.jpsi) : nullptr;
  const auto* photon = hasPhoton ? &genParticles_->at(match.photon) : nullptr;
  const int status =
      (hasJpsi ? 1 : 0) | (hasPhoton ? 2 : 0) | (hasChi ? 4 : 0) | (kind == 1 && match.directionCompatible ? 8 : 0);

  chiGenMatchStatus_.push_back(status);
  chiGenIndex_.push_back(hasChi ? match.chi : -1);
  chiGenPdgId_.push_back(chi ? chi->pdgId() : 0);
  chiGenJpsiIndex_.push_back(hasJpsi ? match.jpsi : -1);
  chiGenPhotonIndex_.push_back(hasPhoton ? match.photon : -1);
  chiGenPhotonDirectionKind_.push_back(hasPhoton ? match.directionKind : -1);
  for (auto pair : {std::make_pair(&chiGenPx_, chi ? float(chi->px()) : kMissing),
                    std::make_pair(&chiGenPy_, chi ? float(chi->py()) : kMissing),
                    std::make_pair(&chiGenPz_, chi ? float(chi->pz()) : kMissing),
                    std::make_pair(&chiGenE_, chi ? float(chi->energy()) : kMissing),
                    std::make_pair(&chiGenMass_, chi ? float(chi->mass()) : kMissing),
                    std::make_pair(&chiGenVX_, chi ? float(chi->vx()) : kMissing),
                    std::make_pair(&chiGenVY_, chi ? float(chi->vy()) : kMissing),
                    std::make_pair(&chiGenVZ_, chi ? float(chi->vz()) : kMissing),
                    std::make_pair(&chiGenJpsiPx_, jpsi ? float(jpsi->px()) : kMissing),
                    std::make_pair(&chiGenJpsiPy_, jpsi ? float(jpsi->py()) : kMissing),
                    std::make_pair(&chiGenJpsiPz_, jpsi ? float(jpsi->pz()) : kMissing),
                    std::make_pair(&chiGenJpsiE_, jpsi ? float(jpsi->energy()) : kMissing),
                    std::make_pair(&chiGenPhotonPx_, photon ? float(photon->px()) : kMissing),
                    std::make_pair(&chiGenPhotonPy_, photon ? float(photon->py()) : kMissing),
                    std::make_pair(&chiGenPhotonPz_, photon ? float(photon->pz()) : kMissing),
                    std::make_pair(&chiGenPhotonE_, photon ? float(photon->energy()) : kMissing),
                    std::make_pair(&chiGenPhotonVX_, photon ? float(photon->vx()) : kMissing),
                    std::make_pair(&chiGenPhotonVY_, photon ? float(photon->vy()) : kMissing),
                    std::make_pair(&chiGenPhotonVZ_, photon ? float(photon->vz()) : kMissing)})
    pair.first->push_back(pair.second);
  chiGenJpsiDeltaR_.push_back(hasJpsi ? match.jpsiDeltaR : kMissing);
  chiGenJpsiPtRatio_.push_back(hasJpsi ? match.jpsiPtRatio : kMissing);
  chiGenPhotonFreeDeltaR_.push_back(hasPhoton ? match.photonFreeDeltaR : kMissing);
  chiGenPhotonPromptDeltaR_.push_back(hasPhoton ? match.photonPromptDeltaR : kMissing);
  chiGenPhotonPtRatio_.push_back(hasPhoton ? match.photonPtRatio : kMissing);

  photonGenMatchStatus_.push_back(status);
  photonGenIndex_.push_back(hasPhoton ? match.photon : -1);
  photonGenPx_.push_back(photon ? photon->px() : kMissing);
  photonGenPy_.push_back(photon ? photon->py() : kMissing);
  photonGenPz_.push_back(photon ? photon->pz() : kMissing);
  photonGenE_.push_back(photon ? photon->energy() : kMissing);
  const float genEnergy = photon ? photon->energy() : kMissing;
  auto closure = [genEnergy](float estimate) {
    return std::isfinite(genEnergy) && std::isfinite(estimate) ? estimate - genEnergy : kMissing;
  };
  photonFreeChiC0EResidualToGen_.push_back(closure(value(candidate, "freeChiC0PhotonEnergy")));
  photonFreeChiC1EResidualToGen_.push_back(closure(value(candidate, "freeChiC1PhotonEnergy")));
  photonFreeChiC2EResidualToGen_.push_back(closure(value(candidate, "freeChiC2PhotonEnergy")));
  photonPromptChiC0EResidualToGen_.push_back(closure(value(candidate, "promptChiC0PhotonEnergy")));
  photonPromptChiC1EResidualToGen_.push_back(closure(value(candidate, "promptChiC1PhotonEnergy")));
  photonPromptChiC2EResidualToGen_.push_back(closure(value(candidate, "promptChiC2PhotonEnergy")));
}

void ChiCCovarianceNtuplizer::appendConventional(const pat::CompositeCandidate& chi) {
  const auto* jpsi = chi.daughter("dimuon");
  const auto* photon = dynamic_cast<const pat::CompositeCandidate*>(chi.daughter("photon"));
  const auto* state = dynamic_cast<const pat::CompositeCandidate*>(jpsi);
  if (jpsi == nullptr || photon == nullptr || state == nullptr)
    return;
  const unsigned int photonIndex = photonKind_.size(), electronBegin = electronPx_.size();
  const int linkedOnia = oniaIndex(jpsi);
  chiKind_.push_back(0);
  chiJpsiIndex_.push_back(appendOniaLinks(linkedOnia));
  chiPVIndex_.push_back(intValue(*state, "pvIndex"));
  chiPhotonBegin_.push_back(photonIndex);
  chiPhotonSize_.push_back(1);
  chiPx_.push_back(chi.px());
  chiPy_.push_back(chi.py());
  chiPz_.push_back(chi.pz());
  chiE_.push_back(chi.energy());
  chiRawPx_.push_back(chi.px());
  chiRawPy_.push_back(chi.py());
  chiRawPz_.push_back(chi.pz());
  chiRawE_.push_back(chi.energy());
  chiFreeX_.push_back(jpsi->vx());
  chiFreeY_.push_back(jpsi->vy());
  chiFreeZ_.push_back(jpsi->vz());
  chiPromptX_.push_back(value(*state, "promptOriginX"));
  chiPromptY_.push_back(value(*state, "promptOriginY"));
  chiPromptZ_.push_back(value(*state, "promptOriginZ"));
  chiFreeVertexStatus_.push_back(1);
  chiPromptVertexStatus_.push_back(intValue(*state, "promptOriginStatus"));
  chiFreeCovBegin_.push_back(chiFreeCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      chiFreeCov_.push_back(value(*state, "jpsiVertexCov" + std::to_string(row) + std::to_string(column)));
  chiFreeCovSize_.push_back(6);
  chiPromptCovBegin_.push_back(chiPromptCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      chiPromptCov_.push_back(value(*state, "promptOriginCov" + std::to_string(row) + std::to_string(column)));
  chiPromptCovSize_.push_back(6);
  photonKind_.push_back(0);
  photonExpectedElectronSize_.push_back(2);
  photonMissingElectronMeasured_.push_back(0);
  photonGeometryRoot_.push_back(-1);
  photonGeometryStatus_.push_back(intValue(chi, "conventionalPhotonStateStatus"));
  photonSourceCollection_.push_back(-1);
  photonSourceKey_.push_back(-1);
  photonElectronBegin_.push_back(electronBegin);
  photonPx_.push_back(value(chi, "conventionalPhotonPx", photon->px()));
  photonPy_.push_back(value(chi, "conventionalPhotonPy", photon->py()));
  photonPz_.push_back(value(chi, "conventionalPhotonPz", photon->pz()));
  photonE_.push_back(value(chi, "conventionalPhotonE", photon->energy()));
  photonVisiblePx_.push_back(value(chi, "conventionalPhotonPx", photon->px()));
  photonVisiblePy_.push_back(value(chi, "conventionalPhotonPy", photon->py()));
  photonVisiblePz_.push_back(value(chi, "conventionalPhotonPz", photon->pz()));
  photonVisibleE_.push_back(value(chi, "conventionalPhotonE", photon->energy()));
  photonConversionX_.push_back(value(chi, "conventionalConversionX", photon->vx()));
  photonConversionY_.push_back(value(chi, "conventionalConversionY", photon->vy()));
  photonConversionZ_.push_back(value(chi, "conventionalConversionZ", photon->vz()));
  photonConversionCovBegin_.push_back(photonConversionCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      photonConversionCov_.push_back(
          value(chi, "conventionalConversionCov" + std::to_string(row) + std::to_string(column)));
  photonConversionCovSize_.push_back(6);
  photonStateCovBegin_.push_back(photonStateCov_.size());
  for (unsigned int row = 0; row != 7; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      photonStateCov_.push_back(
          value(chi, "conventionalPhotonStateCov" + std::to_string(row) + std::to_string(column)));
  photonStateCovSize_.push_back(28);
  photonFreeDirectionX_.push_back(kMissing);
  photonFreeDirectionY_.push_back(kMissing);
  photonFreeDirectionZ_.push_back(kMissing);
  photonFreeDirectionCovBegin_.push_back(photonFreeDirectionCov_.size());
  photonFreeDirectionCovSize_.push_back(0);
  photonPromptDirectionX_.push_back(kMissing);
  photonPromptDirectionY_.push_back(kMissing);
  photonPromptDirectionZ_.push_back(kMissing);
  photonPromptDirectionCovBegin_.push_back(photonPromptDirectionCov_.size());
  photonPromptDirectionCovSize_.push_back(0);
  photonObservedElectronE_.push_back(kMissing);
  photonRawThreeBodyMass_.push_back(chi.mass());
  photonRawDeltaM_.push_back(chi.mass() - jpsi->mass());
  photonVisibleMass_.push_back(chi.mass());
  photonVisibleDeltaM_.push_back(chi.mass() - jpsi->mass());
  photonPromptCompatibilityChi2_.push_back(value(*state, "promptCompatibilityChi2"));
  photonOriginDirectionDelta_.push_back(kMissing);
  photonOriginDirectionDeltaError_.push_back(kMissing);
  photonOriginAngularSignificance_.push_back(kMissing);
  photonOriginDeltaChiC0E_.push_back(kMissing);
  photonOriginDeltaChiC1E_.push_back(kMissing);
  photonOriginDeltaChiC2E_.push_back(kMissing);
  photonFreeStatus_.push_back(0);
  photonPromptStatus_.push_back(0);
  photonDirectionCrossCovStatus_.push_back(-1);
  const reco::Track* track0 = photon->hasUserData("track0") ? photon->userData<reco::Track>("track0") : nullptr;
  const reco::Track* track1 = photon->hasUserData("track1") ? photon->userData<reco::Track>("track1") : nullptr;
  if (track0)
    appendElectronFromTrack(*track0, photonIndex, 2, 0);
  if (track1)
    appendElectronFromTrack(*track1, photonIndex, 2, 1);
  photonElectronSize_.push_back(electronPx_.size() - electronBegin);
  for (auto* destination : {&photonFreeChiC0E_,
                            &photonFreeChiC1E_,
                            &photonFreeChiC2E_,
                            &photonFreeChiC0Err_,
                            &photonFreeChiC1Err_,
                            &photonFreeChiC2Err_,
                            &photonPromptChiC0E_,
                            &photonPromptChiC1E_,
                            &photonPromptChiC2E_,
                            &photonPromptChiC0Err_,
                            &photonPromptChiC1Err_,
                            &photonPromptChiC2Err_})
    destination->push_back(kMissing);
  for (auto* destination : {&photonFreeChiC0Status_,
                            &photonFreeChiC1Status_,
                            &photonFreeChiC2Status_,
                            &photonPromptChiC0Status_,
                            &photonPromptChiC1Status_,
                            &photonPromptChiC2Status_})
    destination->push_back(0);
}

void ChiCCovarianceNtuplizer::appendOneProng(const pat::CompositeCandidate& chi) {
  const auto* jpsi = dynamic_cast<const pat::CompositeCandidate*>(chi.daughter("jpsi"));
  const auto* leg = dynamic_cast<const pat::CompositeCandidate*>(chi.daughter("leg"));
  if (jpsi == nullptr || leg == nullptr)
    return;
  const unsigned int photonIndex = photonKind_.size(), electronBegin = electronPx_.size();
  const int linkedOnia = intValue(*jpsi, "oniaStateIndex");
  chiKind_.push_back(1);
  chiJpsiIndex_.push_back(appendOniaLinks(linkedOnia));
  chiPVIndex_.push_back(intValue(*jpsi, "pvIndex"));
  chiPhotonBegin_.push_back(photonIndex);
  chiPhotonSize_.push_back(1);
  chiPx_.push_back(chi.px());
  chiPy_.push_back(chi.py());
  chiPz_.push_back(chi.pz());
  chiE_.push_back(chi.energy());
  chiRawPx_.push_back(jpsi->px() + leg->px());
  chiRawPy_.push_back(jpsi->py() + leg->py());
  chiRawPz_.push_back(jpsi->pz() + leg->pz());
  chiRawE_.push_back(jpsi->energy() + leg->energy());
  chiFreeX_.push_back(jpsi->vx());
  chiFreeY_.push_back(jpsi->vy());
  chiFreeZ_.push_back(jpsi->vz());
  chiPromptX_.push_back(value(*jpsi, "promptOriginX"));
  chiPromptY_.push_back(value(*jpsi, "promptOriginY"));
  chiPromptZ_.push_back(value(*jpsi, "promptOriginZ"));
  chiFreeVertexStatus_.push_back(intValue(chi, "freeHypothesisStatus"));
  chiPromptVertexStatus_.push_back(intValue(chi, "promptHypothesisStatus"));
  chiFreeCovBegin_.push_back(chiFreeCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      chiFreeCov_.push_back(value(*jpsi, "jpsiVertexCov" + std::to_string(row) + std::to_string(column)));
  chiFreeCovSize_.push_back(6);
  chiPromptCovBegin_.push_back(chiPromptCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      chiPromptCov_.push_back(value(*jpsi, "promptOriginCov" + std::to_string(row) + std::to_string(column)));
  chiPromptCovSize_.push_back(6);
  photonKind_.push_back(1);
  photonExpectedElectronSize_.push_back(2);
  photonMissingElectronMeasured_.push_back(0);
  photonGeometryRoot_.push_back(intValue(chi, "geometryRoot"));
  photonGeometryStatus_.push_back(intValue(chi, "geometryStatus"));
  photonSourceCollection_.push_back(intValue(*leg, "sourceCollection"));
  photonSourceKey_.push_back(intValue(*leg, "sourceKey"));
  photonElectronBegin_.push_back(electronBegin);
  photonPx_.push_back(kMissing);
  photonPy_.push_back(kMissing);
  photonPz_.push_back(kMissing);
  photonE_.push_back(kMissing);
  photonVisiblePx_.push_back(chi.px() - jpsi->px());
  photonVisiblePy_.push_back(chi.py() - jpsi->py());
  photonVisiblePz_.push_back(chi.pz() - jpsi->pz());
  photonVisibleE_.push_back(chi.energy() - jpsi->energy());
  photonConversionX_.push_back(chi.vx());
  photonConversionY_.push_back(chi.vy());
  photonConversionZ_.push_back(chi.vz());
  photonConversionCovBegin_.push_back(photonConversionCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      photonConversionCov_.push_back(value(chi, "conversionCov" + std::to_string(row) + std::to_string(column)));
  photonConversionCovSize_.push_back(6);
  photonStateCovBegin_.push_back(photonStateCov_.size());
  photonStateCovSize_.push_back(0);
  photonFreeDirectionX_.push_back(value(chi, "freePhotonDirectionX"));
  photonFreeDirectionY_.push_back(value(chi, "freePhotonDirectionY"));
  photonFreeDirectionZ_.push_back(value(chi, "freePhotonDirectionZ"));
  photonFreeDirectionCovBegin_.push_back(photonFreeDirectionCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      photonFreeDirectionCov_.push_back(
          value(chi, "freePhotonDirectionCov" + std::to_string(row) + std::to_string(column)));
  photonFreeDirectionCovSize_.push_back(6);
  photonPromptDirectionX_.push_back(value(chi, "promptPhotonDirectionX"));
  photonPromptDirectionY_.push_back(value(chi, "promptPhotonDirectionY"));
  photonPromptDirectionZ_.push_back(value(chi, "promptPhotonDirectionZ"));
  photonPromptDirectionCovBegin_.push_back(photonPromptDirectionCov_.size());
  for (unsigned int row = 0; row != 3; ++row)
    for (unsigned int column = 0; column <= row; ++column)
      photonPromptDirectionCov_.push_back(
          value(chi, "promptPhotonDirectionCov" + std::to_string(row) + std::to_string(column)));
  photonPromptDirectionCovSize_.push_back(6);
  photonObservedElectronE_.push_back(value(chi, "observedElectronEnergy"));
  photonRawThreeBodyMass_.push_back(value(chi, "rawThreeBodyMass"));
  photonRawDeltaM_.push_back(value(chi, "rawPartialMassDifference"));
  photonVisibleMass_.push_back(value(chi, "geometryVisibleMass"));
  photonVisibleDeltaM_.push_back(value(chi, "geometryVisibleMassDifference"));
  photonPromptCompatibilityChi2_.push_back(value(chi, "promptCompatibilityChi2"));
  photonOriginDirectionDelta_.push_back(value(chi, "originDirectionDelta"));
  const float angleError = originAngleError(chi);
  photonOriginDirectionDeltaError_.push_back(angleError);
  photonOriginAngularSignificance_.push_back(
      std::isfinite(angleError) && angleError > 0.f ? value(chi, "originDirectionDelta") / angleError : kMissing);
  photonOriginDeltaChiC0E_.push_back(value(chi, "promptChiC0PhotonEnergy") - value(chi, "freeChiC0PhotonEnergy"));
  photonOriginDeltaChiC1E_.push_back(value(chi, "promptChiC1PhotonEnergy") - value(chi, "freeChiC1PhotonEnergy"));
  photonOriginDeltaChiC2E_.push_back(value(chi, "promptChiC2PhotonEnergy") - value(chi, "freeChiC2PhotonEnergy"));
  photonFreeStatus_.push_back(intValue(chi, "freeHypothesisStatus"));
  photonPromptStatus_.push_back(intValue(chi, "promptHypothesisStatus"));
  photonDirectionCrossCovStatus_.push_back(intValue(chi, "photonDirectionCrossCovarianceStatus"));
  photonFreeChiC0E_.push_back(value(chi, "freeChiC0PhotonEnergy"));
  photonFreeChiC1E_.push_back(value(chi, "freeChiC1PhotonEnergy"));
  photonFreeChiC2E_.push_back(value(chi, "freeChiC2PhotonEnergy"));
  photonFreeChiC0Err_.push_back(value(chi, "freeChiC0PhotonEnergyError"));
  photonFreeChiC1Err_.push_back(value(chi, "freeChiC1PhotonEnergyError"));
  photonFreeChiC2Err_.push_back(value(chi, "freeChiC2PhotonEnergyError"));
  photonPromptChiC0E_.push_back(value(chi, "promptChiC0PhotonEnergy"));
  photonPromptChiC1E_.push_back(value(chi, "promptChiC1PhotonEnergy"));
  photonPromptChiC2E_.push_back(value(chi, "promptChiC2PhotonEnergy"));
  photonPromptChiC0Err_.push_back(value(chi, "promptChiC0PhotonEnergyError"));
  photonPromptChiC1Err_.push_back(value(chi, "promptChiC1PhotonEnergyError"));
  photonPromptChiC2Err_.push_back(value(chi, "promptChiC2PhotonEnergyError"));
  photonFreeChiC0Status_.push_back(intValue(chi, "freeChiC0MassHypothesisStatus"));
  photonFreeChiC1Status_.push_back(intValue(chi, "freeChiC1MassHypothesisStatus"));
  photonFreeChiC2Status_.push_back(intValue(chi, "freeChiC2MassHypothesisStatus"));
  photonPromptChiC0Status_.push_back(intValue(chi, "promptChiC0MassHypothesisStatus"));
  photonPromptChiC1Status_.push_back(intValue(chi, "promptChiC1MassHypothesisStatus"));
  photonPromptChiC2Status_.push_back(intValue(chi, "promptChiC2MassHypothesisStatus"));
  appendElectronFromOneProng(*leg, photonIndex);
  photonElectronSize_.push_back(electronPx_.size() - electronBegin);
}

void ChiCCovarianceNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  clear();
  const auto& vertices = event.get(primaryVerticesToken_);
  const auto& onias = event.get(oniaStatesToken_);
  const auto& conventional = event.get(conventionalToken_);

  const auto oniaMuonsHandle = event.getHandle(oniaMuonsToken_);
  oniaMuons_ = oniaMuonsHandle.isValid() ? oniaMuonsHandle.product() : nullptr;

  const auto packedHandle = event.getHandle(packedCandidatesToken_);
  const auto lostHandle = event.getHandle(lostTracksToken_);
  if (packedHandle.isValid() && lostHandle.isValid()) {
    pvNTrack_ = chic::countAssociatedTracks(vertices, *packedHandle, *lostHandle, trackActivityConfig_);
  }

  const pat::CompositeCandidateCollection* oneProng = nullptr;
  if (doOneProng_) {
    const auto oneProngHandle = event.getHandle(oneProngToken_);
    if (oneProngHandle.isValid()) {
      oneProng = oneProngHandle.product();
    }
  }

  genParticles_ = nullptr;
  genIndex_.clear();
  if (isMC_) {
    const auto genParticlesHandle = event.getHandle(genParticlesToken_);
    if (genParticlesHandle.isValid()) {
      genParticles_ = genParticlesHandle.product();
      for (size_t index = 0; index < genParticles_->size(); ++index) {
        genIndex_[&genParticles_->at(index)] = static_cast<int>(index);
      }
    }
  }

  for (const auto& vertex : vertices) {
    pvX_.push_back(vertex.x());
    pvY_.push_back(vertex.y());
    pvZ_.push_back(vertex.z());
    pvCovBegin_.push_back(pvCov_.size());
    appendLowerTriangle(vertex.covariance(), 3, pvCov_);
    pvCovSize_.push_back(6);
  }

  for (const auto& onia : onias) {
    const unsigned int index = jpsiPx_.size(), muonBegin = muonPx_.size();
    jpsiPx_.push_back(onia.px());
    jpsiPy_.push_back(onia.py());
    jpsiPz_.push_back(onia.pz());
    jpsiE_.push_back(onia.energy());
    jpsiRawPx_.push_back(value(onia, "rawJpsiPx"));
    jpsiRawPy_.push_back(value(onia, "rawJpsiPy"));
    jpsiRawPz_.push_back(value(onia, "rawJpsiPz"));
    jpsiRawE_.push_back(value(onia, "rawJpsiE"));
    jpsiMass_.push_back(value(onia, "jpsiFitMass"));
    jpsiFitChi2_.push_back(value(onia, "jpsiFitChi2"));
    jpsiFitNdof_.push_back(value(onia, "jpsiFitNdof"));
    jpsiPromptChi2_.push_back(value(onia, "promptCompatibilityChi2"));
    jpsiFitStatus_.push_back(intValue(onia, "jpsiFitStatus"));
    jpsiPVIndex_.push_back(intValue(onia, "pvIndex"));
    jpsiPVCrossCovStatus_.push_back(intValue(onia, "pvJpsiCrossCovarianceStatus"));
    jpsiPromptStatus_.push_back(intValue(onia, "promptOriginStatus"));
    jpsiMuonBegin_.push_back(muonBegin);
    jpsiSVX_.push_back(value(onia, "jpsiVertexX"));
    jpsiSVY_.push_back(value(onia, "jpsiVertexY"));
    jpsiSVZ_.push_back(value(onia, "jpsiVertexZ"));
    jpsiPVX_.push_back(value(onia, "primaryVertexX"));
    jpsiPVY_.push_back(value(onia, "primaryVertexY"));
    jpsiPVZ_.push_back(value(onia, "primaryVertexZ"));
    jpsiPVDeltaX_.push_back(value(onia, "jpsiVertexX") - value(onia, "primaryVertexX"));
    jpsiPVDeltaY_.push_back(value(onia, "jpsiVertexY") - value(onia, "primaryVertexY"));
    jpsiPVDeltaZ_.push_back(value(onia, "jpsiVertexZ") - value(onia, "primaryVertexZ"));
    jpsiPVDeltaCovBegin_.push_back(jpsiPVDeltaCov_.size());
    for (unsigned int row = 0; row != 3; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        jpsiPVDeltaCov_.push_back(value(onia, "jpsiVertexCov" + std::to_string(row) + std::to_string(column)) +
                                  value(onia, "primaryVertexCov" + std::to_string(row) + std::to_string(column)));
    jpsiPVDeltaCovSize_.push_back(6);
    jpsiPromptX_.push_back(value(onia, "promptOriginX"));
    jpsiPromptY_.push_back(value(onia, "promptOriginY"));
    jpsiPromptZ_.push_back(value(onia, "promptOriginZ"));
    for (auto tuple : {std::make_pair(&jpsiSVCov_, "jpsiVertex"),
                       std::make_pair(&jpsiPVCov_, "primaryVertex"),
                       std::make_pair(&jpsiPromptCov_, "promptOrigin")}) {
      auto* output = tuple.first;
      std::vector<unsigned int>* begin =
          output == &jpsiSVCov_ ? &jpsiSVCovBegin_ : output == &jpsiPVCov_ ? &jpsiPVCovBegin_ : &jpsiPromptCovBegin_;
      std::vector<unsigned int>* size =
          output == &jpsiSVCov_ ? &jpsiSVCovSize_ : output == &jpsiPVCov_ ? &jpsiPVCovSize_ : &jpsiPromptCovSize_;
      begin->push_back(output->size());
      for (unsigned int row = 0; row != 3; ++row)
        for (unsigned int column = 0; column <= row; ++column)
          output->push_back(value(onia, std::string(tuple.second) + "Cov" + std::to_string(row) + std::to_string(column)));
      size->push_back(6);
    }
    jpsiStateCovBegin_.push_back(jpsiStateCov_.size());
    for (unsigned int row = 0; row != 7; ++row)
      for (unsigned int column = 0; column <= row; ++column)
        jpsiStateCov_.push_back(value(onia, "jpsiStateCov" + std::to_string(row) + std::to_string(column)));
    jpsiStateCovSize_.push_back(28);
    const auto* mu1 = dynamic_cast<const pat::Muon*>(onia.daughter("muon1"));
    const auto* mu2 = dynamic_cast<const pat::Muon*>(onia.daughter("muon2"));
    jpsiOniaIndex_.push_back(intValue(onia, "oniaStateIndex", static_cast<int>(index)));
    const pat::Muon* muPlus =
        (mu1 && mu1->charge() > 0) ? mu1 : ((mu2 && mu2->charge() > 0) ? mu2 : nullptr);
    const pat::Muon* muMinus =
        (mu1 && mu1->charge() < 0) ? mu1 : ((mu2 && mu2->charge() < 0) ? mu2 : nullptr);
    jpsiOniaMuPlusIndex_.push_back(muPlus ? oniaMuonIndex(*muPlus) : -1);
    jpsiOniaMuMinusIndex_.push_back(muMinus ? oniaMuonIndex(*muMinus) : -1);
    if (mu1)
      appendMuon(*mu1, index);
    if (mu2)
      appendMuon(*mu2, index);
    jpsiMuonSize_.push_back(muonPx_.size() - muonBegin);
  }

  for (const auto& candidate : conventional) {
    const size_t before = chiKind_.size();
    appendConventional(candidate);
    if (chiKind_.size() != before)
      appendTruth(candidate, 0);
  }
  if (oneProng) {
    for (const auto& candidate : *oneProng) {
      const size_t before = chiKind_.size();
      appendOneProng(candidate);
      if (chiKind_.size() != before)
        appendTruth(candidate, 1);
    }
  }
  for (const int vertexIndex : chiPVIndex_) {
    chiPVNTrack_.push_back(vertexIndex >= 0 && static_cast<size_t>(vertexIndex) < pvNTrack_.size()
                               ? static_cast<int>(pvNTrack_[vertexIndex])
                               : -1);
  }

  eventRun_.push_back(event.id().run());
  eventLumi_.push_back(event.id().luminosityBlock());
  eventNumber_.push_back(event.id().event());
  eventNPV_.push_back(pvX_.size());
  eventNJpsi_.push_back(jpsiPx_.size());
  eventNMuon_.push_back(muonPx_.size());
  eventNChi_.push_back(chiPx_.size());
  eventNPhoton_.push_back(photonPx_.size());
  eventNElectron_.push_back(electronPx_.size());
  tree_->Fill();
}

void ChiCCovarianceNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("oniaStates", edm::InputTag("chiCSharedOniaStates"));
  description.add<edm::InputTag>("oniaMuons", edm::InputTag("chiCSharedOniaSelectedMuons"));
  description.add<edm::InputTag>("conventionalCandidates", edm::InputTag("chiCConventionalCandidates"));
  description.add<edm::InputTag>("oneProngCandidates", edm::InputTag("chiCSharedOniaOneProngCandidates"));
  description.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  description.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  description.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  description.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
  description.add<double>("trackMinPt", 0.4);
  description.add<double>("trackMaxAbsEta", 2.4);
  description.add<double>("trackMaxRelPtError", 0.1);
  description.add<double>("trackMaxDzSignificance", 3.0);
  description.add<double>("trackMaxDxySignificance", 3.0);
  description.add<bool>("doOneProng", true);
  description.add<bool>("isMC", false);
  descriptions.addWithDefaultLabel(description);
}

DEFINE_FWK_MODULE(ChiCCovarianceNtuplizer);
