// -*- C++ -*-
#ifndef VertexCompositeAnalysis__ChiCFlatNtuplizer_h
#define VertexCompositeAnalysis__ChiCFlatNtuplizer_h

#include <array>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class TTree;

class ChiCFlatNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCFlatNtuplizer(const edm::ParameterSet&);
  ~ChiCFlatNtuplizer() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

private:
  struct SourceConfig {
    edm::InputTag tag;
    edm::EDGetTokenT<pat::CompositeCandidateCollection> token;
    std::string name;
    int pdgId;
    unsigned int index;
  };

  using VertexCollection = reco::VertexCollection;

  void resetBranches();
  void fillCandidate(const SourceConfig&, 
                     const pat::CompositeCandidate&, 
                     const reco::Vertex* primaryVertex,
                     const edm::Handle<edm::View<pat::IsolatedTrack>>& isoTracksHandle);

  edm::EDGetTokenT<VertexCollection> pvToken_;
  edm::EDGetTokenT<edm::View<pat::IsolatedTrack>> isoTracksToken_;
  std::vector<SourceConfig> sources_;
  std::string treeName_;
  bool useDeDx_;

  edm::Service<TFileService> fileService_;
  TTree* tree_;

  // Per-entry content (one entry per candidate)
  unsigned int run_{};
  unsigned int lumi_{};
  unsigned long long event_{};
  int sourceIndex_{};
  int sourcePdgId_{};
  std::string sourceLabel_{};

  int candCharge_{};
  float candMass_{};
  float candPt_{};
  float candEta_{};
  float candPhi_{};
  float candRapidity_{};
  float candVx_{};
  float candVy_{};
  float candVz_{};
  float candD3D_{};

  float candAcoplanarity_{};
  float candSphericity_{};
  float candLambda1_{};
  float candLambda2_{};
  float candLambda3_{};

  std::array<float, 4> dauPt_{};
  std::array<float, 4> dauEta_{};
  std::array<float, 4> dauPhi_{};
  std::array<int, 4> dauCharge_{};
  std::array<float, 4> dauDxy_{};
  std::array<float, 4> dauDz_{};
  std::array<float, 4> dauD3d_{};
  std::array<float, 4> dauMass_{};
  std::array<float, 4> dauDeDx_{};

  std::array<float, 6> pairDca_{};
};

#endif
