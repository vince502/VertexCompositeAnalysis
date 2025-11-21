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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
  using GenParticleCollection = reco::GenParticleCollection;

  void resetBranches();
  void fillCandidate(const SourceConfig&, 
                     const pat::CompositeCandidate&, 
                     const reco::Vertex* primaryVertex,
                     const edm::Handle<GenParticleCollection>* genParticles = nullptr);
  
  const reco::GenParticle* findGenMatch(const pat::CompositeCandidate& cand,
                                        const edm::Handle<GenParticleCollection>& genParticles) const;
  void collectStablePions(const reco::GenParticle& particle,
                          std::vector<const reco::GenParticle*>& pions) const;

  edm::EDGetTokenT<VertexCollection> pvToken_;
  edm::EDGetTokenT<GenParticleCollection> genToken_;
  bool useGenMatching_;
  std::vector<SourceConfig> sources_;
  std::string treeName_;

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
  float candDecayLength3D_{};
  float candDecayLength2D_{};
  float candPointingAngle3D_{};
  float candPointingAngle2D_{};
  float candCosPointingAngle3D_{};
  float candCosPointingAngle2D_{};
  float candDxy_{};
  float candDz_{};
  float candDxySig_{};
  float candDzSig_{};

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
  std::array<float, 4> dauDxySig_{};
  std::array<float, 4> dauDzSig_{};
  std::array<float, 4> dauChi2_{};
  std::array<int, 4> dauNhits_{};
  std::array<int, 4> dauNpixHits_{};
  std::array<float, 4> dauPtErr_{};
  std::array<float, 4> dauEtaErr_{};
  std::array<float, 4> dauPhiErr_{};

  std::array<float, 6> pairDca_{};
  std::array<float, 6> pairMass_{};
  std::array<float, 6> pairPt_{};
  std::array<float, 6> pairEta_{};
  std::array<float, 6> pairPhi_{};

  // Event-level information
  int nPV_{};
  float pvX_{};
  float pvY_{};
  float pvZ_{};
  float pvNdof_{};
  float pvChi2_{};

  // Generator-level and q-vector variables
  int genMatch_{};
  float genMass_{};
  float genPt_{};
  float genEta_{};
  float genPhi_{};
  float genY_{};
  float q2Magnitude_{};
  float q2Phase_{};
};

#endif
