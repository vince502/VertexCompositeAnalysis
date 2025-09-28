#ifndef VertexCompositeAnalysis__ChiCGenNtuplizer_h
#define VertexCompositeAnalysis__ChiCGenNtuplizer_h

#include <array>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class ChiCGenNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCGenNtuplizer(const edm::ParameterSet&);
  ~ChiCGenNtuplizer() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

private:
  using GenParticleCollection = reco::GenParticleCollection;
  using VertexCollection = reco::VertexCollection;
  using XYZPoint = math::XYZPoint;
  using XYZVector = math::XYZVector;

  void resetBranches();
  void collectStablePions(const reco::GenParticle& particle,
                          std::vector<const reco::GenParticle*>& pions) const;
  static float distance3D(const XYZPoint& from, const XYZPoint& to);
  static float pairDca(const reco::GenParticle& first, const reco::GenParticle& second);

  edm::EDGetTokenT<GenParticleCollection> genToken_;
  edm::EDGetTokenT<VertexCollection> pvToken_;

  std::string treeName_;

  edm::Service<TFileService> fileService_;
  TTree* tree_;

  unsigned int run_;
  unsigned int lumi_;
  unsigned long long event_;
  int chic_found_;
  float chic_pt_;
  float chic_eta_;
  float chic_phi_;
  float chic_y_;
  float chic_mass_;
  float chic_d3d_;

  std::array<float, 4> pion_pt_;
  std::array<float, 4> pion_eta_;
  std::array<float, 4> pion_phi_;
  std::array<float, 4> pion_d3d_;

  std::array<float, 6> dca_pairs_;
};

#endif
