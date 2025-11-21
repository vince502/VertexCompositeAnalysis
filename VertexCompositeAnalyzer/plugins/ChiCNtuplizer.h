// -*- C++ -*-
#ifndef VertexCompositeAnalysis__ChiCNtuplizer_h
#define VertexCompositeAnalysis__ChiCNtuplizer_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class ChiCNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ChiCNtuplizer(const edm::ParameterSet&);
  ~ChiCNtuplizer() override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

private:
  struct SourceConfig {
    edm::InputTag tag;
    edm::EDGetTokenT<pat::CompositeCandidateCollection> token;
    std::string name;
    int pdgId;
    unsigned int index;
  };

  void resetEventContent();
  void fillCandidates(const SourceConfig&, const pat::CompositeCandidateCollection&, const reco::Vertex* primaryVertex);

  std::vector<SourceConfig> sources_;
  std::string treeName_;
  bool storeDaughterInfo_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;

  edm::Service<TFileService> fileService_;
  TTree* tree_;

  // Event-level content
  unsigned int run_;
  unsigned int lumi_;
  unsigned long long event_;

  // Candidate-level content
  std::vector<int> cand_type_;
  std::vector<std::string> cand_label_;
  std::vector<int> cand_pdgId_;
  std::vector<float> cand_mass_;
  std::vector<float> cand_pt_;
  std::vector<float> cand_eta_;
  std::vector<float> cand_phi_;
  std::vector<float> cand_y_;
  std::vector<float> cand_vx_;
  std::vector<float> cand_vy_;
  std::vector<float> cand_vz_;
  std::vector<int> cand_charge_;
  std::vector<unsigned int> cand_nDau_;
  
  // Geometry variables
  std::vector<float> cand_d3d_;           // 3D distance from PV
  std::vector<float> cand_decayLength3D_;  // 3D decay length
  std::vector<float> cand_decayLength2D_; // 2D decay length
  std::vector<float> cand_pointingAngle3D_; // 3D pointing angle
  std::vector<float> cand_pointingAngle2D_; // 2D pointing angle
  std::vector<float> cand_cosPointingAngle3D_; // 3D cos pointing angle
  std::vector<float> cand_cosPointingAngle2D_; // 2D cos pointing angle

  // Daughter-level content (optional)
  std::vector<unsigned int> cand_dauStart_;
  std::vector<unsigned int> cand_dauCount_;
  std::vector<float> dau_pt_;
  std::vector<float> dau_eta_;
  std::vector<float> dau_phi_;
  std::vector<float> dau_mass_;
  std::vector<int> dau_charge_;
  std::vector<int> dau_pdgId_;
  
  // Daughter geometry variables
  std::vector<float> dau_dxy_;    // dxy impact parameter
  std::vector<float> dau_dz_;      // dz impact parameter
  std::vector<float> dau_d3d_;     // 3D impact parameter
  
  // Pair DCA (for 4-pion: 6 pairs)
  std::vector<float> pair_dca_;
};

#endif
