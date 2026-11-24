/**
   \file
   Declaration of PoniaOniaPhotonProducer
   \author
   Alberto Sanchez-Hernandez (September 2014)
   Unified for VertexCompositeAnalysis
*/

#ifndef VertexCompositeAnalysis_VertexCompositeProducer_PoniaOniaPhotonProducer_h
#define VertexCompositeAnalysis_VertexCompositeProducer_PoniaOniaPhotonProducer_h

#include <cmath>
#include <limits>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

/**
   Create a Chi(b,c) candidate by matching dimuon and conversion
 */
class PoniaOniaPhotonProducer : public edm::stream::EDProducer<> {
public:
  explicit PoniaOniaPhotonProducer(const edm::ParameterSet& ps);
  ~PoniaOniaPhotonProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event& event, const edm::EventSetup& esetup) override;
  void endStream() override;

  pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate& dimuon,
                                           const pat::CompositeCandidate& photon,
                                           const reco::Candidate::Point* oniaVertexPos) const;

  float Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point& p) const;
  bool cutDeltaMass(const pat::CompositeCandidate& chiCand, const pat::CompositeCandidate& dimuonCand) const;
  bool cutdz(float dz) const { return std::isfinite(dz) && dz < dzMax_; }

  edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> photon_Label_;

  bool pi0OnlineSwitch_;
  std::vector<double> deltaMass_;
  double dzMax_;
  bool triggerMatch_;
  double minPhotonPt_;
  bool applyPairCuts_;
  bool requireOniaVertex_;

  int candidates_;
  int delta_mass_fail_;
  int dz_cut_fail_;
  int pizero_fail_;
  int missing_vertex_fail_;
};

#endif
