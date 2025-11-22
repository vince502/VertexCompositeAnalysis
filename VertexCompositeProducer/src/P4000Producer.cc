// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      P4000Producer
//
// Producer for P(4000) -> J/ψ(μ+μ-) + φ(K+K-)

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/P4000Producer.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <vector>

P4000Producer::P4000Producer(const edm::ParameterSet& cfg)
  : jpsiToken_(consumes<ResonanceCollection>(cfg.getParameter<edm::InputTag>("jpsiCollection"))),
    phiToken_(consumes<ResonanceCollection>(cfg.getParameter<edm::InputTag>("phiCollection"))),
    vertexToken_(cfg.exists("vertexRecoAlgorithm") ? consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertexRecoAlgorithm")) : edm::EDGetTokenT<reco::VertexCollection>()),
    beamSpotToken_(cfg.exists("beamSpot") ? consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot")) : edm::EDGetTokenT<reco::BeamSpot>()),
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    applyMassWindow_(cfg.getParameter<bool>("applyMassWindow")),
    requireUniqueTracks_(cfg.getParameter<bool>("requireUniqueTracks")),
    useVertexFitting_(cfg.getParameter<bool>("useVertexFitting"))
{
  if (cfg.exists("resonanceMassSigmas")) {
    resonanceMassSigmas_ = cfg.getParameter<std::vector<double> >("resonanceMassSigmas");
  } else {
    resonanceMassSigmas_ = {0.003096916, 0.001019445};  // J/ψ and φ mass uncertainties
  }

  const auto& statePsets = cfg.getParameter<std::vector<edm::ParameterSet> >("states");
  states_.reserve(statePsets.size());
  for (const auto& ps : statePsets) {
    P4000StateConfig state;
    state.name = ps.getParameter<std::string>("name");
    state.pdgId = ps.getParameter<int>("pdgId");
    state.mass = ps.getParameter<double>("mass");
    state.massWindow = ps.getParameter<double>("massWindow");
    states_.push_back(state);
    produces<P4000Collection>(state.name);
  }
}

P4000Producer::~P4000Producer() = default;

void P4000Producer::beginJob() {}

bool P4000Producer::shareTracks(const reco::VertexCompositeCandidate& jpsi,
                                 const reco::VertexCompositeCandidate& phi) const {
  std::set<reco::TrackRef> jpsiTracks;
  for (size_t i = 0; i < jpsi.numberOfDaughters(); ++i) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(jpsi.daughter(i));
    if (dau && dau->track().isNonnull()) {
      jpsiTracks.insert(dau->track());
    }
  }

  for (size_t i = 0; i < phi.numberOfDaughters(); ++i) {
    const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(phi.daughter(i));
    if (dau && dau->track().isNonnull()) {
      if (jpsiTracks.find(dau->track()) != jpsiTracks.end()) {
        return true;
      }
    }
  }
  return false;
}

void P4000Producer::produce(edm::Event& event, const edm::EventSetup& iSetup) {
  edm::Handle<ResonanceCollection> jpsiHandle;
  edm::Handle<ResonanceCollection> phiHandle;
  event.getByToken(jpsiToken_, jpsiHandle);
  event.getByToken(phiToken_, phiHandle);

  if (!jpsiHandle.isValid() || !phiHandle.isValid()) return;
  if (jpsiHandle->empty() || phiHandle->empty()) return;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::BeamSpot> beamSpot;
  const MagneticField* magField = nullptr;

  if (useVertexFitting_) {
    if (!vertexToken_.isUninitialized()) {
      event.getByToken(vertexToken_, vertices);
    }
    if (!beamSpotToken_.isUninitialized()) {
      event.getByToken(beamSpotToken_, beamSpot);
    }
    magField = &iSetup.getData(bFieldToken_);
    if (!vertices.isValid() || !beamSpot.isValid()) {
      throw cms::Exception("InvalidInput") << "Vertex or BeamSpot not available but useVertexFitting is enabled";
    }
  }

  const reco::Vertex* primaryVertex = nullptr;
  if (vertices.isValid() && !vertices->empty()) {
    primaryVertex = &vertices->front();
  }

  for (const auto& state : states_) {
    auto output = std::make_unique<P4000Collection>();

    for (const auto& jpsi : *jpsiHandle) {
      for (const auto& phi : *phiHandle) {
        if (requireUniqueTracks_ && shareTracks(jpsi, phi)) continue;

        double mass = (jpsi.p4() + phi.p4()).mass();
        if (applyMassWindow_) {
          if (mass < state.mass - state.massWindow || mass > state.mass + state.massWindow) continue;
        }

        pat::CompositeCandidate p4000;
        
        // Create pat::CompositeCandidate copies from reco::VertexCompositeCandidate
        pat::CompositeCandidate jpsiPat;
        jpsiPat.setP4(jpsi.p4());
        jpsiPat.setVertex(jpsi.vertex());
        jpsiPat.setPdgId(jpsi.pdgId());
        for (size_t i = 0; i < jpsi.numberOfDaughters(); ++i) {
          jpsiPat.addDaughter(*jpsi.daughter(i));
        }
        
        pat::CompositeCandidate phiPat;
        phiPat.setP4(phi.p4());
        phiPat.setVertex(phi.vertex());
        phiPat.setPdgId(phi.pdgId());
        for (size_t i = 0; i < phi.numberOfDaughters(); ++i) {
          phiPat.addDaughter(*phi.daughter(i));
        }
        
        p4000.addDaughter(jpsiPat, "Jpsi");
        p4000.addDaughter(phiPat, "Phi");

        AddFourMomenta addP4;
        addP4.set(p4000);

        if (useVertexFitting_ && magField && primaryVertex) {
          // Kinematic vertex fit
          std::vector<reco::TransientTrack> transTracks;
          
          // Get tracks from J/ψ daughters
          for (size_t i = 0; i < jpsi.numberOfDaughters(); ++i) {
            const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(jpsi.daughter(i));
            if (dau && dau->track().isNonnull()) {
              transTracks.push_back(reco::TransientTrack(*dau->track(), magField));
            }
          }
          
          // Get tracks from φ daughters
          for (size_t i = 0; i < phi.numberOfDaughters(); ++i) {
            const auto* dau = dynamic_cast<const reco::RecoChargedCandidate*>(phi.daughter(i));
            if (dau && dau->track().isNonnull()) {
              transTracks.push_back(reco::TransientTrack(*dau->track(), magField));
            }
          }

          if (transTracks.size() == 4) {
            KinematicParticleFactoryFromTransientTrack pFactory;
            std::vector<RefCountedKinematicParticle> particles;

            // J/ψ daughters (muons)
            float muonMass = 0.1056583745f;
            float muonMassSigma = 0.0001056583745f;
            float chi = 0.0f;
            float ndf = 0.0f;
            particles.push_back(pFactory.particle(transTracks[0], muonMass, chi, ndf, muonMassSigma));
            particles.push_back(pFactory.particle(transTracks[1], muonMass, chi, ndf, muonMassSigma));
            
            // φ daughters (kaons)
            float kaonMass = 0.493677f;
            float kaonMassSigma = 0.000493677f;
            particles.push_back(pFactory.particle(transTracks[2], kaonMass, chi, ndf, kaonMassSigma));
            particles.push_back(pFactory.particle(transTracks[3], kaonMass, chi, ndf, kaonMassSigma));

            // Fit J/ψ
            KinematicParticleVertexFitter jpsiFitter;
            RefCountedKinematicTree jpsiTree = jpsiFitter.fit(std::vector<RefCountedKinematicParticle>(particles.begin(), particles.begin() + 2));
            
            // Fit φ
            RefCountedKinematicTree phiTree = jpsiFitter.fit(std::vector<RefCountedKinematicParticle>(particles.begin() + 2, particles.end()));

            if (jpsiTree->isValid() && phiTree->isValid()) {
              jpsiTree->movePointerToTheTop();
              phiTree->movePointerToTheTop();
              
              RefCountedKinematicParticle jpsiParticle = jpsiTree->currentParticle();
              RefCountedKinematicParticle phiParticle = phiTree->currentParticle();

              // Fit P(4000)
              std::vector<RefCountedKinematicParticle> p4000Particles;
              p4000Particles.push_back(jpsiParticle);
              p4000Particles.push_back(phiParticle);

              RefCountedKinematicTree p4000Tree = jpsiFitter.fit(p4000Particles);

              if (p4000Tree->isValid()) {
                p4000Tree->movePointerToTheTop();
                RefCountedKinematicParticle p4000Particle = p4000Tree->currentParticle();
                
                GlobalVector p4000Mom = p4000Particle->currentState().globalMomentum();
                double p4000E = p4000Particle->currentState().kinematicParameters().energy();
                p4000.setP4(reco::Particle::LorentzVector(p4000Mom.x(), p4000Mom.y(), p4000Mom.z(), p4000E));

                RefCountedKinematicVertex p4000Vtx = p4000Tree->currentDecayVertex();
                if (p4000Vtx->vertexIsValid()) {
                  GlobalPoint vtxPos = p4000Vtx->position();
                  p4000.setVertex(reco::Candidate::Point(vtxPos.x(), vtxPos.y(), vtxPos.z()));
                  
                  double vtxChi2 = p4000Vtx->chiSquared();
                  double vtxNdof = p4000Vtx->degreesOfFreedom();
                  p4000.addUserFloat("VtxChi2", vtxChi2);
                  p4000.addUserFloat("VtxNdof", vtxNdof);
                  p4000.addUserFloat("VtxProb", TMath::Prob(vtxChi2, vtxNdof));
                }
              }
            }
          }
        } else {
          // Simple combination without vertex fit
          reco::Candidate::Point jpsiVtx = jpsi.vertex();
          reco::Candidate::Point phiVtx = phi.vertex();
          reco::Candidate::Point vtxPos(
            (jpsiVtx.x() + phiVtx.x()) * 0.5,
            (jpsiVtx.y() + phiVtx.y()) * 0.5,
            (jpsiVtx.z() + phiVtx.z()) * 0.5
          );
          p4000.setVertex(vtxPos);
        }

        p4000.setPdgId(state.pdgId);
        p4000.addUserFloat("jpsiMass", jpsi.mass());
        p4000.addUserFloat("phiMass", phi.mass());
        p4000.addUserFloat("jpsiPt", jpsi.pt());
        p4000.addUserFloat("phiPt", phi.pt());

        output->push_back(p4000);
      }
    }

    event.put(std::move(output), state.name);
  }
}

void P4000Producer::endJob() {}

#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(P4000Producer);
