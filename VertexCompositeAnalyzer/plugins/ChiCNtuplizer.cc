#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCNtuplizer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <cmath>

ChiCNtuplizer::ChiCNtuplizer(const edm::ParameterSet& cfg)
    : tree_(nullptr) {
  usesResource("TFileService");

  treeName_ = cfg.getUntrackedParameter<std::string>("treeName", "ChiCNtuple");
  storeDaughterInfo_ = cfg.getUntrackedParameter<bool>("storeDaughterInfo", true);

  // Primary vertex collection (optional)
  if (cfg.existsAs<edm::InputTag>("primaryVertices")) {
    const auto pvTag = cfg.getParameter<edm::InputTag>("primaryVertices");
    pvToken_ = consumes<reco::VertexCollection>(pvTag);
  }

  const auto& sourcePsets = cfg.getParameter<std::vector<edm::ParameterSet>>("sources");
  if (sourcePsets.empty()) {
    throw cms::Exception("Configuration")
        << "ChiCNtuplizer requires at least one entry in 'sources'.";
  }

  sources_.reserve(sourcePsets.size());
  for (size_t idx = 0; idx < sourcePsets.size(); ++idx) {
    const auto& ps = sourcePsets[idx];
    SourceConfig sc;
    sc.tag = ps.getParameter<edm::InputTag>("collection");
    sc.name = ps.getParameter<std::string>("name");
    sc.pdgId = ps.getParameter<int>("pdgId");
    sc.index = static_cast<unsigned int>(idx);
    sc.token = consumes<pat::CompositeCandidateCollection>(sc.tag);
    sources_.push_back(sc);
  }
}

ChiCNtuplizer::~ChiCNtuplizer() = default;

void ChiCNtuplizer::beginJob() {
  if (!fileService_.isAvailable()) {
    throw cms::Exception("Configuration")
        << "TFileService is not available. Please add it to the process.";
  }

  tree_ = fileService_->make<TTree>(treeName_.c_str(), treeName_.c_str());

  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");

  tree_->Branch("cand_type", &cand_type_);
  tree_->Branch("cand_label", &cand_label_);
  tree_->Branch("cand_pdgId", &cand_pdgId_);
  tree_->Branch("cand_mass", &cand_mass_);
  tree_->Branch("cand_pt", &cand_pt_);
  tree_->Branch("cand_eta", &cand_eta_);
  tree_->Branch("cand_phi", &cand_phi_);
  tree_->Branch("cand_y", &cand_y_);
  tree_->Branch("cand_vx", &cand_vx_);
  tree_->Branch("cand_vy", &cand_vy_);
  tree_->Branch("cand_vz", &cand_vz_);
  tree_->Branch("cand_charge", &cand_charge_);
  tree_->Branch("cand_nDau", &cand_nDau_);
  
  // Geometry branches
  tree_->Branch("cand_d3d", &cand_d3d_);
  tree_->Branch("cand_decayLength3D", &cand_decayLength3D_);
  tree_->Branch("cand_decayLength2D", &cand_decayLength2D_);
  tree_->Branch("cand_pointingAngle3D", &cand_pointingAngle3D_);
  tree_->Branch("cand_pointingAngle2D", &cand_pointingAngle2D_);
  tree_->Branch("cand_cosPointingAngle3D", &cand_cosPointingAngle3D_);
  tree_->Branch("cand_cosPointingAngle2D", &cand_cosPointingAngle2D_);

  if (storeDaughterInfo_) {
    tree_->Branch("cand_dauStart", &cand_dauStart_);
    tree_->Branch("cand_dauCount", &cand_dauCount_);
    tree_->Branch("dau_pt", &dau_pt_);
    tree_->Branch("dau_eta", &dau_eta_);
    tree_->Branch("dau_phi", &dau_phi_);
    tree_->Branch("dau_mass", &dau_mass_);
    tree_->Branch("dau_charge", &dau_charge_);
    tree_->Branch("dau_pdgId", &dau_pdgId_);
    tree_->Branch("dau_dxy", &dau_dxy_);
    tree_->Branch("dau_dz", &dau_dz_);
    tree_->Branch("dau_d3d", &dau_d3d_);
    tree_->Branch("pair_dca", &pair_dca_);
  }
}

void ChiCNtuplizer::endJob() {}

void ChiCNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  resetEventContent();

  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  // Get primary vertex
  const reco::Vertex* primaryVertex = nullptr;
  if (!pvToken_.isUninitialized()) {
    edm::Handle<reco::VertexCollection> pvHandle;
    if (event.getByToken(pvToken_, pvHandle) && pvHandle.isValid() && !pvHandle->empty()) {
      primaryVertex = &pvHandle->front();
    }
  }

  for (const auto& src : sources_) {
    edm::Handle<pat::CompositeCandidateCollection> handle;
    event.getByToken(src.token, handle);
    if (!handle.isValid()) {
      continue;
    }
    fillCandidates(src, *handle, primaryVertex);
  }

  tree_->Fill();
}

void ChiCNtuplizer::resetEventContent() {
  cand_type_.clear();
  cand_label_.clear();
  cand_pdgId_.clear();
  cand_mass_.clear();
  cand_pt_.clear();
  cand_eta_.clear();
  cand_phi_.clear();
  cand_y_.clear();
  cand_vx_.clear();
  cand_vy_.clear();
  cand_vz_.clear();
  cand_charge_.clear();
  cand_nDau_.clear();
  
  cand_d3d_.clear();
  cand_decayLength3D_.clear();
  cand_decayLength2D_.clear();
  cand_pointingAngle3D_.clear();
  cand_pointingAngle2D_.clear();
  cand_cosPointingAngle3D_.clear();
  cand_cosPointingAngle2D_.clear();

  if (storeDaughterInfo_) {
    cand_dauStart_.clear();
    cand_dauCount_.clear();
    dau_pt_.clear();
    dau_eta_.clear();
    dau_phi_.clear();
    dau_mass_.clear();
    dau_charge_.clear();
    dau_pdgId_.clear();
    dau_dxy_.clear();
    dau_dz_.clear();
    dau_d3d_.clear();
    pair_dca_.clear();
  }
}

void ChiCNtuplizer::fillCandidates(const SourceConfig& src,
                                   const pat::CompositeCandidateCollection& candidates,
                                   const reco::Vertex* primaryVertex) {
  for (const auto& cand : candidates) {
    cand_type_.push_back(static_cast<int>(src.index));
    cand_label_.push_back(src.name);

    const int pdgId = cand.pdgId() != 0 ? cand.pdgId() : src.pdgId;
    cand_pdgId_.push_back(pdgId);

    cand_mass_.push_back(static_cast<float>(cand.mass()));
    cand_pt_.push_back(static_cast<float>(cand.pt()));
    cand_eta_.push_back(static_cast<float>(cand.eta()));
    cand_phi_.push_back(static_cast<float>(cand.phi()));
    cand_y_.push_back(static_cast<float>(cand.rapidity()));

    const auto& vtx = cand.vertex();
    cand_vx_.push_back(static_cast<float>(vtx.x()));
    cand_vy_.push_back(static_cast<float>(vtx.y()));
    cand_vz_.push_back(static_cast<float>(vtx.z()));

    cand_charge_.push_back(cand.charge());

    const auto nDau = static_cast<unsigned int>(cand.numberOfDaughters());
    cand_nDau_.push_back(nDau);

    // Compute geometry variables
    float d3d = -1.f;
    float decayLength3D = -1.f;
    float decayLength2D = -1.f;
    float pointingAngle3D = -99.f;
    float pointingAngle2D = -99.f;
    float cosPointingAngle3D = -99.f;
    float cosPointingAngle2D = -99.f;

    if (primaryVertex) {
      const auto& pvPos = primaryVertex->position();
      const double dx = vtx.x() - pvPos.x();
      const double dy = vtx.y() - pvPos.y();
      const double dz = vtx.z() - pvPos.z();
      
      d3d = static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
      decayLength3D = d3d;
      decayLength2D = static_cast<float>(std::sqrt(dx * dx + dy * dy));
      
      // Pointing angle: angle between candidate momentum and vector from PV to decay vertex
      const math::XYZVector flightDir(dx, dy, dz);
      const math::XYZVector candMom(cand.px(), cand.py(), cand.pz());
      
      const double flightDirMag = std::sqrt(flightDir.mag2());
      const double candMomMag = std::sqrt(candMom.mag2());
      
      if (flightDirMag > 0 && candMomMag > 0) {
        const double cosAngle3D = flightDir.Dot(candMom) / (flightDirMag * candMomMag);
        cosPointingAngle3D = static_cast<float>(cosAngle3D);
        pointingAngle3D = static_cast<float>(std::acos(std::max(-1.0, std::min(1.0, cosAngle3D))));
        
        const math::XYZVector flightDir2D(dx, dy, 0);
        const math::XYZVector candMom2D(cand.px(), cand.py(), 0);
        const double flightDir2DMag = std::sqrt(flightDir2D.mag2());
        const double candMom2DMag = std::sqrt(candMom2D.mag2());
        if (flightDir2DMag > 0 && candMom2DMag > 0) {
          const double cosAngle2D = flightDir2D.Dot(candMom2D) / (flightDir2DMag * candMom2DMag);
          cosPointingAngle2D = static_cast<float>(cosAngle2D);
          pointingAngle2D = static_cast<float>(std::acos(std::max(-1.0, std::min(1.0, cosAngle2D))));
        }
      }
    }
    
    cand_d3d_.push_back(d3d);
    cand_decayLength3D_.push_back(decayLength3D);
    cand_decayLength2D_.push_back(decayLength2D);
    cand_pointingAngle3D_.push_back(pointingAngle3D);
    cand_pointingAngle2D_.push_back(pointingAngle2D);
    cand_cosPointingAngle3D_.push_back(cosPointingAngle3D);
    cand_cosPointingAngle2D_.push_back(cosPointingAngle2D);

    if (storeDaughterInfo_) {
      cand_dauStart_.push_back(static_cast<unsigned int>(dau_pt_.size()));
      cand_dauCount_.push_back(nDau);
      
      // Store daughter indices for DCA calculation
      std::vector<unsigned int> dauIndices;
      dauIndices.reserve(nDau);

      for (unsigned int i = 0; i < nDau; ++i) {
        const auto* dau = cand.daughter(i);
        if (!dau) {
          dau_pt_.push_back(-1.f);
          dau_eta_.push_back(-99.f);
          dau_phi_.push_back(-99.f);
          dau_mass_.push_back(-1.f);
          dau_charge_.push_back(0);
          dau_pdgId_.push_back(0);
          dau_dxy_.push_back(0.f);
          dau_dz_.push_back(0.f);
          dau_d3d_.push_back(0.f);
          continue;
        }

        dau_pt_.push_back(static_cast<float>(dau->pt()));
        dau_eta_.push_back(static_cast<float>(dau->eta()));
        dau_phi_.push_back(static_cast<float>(dau->phi()));
        dau_mass_.push_back(static_cast<float>(dau->mass()));
        dau_charge_.push_back(dau->charge());
        dau_pdgId_.push_back(dau->pdgId());
        
        // Get track for impact parameters
        const reco::Track* trackPtr = nullptr;
        if (const auto* recoDau = dynamic_cast<const reco::RecoChargedCandidate*>(dau)) {
          const auto trackRef = recoDau->track();
          if (!trackRef.isNull())
            trackPtr = trackRef.get();
        }
        if (!trackPtr)
          trackPtr = dau->bestTrack();
        
        if (trackPtr && primaryVertex) {
          const auto& pvPos = primaryVertex->position();
          const double dxy = trackPtr->dxy(pvPos);
          const double dz = trackPtr->dz(pvPos);
          dau_dxy_.push_back(static_cast<float>(dxy));
          dau_dz_.push_back(static_cast<float>(dz));
          dau_d3d_.push_back(static_cast<float>(std::sqrt(dxy * dxy + dz * dz)));
        } else {
          dau_dxy_.push_back(0.f);
          dau_dz_.push_back(0.f);
          dau_d3d_.push_back(0.f);
        }
        
        dauIndices.push_back(static_cast<unsigned int>(dau_pt_.size() - 1));
      }
      
      // Compute pairwise DCA for daughters (for 4-pion: 6 pairs)
      // Only compute if we have exactly 4 daughters
      if (nDau == 4 && dauIndices.size() == 4) {
        const auto computePairDCA = [&](unsigned int i, unsigned int j) -> float {
          if (i >= dauIndices.size() || j >= dauIndices.size())
            return 0.f;
          
          const auto* dau1 = cand.daughter(i);
          const auto* dau2 = cand.daughter(j);
          if (!dau1 || !dau2)
            return 0.f;
          
          const auto& v1 = dau1->vertex();
          const auto& v2 = dau2->vertex();
          const double dx = v1.x() - v2.x();
          const double dy = v1.y() - v2.y();
          const double dz = v1.z() - v2.z();
          return static_cast<float>(std::sqrt(dx * dx + dy * dy + dz * dz));
        };
        
        // Pairs: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
        pair_dca_.push_back(computePairDCA(0, 1));
        pair_dca_.push_back(computePairDCA(0, 2));
        pair_dca_.push_back(computePairDCA(0, 3));
        pair_dca_.push_back(computePairDCA(1, 2));
        pair_dca_.push_back(computePairDCA(1, 3));
        pair_dca_.push_back(computePairDCA(2, 3));
      } else {
        // Not 4 daughters, fill with zeros
        for (int i = 0; i < 6; ++i) {
          pair_dca_.push_back(0.f);
        }
      }
    }
  }
}

DEFINE_FWK_MODULE(ChiCNtuplizer);
