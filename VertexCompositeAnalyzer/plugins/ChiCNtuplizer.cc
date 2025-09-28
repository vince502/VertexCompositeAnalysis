#include "VertexCompositeAnalysis/VertexCompositeAnalyzer/plugins/ChiCNtuplizer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"

ChiCNtuplizer::ChiCNtuplizer(const edm::ParameterSet& cfg)
    : tree_(nullptr) {
  usesResource("TFileService");

  treeName_ = cfg.getUntrackedParameter<std::string>("treeName", "ChiCNtuple");
  storeDaughterInfo_ = cfg.getUntrackedParameter<bool>("storeDaughterInfo", true);

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

  if (storeDaughterInfo_) {
    tree_->Branch("cand_dauStart", &cand_dauStart_);
    tree_->Branch("cand_dauCount", &cand_dauCount_);
    tree_->Branch("dau_pt", &dau_pt_);
    tree_->Branch("dau_eta", &dau_eta_);
    tree_->Branch("dau_phi", &dau_phi_);
    tree_->Branch("dau_mass", &dau_mass_);
    tree_->Branch("dau_charge", &dau_charge_);
    tree_->Branch("dau_pdgId", &dau_pdgId_);
  }
}

void ChiCNtuplizer::endJob() {}

void ChiCNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  resetEventContent();

  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  for (const auto& src : sources_) {
    edm::Handle<pat::CompositeCandidateCollection> handle;
    event.getByToken(src.token, handle);
    if (!handle.isValid()) {
      continue;
    }
    fillCandidates(src, *handle);
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

  if (storeDaughterInfo_) {
    cand_dauStart_.clear();
    cand_dauCount_.clear();
    dau_pt_.clear();
    dau_eta_.clear();
    dau_phi_.clear();
    dau_mass_.clear();
    dau_charge_.clear();
    dau_pdgId_.clear();
  }
}

void ChiCNtuplizer::fillCandidates(const SourceConfig& src,
                                   const pat::CompositeCandidateCollection& candidates) {
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

    if (storeDaughterInfo_) {
      cand_dauStart_.push_back(static_cast<unsigned int>(dau_pt_.size()));
      cand_dauCount_.push_back(nDau);

      for (unsigned int i = 0; i < nDau; ++i) {
        const auto* dau = cand.daughter(i);
        if (!dau) {
          dau_pt_.push_back(-1.f);
          dau_eta_.push_back(-99.f);
          dau_phi_.push_back(-99.f);
          dau_mass_.push_back(-1.f);
          dau_charge_.push_back(0);
          dau_pdgId_.push_back(0);
          continue;
        }

        dau_pt_.push_back(static_cast<float>(dau->pt()));
        dau_eta_.push_back(static_cast<float>(dau->eta()));
        dau_phi_.push_back(static_cast<float>(dau->phi()));
        dau_mass_.push_back(static_cast<float>(dau->mass()));
        dau_charge_.push_back(dau->charge());
        dau_pdgId_.push_back(dau->pdgId());
      }
    }
  }
}

DEFINE_FWK_MODULE(ChiCNtuplizer);
