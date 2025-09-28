#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <iostream>
#include <string>
#include <vector>

// Usage within ROOT (6.32 or newer):
// root -l -q 'flattenChiCNtuple.C("input.root", "ChiCNtuple", "ChiCNtuple_flat.root")'
// The macro writes two trees:
//   - ChiC:     one entry per reconstructed ChiC candidate
//   - ChiC_dau: optional, one entry per daughter (only if daughter info exists)

void flattenChiCNtuple(const char* inputFile,
                       const char* treeName="ChiCNtuple",
                       const char* outputFile="ChiCNtuple_flat.root") {
  TFile inFile(inputFile, "READ");
  if (inFile.IsZombie()) {
    std::cerr << "[flattenChiCNtuple] Failed to open input file " << inputFile << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(inFile.Get(treeName));
  if (!tree) {
    std::cerr << "[flattenChiCNtuple] Tree '" << treeName << "' not found in "
              << inputFile << std::endl;
    return;
  }

  unsigned int run = 0;
  unsigned int lumi = 0;
  unsigned long long event = 0;
  std::vector<int>* cand_type = nullptr;
  std::vector<std::string>* cand_label = nullptr;
  std::vector<int>* cand_pdgId = nullptr;
  std::vector<float>* cand_mass = nullptr;
  std::vector<float>* cand_pt = nullptr;
  std::vector<float>* cand_eta = nullptr;
  std::vector<float>* cand_phi = nullptr;
  std::vector<float>* cand_y = nullptr;
  std::vector<float>* cand_vx = nullptr;
  std::vector<float>* cand_vy = nullptr;
  std::vector<float>* cand_vz = nullptr;
  std::vector<int>* cand_charge = nullptr;
  std::vector<unsigned int>* cand_nDau = nullptr;

  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("lumi", &lumi);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("cand_type", &cand_type);
  tree->SetBranchAddress("cand_label", &cand_label);
  tree->SetBranchAddress("cand_pdgId", &cand_pdgId);
  tree->SetBranchAddress("cand_mass", &cand_mass);
  tree->SetBranchAddress("cand_pt", &cand_pt);
  tree->SetBranchAddress("cand_eta", &cand_eta);
  tree->SetBranchAddress("cand_phi", &cand_phi);
  tree->SetBranchAddress("cand_y", &cand_y);
  tree->SetBranchAddress("cand_vx", &cand_vx);
  tree->SetBranchAddress("cand_vy", &cand_vy);
  tree->SetBranchAddress("cand_vz", &cand_vz);
  tree->SetBranchAddress("cand_charge", &cand_charge);
  tree->SetBranchAddress("cand_nDau", &cand_nDau);

  // Optional daughter information
  bool hasDaughters = (tree->GetBranch("cand_dauStart") != nullptr);
  std::vector<unsigned int>* cand_dauStart = nullptr;
  std::vector<unsigned int>* cand_dauCount = nullptr;
  std::vector<float>* dau_pt = nullptr;
  std::vector<float>* dau_eta = nullptr;
  std::vector<float>* dau_phi = nullptr;
  std::vector<float>* dau_mass = nullptr;
  std::vector<int>* dau_charge = nullptr;
  std::vector<int>* dau_pdgId = nullptr;

  if (hasDaughters) {
    tree->SetBranchAddress("cand_dauStart", &cand_dauStart);
    tree->SetBranchAddress("cand_dauCount", &cand_dauCount);
    tree->SetBranchAddress("dau_pt", &dau_pt);
    tree->SetBranchAddress("dau_eta", &dau_eta);
    tree->SetBranchAddress("dau_phi", &dau_phi);
    tree->SetBranchAddress("dau_mass", &dau_mass);
    tree->SetBranchAddress("dau_charge", &dau_charge);
    tree->SetBranchAddress("dau_pdgId", &dau_pdgId);
  }

  TFile outFile(outputFile, "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "[flattenChiCNtuple] Failed to open output file " << outputFile << std::endl;
    return;
  }

  // Candidate-level flattened tree
  TTree flat("ChiC", "Flattened ChiC candidates");
  unsigned int out_run = 0;
  unsigned int out_lumi = 0;
  unsigned long long out_event = 0;
  int out_type = -1;
  std::string out_label;
  int out_pdgId = 0;
  float out_mass = 0;
  float out_pt = 0;
  float out_eta = 0;
  float out_phi = 0;
  float out_y = 0;
  float out_vx = 0;
  float out_vy = 0;
  float out_vz = 0;
  int out_charge = 0;
  unsigned int out_nDau = 0;
  unsigned int out_candIndex = 0;

  flat.Branch("run", &out_run, "run/i");
  flat.Branch("lumi", &out_lumi, "lumi/i");
  flat.Branch("event", &out_event, "event/l");
  flat.Branch("candIndex", &out_candIndex, "candIndex/i");
  flat.Branch("cand_type", &out_type, "cand_type/I");
  flat.Branch("cand_label", &out_label);
  flat.Branch("cand_pdgId", &out_pdgId, "cand_pdgId/I");
  flat.Branch("cand_mass", &out_mass, "cand_mass/F");
  flat.Branch("cand_pt", &out_pt, "cand_pt/F");
  flat.Branch("cand_eta", &out_eta, "cand_eta/F");
  flat.Branch("cand_phi", &out_phi, "cand_phi/F");
  flat.Branch("cand_y", &out_y, "cand_y/F");
  flat.Branch("cand_vx", &out_vx, "cand_vx/F");
  flat.Branch("cand_vy", &out_vy, "cand_vy/F");
  flat.Branch("cand_vz", &out_vz, "cand_vz/F");
  flat.Branch("cand_charge", &out_charge, "cand_charge/I");
  flat.Branch("cand_nDau", &out_nDau, "cand_nDau/i");

  // Optional daughter-level tree
  TTree* daughterTree = nullptr;
  unsigned int dau_index = 0;
  float dau_pt_out = 0;
  float dau_eta_out = 0;
  float dau_phi_out = 0;
  float dau_mass_out = 0;
  int dau_charge_out = 0;
  int dau_pdgId_out = 0;

  if (hasDaughters) {
    daughterTree = new TTree("ChiC_dau", "ChiC daughter table");
    daughterTree->Branch("run", &out_run, "run/i");
    daughterTree->Branch("event", &out_event, "event/l");
    daughterTree->Branch("candIndex", &out_candIndex, "candIndex/i");
    daughterTree->Branch("dauIndex", &dau_index, "dauIndex/i");
    daughterTree->Branch("dau_pt", &dau_pt_out, "dau_pt/F");
    daughterTree->Branch("dau_eta", &dau_eta_out, "dau_eta/F");
    daughterTree->Branch("dau_phi", &dau_phi_out, "dau_phi/F");
    daughterTree->Branch("dau_mass", &dau_mass_out, "dau_mass/F");
    daughterTree->Branch("dau_charge", &dau_charge_out, "dau_charge/I");
    daughterTree->Branch("dau_pdgId", &dau_pdgId_out, "dau_pdgId/I");
  }

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry);

    out_run = run;
    out_lumi = lumi;
    out_event = event;

    const size_t nCand = cand_type ? cand_type->size() : 0;
    for (size_t icand = 0; icand < nCand; ++icand) {
      out_candIndex = static_cast<unsigned int>(icand);
      out_type = cand_type->at(icand);
      out_label = cand_label ? cand_label->at(icand) : std::string("" );
      out_pdgId = cand_pdgId ? cand_pdgId->at(icand) : 0;
      out_mass = cand_mass ? cand_mass->at(icand) : 0.f;
      out_pt = cand_pt ? cand_pt->at(icand) : 0.f;
      out_eta = cand_eta ? cand_eta->at(icand) : 0.f;
      out_phi = cand_phi ? cand_phi->at(icand) : 0.f;
      out_y = cand_y ? cand_y->at(icand) : 0.f;
      out_vx = cand_vx ? cand_vx->at(icand) : 0.f;
      out_vy = cand_vy ? cand_vy->at(icand) : 0.f;
      out_vz = cand_vz ? cand_vz->at(icand) : 0.f;
      out_charge = cand_charge ? cand_charge->at(icand) : 0;
      out_nDau = cand_nDau ? cand_nDau->at(icand) : 0;

      flat.Fill();

      if (hasDaughters && cand_dauStart && cand_dauCount) {
        const unsigned int start = cand_dauStart->at(icand);
        const unsigned int count = cand_dauCount->at(icand);
        for (unsigned int idau = 0; idau < count; ++idau) {
          const unsigned int offset = start + idau;
          dau_index = idau;
          dau_pt_out = dau_pt ? dau_pt->at(offset) : 0.f;
          dau_eta_out = dau_eta ? dau_eta->at(offset) : 0.f;
          dau_phi_out = dau_phi ? dau_phi->at(offset) : 0.f;
          dau_mass_out = dau_mass ? dau_mass->at(offset) : 0.f;
          dau_charge_out = dau_charge ? dau_charge->at(offset) : 0;
          dau_pdgId_out = dau_pdgId ? dau_pdgId->at(offset) : 0;
          daughterTree->Fill();
        }
      }
    }
  }

  outFile.Write();
  outFile.Close();
  inFile.Close();
}
