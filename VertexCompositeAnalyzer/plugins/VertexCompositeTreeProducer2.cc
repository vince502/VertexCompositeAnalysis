
// #include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
// #include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
// #include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

#include "VertexCompositeTreeProducer2.hxx"

// constants, enums and typedefs
//
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>
    SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

// ------------ method called to for each event  ------------
void VertexCompositeTreeProducer2::analyze(const edm::Event &iEvent,
                                           const edm::EventSetup &iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;

  if (doGenNtuple_)
    fillGEN(iEvent, iSetup);
  if (doRecoNtuple_)
    fillRECO(iEvent, iSetup);

  if (saveTree_)
    VertexCompositeNtuple->Fill();
};

void VertexCompositeTreeProducer2::fillRECO(const edm::Event &iEvent,
                                            const edm::EventSetup &iSetup) {
  // get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tok_generalTrk_, tracks);

  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
  iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,
                    v0candidates);
  const reco::VertexCompositeCandidateCollection *v0candidates_ =
      v0candidates.product();

  edm::Handle<MVACollection> mvavalues;
  edm::Handle<MVACollection> mvavalues2;
  if (useAnyMVA_) {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    assert((*mvavalues).size() == v0candidates->size());
    if (doubleCand_) {
      iEvent.getByToken(MVAValues_Token2_, mvavalues2);
      assert((*mvavalues2).size() == v0candidates->size());
    }
  }
  edm::Handle<std::vector<float>> dcaValues;
  edm::Handle<std::vector<float>> dcaErrors;
  if (useDCA_) {
    iEvent.getByToken(tok_DCAVal_, dcaValues);
    iEvent.getByToken(tok_DCAErr_, dcaErrors);
    assert((*dcaValues).size() == v0candidates->size());
    assert((*dcaErrors).size() == v0candidates->size());
  }

  edm::Handle<reco::GenParticleCollection> genpars;
  if (doGenMatching_ || doGenMatchingTOF_)
    iEvent.getByToken(tok_genParticle_, genpars);

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle1;
  iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle2;
  iEvent.getByToken(Dedx_Token2_, dEdxHandle2);

  centrality = -1;
  if (isCentrality_) {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);

    iEvent.getByToken(tok_centBinLabel_, cbin_);
    centrality = *cbin_;

    HFsumET = cent->EtHFtowerSum();
    Npixel = cent->multiplicityPixel();
    //      int ntrk = cent->Ntracks();
  }

  if (isEventPlane_) {
    edm::Handle<reco::EvtPlaneCollection> eventplanes;
    iEvent.getByToken(tok_eventplaneSrc_, eventplanes);

    const reco::EvtPlane &ephfp1 = (*eventplanes)[0];
    const reco::EvtPlane &ephfm1 = (*eventplanes)[1];
    const reco::EvtPlane &ephfp2 = (*eventplanes)[6];
    const reco::EvtPlane &ephfm2 = (*eventplanes)[7];
    const reco::EvtPlane &ephfp3 = (*eventplanes)[13];
    const reco::EvtPlane &ephfm3 = (*eventplanes)[14];

    ephfpAngle[0] = ephfp1.angle(2);
    ephfpAngle[1] = ephfp2.angle(2);
    ephfpAngle[2] = ephfp3.angle(2);

    ephfmAngle[0] = ephfm1.angle(2);
    ephfmAngle[1] = ephfm2.angle(2);
    ephfmAngle[2] = ephfm3.angle(2);

    ephfpQ[0] = ephfp1.q(2);
    ephfpQ[1] = ephfp2.q(2);
    ephfpQ[2] = ephfp3.q(2);

    ephfmQ[0] = ephfm1.q(2);
    ephfmQ[1] = ephfm2.q(2);
    ephfmQ[2] = ephfm3.q(2);

    ephfpSumW = ephfp2.sumw();
    ephfmSumW = ephfm2.sumw();
  }

  // best vertex
  bestvz = -999.9;
  bestvx = -999.9;
  bestvy = -999.9;
  double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
  const reco::Vertex &vtx = (*vertices)[0];
  bestvz = vtx.z();
  bestvx = vtx.x();
  bestvy = vtx.y();
  bestvzError = vtx.zError();
  bestvxError = vtx.xError();
  bestvyError = vtx.yError();

  // Ntrkoffline
  Ntrkoffline = 0;
  if (multMax_ != -1 && multMin_ != -1) {
    for (unsigned it = 0; it < tracks->size(); ++it) {

      const reco::Track &trk = (*tracks)[it];

      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzvtx = trk.dz(bestvtx);
      double dxyvtx = trk.dxy(bestvtx);
      double dzerror =
          sqrt(trk.dzError() * trk.dzError() + bestvzError * bestvzError);
      double dxyerror =
          sqrt(trk.d0Error() * trk.d0Error() + bestvxError * bestvyError);

      if (!trk.quality(reco::TrackBase::highPurity))
        continue;
      if (fabs(trk.ptError()) / trk.pt() > 0.10)
        continue;
      if (fabs(dzvtx / dzerror) > 3)
        continue;
      if (fabs(dxyvtx / dxyerror) > 3)
        continue;

      double eta = trk.eta();
      double pt = trk.pt();

      if (fabs(eta) > 2.4)
        continue;
      if (pt <= 0.4)
        continue;
      Ntrkoffline++;
    }
  }

  // Gen info for matching
  std::vector<reco::GenParticleRef> genRefs;
  if (doGenMatching_) {
    if (!genpars.isValid()) {
      cout << "Gen matching cannot be done without Gen collection!!" << endl;
      return;
    }
#ifdef DEBUG
    cout << "[RECO-GEN] Gen particle collection size : " << genpars->size();
#endif
    for (unsigned int it = 0; it < genpars->size(); ++it) {
      const reco::GenParticle &trk = (*genpars)[it];
      int id = trk.pdgId();
      if (fabs(id) != PID_)
        continue; // check is target
      if (decayInGen_ && trk.numberOfDaughters() != 2 && !threeProngDecay_)
        continue; // check 2-pron decay if target decays in Gen
      if (decayInGen_ && trk.numberOfDaughters() != 3 && threeProngDecay_)
        continue; // check 2-pron decay if target decays in Gen

      // wrong when considering two layer decay
      int nDau = threeProngDecay_ ? 3 : 2;
      std::vector<unsigned int> idxs;
      std::vector<unsigned int> permutations(nDau);
      std::iota(permutations.begin(), permutations.end(), 0);
      std::sort(permutations.begin(), permutations.end());
      if (!threeProngDecay_) {
        do {
          auto Dd1 = trk.daughter(permutations.at(0));
          auto Dd2 = trk.daughter(permutations.at(1));
          if (abs(Dd1->pdgId()) == PID_dau1_ &&
              abs(Dd2->pdgId()) == PID_dau2_) {
            idxs = permutations;
            break;
          }
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));
      } else {
        do {
          auto Dd1 = trk.daughter(permutations.at(0));
          auto Dd2 = trk.daughter(permutations.at(1));
          auto Dd3 = trk.daughter(permutations.at(2));

          if (abs(Dd1->pdgId()) == PID_dau1_ &&
              abs(Dd2->pdgId()) == PID_dau2_ &&
              abs(Dd3->pdgId() == PID_dau3_)) {
            idxs = permutations;
            break;
          }
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));
      }
      if (decayInGen_ && idxs.empty())
        continue;
      genRefs.push_back(reco::GenParticleRef(genpars, it));
    }
#ifdef DEBUG
    cout << " Passing : " << genRefs.size() << endl;
    for (unsigned int igen = 0; igen < genRefs.size(); igen++) {
      auto const &theGen = genRefs.at(igen);
      const auto *genDau0 = theGen->daughter(0);
      const auto *genDau1 = theGen->daughter(1);
      cout << Form("GEN COL : Gen pT eta phi %.4f, %.4f, %.4f, trk pt eta "
                   "phi's : (%.4f, %.4f, %.4f) (%.4f, %.4f, %.4f)",
                   theGen->pt(), theGen->eta(), theGen->phi(), genDau0->pt(),
                   genDau0->eta(), genDau0->phi(), genDau1->pt(),
                   genDau1->eta(), genDau1->phi())
           << endl;
    }
#endif
    // if (genRefs.size()>1) std::cout << "More than one target of generated
    // particles\n";
  }

  // RECO Candidate info
  candSize = v0candidates_->size();
  if (debug_)
    std::cout << "Reco cand Size : " << candSize << std::endl;
  for (unsigned it = 0; it < v0candidates_->size(); ++it) {

    const reco::VertexCompositeCandidate &trk = (*v0candidates_)[it];

    double secvz = -999.9, secvx = -999.9, secvy = -999.9;
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();
    // cout << "size pt eta phi mass " << candSize << ", "<< trk.pt() << ", " <<
    // trk.eta() << ", " << trk.phi() << ", MM" << trk.mass() << endl;

    eta[it] = trk.eta();
    y[it] = trk.rapidity();
    pt[it] = trk.pt();
    phi[it] = trk.phi();
    flavor[it] = trk.pdgId() / abs(trk.pdgId());

    mva[it] = 0.0;
    if (useAnyMVA_)
      mva[it] = (*mvavalues)[it];

    dca3D[it] = -1.0;
    dcaErr3D[it] = -1.0;
    if (useDCA_) {
      dca3D[it] = dcaValues->at(it);
      dcaErr3D[it] = dcaErrors->at(it);
    }

    double px = trk.px();
    double py = trk.py();
    double pz = trk.pz();
    mass[it] = trk.mass();

    const reco::Candidate *d1 = trk.daughter(0);
    const reco::Candidate *d2 = trk.daughter(1);
    if (doubleCand_) {
      flavor1[it] = d1->pdgId() / abs(d1->pdgId());
      flavor2[it] = d2->pdgId() / abs(d2->pdgId());
      if (debug_ && d1->pt() == d2->pt())
        std::cout << "Two daughter is same" << std::endl;
    }
    const reco::Candidate *d3 = 0;
    if (threeProngDecay_)
      d3 = trk.daughter(2);

    // Gen match
    if (doGenMatching_) {
      const auto nGen = genRefs.size();
      if (!doGenDoubleDecay_) {
        matchGEN[it] = false;
        isSwap[it] = false;
        idmom_reco1[it] = -77;
        idBAnc_reco1[it] = -77;
        matchToGen1[it] = MAXCAN + 1;

        for (unsigned int igen = 0; igen < nGen; igen++) {
          auto const &theGen = genRefs.at(igen);
          // Only works for 2 body two layer decay

          reco::Candidate const *recoDaus[3] = {nullptr, nullptr, nullptr};
          reco::Candidate const *genDaus[3] = {nullptr, nullptr, nullptr};

          if (matchHadron(&trk, *theGen)) {
            matchGEN1[it] = true;
            matchToGen1[it] = igen;
            isSwap[it] = checkSwap(&trk, *theGen);
            auto mom_ref = findMother(theGen);
            if (mom_ref.isNonnull())
              idmom_reco1[it] = mom_ref->pdgId();
            int __count_anc__ = 0;
            auto __ref_anc__ = mom_ref;
            while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
              __ref_anc__ = findMother(__ref_anc__);
              if (__ref_anc__.isNonnull()) {
                if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                  idBAnc_reco2[it] = __ref_anc__->pdgId();
                }
              }
            }
            gen_D1pT_[it] = theGen->pt();
            gen_D1eta_[it] = theGen->eta();
            gen_D1phi_[it] = theGen->phi();
            gen_D1mass_[it] = theGen->mass();
            gen_D1y_[it] = theGen->rapidity();
            gen_D1charge_[it] = theGen->charge();
            gen_D1pdgId_[it] = theGen->pdgId();

            // all done in genDecayLength
            // gen_decayLength3D_;
            // gen_decayLength2D_;
            // gen_angle3D_;
            // gen_angle2D_;
            genDecayLength(*theGen, gen_D1decayLength2D_[it],
                           gen_D1decayLength3D_[it], gen_D1angle2D_[it],
                           gen_D1angle3D_[it]);
            getAncestorId(*theGen, gen_D1ancestorId_[it],
                          gen_D1ancestorFlavor_[it]);

            const auto *genDau0 = theGen->daughter(0);
            const auto *genDau1 = theGen->daughter(1);

            gen_D1pTD1_[it] = genDau0->pt();
            gen_D1etaD1_[it] = genDau0->eta();
            gen_D1phiD1_[it] = genDau0->phi();
            gen_D1massD1_[it] = genDau0->mass();
            gen_D1yD1_[it] = genDau0->rapidity();
            gen_D1chargeD1_[it] = genDau0->charge();
            gen_D1pdgIdD1_[it] = genDau0->pdgId();

            gen_D1pTD2_[it] = genDau1->pt();
            gen_D1etaD2_[it] = genDau1->eta();
            gen_D1phiD2_[it] = genDau1->phi();
            gen_D1massD2_[it] = genDau1->mass();
            gen_D1yD2_[it] = genDau1->rapidity();
            gen_D1chargeD2_[it] = genDau1->charge();
            gen_D1pdgIdD2_[it] = genDau1->pdgId();
          }
        }
        if (!matchGEN[it]) {
          gen_D1pT_[it] = -99;
          gen_D1eta_[it] = -99;
          gen_D1phi_[it] = -99;
          gen_D1mass_[it] = -99;
          gen_D1y_[it] = -99;
          gen_D1decayLength3D_[it] = -99;
          gen_D1decayLength2D_[it] = -99;
          gen_D1angle3D_[it] = -99;
          gen_D1angle2D_[it] = -99;
          gen_D1pTD1_[it] = -99;
          gen_D1etaD1_[it] = -99;
          gen_D1phiD1_[it] = -99;
          gen_D1massD1_[it] = -99;
          gen_D1yD1_[it] = -99;
          gen_D1pTD2_[it] = -99;
          gen_D1etaD2_[it] = -99;
          gen_D1phiD2_[it] = -99;
          gen_D1massD2_[it] = -99;
          gen_D1yD2_[it] = -99;
        }
      }
      if (doGenDoubleDecay_) {
        matchGEN1[it] = false;
        matchGEN2[it] = false;
        isSwap1[it] = false;
        isSwap2[it] = false;
        idmom_reco1[it] = -77;
        idmom_reco2[it] = -77;
        idBAnc_reco1[it] = -77;
        idBAnc_reco2[it] = -77;
        matchToGen1[it] = MAXCAN + 1;
        matchToGen2[it] = MAXCAN + 1;
        if (debug_)
          std::cout << "nGen : " << nGen << std::endl;
        matchGEN[it] = false;
#ifdef DEBUG
// cout << "nGen : " <<  nGen << endl;
#endif
        for (unsigned int igen = 0; igen < nGen; igen++) {
          auto const &theGen = genRefs.at(igen);
          // Only works for 2 body two layer decay
          // reco::Candidate const* recoD1 = trk.daughter(0);
          // reco::Candidate const* recoD2 = trk.daughter(1);

          const auto nGenDau = theGen->numberOfDaughters();
          if (debug_)
            std::cout << "nGenDau: " << nGenDau << std::endl;

          if (matchHadron(d1, *theGen)) {
            matchGEN1[it] = true;
            matchToGen1[it] = igen;
            isSwap1[it] = checkSwap(d1, *theGen);
            auto mom_ref = findMother(theGen);
            if (mom_ref.isNonnull())
              idmom_reco1[it] = mom_ref->pdgId();
            int __count_anc__ = 0;
            auto __ref_anc__ = mom_ref;
            while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
              __ref_anc__ = findMother(__ref_anc__);
              if (__ref_anc__.isNonnull()) {
                if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                  idBAnc_reco2[it] = __ref_anc__->pdgId();
                }
              }
            }
            gen_D1pT_[it] = theGen->pt();
            gen_D1eta_[it] = theGen->eta();
            gen_D1phi_[it] = theGen->phi();
            gen_D1mass_[it] = theGen->mass();
            gen_D1y_[it] = theGen->rapidity();
            gen_D1charge_[it] = theGen->charge();
            gen_D1pdgId_[it] = theGen->pdgId();

            genDecayLength(*theGen, gen_D1decayLength2D_[it],
                           gen_D1decayLength3D_[it], gen_D1angle2D_[it],
                           gen_D1angle3D_[it]);
            getAncestorId(*theGen, gen_D1ancestorId_[it],
                          gen_D1ancestorFlavor_[it]);

            const auto *genDau0 = theGen->daughter(0);
            const auto *genDau1 = theGen->daughter(1);

            gen_D1pTD1_[it] = genDau0->pt();
            gen_D1etaD1_[it] = genDau0->eta();
            gen_D1phiD1_[it] = genDau0->phi();
            gen_D1massD1_[it] = genDau0->mass();
            gen_D1yD1_[it] = genDau0->rapidity();
            gen_D1chargeD1_[it] = genDau0->charge();
            gen_D1pdgIdD1_[it] = genDau0->pdgId();

            gen_D1pTD2_[it] = genDau1->pt();
            gen_D1etaD2_[it] = genDau1->eta();
            gen_D1phiD2_[it] = genDau1->phi();
            gen_D1massD2_[it] = genDau1->mass();
            gen_D1yD2_[it] = genDau1->rapidity();
            gen_D1chargeD2_[it] = genDau1->charge();
            gen_D1pdgIdD2_[it] = genDau1->pdgId();
          }

          if (matchHadron(d2, *theGen)) {
            matchGEN2[it] = true;
            matchToGen2[it] = igen;
            isSwap2[it] = checkSwap(d2, *theGen);
            auto mom_ref = findMother(theGen);
            if (mom_ref.isNonnull())
              idmom_reco2[it] = mom_ref->pdgId();
            int __count_anc__ = 0;
            auto __ref_anc__ = mom_ref;
            while (__ref_anc__.isNonnull() && __count_anc__ < 50) {
              __ref_anc__ = findMother(__ref_anc__);
              if (__ref_anc__.isNonnull()) {
                if (((int)abs(__ref_anc__->pdgId())) % 1000 / 100 == 5) {
                  idBAnc_reco2[it] = __ref_anc__->pdgId();
                }
              }
            }
            gen_D2pT_[it] = theGen->pt();
            gen_D2eta_[it] = theGen->eta();
            gen_D2phi_[it] = theGen->phi();
            gen_D2mass_[it] = theGen->mass();
            gen_D2y_[it] = theGen->rapidity();
            gen_D2charge_[it] = theGen->charge();
            gen_D2pdgId_[it] = theGen->pdgId();

            genDecayLength(*theGen, gen_D2decayLength2D_[it],
                           gen_D2decayLength3D_[it], gen_D2angle2D_[it],
                           gen_D2angle3D_[it]);
            getAncestorId(*theGen, gen_D2ancestorId_[it],
                          gen_D2ancestorFlavor_[it]);

            const auto *genDau0 = theGen->daughter(0);
            const auto *genDau1 = theGen->daughter(1);

            gen_D2pTD1_[it] = genDau0->pt();
            gen_D2etaD1_[it] = genDau0->eta();
            gen_D2phiD1_[it] = genDau0->phi();
            gen_D2massD1_[it] = genDau0->mass();
            gen_D2yD1_[it] = genDau0->rapidity();
            gen_D2chargeD1_[it] = genDau0->charge();
            gen_D2pdgIdD1_[it] = genDau0->pdgId();

            gen_D2pTD2_[it] = genDau1->pt();
            gen_D2etaD2_[it] = genDau1->eta();
            gen_D2phiD2_[it] = genDau1->phi();
            gen_D2massD2_[it] = genDau1->mass();
            gen_D2yD2_[it] = genDau1->rapidity();
            gen_D2chargeD2_[it] = genDau1->charge();
            gen_D2pdgIdD2_[it] = genDau1->pdgId();
          }
        } // END for nGen
        if (!matchGEN1[it]) {
          gen_D1pT_[it] = -99;
          gen_D1eta_[it] = -99;
          gen_D1phi_[it] = -99;
          gen_D1mass_[it] = -99;
          gen_D1y_[it] = -99;
          gen_D1decayLength3D_[it] = -99;
          gen_D1decayLength2D_[it] = -99;
          gen_D1angle3D_[it] = -99;
          gen_D1angle2D_[it] = -99;
          gen_D1pTD1_[it] = -99;
          gen_D1etaD1_[it] = -99;
          gen_D1phiD1_[it] = -99;
          gen_D1massD1_[it] = -99;
          gen_D1yD1_[it] = -99;
          gen_D1pTD2_[it] = -99;
          gen_D1etaD2_[it] = -99;
          gen_D1phiD2_[it] = -99;
          gen_D1massD2_[it] = -99;
          gen_D1yD2_[it] = -99;
        }
        if (!matchGEN2[it]) {
          gen_D2pT_[it] = -99;
          gen_D2eta_[it] = -99;
          gen_D2phi_[it] = -99;
          gen_D2mass_[it] = -99;
          gen_D2y_[it] = -99;
          gen_D2decayLength3D_[it] = -99;
          gen_D2decayLength2D_[it] = -99;
          gen_D2angle3D_[it] = -99;
          gen_D2angle2D_[it] = -99;
          gen_D2pTD1_[it] = -99;
          gen_D2etaD1_[it] = -99;
          gen_D2phiD1_[it] = -99;
          gen_D2massD1_[it] = -99;
          gen_D2yD1_[it] = -99;
          gen_D2pTD2_[it] = -99;
          gen_D2etaD2_[it] = -99;
          gen_D2phiD2_[it] = -99;
          gen_D2massD2_[it] = -99;
          gen_D2yD2_[it] = -99;
        }
      } // END if doGenDoubleDecay_
      matchGEN[it] = (matchGEN1[it] && matchGEN2[it]);
#ifdef DEBUG
      // if(matchGEN1[it]){
      //   const auto* d1trk1 = d1->daughter(0);
      //   const auto* d1trk2 = d1->daughter(1);
      //   cout << Form("RECO pt eta phi 1 : (%.4f, %.4f, %.4f) -> (%.4f %.4f
      //   %.4f), (%.4f %.4f %.4f)", d1->pt(), d1->eta(), d1->phi(),
      //   d1trk1->pt(), d1trk1->eta(), d1trk1->phi(), d1trk2->pt(),
      //   d1trk2->eta(), d1trk2->phi()) << endl;
      // }
      // if(matchGEN2[it]){
      //   const auto* d2trk1 = d2->daughter(0);
      //   const auto* d2trk2 = d2->daughter(1);
      //   cout << Form("RECO pt eta phi 2 : (%.4f, %.4f, %.4f) -> (%.4f %.4f
      //   %.4f), (%.4f %.4f %.4f)", d2->pt(), d2->eta(), d2->phi(),
      //   d2trk1->pt(), d2trk1->eta(), d2trk1->phi(), d2trk2->pt(),
      //   d2trk2->eta(), d2trk2->phi()) << endl;

      //   cout << Form("Result (match, swap, matchTo) Reco DD %d : Reco 1 (%d,
      //   %d, %d),  Reco 2 (%d, %d, %d)", it, matchGEN1[it], isSwap1[it],
      //   matchToGen1[it], matchGEN2[it], isSwap2[it], matchToGen2[it]) <<
      //   endl;
      // }
      if (matchGEN[it]) {
        const auto *d1trk1 = d1->daughter(0);
        const auto *d1trk2 = d1->daughter(1);
        cout << Form("RECO pt eta phi 1 : (%.4f, %.4f, %.4f) -> (%.4f %.4f "
                     "%.4f), (%.4f %.4f %.4f)",
                     d1->pt(), d1->eta(), d1->phi(), d1trk1->pt(),
                     d1trk1->eta(), d1trk1->phi(), d1trk2->pt(), d1trk2->eta(),
                     d1trk2->phi())
             << endl;
        const auto *d2trk1 = d2->daughter(0);
        const auto *d2trk2 = d2->daughter(1);
        cout << Form("RECO pt eta phi 2 : (%.4f, %.4f, %.4f) -> (%.4f %.4f "
                     "%.4f), (%.4f %.4f %.4f)",
                     d2->pt(), d2->eta(), d2->phi(), d2trk1->pt(),
                     d2trk1->eta(), d2trk1->phi(), d2trk2->pt(), d2trk2->eta(),
                     d2trk2->phi())
             << endl;

        cout << Form("Result (match, swap, matchTo) Reco DD %d : Reco 1 (%d, "
                     "%d, %d),  Reco 2 (%d, %d, %d)",
                     it, matchGEN1[it], isSwap1[it], matchToGen1[it],
                     matchGEN2[it], isSwap2[it], matchToGen2[it])
             << endl;
      }
      // cout << Form("Result (match, swap, matchTo) Reco DD %d : Reco 1 (%d,
      // %d, %d),  Reco 2 (%d, %d, %d)", it, matchGEN1[it], isSwap1[it],
      // matchToGen1[it], matchGEN2[it], isSwap2[it], matchToGen2[it]) << endl;
#endif
    }

    double pxd1 = d1->px();
    double pyd1 = d1->py();
    double pzd1 = d1->pz();
    double pxd2 = d2->px();
    double pyd2 = d2->py();
    double pzd2 = d2->pz();

    TVector3 dauvec1(pxd1, pyd1, pzd1);
    TVector3 dauvec2(pxd2, pyd2, pzd2);

    // pt
    pt1[it] = d1->pt();
    pt2[it] = d2->pt();

    // momentum
    p1[it] = d1->p();
    p2[it] = d2->p();

    // eta
    eta1[it] = d1->eta();
    eta2[it] = d2->eta();

    y1[it] = d1->rapidity();
    y2[it] = d2->rapidity();

    // phi
    phi1[it] = d1->phi();
    phi2[it] = d2->phi();

    // charge
    charge1[it] = d1->charge();
    charge2[it] = d2->charge();

    double pxd3 = -999.9;
    double pyd3 = -999.9;
    double pzd3 = -999.9;
    if (threeProngDecay_ && d3) {
      pxd3 = d3->px();
      pyd3 = d3->py();
      pzd3 = d3->pz();
      pt3[it] = d3->pt();
      p3[it] = d3->p();
      eta3[it] = d3->eta();
      phi3[it] = d3->phi();
      charge3[it] = d3->charge();
    }
    TVector3 dauvec3(pxd3, pyd3, pzd3);

    pid1[it] = -99999;
    pid2[it] = -99999;
    if (doGenMatchingTOF_) {
      for (unsigned it = 0; it < genpars->size(); ++it) {

        const reco::GenParticle &trk = (*genpars)[it];

        if (trk.pt() < 0.001)
          continue;

        int id = trk.pdgId();
        TVector3 trkvect(trk.px(), trk.py(), trk.pz());

        if (fabs(id) != PID_ && trk.charge()) {
          // matching daughter 1
          double deltaR = trkvect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt1[it]) / pt1[it]) < 0.5 &&
              trk.charge() == charge1[it] && pid1[it] == -99999) {
            pid1[it] = id;
          }

          // matching daughter 2
          deltaR = trkvect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt2[it]) / pt2[it]) < 0.5 &&
              trk.charge() == charge2[it] && pid2[it] == -99999) {
            pid2[it] = id;
          }
        }

        if (fabs(id) == PID_ && trk.numberOfDaughters() == 2) {
          const reco::Candidate *Dd1 = trk.daughter(0);
          const reco::Candidate *Dd2 = trk.daughter(1);
          TVector3 d1vect(Dd1->px(), Dd1->py(), Dd1->pz());
          TVector3 d2vect(Dd2->px(), Dd2->py(), Dd2->pz());
          int id1 = Dd1->pdgId();
          int id2 = Dd2->pdgId();

          double deltaR = d1vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt1[it]) / pt1[it]) < 0.5 &&
              Dd1->charge() == charge1[it] && pid1[it] == -99999) {
            pid1[it] = id1;
          }
          deltaR = d2vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt1[it]) / pt1[it]) < 0.5 &&
              Dd2->charge() == charge1[it] && pid1[it] == -99999) {
            pid1[it] = id1;
          }

          deltaR = d1vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt2[it]) / pt2[it]) < 0.5 &&
              Dd1->charge() == charge2[it] && pid2[it] == -99999) {
            pid2[it] = id2;
          }
          deltaR = d2vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt2[it]) / pt2[it]) < 0.5 &&
              Dd2->charge() == charge2[it] && pid2[it] == -99999) {
            pid2[it] = id2;
          }
        }

        if (pid1[it] != -99999 && pid2[it] != -99999)
          break;
      }
    }

    // vtxChi2
    vtxChi2[it] = trk.vertexChi2();
    ndf[it] = trk.vertexNdof();
    VtxProb[it] = TMath::Prob(vtxChi2[it], ndf[it]);

    // PAngle
    TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
    TVector3 secvec(px, py, pz);

    TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
    TVector3 secvec2D(px, py, 0);

    agl[it] = cos(secvec.Angle(ptosvec));
    agl_abs[it] = secvec.Angle(ptosvec);

    agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

    // Decay length 3D
    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>
        SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;
    typedef ROOT::Math::SVector<double, 6> SVector6;

    SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
    SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

    dl[it] = ROOT::Math::Mag(distanceVector);
    dlerror[it] =
        sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl[it];

    dlos[it] = dl[it] / dlerror[it];

    // correct way for both DCA and its Error
    // std::cout << "By cur3DIP                " << dca3D[it] << " +/- " <<
    // dcaErr3D[it] <<"\n"; incorrect way for DCA error std::cout << "By decay
    // length and alpha " << std::sin(agl_abs[it])*dl[it]<< " +/- " <<
    // dlerror[it]* std::sin(agl_abs[it]) <<"\n"; std::cout << "\n";

    // Decay length 2D
    SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1),
                vtx.covariance(1, 1), 0, 0, 0);
    SVector6 v2(trk.vertexCovariance(0, 0), trk.vertexCovariance(0, 1),
                trk.vertexCovariance(1, 1), 0, 0, 0);

    SMatrixSym3D sv1(v1);
    SMatrixSym3D sv2(v2);

    SMatrixSym3D totalCov2D = sv1 + sv2;
    SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

    dl2D[it] = ROOT::Math::Mag(distanceVector2D);
    double dl2Derror =
        sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / dl2D[it];

    dlos2D[it] = dl2D[it] / dl2Derror;

    // trk info
    auto dau1 = d1->get<reco::TrackRef>();
    if (!twoLayerDecay_) {
      // trk quality
      trkquality1[it] = dau1->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx1[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx1[it] = dEdxTrack[dau1].dEdx();
      }

      T4dedx1[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx1[it] = dEdxTrack[dau1].dEdx();
      }

      // track Chi2
      trkChi1[it] = dau1->normalizedChi2();

      // track pT error
      ptErr1[it] = dau1->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit1[it] = dau1->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest1 = dau1->dz(bestvtx);
      double dxybest1 = dau1->dxy(bestvtx);
      double dzerror1 =
          sqrt(dau1->dzError() * dau1->dzError() + bestvzError * bestvzError);
      double dxyerror1 =
          sqrt(dau1->d0Error() * dau1->d0Error() + bestvxError * bestvyError);

      dzos1[it] = dzbest1 / dzerror1;
      dxyos1[it] = dxybest1 / dxyerror1;
    }

    if (!doubleCand_) {
      auto dau2 = d2->get<reco::TrackRef>();

      // trk quality
      trkquality2[it] = dau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx2[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx2[it] = dEdxTrack[dau2].dEdx();
      }

      T4dedx2[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx2[it] = dEdxTrack[dau2].dEdx();
      }

      // track Chi2
      trkChi2[it] = dau2->normalizedChi2();

      // track pT error
      ptErr2[it] = dau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit2[it] = dau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest2 = dau2->dz(bestvtx);
      double dxybest2 = dau2->dxy(bestvtx);
      double dzerror2 =
          sqrt(dau2->dzError() * dau2->dzError() + bestvzError * bestvzError);
      double dxyerror2 =
          sqrt(dau2->d0Error() * dau2->d0Error() + bestvxError * bestvyError);

      dzos2[it] = dzbest2 / dzerror2;
      dxyos2[it] = dxybest2 / dxyerror2;

      if (doMuon_) {
        edm::Handle<reco::MuonCollection> theMuonHandle;
        iEvent.getByToken(tok_muon_, theMuonHandle);

        nmatchedch1[it] = -1;
        nmatchedst1[it] = -1;
        matchedenergy1[it] = -1;
        nmatchedch2[it] = -1;
        nmatchedst2[it] = -1;
        matchedenergy2[it] = -1;

        double x_exp = -999.;
        double y_exp = -999.;
        double xerr_exp = -999.;
        double yerr_exp = -999.;
        double dxdz_exp = -999.;
        double dydz_exp = -999.;
        double dxdzerr_exp = -999.;
        double dydzerr_exp = -999.;

        double x_seg = -999.;
        double y_seg = -999.;
        double xerr_seg = -999.;
        double yerr_seg = -999.;
        double dxdz_seg = -999.;
        double dydz_seg = -999.;
        double dxdzerr_seg = -999.;
        double dydzerr_seg = -999.;

        double dx_seg = 999.;
        double dy_seg = 999.;
        double dxerr_seg = 999.;
        double dyerr_seg = 999.;
        double dxSig_seg = 999.;
        double dySig_seg = 999.;
        double ddxdz_seg = 999.;
        double ddydz_seg = 999.;
        double ddxdzerr_seg = 999.;
        double ddydzerr_seg = 999.;
        double ddxdzSig_seg = 999.;
        double ddydzSig_seg = 999.;

        onestmuon1[it] = false;
        pfmuon1[it] = false;
        glbmuon1[it] = false;
        trkmuon1[it] = false;
        calomuon1[it] = false;
        softmuon1[it] = false;
        onestmuon2[it] = false;
        pfmuon2[it] = false;
        glbmuon2[it] = false;
        trkmuon2[it] = false;
        calomuon2[it] = false;
        softmuon2[it] = false;

        const int muId1 = muAssocToTrack(dau1, theMuonHandle);
        const int muId2 = muAssocToTrack(dau2, theMuonHandle);

        if (muId1 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId1];

          onestmuon1[it] = muon::isGoodMuon(
              cand, muon::selectionTypeFromString("TMOneStationTight"));
          pfmuon1[it] = cand.isPFMuon();
          glbmuon1[it] = cand.isGlobalMuon();
          trkmuon1[it] = cand.isTrackerMuon();
          calomuon1[it] = cand.isCaloMuon();

          if (glbmuon1[it] && trkmuon1[it] &&
              cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                  5 &&
              cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() >
                  0 &&
              fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
              fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
            softmuon1[it] = true;
        }

        if (muId2 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId2];

          onestmuon2[it] = muon::isGoodMuon(
              cand, muon::selectionTypeFromString("TMOneStationTight"));
          pfmuon2[it] = cand.isPFMuon();
          glbmuon2[it] = cand.isGlobalMuon();
          trkmuon2[it] = cand.isTrackerMuon();
          calomuon2[it] = cand.isCaloMuon();

          if (glbmuon2[it] && trkmuon2[it] &&
              cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                  5 &&
              cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() >
                  0 &&
              fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
              fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
            softmuon2[it] = true;
        }

        if (doMuonFull_) {

          if (muId1 != -1) {
            const reco::Muon &cand = (*theMuonHandle)[muId1];

            nmatchedch1[it] = cand.numberOfMatches();
            nmatchedst1[it] = cand.numberOfMatchedStations();

            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy1[it] = muenergy.hadMax;

            const std::vector<reco::MuonChamberMatch> &muchmatches =
                cand.matches();

            for (unsigned int ich = 0; ich < muchmatches.size(); ich++) {
              x_exp = muchmatches[ich].x;
              y_exp = muchmatches[ich].y;
              xerr_exp = muchmatches[ich].xErr;
              yerr_exp = muchmatches[ich].yErr;
              dxdz_exp = muchmatches[ich].dXdZ;
              dydz_exp = muchmatches[ich].dYdZ;
              dxdzerr_exp = muchmatches[ich].dXdZErr;
              dydzerr_exp = muchmatches[ich].dYdZErr;

              std::vector<reco::MuonSegmentMatch> musegmatches =
                  muchmatches[ich].segmentMatches;

              if (!musegmatches.size())
                continue;
              for (unsigned int jseg = 0; jseg < musegmatches.size(); jseg++) {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;

                if (sqrt((x_seg - x_exp) * (x_seg - x_exp) +
                         (y_seg - y_exp) * (y_seg - y_exp)) <
                    sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                  dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg * dxdzerr_seg +
                                      dxdzerr_exp * dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg * dydzerr_seg +
                                      dydzerr_exp * dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }

              dx1_seg_[it] = dx_seg;
              dy1_seg_[it] = dy_seg;
              dxSig1_seg_[it] = dxSig_seg;
              dySig1_seg_[it] = dySig_seg;
              ddxdz1_seg_[it] = ddxdz_seg;
              ddydz1_seg_[it] = ddydz_seg;
              ddxdzSig1_seg_[it] = ddxdzSig_seg;
              ddydzSig1_seg_[it] = ddydzSig_seg;
            }
          }

          if (muId2 != -1) {
            const reco::Muon &cand = (*theMuonHandle)[muId2];

            nmatchedch2[it] = cand.numberOfMatches();
            nmatchedst2[it] = cand.numberOfMatchedStations();

            reco::MuonEnergy muenergy = cand.calEnergy();
            matchedenergy2[it] = muenergy.hadMax;

            const std::vector<reco::MuonChamberMatch> &muchmatches =
                cand.matches();
            for (unsigned int ich = 0; ich < muchmatches.size(); ich++)
            //                        for(unsigned int ich=0;ich<1;ich++)
            {
              x_exp = muchmatches[ich].x;
              y_exp = muchmatches[ich].y;
              xerr_exp = muchmatches[ich].xErr;
              yerr_exp = muchmatches[ich].yErr;
              dxdz_exp = muchmatches[ich].dXdZ;
              dydz_exp = muchmatches[ich].dYdZ;
              dxdzerr_exp = muchmatches[ich].dXdZErr;
              dydzerr_exp = muchmatches[ich].dYdZErr;

              std::vector<reco::MuonSegmentMatch> musegmatches =
                  muchmatches[ich].segmentMatches;

              if (!musegmatches.size())
                continue;
              for (unsigned int jseg = 0; jseg < musegmatches.size(); jseg++) {
                x_seg = musegmatches[jseg].x;
                y_seg = musegmatches[jseg].y;
                xerr_seg = musegmatches[jseg].xErr;
                yerr_seg = musegmatches[jseg].yErr;
                dxdz_seg = musegmatches[jseg].dXdZ;
                dydz_seg = musegmatches[jseg].dYdZ;
                dxdzerr_seg = musegmatches[jseg].dXdZErr;
                dydzerr_seg = musegmatches[jseg].dYdZErr;

                if (sqrt((x_seg - x_exp) * (x_seg - x_exp) +
                         (y_seg - y_exp) * (y_seg - y_exp)) <
                    sqrt(dx_seg * dx_seg + dy_seg * dy_seg)) {
                  dx_seg = x_seg - x_exp;
                  dy_seg = y_seg - y_exp;
                  dxerr_seg = sqrt(xerr_seg * xerr_seg + xerr_exp * xerr_exp);
                  dyerr_seg = sqrt(yerr_seg * yerr_seg + yerr_exp * yerr_exp);
                  dxSig_seg = dx_seg / dxerr_seg;
                  dySig_seg = dy_seg / dyerr_seg;
                  ddxdz_seg = dxdz_seg - dxdz_exp;
                  ddydz_seg = dydz_seg - dydz_exp;
                  ddxdzerr_seg = sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                  ddydzerr_seg = sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                  ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                  ddydzSig_seg = ddydz_seg / ddydzerr_seg;
                }
              }

              dx2_seg_[it] = dx_seg;
              dy2_seg_[it] = dy_seg;
              dxSig2_seg_[it] = dxSig_seg;
              dySig2_seg_[it] = dySig_seg;
              ddxdz2_seg_[it] = ddxdz_seg;
              ddydz2_seg_[it] = ddydz_seg;
              ddxdzSig2_seg_[it] = ddxdzSig_seg;
              ddydzSig2_seg_[it] = ddydzSig_seg;
            }
          }
        } // doMuonFull
      }
    }

    if (twoLayerDecay_) {
      grand_mass[it] = d1->mass();
      mva1[it] = (*mvavalues)[it];

      const reco::Candidate *gd1 = d1->daughter(0);
      const reco::Candidate *gd2 = d1->daughter(1);

      double gpxd1 = gd1->px();
      double gpyd1 = gd1->py();
      double gpzd1 = gd1->pz();
      double gpxd2 = gd2->px();
      double gpyd2 = gd2->py();
      double gpzd2 = gd2->pz();

      TVector3 gdauvec1(gpxd1, gpyd1, gpzd1);
      TVector3 gdauvec2(gpxd2, gpyd2, gpzd2);

      auto gdau1 = gd1->get<reco::TrackRef>();
      auto gdau2 = gd2->get<reco::TrackRef>();

      // trk quality

      grand_trkquality1[it] = gdau1->quality(reco::TrackBase::highPurity);
      grand_trkquality2[it] = gdau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      grand_H2dedx1[it] = -999.9;
      grand_H2dedx2[it] = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        grand_H2dedx1[it] = dEdxTrack[gdau1].dEdx();
        grand_H2dedx2[it] = dEdxTrack[gdau2].dEdx();
      }

      grand_T4dedx1[it] = -999.9;
      grand_T4dedx2[it] = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        grand_T4dedx1[it] = dEdxTrack[gdau1].dEdx();
        grand_T4dedx2[it] = dEdxTrack[gdau2].dEdx();
      }

      // track pt
      grand_pt1[it] = gd1->pt();
      grand_pt2[it] = gd2->pt();

      // track momentum
      grand_p1[it] = gd1->p();
      grand_p2[it] = gd2->p();

      // track eta
      grand_eta1[it] = gd1->eta();
      grand_eta2[it] = gd2->eta();

      grand_phi1[it] = gd1->phi();
      grand_phi2[it] = gd2->phi();

      // track charge
      grand_charge1[it] = gd1->charge();
      grand_charge2[it] = gd2->charge();

      // track Chi2
      grand_trkChi1[it] = gdau1->normalizedChi2();
      grand_trkChi2[it] = gdau2->normalizedChi2();

      // track pT error
      grand_ptErr1[it] = gdau1->ptError();
      grand_ptErr2[it] = gdau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = d1->vz();
      secvx = d1->vx();
      secvy = d1->vy();

      // trkNHits
      grand_nhit1[it] = gdau1->numberOfValidHits();
      grand_nhit2[it] = gdau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double gdzbest1 = gdau1->dz(bestvtx);
      double gdxybest1 = gdau1->dxy(bestvtx);
      double gdzerror1 =
          sqrt(gdau1->dzError() * gdau1->dzError() + bestvzError * bestvzError);
      double gdxyerror1 =
          sqrt(gdau1->d0Error() * gdau1->d0Error() + bestvxError * bestvyError);

      grand_dzos1[it] = gdzbest1 / gdzerror1;
      grand_dxyos1[it] = gdxybest1 / gdxyerror1;

      double gdzbest2 = gdau2->dz(bestvtx);
      double gdxybest2 = gdau2->dxy(bestvtx);
      double gdzerror2 =
          sqrt(gdau2->dzError() * gdau2->dzError() + bestvzError * bestvzError);
      double gdxyerror2 =
          sqrt(gdau2->d0Error() * gdau2->d0Error() + bestvxError * bestvyError);

      grand_dzos2[it] = gdzbest2 / gdzerror2;
      grand_dxyos2[it] = gdxybest2 / gdxyerror2;

      // vtxChi2
      grand_vtxChi2[it] = d1->vertexChi2();
      grand_ndf[it] = d1->vertexNdof();
      grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it], grand_ndf[it]);

      // PAngle
      TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
      TVector3 secvec(d1->px(), d1->py(), d1->pz());

      TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
      TVector3 secvec2D(d1->px(), d1->py(), 0);

      grand_agl[it] = cos(secvec.Angle(ptosvec));
      grand_agl_abs[it] = secvec.Angle(ptosvec);

      grand_agl2D[it] = cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

      // Decay length 3D
      typedef ROOT::Math::SMatrix<double, 3, 3,
                                  ROOT::Math::MatRepSym<double, 3>>
          SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;
      typedef ROOT::Math::SVector<double, 6> SVector6;

      SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
      SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

      grand_dl[it] = ROOT::Math::Mag(distanceVector);
      grand_dlerror[it] =
          sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / grand_dl[it];

      grand_dlos[it] = grand_dl[it] / grand_dlerror[it];

      // Decay length 2D
      SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1),
                  vtx.covariance(1, 1), 0, 0, 0);
      SVector6 v2(d1->vertexCovariance(0, 0), d1->vertexCovariance(0, 1),
                  d1->vertexCovariance(1, 1), 0, 0, 0);

      SMatrixSym3D sv1(v1);
      SMatrixSym3D sv2(v2);

      SMatrixSym3D totalCov2D = sv1 + sv2;
      SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

      double gdl2D = ROOT::Math::Mag(distanceVector2D);
      double gdl2Derror =
          sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / gdl2D;

      grand_dlos2D[it] = gdl2D / gdl2Derror;
      if (doubleCand_) {
        grand_mass2[it] = d2->mass();
        mva2[it] = (*mvavalues2)[it];

        const reco::Candidate *gd21 = d2->daughter(0);
        const reco::Candidate *gd22 = d2->daughter(1);

        double gpxd21 = gd21->px();
        double gpyd21 = gd21->py();
        double gpzd21 = gd21->pz();
        double gpxd22 = gd22->px();
        double gpyd22 = gd22->py();
        double gpzd22 = gd22->pz();

        TVector3 gdauvec21(gpxd21, gpyd21, gpzd21);
        TVector3 gdauvec22(gpxd22, gpyd22, gpzd22);

        auto gdau21 = gd21->get<reco::TrackRef>();
        auto gdau22 = gd22->get<reco::TrackRef>();

        // trk quality

        grand_trkquality21[it] = gdau21->quality(reco::TrackBase::highPurity);
        grand_trkquality22[it] = gdau22->quality(reco::TrackBase::highPurity);

        // trk dEdx
        grand_H2dedx21[it] = -999.9;
        grand_H2dedx22[it] = -999.9;

        if (dEdxHandle1.isValid()) {
          const edm::ValueMap<reco::DeDxData> dEdxTrack =
              *dEdxHandle1.product();
          grand_H2dedx21[it] = dEdxTrack[gdau21].dEdx();
          grand_H2dedx22[it] = dEdxTrack[gdau22].dEdx();
        }

        grand_T4dedx21[it] = -999.9;
        grand_T4dedx22[it] = -999.9;

        if (dEdxHandle2.isValid()) {
          const edm::ValueMap<reco::DeDxData> dEdxTrack =
              *dEdxHandle2.product();
          grand_T4dedx21[it] = dEdxTrack[gdau21].dEdx();
          grand_T4dedx22[it] = dEdxTrack[gdau22].dEdx();
        }

        // track pt
        grand_pt21[it] = gd21->pt();
        grand_pt22[it] = gd22->pt();

        // track momentum
        grand_p21[it] = gd21->p();
        grand_p22[it] = gd22->p();

        // track eta
        grand_eta21[it] = gd21->eta();
        grand_eta22[it] = gd22->eta();

        // track phi
        grand_phi21[it] = gd21->phi();
        grand_phi22[it] = gd22->phi();

        // track charge
        grand_charge21[it] = gd21->charge();
        grand_charge22[it] = gd22->charge();

        // track Chi2
        grand_trkChi21[it] = gdau21->normalizedChi2();
        grand_trkChi22[it] = gdau22->normalizedChi2();

        // track pT error
        grand_ptErr21[it] = gdau21->ptError();
        grand_ptErr22[it] = gdau22->ptError();

        // vertexCovariance 00-xError 11-y 22-z
        secvz = d2->vz();
        secvx = d2->vx();
        secvy = d2->vy();

        // trkNHits
        grand_nhit21[it] = gdau21->numberOfValidHits();
        grand_nhit22[it] = gdau22->numberOfValidHits();

        // DCA
        //  math::XYZPoint bestvtx2(bestvx,bestvy,bestvz);

        double gdzbest21 = gdau21->dz(bestvtx);
        double gdxybest21 = gdau21->dxy(bestvtx);
        double gdzerror21 = sqrt(gdau21->dzError() * gdau21->dzError() +
                                 bestvzError * bestvzError);
        double gdxyerror21 = sqrt(gdau21->d0Error() * gdau21->d0Error() +
                                  bestvxError * bestvyError);

        grand_dzos21[it] = gdzbest21 / gdzerror21;
        grand_dxyos21[it] = gdxybest21 / gdxyerror21;

        double gdzbest22 = gdau22->dz(bestvtx);
        double gdxybest22 = gdau22->dxy(bestvtx);
        double gdzerror22 = sqrt(gdau22->dzError() * gdau22->dzError() +
                                 bestvzError * bestvzError);
        double gdxyerror22 = sqrt(gdau22->d0Error() * gdau22->d0Error() +
                                  bestvxError * bestvyError);

        grand_dzos22[it] = gdzbest22 / gdzerror22;
        grand_dxyos22[it] = gdxybest22 / gdxyerror22;

        // vtxChi2
        grand_vtxChi22[it] = d2->vertexChi2();
        grand_ndf2[it] = d2->vertexNdof();
        grand_VtxProb2[it] = TMath::Prob(grand_vtxChi22[it], grand_ndf2[it]);

        // PAngle
        TVector3 ptosvec2(secvx - bestvx, secvy - bestvy, secvz - bestvz);
        TVector3 secvec2(d2->px(), d2->py(), d2->pz());

        TVector3 ptosvec2D2(secvx - bestvx, secvy - bestvy, 0);
        TVector3 secvec2D2(d2->px(), d2->py(), 0);

        grand_agl2[it] = cos(secvec2.Angle(ptosvec2));
        grand_agl_abs2[it] = secvec2.Angle(ptosvec2);

        grand_agl2D2[it] = cos(secvec2D2.Angle(ptosvec2D2));
        grand_agl2D_abs2[it] = secvec2D2.Angle(ptosvec2D2);

        // Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3,
                                    ROOT::Math::MatRepSym<double, 3>>
            SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;

        SMatrixSym3D totalCov2 = vtx.covariance() + d2->vertexCovariance();
        SVector3 distanceVector2(secvx - bestvx, secvy - bestvy,
                                 secvz - bestvz);

        grand_dl2[it] = ROOT::Math::Mag(distanceVector2);
        grand_dlerror2[it] =
            sqrt(ROOT::Math::Similarity(totalCov2, distanceVector2)) /
            grand_dl2[it];

        grand_dlos2[it] = grand_dl2[it] / grand_dlerror2[it];

        // Decay length 2D
        SVector6 v21(vtx.covariance(0, 0), vtx.covariance(0, 1),
                     vtx.covariance(1, 1), 0, 0, 0);
        SVector6 v22(d2->vertexCovariance(0, 0), d2->vertexCovariance(0, 1),
                     d2->vertexCovariance(1, 1), 0, 0, 0);

        SMatrixSym3D sv21(v1);
        SMatrixSym3D sv22(v2);

        SMatrixSym3D totalCov2D2 = sv21 + sv22;
        SVector3 distanceVector2D2(secvx - bestvx, secvy - bestvy, 0);

        double gdl2D2 = ROOT::Math::Mag(distanceVector2D2);
        double gdl2Derror2 =
            sqrt(ROOT::Math::Similarity(totalCov2D2, distanceVector2D2)) /
            gdl2D2;

        grand_dlos2D2[it] = gdl2D2 / gdl2Derror2;
      }
    }

    if (saveHistogram_) {
      for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++)
        for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
          if (pt[it] < pTBins_[ipt + 1] && pt[it] > pTBins_[ipt] &&
              y[it] < yBins_[iy + 1] && y[it] > yBins_[iy]) {
            hMassVsMVA[iy][ipt]->Fill(mva[it], mass[it]);
            //                h3DDCAVsMVA[iy][ipt]->Fill(mva[it],dl[it]*sin(agl_abs[it]));
            //                h2DDCAVsMVA[iy][ipt]->Fill(mva[it],dl2D[it]*sin(agl2D_abs[it]));

            if (saveAllHistogram_) {
              hpTVsMVA[iy][ipt]->Fill(mva[it], pt[it]);
              hetaVsMVA[iy][ipt]->Fill(mva[it], eta[it]);
              hyVsMVA[iy][ipt]->Fill(mva[it], y[it]);
              hVtxProbVsMVA[iy][ipt]->Fill(mva[it], VtxProb[it]);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl[it]);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl_abs[it]);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D[it]);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D_abs[it]);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos[it]);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl[it]);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it],
                                                             dlos2D[it]);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl2D[it]);
              hzDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it],
                                                             dzos1[it]);
              hxyDCASignificancedaughter1VsMVA[iy][ipt]->Fill(mva[it],
                                                              dxyos1[it]);
              hNHitD1VsMVA[iy][ipt]->Fill(mva[it], nhit1[it]);
              hpTD1VsMVA[iy][ipt]->Fill(mva[it], pt1[it]);
              hpTerrD1VsMVA[iy][ipt]->Fill(mva[it], ptErr1[it] / pt1[it]);
              hEtaD1VsMVA[iy][ipt]->Fill(mva[it], eta1[it]);
              hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva[it], H2dedx1[it]);
              hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1[it], H2dedx1[it]);
              hzDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it],
                                                             dzos2[it]);
              hxyDCASignificancedaughter2VsMVA[iy][ipt]->Fill(mva[it],
                                                              dxyos2[it]);
              hNHitD2VsMVA[iy][ipt]->Fill(mva[it], nhit2[it]);
              hpTD2VsMVA[iy][ipt]->Fill(mva[it], pt2[it]);
              hpTerrD2VsMVA[iy][ipt]->Fill(mva[it], ptErr2[it] / pt2[it]);
              hEtaD2VsMVA[iy][ipt]->Fill(mva[it], eta2[it]);
              hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva[it], H2dedx2[it]);
              hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2[it], H2dedx2[it]);
              if (threeProngDecay_) {
                hzDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it],
                                                               dzos3[it]);
                hxyDCASignificancedaughter3VsMVA[iy][ipt]->Fill(mva[it],
                                                                dxyos3[it]);
                hNHitD3VsMVA[iy][ipt]->Fill(mva[it], nhit3[it]);
                hpTD3VsMVA[iy][ipt]->Fill(mva[it], pt3[it]);
                hpTerrD3VsMVA[iy][ipt]->Fill(mva[it], ptErr3[it] / pt3[it]);
                hEtaD3VsMVA[iy][ipt]->Fill(mva[it], eta3[it]);
                hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva[it], H2dedx3[it]);
                hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1[it], H2dedx3[it]);
              }
            }
          }
        }
    }
  }
}

void VertexCompositeTreeProducer2::fillGEN(const edm::Event &iEvent,
                                           const edm::EventSetup &iSetup) {
  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);

  candSize_gen = 0;
#ifdef DEBUG
  cout << " fillGEN: Passing : " << genpars->size() << endl;
#endif
  for (unsigned it = 0; it < genpars->size(); ++it) {

    const reco::GenParticle &trk = (*genpars)[it];

    const reco::Candidate *Dd1 = trk.daughter(0);
    const reco::Candidate *Dd2 = trk.daughter(1);

    int id = trk.pdgId();
    // if( fabs(id)!=PID_) continue; //check is target
    if (fabs(id) != PID_)
      continue; // check is target
    if (Dd1 == nullptr || Dd2 == nullptr)
      continue; // check is target
    // if( doGenDoubleDecay_ && !(abs(Dd1->pdgId()) == 421 && abs(Dd2->pdgId())
    // == 421) ) continue; //check is target if(!(trk.statusFlags().isLastCopy()
    // && (trk.pdgId()==21 || std::abs(trk.pdgId())<=6))) continue;
    if (debug_)
      std::cout << "pass id, decay dau : " << trk.numberOfDaughters()
                << std::endl;

    if (decayInGen_ &&
        (trk.numberOfDaughters() != 2 && trk.numberOfDaughters() != 3))
      continue; // check 2-pron decay if target decays in Gen

    if (debug_)
      std::cout << "pass decay" << std::endl;

    candSize_gen += 1;
#ifdef DEBUG
    cout << Form("fillGEN : Gen pT eta phi %.4f, %.4f, %.4f, trk pt eta phi's "
                 ": (%.4f, %.4f, %.4f) (%.4f, %.4f, %.4f) PDG ID %d, %d",
                 trk.pt(), trk.eta(), trk.phi(), Dd1->pt(), Dd1->eta(),
                 Dd1->phi(), Dd2->pt(), Dd2->eta(), Dd2->phi(), Dd1->pdgId(),
                 Dd2->pdgId())
         << endl;
#endif

    pt_gen[candSize_gen - 1] = trk.pt();
    mass_gen[candSize_gen - 1] = trk.mass();
    eta_gen[candSize_gen - 1] = trk.eta();
    phi_gen[candSize_gen - 1] = trk.phi();
    status_gen[candSize_gen - 1] = trk.status();
    idself[candSize_gen - 1] = trk.pdgId();
    idmom[candSize_gen - 1] = -77;
    y_gen[candSize_gen - 1] = trk.rapidity();
    ptmom[candSize_gen - 1] = -999.0;
    etamom[candSize_gen - 1] = -999.0;
    phimom[candSize_gen - 1] = -999.0;
    ymom[candSize_gen - 1] = -999.0;
    statusmom[candSize_gen - 1] = -999;
    genDecayLength(trk, dl2D_gen[candSize_gen - 1], dl3D_gen[candSize_gen - 1],
                   angle2D_gen[candSize_gen - 1],
                   angle3D_gen[candSize_gen - 1]);

    if (trk.numberOfMothers() != 0) {
      const reco::Candidate *mom = trk.mother();
      idmom[candSize_gen - 1] = mom->pdgId();
      ptmom[candSize_gen - 1] = mom->pt();
      etamom[candSize_gen - 1] = mom->eta();
      phimom[candSize_gen - 1] = mom->phi();
      ymom[candSize_gen - 1] = mom->rapidity();
      statusmom[candSize_gen - 1] = mom->status();
    }

    if (!decayInGen_)
      continue;

    const reco::Candidate *Dd3 = trk.daughter(2);

    iddau1[candSize_gen - 1] = fabs(Dd1->pdgId());
    iddau2[candSize_gen - 1] = fabs(Dd2->pdgId());
    if (Dd3)
      iddau3[candSize_gen - 1] = fabs(Dd3->pdgId());
    pt_gen1[candSize_gen - 1] = Dd1->pt();
    mass_gen1[candSize_gen - 1] = Dd1->mass();
    eta_gen1[candSize_gen - 1] = Dd1->eta();
    phi_gen1[candSize_gen - 1] = Dd1->phi();
    status_gen1[candSize_gen - 1] = Dd1->status();
    idself1[candSize_gen - 1] = Dd1->pdgId();

    pt_gen2[candSize_gen - 1] = Dd2->pt();
    mass_gen2[candSize_gen - 1] = Dd2->mass();
    eta_gen2[candSize_gen - 1] = Dd2->eta();
    phi_gen2[candSize_gen - 1] = Dd2->phi();
    status_gen2[candSize_gen - 1] = Dd2->status();
    idself2[candSize_gen - 1] = Dd2->pdgId();
  }
}

// ------------ method called once each job just before starting event
// loop  ------------
void VertexCompositeTreeProducer2::beginJob() {
  TH1D::SetDefaultSumw2();

  if (!doRecoNtuple_ && !doGenNtuple_) {
    cout << "No output for either RECO or GEN!! Fix config!!" << endl;
    return;
  }

  if (twoLayerDecay_ && doMuon_) {
    cout << "Muons cannot be coming from two layer decay!! Fix config!!"
         << endl;
    return;
  }

  if (saveHistogram_)
    initHistogram();
  if (saveTree_)
    initTree();
};

int VertexCompositeTreeProducer2::muAssocToTrack(
    const reco::TrackRef &trackref,
    const edm::Handle<reco::MuonCollection> &muonh) const {
  auto muon =
      std::find_if(muonh->cbegin(), muonh->cend(), [&](const reco::Muon &m) {
        return (m.track().isNonnull() && m.track() == trackref);
      });
  return (muon != muonh->cend() ? std::distance(muonh->cbegin(), muon) : -1);
};

reco::GenParticleRef VertexCompositeTreeProducer2::findMother(
    const reco::GenParticleRef &genParRef) {
  if (genParRef.isNull())
    return genParRef;
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId();
  const int pdg_OLD = pdg;
  while (pdg == pdg_OLD && genMomRef->numberOfMothers() > 0) {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  if (pdg == pdg_OLD)
    genMomRef = reco::GenParticleRef();
  return genMomRef;
};

void VertexCompositeTreeProducer2::genDecayLength(
    const reco::GenParticle &gCand, float &gen_decayLength2D_,
    float &gen_decayLength3D_, float &gen_angle2D_, float &gen_angle3D_) {
  gen_decayLength2D_ = -99.;
  gen_decayLength3D_ = -99.;
  gen_angle2D_ = -99;
  gen_angle3D_ = -99;

  if (gCand.numberOfDaughters() == 0 || !gCand.daughter(0))
    return;
  const auto &dauVtx = gCand.daughter(0)->vertex();
  TVector3 ptosvec(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(),
                   dauVtx.Z() - genVertex_.z());
  TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
  gen_angle3D_ = secvec.Angle(ptosvec);
  gen_decayLength3D_ = ptosvec.Mag();
  TVector3 ptosvec2D(dauVtx.X() - genVertex_.x(), dauVtx.Y() - genVertex_.y(),
                     0.0);
  TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
  gen_angle2D_ = secvec2D.Angle(ptosvec2D);
  gen_decayLength2D_ = ptosvec2D.Mag();
};

void VertexCompositeTreeProducer2::getAncestorId(const reco::GenParticle &gCand,
                                                 int &gen_ancestorId_,
                                                 int &gen_ancestorFlavor_) {
  gen_ancestorId_ = 0;
  gen_ancestorFlavor_ = 0;
  for (auto mothers = gCand.motherRefVector(); !mothers.empty();) {
    auto mom = mothers.at(0);
    mothers = mom->motherRefVector();
    gen_ancestorId_ = mom->pdgId();
    const auto idstr = std::to_string(std::abs(gen_ancestorId_));
    gen_ancestorFlavor_ =
        std::stoi(std::string{idstr.begin(), idstr.begin() + 1});
    if (idstr[0] == '5') {
      break;
    }
    if (std::abs(gen_ancestorId_) <= 40)
      break;
  }
};
