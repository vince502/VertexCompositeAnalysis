#include "VertexCompositeNtupleProducer2.hxx"

void VertexCompositeNtupleProducer2::fillRECO(const edm::Event &iEvent,
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
  if (useAnyMVA_) {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    assert((*mvavalues).size() == v0candidates->size());
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

    HFsumETPlus = cent->EtHFtowerSumPlus();
    HFsumETMinus = cent->EtHFtowerSumMinus();
    ZDCPlus = cent->zdcSumPlus();
    ZDCMinus = cent->zdcSumMinus();
    Npixel = cent->multiplicityPixel();
    //  int ntrk = cent->Ntracks();
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
    // if (genRefs.size()>1) std::cout << "More than one target of generated
    // particles\n";
  }

  // RECO Candidate info
  for (unsigned int it = 0; it < v0candidates_->size(); ++it) {

    const reco::VertexCompositeCandidate &trk = (*v0candidates_)[it];

    double secvz = -999.9, secvx = -999.9, secvy = -999.9;
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();

    eta = trk.eta();
    phi = trk.phi();
    y = trk.rapidity();
    pt = trk.pt();
    flavor = trk.pdgId() / abs(trk.pdgId());

    mva = 0.0;
    if (useAnyMVA_)
      mva = (*mvavalues)[it];

    dca3D = -1.0;
    dcaErr3D = -1.0;
    if (useDCA_) {
      dca3D = dcaValues->at(it);
      dcaErr3D = dcaErrors->at(it);
    }

    double px = trk.px();
    double py = trk.py();
    double pz = trk.pz();
    mass = trk.mass();

    const reco::Candidate *d1 = trk.daughter(0);
    const reco::Candidate *d2 = trk.daughter(1);
    const reco::Candidate *d3 = 0;
    if (threeProngDecay_)
      d3 = trk.daughter(2);

    // Gen match
    if (doGenMatching_) {
      matchGEN = false;
      isSwap = false;
      idmom_reco = -77;

      const auto nGen = genRefs.size();
      for (unsigned int igen = 0; igen < nGen; igen++) {
        const auto &genRef = genRefs.at(igen);

        reco::Candidate const *recoDaus[3] = {nullptr, nullptr, nullptr};
        reco::Candidate const *genDaus[3] = {nullptr, nullptr, nullptr};

        const auto nGenDau = genRef->numberOfDaughters();
        std::vector<unsigned int> permutations(nGenDau);
        std::iota(permutations.begin(), permutations.end(), 0);
        std::sort(permutations.begin(), permutations.end());

        do {
          matchGEN = false;
          for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
            genDaus[iDau] = genRef->daughter(permutations.at(iDau));
            recoDaus[iDau] = trk.daughter(iDau);
          }

          for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
            const double dR =
                reco::deltaR(genDaus[iDau]->eta(), genDaus[iDau]->phi(),
                             recoDaus[iDau]->eta(), recoDaus[iDau]->phi());
            const double dPt = abs(genDaus[iDau]->pt() - recoDaus[iDau]->pt()) /
                               recoDaus[iDau]->pt();
            const bool unMatchCharge =
                genDaus[iDau]->charge() != recoDaus[iDau]->charge();
            const bool unMatchDR = dR > deltaR_;
            const bool unMatchDPt = dPt > 0.5;
            matchGEN = matchGEN || unMatchCharge || unMatchDR || unMatchDPt;
          }
          matchGEN = !matchGEN;
          if (matchGEN)
            break;
        } while (
            std::next_permutation(permutations.begin(), permutations.end()));

        for (unsigned int iDau = 0; iDau < nGenDau; ++iDau) {
          const double diffMass =
              abs(genDaus[iDau]->mass() - recoDaus[iDau]->mass());
          isSwap = isSwap || diffMass > 0.01;
        }
        if (matchGEN) {
          auto mom_ref = findMother(genRef);
          if (mom_ref.isNonnull())
            idmom_reco = mom_ref->pdgId();

          gen_pT_ = genRef->pt();
          gen_eta_ = genRef->eta();
          gen_phi_ = genRef->phi();
          gen_mass_ = genRef->mass();
          gen_y_ = genRef->rapidity();
          gen_charge_ = genRef->charge();
          gen_pdgId_ = genRef->pdgId();

          // all done in genDecayLength
          // gen_decayLength3D_;
          // gen_decayLength2D_;
          // gen_angle3D_;
          // gen_angle2D_;
          genDecayLength(*genRef);
          getAncestorId(*genRef);

          gen_pTD1_ = genDaus[0]->pt();
          gen_etaD1_ = genDaus[0]->eta();
          gen_phiD1_ = genDaus[0]->phi();
          gen_massD1_ = genDaus[0]->mass();
          gen_yD1_ = genDaus[0]->rapidity();
          gen_chargeD1_ = genDaus[0]->charge();
          gen_pdgIdD1_ = genDaus[0]->pdgId();

          gen_pTD2_ = genDaus[1]->pt();
          gen_etaD2_ = genDaus[1]->eta();
          gen_phiD2_ = genDaus[1]->phi();
          gen_massD2_ = genDaus[1]->mass();
          gen_yD2_ = genDaus[1]->rapidity();
          gen_chargeD2_ = genDaus[1]->charge();
          gen_pdgIdD2_ = genDaus[1]->pdgId();
        } else {
          gen_pT_ = -99;
          gen_eta_ = -99;
          gen_phi_ = -99;
          gen_mass_ = -99;
          gen_y_ = -99;

          gen_decayLength3D_ = -99;
          gen_decayLength2D_ = -99;
          gen_angle3D_ = -99;
          gen_angle2D_ = -99;

          gen_pTD1_ = -99;
          gen_etaD1_ = -99;
          gen_phiD1_ = -99;
          gen_massD1_ = -99;
          gen_yD1_ = -99;

          gen_pTD2_ = -99;
          gen_etaD2_ = -99;
          gen_phiD2_ = -99;
          gen_massD2_ = -99;
          gen_yD2_ = -99;
        }
        if (matchGEN)
          break;
      }
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
    pt1 = d1->pt();
    pt2 = d2->pt();

    // momentum
    p1 = d1->p();
    p2 = d2->p();

    // eta
    eta1 = d1->eta();
    eta2 = d2->eta();

    // phi
    phi1 = d1->phi();
    phi2 = d2->phi();

    // charge
    charge1 = d1->charge();
    charge2 = d2->charge();

    double pxd3 = -999.9;
    double pyd3 = -999.9;
    double pzd3 = -999.9;
    if (threeProngDecay_ && d3) {
      pxd3 = d3->px();
      pyd3 = d3->py();
      pzd3 = d3->pz();
      pt3 = d3->pt();
      p3 = d3->p();
      eta3 = d3->eta();
      phi3 = d3->phi();
      charge3 = d3->charge();
    }
    TVector3 dauvec3(pxd3, pyd3, pzd3);

    // useless
    pid1 = -99999;
    pid2 = -99999;
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
          if (deltaR < deltaR_ && fabs((trk.pt() - pt1) / pt1) < 0.5 &&
              trk.charge() == charge1 && pid1 == -99999) {
            pid1 = id;
            //      tof1 = ;
          }

          // matching daughter 2
          deltaR = trkvect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((trk.pt() - pt2) / pt2) < 0.5 &&
              trk.charge() == charge2 && pid2 == -99999) {
            pid2 = id;
            //      tof2 = ;
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
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt1) / pt1) < 0.5 &&
              Dd1->charge() == charge1 && pid1 == -99999) {
            pid1 = id1;
          }
          deltaR = d2vect.DeltaR(dauvec1);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt1) / pt1) < 0.5 &&
              Dd2->charge() == charge1 && pid1 == -99999) {
            pid1 = id1;
          }

          deltaR = d1vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd1->pt() - pt2) / pt2) < 0.5 &&
              Dd1->charge() == charge2 && pid2 == -99999) {
            pid2 = id2;
          }
          deltaR = d2vect.DeltaR(dauvec2);
          if (deltaR < deltaR_ && fabs((Dd2->pt() - pt2) / pt2) < 0.5 &&
              Dd2->charge() == charge2 && pid2 == -99999) {
            pid2 = id2;
          }
        }

        if (pid1 != -99999 && pid2 != -99999)
          break;
      }
    }

    // vtxChi2
    vtxChi2 = trk.vertexChi2();
    ndf = trk.vertexNdof();
    VtxProb = TMath::Prob(vtxChi2, ndf);

    // PAngle
    TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
    TVector3 secvec(px, py, pz);

    TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
    TVector3 secvec2D(px, py, 0);

    agl = cos(secvec.Angle(ptosvec));
    agl_abs = secvec.Angle(ptosvec);

    agl2D = cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs = secvec2D.Angle(ptosvec2D);

    // Decay length 3D

    SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
    SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

    dl = ROOT::Math::Mag(distanceVector);
    dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl;

    dlos = dl / dlerror;

    // Decay length 2D
    SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1),
                vtx.covariance(1, 1), 0, 0, 0);
    SVector6 v2(trk.vertexCovariance(0, 0), trk.vertexCovariance(0, 1),
                trk.vertexCovariance(1, 1), 0, 0, 0);

    SMatrixSym3D sv1(v1);
    SMatrixSym3D sv2(v2);

    SMatrixSym3D totalCov2D = sv1 + sv2;
    SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

    dl2D = ROOT::Math::Mag(distanceVector2D);
    double dl2Derror =
        sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / dl2D;

    dlos2D = dl2D / dl2Derror;

    // trk info
    auto dau1 = d1->get<reco::TrackRef>();
    if (!twoLayerDecay_) {
      // trk quality
      trkquality1 = dau1->quality(reco::TrackBase::highPurity);

      // trk dEdx
      H2dedx1 = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        H2dedx1 = dEdxTrack[dau1].dEdx();
      }

      T4dedx1 = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        T4dedx1 = dEdxTrack[dau1].dEdx();
      }

      // track Chi2
      trkChi1 = dau1->normalizedChi2();

      // track pT error
      ptErr1 = dau1->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();

      // trkNHits
      nhit1 = dau1->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double dzbest1 = dau1->dz(bestvtx);
      double dxybest1 = dau1->dxy(bestvtx);
      double dzerror1 =
          sqrt(dau1->dzError() * dau1->dzError() + bestvzError * bestvzError);
      double dxyerror1 =
          sqrt(dau1->d0Error() * dau1->d0Error() + bestvxError * bestvyError);

      dzos1 = dzbest1 / dzerror1;
      dxyos1 = dxybest1 / dxyerror1;
    }

    auto dau2 = d2->get<reco::TrackRef>();

    // trk quality
    trkquality2 = dau2->quality(reco::TrackBase::highPurity);

    // trk dEdx
    H2dedx2 = -999.9;

    if (dEdxHandle1.isValid()) {
      const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
      H2dedx2 = dEdxTrack[dau2].dEdx();
    }

    T4dedx2 = -999.9;

    if (dEdxHandle2.isValid()) {
      const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
      T4dedx2 = dEdxTrack[dau2].dEdx();
    }

    // track Chi2
    trkChi2 = dau2->normalizedChi2();

    // track pT error
    ptErr2 = dau2->ptError();

    // vertexCovariance 00-xError 11-y 22-z
    secvz = trk.vz();
    secvx = trk.vx();
    secvy = trk.vy();

    // trkNHits
    nhit2 = dau2->numberOfValidHits();

    // DCA
    math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

    double dzbest2 = dau2->dz(bestvtx);
    double dxybest2 = dau2->dxy(bestvtx);
    double dzerror2 =
        sqrt(dau2->dzError() * dau2->dzError() + bestvzError * bestvzError);
    double dxyerror2 =
        sqrt(dau2->d0Error() * dau2->d0Error() + bestvxError * bestvyError);

    dzos2 = dzbest2 / dzerror2;
    dxyos2 = dxybest2 / dxyerror2;

    if (doMuon_) {
      edm::Handle<reco::MuonCollection> theMuonHandle;
      iEvent.getByToken(tok_muon_, theMuonHandle);

      nmatchedch1 = -1;
      nmatchedst1 = -1;
      matchedenergy1 = -1;
      nmatchedch2 = -1;
      nmatchedst2 = -1;
      matchedenergy2 = -1;

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

      onestmuon1 = false;
      pfmuon1 = false;
      glbmuon1 = false;
      trkmuon1 = false;
      calomuon1 = false;
      softmuon1 = false;
      onestmuon2 = false;
      pfmuon2 = false;
      glbmuon2 = false;
      trkmuon2 = false;
      calomuon2 = false;
      softmuon2 = false;

      const int muId1 = muAssocToTrack(dau1, theMuonHandle);
      const int muId2 = muAssocToTrack(dau2, theMuonHandle);

      if (muId1 != -1) {
        const reco::Muon &cand = (*theMuonHandle)[muId1];

        onestmuon1 = muon::isGoodMuon(
            cand, muon::selectionTypeFromString("TMOneStationTight"));
        pfmuon1 = cand.isPFMuon();
        glbmuon1 = cand.isGlobalMuon();
        trkmuon1 = cand.isTrackerMuon();
        calomuon1 = cand.isCaloMuon();

        if (
            // glbmuon1 &&
            trkmuon1 &&
            cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                5 &&
            cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
            fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
            fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
          softmuon1 = true;
      }

      if (muId2 != -1) {
        const reco::Muon &cand = (*theMuonHandle)[muId2];

        onestmuon2 = muon::isGoodMuon(
            cand, muon::selectionTypeFromString("TMOneStationTight"));
        pfmuon2 = cand.isPFMuon();
        glbmuon2 = cand.isGlobalMuon();
        trkmuon2 = cand.isTrackerMuon();
        calomuon2 = cand.isCaloMuon();

        if (
            // glbmuon2 &&
            trkmuon2 &&
            cand.innerTrack()->hitPattern().trackerLayersWithMeasurement() >
                5 &&
            cand.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
            fabs(cand.innerTrack()->dxy(vtx.position())) < 0.3 &&
            fabs(cand.innerTrack()->dz(vtx.position())) < 20.)
          softmuon2 = true;
      }

      if (doMuonFull_) {

        if (muId1 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId1];

          nmatchedch1 = cand.numberOfMatches();
          nmatchedst1 = cand.numberOfMatchedStations();

          reco::MuonEnergy muenergy = cand.calEnergy();
          matchedenergy1 = muenergy.hadMax;

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
                ddxdzerr_seg =
                    sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                ddydzerr_seg =
                    sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                ddydzSig_seg = ddydz_seg / ddydzerr_seg;
              }
            }

            dx1_seg_ = dx_seg;
            dy1_seg_ = dy_seg;
            dxSig1_seg_ = dxSig_seg;
            dySig1_seg_ = dySig_seg;
            ddxdz1_seg_ = ddxdz_seg;
            ddydz1_seg_ = ddydz_seg;
            ddxdzSig1_seg_ = ddxdzSig_seg;
            ddydzSig1_seg_ = ddydzSig_seg;
          }
        }

        if (muId2 != -1) {
          const reco::Muon &cand = (*theMuonHandle)[muId2];

          nmatchedch2 = cand.numberOfMatches();
          nmatchedst2 = cand.numberOfMatchedStations();

          reco::MuonEnergy muenergy = cand.calEnergy();
          matchedenergy2 = muenergy.hadMax;

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
                ddxdzerr_seg =
                    sqrt(dxdzerr_seg * dxdzerr_seg + dxdzerr_exp * dxdzerr_exp);
                ddydzerr_seg =
                    sqrt(dydzerr_seg * dydzerr_seg + dydzerr_exp * dydzerr_exp);
                ddxdzSig_seg = ddxdz_seg / ddxdzerr_seg;
                ddydzSig_seg = ddydz_seg / ddydzerr_seg;
              }
            }

            dx2_seg_ = dx_seg;
            dy2_seg_ = dy_seg;
            dxSig2_seg_ = dxSig_seg;
            dySig2_seg_ = dySig_seg;
            ddxdz2_seg_ = ddxdz_seg;
            ddydz2_seg_ = ddydz_seg;
            ddxdzSig2_seg_ = ddxdzSig_seg;
            ddydzSig2_seg_ = ddydzSig_seg;
          }
        }
      } // doMuonFull
    }

    if (twoLayerDecay_) {
      grand_mass = d1->mass();

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

      grand_trkquality1 = gdau1->quality(reco::TrackBase::highPurity);
      grand_trkquality2 = gdau2->quality(reco::TrackBase::highPurity);

      // trk dEdx
      grand_H2dedx1 = -999.9;
      grand_H2dedx2 = -999.9;

      if (dEdxHandle1.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
        grand_H2dedx1 = dEdxTrack[gdau1].dEdx();
        grand_H2dedx2 = dEdxTrack[gdau2].dEdx();
      }

      grand_T4dedx1 = -999.9;
      grand_T4dedx2 = -999.9;

      if (dEdxHandle2.isValid()) {
        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
        grand_T4dedx1 = dEdxTrack[gdau1].dEdx();
        grand_T4dedx2 = dEdxTrack[gdau2].dEdx();
      }

      // track pt
      grand_pt1 = gd1->pt();
      grand_pt2 = gd2->pt();

      // track momentum
      grand_p1 = gd1->p();
      grand_p2 = gd2->p();

      // track eta
      grand_eta1 = gd1->eta();
      grand_eta2 = gd2->eta();

      // track charge
      grand_charge1 = gd1->charge();
      grand_charge2 = gd2->charge();

      // track Chi2
      grand_trkChi1 = gdau1->normalizedChi2();
      grand_trkChi2 = gdau2->normalizedChi2();

      // track pT error
      grand_ptErr1 = gdau1->ptError();
      grand_ptErr2 = gdau2->ptError();

      // vertexCovariance 00-xError 11-y 22-z
      secvz = d1->vz();
      secvx = d1->vx();
      secvy = d1->vy();

      // trkNHits
      grand_nhit1 = gdau1->numberOfValidHits();
      grand_nhit2 = gdau2->numberOfValidHits();

      // DCA
      math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

      double gdzbest1 = gdau1->dz(bestvtx);
      double gdxybest1 = gdau1->dxy(bestvtx);
      double gdzerror1 =
          sqrt(gdau1->dzError() * gdau1->dzError() + bestvzError * bestvzError);
      double gdxyerror1 =
          sqrt(gdau1->d0Error() * gdau1->d0Error() + bestvxError * bestvyError);

      grand_dzos1 = gdzbest1 / gdzerror1;
      grand_dxyos1 = gdxybest1 / gdxyerror1;

      double gdzbest2 = gdau2->dz(bestvtx);
      double gdxybest2 = gdau2->dxy(bestvtx);
      double gdzerror2 =
          sqrt(gdau2->dzError() * gdau2->dzError() + bestvzError * bestvzError);
      double gdxyerror2 =
          sqrt(gdau2->d0Error() * gdau2->d0Error() + bestvxError * bestvyError);

      grand_dzos2 = gdzbest2 / gdzerror2;
      grand_dxyos2 = gdxybest2 / gdxyerror2;

      // vtxChi2
      grand_vtxChi2 = d1->vertexChi2();
      grand_ndf = d1->vertexNdof();
      grand_VtxProb = TMath::Prob(grand_vtxChi2, grand_ndf);

      // PAngle
      TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
      TVector3 secvec(d1->px(), d1->py(), d1->pz());

      TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
      TVector3 secvec2D(d1->px(), d1->py(), 0);

      grand_agl = cos(secvec.Angle(ptosvec));
      grand_agl_abs = secvec.Angle(ptosvec);

      grand_agl2D = cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs = secvec2D.Angle(ptosvec2D);

      // Decay length 3D
      SMatrixSym3D totalCov = vtx.covariance() + d1->vertexCovariance();
      SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

      grand_dl = ROOT::Math::Mag(distanceVector);
      grand_dlerror =
          sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / grand_dl;

      grand_dlos = grand_dl / grand_dlerror;

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

      grand_dlos2D = gdl2D / gdl2Derror;
    }

    if (saveTree_)
      VertexCompositeNtuple->Fill();
    if (saveHistogram_) {
      for (unsigned int ipt = 0; ipt < pTBins_.size() - 1; ipt++)
        for (unsigned int iy = 0; iy < yBins_.size() - 1; iy++) {
          if (pt < pTBins_[ipt + 1] && pt > pTBins_[ipt] &&
              y < yBins_[iy + 1] && y > yBins_[iy]) {
            hMassVsMVA[iy][ipt]->Fill(mva, mass);
            //  h3DDCAVsMVA[iy][ipt]->Fill(mva,dl*sin(agl_abs));
            //  h2DDCAVsMVA[iy][ipt]->Fill(mva,dl2D*sin(agl2D_abs));

            if (saveAllHistogram_) {
              hpTVsMVA[iy][ipt]->Fill(mva, pt);
              hetaVsMVA[iy][ipt]->Fill(mva, eta);
              hyVsMVA[iy][ipt]->Fill(mva, y);
              hVtxProbVsMVA[iy][ipt]->Fill(mva, VtxProb);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl_abs);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva, agl2D_abs);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva, dlos2D);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva, dl2D);
              hzDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva, dzos1);
              hxyDCASignificanceDaugther1VsMVA[iy][ipt]->Fill(mva, dxyos1);
              hNHitD1VsMVA[iy][ipt]->Fill(mva, nhit1);
              hpTD1VsMVA[iy][ipt]->Fill(mva, pt1);
              hpTerrD1VsMVA[iy][ipt]->Fill(mva, ptErr1 / pt1);
              hEtaD1VsMVA[iy][ipt]->Fill(mva, eta1);
              hdedxHarmonic2D1VsMVA[iy][ipt]->Fill(mva, H2dedx1);
              hdedxHarmonic2D1VsP[iy][ipt]->Fill(p1, H2dedx1);
              hzDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva, dzos2);
              hxyDCASignificanceDaugther2VsMVA[iy][ipt]->Fill(mva, dxyos2);
              hNHitD2VsMVA[iy][ipt]->Fill(mva, nhit2);
              hpTD2VsMVA[iy][ipt]->Fill(mva, pt2);
              hpTerrD2VsMVA[iy][ipt]->Fill(mva, ptErr2 / pt2);
              hEtaD2VsMVA[iy][ipt]->Fill(mva, eta2);
              hdedxHarmonic2D2VsMVA[iy][ipt]->Fill(mva, H2dedx2);
              hdedxHarmonic2D2VsP[iy][ipt]->Fill(p2, H2dedx2);
              if (threeProngDecay_) {
                hzDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva, dzos3);
                hxyDCASignificanceDaugther3VsMVA[iy][ipt]->Fill(mva, dxyos3);
                hNHitD3VsMVA[iy][ipt]->Fill(mva, nhit3);
                hpTD3VsMVA[iy][ipt]->Fill(mva, pt3);
                hpTerrD3VsMVA[iy][ipt]->Fill(mva, ptErr3 / pt3);
                hEtaD3VsMVA[iy][ipt]->Fill(mva, eta3);
                hdedxHarmonic2D3VsMVA[iy][ipt]->Fill(mva, H2dedx3);
                hdedxHarmonic2D3VsP[iy][ipt]->Fill(p1, H2dedx3);
              }
            }
          }
        }
    }
  }
}

void VertexCompositeNtupleProducer2::fillGEN(const edm::Event &iEvent,
                                             const edm::EventSetup &iSetup) {
  genVertex_ = reco::Vertex();

  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);
  // generated primary vertex information
  for (const auto &p : *genpars) {
    if (p.statusFlags().isLastCopy() &&
        (p.pdgId() == 21 || std::abs(p.pdgId()) <= 6)) {
      genVertex_ = p.vertex();
      break;
    }
  }
  gen_PVx_ = genVertex_.x();
  gen_PVy_ = genVertex_.y();
  gen_PVz_ = genVertex_.z();

  for (unsigned it = 0; it < genpars->size(); ++it) {
    const reco::GenParticle &trk = (*genpars)[it];
    int id = trk.pdgId();
    if (fabs(id) != PID_)
      continue; // check is target
    if (decayInGen_ && trk.numberOfDaughters() != 2)
      continue; // check 2-pron decay if target decays in Gen

    pt_gen = trk.pt();
    eta_gen = trk.eta();
    phi_gen = trk.phi();
    status_gen = trk.status();
    idmom = -77;
    y_gen = trk.rapidity();
    ptmom = -999.0;
    etamom = -999.0;
    phimom = -999.0;
    ymom = -999.0;
    statusmom = -999;

    if (trk.numberOfMothers() != 0) {
      const reco::Candidate *mom = trk.mother();
      idmom = mom->pdgId();
      ptmom = mom->pt();
      etamom = mom->eta();
      phimom = mom->phi();
      ymom = mom->rapidity();
      statusmom = mom->status();
    }
    if (!decayInGen_)
      continue;

    const reco::Candidate *Dd1 = trk.daughter(0);
    const reco::Candidate *Dd2 = trk.daughter(1);
    const reco::Candidate *Dd3 = trk.daughter(2);

    iddau1 = fabs(Dd1->pdgId());
    iddau2 = fabs(Dd2->pdgId());
    if (Dd3)
      iddau3 = fabs(Dd3->pdgId());

    getAncestorId(trk);
    genDecayLength(trk);

    gen_pT_ = trk.pt();
    gen_eta_ = trk.eta();
    gen_phi_ = trk.phi();
    gen_mass_ = trk.mass();
    gen_y_ = trk.rapidity();
    gen_charge_ = trk.charge();
    gen_pdgId_ = trk.pdgId();

    gen_pTD1_ = Dd1->pt();
    gen_etaD1_ = Dd1->eta();
    gen_phiD1_ = Dd1->phi();
    gen_massD1_ = Dd1->mass();
    gen_yD1_ = Dd1->rapidity();
    gen_chargeD1_ = Dd1->charge();
    gen_pdgIdD1_ = Dd1->pdgId();

    gen_pTD2_ = Dd2->pt();
    gen_etaD2_ = Dd2->eta();
    gen_phiD2_ = Dd2->phi();
    gen_massD2_ = Dd2->mass();
    gen_yD2_ = Dd2->rapidity();
    gen_chargeD2_ = Dd2->charge();
    gen_pdgIdD2_ = Dd2->pdgId();

    genCandidateNtuple->Fill();
  }
}

// ------------ method called once each job just before starting event
// loop  ------------

// ------------ method called once each job just after ending the event
// loop  ------------
