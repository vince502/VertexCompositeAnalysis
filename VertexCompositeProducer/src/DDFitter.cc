// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DDFitter
// 
/**\class DDFitter DDFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/DDFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DDFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// for DCA
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

const float piMassDD = 0.13957018;
const float piMassDDSquared = piMassDD*piMassDD;
const float dStarMassDD = 2.010000;
float piMassDD_sigma = 3.5E-7f;
float DDMassD0_sigma = 1.6E-4f;
float dStarMassDD_sigma = dStarMassDD*1.e-6;


// Constructor and (empty) destructor
DDFitter::DDFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_d0cand = iC.consumes<reco::VertexCompositeCandidateCollection>(theParameters.getParameter<edm::InputTag>("d0Collection"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));

  // Second, initialize post-fit cuts
  mPiKCutMin = theParameters.getParameter<double>(string("mPiKCutMin"));
  mPiKCutMax = theParameters.getParameter<double>(string("mPiKCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkPtErrCut = theParameters.getParameter<double>(string("tkPtErrCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
  tkPtSumCut = theParameters.getParameter<double>(string("tkPtSumCut"));
  tkEtaDiffCut = theParameters.getParameter<double>(string("tkEtaDiffCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut2D = theParameters.getParameter<double>(string("collinearityCut2D"));
  collinCut3D = theParameters.getParameter<double>(string("collinearityCut3D"));
  dStarMassCut = theParameters.getParameter<double>(string("dStarMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPtCut = theParameters.getParameter<double>(string("dPtCut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  alpha2DCut = theParameters.getParameter<double>(string("alpha2DCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));


  useAnyMVA_ = false;
  forestLabel_ = "D0InpPb";
  std::string type = "BDT";
  useForestFromDB_ = true;
  dbFileName_ = "";

  forest_ = nullptr;

  if(theParameters.exists("useAnyMVA")) useAnyMVA_ = theParameters.getParameter<bool>("useAnyMVA");

  if(useAnyMVA_){
    if(theParameters.exists("mvaType"))type = theParameters.getParameter<std::string>("mvaType");
    if(theParameters.exists("GBRForestLabel"))forestLabel_ = theParameters.getParameter<std::string>("GBRForestLabel");
    if(theParameters.exists("GBRForestFileName")){
      dbFileName_ = theParameters.getParameter<std::string>("GBRForestFileName");
      useForestFromDB_ = false;
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeProducer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;
  }

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

DDFitter::~DDFitter() {
  delete forest_;
}

// Method containing the algorithm for vertex reconstruction
void DDFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;
  using namespace std; 

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs;
  std::vector<TransientTrack> theTransTracks;
  std::vector<pat::GenericParticleRef> theD0CandRefs;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::VertexCompositeCandidateCollection> theD0Handle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_d0cand, theD0Handle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
  iEvent.getByToken(token_dedx, dEdxHandle);

  std::cout << "Hi " << std::endl;

  if( !theTrackHandle->size() ) return;

  std::cout << "Bye " << std::endl;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  magField = bFieldHandle.product();

  //needed for IP error
  AnalyticalImpactPointExtrapolator extrapolator(magField);
  TrajectoryStateOnSurface tsos;

  // Setup TMVA
//  mvaValValueMap = auto_ptr<edm::ValueMap<float> >(new edm::ValueMap<float>);
//  edm::ValueMap<float>::Filler mvaFiller(*mvaValValueMap);

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2)
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	      if (tmpRef->quality(qualities[ndx_])){
	        quality_ok = true;
	        break;          
	      }
      }
    }
    if( !quality_ok ) continue;

    if( tmpRef->normalizedChi2() < tkChi2Cut &&
        tmpRef->numberOfValidHits() >= tkNhitsCut &&
        tmpRef->ptError() / tmpRef->pt() < tkPtErrCut &&
        tmpRef->pt() > tkPtCut && fabs(tmpRef->eta()) < tkEtaCut ) {
      TransientTrack tmpTk( *tmpRef, magField );

      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }
  // for(unsigned int idx = 0; idx < theD0Handle->size(); idx ++){
  //   pat::GenericParticleRef tmpRef ( theD0Handle, idx);
  //   theD0CandRefs.push_back( tmpRef );
  // }

  //float posCandMass[2] = {piMassDD, kaonMassD0};
  //float negCandMass[2] = {kaonMassD0, piMassDD};
  //float posCandMass_sigma[2] = {piMassDD_sigma, kaonMassD0_sigma};
  //float negCandMass_sigma[2] = {kaonMassD0_sigma, piMassDD_sigma};
  //int   pdg_id[2] = {421, -421};

  // Loop over tracks and vertex good charged track pairs
std::cout << "Loop1" << std::endl;
  for(unsigned int didx1 = 0; didx1 < theD0Handle->size(); didx1++) {

std::cout << "Loop2" << std::endl;
    for(unsigned int didx2 = didx1 + 1; didx2 < theD0Handle->size(); didx2++) {

      // Not using this on Dstar fit (1)
      // if( (theTrackRefs[didx1]->pt() + theTrackRefs[trdx2]->pt()) < tkPtSumCut) continue;
      // if( abs(theTrackRefs[didx1]->eta() - theTrackRefs[trdx2]->eta()) > tkEtaDiffCut) continue;

      //This vector holds the 3 tracks (K + pi) +pi to be vertexed
      // std::vector<TransientTrack> transTracks;

      // TrackRef pionTrackRef = theTrackRefs[trdx1];
      // TransientTrack* pionTransTkPtr = 0;
      // pionTransTkPtr = &theTransTracks[trdx1];
      VertexCompositeCandidate theD01 = (*theD0Handle)[didx1];
      VertexCompositeCandidate theD02 = (*theD0Handle)[didx2];

      // if( !pionTransTkPtr->impactPointStateAvailable()) continue;
      // const auto& D0Vec1 = theD01.p4();
      // const auto& D0Vec2 = theD02.p4();
      // const reco::Track& thePiTrack = pionTransTkPtr->track();
      // math::PtEtaPhiMLorentzVector pPi(thePiTrack.pt(), thePiTrack.eta(), thePiTrack.phi(), piMassDD);
      // double theDDcandMass = (D0Vec1 + D0Vec2).M();
      // if(theDDcandMass - D0Vec.M() >0.16) continue;


      

      // Calculate DCA of two daughters
//      double dzvtx_pos = positiveTrackRef->dz(bestvtx);
//      double dxyvtx_pos = positiveTrackRef->dxy(bestvtx);
//      double dzerror_pos = sqrt(positiveTrackRef->dzError()*positiveTrackRef->dzError()+zVtxError*zVtxError);
//      double dxyerror_pos = sqrt(positiveTrackRef->d0Error()*positiveTrackRef->d0Error()+xVtxError*yVtxError);
//      double dauLongImpactSig_pos = dzvtx_pos/dzerror_pos;
//      double dauTransImpactSig_pos = dxyvtx_pos/dxyerror_pos;
//
//      double dzvtx_neg = negativeTrackRef->dz(bestvtx);
//      double dxyvtx_neg = negativeTrackRef->dxy(bestvtx);
//      double dzerror_neg = sqrt(negativeTrackRef->dzError()*negativeTrackRef->dzError()+zVtxError*zVtxError);
//      double dxyerror_neg = sqrt(negativeTrackRef->d0Error()*negativeTrackRef->d0Error()+xVtxError*yVtxError);
//      double dauLongImpactSig_neg = dzvtx_neg/dzerror_neg;
//      double dauTransImpactSig_neg = dxyvtx_neg/dxyerror_neg;
//
//      double nhits_pos = positiveTrackRef->numberOfValidHits();
//      double nhits_neg = negativeTrackRef->numberOfValidHits(); 
//    
//      double ptErr_pos = positiveTrackRef->ptError();
//      double ptErr_neg = negativeTrackRef->ptError();
//
//      double dedx_pos=-999.;
//      double dedx_neg=-999.;
//      // Extract dEdx
//      if(dEdxHandle.isValid()){
//        const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
//        dedx_pos = dEdxTrack[positiveTrackRef].dEdx();
//        dedx_neg = dEdxTrack[negativeTrackRef].dEdx();
//      } 
//      dedx_pos = dedx_pos;
//      dedx_neg = dedx_neg;

//      // Fill the vector of TransientTracks to send to KVF
//      transTracks.push_back(*posTransTkPtr);
//      transTracks.push_back(*negTransTkPtr);

      // Trajectory states to calculate DCA for the 2 tracks
//      FreeTrajectoryState posState = posTransTkPtr->impactPointTSCP().theState();
//      FreeTrajectoryState negState = pionTransTkPtr->impactPointTSCP().theState();
//
//      if( !posTransTkPtr->impactPointTSCP().isValid() || !negTransTkPtr->impactPointTSCP().isValid() ) continue;
//
//      // Measure distance between tracks at their closest approach
//      ClosestApproachInRPhi cApp;
//      cApp.calculate(posState, negState);
//      if( !cApp.status() ) continue;
//      float dca = fabs( cApp.distance() );
//      GlobalPoint cxPt = cApp.crossingPoint();
//
//      if (dca < 0. || dca > tkDCACut) continue;
//
//      // Get trajectory states for the tracks at POCA for later cuts
//      TrajectoryStateClosestToPoint posTSCP = posTransTkPtr->trajectoryStateClosestToPoint( cxPt );
//      TrajectoryStateClosestToPoint negTSCP = negTransTkPtr->trajectoryStateClosestToPoint( cxPt );
//
//      if( !posTSCP.isValid() || !negTSCP.isValid() ) continue;
//
//      if( (mass1 > mPiKCutMax || mass1 < mPiKCutMin) && (mass2 > mPiKCutMax || mass2 < mPiKCutMin)) continue;
//      if( totalPt < dPtCut ) continue;


       float chi = 0.0;
       float ndf = 0.0;

       //Creating a KinematicParticleFactory
       KinematicParticleFactoryFromTransientTrack pFactory;
       vector<RefCountedKinematicParticle> d01Daus;
       vector<RefCountedKinematicParticle> d02Daus;
       reco::Candidate* dau10 = theD01.daughter(0);
       reco::Candidate* dau11 = theD01.daughter(1);
       reco::TransientTrack ttk10(*dau10->bestTrack(), magField);
       reco::TransientTrack ttk11(*dau11->bestTrack(), magField);
       float dau10mass =  dau10->mass();
       float dau11mass =  dau11->mass();
       d01Daus.push_back(pFactory.particle(ttk10,dau10mass,chi,ndf,DDMassD0_sigma));
       d01Daus.push_back(pFactory.particle(ttk11,dau11mass,chi,ndf,DDMassD0_sigma));

       reco::Candidate* dau20 = theD02.daughter(0);
       reco::Candidate* dau21 = theD02.daughter(1);
       reco::TransientTrack ttk20(*dau20->bestTrack(), magField);
       reco::TransientTrack ttk21(*dau21->bestTrack(), magField);
       float dau20mass =  dau20->mass();
       float dau21mass =  dau21->mass();
       d02Daus.push_back(pFactory.particle(ttk20,dau20mass,chi,ndf,DDMassD0_sigma));
       d02Daus.push_back(pFactory.particle(ttk21,dau21mass,chi,ndf,DDMassD0_sigma));

       if( ttk10 == ttk20 || ttk10 == ttk21){
std::cout << "Duplicate daughter, continue" << std::endl;
 continue;}
       if( ttk11 == ttk20 || ttk11 == ttk21){
std::cout << "Duplicate daughter, continue" << std::endl;
 continue;}

       KinematicParticleVertexFitter kpvFitter;
       RefCountedKinematicTree d01Tree =  kpvFitter.fit(d01Daus);
       RefCountedKinematicTree d02Tree =  kpvFitter.fit(d02Daus);

       d01Tree->movePointerToTheTop();
       d02Tree->movePointerToTheTop();

       vector<RefCountedKinematicParticle> ddParticles;
       ddParticles.push_back(d01Tree->currentParticle());
       ddParticles.push_back(d02Tree->currentParticle());

       KinematicParticleVertexFitter ddFitter;
       RefCountedKinematicTree ddVertex;
       ddVertex = ddFitter.fit(ddParticles);

       if( !ddVertex->isValid() ) continue;

       ddVertex->movePointerToTheTop();
       RefCountedKinematicParticle ddCand = ddVertex->currentParticle();
       if (!ddCand->currentState().isValid()) continue;

       RefCountedKinematicVertex ddDecayVertex = ddVertex->currentDecayVertex();
       if (!ddDecayVertex->vertexIsValid()) continue;

	     float ddC2Prob = TMath::Prob(ddDecayVertex->chiSquared(),ddDecayVertex->degreesOfFreedom());
	     if (ddC2Prob < VtxChiProbCut) continue;

       ddVertex->movePointerToTheFirstChild();
       RefCountedKinematicParticle posCand = ddVertex->currentParticle();
       ddVertex->movePointerToTheNextChild();
       RefCountedKinematicParticle negCand = ddVertex->currentParticle();

       if(!posCand->currentState().isValid() || !negCand->currentState().isValid()) continue;

       KinematicParameters posCandKP = posCand->currentState().kinematicParameters();
       KinematicParameters negCandKP = negCand->currentState().kinematicParameters();

       GlobalVector ddTotalP = GlobalVector (ddCand->currentState().globalMomentum().x(),
                                                ddCand->currentState().globalMomentum().y(),
                                                ddCand->currentState().globalMomentum().z());

       GlobalVector posCandTotalP = GlobalVector(posCandKP.momentum().x(),posCandKP.momentum().y(),posCandKP.momentum().z());
       GlobalVector negCandTotalP = GlobalVector(negCandKP.momentum().x(),negCandKP.momentum().y(),negCandKP.momentum().z());

       float posCandTotalE = sqrt( posCandTotalP.mag2() + theD01.mass()*theD01.mass() );
       float negCandTotalE = sqrt( negCandTotalP.mag2() + theD02.mass()*theD02.mass() );
       float ddTotalE = posCandTotalE + negCandTotalE;

       const Particle::LorentzVector ddP4(ddTotalP.x(), ddTotalP.y(), ddTotalP.z(), ddTotalE);

       Particle::Point ddVtx((*ddDecayVertex).position().x(), (*ddDecayVertex).position().y(), (*ddDecayVertex).position().z());
       std::vector<double> ddVtxEVec;
       ddVtxEVec.push_back( ddDecayVertex->error().cxx() );
       ddVtxEVec.push_back( ddDecayVertex->error().cyx() );
       ddVtxEVec.push_back( ddDecayVertex->error().cyy() );
       ddVtxEVec.push_back( ddDecayVertex->error().czx() );
       ddVtxEVec.push_back( ddDecayVertex->error().czy() );
       ddVtxEVec.push_back( ddDecayVertex->error().czz() );
       SMatrixSym3D ddVtxCovMatrix(ddVtxEVec.begin(), ddVtxEVec.end());
       const Vertex::CovarianceMatrix ddVtxCov(ddVtxCovMatrix);
       double ddVtxChi2(ddDecayVertex->chiSquared());
       double ddVtxNdof(ddDecayVertex->degreesOfFreedom());
       double ddNormalizedChi2 = ddVtxChi2/ddVtxNdof;

       double rVtxMag = 99999.0; 
       double lVtxMag = 99999.0;
       double sigmaRvtxMag = 999.0;
       double sigmaLvtxMag = 999.0;
       double ddAngle3D = -100.0;
       double ddAngle2D = -100.0;

       GlobalVector ddLineOfFlight = GlobalVector (ddVtx.x() - xVtx,
                                                   ddVtx.y() - yVtx,
                                                   ddVtx.z() - zVtx);

       SMatrixSym3D ddTotalCov;
       if(isVtxPV) ddTotalCov = ddVtxCovMatrix + vtxPrimary->covariance();
       else ddTotalCov = ddVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

       SVector3 distanceVector3D(ddLineOfFlight.x(), ddLineOfFlight.y(), ddLineOfFlight.z());
       SVector3 distanceVector2D(ddLineOfFlight.x(), ddLineOfFlight.y(), 0.0);

       ddAngle3D = angle(ddLineOfFlight.x(), ddLineOfFlight.y(), ddLineOfFlight.z(),
                       ddTotalP.x(), ddTotalP.y(), ddTotalP.z());
       ddAngle2D = angle(ddLineOfFlight.x(), ddLineOfFlight.y(), (float)0.0,
                       ddTotalP.x(), ddTotalP.y(), (float)0.0);

       lVtxMag = ddLineOfFlight.mag();
       rVtxMag = ddLineOfFlight.perp();
       sigmaLvtxMag = sqrt(ROOT::Math::Similarity(ddTotalCov, distanceVector3D)) / lVtxMag;
       sigmaRvtxMag = sqrt(ROOT::Math::Similarity(ddTotalCov, distanceVector2D)) / rVtxMag;

       // DCA error
       tsos = extrapolator.extrapolate(ddCand->currentState().freeTrajectoryState(), RecoVertex::convertPos(vtxPrimary->position()));
       Measurement1D cur3DIP;
       VertexDistance3D a3d;
       GlobalPoint refPoint          = tsos.globalPosition();
       GlobalError refPointErr       = tsos.cartesianError().position();
       GlobalPoint vertexPosition    = RecoVertex::convertPos(vtxPrimary->position());
       GlobalError vertexPositionErr = RecoVertex::convertError(vtxPrimary->error());
       cur3DIP =  (a3d.distance(VertexState(vertexPosition,vertexPositionErr), VertexState(refPoint, refPointErr)));

       if( ddNormalizedChi2 > chi2Cut ||
           rVtxMag < rVtxCut ||
           rVtxMag / sigmaRvtxMag < rVtxSigCut ||
           lVtxMag < lVtxCut ||
           lVtxMag / sigmaLvtxMag < lVtxSigCut ||
           cos(ddAngle3D) < collinCut3D || cos(ddAngle2D) < collinCut2D || ddAngle3D > alphaCut || ddAngle2D > alpha2DCut
       ) continue;

       VertexCompositeCandidate* theDD = 0;
       theDD = new VertexCompositeCandidate(0, ddP4, ddVtx, ddVtxCov, ddVtxChi2, ddVtxNdof);

      //  RecoChargedCandidate
      //    theNegCand(theTrackRefs[trdx1]->charge(), Particle::LorentzVector(negCandTotalP.x(),
      //                                             negCandTotalP.y(), negCandTotalP.z(),
      //                                             negCandTotalE), ddVtx);
      //  theNegCand.setTrack(pionTrackRef);

       AddFourMomenta addp4;
       theDD->addDaughter(theD01);
       theDD->addDaughter(theD02);
      //  int pdgId = (int) theTrackRefs[trdx1]->charge() * 413;
      int pdgId = 30443;
       theDD->setPdgId(pdgId);
       addp4.set( *theDD );
       //if( theDD->mass() < ddMassDD + ddMassCut &&
       //    theDD->mass() > ddMassDD - ddMassCut ) 
       //{
         theDDs.push_back( *theDD );
         dcaVals_.push_back(cur3DIP.value());
         dcaErrs_.push_back(cur3DIP.error());

// per//form MVA evaluation
         if(useAnyMVA_)
         {
      //    //   float gbrVals_[20];
      //    //   gbrVals_[0] = d0P4.Pt();
      //    //   gbrVals_[1] = d0P4.Eta();
      //    //   gbrVals_[2] = d0C2Prob;
      //    //   gbrVals_[3] = lVtxMag / sigmaLvtxMag;
      //    //   gbrVals_[4] = rVtxMag / sigmaRvtxMag;
      //    //   gbrVals_[5] = lVtxMag;
      //    //   gbrVals_[6] = d0Angle3D;
      //    //   gbrVals_[7] = d0Angle2D;
      //    //   gbrVals_[8] = dauLongImpactSig_pos;
      //    //   gbrVals_[9] = dauLongImpactSig_neg;
      //    //   gbrVals_[10] = dauTransImpactSig_pos;
      //    //   gbrVals_[11] = dauTransImpactSig_neg;
      //    //   gbrVals_[12] = nhits_pos;
      //    //   gbrVals_[13] = nhits_neg;
      //    //   gbrVals_[14] = ptErr_pos;
      //    //   gbrVals_[15] = ptErr_neg;
      //    //   gbrVals_[16] = posCandTotalP.perp();
      //    //   gbrVals_[17] = negCandTotalP.perp();
      //    //   gbrVals_[18] = posCandTotalP.eta();
      //    //   gbrVals_[19] = negCandTotalP.eta();

      //    //   GBRForest const * forest = forest_;
      //    //   if(useForestFromDB_){
      //    //     edm::ESHandle<GBRForest> forestHandle;
      //    //     iSetup.get<GBRWrapperRcd>().get(forestLabel_,forestHandle);
      //    //     forest = forestHandle.product();
      //    //   }

      //    //   auto gbrVal = forest->GetClassifier(gbrVals_);
      //    //   mvaVals_.push_back(gbrVal);
         }
       //}

       if(theDD) delete theDD;
      }
  }

//  mvaFiller.insert(theDDs,mvaVals_.begin(),mvaVals_.end());
//  mvaFiller.fill();
//  mvas = std::make_unique<MVACollection>(mvaVals_.begin(),mvaVals_.end());

}
// Get methods

const reco::VertexCompositeCandidateCollection& DDFitter::getDD() const {
  return theDDs;
}

const std::vector<float>& DDFitter::getDCAVals() const{
  return dcaVals_;
}

const std::vector<float>& DDFitter::getDCAErrs() const{
  return dcaErrs_;
}

const std::vector<float>& DDFitter::getMVAVals() const {
  return mvaVals_;
}

/*
auto_ptr<edm::ValueMap<float> > DDFitter::getMVAMap() const {
  return mvaValValueMap;
}
*/

void DDFitter::resetAll() {
    theDDs.clear();
    mvaVals_.clear();
    dcaVals_.clear();
    dcaErrs_.clear();
}
