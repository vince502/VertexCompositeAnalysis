// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      DStar5PFitter
// 
/**\class DStar5PFitter DStar5PFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/DStar5PFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DStar5PFitter.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/commonTools.h"
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

using CC = DStar5PFitter::CC;
using CCC = DStar5PFitter::CCC;

static const float piMassDStar = 0.13957018;
static const float piMassDStarSquared = piMassDStar*piMassDStar;
static const float dStarMassDStar = 2.010000;
static float piMassDStar_sigma = 3.5E-7f;
static float D0MassD0_sigma = 1.6E-4f;
static float dStarMassDStar_sigma = dStarMassDStar*1.e-6;


// Constructor and (empty) destructor
DStar5PFitter::DStar5PFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) :
    bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_d0cand = iC.consumes<CCC>(theParameters.getParameter<edm::InputTag>("d0Collection"));
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
  combineAllTracks = theParameters.getParameter<bool>(string("combineAllTracks"));


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

DStar5PFitter::~DStar5PFitter() {
  delete forest_;
}

// Method containing the algorithm for vertex reconstruction
void DStar5PFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
  Handle<CCC> theD0Handle;
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


  if( !theTrackHandle->size() ) return;
  bFieldHandle = iSetup.getHandle(bField_esToken_);

  magField = bFieldHandle.product();

  //needed for IP error
  AnalyticalImpactPointExtrapolator extrapolator(magField);
  TrajectoryStateOnSurface tsos;

  // Setup TMVA
//  mvaValValueMap = auto_ptr<edm::ValueMap<float> >(new edm::ValueMap<float>);
//  edm::ValueMap<float>::Filler mvaFiller(*mvaValValueMap);

  // Use commonTools to get best vertex for reference
  using namespace VertexCompositeProducerCommonTools;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  
  // Get first valid vertex for reference (used for legacy compatibility)
  auto [bestvtx, vtxIdx] = getBestVertex(vtxCollection, *theBeamSpotHandle, 5);
  bool isVtxPV = (vtxCollection.size() > 0 && vtxIdx < vtxCollection.size());
  
  // Store vertex position and errors for later use
  double xVtx = bestvtx.x();
  double yVtx = bestvtx.y();
  double zVtx = bestvtx.z();
  double xVtxError = 0.0, yVtxError = 0.0, zVtxError = 0.0;
  const reco::Vertex* vtxPrimary = nullptr;
  
  if (isVtxPV) {
    vtxPrimary = &(vtxCollection[vtxIdx]);
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  } else {
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }

  // Vector to store best vertex index for each track
  std::vector<unsigned int> trackVertexIndices;

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  if(theTrackHandle->size() < 5 ) return;
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

      // Find best vertex for this track
      auto [trackVtx, trackVtxIdx] = getBestVertex(vtxCollection, *theBeamSpotHandle, tmpRef.get(), 
                                                     VertexSelectionCriteria::CLOSEST_DZ, 5);
      
      double dzvtx = tmpRef->dz(trackVtx);
      double dxyvtx = tmpRef->dxy(trackVtx);
      
      double trackZVtxError = 0.0, trackXYVtxError = 0.0;
      if (trackVtxIdx < vtxCollection.size()) {
        const reco::Vertex& trackVertex = vtxCollection[trackVtxIdx];
        trackZVtxError = trackVertex.zError();
        trackXYVtxError = trackVertex.xError() * trackVertex.yError();
      } else {
        trackZVtxError = 0.0;
        trackXYVtxError = theBeamSpotHandle->BeamWidthX() * theBeamSpotHandle->BeamWidthY();
      }
      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+trackZVtxError*trackZVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+trackXYVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
        trackVertexIndices.push_back( trackVtxIdx );
      }
    }
  }
  // for(unsigned int idx = 0; idx < theD0Handle->size(); idx ++){
  //   pat::GenericParticleRef tmpRef ( theD0Handle, idx);
  //   theD0CandRefs.push_back( tmpRef );
  // }

  //float posCandMass[2] = {piMassDStar, kaonMassD0};
  //float negCandMass[2] = {kaonMassD0, piMassDStar};
  //float posCandMass_sigma[2] = {piMassDStar_sigma, kaonMassD0_sigma};
  //float negCandMass_sigma[2] = {kaonMassD0_sigma, piMassDStar_sigma};
  //int   pdg_id[2] = {421, -421};

  // Loop over tracks and vertex good charged track pairs
  for(unsigned int didx1 = 0; didx1 < theD0Handle->size(); didx1++) {

    for(unsigned int trdx1 = 0; trdx1 < theTrackRefs.size(); trdx1++) {

      //This vector holds the 5 tracks (4-prong D0) + pi to be vertexed
      std::vector<TransientTrack> transTracks;

      TrackRef pionTrackRef = theTrackRefs[trdx1];
      TransientTrack* pionTransTkPtr = 0;
      pionTransTkPtr = &theTransTracks[trdx1];
      const CC& theD0 = theD0Handle->at(didx1);

      // Check vertex matching between D0 and slow pion if both have valid vertex association
      if (!combineAllTracks && theD0.hasUserFloat("assocVtxIndex")) {
        float d0VtxIndex = theD0.userFloat("assocVtxIndex");
        float slowPiVtxIndex = static_cast<float>(trackVertexIndices[trdx1]);
        // Skip if D0 has a valid vertex index (>= 0) and it doesn't match the slow pion's
        if (d0VtxIndex >= 0.0f && std::abs(d0VtxIndex - slowPiVtxIndex) > 0.1f) continue;
      }

      // if( !pionTransTkPtr->impactPointStateAvailable()) continue;
      const auto& D0Vec = theD0.p4();
      const reco::Track& thePiTrack = pionTransTkPtr->track();
      math::PtEtaPhiMLorentzVector pPi(thePiTrack.pt(), thePiTrack.eta(), thePiTrack.phi(), piMassDStar);
      double theDStarcandMass = (D0Vec + pPi).M();
      // std::cout << "D* - D0 mass : " << theDStarcandMass << ", " << D0Vec.M() << std::endl;
      if(theDStarcandMass - D0Vec.M() >0.16) continue;
      
       float chi = 0.0;
       float ndf = 0.0;

       //Creating a KinematicParticleFactory
       KinematicParticleFactoryFromTransientTrack pFactory;
       vector<RefCountedKinematicParticle> d0Daus;
       const reco::Candidate* dau0 = theD0.daughter(0);
       const reco::Candidate* dau1 = theD0.daughter(1);
       const reco::Candidate* dau2 = theD0.daughter(2);
       const reco::Candidate* dau3 = theD0.daughter(3);
       if(!dau0 || !dau1 || !dau2 || !dau3) continue;
       const reco::Track* trk0 = dau0->bestTrack();
       const reco::Track* trk1 = dau1->bestTrack();
       const reco::Track* trk2 = dau2->bestTrack();
       const reco::Track* trk3 = dau3->bestTrack();
       if(!trk0 || !trk1 || !trk2 || !trk3) continue;
       reco::TransientTrack ttk0(*trk0, magField);
       reco::TransientTrack ttk1(*trk1, magField);
       reco::TransientTrack ttk2(*trk2, magField);
       reco::TransientTrack ttk3(*trk3, magField);
       if(fabs(thePiTrack.eta() - trk0->eta()) < 0.03) continue;
       if(fabs(thePiTrack.eta() - trk1->eta()) < 0.03) continue;
       if(fabs(thePiTrack.eta() - trk2->eta()) < 0.03) continue;
       if(fabs(thePiTrack.eta() - trk3->eta()) < 0.03) continue;
       d0Daus.push_back(pFactory.particle(ttk0,dau0->mass(),chi,ndf,D0MassD0_sigma));
       d0Daus.push_back(pFactory.particle(ttk1,dau1->mass(),chi,ndf,D0MassD0_sigma));
       d0Daus.push_back(pFactory.particle(ttk2,dau2->mass(),chi,ndf,D0MassD0_sigma));
       d0Daus.push_back(pFactory.particle(ttk3,dau3->mass(),chi,ndf,D0MassD0_sigma));

       KinematicParticleVertexFitter kpvFitter;
       RefCountedKinematicTree d0Tree =  kpvFitter.fit(d0Daus);
      if( !d0Tree->isValid() ) continue;

       d0Tree->movePointerToTheTop();

       vector<RefCountedKinematicParticle> dStarParticles;
       dStarParticles.push_back(d0Tree->currentParticle());
       dStarParticles.push_back(pFactory.particle(*pionTransTkPtr,piMassDStar,chi,ndf,piMassDStar_sigma));

       KinematicParticleVertexFitter dStarFitter;
       RefCountedKinematicTree dStarVertex;
       dStarVertex = dStarFitter.fit(dStarParticles);
       if( !dStarVertex->isValid() ) continue;

       dStarVertex->movePointerToTheTop();
       RefCountedKinematicParticle dStarCand = dStarVertex->currentParticle();
       if (!dStarCand->currentState().isValid()) continue;

       RefCountedKinematicVertex dStarDecayVertex = dStarVertex->currentDecayVertex();
       if (!dStarDecayVertex->vertexIsValid()) continue;

	     float dStarC2Prob = TMath::Prob(dStarDecayVertex->chiSquared(),dStarDecayVertex->degreesOfFreedom());
	     if (dStarC2Prob < VtxChiProbCut) continue;

       dStarVertex->movePointerToTheFirstChild();
       RefCountedKinematicParticle posCand = dStarVertex->currentParticle();
       dStarVertex->movePointerToTheNextChild();
       RefCountedKinematicParticle negCand = dStarVertex->currentParticle();

       if(!posCand->currentState().isValid() || !negCand->currentState().isValid()) continue;

       KinematicParameters posCandKP = posCand->currentState().kinematicParameters();
       KinematicParameters negCandKP = negCand->currentState().kinematicParameters();

       GlobalVector dStarTotalP = GlobalVector (dStarCand->currentState().globalMomentum().x(),
                                                dStarCand->currentState().globalMomentum().y(),
                                                dStarCand->currentState().globalMomentum().z());

       GlobalVector posCandTotalP = GlobalVector(posCandKP.momentum().x(),posCandKP.momentum().y(),posCandKP.momentum().z());
       GlobalVector negCandTotalP = GlobalVector(negCandKP.momentum().x(),negCandKP.momentum().y(),negCandKP.momentum().z());

       float posCandTotalE = sqrt( posCandTotalP.mag2() + theD0.mass()*theD0.mass() );
       float negCandTotalE = sqrt( negCandTotalP.mag2() + piMassDStar*piMassDStar );
       float dStarTotalE = posCandTotalE + negCandTotalE;

       const Particle::LorentzVector dStarP4(dStarTotalP.x(), dStarTotalP.y(), dStarTotalP.z(), dStarTotalE);

       Particle::Point dStarVtx((*dStarDecayVertex).position().x(), (*dStarDecayVertex).position().y(), (*dStarDecayVertex).position().z());
       std::vector<double> dStarVtxEVec;
       dStarVtxEVec.push_back( dStarDecayVertex->error().cxx() );
       dStarVtxEVec.push_back( dStarDecayVertex->error().cyx() );
       dStarVtxEVec.push_back( dStarDecayVertex->error().cyy() );
       dStarVtxEVec.push_back( dStarDecayVertex->error().czx() );
       dStarVtxEVec.push_back( dStarDecayVertex->error().czy() );
       dStarVtxEVec.push_back( dStarDecayVertex->error().czz() );
       SMatrixSym3D dStarVtxCovMatrix(dStarVtxEVec.begin(), dStarVtxEVec.end());
       const Vertex::CovarianceMatrix dStarVtxCov(dStarVtxCovMatrix);
       double dStarVtxChi2(dStarDecayVertex->chiSquared());
       double dStarVtxNdof(dStarDecayVertex->degreesOfFreedom());
       double dStarNormalizedChi2 = (dStarVtxNdof > 0.) ? dStarVtxChi2 / dStarVtxNdof : -1.;

       double rVtxMag = 99999.0; 
       double lVtxMag = 99999.0;
       double sigmaRvtxMag = 999.0;
       double sigmaLvtxMag = 999.0;
       double dStarAngle3D = -100.0;
       double dStarAngle2D = -100.0;

       GlobalVector dStarLineOfFlight = GlobalVector (dStarVtx.x() - xVtx,
                                                   dStarVtx.y() - yVtx,
                                                   dStarVtx.z() - zVtx);

       SMatrixSym3D dStarTotalCov;
       if(isVtxPV) dStarTotalCov = dStarVtxCovMatrix + vtxPrimary->covariance();
       else dStarTotalCov = dStarVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

       SVector3 distanceVector3D(dStarLineOfFlight.x(), dStarLineOfFlight.y(), dStarLineOfFlight.z());
       SVector3 distanceVector2D(dStarLineOfFlight.x(), dStarLineOfFlight.y(), 0.0);

       dStarAngle3D = angle(dStarLineOfFlight.x(), dStarLineOfFlight.y(), dStarLineOfFlight.z(),
                       dStarTotalP.x(), dStarTotalP.y(), dStarTotalP.z());
       dStarAngle2D = angle(dStarLineOfFlight.x(), dStarLineOfFlight.y(), (float)0.0,
                       dStarTotalP.x(), dStarTotalP.y(), (float)0.0);

       lVtxMag = dStarLineOfFlight.mag();
       rVtxMag = dStarLineOfFlight.perp();
       sigmaLvtxMag = sqrt(ROOT::Math::Similarity(dStarTotalCov, distanceVector3D)) / lVtxMag;
       sigmaRvtxMag = sqrt(ROOT::Math::Similarity(dStarTotalCov, distanceVector2D)) / rVtxMag;

       // DCA error
       tsos = extrapolator.extrapolate(dStarCand->currentState().freeTrajectoryState(), RecoVertex::convertPos(vtxPrimary->position()));
       Measurement1D cur3DIP;
       VertexDistance3D a3d;
       GlobalPoint refPoint          = tsos.globalPosition();
       GlobalError refPointErr       = tsos.cartesianError().position();
       GlobalPoint vertexPosition    = RecoVertex::convertPos(vtxPrimary->position());
       GlobalError vertexPositionErr = RecoVertex::convertError(vtxPrimary->error());
       cur3DIP =  (a3d.distance(VertexState(vertexPosition,vertexPositionErr), VertexState(refPoint, refPointErr)));

       if( dStarNormalizedChi2 > chi2Cut ||
           rVtxMag < rVtxCut ||
           rVtxMag / sigmaRvtxMag < rVtxSigCut ||
           lVtxMag < lVtxCut ||
           lVtxMag / sigmaLvtxMag < lVtxSigCut ||
           cos(dStarAngle3D) < collinCut3D || cos(dStarAngle2D) < collinCut2D || dStarAngle3D > alphaCut || dStarAngle2D > alpha2DCut
       ) continue;

       auto theDStar = std::make_unique<CC>();
       const int charge = theTrackRefs[trdx1]->charge();
       theDStar->setP4(dStarP4);
       theDStar->setCharge(charge);
       theDStar->setPdgId(charge * 413);
       theDStar->setVertex(reco::Candidate::Point(dStarVtx.x(), dStarVtx.y(), dStarVtx.z()));

       RecoChargedCandidate
         theNegCand(charge, Particle::LorentzVector(negCandTotalP.x(),
                                                  negCandTotalP.y(), negCandTotalP.z(),
                                                  negCandTotalE), dStarVtx);
       theNegCand.setTrack(pionTrackRef);

       AddFourMomenta addp4;
       theDStar->addDaughter(theD0, "D0");
       theDStar->addDaughter(theNegCand, "slowPi");
       addp4.set(*theDStar);

       if( std::abs(theDStar->mass() - dStarMassDStar) <= dStarMassCut )
       {
         const math::XYZPoint dStarXYZ(dStarVtx.x(), dStarVtx.y(), dStarVtx.z());
         reco::Vertex dStarVtxObj(dStarXYZ, dStarVtxCov, dStarVtxChi2, dStarVtxNdof, theDStar->numberOfDaughters());
         theDStar->addUserData("Vtx", dStarVtxObj);
         theDStar->addUserFloat("VtxChi2", dStarVtxChi2);
         theDStar->addUserFloat("VtxNdof", dStarVtxNdof);
         theDStar->addUserFloat("vertexChi2", dStarVtxChi2);
         theDStar->addUserFloat("vertexNdof", dStarVtxNdof);
         theDStar->addUserFloat("vertexNormalizedChi2", dStarNormalizedChi2);
         theDStar->addUserFloat("alpha3D", dStarAngle3D);
         theDStar->addUserFloat("alpha2D", dStarAngle2D);
         theDStar->addUserFloat("decaylength3D", lVtxMag);
         theDStar->addUserFloat("decaylength2D", rVtxMag);
         theDStar->addUserFloat("decaylengthsignif3D", (sigmaLvtxMag > 0.) ? lVtxMag / sigmaLvtxMag : -1.f);
         theDStar->addUserFloat("decaylengthsignif2D", (sigmaRvtxMag > 0.) ? rVtxMag / sigmaRvtxMag : -1.f);
         theDStar->addUserFloat("dca3D", cur3DIP.value());
         theDStar->addUserFloat("dca3DErr", cur3DIP.error());
         theDStar->addUserFloat("deltaM", theDStar->mass() - D0Vec.M());

         theDStars.emplace_back(std::move(*theDStar));
         dcaVals_.push_back(cur3DIP.value());
         dcaErrs_.push_back(cur3DIP.error());
         detlaM_.push_back( theDStars.back().mass() - D0Vec.M());

         if(useAnyMVA_)
         {
           mvaVals_.push_back(0.f);
         }
       }
      }
  }

//  mvaFiller.insert(theDStars,mvaVals_.begin(),mvaVals_.end());
//  mvaFiller.fill();
//  mvas = std::make_unique<MVACollection>(mvaVals_.begin(),mvaVals_.end());

}
// Get methods

const DStar5PFitter::CCC& DStar5PFitter::getDStar() const {
  return theDStars;
}

const std::vector<float>& DStar5PFitter::getDCAVals() const{
  return dcaVals_;
}

const std::vector<float>& DStar5PFitter::getDCAErrs() const{
  return dcaErrs_;
}

const std::vector<float>& DStar5PFitter::getMVAVals() const {
  return mvaVals_;
}

const std::vector<float>& DStar5PFitter::getDeltaM() const {
  return detlaM_;
}

/*
auto_ptr<edm::ValueMap<float> > DStar5PFitter::getMVAMap() const {
  return mvaValValueMap;
}
*/

void DStar5PFitter::resetAll() {
    theDStars.clear();
    mvaVals_.clear();
    dcaVals_.clear();
    dcaErrs_.clear();
    detlaM_.clear();
}
