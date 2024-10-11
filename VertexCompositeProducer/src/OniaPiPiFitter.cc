// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      OniapipiFitter
// 
/**\class OniapipiFitter OniapipiFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/OniapipiFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/OniapipiFitter.h"
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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

const float piMassB = 0.13957018;
const float piMassBSquared = piMassB*piMassB;
const float kaonMassB = 0.493677;
const float kaonMassBSquared = kaonMassB*kaonMassB;
const float oniaMass = 3.094;
const float bMassB = 3.872;
float piMassB_sigma = 3.5E-7f;
float kaonMassB_sigma = 1.6E-5f;
float d0MassB_sigma = bMassB*1.e-6;

// Constructor and (empty) destructor
OniapipiFitter::OniapipiFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) :
    bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  using std::string;

  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_onias = iC.consumes<pat::CompositeCandidateCollection>(theParameters.getParameter<edm::InputTag>("dimuons"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));

  // Second, initialize post-fit cuts
  mPiDCutMin = theParameters.getParameter<double>(string("mPiDCutMin"));
  mPiDCutMax = theParameters.getParameter<double>(string("mPiDCutMax"));
  batTkDCACut = theParameters.getParameter<double>(string("batTkDCACut"));
  batTkChi2Cut = theParameters.getParameter<double>(string("batTkChi2Cut"));
  batTkNhitsCut = theParameters.getParameter<int>(string("batTkNhitsCut"));
  batTkPtCut = theParameters.getParameter<double>(string("batTkPtCut"));
  batTkPtErrCut = theParameters.getParameter<double>(string("batTkPtErrCut"));
  batTkEtaCut = theParameters.getParameter<double>(string("batTkEtaCut"));
  bVtxChi2Cut = theParameters.getParameter<double>(string("bVtxChi2Cut"));
  bRVtxCut = theParameters.getParameter<double>(string("bVtx2DCut"));
  bRVtxSigCut = theParameters.getParameter<double>(string("bVtxSignificance2DCut"));
  bLVtxCut = theParameters.getParameter<double>(string("bVtx3DCut"));
  bLVtxSigCut = theParameters.getParameter<double>(string("bVtxSignificance3DCut"));
  bCollinCut2D = theParameters.getParameter<double>(string("bCollinCut2D"));
  bCollinCut3D = theParameters.getParameter<double>(string("bCollinCut3D"));
  bMassCut = theParameters.getParameter<double>(string("bMassCut"));
  batDauTransImpactSigCut = theParameters.getParameter<double>(string("batDauTransImpactSigCut"));
  batDauLongImpactSigCut = theParameters.getParameter<double>(string("batDauLongImpactSigCut"));
  bVtxChiProbCut = theParameters.getParameter<double>(string("bVtxChiProbCut"));
  bPtCut = theParameters.getParameter<double>(string("bPtCut"));
  bAlphaCut = theParameters.getParameter<double>(string("bAlphaCut"));
  bAlpha2DCut = theParameters.getParameter<double>(string("bAlpha2DCut"));
  isWrongSignB = theParameters.getParameter<bool>(string("isWrongSignB"));
}

OniapipiFitter::~OniapipiFitter() {
}

// Method containing the algorithm for vertex reconstruction
void OniapipiFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<pat::CompositeCandidateCollection> theOniaHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_onias, theOniaHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
  iEvent.getByToken(token_dedx, dEdxHandle);

  if( !theTrackHandle->size() ) return;
  bFieldHandle = iSetup.getHandle(bField_esToken_);

  magField = bFieldHandle.product();

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
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2) {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  } else {
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

    if( tmpRef->normalizedChi2() < batTkChi2Cut &&
        tmpRef->numberOfValidHits() >= batTkNhitsCut &&
        tmpRef->ptError() / tmpRef->pt() < batTkPtErrCut &&
        tmpRef->pt() > batTkPtCut && fabs(tmpRef->eta()) < batTkEtaCut ) {
      TransientTrack tmpTk( *tmpRef, magField );

      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > batDauTransImpactSigCut && fabs(dauLongImpactSig) > batDauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

  const pat::CompositeCandidateCollection theOnias = *(theOniaHandle.product());
  for(unsigned it=0; it<theOnias.size(); ++it){

    const pat::CompositeCandidate & theOnia = theOnias[it];

    float massWindow = 2.040;
    // TODO Rename variables
    // if(theOnia.mass() > d0MassB + massWindow || theOnia.mass() < d0MassB - massWindow) continue;
    if(theOnia.mass() > 2.8 || theOnia.mass() < 3.7) continue;

    for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {
      if ( theTrackRefs[trdx].isNull() ) continue;
      for(unsigned int trdx2 = trdx+1; trdx2 < theTrackRefs.size(); trdx2++) {
        if ( theTrackRefs[trdx2].isNull() ) continue;
        if( theTrackRefs[trdx]->charge() == theTrackRefs[trdx2]->charge()) continue;
        TransientTrack* trk1TransTkPtr = &theTransTracks[trdx];
        TransientTrack* trk2TransTkPtr = &theTransTracks[trdx2];

        //Creating a KinematicParticleFactory
        float chi = 0.;
        float ndf = 0.;
        KinematicParticleFactoryFromTransientTrack pFactory;
        vector<RefCountedKinematicParticle> oniaDaus;
        const reco::Candidate* dau0 = theOnia.daughter(0);
        const reco::Candidate* dau1 = theOnia.daughter(1);
        reco::TransientTrack ttk0(*dau0->bestTrack(), magField);
        reco::TransientTrack ttk1(*dau1->bestTrack(), magField);
        float dau0mass =  dau0->mass();
        float dau1mass =  dau1->mass();
        //if(debug_) 
          cout << "mass of two muons : " << dau0mass << ", " <<dau1mass << endl;
        oniaDaus.push_back(pFactory.particle(ttk0,dau0mass,chi,ndf,piMassB_sigma));
        oniaDaus.push_back(pFactory.particle(ttk1,dau1mass,chi,ndf,piMassB_sigma));

        KinematicParticleVertexFitter kpvFitter;
        RefCountedKinematicTree oniaTree =  kpvFitter.fit(oniaDaus);
        oniaTree->movePointerToTheTop();

        // Onia + trk + trk fit
        float chi2 = 0.;
        float ndf2 = 0.;
        vector<RefCountedKinematicParticle> ottParticles;
        ottParticles.push_back(oniaTree->currentParticle());
        ottParticles.push_back(pFactory.particle(*trk1TransTkPtr, piMassB, chi2, ndf2, piMassB_sigma));
        ottParticles.push_back(pFactory.particle(*trk2TransTkPtr, piMassB, chi2, ndf2, piMassB_sigma));

        KinematicParticleVertexFitter ottFitter;
        RefCountedKinematicTree ottVertex;
        ottVertex = ottFitter.fit(ottParticles);
        if(!ottVertex->isValid()) continue;

        ottVertex->movePointerToTheTop();

        RefCountedKinematicParticle ottCand = ottVertex->currentParticle();
        if(!ottCand->currentState().isValid()) continue;

        RefCountedKinematicVertex ottTopVertex = ottVertex->currentDecayVertex();

	      float ottC2Prob = TMath::Prob(ottTopVertex->chiSquared(),ottTopVertex->degreesOfFreedom());
	      if (ottC2Prob < bVtxChiProbCut) continue;

        // get children from final B fit
        ottVertex->movePointerToTheFirstChild();
        RefCountedKinematicParticle ottOniaCand = ottVertex->currentParticle();
        ottVertex->movePointerToTheNextChild();
        RefCountedKinematicParticle ottTrkCand1 = ottVertex->currentParticle();
        ottVertex->movePointerToTheNextChild();
        RefCountedKinematicParticle ottTrkCand2 = ottVertex->currentParticle();

        if(!ottOniaCand->currentState().isValid() || !ottTrkCand1->currentState().isValid()|| !ottTrkCand2->currentState().isValid()) continue;

        // get batchlor pion and D0 parameters from B fit
        KinematicParameters ottOniaCandKP = ottOniaCand->currentState().kinematicParameters();
        KinematicParameters ottTrkCand1KP = ottTrkCand1->currentState().kinematicParameters();
        KinematicParameters ottTrkCand2KP = ottTrkCand2->currentState().kinematicParameters();

        GlobalVector ottTotalP = GlobalVector (ottCand->currentState().globalMomentum().x(),
                                             ottCand->currentState().globalMomentum().y(),
                                             ottCand->currentState().globalMomentum().z());

        GlobalVector ottOniaTotalP = GlobalVector(ottOniaCandKP.momentum().x(),ottOniaCandKP.momentum().y(),ottOniaCandKP.momentum().z());
        GlobalVector ottTrk1TotalP = GlobalVector(ottTrkCand1KP.momentum().x(),ottTrkCand1KP.momentum().y(),ottTrkCand1KP.momentum().z());
        GlobalVector ottTrk2TotalP = GlobalVector(ottTrkCand2KP.momentum().x(),ottTrkCand2KP.momentum().y(),ottTrkCand2KP.momentum().z());

        double ottOniaTotalE = sqrt( ottOniaTotalP.mag2() + oniaMass*oniaMass  );
        double trk1TotalE = sqrt( ottTrk1TotalP.mag2() + piMassBSquared );
        double trk2TotalE = sqrt( ottTrk2TotalP.mag2() + piMassBSquared );
        double ottTotalE = ottOniaTotalE + trk1TotalE + trk2TotalE ;

        const Particle::LorentzVector ottP4(ottTotalP.x(), ottTotalP.y(), ottTotalP.z(), ottTotalE);

        Particle::Point ottVtx((*ottTopVertex).position().x(), (*ottTopVertex).position().y(), (*ottTopVertex).position().z());
        std::vector<double> bVtxEVec;
        bVtxEVec.push_back( ottTopVertex->error().cxx() );
        bVtxEVec.push_back( ottTopVertex->error().cyx() );
        bVtxEVec.push_back( ottTopVertex->error().cyy() );
        bVtxEVec.push_back( ottTopVertex->error().czx() );
        bVtxEVec.push_back( ottTopVertex->error().czy() );
        bVtxEVec.push_back( ottTopVertex->error().czz() );
        SMatrixSym3D bVtxCovMatrix(bVtxEVec.begin(), bVtxEVec.end());
        const Vertex::CovarianceMatrix bVtxCov(bVtxCovMatrix);
        double bVtxChi2(ottTopVertex->chiSquared());
        double bVtxNdof(ottTopVertex->degreesOfFreedom());
        double bNormalizedChi2 = bVtxChi2/bVtxNdof;

        double bRVtxMag = 99999.0;
        double bLVtxMag = 99999.0;
        double bSigmaRvtxMag = 999.0;
        double bSigmaLvtxMag = 999.0;

        GlobalVector bLineOfFlight = GlobalVector (ottVtx.x() - xVtx,
                                                     ottVtx.y() - yVtx,
                                                     ottVtx.z() - zVtx);

        SMatrixSym3D bTotalCov;
        if(isVtxPV) bTotalCov = bVtxCovMatrix + vtxPrimary->covariance();
        else bTotalCov = bVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

        SVector3 distanceVector3D(bLineOfFlight.x(), bLineOfFlight.y(), bLineOfFlight.z());
        SVector3 distanceVector2D(bLineOfFlight.x(), bLineOfFlight.y(), 0.0);

        double bAngle3D = angle(bLineOfFlight.x(), bLineOfFlight.y(), bLineOfFlight.z(),
                          ottTotalP.x(), ottTotalP.y(), ottTotalP.z());
        double bAngle2D = angle(bLineOfFlight.x(), bLineOfFlight.y(), (float)0.0,
                         ottTotalP.x(), ottTotalP.y(), (float)0.0);
        bLVtxMag = bLineOfFlight.mag();
        bRVtxMag = bLineOfFlight.perp();
        bSigmaLvtxMag = sqrt(ROOT::Math::Similarity(bTotalCov, distanceVector3D)) / bLVtxMag;
        bSigmaRvtxMag = sqrt(ROOT::Math::Similarity(bTotalCov, distanceVector2D)) / bRVtxMag;

        if( bNormalizedChi2 > bVtxChi2Cut ||
            bRVtxMag < bRVtxCut ||
            bRVtxMag / bSigmaRvtxMag < bRVtxSigCut ||
            bLVtxMag < bLVtxCut ||
            bLVtxMag / bSigmaLvtxMag < bLVtxSigCut ||
            cos(bAngle3D) < bCollinCut3D || cos(bAngle2D) < bCollinCut2D || bAngle3D > bAlphaCut || bAngle2D > bAlpha2DCut
        ) continue;

        // GlobalVector ottOniaTotalP = GlobalVector(ottOniaCandKP.momentum().x(),ottOniaCandKP.momentum().y(),ottOniaCandKP.momentum().z());
        // GlobalVector ottTrk1TotalP = GlobalVector(ottTrkCand1KP.momentum().x(),ottTrkCand1KP.momentum().y(),ottTrkCand1KP.momentum().z());
        // GlobalVector ottTrk2TotalP = GlobalVector(ottTrkCand2KP.momentum().x(),ottTrkCand2KP.momentum().y(),ottTrkCand2KP.momentum().z());

        RecoChargedCandidate pion1candidate(theTrackRefs[trdx]->charge(), Particle::LorentzVector(ottTrk1TotalP.x(), ottTrk1TotalP.y(), ottTrk1TotalP.z(), trk1TotalE), ottVtx);
        RecoChargedCandidate pion2candidate(theTrackRefs[trdx2]->charge(), Particle::LorentzVector(ottTrk2TotalP.x(), ottTrk2TotalP.y(), ottTrk2TotalP.z(), trk2TotalE), ottVtx);
        pion1candidate.setTrack(theTrackRefs[trdx]);
        pion2candidate.setTrack(theTrackRefs[trdx2]);

        AddFourMomenta addp4;

        VertexCompositeCandidate* theB = 0;
        theB = new VertexCompositeCandidate(0, ottP4, ottVtx, bVtxCov, bVtxChi2, bVtxNdof);
        theB->addDaughter(theOnia);
        theB->addDaughter(PionCand);
        theB->addDaughter(PionCand);

        // if(theOnia.pdgId()<0) {theB->setPdgId(521); theB->setCharge(theTrackRefs[trdx]->charge());}
        // else if(theOnia.pdgId()>0) {theB->setPdgId(-521); theB->setCharge(theTrackRefs[trdx]->charge());}

        addp4.set( *theB );

        if( theB->mass() < bMassB + bMassCut &&
            theB->mass() > bMassB - bMassCut ) theBs.push_back( *theB );
        if(theB) delete theB;
           theB = 0;
      }
  }
}
// Get methods

const reco::VertexCompositeCandidateCollection& OniapipiFitter::getB() const {
  return theBs;
}

/*
const std::vector<float>& OniapipiFitter::getMVAVals() const {
  return mvaVals_;
}
*/
/*
auto_ptr<edm::ValueMap<float> > OniapipiFitter::getMVAMap() const {
  return mvaValValueMap;
}
*/

void OniapipiFitter::resetAll() {
    theBs.clear();
//    mvaVals_.clear();
}
