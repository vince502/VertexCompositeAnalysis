
#include "PrimaryVertexAnalyzerMultiFS.hxx"


void PrimaryVertexAnalyzerMultiFS::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup){
    using std::vector;
    using namespace edm;
    using namespace reco;

    using std::cout;
    using std::endl;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_, vertices);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_, v0candidates);
    const reco::VertexCompositeCandidateCollection *v0candidates_ = v0candidates.product();

    VertexCompositeCandidate v01;
    if (v0candidates_->size() > 0){
        v01 = (*v0candidates_)[0];
        cout << "V0 candidate : " <<  v01.daughter(0)->pt() << ", " <<  v01.daughter(1)->pt() << endl;
    } 
    const reco::VertexCollection vtxCollection = *(vertices.product());
    reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
    if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2){
    //   xVtx = vtxPrimary->x();
    //   yVtx = vtxPrimary->y();
    //   zVtx = vtxPrimary->z();
    //   xVtxError = vtxPrimary->xError();
    //   yVtxError = vtxPrimary->yError();
    //   zVtxError = vtxPrimary->zError();
        cout << "Analyzing First (best) vertex candidate with assoc Ntracks " << vtxPrimary->tracksSize() << "in pT(weight)" << endl;
        for( auto vtxTrk = vtxPrimary->tracks_begin(); vtxTrk != vtxPrimary->tracks_end(); vtxTrk ++){
            auto& tk = *vtxTrk;
            cout << tk->pt() << "(" << vtxPrimary->trackWeight(*vtxTrk) << ") " ;
        }
        cout << endl;
    }

};