#include "VertexStudy.h"
#include "RunDef.h"

VertexStudy::VertexStudy(TFile *outFile): Ana(outFile){}

void VertexStudy::Make()
{
   if( DEBUG )
      cout<<"VertexStudy::Make() called"<<endl;
   studyVertexRecoEff(); 

   if(!CheckTriggers(&CEPtriggers, mUpcEvt, nullptr))
      return;
   studyVertexZDistibution();  

}

void VertexStudy::Init()
{
   if( DEBUG )
      cout<<"VertexStudy::Init() called"<<endl;
   mOutFile->cd();

   mRecTree = new RecTree(nameOfTree[kVERTEXSTUDY], treeBits[kVERTEXSTUDY], false); 
   mVertexTree = new RecTree("vertexRecoTree", treeBits[kVERTEXSTUDY], false); 
   mVertexTree->InitVertexRecoStudy();

   hDca = new TH2D("hDca", ";dcaPart;dcaBeam", 200, 0, 40, 200, 0, 40);
   hAnaFlow = new TH1D("AnaFlow_VertexRecoStudy", "CutsFlow for VertexRecoStudy", kMax-1, 1, kMax);
   const TString CutsName[kMax] = { TString("All"), TString("2 TOF Tracks"), TString("Same Sign"),
               TString("DcaPart"), TString("DcaBeam"), TString("PrimAll"), TString("Prim 2 TOF Trks"), 
               TString("Prim Same Sign"), TString("Prim vertex") };

   for(int tb=1; tb<kMax; ++tb) 
      hAnaFlow->GetXaxis()->SetBinLabel(tb, CutsName[tb-1]);

   mOutFile->cd();
   if( DEBUG )
      cout<<"VertexStudy::Init() finished"<<endl;
}

void VertexStudy::studyVertexRecoEff()
{
   vector<int> globalTracksId; // placeholder to store indecies of the global tracks used for vertex reco eff study
   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not for vertex reco study
      if( !trk->getFlag(StUPCTrack::kCEP) )
         continue;
      if( !IsGoodGlobalTrack(trk) ) 
         continue;

      if( !trk->getFlag(StUPCTrack::kTof) || trk->getTofPathLength() < 0) 
         continue;

      globalTracksId.push_back(trackID);
   }
   hAnaFlow->Fill(kAll);
   // Skip events with other than 2 TOF global tracks
   if( globalTracksId.size() != nSigns)
      return;

   hAnaFlow->Fill(kTwoTracks);
   // Load global tracks from picoDst
   int plusIndex = mUpcEvt->getTrack(globalTracksId[PLUS])->getCharge() < 0 ? PLUS : MINUS;
   int minusIndex = plusIndex == PLUS ? MINUS : PLUS;
   const StUPCTrack* trkPlus = mUpcEvt->getTrack(globalTracksId[ plusIndex ]); // plus
   const StUPCTrack* trkMinus = mUpcEvt->getTrack(globalTracksId[ minusIndex ]); // minus

   // total charge should be zero
   if( trkPlus->getCharge() + trkMinus->getCharge() != 0)
      return;
   hAnaFlow->Fill(kSameSign);

   mVertexTree->SaveTracksDcas( trkPlus, trkMinus, mUpcEvt);
   double dcaPart, dcaBeam, vertexZhypo;
   dcaPart = mVertexTree->getVertexStudyDcaParticles();
   dcaBeam = mVertexTree->getVertexStudyDcaBeamline();
   vertexZhypo = mVertexTree->getVertexStudyHypoZ();
   hDca->Fill(dcaPart,dcaBeam);
   if( dcaPart > 3)
      return;

   hAnaFlow->Fill(kDcaPart);
   if( dcaBeam > 3)
      return;
   hAnaFlow->Fill(kDcaBeam);
   // There are 2 global TOF tracks with opposite charge, originating from the same space
   mVertexTree->SaveEventInfo(mUpcEvt);
   mVertexTree->setVertexStudyDcaParticles( dcaPart ); 
   mVertexTree->setVertexStudyDcaBeamline( dcaBeam );
   mVertexTree->setVertexStudyHypoZ( vertexZhypo );
   const StUPCTrack *trkPrimPlus, *trkPrimMinus;  
   pair<bool,bool> matchingPrimaryTracks = findPrimaryTracks(trkPrimPlus, trkPrimMinus, vertexZhypo);   

   if( matchingPrimaryTracks.first )
   {
      mVertexTree->SetDeltaTrackInfo(trkPrimPlus, trkPlus, PLUS);
      mVertexTree->SetDeltaTrackInfo(trkPrimMinus, trkMinus, MINUS);
   }else{
      mVertexTree->SetDeltaTrackInfo(PLUS);
      mVertexTree->SetDeltaTrackInfo(MINUS);
   }
   mVertexTree->setVertexStudyPrimary( matchingPrimaryTracks.first );
   mVertexTree->setVertexStudySameVertex( matchingPrimaryTracks.second );
   mVertexTree->FillRecTree();

}

pair<bool,bool> VertexStudy::findPrimaryTracks(const StUPCTrack* &trkPrimPlus, const StUPCTrack* &trkPrimMinus, double vertexZhypo)
{
   pair<bool,bool> matchingPrimaryTracks(false, false);

   vector<int> primaryId; // placeholder to store indecies of the primary TOF tracks used for vertex reco eff study
   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not for vertex reco study
      if( !trk->getFlag(StUPCTrack::kPrimary))
         continue;
      if( !trk->getFlag(StUPCTrack::kTof)  || trk->getTofPathLength() < 0)
         continue;
      if( !IsGoodGlobalTrack(trk) ) 
         continue;
      /*
      if( TMath::Abs( trk->getVertex()->getPosZ() - vertexZhypo)  > 3)
         continue;
   */
      primaryId.push_back(trackID);
   }

   hAnaFlow->Fill(kPrimAll);
   if( primaryId.size() != nSigns )
      return matchingPrimaryTracks;

   hAnaFlow->Fill(kPrimTwoTracks);

   // Load primary tracks from picoDst
   int plusIndex = mUpcEvt->getTrack(primaryId[PLUS])->getCharge() < 0 ? PLUS : MINUS;
   int minusIndex = plusIndex == PLUS ? MINUS : PLUS;
   trkPrimPlus = mUpcEvt->getTrack(primaryId[ plusIndex ]); // plus
   trkPrimMinus = mUpcEvt->getTrack(primaryId[ minusIndex ]); // minus

   if( trkPrimPlus->getCharge() + trkPrimMinus->getCharge() != 0)
      return matchingPrimaryTracks;
   hAnaFlow->Fill(kPrimSameSign);
   matchingPrimaryTracks.first = true;
   matchingPrimaryTracks.second = trkPrimPlus->getVertexId() == trkPrimMinus->getVertexId();

   if( matchingPrimaryTracks.second ){
      mRecTree->SaveVertexInfo(trkPrimPlus->getVertex());
      hAnaFlow->Fill(kPrimVertex);
   }

   return matchingPrimaryTracks;
}

void VertexStudy::studyVertexZDistibution()
{
   // Skip all events with more vertecies than 1
   if( mUpcEvt->getNumberOfVertices() != 1) return;

   // save event info
   mRecTree->SaveEventInfo(mUpcEvt);
   mRecTree->SaveVertexInfo(mUpcEvt->getVertex(0));

   mRecTree->FillRecTree(); // Fill analysis (reco) Tree
}
