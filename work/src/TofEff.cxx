#include "TofEff.h"
#include "RunDef.h"

//_____________________________________________________________________________
TofEff::TofEff(TFile *outFile): Ana(outFile){}

TofEff::~TofEff(){
   //if(mUtil) delete mUtil;
}


void TofEff::Init(){
   if( DEBUG )
      cout<<"TofEff::Init() called"<<endl;
   mOutFile->cd();

   mRecTree = new RecTree(nameOfTree[kMAINANA] + "_tofEff", treeBits[kMAINANA], false); 

   tofEffTree = new TTree("tofEffTree", "tofEffTree");
   tofEffTree->Branch("probePt", &probePt);
   tofEffTree->Branch("probeEta", &probeEta);
   tofEffTree->Branch("probePhi", &probePhi);
   tofEffTree->Branch("probeTof", &probeTof);
   tofEffTree->Branch("probeTofHit", &probeTofHit);

   tofEffTree->Branch("tagPt", &tagPt);
   tofEffTree->Branch("tagEta", &tagEta);
   tofEffTree->Branch("tagPhi", &tagPhi);

   tofEffTree->Branch("pTMiss", &pTMiss); 
   tofEffTree->Branch("invMass", &invMass);
   tofEffTree->Branch("deltaPhiProtons", &deltaPhiProtons); 
   tofEffTree->Branch("pTState", &pTState);
   tofEffTree->Branch("phiState", &phiState); 
   tofEffTree->Branch("thetaState", &thetaState); 

   mOutFile->cd();
}

void TofEff::Make(){
   if(!CheckTriggers(&CEPtriggers, mUpcEvt, nullptr))
      return;

   // Skip all events with more vertecies than 1
   if( mUpcEvt->getNumberOfVertices() != 1) 
      return;

   // save event info
   mRecTree->SaveVertexInfo(mUpcEvt->getVertex(0));

   if( abs(mRecTree->getVertexZInCm()) > vertexRange ) 
      return;


   AnaRpTracks(mRpEvt);

   vector<unsigned int> mRpTrackPerSide[nSides];
   for (unsigned int iSide = 0; iSide < nSides; ++iSide)
   {
      for (unsigned int iTrck = 0; iTrck < mRpTrackIdVec_perSide[iSide].size(); ++iTrck)
      {
         StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][iTrck]);
         if(!trackRP) 
            return;

         if( !RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y(), trackRP->branch()) )
            continue;

         mRpTrackPerSide[iSide].push_back(mRpTrackIdVec_perSide[iSide][iTrck]);
         mRecTree->SaveRPinfo(trackRP, iSide);
      }
   }

   if( !(mRpTrackPerSide[E].size()==1 && mRpTrackPerSide[W].size()==1))
      return;

   mRecTree->SaveRPConfigInfo();
   if( mRecTree->getIsElastic())
      return;

   mRPpTBalance = mRpEvt->getTrack( mRpTrackPerSide[E][0] )->pVec() + mRpEvt->getTrack( mRpTrackPerSide[W][0] )->pVec();


   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   int totalCharge = 0;
   vector<int> hadronId;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // Or they are originating from different vertex than selected
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID )
         continue;

      if( !trk->getFlag(StUPCTrack::kTof) || !IsGoodTofTrack(trk)) 
         continue;

      if(!IsGoodTrack(trk))
         continue; 

      hadronId.push_back(trackID);
      totalCharge += static_cast<int>( trk->getCharge() );
   } 

   if( hadronId.size() == 0 || hadronId.size() > 2)
      return;

   CalculateTOFEff( hadronId[0] );
   if( hadronId.size() == 2){
      CalculateTOFEff( hadronId[1] );
   }
}


void TofEff::CalculateTOFEff(unsigned int tagID)
{
   // eta cut and PID not implemented
   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   vector<int> hadronId;

   TLorentzVector tag, probe, state;
   const StUPCTrack* tagTrk = mUpcEvt->getTrack(tagID);
   int tagCharge = static_cast<int>( tagTrk->getCharge() );
   tagTrk->getLorentzVector(tag, mUtil->mass(PION));


   if( abs(tagTrk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
      return;
   if( !IsGoodEtaTrack(tagTrk) )
      return;

   double missPt, minMissPt;
   unsigned int probeID;
   minMissPt = 999;
   for(unsigned int trackID = 0; trackID < (unsigned int)mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      if( trackID == tagID)
         continue;

      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      if( (trk->getCharge() + tagCharge) != 0 )
         continue;

      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID || !IsGoodTrack(trk)) 
         continue;
      if( abs(trk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
         continue;
      if( !IsGoodEtaTrack(trk) )
         continue;
/*
      if( abs(tagTrk->getPhi() - trk->getPhi()) > 2.0)
         continue;
*/
      trk->getLorentzVector(probe, mUtil->mass(PION));
      state = tag + probe;

/*
      if( abs(state.Phi()) > 2.0 || abs(state.Phi()) < 1.3)
         continue;
      if( state.Pt() < 0.8)
         continue;

      if( state.Theta() < 0.8 || state.Theta() > 2.2)
         continue;
*/
      missPt = (state.Vect() + mRPpTBalance).Pt();
      if( missPt < minMissPt){
         minMissPt = missPt;
         probeID = trackID;
      }
   } 


   if( minMissPt > 1 )
      return;
   
   const StUPCTrack* probeTrk = mUpcEvt->getTrack(probeID);
   probeTrk->getLorentzVector(probe, mUtil->mass(PION));
   state = tag + probe;

   pTMiss = minMissPt;
   // Fill sample total
   invMass = state.M();

   tagPt = tagTrk->getPt();
   tagEta = tagTrk->getEta();
   tagPhi = tagTrk->getPhi();

   deltaPhiProtons= TMath::Abs(mRecTree->getPhiRp(East) - mRecTree->getPhiRp(West))*convertToDegree;
   if(deltaPhiProtons > 180)
      deltaPhiProtons = 360 - deltaPhiProtons;

   pTState = state.Pt();
   phiState = state.Phi();
   thetaState = state.Theta();  

   probePt = probeTrk->getPt();
   probeEta = probeTrk->getEta();
   probePhi = probeTrk->getPhi();
   probeTof = probeTrk->getFlag(StUPCTrack::kTof);
   probeTofHit = IsGoodTofTrack(probeTrk);
   tofEffTree->Fill();
}

