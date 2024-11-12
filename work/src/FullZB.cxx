#include "FullZB.h"
#include "RunDef.h"

//_____________________________________________________________________________
FullZB::FullZB(TFile *outFile): Ana(outFile){}

FullZB::~FullZB(){
   //if(mUtil) delete mUtil;
}


void FullZB::Init(){
   if( DEBUG )
      cout<<"FullZB::Init() called"<<endl;
   mOutFile->cd();

   mRecTree = new RecTree(nameOfTree[kFULLZB], treeBits[kFULLZB], false); 

   mOutFile->cd();
}

void FullZB::Make(){
   if(!CheckTriggers(&ZBtriggers, mUpcEvt, nullptr))
      return;

   mRecTree->SaveTriggerInfo(mUpcEvt, mRpEvt);
   mRecTree->SaveEventInfo(mUpcEvt);

   AnaRpTracks(mRpEvt);
   for (int i = 0; i < nBranches; ++i)
      mRecTree->setNRpTracksPerBranch( mRpTrackIdVec_perBranch[i].size(), i );

   unsigned int nTofGoodTracks = 0;
   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      if( !trk->getFlag(StUPCTrack::kPrimary) || !IsGoodTrack(trk)) 
         continue;

      if( !trk->getFlag(StUPCTrack::kTof) || !IsGoodTofTrack(trk)) 
         continue;

      nTofGoodTracks++;
   } 
   
   mRecTree->setNGoodTofTracks( nTofGoodTracks );

   mRecTree->FillRecTree();
}



