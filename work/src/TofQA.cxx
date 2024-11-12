#include "TofQA.h"
#include "RunDef.h"

TofQA::TofQA(TFile *outFile): Ana(outFile){}

void TofQA::Make()
{
   if(!CheckTriggers(&CEPtriggers, mUpcEvt, nullptr))
      return;

   unsigned int newRunNumber = mUpcEvt->getRunNumber();
   if(newRunNumber != mRunNumber){
      mRunNumber = newRunNumber;
      hEtaPhi->SetName(Form("hEtaPhi_%i",mRunNumber));
   }

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof) || !IsGoodTrack(trk) || !IsGoodTofTrack(trk)) 
         continue;

      hEtaPhi->Fill(trk->getEta(), trk->getPhi());
   } 

}

void TofQA::Init()
{
   if( DEBUG )
      cout<<"TofQA::Init() called"<<endl;

   mRunNumber = 999;
   mOutFile->mkdir("TofQA")->cd();
   hEtaPhi = new TH2D(Form("hEtaPhi_%i",mRunNumber),"",100,-2,2,100,-4,4);
   mOutFile->cd();
}
