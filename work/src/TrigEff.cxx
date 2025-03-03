#include "TrigEff.h"
#include "RunDef.h"

//_____________________________________________________________________________
TrigEff::TrigEff(TFile *outFile): Ana(outFile){}

TrigEff::~TrigEff(){
   //if(mUtil) delete mUtil;
}


void TrigEff::Init(){
   if( DEBUG )
      cout<<"TrigEff::Init() called"<<endl;
   mOutFile->cd();

   hAnalysisFlow = new TH1D("AnaFlow_TrigEff", "CutsFlow for trigger efficiency study", nCuts-1, 1, nCuts);
   for(int tb=1; tb<nCuts; ++tb) 
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mCutName[tb-1]);

   mRecTree = new RecTree(nameOfTree[kTRIGEFF], treeBits[kTRIGEFF], false); 

   mOutFile->mkdir("Zerobias")->cd();
   for (int iRp = 0; iRp < 2*nRomanPots; ++iRp)
   {
      hRpAdc[iRp]= new TH1D( mUtil->rpName(iRp/2) + Form("_%i_ADC",iRp%2), "ADC", 100, 0, 600);
      hRpAdcInWindow[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_ADCinTAC",iRp%2), "ADC in TAC window", 100, 0, 600);
      hRpTac[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_TAC",iRp%2), "TAC", 100, 0, 2000);
   }
   hNTpPerRP = new TH2D("hNTpPerRP","Number of TPs per RP",8,-0.5,7.5,30,-0.5,29.5);
   hNClustersPerPlane = new TH2D("hNClustersPerPlane","Number of clusters per plane",32,-0.5,31.5,30,-0.5,29.5);
   for(int iRp=0; iRp<nRomanPots; ++iRp)
   {
      hNTpPerRP->GetXaxis()->SetBinLabel(iRp+1, mUtil->rpName(iRp));
      for(int iPl=0; iPl<nPlanes; ++iPl) 
         hNClustersPerPlane->GetXaxis()->SetBinLabel(iRp*nPlanes + iPl + 1, mUtil->rpName(iRp)+"_"+mUtil->planeName(iPl));
   }
   mOutFile->cd();

   mOutFile->mkdir("rpTrigEff")->cd();
   TString cutName[] ={ TString(""), TString("_PxPyCut_")};
   for (int iConf = 0; iConf < nRpConfigurations; ++iConf)
   {
      for (int iCut = 0; iCut < 2; ++iCut)
      {
         hRpTrigEffPtMiss[0][iCut][iConf] = new TH1D("hRpTrigEffPtMiss_Total_" + cutName[iCut] + mUtil->rpConfigName(iConf),"Sample total",500,0,5);
         hRpTrigEffPtMiss[1][iCut][iConf] = new TH1D("hRpTrigEffPtMiss_Passed_" + cutName[iCut] + mUtil->rpConfigName(iConf),"Sample passed",500,0,5);

         for (int iBr = 0; iBr < nBranches; ++iBr){
            hRpTrigEffPxPy[0][iCut][iBr][iConf] = new TH2D("hRpTrigEffPxPy_" + cutName[iCut]+mUtil->branchName(iBr) + "_Total_" + mUtil->rpConfigName(iConf),"Sample total",30,-1,1,30,-1,1);
            hRpTrigEffPxPy[1][iCut][iBr][iConf] = new TH2D("hRpTrigEffPxPy_" + cutName[iCut]+mUtil->branchName(iBr) + "_Passed_" + mUtil->rpConfigName(iConf),"Sample passed",30,-1,1,30,-1,1);
         }
      }    
      /*
      for (int iRp = 0; iRp < nRomanPots; ++iRp){
         hRpTrigEffXY[0][iRp][iConf] = new TH2D("hRpTrigEffXY_"+mUtil->rpName(iRp) + "_Total_" + mUtil->rpConfigName(iConf),"Sample total",30,-0.05,0.05,45,-0.075,0.075);
         hRpTrigEffXY[1][iRp][iConf] = new TH2D("hRpTrigEffXY_"+mUtil->rpName(iRp) + "_Passed_" + mUtil->rpConfigName(iConf),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
         hRpTrigEffXYOffSub[0][iRp][iConf] = new TH2D("hRpTrigEffXYOffSub_"+mUtil->rpName(iRp) + "_Total_" + mUtil->rpConfigName(iConf),"Sample total",30,-0.05,0.05,45,-0.075,0.075);
         hRpTrigEffXYOffSub[1][iRp][iConf] = new TH2D("hRpTrigEffXYOffSub_"+mUtil->rpName(iRp) + "_Passed_" + mUtil->rpConfigName(iConf),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
      }*/
   }
   mOutFile->cd();
}

void TrigEff::Make(){
   hAnalysisFlow->Fill(ALL);
   if(!CheckTriggers(&ZBtriggers, mUpcEvt, nullptr))
      return;
   hAnalysisFlow->Fill(TRIG);

   if( mRpEvt == nullptr){
      cerr<<"Error RP event not set"<<endl;
      return;
   }

   // save event info
   mRecTree->SaveEventInfo(mUpcEvt);
   mRecTree->SaveTriggerInfo(mUpcEvt, mRpEvt);
   
   // Analyze RP tracks and fill RP info
   AnaRpTracks(mRpEvt);
   FillRPInfo();
   CheckRPTrigEff();

   // Skip all events with more vertecies than 1
   if( mUpcEvt->getNumberOfVertices() != 1) return;
   hAnalysisFlow->Fill(ONETOFVX);
   mRecTree->SaveVertexInfo(mUpcEvt->getVertex(0));

   vector<int> hadronId;
   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof) || !IsGoodTrack(trk) || !IsGoodTofTrack(trk)) 
         continue;

      hadronId.push_back(trackID);
   } 

   // Skip events with less than 2 TOF tracks
   if(hadronId.size() != nSigns)
      return;

   hAnalysisFlow->Fill(TWOTOFTRKS);
   // save tracks info
   for (unsigned int id = 0; id < hadronId.size(); ++id)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);
      mRecTree->SaveTrackInfo(trk, id);
   }
   mRecTree->FillRecTree();
}

void TrigEff::FillRPInfo()
{
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      for (int iPmt = 0; iPmt < 2; ++iPmt)
      {
         hRpAdc[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
         hRpTac[2*iRp+iPmt]->Fill(mRpEvt->tac(iRp, iPmt));      
         if(mRpEvt->tac(iRp, iPmt) > 200 && mRpEvt->tac(iRp, iPmt) < 1750)
            hRpAdcInWindow[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
      }

      hNTpPerRP->Fill(iRp,mTrackPointIdVec[iRp].size());
      for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
         hNClustersPerPlane->Fill(iRp*nPlanes + iPlane, mClusterIdVec[iRp*nPlanes + iPlane].size() );
   }
}

void TrigEff::CheckRPTrigEff()
{
   if( !(mRpTrackIdVec_perSide[E].size()==1 && mRpTrackIdVec_perSide[W].size()==1))
      return;

   bool protonInRange[nSides];
   int rpOri[nSides];
   int branch[nSides];
   double pX[nSides], pY[nSides];
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]);
      if(!trackRP) return;
      branch[iSide] = trackRP->branch();
      pX[iSide] = trackRP->pVec().X(); 
      pY[iSide] = trackRP->pVec().Y(); 
      protonInRange[iSide] = RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y(), branch[iSide]);
      rpOri[iSide] = trackRP->pVec().Y() > 0 ? Up : Down; 
   }


   // There are two good RP track = sample total
   int rpConfig = mUtil->mapBranchConfiguration(rpOri[E], rpOri[W]) / nRpConfigurations;
   double missPt = (mRpEvt->getTrack( mRpTrackIdVec_perSide[E][0] )->pVec() + mRpEvt->getTrack( mRpTrackIdVec_perSide[W][0] )->pVec()).Pt();
   hRpTrigEffPtMiss[0][0][rpConfig]->Fill(missPt);
   if(protonInRange[0] && protonInRange[1]) 
      hRpTrigEffPtMiss[0][1][rpConfig]->Fill(missPt);

   for (int iSide = 0; iSide < nSides; ++iSide){
      hRpTrigEffPxPy[0][0][branch[iSide]][rpConfig]->Fill(pX[iSide],pY[iSide]);
      if(protonInRange[0] && protonInRange[1]) 
         hRpTrigEffPxPy[0][1][branch[iSide]][rpConfig]->Fill(pX[iSide],pY[iSide]);
   }

   // Check if there is corresponding RP_IT or RP_ET trigger
   if( (rpConfig == IT && !mRecTree->getRpItDsmBit()) || (rpConfig == ET && !mRecTree->getRpEtDsmBit()))
      return;

   hRpTrigEffPtMiss[1][0][rpConfig]->Fill(missPt);
   if(protonInRange[0] && protonInRange[1]) 
      hRpTrigEffPtMiss[1][1][rpConfig]->Fill(missPt);

   for (int iSide = 0; iSide < nSides; ++iSide){
      hRpTrigEffPxPy[1][0][branch[iSide]][rpConfig]->Fill(pX[iSide],pY[iSide]);
      if(protonInRange[0] && protonInRange[1]) 
         hRpTrigEffPxPy[1][1][branch[iSide]][rpConfig]->Fill(pX[iSide],pY[iSide]);
   }
}


const TString TrigEff::mCutName[nCuts] = { TString("All"), TString("Zerobias"), TString("1 TOF vtx"), 
      TString("2 good TOF trks")};
