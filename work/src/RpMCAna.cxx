#include "UpcDstLibreries.h"
#include "RunDef.h"
#include "RpMCAna.h"

RpMCAna::RpMCAna(TFile *outFile): Ana(outFile){
   embedMaker = new EmbedMaker();
   embedMaker->Init(outFile);
   embedMaker->setAfterburner(AFTERBURNER);
}

RpMCAna::~RpMCAna(){
   FillAverageEffInFV();
   for (int iSet = 0; iSet < DATA; ++iSet)
      if(mElAna[iSet]) delete mElAna[iSet];

   if(embedMaker) delete embedMaker;
}

void RpMCAna::Make()
{
   // check RP trigger
   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      runEmbedding(iSet);
      mElAna[iSet]->SetRpPosition(mCorrection, mOffSet); 
      mElAna[iSet]->SetEvent(mUpcEvt, mEmbedEvt, mMcEvt);
      mElAna[iSet]->SetRPMCInfo(mVertex,mMomentum);
      mElAna[iSet]->Make();
      runMCEff(iSet);
      if(mEmbedEvt) delete mEmbedEvt;
   }
}

void RpMCAna::Init()
{
   if( DEBUG )
      cout<<"RpMCAna::Init() called"<<endl;

   mRecTree = new RecTree(nameOfTree[kRPMCANA],treeBits[kRPMCANA], false); 

   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      mElAna[iSet] = new ElasticAna(mOutFile);
      mElAna[iSet]->SetAnaName(mUtil->dataSetName(iSet));
      mElAna[iSet]->SetTriggers(&noTriggers);  
      mElAna[iSet]->Init();
      mElAna[iSet]->InitRPMCInfo();

      // RP efficiency
      mOutFile->mkdir("RPMC_"+mUtil->dataSetName(iSet))->cd();
      hMCEffFlow[iSet] = new TH1D("AnaFlow_" + anaName, "CutsFlow", kMax-1, 1, kMax);
      for(int tb=1; tb<kMax; ++tb) 
         hMCEffFlow[iSet]->GetXaxis()->SetBinLabel(tb, mMCEffFlowCutsName[tb-1]);

      hMCtoRecoDiff[iSet] = new TH1D("hMCtoRecoDiff_"+mUtil->dataSetName(iSet),"hMCtoRecoDiff",1000,0,1);
      for (int iBr = 0; iBr < nBranches; ++iBr){
         hEffPxPy[0][iBr][iSet] = new TH2D("hEffPxPy_MC_"+mUtil->branchName(iBr) + "_Total","Sample total",55,-0.55,0.55,120,-1.2,1.2);
         hEffPxPy[1][iBr][iSet] = new TH2D("hEffPxPy_MC_"+mUtil->branchName(iBr) + "_PassedBefore","Sample passed",55,-0.55,0.55,120,-1.2,1.2);
         hEffPxPy[2][iBr][iSet] = new TH2D("hEffPxPy_MC_"+mUtil->branchName(iBr) + "_Passed","Sample passed",55,-0.55,0.55,120,-1.2,1.2);         
         hEffRpCheckedPxPy[0][iBr][iSet] = new TH2D("hEffRpCheckedPxPy_MC_"+mUtil->branchName(iBr) + "_Total","Sample total",55,-0.55,0.55,120,-1.2,1.2);
         hEffRpCheckedPxPy[1][iBr][iSet] = new TH2D("hEffRpCheckedPxPy_MC_"+mUtil->branchName(iBr) + "_PassedBefore","Sample passed",55,-0.55,0.55,120,-1.2,1.2);
         hEffRpCheckedPxPy[2][iBr][iSet] = new TH2D("hEffRpCheckedPxPy_MC_"+mUtil->branchName(iBr) + "_Passed","Sample passed",55,-0.55,0.55,120,-1.2,1.2);              
      }
      for (int iRp = 0; iRp < nRomanPots; ++iRp){ 
         hEffXY[0][iRp][iSet] = new TH2D("hEffXY_MC_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
         hEffXY[1][iRp][iSet] = new TH2D("hEffXY_MC_"+mUtil->rpName(iRp) + "_PassedBefore","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffXY[2][iRp][iSet] = new TH2D("hEffXY_MC_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffXYOffSub[0][iRp][iSet] = new TH2D("hEffXYOffSub_MC_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
         hEffXYOffSub[1][iRp][iSet] = new TH2D("hEffXYOffSub_MC_"+mUtil->rpName(iRp) + "_PassedBefore","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffXYOffSub[2][iRp][iSet] = new TH2D("hEffXYOffSub_MC_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXY[0][iRp][iSet] = new TH2D("hEffRpCheckedXY_MC_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXY[1][iRp][iSet] = new TH2D("hEffRpCheckedXY_MC_"+mUtil->rpName(iRp) + "_PassedBefore","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXY[2][iRp][iSet] = new TH2D("hEffRpCheckedXY_MC_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXYOffSub[0][iRp][iSet] = new TH2D("hEffRpCheckedXYOffSub_MC_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXYOffSub[1][iRp][iSet] = new TH2D("hEffRpCheckedXYOffSub_MC_"+mUtil->rpName(iRp) + "_PassedBefore","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
         hEffRpCheckedXYOffSub[2][iRp][iSet] = new TH2D("hEffRpCheckedXYOffSub_MC_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
      }
      hNTpPerRP[iSet] = new TH2D("hNTpPerRP_"+mUtil->dataSetName(iSet),"Number of TPs per RP",8,-0.5,7.5,30,-0.5,29.5);
      hNClustersPerPlane[iSet] = new TH2D("hNClustersPerPlane_"+mUtil->dataSetName(iSet),"Number of clusters per plane",32,-0.5,31.5,30,-0.5,29.5);
      for(int iRp=0; iRp<nRomanPots; ++iRp)
      {
         hNTpPerRP[iSet]->GetXaxis()->SetBinLabel(iRp+1, mUtil->rpName(iRp));
         for(int iPl=0; iPl<nPlanes; ++iPl) 
            hNClustersPerPlane[iSet]->GetXaxis()->SetBinLabel(iRp*nPlanes + iPl + 1, mUtil->rpName(iRp)+"_"+mUtil->planeName(iPl));
      }
      mOutFile->cd();

      mOutFile->mkdir("RPMC_"+mUtil->dataSetName(iSet)+"_EffStudy")->cd();
      mRunNumber = 999;
      for (int iBr = 0; iBr < nBranches; ++iBr){
         hAverageEff[iBr][iSet] = new TH1D("hAverageEff_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet) + Form("_%i",mRunNumber),";<Eff>;",400,0.0,2.0);
         hAverageEffRpChecked[iBr][iSet] = new TH1D("hAverageEffRpChecked_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet) + Form("_%i",mRunNumber),";<Eff>;",400,0.0,2.0);
         hAvrgEffperBr[iBr][iSet] = new TH1D("hAverageEffPerBr_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet),";<Eff>;",400,0.0,2.0);
      }
      mOutFile->cd();
   }
}


//void RpMCAna::SetRpPosition(TVector3 (&corr)[nRomanPots], TVector3 (&offsets)[nRomanPots])
void RpMCAna::SetRpPosition(TVector3* corr, TVector3 *offsets)
{
   mCorrection = corr;
   mOffSet = offsets;
   embedMaker->setOffsetCorrections(corr);
}

void RpMCAna::SetMCInfo(double (&mc_vtx)[nCoordinates], double (&mc_p)[nCoordinates][nSides])
{
   mVertex = mc_vtx;
   mMomentum = mc_p;
}

bool RpMCAna::IsTrackInRp(int side)
{
   if( abs(mVertex[Z]) > vertexRange/100 )
      return false;

   return true;

   int branch = (mMomentum[Z][side] < 0) ? EU : WU; 
   if( mMomentum[Y][side] < 0)
      branch++;

   for (int iRp = 0; iRp < nStationPerSide; ++iRp)
   {
      int RP = mUtil->rpPerBranchStationOrder(branch, iRp); 
      double x, y;
      ProjectToRP(x,y,RP,side);

      if(!IsInRpRange(x,y,RP,mOffSet[RP]))
         return false;
   }

  return true;
}

void RpMCAna::ProjectToRP(double& x, double& y, int RP, int side)
{
   x = mVertex[X]/1000 + (mUtil->rpZPosition(RP) - mVertex[Z])*(mMomentum[X][side]/mMomentum[Z][side]);
   y = mVertex[Y]/1000 + (mUtil->rpZPosition(RP) - mVertex[Z])*(mMomentum[Y][side]/mMomentum[Z][side]);
}


bool RpMCAna::AreTracksInRp()
{
  return IsTrackInRp(East) && IsTrackInRp(West);
}

void RpMCAna::runEmbedding(bool embedding)
{
   embedMaker->setMCEvent(mMcEvt);
   embedMaker->setZBEvent(mRpEvt);
   mEmbedEvt = new StRPEvent(*mRpEvt);
   mEmbedEvt->clearEvent();
   embedMaker->setRPEvent(mEmbedEvt);
   embedMaker->setEmbedding(embedding);
   embedMaker->setVertex( mVertex[X], mVertex[Y], mVertex[Z]);
   embedMaker->MakeTracks(mUtil->p0(),mUtil->p0());
}

void RpMCAna::runMCEff(int set)
{
   int branch;
   AnaRpTracks(mEmbedEvt);
   // plot efficiency as function of momenta
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      hMCEffFlow[set]->Fill(kAll);
      if(!IsTrackInRp(iSide))
         continue;

      hMCEffFlow[set]->Fill(kInside);
      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;

      hMCEffFlow[set]->Fill(kTrigger);
      // fill sample total 
      hEffPxPy[0][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);

      if( mRpTrackIdVec_perBranch[branch].size()!=1 )
         continue;
      hMCEffFlow[set]->Fill(kOneTrack);      
      StUPCRpsTrack* trk = mEmbedEvt->getTrack(mRpTrackIdVec_perBranch[branch][0]);
      if(!trk) continue;

      double MCtoRecoDiff = sqrt( pow(trk->pVec().X() - mMomentum[X][iSide], 2) + pow(trk->pVec().Y() - mMomentum[Y][iSide], 2));
      // plot difference between MC and reconstructed hit
      if(RPInFidRange(mMomentum[X][iSide], mMomentum[Y][iSide], branch))
         hMCtoRecoDiff[set]->Fill( MCtoRecoDiff);

      hEffPxPy[1][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
      if(MCtoRecoDiff > maxMcToRecoDist)
         continue;

      hMCEffFlow[set]->Fill(kDist);
      hEffPxPy[2][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
   }
   // plot efficiency as function of momenta, checked signal in RP_paired
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      if(!IsTrackInRp(iSide))
         continue;

      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;

      for (int iSt = 0; iSt < nStationPerSide; ++iSt)
      {
         if( mTrackPointIdVec[ mUtil->rpPerBranchStationOrder(branch, iSt) ].size() !=1 )
            continue;

         // fill sample total 
         hEffRpCheckedPxPy[0][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);

         if( mRpTrackIdVec_perBranch[branch].size()!=1 )
            continue;

         StUPCRpsTrack* trk = mEmbedEvt->getTrack(mRpTrackIdVec_perBranch[branch][0]);
         if(!trk) continue;

         double MCtoRecoDiff = sqrt( pow(trk->pVec().X() - mMomentum[X][iSide], 2) + pow(trk->pVec().Y() - mMomentum[Y][iSide], 2));

         hEffRpCheckedPxPy[1][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
         if(MCtoRecoDiff > maxMcToRecoDist)
            continue;

         hEffRpCheckedPxPy[2][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
      }
   }

   // plot efficiency as function of RP coordinates (x,y)
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;
      // fill sample total 
      double xMC[nStationPerSide], yMC[nStationPerSide];
      for (int iRp = 0; iRp < nStationPerSide; ++iRp)
      {
         int RP = mUtil->rpPerBranchStationOrder(branch, iRp); 
         ProjectToRP(xMC[iRp],yMC[iRp],RP,iSide);
         hEffXY[0][RP][set]->Fill(xMC[iRp]*100,yMC[iRp]*100);
         hEffXYOffSub[0][RP][set]->Fill((xMC[iRp])*100,(yMC[iRp]-mOffSet[RP][Y])*100);
         if( IsInRpRange( xMC[iRp], yMC[iRp], RP, mOffSet[RP]))
         { // plot # hits in plane (x vs y) and #TPs in RP in general
            hNTpPerRP[set]->Fill(RP,mTrackPointIdVec[RP].size());
            for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
               hNClustersPerPlane[set]->Fill(RP*nPlanes + iPlane, mClusterIdVec[RP*nPlanes + iPlane].size() );
         }

         if( mTrackPointIdVec[RP].size() !=1 )
            continue;

         StUPCRpsTrackPoint *trkPoint = mEmbedEvt->getTrackPoint(mTrackPointIdVec[RP][0]);
         if(!trkPoint) continue;
         double x = trkPoint->x();
         double y = trkPoint->y();
         hEffXY[1][RP][set]->Fill(xMC[iRp]*100,yMC[iRp]*100);
         hEffXYOffSub[1][RP][set]->Fill((xMC[iRp])*100,(yMC[iRp]-mOffSet[RP][Y])*100);

         double MCtoRecoDiff = sqrt( pow(x - xMC[iRp], 2) + pow(y - yMC[iRp], 2));
         if(MCtoRecoDiff > maxDiffBetweenRpHitMCAndReco)
            continue;
         hEffXY[2][RP][set]->Fill(xMC[iRp]*100,yMC[iRp]*100);
         hEffXYOffSub[2][RP][set]->Fill((xMC[iRp])*100,(yMC[iRp]-mOffSet[RP][Y])*100);      }      
   }
   // plot efficiency as function of RP coordinates (x,y), checked signal in RP_paired
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;
      // fill sample total 
      double xMC[nStationPerSide], yMC[nStationPerSide];
      for (int iSt = 0; iSt < nStationPerSide; ++iSt)
      {
         int rpPaired = mUtil->rpPerBranchStationOrder(branch, (iSt+1)%nStationPerSide );
         if( mTrackPointIdVec[rpPaired].size() !=1 )
            continue;

         int RP = mUtil->rpPerBranchStationOrder(branch, iSt); 
         ProjectToRP(xMC[iSt],yMC[iSt],RP,iSide);
         hEffRpCheckedXY[0][RP][set]->Fill(xMC[iSt]*100,yMC[iSt]*100);
         hEffRpCheckedXYOffSub[0][RP][set]->Fill((xMC[iSt])*100,(yMC[iSt]-mOffSet[RP][Y])*100);

         if( mTrackPointIdVec[RP].size() !=1 )
            continue;

         StUPCRpsTrackPoint *trkPoint = mEmbedEvt->getTrackPoint(mTrackPointIdVec[RP][0]);
         if(!trkPoint) continue;
         double x = trkPoint->x();
         double y = trkPoint->y();
         hEffRpCheckedXY[1][RP][set]->Fill(xMC[iSt]*100,yMC[iSt]*100);
         hEffRpCheckedXYOffSub[1][RP][set]->Fill((xMC[iSt])*100,(yMC[iSt]-mOffSet[RP][Y])*100);

         double MCtoRecoDiff = sqrt( pow(x - xMC[iSt], 2) + pow(y - yMC[iSt], 2));
         if(MCtoRecoDiff > maxDiffBetweenRpHitMCAndReco)
            continue;
         hEffRpCheckedXY[2][RP][set]->Fill(xMC[iSt]*100,yMC[iSt]*100);
         hEffRpCheckedXYOffSub[2][RP][set]->Fill((xMC[iSt])*100,(yMC[iSt]-mOffSet[RP][Y])*100);      }      
   }
}

void RpMCAna::FillAverageEffInFV()
{
   auto fillAverageEff = [&](TH2D* total, TH2D* passed, TH1D* hAvrgEff, TH1D* hEffperBr, int br) 
   {
      auto hEff = (TH2D*)passed->Clone();
      hEff->Divide(total);
      double nPassed, nTotal;
      nPassed = 0;
      nTotal = 0;
      for (int i = 2; i <= hEff->GetNbinsX(); ++i) { 
         double px = hEff->GetXaxis()->GetBinCenter(i);
         for (int j = 2; j <= hEff->GetNbinsY(); ++j) {
            double py = hEff->GetYaxis()->GetBinCenter(j);
            if( !RPInFidRange( px, py, br) )
               continue;
            //if( hEff->GetBinContent(i, j) > 0)
            hAvrgEff->Fill( hEff->GetBinContent(i, j) );
            nPassed += passed->GetBinContent(i, j);
            nTotal += total->GetBinContent(i, j);
         }
      }
      double rartio = nTotal == 0 ? 0.0 : nPassed/ nTotal;
      if(rartio > 1.0)
         cout<<"Warrning: in FillAverageEffInFV the ff is bigger than 1. The eff = "<<rartio<<" for br "<<mUtil->branchName(br)<<" for run "<<mRunNumber<<endl;

      hEffperBr->Fill( rartio );
   };

   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      for (int iBr = 0; iBr < nBranches; ++iBr)
      { 
         hAverageEff[iBr][iSet]->SetName("hAverageEff_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet) + Form("_%i",mRunNumber));
         hAverageEffRpChecked[iBr][iSet]->SetName("hAverageEffRpChecked_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet) + Form("_%i",mRunNumber));

         fillAverageEff(hEffPxPy[0][iBr][iSet], hEffPxPy[2][iBr][iSet], hAverageEff[iBr][iSet], hAvrgEffperBr[iBr][iSet], iBr);
         //fillAverageEff(hEffRpCheckedPxPy[0][iBr][iSet], hEffRpCheckedPxPy[2][iBr][iSet], hAverageEffRpChecked[iBr][iSet], iBr);
      }
   }

}//FillAverageEffInFV


bool RpMCAna::Veto(const StUPCEvent *upcEvt, const StRPEvent *rpEvt) const
{
   mRecTree->SaveEventInfo(upcEvt);
   mRecTree->SaveTriggerInfo(upcEvt, rpEvt);

   if( mRecTree->getBbcDsmBit() ) return false; // trigger signal in BBC
   if( mRecTree->getZdcDsmBit() ) return false; // trigger signal in ZDC
   if( mRecTree->getNVertecies() > 0 ) return false; // more TPC-TOF vertecies than 0
   if( mRecTree->getTofMult() > 8 ) return false; 

   bool signalInBrunch[nBranches] = { false, false, false, false};
   for (int iBr = 0; iBr < nBranches; ++iBr)
      for (int iSt = 0; iSt < nStationPerSide; ++iSt)
         if(mRecTree->getRpTrigBits(mUtil->rpPerBranchStationOrder(iBr, iSt)))
            signalInBrunch[iBr] = true;

   for (int iSide = 0; iSide < nSides; ++iSide)
      if( signalInBrunch[nSides*iSide] && signalInBrunch[nSides*iSide + 1]) // the same as !(brunchUp != brunchDown)
         return false; 

   return true;
}//Veto


const TString RpMCAna::mMCEffFlowCutsName[kMax] = { TString("All"), TString("Inside RP"), TString("Trigger"),
               TString("1 track"), TString("Distance")};
