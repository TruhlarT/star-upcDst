#include "MainAna.h"
#include "RunDef.h"

MainAna::MainAna(TFile *outFile): Ana(outFile){}

void MainAna::Make()
{
   runMainAna(mRecTree, NOMINAL, NOMINAL); // run standard annalysis
   
   // Vary nHits fit and dEdx cuts
   for (unsigned int i = 0; i < nTPCnHitsStudies; ++i){
      std::pair<VARIATION, VARIATION> nHitsVar = mUtil->varyNHits(i);  
      runMainAna(mHitsTree[i], nHitsVar.first, nHitsVar.second);
   }

   // Vary pid cuts
   VARIATION pidStatus[] = { LOOSE, TIGHT};
   for (unsigned int i = 0; i < nPidVariation; ++i)
      runMainAna(mPidTree[i], NOMINAL, NOMINAL, pidStatus[i]);

   for (unsigned int i = 1; i <= nDcaVariation; ++i)
      runMainAna(mDcaTree[i-1], NOMINAL, NOMINAL, NOMINAL, i);

}//Make

void MainAna::runMainAna(RecTree* recTree, VARIATION nHitFitStatus, VARIATION nHitDEdxStatus, VARIATION pid, unsigned int dca)
{
   const VARIATION dcaZStatus[] = { NOMINAL, LOOSE, TIGHT, NOMINAL, NOMINAL};
   const VARIATION dcaXYStatus[] = { NOMINAL, NOMINAL, NOMINAL, LOOSE, TIGHT};

   VARIATION dcaZ = dcaZStatus[dca]; 
   VARIATION dcaXY = dcaXYStatus[dca];



   mCurrentTree = recTree;
   bool fillPlots = nHitFitStatus == NOMINAL && nHitDEdxStatus == NOMINAL && pid == NOMINAL && dca == NOMINAL;

   if( fillPlots ) 
      hAnalysisFlow->Fill(ALL);
   if(!CheckTriggers(&CEPtriggers, mUpcEvt, hTriggerBits))
      return;

   if( fillPlots ) 
      hAnalysisFlow->Fill(TRIG);
/*
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      for (int iPmt = 0; iPmt < 2; ++iPmt)
      {
         hRpAdc[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
         hRpTac[2*iRp+iPmt]->Fill(mRpEvt->tac(iRp, iPmt));      
         if(mRpEvt->tac(iRp, iPmt) > 200 && mRpEvt->tac(iRp, iPmt) < 1750)
            hRpAdcInWindow[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
      }
*/
   AnaRpTracks(mRpEvt);

   vector<unsigned int> mRpTrackPerSide[nSides];
   unsigned int nRpTracksInFV[] = {0,0,0,0};
   unsigned int nRpTracksTotal[] = {0,0,0,0};
   for (unsigned int iSide = 0; iSide < nSides; ++iSide)
   {
      for (unsigned int iTrck = 0; iTrck < mRpTrackIdVec_perSide[iSide].size(); ++iTrck)
      {
         StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][iTrck]);
         if(!trackRP) 
            return;
         nRpTracksTotal[trackRP->branch()]++;
         if( !RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y(), trackRP->branch()) )
            continue;
         nRpTracksInFV[trackRP->branch()]++;
         mRpTrackPerSide[iSide].push_back(mRpTrackIdVec_perSide[iSide][iTrck]);
         mCurrentTree->SaveRPinfo(trackRP, iSide);
      }
   }

   for (int i = 0; i < nBranches; ++i){
      hNRPTracksInFV[i]->Fill( nRpTracksInFV[i] );
      recTree->setNRpTracksPerBranch(nRpTracksTotal[i], i);
   }

   if( !(mRpTrackPerSide[E].size()==1 && mRpTrackPerSide[W].size()==1))
      return;
   if( fillPlots ) hAnalysisFlow->Fill(TWOINFID);

   mRPpTBalance = mRpEvt->getTrack( mRpTrackPerSide[E][0] )->pVec() + mRpEvt->getTrack( mRpTrackPerSide[W][0] )->pVec();

   // Skip all events with more vertecies than 1
   hNTOFVerticies->Fill( mUpcEvt->getNumberOfVertices() );
   if( mUpcEvt->getNumberOfVertices() != 1) 
      return;
   if( fillPlots ) 
      hAnalysisFlow->Fill(ONETOFVX);

   // save event info
   mRunNumber = mUpcEvt->getRunNumber();
   mCurrentTree->SaveEventInfo(mUpcEvt);
   mCurrentTree->SaveVertexInfo(mUpcEvt->getVertex(0));

   hZVertex->Fill( mCurrentTree->getVertexZInCm() );
   if( abs(mCurrentTree->getVertexZInCm()) > vertexRange ) 
      return;
   if( fillPlots ) 
      hAnalysisFlow->Fill(ZVERTEX);

   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   int totalCharge = 0;
   unsigned int nTpcGoodTracks = 0;
   vector<int> hadronId;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // Or they are originating from different vertex than selected
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID )
         continue;

      // count number of good quality TPC tracks in the vertex
      if(IsGoodTrack(trk, nHitFitStatus, nHitDEdxStatus, dcaZ, dcaXY))
         nTpcGoodTracks++;

      if( !trk->getFlag(StUPCTrack::kTof) || !IsGoodTofTrack(trk)) 
         continue;

      hNHitFit->Fill(trk->getNhitsFit()); 
      hNHitDedx->Fill( trk->getNhitsDEdx()); 
      hDcaXY->Fill(trk->getDcaXY()); 
      hDCAZ->Fill(trk->getDcaZ()); 

      if(!IsGoodTrack(trk, nHitFitStatus, nHitDEdxStatus))
         continue; 

      hadronId.push_back(trackID);
      totalCharge += static_cast<int>( trk->getCharge() );
   } 
   mCurrentTree->setNGoodTpcTrks(nTpcGoodTracks);
   hNTOFGoodTracks->Fill( hadronId.size() );

/*
   const int nVertecies = mUpcEvt->getNumberOfVertices();

   vector<int> tracksInVertex[nVertecies]; // save the tracks Id of good TOF tracks in that vertex
   int chargeInVertex[nVertecies]; 
   for (int iVertex = 0; iVertex < nVertecies; ++iVertex)
      chargeInVertex[iVertex] = 0;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof) || !mUtil->IsGoodTrack(trk)) 
         continue;

      tracksInVertex[trk->getVertexId()].push_back(trackID);
      chargeInVertex[trk->getVertexId()] += static_cast<int>( trk->getCharge() );
   } 

   for (int iVertex = 0; iVertex < nVertecies; ++iVertex)
   {
      if( tracksInVertex[iVertex].size() != 2)
         continue;
      if( chargeInVertex[iVertex] )
         continue;
      
   }
*/

   // Skip events with other than 2 TOF tracks
   if(hadronId.size() != nSigns)
      return;

   if( fillPlots ) 
      hAnalysisFlow->Fill(TWOTOFTRKS);

   TLorentzVector hadron[nParticles], state[nParticles];
   bool inEtaRange[nSigns];
   for (unsigned int id = 0; id < nSigns; ++id)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);
      unsigned int trackOrder = id;
      if(!totalCharge)
         trackOrder = trk->getCharge() > 0 ? PLUS : MINUS;

      mCurrentTree->SaveTrackInfo(trk, trackOrder);
      for (int iPart = 0; iPart < nParticles; ++iPart){
         trk->getLorentzVector(hadron[iPart], mUtil->mass(iPart));
         state[iPart] += hadron[iPart];
      }
      inEtaRange[id] = IsGoodEtaTrack(trk);
      hEtaZVertex->Fill( trk->getEta(), mCurrentTree->getVertexZInCm());
   }

   if( !inEtaRange[PLUS] || !inEtaRange[MINUS])
      return;

   if( fillPlots ) 
      hAnalysisFlow->Fill(ETA);

   mCurrentTree->SaveTriggerInfo(mUpcEvt, mRpEvt);
   mCurrentTree->CalculatePID(true, true, -1, pid);

   hPIDStats[0]->Fill(mCurrentTree->getPairID());

   TVector3 missingMomenta = state[mCurrentTree->getPairID()].Vect() + mRPpTBalance;
   SaveMissingMomenta(missingMomenta);
   mCurrentTree->SaveRPConfigInfo();

   if( pid == NOMINAL)
      FillMSquared();
   mCurrentTree->CalculatePID(true, true, -1, pid);
   if( mCurrentTree->getPairID() == -1)
      return;

   // Fill Ana flow
   mCurrentTree->SaveStateInfo(state[mCurrentTree->getPairID()]);
   if(totalCharge) // Total charge is not zero => background event
   {
      if( fillPlots ) 
         mCurrentTree->FillBcgTree(); // Fill background Tree

      return;
   }

   if( fillPlots ) 
      hAnalysisFlow->Fill(OPPOSITE);

   hPIDStats[1]->Fill(mCurrentTree->getPairID());
   bool exclusive = false;

   if( missingMomenta.Pt() < exclusivityCut)
   {
      exclusive = true;
      if( fillPlots ) 
         hAnalysisFlow->Fill(EXCLUSIVE);
   }
   if( !keepNonExclusive && !exclusive )
      return;

   // Save their DCA to later study vertex reco eff
   //TVector3 vertex(mUpcEvt->getVertex(0)->getPosX(), mUpcEvt->getVertex(0)->getPosY(), mUpcEvt->getVertex(0)->getPosZ());
   mCurrentTree->SaveTracksDcas( mUpcEvt->getTrack(hadronId[PLUS]), mUpcEvt->getTrack(hadronId[MINUS]), mUpcEvt);

   mCurrentTree->FillRecTree(); // Fill analysis (reco) Tree

   if( mCurrentTree->getPairID() == PION && exclusive && fillPlots )
      hAnalysisFlow->Fill(PIPI);
   if( mCurrentTree->getPairID() == KAON && exclusive && fillPlots )
      hAnalysisFlow->Fill(KK);
   if( mCurrentTree->getPairID() == PROTON && exclusive && fillPlots )
      hAnalysisFlow->Fill(PPBAR);

   if( !fillPlots || !exclusive)
      return;

   for (unsigned int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackPerSide[iSide][0]);
      hRPFV[trackRP->branch()]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   }

}//runMainAna

void MainAna::Init()
{
   if( DEBUG )
      cout<<"MainAna::Init() called"<<endl;
   mOutFile->cd();

   mOutFile->mkdir("PID")->cd();
   hPIDStats[0] = new TH1D("PIDStatsTotCh", "PID stats before tot. charge cut", 4, -1.5, 2.5);
   hPIDStats[1] = new TH1D("PIDStatsExlsv", "PID stats before pT exclusive cut", 4, -1.5, 2.5);
   for (int pid = 0; pid < 4; ++pid)
      for (int stat = 0; stat < 3; ++stat)
         for (int part = 0; part < 4; ++part)
            hMSquared[pid][stat][part] = new TH1D( Form("hMSquared_%i_%i_%i",pid, stat, part), Form("hMSquared_%i_%i_%i",pid, stat, part), 200, -0.5, 1.5);
   mOutFile->cd();

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfTree[kMAINANA] , treeBits[kMAINANA], true); 
   for (unsigned int i = 0; i < nTPCnHitsStudies; ++i)
      mHitsTree[i] = new RecTree(nameOfTree[kMAINANA] +"_"+mUtil->nHitVaryName(i), treeBits[kMAINANA], false); 
   for (unsigned int i = 0; i < nPidVariation; ++i)
      mPidTree[i] = new RecTree(nameOfTree[kMAINANA] +"_"+mUtil->pidVaryName(i), treeBits[kMAINANA], false);
   for (unsigned int i = 0; i < nDcaVariation; ++i)
      mDcaTree[i] = new RecTree(nameOfTree[kMAINANA] +"_"+mUtil->dcaVaryName(i), treeBits[kMAINANA], false);

   mOutFile->mkdir("CPT2noBBCL")->cd();
   for (int iRp = 0; iRp < 2*nRomanPots; ++iRp)
   {
      hRpAdc[iRp]= new TH1D( mUtil->rpName(iRp/2) + Form("_%i_ADC",iRp%2), "ADC", 100, 0, 600);
      hRpAdcInWindow[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_ADCinTAC",iRp%2), "ADC in TAC window", 100, 0, 600);
      hRpTac[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_TAC",iRp%2), "TAC", 100, 0, 2000);
   }

   mOutFile->mkdir("AnaPlots")->cd();
   hAnalysisFlow = new TH1D("hAnaFlow", "CutsFlow", nAnalysisCuts-1, 1, nAnalysisCuts);
   for(int tb=1; tb<nAnalysisCuts; ++tb) 
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisCutName(tb));

   for (int i = 0; i < nBranches; ++i){
      hRPFV[i] = new TH2D("hRPFV"+mUtil->branchName(i),"", 220,-1.1,1.1, 220,-1.1,1.1);
      hNRPTracksInFV[i] = new TH1D("hNRPTracksInFV"+mUtil->branchName(i),"", 50,-0.5,49.5);
   }

   hEtaZVertex = new TH2D("hEtaZVertex","", 140,-1.2,1.2, 40,-100,100);

   hNTOFVerticies = new TH1D("hNTOFVerticies","", 10,-0.5,9.5);
   hZVertex = new TH1D("hZVertex","", 80,-200,200);
   hNTOFGoodTracks = new TH1D("hNTOFGoodTracks","", 15,-0.5,14.5);
   hNHitFit = new TH1D("hNHitFit","", 50, -0.5, 49.5);
   hNHitDedx = new TH1D("hNHitDedx","", 50, -0.5, 49.5);
   hDcaXY = new TH1D("hDcaXY","", 60,0,3);
   hDCAZ = new TH1D("hDCAZ","", 120,-3.0,3.0);

   mOutFile->cd();
   
}


bool MainAna::IsPairOf(int type)
{
   // 0 - PION
   // 1 - KAON
   // 2 - PROTON
   if( type < 0 || type > 2)
      return false;

   return (mCurrentTree->getNSigmaTPC(PLUS, type) < 3 && mCurrentTree->getNSigmaTPC(MINUS, type) < 3);

}

void MainAna::SaveMissingMomenta(TVector3 missP)
{
   mCurrentTree->setPtMissing( missP.Pt() );
   mCurrentTree->setPxMissing( missP.X() );
   mCurrentTree->setPyMissing( missP.Y() );  
}

void MainAna::FillMSquared()
{
   int totCharge = mCurrentTree->getCharge(PLUS) + mCurrentTree->getCharge(MINUS);
   double pTMissing = mCurrentTree->getPtMissing();
   double mSquared = mCurrentTree->getMSquared();
   
   // hMSquared[strategy][state][particle]
   // particle: All Pion Kaon Proton
   // state: before tot charge | after tot charge | exclusive
   // strategy: standard PID no pT no m2 | standard PID no m2 | selective dEdx | complete standard PID

   // fill hMSquared for "before tot charge cut"
   for (int iState = 0; iState < 3; ++iState)
   {
      // iState == 0 -> fill hMSquared for "before tot charge cut"
      // iState == 1 -> fill hMSquared for "after tot charge cut"
      // iState == 2 -> fill hMSquared for "after exclusive cut"

      if( iState == 1 && totCharge != 0 )
         continue;

      if( iState == 2 && pTMissing > exclusivityCut )
         continue;

      for (int iStrategy = 0; iStrategy < 4; ++iStrategy) // for each strategy
         hMSquared[iStrategy][iState][0]->Fill(mSquared);

      
      mCurrentTree->CalculatePID(false, false);
      int partType = mCurrentTree->getPairID() + 1;
      if( partType > 0 )
         hMSquared[0][iState][partType]->Fill(mSquared);
      
      mCurrentTree->CalculatePID(true, false);
      partType = mCurrentTree->getPairID() + 1;
      if( partType > 0 )
         hMSquared[1][iState][partType]->Fill(mSquared);

      for(int iPart = 0; iPart < nParticles; ++iPart)
         if( IsPairOf(iPart))
            hMSquared[2][iState][iPart+1]->Fill(mSquared);

      mCurrentTree->CalculatePID();
      partType = mCurrentTree->getPairID() + 1;
      if( partType > 0 )
         hMSquared[3][iState][partType]->Fill(mSquared);
   }
}
