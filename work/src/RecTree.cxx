#include "RecTree.h"

//_____________________________________________________________________________
RecTree::RecTree(TString treeName, bitset<8> treeVersion, bool isBcgTree) {
   // treeVersion bits: 11111111 = all information is stored
   // 0000 0001 - event info
   // 0000 0010 - Vertex info
   // ...

   mUtil = new Util();
   //standard reconstructed tree
   mRecTree = new TTree(treeName, treeName);

   // event info
   mRecTree->Branch("runNumber", &mRunNumber);
   if( treeVersion.test(0) )
   {
      mRecTree->Branch("nVertecies", &mNVertecies);
      mRecTree->Branch("nGoodTpcTrks", &mNGoodTpcTrks);
      mRecTree->Branch("tofMult", &mTofMult);
      mRecTree->Branch("nTofGoodTracks", &mNTofGoodTracks);
      for (int iBr = 0; iBr < nBranches; ++iBr)
         mRecTree->Branch("nRpTracksPer" + mUtil->branchName(iBr), &mNRpTracksPerBranch[iBr]);

      mRecTree->Branch("eventNumber", &mEventNumber);
      mRecTree->Branch("fillNumber", &mFillNumber);
      mRecTree->Branch("bunchCrossId", &mBunchCrossId);
      mRecTree->Branch("bunchCrossId7bit", &mBunchCrossId7bit);
   }

   // Vertex info
   mVertexIsSet = false;
   if( treeVersion.test(1) )
   {
      for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
      {
         TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
         mRecTree->Branch("vertexZInCm" + tag, &mVertexZInCm[iTag]);
         mRecTree->Branch("vertexYInCm" + tag, &mVertexYInCm[iTag]);
         mRecTree->Branch("vertexXInCm" + tag, &mVertexXInCm[iTag]);
         if( iTag == TRUEMC )
            mRecTree->Branch("vertexIdTruth" + tag, &mVertexIdTruth);
      }
   }


   // State info
   if( treeVersion.test(2) )
   {
      mRecTree->Branch("invMass", &mInvMass);
      mRecTree->Branch("theta", &mTheta);
      mRecTree->Branch("phi", &mPhi);
      mRecTree->Branch("p", &mP);
      mRecTree->Branch("pt", &mPt);
      mRecTree->Branch("pairRapidity", &mRap);
      mRecTree->Branch("mSquared", &mMSquared);
      mRecTree->Branch("cEq",&cEq);
      mRecTree->Branch("bEq",&bEq);
      mRecTree->Branch("aEq",&aEq);
      mRecTree->Branch("pairID", &mPairID);
   }

   // Central hadrons info
   if( treeVersion.test(3) )
   {
      for (int i = 0; i < nSigns; ++i)
      { 
         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            mRecTree->Branch(Form("momentumInGev%i",i) + tag, &mMomentumInGev[i][iTag]);
            mRecTree->Branch(Form("pTInGev%i",i) + tag, &mPtInGev[i][iTag]);
            mRecTree->Branch(Form("pXInGev%i",i) + tag, &mPxInGev[i][iTag]);
            mRecTree->Branch(Form("pYInGev%i",i) + tag, &mPyInGev[i][iTag]);
            mRecTree->Branch(Form("pZInGev%i",i) + tag, &mPzInGev[i][iTag]);
            mRecTree->Branch(Form("charge%i",i) + tag, &mCharge[i][iTag]);
            mRecTree->Branch(Form("eta%i",i) + tag, &mEtaHadrons[i][iTag]);
            mRecTree->Branch(Form("phi%i",i) + tag, &mPhiHadrons[i][iTag]);
            if( iTag == TRUEMC )
               mRecTree->Branch(Form("QA%i",i) + tag, &mQATruth[i]);
         } 

         mRecTree->Branch(Form("dEdxInKevCm%i",i), &mDEdxInKevCm[i]);
         mRecTree->Branch(Form("tofTimeInNs%i",i), &mTofTimeInNs[i]);
         mRecTree->Branch(Form("tofLengthInCm%i",i), &mTofLengthInCm[i]);
         mRecTree->Branch(Form("dcaXYInCm%i",i), &mDcaXYInCm[i]);
         mRecTree->Branch(Form("dcaZInCm%i",i), &mDcaZInCm[i]);
         mRecTree->Branch(Form("nHitsFit%i",i), &mNHitsFit[i]);
         mRecTree->Branch(Form("nHitsDEdx%i",i), &mNHitsDEdx[i]);
         for (int iPart = 0; iPart < nParticles; ++iPart)
            mRecTree->Branch("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(i), &mNSigmaTPC[i][iPart]);
      }      
   }

   // RP track info 
   if( treeVersion.test(4) )
   {
      mRecTree->Branch("isElastic", &mIsElastic);
      mRecTree->Branch("RPDeltaPhi", &mRPDeltaPhi);
      for (int i = 0; i < nSides; ++i)
      {
         mRecTree->Branch("rpX" + mUtil->sideName(i), &mRpX[i]);
         mRecTree->Branch("rpY" + mUtil->sideName(i), &mRpY[i]);
         mRecTree->Branch("rpZ" + mUtil->sideName(i), &mRpZ[i]);
         mRecTree->Branch("thetaRp" + mUtil->sideName(i), &mThetaRp[i]);
         mRecTree->Branch("phiRp" + mUtil->sideName(i), &mPhiRp[i]);
         mRecTree->Branch("timeRp" + mUtil->sideName(i), &mTimeRp[i]);
         mRecTree->Branch("pRp" + mUtil->sideName(i), &mPRp[i]);
         mRecTree->Branch("ptRp" + mUtil->sideName(i), &mPtRp[i]);
         mRecTree->Branch("etaRp" + mUtil->sideName(i), &mEtaRp[i]);
         mRecTree->Branch("pXRp" + mUtil->sideName(i), &mPx[i]);
         mRecTree->Branch("pYRp" + mUtil->sideName(i), &mPy[i]);
         mRecTree->Branch("pZRp" + mUtil->sideName(i), &mPz[i]);
         mRecTree->Branch("t" + mUtil->sideName(i), &mT[i]);
         mRecTree->Branch("xi" + mUtil->sideName(i), &mXi[i]);
      }
   }

   // pT Missing info (if cental hadrons and RP track info)
   if( treeVersion.test(3) && treeVersion.test(4))
   {
      mRecTree->Branch("pTMissing", &mPtMissing);
      mRecTree->Branch("pXMissing", &mPxMissing);
      mRecTree->Branch("pYMissing", &mPyMissing);
   }

   // ZDC and BBC info
   if( treeVersion.test(5) )
   {
      mRecTree->Branch("BbcSmallEast", &mBbcSmallEast);
      mRecTree->Branch("BbcSmallWest", &mBbcSmallWest);
      mRecTree->Branch("BbcLargeEast", &mBbcLargeEast);
      mRecTree->Branch("BbcLargeWest", &mBbcLargeWest);
      for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
      {
         mRecTree->Branch(Form("ZdcAdcEastPmt%i",iPmt), &mZdcAdcEastPmt[iPmt]);
         mRecTree->Branch(Form("ZdcAdcWestPmt%i",iPmt), &mZdcAdcWestPmt[iPmt]);
      }
      mRecTree->Branch("ZdcTdcEast", &mZdcTdcEast);
      mRecTree->Branch("ZdcTdcWest", &mZdcTdcWest);
      mRecTree->Branch("ZdcTimeDiff", &mZdcTimeDiff);
      mRecTree->Branch("ZdcVertexZ", &mZdcVertexZ);
      mRecTree->Branch("ZdcEastRate", &mZdcEastRate); 
      mRecTree->Branch("ZdcWestRate", &mZdcWestRate); 
      mRecTree->Branch("ZdcEastUA", &mZdcEastUA); 
      mRecTree->Branch("ZdcWestUA", &mZdcWestUA); 
   }

   // trigger efficiency info
   if( treeVersion.test(6) )
   {
      mRecTree->Branch("BbcDsmBit", &mBbcDsmBit);  
      mRecTree->Branch("BbcSmallEDsmBit", &mBbcSmallEDsmBit);
      mRecTree->Branch("BbcSmallWDsmBit", &mBbcSmallWDsmBit);
      mRecTree->Branch("BbcLargeEDsmBit", &mBbcLargeEDsmBit);
      mRecTree->Branch("BbcLargeWDsmBit", &mBbcLargeWDsmBit);
      mRecTree->Branch("ZdcDsmBit", &mZdcDsmBit);
      mRecTree->Branch("ZdcEDsmBit", &mZdcEDsmBit);
      mRecTree->Branch("ZdcWDsmBit", &mZdcWDsmBit);

      mRecTree->Branch("BbcTrigBit", &mBbcTrigBit);  
      mRecTree->Branch("ZdcTrigBit", &mZdcTrigBit);
      mRecTree->Branch("ZdcETrigBit", &mZdcETrigBit);
      mRecTree->Branch("ZdcWTrigBit", &mZdcWTrigBit);  

      mRecTree->Branch("tofTrigBit", &mTofTrigBit);
      mRecTree->Branch("tofDsmBit", &mTofDsmBit);
      mRecTree->Branch("tofDsmABit", &mTofDsmABit);
      mRecTree->Branch("tofDsmBBit", &mTofDsmBBit);

      mRecTree->Branch("RpTrigBit", &mRpTrigBit);  
      mRecTree->Branch("rpEtTrigBit", &mRpEtTrigBit);  
      mRecTree->Branch("rpItTrigBit", &mRpItTrigBit);  
      mRecTree->Branch("RpDsmBit", &mRpDsmBit);  
      mRecTree->Branch("rpEtDsmBit", &mRpEtDsmBit);
      mRecTree->Branch("rpItDsmBit", &mRpItDsmBit); 
      for (int iRp = 0; iRp < nRomanPots; ++iRp)
      {
         mRecTree->Branch("rp" + mUtil->rpName(iRp) + "DsmBit", &mRpTrigBits[iRp]); 
      }
   }

   // Save vertex reco info
   if( treeVersion.test(7) )
   {
      mRecTree->Branch("dcaParticles", &mDcaParticles);
      mRecTree->Branch("dcaBeamline", &mDcaBeamline);
      mRecTree->Branch("vertexZHypo", &mVertexZHypo);
   }

   // Setting background Tree
   if(isBcgTree){
      mBcgTree = mRecTree->CloneTree(0);
      mBcgTree->SetName(treeName + "_Bcg");
   }
}//RecTree::CreateRecTree


//_____________________________________________________________________________
RecTree::RecTree(TTree* tree, bitset<8> treeVersion) {

   mUtil = new Util();
   // event info
   tree->SetBranchAddress("runNumber", &mRunNumber);
   if( treeVersion.test(0) )
   {
      tree->SetBranchAddress("nVertecies", &mNVertecies);
      tree->SetBranchAddress("nGoodTpcTrks", &mNGoodTpcTrks);
      tree->SetBranchAddress("tofMult", &mTofMult);
      tree->SetBranchAddress("nTofGoodTracks", &mNTofGoodTracks);
      for (int iBr = 0; iBr < nBranches; ++iBr)
         tree->SetBranchAddress("nRpTracksPer" + mUtil->branchName(iBr), &mNRpTracksPerBranch[iBr]);
      
      tree->SetBranchAddress("eventNumber", &mEventNumber);
      tree->SetBranchAddress("fillNumber", &mFillNumber);
      tree->SetBranchAddress("bunchCrossId", &mBunchCrossId);
      tree->SetBranchAddress("bunchCrossId7bit", &mBunchCrossId7bit);
   }

   // Vertex info
   if( treeVersion.test(1) )
   {
      for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
      {
         TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
         tree->SetBranchAddress("vertexZInCm" + tag, &mVertexZInCm[iTag]);
         tree->SetBranchAddress("vertexYInCm" + tag, &mVertexYInCm[iTag]);
         tree->SetBranchAddress("vertexXInCm" + tag, &mVertexXInCm[iTag]);
         if( iTag == TRUEMC )
            tree->SetBranchAddress("vertexIdTruth" + tag, &mVertexIdTruth);
      }
   }


   // State info
   if( treeVersion.test(2) )
   {
      tree->SetBranchAddress("invMass", &mInvMass);
      tree->SetBranchAddress("theta", &mTheta);
      tree->SetBranchAddress("phi", &mPhi);
      tree->SetBranchAddress("p", &mP);
      tree->SetBranchAddress("pt", &mPt);
      tree->SetBranchAddress("pairRapidity", &mRap);
      tree->SetBranchAddress("mSquared", &mMSquared);
      tree->SetBranchAddress("cEq",&cEq);
      tree->SetBranchAddress("bEq",&bEq);
      tree->SetBranchAddress("aEq",&aEq);
      tree->SetBranchAddress("pairID", &mPairID);
   }

   // Central hadrons info
   if( treeVersion.test(3) )
   {
      for (int i = 0; i < nSigns; ++i)
      { 
         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            tree->SetBranchAddress(Form("momentumInGev%i",i) + tag, &mMomentumInGev[i][iTag]);
            tree->SetBranchAddress(Form("pTInGev%i",i) + tag, &mPtInGev[i][iTag]);
            tree->SetBranchAddress(Form("pXInGev%i",i) + tag, &mPxInGev[i][iTag]);
            tree->SetBranchAddress(Form("pYInGev%i",i) + tag, &mPyInGev[i][iTag]);
            tree->SetBranchAddress(Form("pZInGev%i",i) + tag, &mPzInGev[i][iTag]);
            tree->SetBranchAddress(Form("charge%i",i) + tag, &mCharge[i][iTag]);
            tree->SetBranchAddress(Form("eta%i",i) + tag, &mEtaHadrons[i][iTag]);
            tree->SetBranchAddress(Form("phi%i",i) + tag, &mPhiHadrons[i][iTag]);
            if( iTag == TRUEMC )
               tree->SetBranchAddress(Form("QA%i",i) + tag, &mQATruth[i]);
         } 

         tree->SetBranchAddress(Form("dEdxInKevCm%i",i), &mDEdxInKevCm[i]);
         tree->SetBranchAddress(Form("tofTimeInNs%i",i), &mTofTimeInNs[i]);
         tree->SetBranchAddress(Form("tofLengthInCm%i",i), &mTofLengthInCm[i]);
         tree->SetBranchAddress(Form("dcaXYInCm%i",i), &mDcaXYInCm[i]);
         tree->SetBranchAddress(Form("dcaZInCm%i",i), &mDcaZInCm[i]);
         tree->SetBranchAddress(Form("nHitsFit%i",i), &mNHitsFit[i]);
         tree->SetBranchAddress(Form("nHitsDEdx%i",i), &mNHitsDEdx[i]);
         for (int iPart = 0; iPart < nParticles; ++iPart)
            tree->SetBranchAddress("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(i), &mNSigmaTPC[i][iPart]);
      }      
   }

   // RP track info 
   if( treeVersion.test(4) )
   {
      tree->SetBranchAddress("isElastic", &mIsElastic); 
      tree->SetBranchAddress("RPDeltaPhi", &mRPDeltaPhi);
      for (int i = 0; i < nSides; ++i)
      {
         tree->SetBranchAddress("rpX" + mUtil->sideName(i), &mRpX[i]);
         tree->SetBranchAddress("rpY" + mUtil->sideName(i), &mRpY[i]);
         tree->SetBranchAddress("rpZ" + mUtil->sideName(i), &mRpZ[i]);
         tree->SetBranchAddress("thetaRp" + mUtil->sideName(i), &mThetaRp[i]);
         tree->SetBranchAddress("phiRp" + mUtil->sideName(i), &mPhiRp[i]);
         tree->SetBranchAddress("timeRp" + mUtil->sideName(i), &mTimeRp[i]);
         tree->SetBranchAddress("pRp" + mUtil->sideName(i), &mPRp[i]);
         tree->SetBranchAddress("ptRp" + mUtil->sideName(i), &mPtRp[i]);
         tree->SetBranchAddress("etaRp" + mUtil->sideName(i), &mEtaRp[i]);
         tree->SetBranchAddress("pXRp" + mUtil->sideName(i), &mPx[i]);
         tree->SetBranchAddress("pYRp" + mUtil->sideName(i), &mPy[i]);
         tree->SetBranchAddress("pZRp" + mUtil->sideName(i), &mPz[i]);
         tree->SetBranchAddress("t" + mUtil->sideName(i), &mT[i]);
         tree->SetBranchAddress("xi" + mUtil->sideName(i), &mXi[i]);
      }
   }

   // pT Missing info (if cental hadrons and RP track info)
   if( treeVersion.test(3) && treeVersion.test(4))
   {
      tree->SetBranchAddress("pTMissing", &mPtMissing);
      tree->SetBranchAddress("pXMissing", &mPxMissing);
      tree->SetBranchAddress("pYMissing", &mPyMissing);
   }

   // ZDC and BBC info
   if( treeVersion.test(5) )
   {
      tree->SetBranchAddress("BbcSmallEast", &mBbcSmallEast);
      tree->SetBranchAddress("BbcSmallWest", &mBbcSmallWest);
      tree->SetBranchAddress("BbcLargeEast", &mBbcLargeEast);
      tree->SetBranchAddress("BbcLargeWest", &mBbcLargeWest);
      for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
      {
         tree->SetBranchAddress(Form("ZdcAdcEastPmt%i",iPmt), &mZdcAdcEastPmt[iPmt]);
         tree->SetBranchAddress(Form("ZdcAdcWestPmt%i",iPmt), &mZdcAdcWestPmt[iPmt]);
      }
      tree->SetBranchAddress("ZdcTdcEast", &mZdcTdcEast);
      tree->SetBranchAddress("ZdcTdcWest", &mZdcTdcWest);
      tree->SetBranchAddress("ZdcTimeDiff", &mZdcTimeDiff);
      tree->SetBranchAddress("ZdcVertexZ", &mZdcVertexZ);
      tree->SetBranchAddress("ZdcEastRate", &mZdcEastRate); 
      tree->SetBranchAddress("ZdcWestRate", &mZdcWestRate); 
      tree->SetBranchAddress("ZdcEastUA", &mZdcEastUA); 
      tree->SetBranchAddress("ZdcWestUA", &mZdcWestUA); 
   }

   // trigger efficiency info
   if( treeVersion.test(6) )
   {
      tree->SetBranchAddress("BbcDsmBit", &mBbcDsmBit);  
      tree->SetBranchAddress("BbcSmallEDsmBit", &mBbcSmallEDsmBit);
      tree->SetBranchAddress("BbcSmallWDsmBit", &mBbcSmallWDsmBit);
      tree->SetBranchAddress("BbcLargeEDsmBit", &mBbcLargeEDsmBit);
      tree->SetBranchAddress("BbcLargeWDsmBit", &mBbcLargeWDsmBit);
      tree->SetBranchAddress("ZdcDsmBit", &mZdcDsmBit);
      tree->SetBranchAddress("ZdcEDsmBit", &mZdcEDsmBit);
      tree->SetBranchAddress("ZdcWDsmBit", &mZdcWDsmBit);

      tree->SetBranchAddress("BbcTrigBit", &mBbcTrigBit);  
      tree->SetBranchAddress("ZdcTrigBit", &mZdcTrigBit);
      tree->SetBranchAddress("ZdcETrigBit", &mZdcETrigBit);
      tree->SetBranchAddress("ZdcWTrigBit", &mZdcWTrigBit);  

      tree->SetBranchAddress("tofTrigBit", &mTofTrigBit);
      tree->SetBranchAddress("tofDsmBit", &mTofDsmBit);
      tree->SetBranchAddress("tofDsmABit", &mTofDsmABit);
      tree->SetBranchAddress("tofDsmBBit", &mTofDsmBBit);

      tree->SetBranchAddress("RpTrigBit", &mRpTrigBit);  
      tree->SetBranchAddress("rpEtTrigBit", &mRpEtTrigBit);  
      tree->SetBranchAddress("rpItTrigBit", &mRpItTrigBit);  
      tree->SetBranchAddress("RpDsmBit", &mRpDsmBit);  
      tree->SetBranchAddress("rpEtDsmBit", &mRpEtDsmBit);
      tree->SetBranchAddress("rpItDsmBit", &mRpItDsmBit); 
      for (int iRp = 0; iRp < nRomanPots; ++iRp)
      {
         tree->SetBranchAddress("rp" + mUtil->rpName(iRp) + "DsmBit", &mRpTrigBits[iRp]); 
      }
   }

   // Save vertex reco info
   if( treeVersion.test(7) )
   {
      tree->SetBranchAddress("dcaParticles", &mDcaParticles);
      tree->SetBranchAddress("dcaBeamline", &mDcaBeamline);
      tree->SetBranchAddress("vertexZHypo", &mVertexZHypo);
   }
   mRecTree = tree;
}//RecTree::LoadTTree

RecTree::~RecTree(){
   if(mUtil) delete mUtil;
   if(mRecTree) delete mRecTree;
}



void RecTree::InitRPMCInfo()
{
   mRecTree->Branch("rpTrueVertexX", &mTrueVertexX);
   mRecTree->Branch("rpTrueVertexY", &mTrueVertexY);
   mRecTree->Branch("rpTrueVertexZ", &mTrueVertexZ);

   for (int i = 0; i < nSides; ++i)
   {
      mRecTree->Branch("rpTruePx" + mUtil->sideName(i), &mTruePx[i]);
      mRecTree->Branch("rpTruePy" + mUtil->sideName(i), &mTruePy[i]);
      mRecTree->Branch("rpTruePz" + mUtil->sideName(i), &mTruePz[i]);
   }
}

void RecTree::InitRPMCInfo(TTree* tree)
{
   tree->SetBranchAddress("rpTrueVertexX", &mTrueVertexX);
   tree->SetBranchAddress("rpTrueVertexY", &mTrueVertexY);
   tree->SetBranchAddress("rpTrueVertexZ", &mTrueVertexZ);

   for (int i = 0; i < nSides; ++i)
   {
      tree->SetBranchAddress("rpTruePx" + mUtil->sideName(i), &mTruePx[i]);
      tree->SetBranchAddress("rpTruePy" + mUtil->sideName(i), &mTruePy[i]);
      tree->SetBranchAddress("rpTruePz" + mUtil->sideName(i), &mTruePz[i]);
   }
}

void RecTree::InitVertexRecoStudy()
{

   //cout<<"RecTree::InitVertexRecoStudy() called"<<endl;   
   mRecTree->Branch("primary", &mPrimary);
   mRecTree->Branch("sameVertex", &mSameVertex);

   for (int i = 0; i < nSigns; ++i)
   { 
      mRecTree->Branch(Form("deltaPt%i",i), &mDeltaPt[i]);
      mRecTree->Branch(Form("deltaEta%i",i), &mDeltaEta[i]);
      mRecTree->Branch(Form("deltaPhi%i",i), &mDeltaPhi[i]);
      mRecTree->Branch(Form("deltaCurvature%i",i), &mDeltaCurvature[i]);
      mRecTree->Branch(Form("deltaDipAngle%i",i), &mDeltaDipAngle[i]);
      mRecTree->Branch(Form("deltaPhase%i",i), &mDeltaPhase[i]); 
   }    

}

void RecTree::InitVertexRecoStudy(TTree* tree)
{

   tree->SetBranchAddress("primary", &mPrimary);
   tree->SetBranchAddress("sameVertex", &mSameVertex);

   for (int i = 0; i < nSigns; ++i)
   { 
      tree->SetBranchAddress(Form("deltaPt%i",i), &mDeltaPt[i]);
      tree->SetBranchAddress(Form("deltaEta%i",i), &mDeltaEta[i]);
      tree->SetBranchAddress(Form("deltaPhi%i",i), &mDeltaPhi[i]);
      tree->SetBranchAddress(Form("deltaCurvature%i",i), &mDeltaCurvature[i]);
      tree->SetBranchAddress(Form("deltaDipAngle%i",i), &mDeltaDipAngle[i]);
      tree->SetBranchAddress(Form("deltaPhase%i",i), &mDeltaPhase[i]); 
   } 
}

void RecTree::CalculatePID(bool pTSpace, bool m2, bool smearing)
{

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = this->getNSigmaTPC(PLUS, iPart)*this->getNSigmaTPC(PLUS, iPart) + 
                     this->getNSigmaTPC(MINUS, iPart)*this->getNSigmaTPC(MINUS, iPart);

   TRandom3 gen(0);
   // smearing: https://pdf.sciencedirectassets.com/271580/1-s2.0-S0168900206X06836/1-s2.0-S0168900205024411/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEN%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIF%2Bk2nayCi0a%2BhnGW3aeik0DOu7dYBhf4p2pTe7EC%2F2cAiAyom6s2fbSxaeZS1pbQSrauTV9lyJ7VgCUwlIwdMXB2Sq8BQjY%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMEUcqst4GIoIr50JjKpAFesanjys0FoPa4RlzLXQePkRdKE12Dwvdr6MPoj4Lwou7rERi6rn167Y34XkiCuoMiwqUnZ9wrTSn2ykasj%2BSpIkMn9N4NtyPcyZ4ZHAZ%2FHdSv0BUUP96F1ru9iS4IfT3ii%2B1N9mdwMvTAyGGAmhDk1PTFQBOY4hMGDIpGwuEHFeew6A2mLsnEAa0f8OdnoWO0KHAHJkeXsLQDjRsNgBrWPfedAHqqXdQ%2BnKHKgqNxdXkDETFGaijMho4kOrnSmye2YhsPJXrNBErcrchnEEQmhf403f89bmThRgIVf88RhOMxSS51OW0Edft6kyY2sqyCuE7ARkexR4qlYO1yrQIeSIHdyrk1Vj9pXEpLlpnrJH%2FcrqZvBDBGsUCoc6TEmB4LN09hLGDOuHDzl0j7KGENjobhw6FQKVwhe8fxvZ1kzIBUAaPpht36%2BEmRw1ymUKTxihSSdLd5%2FgZ7PBRTYv7CuZTs2q5vmd%2B9a7IUy9wijB1s4lWAM5d9nfcCoplHV4rJj9DRTH8IncT3SND15mONnFk4%2FkQbtmWRZ65h0FO0O3myLkESQVN%2BDKY%2BIuP9keZoo%2BxrBLtnokZv90J9eknDblb9tlc0osDdaz%2BSZUjaqHxk5O8Ig6u9LXDMo8jzyuV7reAe2NIgT0kb4cxNKMfRUL84tukhnxBQ9%2BRND23Q%2Fhoxnc8qEuMmHpNX7RGwGTK9Ws3DjtX%2BXk6PGzrbEo%2FcICzSedfL3vA1aCFkQepgZ9y%2BYiWeXWDIJ002FUTYFySYntK9RpZV%2Fry19NUUKS2%2FVuUV6nAYodvu1egFzGrdxK9Vl1UwZBCxQG66XGotuLiwY%2BpdhfLwad590Y2GgoIpav5ZnHInvNcB3nEszpK%2B2IwzKf4tQY6sgFbj5NoS1TvJQVVJ0CXdmZmaRQQt5%2BNSlJVywKy5O7a8f4z8eQTqtHYQ2AfmwkE9lI2iF3dWwcnFUyYdYj3vg3wlD09UHbbddQHo5BbsD3aEiCrNXmFCb0EvHhBvtB4NRf4mMwQ%2FTW4uLEtNNSRnNWItQDhwOL%2Bloj15eu67SYZ6aumK0Zo%2B%2BqytqPVd%2BJSfglRYggAV08zGpPBUNMfOgP5kv4sALPQaEDGNDQ9BRCTJp7I&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240815T150733Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYRQNPFS6N%2F20240815%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=79452b3bff35fb8132a52be2405c1ace5a1fb16542e8ee63ea5d8de7e8f44afe&hash=db8555accce1201a8f1b728e0d8bfc773dc3cde54b789515cc18a565004e8c8d&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0168900205024411&tid=spdf-911c6aea-15e8-49e6-98a0-afecd0542f9c&sid=ed83bb835c152243ba39aa9488c886a4dacdgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=191f590d5d025e585701&rr=8b3a184e6d0e4114&cc=cz
   double deltaTOF = this->getTofTimeInNs(MINUS) - this->getTofTimeInNs(PLUS);
   if(smearing)
      deltaTOF+= gen.Gaus(0.0, 0.11)+gen.Gaus(0.0, 0.11);
   double speedOfLight2 = mUtil->c()*mUtil->c();
   double speedOfLight4 = speedOfLight2*speedOfLight2; 
   double length1Squared = smearing ? pow((this->getTofLengthInCm(PLUS)+gen.Gaus(0.0, 0.5) )/100,2) : this->getTofLengthInCm(PLUS)*this->getTofLengthInCm(PLUS)/(100*100); // convert TOFlength from cm to m
   double length2Squared = smearing ? pow((this->getTofLengthInCm(MINUS)+gen.Gaus(0.0, 0.5) )/100,2) : this->getTofLengthInCm(MINUS)*this->getTofLengthInCm(MINUS)/(100*100); // convert TOFlength from cm to m

   double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
   double deltaTime4 = deltaTime2*deltaTime2;

   double momentum1 = smearing ? this->getMomentumInGeV(PLUS,TRUEMC) + gen.Gaus(0.0, this->getPtInGev(PLUS,TRUEMC) < 0.3 ? 0.006 : 0.0024 + 0.012*this->getPtInGev(PLUS,TRUEMC)) : this->getMomentumInGeV(PLUS);
   double momentum2 = smearing ? this->getMomentumInGeV(MINUS,TRUEMC) + gen.Gaus(0.0, this->getPtInGev(MINUS,TRUEMC) < 0.3 ? 0.006 : 0.0024 + 0.012*this->getPtInGev(MINUS,TRUEMC)) : this->getMomentumInGeV(MINUS);

   double oneOverMomentum1sq = 1.0/pow(momentum1,2);
   double oneOverMomentum2sq = 1.0/pow(momentum2,2);
   cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
   bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
   aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
   double mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
   double pT[] = { this->getPtInGev(PLUS), this->getPtInGev(MINUS)};

   if(chiPair[PION] > 9 && chiPair[KAON] > 9 && chiPair[PROTON] < 9 && ( m2 ? mSquared > m2minProtons : true)
      && (pTSpace ? (pT[PLUS] > minPt[PROTON] && pT[MINUS] > minPt[PROTON] && min(pT[PLUS], pT[MINUS]) < minPtPair[PROTON]) : true))
   { // it is... proton!
         pairID = PROTON;
   }
   else if(chiPair[PION] > 9 && chiPair[KAON] < 9 && chiPair[PROTON] > 9 && ( m2 ? mSquared > m2minKaons : true)
      && (pTSpace ? (pT[PLUS] > minPt[KAON] && pT[MINUS] > minPt[KAON] && min(pT[PLUS], pT[MINUS]) < minPtPair[KAON]) : true))
   { // it is... kaon!
         pairID = KAON;
   } 
   else if( chiPair[PION] < 12 && (pTSpace ? (pT[PLUS] > minPt[PION] && pT[MINUS] > minPt[PION]) : true) ) // it is... pion!
   {
      pairID = PION;
   }
   
   this->setPairID(pairID);
   this->setMSquared(mSquared);
}

unsigned int RecTree::getBranch(unsigned int side)
{
   int branch = (mPz[side] < 0) ? EU : WU; 
   if( mPy[side] < 0)
      branch++;
   return branch;
}//getBranch