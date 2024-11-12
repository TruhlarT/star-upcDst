#include "UpcDstLibreries.h"
#include "RunDef.h"
#include "ElasticAna.h"

//_____________________________________________________________________________
ElasticAna::ElasticAna(TFile *outFile): Ana(outFile){
   SetAnaName(studyName[kELASTICANA]);
   SetTriggers(&ElasticTriggers);
}

ElasticAna::~ElasticAna(){
   //if(mUtil) delete mUtil;
}

void ElasticAna::Init(){
   if( DEBUG )
      cout<<"ElasticAna::Init() called"<<endl;
   mOutFile->cd();
   hAnalysisFlow = new TH1D("SimpleAnaFlow_" + anaName, "CutsFlow for simple elastic ana", nCuts-1, 1, nCuts);
   for(int tb=1; tb<nCuts; ++tb) 
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mCutName[tb-1]);

   hElasticAnaFlow = new TH1D("AnaFlow_" + anaName, "CutsFlow for elastic ana", kMax-1, 1, kMax);
   for(int tb=1; tb<kMax; ++tb) 
      hElasticAnaFlow->GetXaxis()->SetBinLabel(tb, mElAnaCutsName[tb-1]);

   mRecTree = new RecTree(anaName + "Tree",treeBits[kELASTICANA], false); 
   mOutFile->mkdir(anaName)->cd();

   hTheta[E][X] = new TH1D("hThetaEX", "hThetaEX", 3000, -0.03, 0.03);
   hTheta[E][Y] = new TH1D("hThetaEY", "hThetaEY", 3000, -0.03, 0.03);
   hTheta[W][X] = new TH1D("hThetaWX", "hThetaWX", 3000, -0.03, 0.03);
   hTheta[W][Y] = new TH1D("hThetaWY", "hThetaWY", 3000, -0.03, 0.03);

   hDeltaTheta[0] = new TH1D("hDeltaTheta", "hDeltaTheta", 3000, -0.03, 0.03);
   hDeltaTheta[1] = new TH1D("hDeltaTheta_ET", "hDeltaTheta_ET", 3000, -0.03, 0.03); 
   hThetaCorr = new TH2D("hThetaCorr", "hThetaCorr", 200, -0.002, 0.002, 200, -0.002, 0.002);// bin size 0.000 040 rad

   hT = new TH1D("hT", "hT", 50, 0.2, 1.2);

   hDCutR = new TH1D("hDCutR", "hDCutR", 200, 0, 0.1);
   hDCut = new TH2D("hDCut", "hDCut", 200, -0.005, 0.005, 200, -0.005, 0.005);

   for (int iRp = 0; iRp < 2*nRomanPots; ++iRp)
   {
      hRpAdc[iRp]= new TH1D( mUtil->rpName(iRp/2) + Form("_%i_ADC",iRp%2), "ADC", 100, 0, 600);
      hRpAdcInWindow[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_ADCinTAC",iRp%2), "ADC in TAC window", 100, 0, 600);
      hRpTac[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_TAC",iRp%2), "TAC", 100, 0, 2000);
   }

   hFindHitDiff = new TH1D("hFindHitDiff", "hFindHitDiff", 100, .0, 0.1);
   
   for (int iBr = 0; iBr < nBranches; ++iBr){
      hEffPxPy[0][iBr] = new TH2D("hEffPxPy_"+mUtil->branchName(iBr) + "_Total","Sample total",55,-0.55,0.55,120,-1.2,1.2);
      hEffPxPy[1][iBr] = new TH2D("hEffPxPy_"+mUtil->branchName(iBr) + "_Passed","Sample passed",55,-0.55,0.55,120,-1.2,1.2);
   }
   for (int iRp = 0; iRp < nRomanPots; ++iRp){ 
      hEffXY[0][iRp] = new TH2D("hEffXY_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
      hEffXY[1][iRp] = new TH2D("hEffXY_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
      hEffXYOffSub[0][iRp] = new TH2D("hEffXYOffSub_"+mUtil->rpName(iRp) + "_Total","Sample total",27,-4.5,4.5,51,-8.5,8.5);
      hEffXYOffSub[1][iRp] = new TH2D("hEffXYOffSub_"+mUtil->rpName(iRp) + "_Passed","Sample passed",27,-4.5,4.5,51,-8.5,8.5);
   }

   if( anaName == "ElasticAna")
   {
      mRunNumber = 999;
      for (int iArm = 0; iArm < nArms; ++iArm)
      {
         for (int iCoor = 0; iCoor < nCoordinates - 1; ++iCoor)
            hVertexExtraction[iArm][iCoor] = new TH2D("hVertexExtraction_" + mUtil->armName(iArm) + "_" + mUtil->coordinateName(iCoor), "hVertexExtraction_" + mUtil->armName(iArm) + "_" + mUtil->coordinateName(iCoor), 200, -5, 5, 200, -80, 80);
         hVertexZExtraction[iArm] = new TH1D("hVertexExtraction_" + mUtil->armName(iArm) + "_Z", "hVertexExtraction_" + mUtil->armName(iArm) + "_Z", 100, -300, 300);
      }
   }

   for (int iRP = 0; iRP < nRomanPots; ++iRP)
   {
      hClusterLength[iRP] = new TH1D("hClusterLength_"+mUtil->rpName(iRP), "hClusterLength_"+mUtil->rpName(iRP), 80, -0.5, 79.5);
      hClusterEnergy[iRP] = new TH1D("hClusterEnergy_"+mUtil->rpName(iRP), "hClusterEnergy_"+mUtil->rpName(iRP), 250, -0.5, 249.5);
      hNClusterPerRP[iRP] = new TH1D("hNClusterPerRP_"+mUtil->rpName(iRP),"hNClusterPerRP_"+mUtil->rpName(iRP), 100, -0.5, 99.5);
      for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
         hClusterPerPlain[iRP][iPlane] = new TH1D("hClusterPerPlain_"+mUtil->planeName(iPlane)+"_"+mUtil->rpName(iRP), "hClusterPerPlain_"+mUtil->rpName(iRP)+"_"+mUtil->planeName(iPlane), 80, -0.5, 79.5);
   }

   for (int iBr = 0; iBr < nBranches; ++iBr)
      for (int iCoor = 0; iCoor < nCoordinates - 1; ++iCoor)
         hDeltaProjection[iCoor][iBr] = new TH1D("hDeltaProj_" + mUtil->coordinateName(iCoor) + "_" + mUtil->branchName(iBr), "hDeltaProj_" + mUtil->coordinateName(iCoor) + "_" + mUtil->branchName(iBr), 24, -0.6, 0.6);

   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      for (int iCoor = 0; iCoor < nCoordinates - 1; ++iCoor)
         hDeltaFitProjection[iRp][iCoor] = new TH1D("hDeltaFitProj_" + mUtil->coordinateName(iCoor) + "_" + mUtil->rpName(iRp), "hDeltaFitProj_" + mUtil->rpName(iRp) + "_" + mUtil->coordinateName(iCoor), 60, -0.8, 0.8);

   for (int iBr = 0; iBr < nBranches; ++iBr)
      hAcceptancePxPy[iBr] = new TH2D("hAcceptancePxPy_"+mUtil->branchName(iBr),"Px, Py acceptance",100,-1,1,100,-1,1);
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      hAcceptanceXY[iRp] = new TH2D("hAcceptanceXY_"+mUtil->rpName(iRp),"x, y acceptance",30,-0.05,0.05,45,-0.075,0.075);
   mOutFile->cd();

   if( runALIGNMENT ){
      mOutFile->mkdir("Alignment")->cd();
      for (unsigned int i = 0; i < nAligIteration; ++i)
      {
         for (unsigned int iRp = 0; iRp < nRomanPots; ++iRp)
         {
            hAligEvents[iRp][i] = new TH1D("hAligEvents_" + mUtil->rpName(iRp) + Form("_%i",i), "", 100, 0, 1000);
            hAligXCorr[iRp][i] = new TH1D("hAligXCorr_" + mUtil->rpName(iRp) + Form("_%i",i), "", 100, -0.2, 0.2);
            hAligYCorr[iRp][i] = new TH1D("hAligYCorr_" + mUtil->rpName(iRp) + Form("_%i",i), "", 100, -0.2, 0.2);         
         }
      }
      mOutFile->cd();      
   }
}

void ElasticAna::Make()
{
   hElasticAnaFlow->Fill(ALL);

   if(!CheckTriggers(trigger, mUpcEvt, nullptr))
      return;

   for (int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      for (int iPmt = 0; iPmt < 2; ++iPmt)
      {
         hRpAdc[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
         hRpTac[2*iRp+iPmt]->Fill(mRpEvt->tac(iRp, iPmt));      
         if(mRpEvt->tac(iRp, iPmt) > 200 && mRpEvt->tac(iRp, iPmt) < 1750)
            hRpAdcInWindow[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
      }
   }
   hElasticAnaFlow->Fill(TRIG);
   AnaRpTracks(mRpEvt);
   runRpDataDrivenEffStudy();
   FillAccaptancePlots();

   bool isElastic = IsElasticEventCandidate();

   FillElasticPlots();
   if( !isElastic )
      return;

   hElasticAnaFlow->Fill(kElastic);
   FillClusterInfo();
   runMCValidation();

   if( runALIGNMENT )
      runAlignment();

   if( trigger->size() )    
      runVertexExtraction();

   // save event info
   mRecTree->SaveEventInfo(mUpcEvt);
   for (int iSide = 0; iSide < nSides; ++iSide)
      mRecTree->SaveRPinfo(mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]), iSide);

   mRecTree->SaveTriggerInfo(mUpcEvt, mRpEvt);
   mRecTree->FillRecTree();
}

bool  ElasticAna::IsElasticEventCandidate(int RPStudy)
{
   bitset<8> rpsBits, rpsBitsEast, rpsBitsWest; // [0,0,0,0,0,0,0,0]
   // in cout -> [W2D, W2U, W1D, W1U, E2D, E2U, E1D, E1U] 

   static const std::bitset<8> EUWD( std::string("10100101") );
   static const std::bitset<8> EDWU( std::string("01011010") );

   // init trackEW
   trackEW.t = -99.0;
   trackEW.nPoints = 0;

   if( RPStudy != -1)
      rpsBits.set(RPStudy);

   unsigned int nPoints = 0;
   for (int iRP = 0; iRP < nRomanPots; ++iRP){
      nPoints += iRP == RPStudy ? 1 : mTrackPointIdVec[iRP].size();
      if(mTrackPointIdVec[iRP].size() >= 1)
         rpsBits.set(iRP);
   }

   if(rpsBits != EUWD && rpsBits != EDWU)
      return false;

   // init trackEW
   trackEW.nPoints = -1;

   if(nPoints != nStations)
      return false;

   trackEW = makeTrack(rpsBits, RPStudy);

   rpsBitsEast = rpsBits;
   rpsBitsWest = rpsBits;
   for (unsigned int iRP = 0; iRP < nStations; ++iRP)
   {
      rpsBitsWest.reset(iRP);
      rpsBitsEast.reset(iRP + 4);
   }
   trackW = makeTrack(rpsBitsWest, RPStudy);
   trackE = makeTrack(rpsBitsEast, RPStudy);

   return (IsColinear() && IsInDCut() && true); 
   //return (IsColinear() && IsInDCut() && InFiducial()); // The same fiducial region as for RHICf elastic analysis
}

bool ElasticAna::InFiducial()
{
   double phiE = trackE.phi;
   double phiW = trackW.phi;  
   return (IsInGeoWindow( phiE ) && IsInGeoWindow( phiW ));
}//InFiducial

bool ElasticAna::IsInDCut()
{
   double dx0 = trackE.X0 - trackW.X0;
   double dy0 = trackE.Y0 - trackW.Y0;
   return sqrt( dx0*dx0 + dy0*dy0 ) <= DCUT; // meters
}//IsInDCut

bool ElasticAna::IsColinear()
{
   double thetaEast = sqrt( trackE.thX*trackE.thX + trackE.thY*trackE.thY );
   double thetaWest = sqrt( trackW.thX*trackW.thX + trackW.thY*trackW.thY );
   double deltaTheta = thetaEast - thetaWest;
   
   return abs(deltaTheta) <= nSigma*sigma;
}//IsColinear


ElasticAna::track_t ElasticAna::makeTrack(bitset<8> rpsBits, int copyRP){
   double vertex[3] = {vertex_X, vertex_Y, vertex_Z}; // to reconstruct local tracks = with only 1 RP 
   track_t track;
   double psqr = mUtil->p0()*mUtil->p0(); 
   const unsigned int mirrorRP[nRomanPots] = { W1D, W1U, W2D, W2U, E1D, E1U, E2D, E2U};

   unsigned int nPoints = 0;
   StUPCRpsTrackPoint *trkPoint;
   for (int iRP = 0; iRP < nRomanPots; ++iRP)
   {
      if(!rpsBits.test(iRP))
         continue;

      trkPoint = nullptr;
      if(iRP != copyRP)
         trkPoint = mRpEvt->getTrackPoint(mTrackPointIdVec[iRP][0]);
      else 
         trkPoint = mRpEvt->getTrackPoint(mTrackPointIdVec[mirrorRP[iRP]][0]);

      track.posX[nPoints] = trkPoint->x();
      track.posY[nPoints] = trkPoint->y();
      track.posZ[nPoints] = trkPoint->z();
      nPoints++;
   }


   double xT, yT, thx, thy;
   double sxz=0., sx=0., sz=0., szz=0.; 
   double syz=0., sy=0.;
   for(unsigned int iP = 0; iP < nPoints; ++iP)
   {
      sx += track.posX[iP];
      sy += track.posY[iP];
      sz += track.posZ[iP];
      sxz += track.posX[iP]*track.posZ[iP];
      syz += track.posY[iP]*track.posZ[iP];
      szz += track.posZ[iP]*track.posZ[iP];
   }     
   double fp= (double)nPoints;
   if(nPoints == 1){
      thx = (track.posX[0]-vertex[0])/(track.posZ[0]-vertex[2]);
      thy = (track.posY[0]-vertex[1])/(track.posZ[0]-vertex[2]);
      xT = vertex[0]; 
      yT = vertex[1];
   }else{
      thx = (fp*sxz-sx*sz)/(fp*szz-sz*sz);
      thy = (fp*syz-sy*sz)/(fp*szz-sz*sz);    
      xT = sx/fp - thx*sz/fp;
      yT = sy/fp - thy*sz/fp;
   }

   track.thX = thx;
   track.thY = thy;
   track.phi = atan2(thy, thx);
   track.t = psqr*(thx*thx+thy*thy);
   track.X0  = xT;
   track.Y0  = yT;
   track.nPoints = nPoints;

   return track;
}


void ElasticAna::FillClusterInfo()
{
   // Fill cluster information
   int nClusters[] = {0, 0, 0, 0, 0, 0, 0, 0};
   int nClustersPerPlain[nRomanPots][nPlanes];
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
         nClustersPerPlain[iRp][iPlane] = 0;

   for (unsigned int i = 0; i < mRpEvt->getNumberOfClusters(); ++i)
   {
      StUPCRpsCluster *clstr = mRpEvt->getCluster(i);
      int clusterRP = clstr->romanPotId();
      int clusterPlane = clstr->planeId();
      hClusterEnergy[clusterRP]->Fill( clstr->energy() );
      hClusterLength[clusterRP]->Fill( clstr->length() );
      nClustersPerPlain[clusterRP][clusterPlane]++;
      nClusters[clusterRP]++;
   }

   for (int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      hNClusterPerRP[iRp]->Fill(nClusters[iRp]);
      for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
         hClusterPerPlain[iRp][iPlane]->Fill(nClustersPerPlain[iRp][iPlane]);
   }
}


void ElasticAna::closerTest(const vector<TVector3> trackPoints, int arm)
{
   const int ArmToRP[2][4] = { {0, 2, 5, 7}, {1, 3, 4, 6}};
   int rpStudy;
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      for (int iRP = 0; iRP < 2; ++iRP)
      {
         /// Make fit
         rpStudy = iRP + 2*iSide;
         TVector3 TPProjection = mUtil->fitLine(trackPoints, trackPoints[rpStudy].Z());

         hDeltaFitProjection[ArmToRP[arm][rpStudy]][X]->Fill((trackPoints[rpStudy].X() - TPProjection.X())*1000);
         hDeltaFitProjection[ArmToRP[arm][rpStudy]][Y]->Fill((trackPoints[rpStudy].Y() - TPProjection.Y())*1000);        
      }
   }
}


void ElasticAna::runVertexExtraction()
{
   unsigned int newRunNumber = mUpcEvt->getRunNumber();
   if(newRunNumber != mRunNumber){
      mRunNumber = newRunNumber;
      for (int iArm = 0; iArm < nArms; ++iArm)
      {
         for (int iCoor = 0; iCoor < nCoordinates - 1; ++iCoor)
            hVertexExtraction[iArm][iCoor]->SetName("hVertexExtraction_" + mUtil->armName(iArm) + "_" + mUtil->coordinateName(iCoor) + Form("_%i",mRunNumber) );
         hVertexZExtraction[iArm]->SetName("hVertexExtraction_" + mUtil->armName(iArm) + "_Z"  + Form("_%i",mRunNumber) );
      }
   }  


   int arm;
   double yHit[nSides];
   const int BranchToArm[] = { 0, 1, 1, 0};

   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack *trk = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]);
      StUPCRpsTrackPoint *trkPointFirst = trk->getTrackPoint(0);

      if( !trkPointFirst)
      {
         cout<<"ElasticAna::runVertexExtraction loaded nullptr trackpoint...  ERROR"<<endl;
         return;
      }
      arm = BranchToArm[trk->branch()]; 

      yHit[iSide] = trkPointFirst->y()*100;

      hVertexExtraction[arm][X]->Fill(trk->thetaRp(X)*1000, trkPointFirst->x()*1000);
      hVertexExtraction[arm][Y]->Fill(trk->thetaRp(Y)*1000, trkPointFirst->y()*1000);
   }
   hVertexZExtraction[arm]->Fill(1578*(0.2 - (yHit[W] + yHit[E]))/(yHit[W] - yHit[E]));
}

void ElasticAna::runMCValidation()
{
   int arm;
   vector<TVector3> trackPoints;
   const int BranchToArm[] = { 0, 1, 1, 0};
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack *trk = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]);
      for (int iRp = 0; iRp < nStationPerSide; ++iRp)
      {
         StUPCRpsTrackPoint *trkPoint = trk->getTrackPoint(iRp);
         if( trkPoint == nullptr)
         {
            cout<<"ElasticAna::runMCValidation loaded nullptr trackpoint...  ERROR"<<endl;
            return;
         }
         trackPoints.push_back( trkPoint->positionVec() );
      }
      arm = BranchToArm[trk->branch()]; 

      // vertex assumed to be in [0,0,0]
      double projX = trackPoints[0 + 2*iSide].x() + abs(trackPoints[1 + 2*iSide].z() - trackPoints[0 + 2*iSide].z())*((trackPoints[0 + 2*iSide].x() + 0.001)/abs(trackPoints[0 + 2*iSide].z()));
      double projY = trackPoints[0 + 2*iSide].y() + abs(trackPoints[1 + 2*iSide].z() - trackPoints[0 + 2*iSide].z())*((trackPoints[0 + 2*iSide].y() - 0.0011)/abs(trackPoints[0 + 2*iSide].z()));        

      hDeltaProjection[X][trk->branch()]->Fill((trackPoints[1 + 2*iSide].x() - projX)*1000);
      hDeltaProjection[Y][trk->branch()]->Fill((trackPoints[1 + 2*iSide].y() - projY)*1000);
   }
   closerTest(trackPoints, arm);    
}

void ElasticAna::FillAccaptancePlots()
{
   for (int iBr = 0; iBr < nBranches; ++iBr)
   {
      if( mRpTrackIdVec_perBranch[iBr].size()!=1 )
         continue;

      StUPCRpsTrack* trk = mRpEvt->getTrack(mRpTrackIdVec_perBranch[iBr][0]);
      if(!trk)    
      {
         cout<<"ElasticAna::FillAccaptancePlots loaded nullptr track...  ERROR"<<endl;
         return;
      }
      hAcceptancePxPy[iBr]->Fill(trk->pVec().X(), trk->pVec().Y());
   }

   for (int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      if( mTrackPointIdVec[iRp].size() !=1 )
         continue;

      StUPCRpsTrackPoint *trkPoint = mRpEvt->getTrackPoint(mTrackPointIdVec[iRp][0]);
      if(!trkPoint) 
      {
         cout<<"ElasticAna::FillAccaptancePlots loaded nullptr trackpoint...  ERROR"<<endl;
         return;
      }
      hAcceptanceXY[iRp]->Fill(trkPoint->x(),trkPoint->y());   
   }
}

const TString ElasticAna::mCutName[nCuts] = { TString("All"), TString("ET"), TString("2 RP trks"),
               TString("In RP fiducial"), TString("Collinearity")};

const TString ElasticAna::mElAnaCutsName[kMax] = { TString("All"), TString("El trig"), TString("ET pattern"),
               TString("4 points"), TString("t-range") ,TString("Collinear"), TString("Geo Window"), 
               TString("dCut"), TString("elastic")};


void ElasticAna::runAlignment()
{
   vector<TVector3> trackPoints;
   vector<unsigned int> rpIdVec;
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      if( mTrackPointIdVec[iRp].size() !=1 )
         continue;

      StUPCRpsTrackPoint *trkPoint = mRpEvt->getTrackPoint(mTrackPointIdVec[iRp][0]);
      if(!trkPoint)
      {
         cout<<"ElasticAna::runAlignment loaded nullptr trackpoint...  ERROR"<<endl;
         return;
      }
      rpIdVec.push_back(iRp); 
      trackPoints.push_back( trkPoint->positionVec() ); 
   }

   if( trackPoints.size() != nStations)
   {
      cout<<"ERROR in ElasticAna::runAlignment not elastic event"<<endl;
      return;
   }

   for (unsigned int iTp = 0; iTp < trackPoints.size(); ++iTp)
   {
      TVector3 TP = mUtil->fitLine(trackPoints, trackPoints[iTp].Z());

      alignment[ rpIdVec[iTp] ][X].push_back(trackPoints[iTp].X() - TP.X());
      alignment[ rpIdVec[iTp] ][Y].push_back(trackPoints[iTp].Y() - TP.Y()); 
   }  
}

TVector3 ElasticAna::CalculateOffsetCorrForRP(unsigned int Rp, unsigned int iter)
{
   TVector3 rpCorrection(99,99,0);
   unsigned int nCorrevtiveEvents = alignment[Rp][X].size();
   if( nCorrevtiveEvents == 0)
      return rpCorrection;

   rpCorrection.SetXYZ(0.,0.,0.);
   // calculate the average correction per run
   for (unsigned int id = 0; id < nCorrevtiveEvents; ++id)
   {
      rpCorrection[X] += alignment[Rp][X][id];
      rpCorrection[Y] += alignment[Rp][Y][id];
   }
   rpCorrection[X] = rpCorrection[X]/nCorrevtiveEvents;
   rpCorrection[Y] = rpCorrection[Y]/nCorrevtiveEvents;
   rpCorrection[Z] = nCorrevtiveEvents;

   hAligXCorr[Rp][iter]->Fill(rpCorrection[X]*1000);
   hAligYCorr[Rp][iter]->Fill(rpCorrection[Y]*1000);
   hAligEvents[Rp][iter]->Fill(nCorrevtiveEvents);

   alignment[Rp][X].clear();
   alignment[Rp][Y].clear();

   return rpCorrection;
}


void ElasticAna::InitRPMCInfo()
{
   mRecTree->InitRPMCInfo();
}

void ElasticAna::SetRPMCInfo(double *mc_vrtx, double (*mc_p)[nSides])
{
   mRecTree->setTrueRPVertexZInCm( mc_vrtx[Z] );
   mRecTree->setTrueRPVertexXInCm( mc_vrtx[X] );
   mRecTree->setTrueRPVertexYInCm( mc_vrtx[Y] );
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      mRecTree->setRPTruePx( mc_p[X][iSide], iSide);
      mRecTree->setRPTruePy( mc_p[Y][iSide], iSide);
      mRecTree->setRPTruePz( mc_p[Z][iSide], iSide);
   }

}

void ElasticAna::FillElasticPlots()
{
   if(trackEW.nPoints == 0)
      return;

   hElasticAnaFlow->Fill(kETPattern);
   if(trackEW.nPoints < 0)
      return;

   hElasticAnaFlow->Fill(kFourPoints);
   if( trackEW.t < 0.23 )//|| trackEW.t > 0.67) // not use the cut and plot the N(t)
      return;
   hElasticAnaFlow->Fill(kTRange);

   bool inFiducial = InFiducial();
   bool isColinear = IsColinear();
   bool isInDCut = IsInDCut();
   if(inFiducial)
      hElasticAnaFlow->Fill(kFiducial);
   if(isInDCut)
      hElasticAnaFlow->Fill(kDCut);
   if(isColinear)
      hElasticAnaFlow->Fill(kColinear);

   hThetaCorr->Fill(trackE.thX - trackW.thX,  trackE.thY - trackW.thY);
   hTheta[E][X]->Fill(trackE.thX);
   hTheta[E][Y]->Fill(trackE.thY);
   hTheta[W][X]->Fill(trackW.thX);
   hTheta[W][Y]->Fill(trackW.thY);


   double thetaEast = sqrt( trackE.thX*trackE.thX + trackE.thY*trackE.thY );
   double thetaWest = sqrt( trackW.thX*trackW.thX + trackW.thY*trackW.thY );
   double deltaTheta = thetaEast - thetaWest;
   hDeltaTheta[0]->Fill( deltaTheta );
   
   double dx0 = trackE.X0 - trackW.X0;
   double dy0 = trackE.Y0 - trackW.Y0;

   if(isColinear && inFiducial)
   {
      hDCut->Fill( dx0, dy0 );
      hDCutR->Fill(sqrt( dx0*dx0 + dy0*dy0 ));
   }
   if(isColinear && isInDCut && inFiducial)
   {
      hDeltaTheta[1]->Fill( deltaTheta ); 
      hT->Fill(trackEW.t);  
   }


}

void ElasticAna::runRpDataDrivenEffStudy()
{
   for (unsigned int iRP = 0; iRP < nRomanPots; ++iRP)
   {
      if(!IsElasticEventCandidate(iRP))
         continue;
      //if( trackEW.t < 0.23 || trackEW.t > 0.67)
      //   continue;
      FillEffPlots( iRP, checkRP(iRP)); 
   }
}//runRpDataDrivenEffStudy

bool ElasticAna::checkRP(unsigned int iRP)
{
   if(mTrackPointIdVec[iRP].size() != 1)
      return false;

   StUPCRpsTrackPoint *trkPoint = mRpEvt->getTrackPoint(mTrackPointIdVec[iRP][0]);
   double trueX = trkPoint->x();
   double trueY = trkPoint->y();

   double fitX, fitY;
   getFitPoint(fitX, fitY, trackEW, iRP );

   double distInX = trueX - fitX;
   double distInY = trueY - fitY;
   double dist = sqrt(distInX*distInX+distInY*distInY );
   hFindHitDiff->Fill( dist);
   return dist <= maxDistBetweenFitAndHit; //true;
}//checkRP

void ElasticAna::getFitPoint(double& x, double& y, track_t track, unsigned int rp)
{
  x = track.X0 +  track.thX*mUtil->rpZPosition(rp);
  y = track.Y0 +  track.thY*mUtil->rpZPosition(rp);
}//getFitPoint

void ElasticAna::FillEffPlots( unsigned int iRP, bool passed)
{
   TVector3 momentum;

   track_t trk = iRP < nStations ? trackE : trackW;

   double theta = sqrt( trk.thX*trk.thX + trk.thY*trk.thY );
   momentum.SetMagThetaPhi(mUtil->p0(), theta, trk.phi);

   double fitX, fitY;
   getFitPoint(fitX, fitY, trackEW, iRP );
   fitX*=100; fitY*=100; // conevrt from meters to cm
   int branch = mUtil->branchPerRp( iRP );
   for (int i = 0; i < (passed ? 2 : 1);++i)
   {
      hEffPxPy[i][branch]->Fill(momentum.x(), momentum.y());
      hEffXY[i][iRP]->Fill(fitX, fitY);
      hEffXYOffSub[i][iRP]->Fill(fitX, fitY-((mOffSet[iRP][Y]+mCorrection[iRP][Y])*100));
   }
}//FillEffPlots
