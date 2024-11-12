#include "PlotManager.h"

void PlotManager::runDsmEffStudy(int study)
{
   TString dataSetName = studyName[study];
   changeDir(study, dataSetName);

   TString dsmEffName[nDSMEFF] = { TString("RP ET"), TString("RP IT"), TString("RP ET Veto"), TString("RP IT Veto"), TString("TOF(1)"), 
      TString("TOF(10)"), TString("BBC East"), TString("BBC West"), TString("BBCL East"), TString("BBCL West"), TString("ZDC East"), TString("ZDC West")};;

   hDSMEff[0] = new TH1D("hDSMEff_total_"+dataSetName, "hDSMEff_total", nDSMEFF-1, 1, nDSMEFF);
   hDSMEff[1] = new TH1D("hDSMEff_passed_"+dataSetName, "hDSMEff_passed", nDSMEFF-1, 1, nDSMEFF);

   CheckDSMEff( mRecTree[study] );
   if(!TEfficiency::CheckConsistency(*hDSMEff[1],*hDSMEff[0])){
      cout<<"Error in printEfficiencies: CheckConsistency()"<<endl;
      return;
   }

   TEfficiency* pEff = new TEfficiency(*hDSMEff[1],*hDSMEff[0]);

   CreateCanvas(&canvas,"DSMEff_"+dataSetName);
   gPad->SetMargin(0.13,0.03,0.23,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   TGraphAsymmErrors *gr = pEff->CreateGraph();
   SetGraphStyle(gr);
   gr->SetTitle(";;DSM bit eff");
   gr->GetYaxis()->SetTitleOffset(2.1);
   for(int tb=1; tb<=nDSMEFF; ++tb)
      gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(tb+0.5) , dsmEffName[tb-1]);

   gr->Draw("ap");

   DrawSTARInternal(0.17, 0.79, 0.4, 0.84);
   DrawDataSetInfo(dataSetName);

   WriteCanvas("DSMEff_"+dataSetName);
   canvas->Close();   
}//runDsmEffStudy

void PlotManager::DrawDataSetInfo( TString dataSet )
{

   if( dataSet == studyName[kTRIGEFF]){
      CreateText(0.17,0.74,0.4,0.79);
      text -> AddText("Zero-bias data");   
      text -> Draw("same");
   }

   if( dataSet == studyName[kFULLZB]){
      CreateText(0.17,0.64,0.4,0.79);
      text -> AddText("Zero-bias data");   
      text -> AddText("1 TOF prim. vrtx");
      text -> AddText("2 TOF prim. trks");
      text -> Draw("same");
   }
}//DrawDataSetInfo

void PlotManager::CheckDSMEff(RecTree* recTree)
{
   TTree *tree = recTree->getTTree();

   for(Long64_t iev=0; iev<tree->GetEntries(); ++iev)
   { //get the event
      tree->GetEntry(iev); 

      if( recTree->getRpEtTrigBit())
         hDSMEff[0]->Fill( kRPET );
      if( recTree->getRpItTrigBit())
         hDSMEff[0]->Fill( kRPIT );
      if( !recTree->getRpEtTrigBit())
         hDSMEff[0]->Fill( kRPETVETO );
      if( !recTree->getRpItTrigBit())
         hDSMEff[0]->Fill( kRPITVETO );
      if( recTree->getTofMult() > 1)
         hDSMEff[0]->Fill( kTOFA );
      if( recTree->getTofMult() <= 10)
         hDSMEff[0]->Fill( kTOFB );
      if( recTree->getBbcSmallEast() <= 20)
         hDSMEff[0]->Fill( kBBCE );
      if( recTree->getBbcSmallWest() <= 20)
         hDSMEff[0]->Fill( kBBCW );
      if( recTree->getBbcLargeEast() <= 50)
         hDSMEff[0]->Fill( kBBCLE );
      if( recTree->getBbcLargeWest() <= 50)
         hDSMEff[0]->Fill( kBBCLW );
      if( !recTree->getZdcETrigBit() )
         hDSMEff[0]->Fill( kZDCE );
      if( !recTree->getZdcWTrigBit() )
         hDSMEff[0]->Fill( kZDCW );

      if( recTree->getRpEtDsmBit() && recTree->getRpEtTrigBit())
         hDSMEff[1]->Fill( kRPET );
      if( recTree->getRpItDsmBit() && recTree->getRpItTrigBit())
         hDSMEff[1]->Fill( kRPIT );
      if( !recTree->getRpEtDsmBit() && !recTree->getRpEtTrigBit())
         hDSMEff[1]->Fill( kRPETVETO );
      if( !recTree->getRpItDsmBit() && !recTree->getRpItTrigBit())
         hDSMEff[1]->Fill( kRPITVETO );
      if( recTree->getTofDsmABit() && recTree->getTofMult() > 1)
         hDSMEff[1]->Fill( kTOFA );
      if( !recTree->getTofDsmBBit() && recTree->getTofMult() <= 10)
         hDSMEff[1]->Fill( kTOFB );
      if( !recTree->getBbcSmallEDsmBit() && recTree->getBbcSmallEast() <= 20)
         hDSMEff[1]->Fill( kBBCE );
      if( !recTree->getBbcSmallWDsmBit() && recTree->getBbcSmallWest() <= 20)
         hDSMEff[1]->Fill( kBBCW );
      if( !recTree->getBbcLargeEDsmBit() && recTree->getBbcLargeEast() <= 50)
         hDSMEff[1]->Fill( kBBCLE );
      if( !recTree->getBbcLargeWDsmBit() && recTree->getBbcLargeWest() <= 50)
         hDSMEff[1]->Fill( kBBCLW );
      if( !recTree->getZdcEDsmBit() && !recTree->getZdcETrigBit())
         hDSMEff[1]->Fill( kZDCE );
      if( !recTree->getZdcWDsmBit() && !recTree->getZdcWTrigBit())
         hDSMEff[1]->Fill( kZDCW );
   }
}//CheckDSMEff

void PlotManager::runTofTrigStudy()
{
   unsigned int nPassed, nTotal;      
   nTotal = mTree[kTRIGEFF]->GetEntries(); 
   nPassed = -1;
   mTree[kTRIGEFF]->Draw("tofMult","tofDsmABit");
   TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
   if(htemp) nPassed = htemp->GetEntries();
   
   CalculateEfficiency("TOF trig eff = ", nPassed, nTotal);
}//runTofTrigStudy

void PlotManager::runRpTrigStudy()
{
   changeDir(kTRIGEFF, "rpTrigStudy");
   TH1D *hist[2];
   TH2D *hist2D[2][nBranches];
   TH2D *hEff[2];

   TString histDir = "rpTrigEff";
   TString histName;
   TString histLabel[2] = { TString("Total"), TString("Passed")};
   TString histTotalName;

   TCanvas* effCanvas[nSides];

   const int nSamples = 2;
   TString sampleTag[] = { TString(""), TString("_PxPyCut_")};

   for (int iSample = 0; iSample < nSamples; ++iSample)
   {
      CreateCanvas(&canvas,"rpTrigEffPtMiss"+sampleTag[iSample]);
      SetGPad(.0,.0,.0,.0);
      canvas->Divide(nRpConfigurations,1,.0,.0,0);
      CreateCanvas(&effCanvas[East],"rpTrigEffPxPy_"+mUtil->sideName(E)+sampleTag[iSample]);
      SetGPad(.0,.0,.0,.0);
      effCanvas[East]->Divide(nRpConfigurations,1,.0,.0,0);
      CreateCanvas(&effCanvas[West],"rpTrigEffPxPy_"+mUtil->sideName(W)+sampleTag[iSample]);
      SetGPad(.0,.0,.0,.0);
      effCanvas[West]->Divide(nRpConfigurations,1,.0,.0,0);

      for (int iConf = 0; iConf < nRpConfigurations; ++iConf)
      {
         // Get histograms
         histName = "hRpTrigEffPtMiss";
         for (int i = 0; i < 2; ++i){
            histTotalName = histDir + "/" + histName + "_" + histLabel[i] + "_" + sampleTag[iSample] + mUtil->rpConfigName(iConf);
            hist[i] = (TH1D*)inFile->Get( histTotalName );
            if(!hist[i]){
               cout<<"Error in runRpTrigStudy: cannot loaded "<<histTotalName<<endl; 
               return;
            }
            hist[i]->Rebin(5);
         }
         if(!TEfficiency::CheckConsistency(*hist[1],*hist[0])){
            cout<<"Error in runRpTrigStudy: CheckConsistency()"<<endl;
            return;
         }
         TEfficiency* pEff = new TEfficiency(*hist[1],*hist[0]);
         canvas->cd(iConf+1);
         SetGPad(0.08,0.01,0.13,0.01); //(Float_t left, Float_t right, Float_t bottom, Float_t top)
         TGraphAsymmErrors *gr = pEff->CreateGraph();
         SetGraphStyle(gr);
         gr->SetTitle(";p_{T}^{RP,miss} [GeV];#varepsilon_{RP}^{trig}");
         gr->GetYaxis()->SetTitleOffset(0.6);
         gr->GetXaxis()->SetTitleOffset(1.1);
         gr->Draw("ap");
         if( iSample == 1)
            CalculateEfficiency("RP trig eff for " + mUtil->rpConfigName(iConf) + " comb = ", hist[1]->GetEntries(), hist[0]->GetEntries());

         for (int iBr = 0; iBr < nBranches; ++iBr){
            // Get histograms
            histName = "hRpTrigEffPxPy";
            for (int i = 0; i < 2; ++i){
               histTotalName = histDir + "/" + histName + "_" + sampleTag[iSample] + mUtil->branchName(iBr) + "_" + histLabel[i] + "_" + mUtil->rpConfigName(iConf);
               hist2D[i][iBr] = (TH2D*)inFile->Get( histTotalName );
               if(!hist2D[i][iBr]){
                  cout<<"Error in runRpTrigStudy: cannot loaded "<<histTotalName<<endl; 
                  return;
               }
            }
         }
         for (int iSide = 0; iSide < nSides; ++iSide)
         {
            hEff[0] = (TH2D*)hist2D[1][iSide*2]->Clone("hRpTrigEff_" + sampleTag[iSample] + mUtil->sideName(iSide) + mUtil->rpConfigName(iConf));
            hEff[0]->Add(hist2D[1][iSide*2 + 1]);
            hEff[1] = (TH2D*)hist2D[0][iSide*2]->Clone("hRpTrigEff_total_" + sampleTag[iSample] + mUtil->sideName(iSide) + mUtil->rpConfigName(iConf));
            hEff[1]->Add(hist2D[0][iSide*2 + 1]);

            hEff[0]->Divide(hEff[1]);
            effCanvas[iSide]->cd(iConf+1); 
            SetGPad(0.14,0.16,0.12,0.02);
            SetHistStyle(hEff[0]);
            hEff[0]->GetXaxis()->SetTitleOffset(1.0);
            hEff[0]->GetYaxis()->SetTitleOffset(1.1);
            hEff[0]->GetZaxis()->SetTitleOffset(0.8);
            hEff[0]->SetTitle(";p_{x} [GeV];p_{y} [GeV]; #varepsilon_{RP}^{trig}");
            hEff[0]->SetNdivisions(508, "X");
            hEff[0]->Draw("colz");
            CreateText(0.55,0.55,0.55,0.55);
            text -> AddText(mUtil->stationName(2*iSide + iConf)); 
            text -> SetTextSize(textSize*2);  
            text -> Draw("same");
            DrawFiducial(iSide);
        
         }
      }    
      WriteCanvas("rpTrigEffPtMiss"+sampleTag[iSample]);
      canvas->Close();

      for (int iSide = 0; iSide < nSides; ++iSide)
      {
         WriteCanvas("rpTrigEffPxPy_"+mUtil->sideName(iSide)+sampleTag[iSample], effCanvas[iSide]);
         effCanvas[iSide]->Close();
      }
   }
}//runRpTrigStudy

void PlotManager::CalculateEfficiency(TString effString, unsigned int nPassed, unsigned int nTotal)
{
   if( DEBUG )
      cerr<<"PlotManager::CalculateEfficiency() just started"<<endl;

   if( nTotal == 0)
      return;

   TH1D *hTOFEff[2];  
   hTOFEff[0] = new TH1D("hTOFEff_total", "hTOFEff_total", 1, 1, 2);
   hTOFEff[1] = new TH1D("hTOFEff_passed", "hTOFEff_passed", 1, 1, 2);

   hTOFEff[0]->SetBinContent(1, nTotal);
   hTOFEff[1]->SetBinContent(1, nPassed);

   if(!TEfficiency::CheckConsistency(*hTOFEff[1],*hTOFEff[0])){
      cout<<"Error in CalculateEfficiency: CheckConsistency()"<<endl;
      return;
   }
   TEfficiency* pEff = new TEfficiency(*hTOFEff[1],*hTOFEff[0]);
   TGraphAsymmErrors *gr = pEff->CreateGraph();

   cout<<effString<<nPassed<<" / "<<nTotal<<" = "<< gr->GetY()[0]<<" + "<<gr->GetEYhigh()[0]<<" - "<<gr->GetEYlow()[0]<<endl;

   delete hTOFEff[0];
   delete hTOFEff[1];
}//CalculateEfficiency