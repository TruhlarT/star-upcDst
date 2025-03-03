#include "PlotManager.h"

void PlotManager::runMainAnaPlots()
{
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to runMainAnaPlots"<<endl;
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots call initMassPlots()"<<endl;
   initMainAnaPlots();
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to runM2Plots"<<endl;
   changeDir(kMAINANA, "PIDplots");
   runM2Plots();
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to drawMissIdProbability"<<endl;
   drawMissIdProbability();
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to runAnaPlots"<<endl;
   runAnaPlots(); 
   //return;
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to runMissingPt"<<endl;
   runMissingPt();
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots going to loop through tree"<<endl;

   mCurrentTree = mRecTree[kMAINANA];
   for(Long64_t iev=0; iev<mTree[kMAINANA]->GetEntries(); ++iev) // Nominal analysis
   { //get the event
      mTree[kMAINANA]->GetEntry(iev); 
      loopThroughTree(NOMINAL, NOMINAL, NOMINAL, NOMINAL, NOMINAL, NOMINAL, 0); 
      for (unsigned int i = 1; i <= nTPCAppSysStudies; ++i)
         loopThroughTree(i, NOMINAL, NOMINAL, NOMINAL, NOMINAL, NOMINAL, i);
      for (unsigned int i = 1; i <= nTOFAppSysStudies; ++i)
         loopThroughTree(NOMINAL, i, NOMINAL, NOMINAL, NOMINAL, NOMINAL, nTPCAppSysStudies + i);
   }
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots() going to runSysStudy"<<endl;
   runSysStudy();

   //return;

   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots call correctMainPlotsForBinWidth()"<<endl;
   correctMainPlotsForBinWidth();
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots call subtractBackgroundFromMainPlots()"<<endl;
   subtractBackgroundFromMainPlots();
/*
   changeDir(kMAINANA, "histos");
   for (size_t iGroup = 0; iGroup < mainAnaHists.size(); ++iGroup)
      for (unsigned int i = 0; i <= nTotalTPCTOFSysStudies; ++i){
         mainAnaHists[iGroup].hMainBcgSubtracted[i]->Write();
         mainAnaHists[iGroup].hMain[i]->Write();
      }
*/
   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots call PlotSysStudy()"<<endl;
   changeDir(kMAINANA, "SysPlots");
   PlotSysStudy();

   if( DEBUG )
      cerr<<"PlotManager::runMainAnaPlots call PlotMassPlots()"<<endl;
   PlotMainAnaPlots();

   if( DEBUG )
      cout<<"PlotManager::runMainAnaPlots() finnished"<<endl;
}


void PlotManager::initMainAnaPlots()
{
   double *pTMissBins = setBinArray(91, 0.0, 0.9); 
   for (int iPart = 0; iPart < nParticles; ++iPart)
   {
      for (int iRpCon = 0; iRpCon <= nRpConfigurations; ++iRpCon)
      {
         TString rpConName = iRpCon < nRpConfigurations ? mUtil->rpConfigTag(iRpCon) : "";
         TString hNameSuffix = mUtil->particleName(iPart);
         if(iRpCon < nRpConfigurations)
            hNameSuffix += "_" + mUtil->rpConfigName(iRpCon);


         hGroup group;
         group.part = iPart;
         group.rpCon = iRpCon;

         for (unsigned int iSys = 0; iSys <= nTotalTPCTOFSysStudies; ++iSys)
         {
            TString hSuffix = iSys == 0 ? "" : Form("_s%i",iSys);
            TString unit = iPart == PROTON ? "pb" : "nb";
            hInvMassCorr[iPart][iRpCon][iSys] = new TH1F("hInvMass"+hNameSuffix+hSuffix, ";m(" + mUtil->pairLabel(iPart) + ") [GeV];d#sigma/dm(" + mUtil->pairLabel(iPart) + ") [" + unit +"/GeV]", massNBins[iPart]-1,massBins[iPart]);//, int(mRuns.size()-1), mRuns.data());
            hInvMassCorr[iPart][iRpCon][iSys]->Sumw2();
            hMissingPtVsInvMass[iPart][iRpCon][iSys] = new TH2F("hMissingPtVsInvMass"+hNameSuffix+hSuffix, rpConName + ";m(" + mUtil->pairLabel(iPart) + ") [GeV];p_{T}^{miss} GeV;d#sigma/dm(" + mUtil->pairLabel(iPart) + ") [" + unit +"/GeV]", massNBins[iPart]-1,massBins[iPart], 90, pTMissBins);//, int(mRuns.size()-1), mRuns.data());
            hMissingPtVsInvMass[iPart][iRpCon][iSys]->Sumw2();

            hDeltaPhi[iPart][iRpCon][iSys] = new TH1F("hDeltaPhi"+hNameSuffix+hSuffix,";#Delta#varphi [deg];d#sigma/d#Delta#varphi [pb/deg]", 18, 0.1, 180);
            hDeltaPhi[iPart][iRpCon][iSys]->Sumw2();    
            hMissingPtVsDeltaPhi[iPart][iRpCon][iSys] = new TH2F("hMissingPtVsDeltaPhi"+hNameSuffix+hSuffix, rpConName + ";#Delta#varphi [deg];p_{T}^{miss} GeV;d#sigma/d#Delta#varphi [pb/deg]", 18, 0, 180, 91, 0.0, 0.9);//, int(mRuns.size()-1), mRuns.data());
            hMissingPtVsDeltaPhi[iPart][iRpCon][iSys]->Sumw2();

         }
         group.hMain = hDeltaPhi[iPart][iRpCon];
         group.hVsPt = hMissingPtVsDeltaPhi[iPart][iRpCon];
         group.hMainBcgSubtracted = new TH1F*[nTotalTPCTOFSysStudies+1];
         group.hName = "DeltaPhi";
         group.xLabel = "#Delta#varphi [deg]";
         group.yLabel = "#Delta#varphi";
         mainAnaHists.push_back(group);

         group.hMain = hInvMassCorr[iPart][iRpCon];
         group.hVsPt = hMissingPtVsInvMass[iPart][iRpCon];
         group.hMainBcgSubtracted = new TH1F*[nTotalTPCTOFSysStudies+1];
         group.hName = "InvMass";
         group.xLabel = "m(" + mUtil->pairLabel(iPart) + ") [GeV]";
         group.yLabel = "m(" + mUtil->pairLabel(iPart) + ")";
         mainAnaHists.push_back(group);
      }
   }

   initGraniittiPlots();
   mCurrentTree = mRecTree[kMAINANA];
}//initMainAnaPlots

void PlotManager::loopThroughTree(UInt_t tpcAppMethod, UInt_t tofAppMethod, UInt_t nHitsFit, UInt_t nHitsDeDx, UInt_t pidMethod, UInt_t dcaMethod, int id)
{
   mCurrentTree->CalculatePID(true, true, -1, pidMethod);
   int pairID = mCurrentTree->getPairID(); 
   if( pairID < 0)
      return;

   double totalEff = 1;

   //totalEff *= getRetainCorrection(); // already on luminosity correction
   totalEff *= getRpEff();
   totalEff *= getVertexCutEff();
   totalEff *= getVertexReconstructionEff(dcaMethod);
   totalEff *= getPtMissEff();
   totalEff *= getPIDEff(pidMethod);
   totalEff *= getLuminosity();
   totalEff *= getTpcEff(tpcAppMethod, nHitsFit, nHitsDeDx);
   totalEff *= getTofEff(tofAppMethod);

   //cout<<Form("Eff: %f = %f*%f*%f*%f*%f*%f*%f*%f",id, totalEff, getRpEff(), getVertexCutEff(), getVertexReconstructionEff(), getPtMissEff(), getPIDEff(pidMethod), getLuminosity(), getTpcEff(tpcAppMethod, nHitsFit, nHitsDeDx), getTofEff(tofAppMethod))<<endl;

   double effpb = totalEff/1000; // convert from nb to pb
   double eff = pairID == PROTON ? effpb : totalEff; // convert from nb to pb for Protons

   double missingPt = mCurrentTree->getPtMissing();
   int rpComb = getRpCombination();
   double invMass = mCurrentTree->getInvMass();

   hMissingPtVsInvMass[pairID][rpComb][id]->Fill(invMass,missingPt, 1./eff);
   hMissingPtVsInvMass[pairID][nRpConfigurations][id]->Fill(invMass,missingPt, 1./eff);
   double deltaPhi = getRpDeltaPhi();
   hMissingPtVsDeltaPhi[pairID][rpComb][id]->Fill(deltaPhi,missingPt, 1./effpb);
   hMissingPtVsDeltaPhi[pairID][nRpConfigurations][id]->Fill(deltaPhi,missingPt, 1./effpb);

   if(missingPt > exclusivityCut)
      return;

   hInvMassCorr[pairID][rpComb][id]->Fill(invMass, 1./eff);
   hInvMassCorr[pairID][nRpConfigurations][id]->Fill(invMass, 1./eff);

   hDeltaPhi[pairID][rpComb][id]->Fill( deltaPhi, 1./effpb );
   hDeltaPhi[pairID][nRpConfigurations][id]->Fill( deltaPhi, 1./effpb );

}//loopThroughTree


void PlotManager::PlotMainAnaPlots()
{
   double scaleFactor[nParticles][nRpConfigurations + 1] = { {1.3, 1.6, 1.3} , {1.35, 2.2, 1.5}, {1.85, 1.8, 1.3}};
   for (size_t iGroup = 0; iGroup < mainAnaHists.size(); ++iGroup) 
   {
      int part = mainAnaHists[iGroup].part;
      int rpCon = mainAnaHists[iGroup].rpCon;
      TString hName = mainAnaHists[iGroup].hName;

      TString nameSuffix = mUtil->particleName(part);
      if(rpCon < nRpConfigurations)
         nameSuffix += mUtil->rpConfigName(rpCon);

      TH1F *hist = (TH1F*)mainAnaHists[iGroup].hMain[0]->Clone( "h" + hName + "_" + nameSuffix );
      TH2F *hMissPtVsX = (TH2F*)mainAnaHists[iGroup].hVsPt[0]->Clone( "h" + hName + "VsPtMiss_" + nameSuffix );

      TH1D *hBkgdHistogram = (TH1D*)backgroundStudy(hMissPtVsX, part, rpCon);

      changeDir(kMAINANA, hName);
      changeSubDir(nameSuffix);

      CreateCanvas(&canvas, hName + "_" + nameSuffix);
      SetHistStyle(hist);

      if(part == PROTON)
         hist->GetYaxis()->SetTitleOffset(1.45);

      hist->GetYaxis()->SetRangeUser(0, scaleFactor[part][rpCon]*hist->GetMaximum());
      hist->Draw("hist E1X0");

      SetHistBcgStyle(hBkgdHistogram);
      hBkgdHistogram->Draw("same hist E1X0");

      DrawSTAR(0.17,0.89,0.33,0.95);
      DrawMainText(part, rpCon);

      CreateLegend(0.64, 0.45, 0.9, 0.59);
      legend->AddEntry(hist, "Data","ple");
      legend->AddEntry(hBkgdHistogram, "Non-Exclusive background","ple");
      legend->Draw("same");

      mCurrDir->cd();
      WriteCanvas(hName + "_" + nameSuffix);
      integrateCrossSection(mainAnaHists[iGroup], hName + "_" + nameSuffix);
      /////////////////////////////////////
      TH1F *histSubtracted = (TH1F*)mainAnaHists[iGroup].hMainBcgSubtracted[0]->Clone( "h" + hName + "Sub_" + nameSuffix );
      double integral = histSubtracted->Integral("width");
      SetHistStyle(histSubtracted);
      //subtractBackground( hist, hBkgdHistogram );
      histSubtracted->GetYaxis()->SetRangeUser(0, scaleFactor[part][rpCon]*histSubtracted->GetMaximum());
      histSubtracted->Draw("AXIS");

      TH1F *histMC = GetGraniittiPlot(part, rpCon, hName, integral, false);
      histMC->Draw("same hist");

      DrawSTAR(0.17,0.89,0.33,0.95);
      DrawSystemDescription(part);
      //DrawCentralKin(part);
      DrawForwardProtonKin(rpCon);

      DrawSysUncertainty(mainAnaHists[iGroup], true);
      histSubtracted->Draw("same E1X0");

      if( hName == "DeltaPhi")
         CreateLegend(0.3, 0.65, 0.65, 0.8);
      else
         CreateLegend(0.64, 0.45, 0.9, 0.59);
      legend->AddEntry(histSubtracted, "Data","ple");
      legend->AddEntry(systUncertainty, "Systematic uncertainty", "f");
      legend->AddEntry(histMC, "GRANIITTI not scaled","l");
      //legend->AddEntry(systUncertainty[bcg], "Background sys. uncert.", "f");
      legend->Draw("same");

      histSubtracted->Draw("AXIS same");
      WriteCanvas(hName + nameSuffix + "bckgSubUnscaled");    
      /////////////////////////////////////
      TString canTag[] = { TString("bckgSubScaled"), TString("sq")};
      for (int i = 0; i < 2; ++i)
      {
         if( i == 1)
         {
            CreateCanvas(&canvas, hName + "_" + nameSuffix, 1165.0, 980.0);
            gPad->SetMargin( hName == "DeltaPhi" ? 0.11 : 0.09, 0.03, 0.11, 0.01); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
            
            if( hName != "DeltaPhi")
               histSubtracted->GetYaxis()->SetTitleOffset(1.13);
            else
               histSubtracted->GetYaxis()->SetTitleOffset(1.4);
         }
         histSubtracted->Draw("AXIS");
         systUncertainty->Draw("E5 same");
         histSubtracted->Draw("E1X0 same");
         
         TH1F *histMC2; 
         TString granTag = " not scaled";
         TString granTag2 = " not scaled";

         TH1F *histContinuum; 
         TH1F *histContinuum2;
         if( hName == "DeltaPhi" && rpCon == nRpConfigurations)
         {
            histMC = GetGraniittiPlot(part, IT, hName, integral, true);
            if(part != PROTON)
               histMC->Draw("same hist"); 
            histContinuum = GetGraniittiPlot(part, IT, hName, integral, true, true);
            histContinuum->Draw("same hist"); 

            granTag = getGraniittiSF(part,IT) == 1.0 ? granTag : TString::Format(" #times%.2f",getGraniittiSF(part,IT));

            histMC2 = GetGraniittiPlot(part, ET, hName, integral, true);
            if(part != PROTON)
               histMC2->Draw("same hist");  

            histContinuum2 = GetGraniittiPlot(part, ET, hName, integral, true, true);
            histContinuum2->Draw("same hist");  

            granTag2 = getGraniittiSF(part,ET) == 1.0 ? granTag2 : TString::Format(" #times%.2f",getGraniittiSF(part,ET));
         }else{
            histMC = GetGraniittiPlot(part, rpCon, hName, integral, true);
            if(part != PROTON)
               histMC->Draw("same hist");

            histContinuum = GetGraniittiPlot(part, rpCon, hName, integral, true, true);
            histContinuum->Draw("same hist");

            if( rpCon == nRpConfigurations && getGraniittiSF(part,ET) != getGraniittiSF(part,IT))
               granTag = " scaled";
            else
               granTag = getGraniittiSF(part,rpCon) == 1.0 ? granTag : TString::Format(" #times%.2f",getGraniittiSF(part,rpCon));
         }

         DrawSTAR( 0.14, 0.89, 0.24, 0.96);
         DrawSystemDescription(part, 0.24, 0.89, 0.88, 0.96);
         //DrawCentralKin(part);
         if( i == 0)
            DrawForwardProtonKin(rpCon);
         else
            DrawForwardProtonKin(rpCon, 0.33, 0.8, 0.87, 0.87);


         double xshift = i == 0 ? 0.0 : 0.17;
         if( hName == "DeltaPhi")
            CreateLegend(0.3, part != PROTON ? 0.45 : 0.55, 0.65, 0.87);
         else
            CreateLegend( (part == PION ? 0.5 : 0.6) - xshift, 0.55, 0.9, 0.8);
         legend->AddEntry(histSubtracted, "Data","ple");
         legend->AddEntry(systUncertainty, "Systematic uncertainty", "f");
         if(part != PROTON)
            legend->AddEntry(histMC, "GRANIITTI Res.+Cont." + granTag,"l");
         legend->AddEntry(histContinuum, "GRANIITTI Cont." + granTag,"l");
         if( hName == "DeltaPhi" && rpCon == nRpConfigurations){
            if(part != PROTON)
               legend->AddEntry(histMC2, "GRANIITTI Res.+Cont." + granTag2,"l");
            legend->AddEntry(histContinuum2, "GRANIITTI Cont." + granTag2,"l");
         }
         //legend->AddEntry(systUncertainty[bcg], "Background sys. uncert.", "f");
         legend->Draw("same");

         histSubtracted->Draw("AXIS same");
         WriteCanvas(hName + nameSuffix + canTag[i]); 
         canvas->Close();
      }
      /////////////////////////////////////  
      CreateCanvas(&canvas, hName + "_" + nameSuffix);      
      gPad->SetMargin(0.07,0.13,0.105,0.03);
      SetHistStyle(hMissPtVsX);
      hMissPtVsX->GetZaxis()->SetTitleOffset(1.4);
      hMissPtVsX->Draw("colz");
      SetPalletRange(hMissPtVsX);


      double xMin = hName == "DeltaPhi" ? 0.0 : massBins[part][2];
      double xMax = hName == "DeltaPhi" ? 180.0 : massBins[part][massNBins[part]-2];

      DrawExclusiveLine(xMin, exclusivityCut, xMax, exclusivityCut);

      CreateLine(xMin, nonExclFitRangeMin[part], xMax, nonExclFitRangeMin[part] );
      line->SetLineStyle(10);
      line->SetLineWidth(8);
      line->SetLineColor(2);
      line->Draw("same");
      CreateLine(xMin, nonExclFitRangeMax[part], xMax, nonExclFitRangeMax[part] );
      line->SetLineStyle(10);
      line->SetLineWidth(8);   
      line->SetLineColor(2);   
      line->Draw("same");

      WriteCanvas("pTMissVs" + hName + "_"+ nameSuffix);
      canvas->Close();
   }
}

void PlotManager::runAnaPlots()
{
   changeDir(kMAINANA, "AnaPlots");
   // load histogram from root file, edit and draw them
   auto drawHist = [&](TString histName, TString title, int ID = -1, bool setLogy = false, TString textInput = "")
   {
      TString dir = "AnaPlots";
      TString histFullName = dir + "/" + histName;

      bool twoDim = ID == -1 ? true : false;

      TH1 *hist = nullptr;
      inFile->GetObject(histFullName, hist);

      if (!hist) {
         cerr<<"Error: cannot loaded histogram "<<histName<<"in PlotMainAna::runAnaPlots()"<<endl; 
         return;
      }

      CreateCanvas(&canvas, histName);
      if(twoDim)
         gPad->SetMargin(0.07,0.13,0.105,0.03); // 2D Margin
      else if(ID == 8)
         gPad->SetMargin(0.07,0.01,0.155,0.01);
      else
         gPad->SetMargin(0.07,0.02,0.105,0.05);

      if(setLogy)
         canvas->SetLogy();
      SetHistStyle(hist);
      hist->SetTitle(title);

      TString opt = twoDim ? "colz" : "E1";
      if(twoDim)
         hist->GetZaxis()->SetTitleOffset(1.6);
      else
         hist->GetYaxis()->SetRangeUser(1, setLogy ? 10*hist->GetMaximum() : 1.2*hist->GetMaximum());

      if( ID >= 3 && ID <= 6) 
         hist->GetYaxis()->SetTitleOffset(1.5);
      else if( ID == 8){
         hist->SetLineColor(1);
         hist->SetMarkerColor(1);
         hist->LabelsOption("d");
         opt = "hist TEXT15";
      }

      hist->Draw(opt);

      DrawSTARInternal();
      DrawSystemDescription();

      if( textInput != ""){
         CreateText(0.75, 0.8, 0.88, 0.88);
         text->AddText( textInput );
         text->Draw("same");
      }

      if(twoDim)
         drawVertexEtaSpace();

      if( ID == 1)
         HighlightBin(hist,2);
      else if( ID == 2)
         HighlightBin(hist,3);
      else if( ID == 3)
         DrawLine(hist,vertexRange, true);
      else if( ID == 4)
         DrawLine(hist, minNHitsFit[NOMINAL]);
      else if( ID == 5)
         DrawLine(hist, minNHitsDEdx[NOMINAL]);
      else if( ID == 6)
         DrawLine(hist,maxDcaXY[NOMINAL]);
      else if( ID == 7)
         DrawLine(hist,maxDcaZ[NOMINAL], true);

      WriteCanvas(histName);
      canvas->Close();

   };

   for (int i = 0; i < nBranches; ++i){
      drawHist("hNRPTracksInFV"+mUtil->branchName(i), ";Number of RP tracks in FV;Number of events", 1, true, "Branch: "+ mUtil->branchName(i));
   }

   drawHist("hEtaZVertex", ";#eta;z_{vtx} [cm];Number of tracks", -1, false); //TH2D
   drawHist("hNTOFVerticies", ";Number of TOF vertices;Number of events",1,true);
   drawHist("hNTOFGoodTracks", ";Number of good TOF tracks in TOF vertex;Number of events",2,true);
   drawHist("hZVertex", ";z_{vtx} [cm];Number of events",3);
   drawHist("hNHitFit", ";N^{fit}_{hits};Number of tracks",4);
   drawHist("hNHitDedx", ";N^{dE/dx}_{hits};Number of tracks",5);
   drawHist("hDcaXY", ";DCA_{xy} [cm];Number of tracks",6);
   drawHist("hDCAZ", ";DCA_{z} [cm];Number of tracks",7,true);
   drawHist("hAnaFlow", ";;Number of events",8, true);

      TString histName = "hRPFV";
   CreateCanvas(&canvas, histName);
   canvas->cd();
   gPad->SetMargin(0.0, 0.0, 0.1, 0.01);
   canvas->Divide(2, 1, .0, .0, 0);

   TH2D* hist[nBranches];

   for (int iBr = 0; iBr < nBranches; ++iBr) 
   {
      
      TString histFullName = "AnaPlots/" + histName + mUtil->branchName(iBr);

      inFile->GetObject(histFullName, hist[iBr]);

      if (!hist[iBr]) {
         cerr<<"Error: cannot loaded histogram "<<histName<<"in PlotMainAna::runAnaPlots()"<<endl; 
         return;
      }

      if(iBr % 2 == 0) 
         continue;

      hist[iBr-1]->Add(hist[iBr]);

      bool Last = iBr == nBranches - 1;
      int side = iBr/2; 

      canvas->cd(side+1);
      SetHistStyle(hist[iBr-1]);
      hist[iBr-1]->SetTitle(";p_{x} [GeV];p_{y} [GeV];Number of tracks");
      hist[iBr-1]->GetYaxis()->SetTitleOffset(1.0);

      if ( Last ) 
         hist[iBr-1]->GetZaxis()->SetTitleOffset(1.3);
      else 
         hist[iBr-1]->GetZaxis()->SetTickLength(0);
      gPad->SetMargin(  side == East ? 0.12 : 0.0, Last ? 0.19 : 0.0, 0.12, 0.01);
      
      gPad->SetLogz();

      hist[iBr-1]->Draw("colz");
      DrawFiducial(side);

      CreateText( side == East ? 0.15 : 0.05, 0.9, side == East ? 0.25 : 0.15, 0.95);
      text->AddText( mUtil->sideName(side) );
      text->Draw("same");
   }
   
   CreateLegend(0.05, 0.52, 0.7, 0.58);
   legend->AddEntry(line, "Forward proton fiducial region", "f");
   legend->SetFillStyle(1001);
   legend->SetFillColor(kWhite);
   legend->Draw("same");
   
   WriteCanvas(histName);
   canvas->Close();
}



void PlotManager::correctMainPlotsForBinWidth()
{
   // correct bin values for bin width
   auto scaledHist1D = [&](TH1* hist)
   {
      for (int i = 1; i <= hist->GetNbinsX(); ++i){
         hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i) );
         hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i) );
      }
   };   

   auto scaledHist2D = [&](TH2* hist)
   {
      for (int i = 1; i <= hist->GetNbinsX(); ++i)
         for (int j = 1; j <= hist->GetNbinsY(); ++j){
            hist->SetBinContent(i,j, hist->GetBinContent(i,j) / hist->GetXaxis()->GetBinWidth(i) );
            hist->SetBinError(i,j, hist->GetBinError(i,j) / hist->GetXaxis()->GetBinWidth(i) );
         }
   }; 

   for (size_t iGroup = 0; iGroup < mainAnaHists.size(); ++iGroup) 
   {
      for (unsigned int i = 0; i <= nTotalTPCTOFSysStudies; ++i)
      {
         scaledHist1D(mainAnaHists[iGroup].hMain[i]);
         scaledHist2D(mainAnaHists[iGroup].hVsPt[i]);
      }    
   }
}//correctMainPlotsForBinWidth

void PlotManager::subtractBackgroundFromMainPlots()
{
   for (size_t iGroup = 0; iGroup < mainAnaHists.size(); ++iGroup) 
   {
      for (unsigned int i = 0; i <= nTotalTPCTOFSysStudies; ++i)
      {
         mainAnaHists[iGroup].hMainBcgSubtracted[i] = (TH1F*) mainAnaHists[iGroup].hMain[i]->Clone( Form("%s_subtracted",mainAnaHists[iGroup].hMain[i]->GetName()) );
         TH1D *hBkgdHistogram = (TH1D*)backgroundStudy(mainAnaHists[iGroup].hVsPt[i], mainAnaHists[iGroup].part, mainAnaHists[iGroup].rpCon);
         subtractBackground( mainAnaHists[iGroup].hMainBcgSubtracted[i], hBkgdHistogram );
         //if( mainAnaHists[iGroup].hName == "DeltaPhi")
         //   mainAnaHists[iGroup].hMainBcgSubtracted[i]->Scale(1000); // convert from nb to pb
      }    
   }
}//subtractBackgroundFromMainPlots

void PlotManager::integrateCrossSection(hGroup mainHists, TString hName)
{
   TH1F *hist = (TH1F*)mainHists.hMain[0]->Clone( hName + "_Integrated");

   string yTitle = hist->GetYaxis()->GetTitle();
   TString unit = " [" + yTitle.substr(yTitle.size() - 7, 2) + "]";
   
   int part = mainHists.part;
   int rpCon = mainHists.rpCon;

   hist->Scale( 1 - mBcgFraction[part][rpCon] );
   double statError = 0.0;  // For statistical uncertainty
   double integral = hist->IntegralAndError(1, hist->GetNbinsX(), statError, "width");

   double systErrorLow = integral*sqrt( pow(relLumUncertainty,2) + pow(relativEffSysLow[part],2) + pow(mBcgFracError[part][rpCon],2));
   double systErrorHigh = integral*sqrt( pow(relLumUncertainty,2) + pow(relativEffSysHigh[part],2) + pow(mBcgFracError[part][rpCon],2));

   cout<<"Integrated cross-section for "<< hist->GetName() <<" is "<<integral<<" +- "<<statError<<" + "<<systErrorHigh<<" - "<<systErrorLow<<unit<<endl;
}//integrateCrossSection