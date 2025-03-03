#include "PlotManager.h"
#include "Ana.h"


//_____________________________________________________________________________
void PlotManager::runRPMCPlots(){
   if( runELASTICANA )
   {
      changeDir(kRPMCANA, "rpMcComparison");
      if( DEBUG )
         cerr<<"PlotManager::runRPMCPlots() going to run runRPMCComparison()"<<endl;
      runRPMCComparison();
   }
   if( DEBUG )
      cerr<<"PlotManager::runRPMCPlots() going to run runRPEff()"<<endl;
   runRPEff();

   changeDir(kRPMCANA, "RPEffStudy");
   runRPEffStudy();

}//runRPMCPlots

//_____________________________________________________________________________
void PlotManager::runRPMCComparison(){
   vector<TString> anaName = { mUtil->dataSetName(MC), mUtil->dataSetName(MCZB), studyName[kELASTICANA] };
   vector<int> colorSet = { mainColor, bckgColor, 1};
   vector<int> markerStyle = { mainMarker, bckgMarker, 20};
   TH1D *hist[anaName.size()];


   vector<TString> histList;
   vector<unsigned int> flags; 
   vector<TString> labelList;
   vector<unsigned int> nBins;
   vector<double> binLow; 
   vector<double> binMax; 

   AddHist("hDeltaTheta","#Delta#theta [rad]",0,-0.001,0.001,0,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("hDCutR","DCut [cm?]",0,0.0,0.012,0,histList,labelList,nBins,binLow,binMax,flags);

   AddHist("hClusterLength","Cluster length",0,0.0,25.0,1,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("hClusterEnergy","E [ADC]",0,0.0,80.0,1,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("hNClusterPerRP","N_{Clusters} per RP",0,0.0,25.0,1,histList,labelList,nBins,binLow,binMax,flags);

   for (int iCoor = 0; iCoor < nCoordinates-1; ++iCoor)
      AddHist("hDeltaFitProj_"+mUtil->coordinateName(iCoor),"#Delta d_{" + mUtil->coordinateName(iCoor) +"} [mm]",0,-0.8,0.8,1,histList,labelList,nBins,binLow,binMax,flags);

   for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
      AddHist("hClusterPerPlain_"+mUtil->planeName(iPlane),"N_{Clusters}",0,0.0,10.0,1,histList,labelList,nBins,binLow,binMax,flags);

   AddHist("hDeltaProj","#Delta d [mm]",0,-0.6,0.6,2,histList,labelList,nBins,binLow,binMax,flags);

   TString histNameInFile; 
   TString yTitle = "Normalized counts";
   setTextSize(30,35);

   vector<unsigned int> xDiv = {1, 4, 4};
   vector<unsigned int> yDiv = {1, 2, 2};
   vector<unsigned int> nDiv;
   for( unsigned int i = 0; i < xDiv.size(); ++i)
      nDiv.push_back(xDiv[i]*yDiv[i]);

   vector<TString> flagName[xDiv.size()];
   int RPordered[] = { E2U, E1U, W1U, W2U, E2D, E1D, W1D, W2D};
   for (int i = 0; i < nRomanPots; ++i)
      flagName[1].push_back( mUtil->rpName( RPordered[i] ) );

   int brOrder[] = { EU, WU, ED, WD};
   for (int iCoor = 0; iCoor < nCoordinates-1; ++iCoor)
      for (int iBr = 0; iBr < nBranches; ++iBr)
         flagName[2].push_back( mUtil->coordinateName(iCoor) + "_" + mUtil->branchName(brOrder[iBr]) );

   for (unsigned int iHist = 0; iHist < histList.size(); ++iHist)
   {
      CreateCanvas(&canvas,histList[iHist]);
      if(nDiv[flags[iHist]] > 1)
      {
         gPad->SetMargin(0.15,0.02,0.12,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
         canvas->Divide(xDiv[flags[iHist]],yDiv[flags[iHist]],.0,.0,0);
      }
      for (unsigned int i = 0; i < nDiv[flags[iHist]]; ++i)
      {
         if( i == 0){
            if( histList[iHist] == TString("hClusterEnergy") )
               CreateLegend(0.7,0.2,0.97,0.4);
            else
               CreateLegend(0.7,0.7,0.97,0.9);
         }
         if(nDiv[flags[iHist]] > 1)
            canvas->cd(i+1);
         gPad->SetLogy();
         for (unsigned int iSet = 0; iSet < anaName.size(); ++iSet)
         {
            histNameInFile = anaName[iSet] + "/" + histList[iHist];
            if(flags[iHist] > 0)
               histNameInFile += "_" + flagName[flags[iHist]][i];

            hist[iSet] = (TH1D*)inFile->Get( histNameInFile );
            if(!hist[iSet]){
               cerr<<"Error: PlotManager::runRPMCComparison() cannot loaded hist "<<histNameInFile<<endl; 
               continue;
            }
            hist[iSet] = (TH1D*)hist[iSet]->Clone("hist" + anaName[iSet]);
            hist[iSet]->GetXaxis()->SetRangeUser( binLow[iHist], binMax[iHist]);
            hist[iSet]->Scale(1.0/double(hist[iSet]->Integral(hist[iSet]->FindFixBin(binLow[iHist]),hist[iSet]->FindFixBin(binMax[iHist]))));
            SetHistStyle(hist[iSet], colorSet[iSet], markerStyle[iSet]);
            hist[iSet]->SetTitle(";" + labelList[iHist] +";"+yTitle); 
            if( iSet == 0)
               hist[iSet]->Draw("hist");
            else
               hist[iSet]->Draw("same hist");
         
            if( i == 0)
               legend->AddEntry(hist[iSet], mUtil->dataSetName(iSet),"ple");
         }

         if(nDiv[flags[iHist]] > 1)
         {
            CreateText(0.4,0.85,0.5,0.95);
            if( flags[iHist] == 1)
               text->AddText(mUtil->rpName( RPordered[i] ));
            else if( flags[iHist] == 2 ) 
               text->AddText( flagName[2][i] );
            text -> Draw("same");
         }
      }
      legend->Draw("same");
      WriteCanvas(histList[iHist]);
   }

//////////////////////////////////////////////////////////////////
   // plot TH2D
   TH2D *hist2D[2];
   // Plot pX, pY accatpance
   for (int iSet = 0; iSet < nDataSets; ++iSet)
   {
      CreateCanvas(&canvas,"pXpYAcceptance"+mUtil->dataSetName(iSet));
      canvas->Divide(2,1);
      for (int iBr = 0; iBr < nBranches; ++iBr)
      {         
         canvas->cd((iBr/2)+1);
         hist2D[iBr%2] = (TH2D*)inFile->Get(anaName[iSet]+"/hAcceptancePxPy_"+mUtil->branchName(iBr));

         if( iBr%2 == 0)
            continue;

         hist2D[0]->Add(hist2D[1]);
         hist2D[0]->Add(hist2D[1]);

         SetGPad();
         hist2D[0]->SetTitle(mUtil->sideName(iBr/2) +";p_{x} [GeV];p_{y} [GeV];Events");
         SetHistStyle(hist2D[0]);
         hist2D[0]->Draw("colz");
         DrawFiducial(iBr/2);
      }
      WriteCanvas("pXpYAcceptance"+mUtil->dataSetName(iSet));
   }
   // Plot x,y accatpance
   /*
   for (int iSet = 0; iSet < nDataSets; ++iSet)
   {
      CreateCanvas(&canvas,"xyAcceptance"+mUtil->dataSetName(iSet));
      canvas->Divide(4,1);
      for (int iSt = 0; iSt < nStations; ++iSt)
      {
         for (int iRp = 0; iRp < nRpPerStation; ++iRp)
         {
            canvas->cd(iSt+1);
            hist2D[iRp] = (TH2D*)inFile->Get(anaName[iSet]+"/hAcceptanceXY_"+mUtil->rpName(2*iSt + iRp));

            if( iRp == 0)
               continue;

            hist2D[0]->Add(hist2D[1]);
            hist2D[0]->Add(hist2D[1]);

            SetGPad();
            hist2D[0]->SetTitle(mUtil->stationName(iSt) +";x [mm];y [mm];Events");
            SetHistStyle(hist2D[0]);
            hist2D[0]->Draw("colz");
         }
      }
      WriteCanvas("xyAcceptance"+mUtil->dataSetName(iSet));
   }
   */
   setTextSize();
}//runRPMCComparison

void PlotManager::runRPEff()
{
   setTextSize(30,35);


   auto plotAverageEff = [&](TH2D* hEff, TString histName, int br, bool elasticFV) 
   {

      TH1D *hAverageEff = new TH1D(histName + "_AverageEff",";<Eff>;",100,0.7,1.0);
      //TH2D *hAverageEff = new TH2D(histName + "_AverageEff",";<Eff>;",55,-0.55,0.55,120,-1.2,1.2);

      TCanvas *canvasEff;
      CreateCanvas(&canvasEff,histName + "_AverageEff"); 
      SetHistStyle(hAverageEff);
      gStyle->SetOptStat(111111);
      hAverageEff->SetStats(true);

      double eff;
      Ana anaObj(nullptr);
      // Loop through the bins and get the bin content
      for (int i = 1; i <= hEff->GetNbinsX(); ++i) { 
         double px = hEff->GetXaxis()->GetBinCenter(i);
         for (int j = 1; j <= hEff->GetNbinsY(); ++j) {
            double py = hEff->GetYaxis()->GetBinCenter(j);
            if( !(elasticFV ? anaObj.RPInElFidRange( px, py) : anaObj.RPInFidRange( px, py, br)) )
               continue;
            eff = hEff->GetBinContent(i, j);
            if( eff <= 0)
               continue;
            //hAverageEff->SetBinContent(i, j, eff);
            hAverageEff->Fill(eff);
         }
      }

      hAverageEff->Draw("colz");
      WriteCanvas(histName + "_AverageEff", canvasEff);
   };



   auto createAndSetupCanvas = [&](TCanvas** canvas, TString name, int divX, int divY, float left, float right, float bottom, float top) 
   {
      CreateCanvas(canvas, name);
      (*canvas)->cd();
      gPad->SetMargin(left, right, bottom, top);
      (*canvas)->Divide(divX, divY, .0, .0, 0);
   };

   auto loadHistograms = [&](TString prefix) 
   {
      return std::vector<TH2D*>{
         (TH2D*)inFile->Get(prefix + "_Total"),
         (TH2D*)inFile->Get(prefix + "_Passed"),
         (TH2D*)inFile->Get(prefix + "_PassedBefore")
      };
   };

   // Calculate efficiency before and after the distance cut
   auto computeEfficiency = [&](TH2D* total, TH2D* passed) 
   {
      auto efficiency = (TH2D*)passed->Clone();
      efficiency->Divide(total);
      return efficiency;
   };

   auto drawHistograms = [&](TCanvas* canvas, TH2D* hist, int totDiv, int div, TString textInput, const std::string& title, bool elasticFV) 
   {
      bool Last = div == totDiv - 1;
      div = div/2 + 1; 

      canvas->cd(div);
      hist->SetTitle(title.c_str());
      SetHistStyle(hist);
      hist->SetMinimum(0.825);
      //hist->GetXaxis()->SetRangeUser(-0.55, 0.65);
      if ( Last ) hist->GetZaxis()->SetTitleOffset(1.3);
      else hist->GetZaxis()->SetTickLength(0);
      gPad->SetMargin(  div == 1 ? 0.1 : 0.0, Last ? 0.16 : 0.0, 0.08, 0.01);
      hist->Draw("colz");
      if( totDiv == nBranches)
         DrawFiducial(div-1, elasticFV);

      CreateText(0.4, 0.45, 0.5, 0.55);
      text->AddText( textInput );
      text->Draw("same");
   };

   auto processEfficiency = [&](TString histName, bool isMC, bool doAvrgEff, bool elasticFV, const std::string& title, const int division, float left, float right, float bottom, float top) 
   {   
      int divX, divY;
      if( division == nRomanPots){
         divX = 4; divY = 1;
      }else if( division == nBranches){
         divX = 2; divY = 1;
      }else{
         return;
      }

      if( isMC ) histName += "MC_";
      const int nCanvases = isMC ? 2 : 1;

      for (int iSet = 0; iSet < (isMC ? DATA : nDataSets); ++iSet) 
      {
         TCanvas* effCanvas[nCanvases];
         TString canvasSuffix[] = { TString(""), TString("Before") }; 

         for (int i = 0; i < nCanvases; ++i)
            createAndSetupCanvas( &effCanvas[i], histName + canvasSuffix[i]+ mUtil->dataSetName(iSet), divX, divY, left, right, bottom, top);

         TH2D* hEff[division][nCanvases];
         for (int iDiv = 0; iDiv < division; ++iDiv) 
         {
            TString divName = division == nRomanPots ? mUtil->rpName(iDiv) : mUtil->branchName(iDiv);
            TString prefix = isMC ? "RPMC_" : ""; 
            prefix += iSet == DATA ? studyName[kELASTICANA] : mUtil->dataSetName(iSet);
            prefix += "/"+histName + divName; 

            auto histograms = loadHistograms(prefix);
            for (int i = 0; i < nCanvases; ++i)
               hEff[iDiv][i] = computeEfficiency(histograms[0], histograms[i+1]);

            if( doAvrgEff )
              plotAverageEff( hEff[iDiv][0], prefix, iDiv, elasticFV);

            if(iDiv % 2 == 0) continue;

            for (int i = 0; i < nCanvases; ++i)
               hEff[iDiv-1][i]->Add(hEff[iDiv][i]);


            TString textInput = division == nRomanPots ? mUtil->branchName(iDiv/2) : mUtil->sideName(iDiv/2);
            for (int i = 0; i < nCanvases; ++i){
               drawHistograms(effCanvas[i], hEff[iDiv-1][i], division, iDiv, textInput, title, elasticFV);
               if(i == 0 && histName.Contains("hEffPxPy") && isMC)
                  plotGradientAndFindStableRegions(hEff[iDiv-1][i], prefix, iDiv/2);
            }
         }
         for (int i = 0; i < nCanvases; ++i)
            WriteCanvas(histName + canvasSuffix[i] + mUtil->dataSetName(iSet), effCanvas[i]);
      }
   };
 
   changeDir(kRPMCANA, "rpEffMC");
   processEfficiency("hEffPxPy_", true, true,  false,";p_{x} [GeV];p_{y} [GeV]; Eff / (10 MeV #times 10 MeV)", nBranches, 0.0, 0.0, 0.0, 0.0);
   processEfficiency("hEffXY_", true, false,  false,";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
   processEfficiency("hEffXYOffSub_", true, false,  false,";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
   processEfficiency("hEffRpCheckedPxPy_", true, true,  false,";p_{x} [GeV];p_{y} [GeV]; Eff / (10 MeV #times 10 MeV)", nBranches, 0.0, 0.0, 0.0, 0.0);
   processEfficiency("hEffRpCheckedXY_", true, false,  false,";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
   processEfficiency("hEffRpCheckedXYOffSub_", true, false, false,";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
   if( runELASTICANA )
   {
      changeDir(kRPMCANA, "rpEff");
      processEfficiency("hEffPxPy_", false, true, true,";p_{x} [GeV];p_{y} [GeV]; Eff / (10 MeV #times 10 MeV)", nBranches, 0.0, 0.0, 0.0, 0.0);
      processEfficiency("hEffXY_", false, true,  false,";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
      //processEfficiency("hEffXYOffSub_MC_", false, ";x [cm];y [cm]; Eff", nRomanPots, 0.15, 0.02, 0.1, 0.0);
   }
   setTextSize();

}//runRPEff


void PlotManager::plotGradientAndFindStableRegions(TH2D *hist, TString hName, int side) {
    if (!hist) {
        std::cerr << "Error: PlotManager::plotGradientAndFindStableRegions Unable to retrieve efficiency histogram"<< std::endl;
        return;
    }

    // Calculate gradients
    TH2D *gradX = (TH2D*)hist->Clone("gradX");
    TH2D *gradY = (TH2D*)hist->Clone("gradY");

    gradX->Reset();
    gradY->Reset();

    for (int i = 2; i < hist->GetNbinsX(); ++i) {
        for (int j = 2; j < hist->GetNbinsY(); ++j) {
            double dx = hist->GetBinContent(i + 1, j) - hist->GetBinContent(i - 1, j);
            double dy = hist->GetBinContent(i, j + 1) - hist->GetBinContent(i, j - 1);

            gradX->SetBinContent(i, j, dx);
            gradY->SetBinContent(i, j, dy);
        }
    }

    // Create a canvas
    TCanvas *canvasGrad;
    CreateCanvas(&canvasGrad,hName + "_gradients"); //Form("%s_gradients",hName));
    //canvasGrad->Divide(2, 1);
    gradX->Add(gradY);
    SetHistStyle(gradX);
    // Plot gradient along X axis
    canvasGrad->cd(1);
    gradX->Draw("colz");
    DrawFiducial(side);
    WriteCanvas(hName + "_gradients", canvasGrad);
    return;
/*    

    // Plot gradient along Y axis
    canvasGrad->cd(2);
    SetHistStyle(gradY);
    gradY->Draw("colz");

    // Identify stable regions
    double threshold = 0.1; // Adjust as needed
    std::vector<std::pair<int, int>> stableRegions;

    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            if (std::abs(gradX->GetBinContent(i, j)) < threshold && std::abs(gradY->GetBinContent(i, j)) < threshold) {
                stableRegions.push_back(std::make_pair(i, j));
            }
        }
    }

    // Mark stable regions on the plot
    for (const auto& region : stableRegions) {
        double x = hist->GetXaxis()->GetBinCenter(region.first);
        double y = hist->GetYaxis()->GetBinCenter(region.second);
        TMarker marker(x, y, 20);
        marker.SetMarkerSize(0.5);
        marker.SetMarkerColor(kRed);
        marker.Draw();
    }
*/
}//plotGradientAndFindStableRegions

void PlotManager::runRPEffStudy()
{
   TH1D *hist;
   for (int iBr = 0; iBr < nBranches; ++iBr)
   {
      
      TString histNameInFile = "hAverageEffPerBr_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(MCZB);
      TString dirNameInFile = "RPMC_"+mUtil->dataSetName(MCZB)+"_EffStudy/";
      hist = (TH1D*)inFile->Get( dirNameInFile  + histNameInFile );
      if(!hist){
         cout<<"Warning in PlotManager::runRPEffStudy(): the hist "<<dirNameInFile  + histNameInFile<<" was not loaded..."<<endl;
         continue;
      }

      CreateCanvas(&canvas,histNameInFile);
      gPad->SetLogy();
      hist->GetXaxis()->SetRangeUser( 0.7, 1.0);
      SetHistStyle(hist);
      hist->SetTitle(";<Eff>;Runs"); 
      hist->Draw("hist");

      CreateText(0.1, 0.85, 0.2, 0.95);
      text->AddText( Form("#mu = %0.3f",hist->GetMean()) );
      text->AddText( Form("#sigma = %0.3f",hist->GetRMS()) );
      text->Draw("same");

      WriteCanvas(histNameInFile);
   }

   return;
   //Study average of an everage eff in single bin in fiducial volume

   const int nStat = 1; // 1 for hAverageEff and 2 for hAverageEff and hAverageEffRpChecked
   TString statName[] = { TString("hAverageEff"), TString("hAverageEffRpChecked") };

   TH1D *hAverageEffMean[nStat][nBranches][DATA], *hAverageEffSigma[nStat][nBranches][DATA];
   TH2D *hAverageEffMeanPerLum[nStat][nBranches][DATA];
   for (int iStat = 0; iStat < nStat; ++iStat){
      for (int iSet = 0; iSet < DATA; ++iSet){
         for (int iBr = 0; iBr < nBranches; ++iBr){
            hAverageEffMean[iStat][iBr][iSet] = new TH1D(statName[iStat] + "Mean_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet),";<Eff>;",200,0.7,1.2);
            hAverageEffSigma[iStat][iBr][iSet] = new TH1D(statName[iStat]  + "Sigma_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet),";<Eff>;",200,0.0,0.2);

            hAverageEffMeanPerLum[iStat][iBr][iSet] = new TH2D(statName[iStat] + "MeanPerLum_" + mUtil->branchName(iBr) + "_" + mUtil->dataSetName(iSet),";L [#mub^{-1}s^{-1}];<Eff>",100,60,160,200,0.7,1.2);
         }
      }
   }
   // Convert branch name to integer 
   auto brInteger = [&](TString brName) 
   {
      int branch = EU;
      if( brName == "ED")
         branch = ED;
      if( brName == "WU")
         branch = WU;
      if( brName == "WD")
         branch = WD;        

      return branch;
   };

   int runNumber;
   TKey *key;
   TClass *cl;
   TH1 *h;
   TString histName, brName;
   int stat;
   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      inFile->cd("RPMC_"+mUtil->dataSetName(iSet)+"_EffStudy");

      TIter histkeyList(gDirectory->GetListOfKeys());
      while ((key = (TKey*)histkeyList()))
      {
         cl = gROOT->GetClass(key->GetClassName());
         if (!cl->InheritsFrom("TH1")) 
            continue;
         h = (TH2*)key->ReadObj();
         histName = h->GetName(); 
         if( !histName.Contains("hAverageEff") )
            continue;

         stat = histName.Contains("hAverageEffRpChecked_") ? 1 : 0;
         if( stat == nStat)
            continue;
         // Remove the part of the string up to and including the underscore
         histName.Remove(0, histName.Index("_") + 1); // Remove hAverageEff_ / hAverageEffRpChecked_
         brName = histName(0, histName.Index("_")); // save branch name
         histName.Remove(0,histName.Index("_") + 1); // remove branch name
         histName.Remove(0,histName.Index("_") + 1); // remove data Set Name
         runNumber = histName.Atoi(); // convert the rest to integer
         if( runNumber == 999)
            continue;

         hAverageEffMean[stat][ brInteger(brName) ][iSet]->Fill(h->GetMean());
         hAverageEffSigma[stat][ brInteger(brName) ][iSet]->Fill(h->GetRMS());

         hAverageEffMeanPerLum[stat][ brInteger(brName) ][iSet]->Fill( mInstLumiPerRun[runNumber], h->GetMean());
      }
   }
   mDir->cd();
   for (int iStat = 0; iStat < nStat; ++iStat){
      for (int iSet = 0; iSet < DATA; ++iSet){
         for (int iBr = 0; iBr < nBranches; ++iBr){
            hAverageEffMean[iStat][iBr][iSet]->Write();
            hAverageEffSigma[iStat][iBr][iSet]->Write();

            hAverageEffMeanPerLum[iStat][iBr][iSet]->Write();
         }
      }
   }
}
