#include "PlotManager.h"


void PlotManager::PlotSysStudy()
{
    vector<int> colorSet = { 1, 2, 4, 6, 8};
    vector<int> markerStyle = { mainMarker, 24, 25, 26, 28};
/*
    auto GetAverageYValue = [&](TH1* hist) 
    {
        double sum = 0;
        double sumWeights = 0;

        // Loop over bins
        for (int i = 1; i <= hist->GetNbinsX(); i++) {
            sum += hist->GetBinContent(i) * hist->GetBinWidth(i);  // Weighted sum of the bin contents
            sumWeights += hist->GetBinWidth(i);        // Sum of the bin widths
        }
        
        if (sumWeights != 0)
            return sum / sumWeights;
        else
            return 1.0;  // Avoid division by zero
    };
*/

    auto GetAverageYValue = [&](TH1* hist) 
    {
        double weightedSum = 0;
        double sumWeights = 0;

        // Loop over bins
        for (int i = 1; i <= hist->GetNbinsX(); i++) {
            double content = hist->GetBinContent(i);
            double binWidth = hist->GetBinWidth(i); 
            double uncertainty = hist->GetBinError(i);

            if (uncertainty > 0) {
                double weight = binWidth / (uncertainty * uncertainty); // Weight by inverse variance
                weightedSum += content * weight;
                sumWeights += weight;
            }
        }

        if (sumWeights != 0)
            return weightedSum / sumWeights; // Return the weighted average
        else
            return 1.0; // Avoid division by zero
    };


    auto drawSysPlot = [&](TString canvasName, vector<TH1*> &hists,const vector<TString> &histLegends, bool plotAverage = true)
    {
        CreateCanvas(&canvas,canvasName);
        CreateLegend(0.08,0.85,0.85,0.95);
        legend->SetNColumns(hists.size());
        hists[0]->GetYaxis()->SetRangeUser(0.9,1.1);
        for (unsigned int i = 0; i < hists.size(); ++i) 
        {
            TH1F *hist = (TH1F*) hists[i]->Clone( canvasName + Form("tmp%i_",i));
            //hist->Draw("hist");
            //WriteCanvas(canvasName + Form("tmp%i_",i));
            // TEfficiency is not suitable since the hist are filled with weight
            hist->Divide(hists[0]);
            SetHistStyle(hist, colorSet[i], markerStyle[i]);

            if( i == 0)
                hist->Draw("hist");
            else
                hist->Draw("same hist");

            legend->AddEntry(hist, histLegends[i],"l");
            if(plotAverage && i!=0)
            { 
                CreateLine(hist->GetXaxis()->GetXmin(),GetAverageYValue(hist),hist->GetXaxis()->GetXmax(),GetAverageYValue(hist));
                line->SetLineWidth( lineWidth );
                line->SetLineStyle(10);
                line->SetLineColor(colorSet[i]);
                line->Draw("same"); 
            }
        }
        mCurrDir->cd();
        legend->Draw("same");
        WriteCanvas(canvasName);
        canvas->Close();
    };

    const vector<TString> histLegTPCApp = { TString("TPC averaged"), TString("TPC(#phi)"), TString("TPC(#eta,z)"), TString("TPC(#eta,#phi)") }; //TString("TPC(pT)"), TString("TPC(#eta,#phi)") };
    const vector<TString> histLegTOFApp = { TString("TOF averaged"), TString("TOF(#phi)"), TString("TOF(#eta,z)"), TString("TOF(#eta,#phi)") }; //TString("TPC(pT)"), TString("TPC(#eta,#phi)") };
    const vector<TString> histLegTPCnHits = { TString("Nominal"), TString("N^{fit}_{hits} loose"), TString("N^{fit}_{hits} tight"), TString("N^{dE/dx}_{hits} loose"), TString("N^{dE/dx}_{hits} tight") };
    const vector<TString> histLegPidSys = { TString("Nominal"), TString("Loose"), TString("Tight") };
    const vector<TString> histLegDcaSys = { TString("Nominal"), TString("DCA_{z} loose"), TString("DCA_{z} tight"), TString("DCA_{xy} loose"), TString("DCA_{xy} tight") };

    const vector<TString> histLegBcgSys = { TString("Nominal"), TString("Smaller range"), TString("Wider range") }; //, TString("Bin-by-bin fits") };
    vector<TH1*> hists; 

    for (size_t iGroup = 0; iGroup < mainAnaHists.size(); ++iGroup) 
    {
        int part = mainAnaHists[iGroup].part;
        int rpConfig = mainAnaHists[iGroup].rpCon;
        TString name = mainAnaHists[iGroup].hName;

        TString nameSuffix = mUtil->particleName(part);
        if(rpConfig < nRpConfigurations)
            nameSuffix += mUtil->rpConfigName(rpConfig);

        changeSubDir(nameSuffix);

        nameSuffix += mainAnaHists[iGroup].hName;
        // Prepare histograms
        // Get the main histogram, the "nominal one"
        TH1F *hSignal = (TH1F*) mainAnaHists[iGroup].hMainBcgSubtracted[0]->Clone( name + "InSysStudy_" + nameSuffix );
        hSignal->SetTitle(";" + mainAnaHists[iGroup].xLabel + ";Ratio to nominal");
      
        // TPCAppSysStudies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = 1; i <= nTPCAppSysStudies; ++i)
                hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_TPCAppSysStudy",hists,histLegTPCApp);

        // TOFAppSysStudies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = nTPCAppSysStudies + 1; i <= nTPCAppSysStudies + nTOFAppSysStudies; ++i)
                hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_TOFAppSysStudy",hists,histLegTOFApp);

        // nTPCnHitsStudies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = nTPCAppSysStudies + nTOFAppSysStudies + 1; i <= nTPCAppSysStudies + nTOFAppSysStudies + nTPCnHitsStudies; ++i)
            hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_TPCnHitsStudy",hists,histLegTPCnHits);

        // PID Sys studies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = nTPCAppSysStudies + nTOFAppSysStudies + nTPCnHitsStudies + 1; i <= nTPCAppSysStudies + nTOFAppSysStudies + nTPCnHitsStudies + nPidVariation; ++i)
            hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_PidStudy",hists,histLegPidSys);

        // DCA Sys studies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = nTPCAppSysStudies + nTOFAppSysStudies + nTPCnHitsStudies + nPidVariation + 1; i <= nTotalTPCTOFSysStudies; ++i)
            hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_DcaStudy",hists,histLegDcaSys);

        // background part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = 1; i <= histLegBcgSys.size() - 1; ++i)
        {
            TH1D *hBkgdHist;
            TH1F *hSigComparison = (TH1F*)mainAnaHists[iGroup].hMain[0]->Clone( name + Form("InSysStudy%i_",i) + nameSuffix );
            if( i < 3)
                hBkgdHist = (TH1D*)backgroundStudy(mainAnaHists[iGroup].hVsPt[0], part, rpConfig, "SysBcg", i);
            else // bin-by-method
                hBkgdHist = (TH1D*)bkgdHistogramComparison(mainAnaHists[iGroup].hVsPt[0], part, rpConfig, false);


            subtractBackground( hSigComparison, hBkgdHist );
            hists.push_back(hSigComparison);

            TH1F* tmp = (TH1F*)hSigComparison->Clone("tmp"); 
            tmp->Divide(hSignal);
            if(i == 1){
                mainAnaHists[iGroup].hBcgSysStudy = (TH1F*)tmp->Clone(name + "hBcgSysStudy" + nameSuffix); 
            }else if( i == 2){
                mainAnaHists[iGroup].hBcgSysStudy->Add(tmp);
                mainAnaHists[iGroup].hBcgSysStudy->Scale(0.5);
            }
        }
        drawSysPlot(nameSuffix + "_backgroundSysStudy",hists,histLegBcgSys, false);
    }
}//PlotSysStudy

void PlotManager::DrawSysUncertainty(hGroup mainHists, bool drawCanvas)
{
    int part = mainHists.part;
    int rpCon = mainHists.rpCon;

    TString nameSuffix = mUtil->particleName(part);
    if(rpCon < nRpConfigurations)
        nameSuffix += mUtil->rpConfigName(rpCon);

    TH1* hist = mainHists.hMainBcgSubtracted[0];

    // Prepare arrays for TGraphAsymmErrors
    int nBins = hist->GetNbinsX();
    double x[nBins], y[nBins], exl[nBins], exh[nBins], zeros[nBins];
    double eBcg[nBins], eLum[nBins], eCorrLow[nBins], eCorrHigh[nBins], eTotLow[nBins], eTotHigh[nBins];

    const std::vector<int> binsToDrawSystUnc[nParticles] = { {9, 13, 20, 30, 32, 43, 50, 60, 70}, {3, 9, 13, 20, 24}, {3, 6, 8}  };
    const std::vector<int> binsToDrawSystUncForDeltaPhi = { 1, nBins-1};

    int nSysBins = mainHists.hName == "DeltaPhi" ? binsToDrawSystUncForDeltaPhi.size() : binsToDrawSystUnc[part].size();
    const std::vector<int>* vec = mainHists.hName == "DeltaPhi" ? &binsToDrawSystUncForDeltaPhi : &(binsToDrawSystUnc[part]);

    // Fill the arrays with bin content and systematic uncertainties
    //for (int i = 0; i < nSysBins; ++i) {
    //    int binId = (*vec)[i];
    for (int i = 0; i < nBins; ++i) {
        int binId = i + 1;
        x[i] = hist->GetBinCenter(binId);
        y[i] = hist->GetBinContent(binId);
        exl[i] = hist->GetBinWidth(binId)*0.5;
        exh[i] = hist->GetBinWidth(binId)*0.5;
        double sysBcgErr = mainHists.hBcgSysStudy->GetBinContent(binId)/100;
        // add background normalisation uncertainity
        sysBcgErr = sqrt( pow(sysBcgErr,2) + pow(mBcgFracError[part][rpCon],2) );  
        zeros[i] = 0.0;
        eBcg[i] = y[i] * sysBcgErr;
        eLum[i] = y[i] * relLumUncertainty;
        eCorrLow[i] = y[i] * relativEffSysLow[part];
        eCorrHigh[i] = y[i] * relativEffSysHigh[part];
        eTotLow[i] = y[i] * sqrt(relLumUncertainty*relLumUncertainty + sysBcgErr*sysBcgErr + relativEffSysLow[part]*relativEffSysLow[part]);
        eTotHigh[i] = y[i] * sqrt(relLumUncertainty*relLumUncertainty + sysBcgErr*sysBcgErr + relativEffSysHigh[part]*relativEffSysHigh[part]);
    }
    // Create a TGraphAsymmErrors for systematic uncertainties
    if(systUncertainty)
        delete systUncertainty;

    //systUncertainty = new TGraphAsymmErrors(nSysBins, x, y, exl, exh, eTotLow, eTotHigh);
    systUncertainty = new TGraphAsymmErrors(nBins, x, y, exl, exh, eTotLow, eTotHigh);
    systUncertainty->SetFillColor(18);
    systUncertainty->SetLineColor(16);
    systUncertainty->SetLineWidth(2);
    // Draw the histogram with statistical error bars and systematic uncertainties
    systUncertainty->Draw("E5 same"); // Draw systematic uncertainties as a filled area

    if( !drawCanvas )
        return; 

    enum SYSTUNC { bcg=0, lum, effCorr, total, nSysUnct }; 
    TGraphAsymmErrors* sysUnct[nSysUnct];


    // Fill the arrays with bin content and systematic uncertainties
    for (int i = 0; i < nBins; ++i) {
        int binId = i + 1;
        x[i] = hist->GetBinCenter(binId);
        y[i] = hist->GetBinContent(binId);
        exl[i] = hist->GetBinWidth(binId)*0.5;
        exh[i] = hist->GetBinWidth(binId)*0.5;
        double sysBcgErr = mainHists.hBcgSysStudy->GetBinContent(binId)/100;
        // add background normalisation uncertainity
        sysBcgErr = sqrt( pow(sysBcgErr,2) + pow(mBcgFracError[part][rpCon],2) );  
        zeros[i] = 0.0;
        eBcg[i] = y[i] * sysBcgErr;
        eLum[i] = y[i] * relLumUncertainty;
        eCorrLow[i] = y[i] * relativEffSysLow[part];
        eCorrHigh[i] = y[i] * relativEffSysHigh[part];
        eTotLow[i] = y[i] * sqrt(relLumUncertainty*relLumUncertainty + sysBcgErr*sysBcgErr + relativEffSysLow[part]*relativEffSysLow[part]);
        eTotHigh[i] = y[i] * sqrt(relLumUncertainty*relLumUncertainty + sysBcgErr*sysBcgErr + relativEffSysHigh[part]*relativEffSysHigh[part]);
    }

    TCanvas *cSysStudy; 
    CreateCanvas(&cSysStudy,"sysStudy_"+nameSuffix);
    sysUnct[total] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eTotLow, eTotHigh);
    SetGraphStyle(sysUnct[total]);
    sysUnct[total]->SetTitle(";m(" + mUtil->pairLabel(part) + ") [GeV];d#sigma/dm(" + mUtil->pairLabel(part) + ") [nb/GeV]");
    sysUnct[total]->SetFillColor(18);
    sysUnct[total]->SetLineColor(16);
    sysUnct[total]->SetLineWidth(2);
    sysUnct[total]->GetXaxis()->SetRangeUser(x[0] - exl[0], x[nBins - 1] + exh[nBins - 1]);
    sysUnct[total]->GetYaxis()->SetRangeUser(sysUnct[total]->GetYaxis()->GetXmin() , 1.4*sysUnct[total]->GetYaxis()->GetXmax() );
    sysUnct[total]->Draw("");
    sysUnct[total]->Draw("E5 same");

    sysUnct[lum] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eLum, eLum);
    sysUnct[lum]->SetFillColor(3);
    sysUnct[lum]->SetFillStyle(3008);
    sysUnct[lum]->Draw("E5 same");

    sysUnct[effCorr] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eCorrLow, eCorrHigh);
    sysUnct[effCorr]->SetFillColor(4);
    sysUnct[effCorr]->SetFillStyle(3009);
    sysUnct[effCorr]->Draw("E5 same");  

    sysUnct[bcg] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eBcg, eBcg);
    sysUnct[bcg]->SetFillColor(2);
    sysUnct[bcg]->SetFillStyle(1001); 
    sysUnct[bcg]->Draw("E5 same"); 

    CreateLegend(0.64, 0.6, 0.9, 0.95);
    legend->AddEntry(sysUnct[total], "Systematic uncertainty", "f");
    legend->AddEntry(sysUnct[bcg], "Background sys. uncert.", "f");
    legend->AddEntry(sysUnct[lum], "Luminosity sys. uncert.", "f");
    legend->AddEntry(sysUnct[effCorr], "Eff. corr. sys. uncert.", "f");
    legend->Draw("same");

    WriteCanvas("sysStudy_"+nameSuffix,cSysStudy);
    cSysStudy->Close();
}//DrawSysUncertainty
     

void PlotManager::runSysStudy()
{
    unsigned int id = nTPCAppSysStudies + nTOFAppSysStudies;
    for (unsigned int i = 1; i <= nTPCnHitsStudies + nPidVariation + nDcaVariation; ++i)
    {
        TString treeName = nameOfTree[kMAINANA] + "_";
        if( i <= nTPCnHitsStudies )
            treeName += mUtil->nHitVaryName(i -1);
        else if(i <= nTPCnHitsStudies + nPidVariation)
            treeName += mUtil->pidVaryName(i -1 - nTPCnHitsStudies);
        else
            treeName += mUtil->dcaVaryName(i -1 - nTPCnHitsStudies - nPidVariation);

        TTree *tree = dynamic_cast<TTree*>( inFile->Get( treeName ) );
        if (!tree){ 
            cerr<<"Error: cannot loaded tree "<<treeName<<"in PlotMainAna::runSysStudy()"<<endl; 
            return;
        }
        mCurrentTree = new RecTree(tree, treeBits[kMAINANA] );
        if (!mCurrentTree){ 
            cerr<<"Error: cannot loaded recTree in PlotMainAna::runSysStudy"<<endl; 
            return;
        }
        for(Long64_t iev=0; iev<tree->GetEntries(); ++iev)
        { //get the event
            tree->GetEntry(iev); 
            if( i <= nTPCnHitsStudies ){
                std::pair<VARIATION, VARIATION> nHitsVar = mUtil->varyNHits(i-1);  
                loopThroughTree(NOMINAL, NOMINAL, nHitsVar.first, nHitsVar.second, NOMINAL, NOMINAL, id + i);
            }else if( i <= nTPCnHitsStudies + nPidVariation ){
                loopThroughTree(NOMINAL, NOMINAL, NOMINAL, NOMINAL, i - nTPCnHitsStudies, NOMINAL, id + i);
            }else{
                loopThroughTree(NOMINAL, NOMINAL, NOMINAL, NOMINAL, NOMINAL, i - nTPCnHitsStudies - nPidVariation, id + i);
            }
        } 
    }
}//runSysStudy


