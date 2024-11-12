#include "PlotManager.h"


void PlotManager::PlotSysStudy()
{
    vector<int> colorSet = { 1, 2, 4, 6, 8};
    vector<int> markerStyle = { mainMarker, 24, 25, 26, 28};

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


    auto drawSysPlot = [&](TString canvasName, vector<TH1*> &hists,const vector<TString> &histLegends, bool plotAverage = true)
    {
        CreateCanvas(&canvas,canvasName);
        CreateLegend(0.08,0.85,0.85,0.95);
        legend->SetNColumns(hists.size());
        hists[0]->GetYaxis()->SetRangeUser(0.9,1.1);
        for (unsigned int i = 0; i < hists.size(); ++i) 
        {
            TH1F *hist = (TH1F*) hists[i]->Clone( canvasName + Form("tmp%i_",i));
            //cerr<<"Going to plot: "<<canvasName<<" by dividing: "<<hists[i]->GetName()<<" and "<<hists[0]->GetName()<<endl;
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
    const vector<TString> histLegTPCnHits = { TString("Nominal"), TString("N^{dE/dx}_{hits} loose"), TString("N^{fit}_{hits} loose"), TString("N^{dE/dx}_{hits} tight"), TString("N^{fit}_{hits} tight") };
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
            if( i != 3)
                hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_TPCAppSysStudy",hists,histLegTPCApp);

        // nTPCnHitsStudies part
        hists.clear();
        hists.push_back(hSignal);
        for (unsigned int i = nTPCAppSysStudies + 1; i <= nTPCAppSysStudies + nTPCnHitsStudies; ++i)
            hists.push_back(mainAnaHists[iGroup].hMainBcgSubtracted[i]);
        drawSysPlot(nameSuffix + "_TPCnHitsStudy",hists,histLegTPCnHits);

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
    double eBcg[nBins], eLum[nBins], eCorr[nBins], eTot[nBins];

    // Fill the arrays with bin content and systematic uncertainties
    for (int i = 0; i < nBins; ++i) {
        x[i] = hist->GetBinCenter(i + 1);
        y[i] = hist->GetBinContent(i + 1);
        exl[i] = hist->GetBinWidth(i + 1)*0.5;
        exh[i] = hist->GetBinWidth(i + 1)*0.5;
        double sysBcgErr = mainHists.hBcgSysStudy->GetBinContent(i + 1)/100;
        // add background normalisation uncertainity
        sysBcgErr = sqrt( pow(sysBcgErr,2) + pow(mBcgFracError[part][rpCon],2) );  
        zeros[i] = 0.0;
        eBcg[i] = y[i] * sysBcgErr;
        eLum[i] = y[i] * relLumUncertainty;
        eCorr[i] = y[i] * relativEffSys[part];
        eTot[i] = y[i] * sqrt(relLumUncertainty*relLumUncertainty + sysBcgErr*sysBcgErr + relativEffSys[part]*relativEffSys[part]);
    }
    // Create a TGraphAsymmErrors for systematic uncertainties
    systUncertainty = new TGraphAsymmErrors(nBins, x, y, exl, exh, eTot, eTot);
    systUncertainty->SetFillColor(14);
    systUncertainty->SetFillStyle(3006); // Set fill style (e.g., hatched)
    // Draw the histogram with statistical error bars and systematic uncertainties
    systUncertainty->Draw("E5 same"); // Draw systematic uncertainties as a filled area

    if( !drawCanvas )
        return; 


    enum SYSTUNC { bcg=0, lum, effCorr, total, nSysUnct }; 
    TGraphAsymmErrors* sysUnct[nSysUnct];

    TCanvas *cSysStudy; 
    CreateCanvas(&cSysStudy,"sysStudy_"+nameSuffix);
    sysUnct[total] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eTot, eTot);
    SetGraphStyle(sysUnct[total]);
    sysUnct[total]->SetTitle(";m(" + mUtil->pairLabel(part) + ") [GeV];d#sigma/dm(" + mUtil->pairLabel(part) + ") [nb/GeV]");
    sysUnct[total]->SetFillColor(14);
    sysUnct[total]->SetFillStyle(3006); // Set fill style (e.g., hatched)
    sysUnct[total]->Draw("");
    sysUnct[total]->Draw("E5 same");

    sysUnct[lum] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eLum, eLum);
    sysUnct[lum]->SetFillColor(3);
    sysUnct[lum]->SetFillStyle(3008);
    sysUnct[lum]->Draw("E5 same");

    sysUnct[effCorr] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eCorr, eCorr);
    sysUnct[effCorr]->SetFillColor(4);
    sysUnct[effCorr]->SetFillStyle(3009);
    sysUnct[effCorr]->Draw("E5 same");  

    sysUnct[bcg] = new TGraphAsymmErrors(nBins, x, zeros, exl, exh, eBcg, eBcg);
    sysUnct[bcg]->SetFillColor(2);
    sysUnct[bcg]->SetFillStyle(1001); 
    sysUnct[bcg]->Draw("E5 same"); 

    CreateLegend(0.64, 0.55, 0.9, 0.95);
    legend->AddEntry(sysUnct[total], "Systematic uncertainty", "f");
    legend->AddEntry(sysUnct[bcg], "Background sys. uncert.", "f");
    legend->AddEntry(sysUnct[lum], "Luminosity sys. uncert.", "f");
    legend->AddEntry(sysUnct[effCorr], "Eff. corr. sys. uncert.", "f");
    legend->Draw("same");

    WriteCanvas("sysStudy_"+nameSuffix,cSysStudy);
    cSysStudy->Close();
}//DrawSysUncertainty
     

void PlotManager::runNHitsVarStudy()
{
    for (unsigned int iVar = 0; iVar < nHitsVariation; ++iVar)
    {
        TTree *tree = dynamic_cast<TTree*>( inFile->Get( nameOfTree[kMAINANA] + "_" + mUtil->nHitVaryName(iVar) ) );
        if (!tree){ 
            cerr<<"Error: cannot loaded tree "<<nameOfTree[kMAINANA] << "_" << mUtil->nHitVaryName(iVar)<<"in PlotMainAna::runNHitsVarStudy()"<<endl; 
            return;
        }
        mCurrentTree = new RecTree(tree, treeBits[kMAINANA] );
        if (!mCurrentTree){ 
            cerr<<"Error: cannot loaded recTree in PlotMainAna::runNHitsVarStudy"<<endl; 
            return;
        }
        for(Long64_t iev=0; iev<tree->GetEntries(); ++iev)
        { //get the event
            tree->GetEntry(iev); 
            loopThroughTree(iVar);
        } 
    }
}//runNHitsVarStudy