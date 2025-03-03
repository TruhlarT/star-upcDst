#include "PlotManager.h"

using namespace std;
//using namespace RooFit;

void PlotManager::runMissingPt()
{
    changeDir(kMAINANA, "MissingPt");
    for (int iPart = 0; iPart < nParticles; ++iPart)
        for (int iRpCon = 0; iRpCon < nRpConfigurations +1; ++iRpCon)
            PlotMissingPt(iPart, iRpCon);
}

void PlotManager::PlotMissingPt(int part, int rpCon)
{
    if( DEBUG )
       cerr<<"PlotManager::PlotMissingPt() running for "<< mUtil->particleName(part)<<endl;
    const TString name[] = { TString("x"), TString("y"), TString("T") };
    const TString varName[] = { TString("pXMissing"), TString("pYMissing"), TString("pTMissing") };
    const TString configName = rpCon < nRpConfigurations ? mUtil->rpConfigName(rpCon) : ""; 
    TString nameSuffix = mUtil->particleName(part);
    if(rpCon < nRpConfigurations)
        nameSuffix += mUtil->rpConfigName(rpCon);

    changeSubDir(nameSuffix);

    // set number of bins for pX, pY and pT - missing
    // and indepenetly for pions, kaons and protons 
    int nBins[] = { 66, 66, 40};
    int nBinsDraw[] = { nBins[part], nBinsDraw[0], nBins[part]};

    // ranges for pX, pY and pT - missing
    double nBinsMaxPt[] = { 0.990, 0.990, 1.2};
    double binMaxDraw[] = {nBinsMaxPt[part]/2,binMaxDraw[0], nBinsMaxPt[part]};
    double binMinDraw[] = {-binMaxDraw[0],binMinDraw[0], 0.0};

    const int pTCutBins = 50;
    const double pTCutMinBin = 60.0;
    const double pTCutMaxBin = 2*pTCutBins + pTCutMinBin - 1;


    TH1F *hEff[12];
    // 0 - data passed  1 - data total
    // 2 - MC passed    3 - MC total
    // The first 4 = efficiency, the second signal to background and the last ones are the significance

    for (int i = 0; i < 12; ++i)
        hEff[i] = new TH1F(Form("hEff_%i",i) + nameSuffix, ";p_{T}^{miss} cut [MeV];", pTCutBins, pTCutMinBin, pTCutMaxBin);

    // placeholders for 
    double fitParam[4]; // pX/pY-missing gaussian fit parameters
    int cepSignal[2]; // integral of signail in pX/pY-missing obtained from the fit
    TString cuts = Form("pairID==%i", part);
    if( rpCon == ET)
        cuts += " && isElastic";
    else if( rpCon == IT)
        cuts += " && !isElastic";
    for (int i = 0; i < 3; ++i) // pX, pY, pT
    {

        CreateCanvas(&canvas,varName[i]+nameSuffix);
        if( part == PION)
            gPad->SetLeftMargin(0.1);
        mTree[kMAINANA]->Draw(varName[i] +Form(">>h(%i, %f, %f)",nBinsDraw[i],binMinDraw[i], binMaxDraw[i]),cuts);
        TH1F* hist = (TH1F*)gPad->GetPrimitive("h");

        hist->SetTitle(";p_{" + name[i] +Form("}^{miss} [GeV] ;Events / %.f MeV",hist->GetXaxis()->GetBinWidth(2)*1000 ));
        SetHistStyle(hist);
        //hist->GetXaxis()->SetTitleOffset(1.3);
        if( name[i] != "T")
            hist->GetYaxis()->SetRangeUser(0, 2*hist->GetMaximum() );   
        else
            hist->GetYaxis()->SetRangeUser(0, (part == PROTON ? 1.9 : 1.3)*hist->GetMaximum() ); 

        if( part == PION)
            hist->GetYaxis()->SetTitleOffset(1.9);

        hist->Draw("E");
        if( i < 2)
        {
            std::tuple<double, double, double> tmpTuple = FitpXYMissing(name[i], part, hist);

            fitParam[2*i] = std::get<0>(tmpTuple);
            fitParam[2*i + 1] = std::get<1>(tmpTuple);
            cepSignal[i] = std::get<2>(tmpTuple);
        }else{
            DrawExclusiveLine(exclusivityCut,0.0,exclusivityCut,hist->GetMaximum()*0.8);
            GenereateMCpTDist(hist, fitParam, part, rpCon, cepSignal, hEff);
        }

        DrawMainText(part, rpCon);

        mCurrDir->cd();
        WriteCanvas("missingP" + name[i]+nameSuffix);
        canvas->Close();
    }
    DrawEffStudy(part, rpCon, hEff);
    if( DEBUG )
       cerr<<"PlotManager::PlotMissingPt() is done..."<<endl;
}// PlotMissingPt

void PlotManager::DrawEffStudy(int part, int rpCon, TH1F* hEff[])
{
    TString nameSuffix = "_" + mUtil->particleName(part);
    if(rpCon < nRpConfigurations)
        nameSuffix += "_" + mUtil->rpConfigName(rpCon);
    if( DEBUG )
       cerr<<"PlotManager::DrawEffStudy() running for "<< nameSuffix<<endl;


    CreateCanvas(&canvas,"pTCutEffStudy"+nameSuffix);
    gPad->SetTicky(0);
    gPad->SetRightMargin(0.09);
    TGraphAsymmErrors *gEff = new TGraphAsymmErrors(hEff[0],hEff[1]);
    TGraphAsymmErrors *gEffMC = new TGraphAsymmErrors(hEff[2],hEff[3]);
    SetGraphStyle(gEff,4,20);
    SetGraphStyle(gEffMC,2,24);

    double yMax = 1.1;
    gEff->SetTitle(";p_{T}^{miss} cut [MeV];Eff,S/(S+B)");
    gEff->SetFillStyle(3003);
    gEff->SetFillColor(4);
    double xMax = gEff->GetXaxis()->GetXmax();
    gEff->GetYaxis()->SetRangeUser(0, yMax); 

    gEff->Draw("ap");
    gEffMC->SetFillStyle(3003);
    gEffMC->SetFillColor(2);
    gEffMC->Draw("p");

    TGraphAsymmErrors *gSigToBcg = new TGraphAsymmErrors(hEff[4], hEff[5]);
    TGraphAsymmErrors *gSigToBcgMC = new TGraphAsymmErrors(hEff[6], hEff[7]);
    SetGraphStyle(gSigToBcg,6,21);
    SetGraphStyle(gSigToBcgMC,8,25);

    gSigToBcg->SetFillStyle(3090);
    gSigToBcg->SetFillColor(6);
    gSigToBcg->Draw("p");

    gSigToBcgMC->SetFillStyle(3090);
    gSigToBcgMC->SetFillColor(8);
    gSigToBcgMC->Draw("p");

    hEff[8]->Divide(hEff[9]);
    SetHistStyle(hEff[8],4,22);
    hEff[8]->SetTitle(";;S/#sqrt{S+B}");


    float rightmax = 2*hEff[8]->GetMaximum();
    float scale = yMax/rightmax;

    hEff[8]->Scale(scale);
    TGraph *g = new TGraph(hEff[8]); 
    g->Draw("p");

    hEff[8]->Draw("same E");
    SetHistStyle(hEff[10],2,26);
    hEff[10]->Divide(hEff[11]);
    hEff[10]->Scale(scale);
    hEff[10]->Draw("same E");


    // draw an axis on the right side
    TGaxis *axis = new TGaxis(xMax,0,xMax, yMax,0,rightmax*yMax,510,"+L");
    axis->SetTitle("S/#sqrt{S+B}");
    SetAxisStyle(axis);
    axis->SetTitleOffset(1.5);
    axis->Draw();

    CreateLegend( 0.2, 0.15, 0.55, 0.35);
    legend->SetNColumns(2);
    legend->AddEntry(gEff, "Eff Data","ple");
    legend->AddEntry(gEffMC, "Eff MC","ple");
    legend->AddEntry(gSigToBcg, "S/(S+B) Data","ple");
    legend->AddEntry(gSigToBcgMC, "S/(S+B) MC","ple");
    legend->AddEntry(hEff[8], "S/#sqrt{S+B} Data","ple");
    legend->AddEntry(hEff[10], "S/#sqrt{S+B} MC","ple");
    legend->Draw("same");

    DrawExclusiveLine(exclusivityCut*1000,0.0,exclusivityCut*1000,1.1);

    mCurrDir->cd();
    WriteCanvas("effStudy"+ nameSuffix);
    canvas->Close();
    if( DEBUG )
       cerr<<"PlotManager::DrawEffStudy() is done..."<<endl;
}// DrawEffStudy

void PlotManager::DrawExclusiveLine(double xl, double yl, double xr, double yr)
{
    CreateLine(xl,yl, xr, yr);
    line->SetLineStyle(10);
    line->SetLineWidth(8);
    line->Draw("same");
}

std::tuple<double, double, double> PlotManager::FitpXYMissing( TString name, int part, TH1 *hist)
{
    if( DEBUG )
       cerr<<"PlotManager::FitpXYMissing() running for "<< mUtil->particleName(part)<<endl;
    double minValForPart[] = {0.3,0.4,0.6};
    double binMin = -minValForPart[part];
    double binMax = minValForPart[part];

    double sigWidthRange[] = { 0.04, 0.045, 0.04};

    // [x or y][part][ default value, min, max]
    double sigMeans[2][3][3] = {{{-0.04, -0.05, 0.03}, {-0.04, -0.05, 0.03}, {-0.02, -0.05, 0.03}}, // X
                         {{0.008, 0.0, 0.02}, {0.008, 0.0, 0.02}, {0.01, -0.01, 0.025}}}; // Y

    double sigMeanDef = sigMeans[ name == "x" ? 0 : 1][part][0];
    double sigMeanMin = sigMeans[ name == "x" ? 0 : 1][part][1];
    double sigMeanMax = sigMeans[ name == "x" ? 0 : 1][part][2];

    // Declare observable x
    RooRealVar x("x","p_{"+ name +"}^{miss} [GeV]",binMin,binMax) ;
    // Turn off RooFit info messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,RooFit::Import(*hist));
    RooPlot* frame = x.frame(RooFit::Title("")) ;
    dh.plotOn(frame,RooFit::DataError(RooAbsData::SumW2),RooFit::Name("data") ); 

    //define function to fit background - polynomial 2 degree
    RooRealVar a0("a0","a0",0.16,0.0,20.0);
    RooRealVar a1("a1","a1", -1.7,-50.0,0.0);
    RooPolynomial bkg("bkg","bkg",x,RooArgList(a0,a1));

    // define function to fit peak in data|
    RooRealVar sigmean("s1","s1",sigMeanDef,sigMeanMin,sigMeanMax);
    RooRealVar sigwidth("s2","s2",3.7e-02,0.025,sigWidthRange[part]);
    RooGaussian sig{"sig", "sig PDF", x, sigmean, sigwidth};

    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.7,0.,1.);
    RooRealVar nsig("nsig","signal events",60000,0,80000);
    RooRealVar nbkg("nbkg","background events",350000,0,450000);

    RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
    x.setRange("fitRange", binMin,binMax);

    //fit data with model
    model.fitTo(dh, RooFit::Save(true), RooFit::Range("fitRange"), RooFit::PrintLevel(-1));

    //plotting
    model.plotOn(frame, RooFit::Name("model") );
    CreateText(0.65,0.44,0.95,0.59);
    text -> SetTextAlign(12);
    text -> AddText(Form("#mu = %.1f #pm %.1f MeV",sigmean.getVal()*1000,sigmean.getError()*1000));
    text -> AddText(Form("#sigma = %.1f #pm %.1f MeV",sigwidth.getVal()*1000,sigwidth.getError()*1000));
    text -> AddText(Form("#chi^{2}/NDF = %.2f",frame->chiSquare()));
    text -> Draw("same");

    model.plotOn(frame, RooFit::Components(bkg),RooFit::LineStyle(kDashed), RooFit::LineColor(kRed), RooFit::Name("background") ); 
    frame->Draw("same hist E");

    CreateLegend(0.2, 0.44, 0.5, 0.59);
    legend->AddEntry(hist, "Data (unlike-sign pairs)","ple");
    legend->AddEntry("model", "Signal + background","l");
    legend->AddEntry("background", "Non-excl. bkgd","l");
    legend->Draw("same");
    if( DEBUG )
       cerr<<"PlotManager::FitpXYMissing() is done..."<<endl;
    return std::make_tuple(sigmean.getVal(), sigwidth.getVal(), nsig.getVal());
}//FitpXYMissing

void PlotManager::GenereateMCpTDist(TH1 *hist, double* fitParam, int part, int rpCon, int* cepSignal, TH1F* histToFill[])
{
    TString nameSuffix = "_" + mUtil->particleName(part);
    if(rpCon < nRpConfigurations)
        nameSuffix += "_" + mUtil->rpConfigName(rpCon);
    if( DEBUG )
       cerr<<"PlotManager::GenereateMCpTDist() running for "<< nameSuffix<<endl;
    Long64_t nev = nMCGenerratedExclusiveEvents;

    double pTMax = hist->GetXaxis()->GetBinLowEdge( hist->FindFixBin(nonExclusiveFitRangeMin[part]) + 1 );

    double fitLimitMin = pTMax;
    double fitLimitMax = nonExclusiveFitRangeMax[part];

    TH1F *hSig = new TH1F("hSig"+nameSuffix, "", int(pTMax/hist->GetXaxis()->GetBinWidth(2)), 0, pTMax);

    // Generate pTMissing distributions based on th pX/pY-Missing fit parametrs obtained by FitpXYMissing
    TRandom3 *rand = new TRandom3();
    for(Long64_t iev=0; iev<nev; ++iev) 
    { 
        if(iev%100000==0 && DEBUG)
            cout<<"Generating "<<iev<<". event..."<<endl;

        int seed = 1000000*(10*part + rpCon) + iev;
        rand->SetSeed(seed);
        double px = rand->Gaus(fitParam[0],fitParam[1]);
        double py = rand->Gaus(fitParam[2],fitParam[3]);
        hSig->Fill( sqrt( px*px + py*py) );

    }

    // Fit the non-exclusive background
    TF1* func = new TF1("fBcg", "[0]*x-[1]*x*x", 0, fitLimitMax);
    func->SetNpx(1e3);
    func->SetParameter(0, 1000000);
    func->SetParameter(1, 1000);
    func->SetParLimits(1,1,1000000);
    hist->Fit( func, "MQNI", "", fitLimitMin, fitLimitMax);
    func->SetLineColor(kRed);
    func->Draw("same");

    // Calculated the CEP sample = hist->Integral - nonExclusiveBackground->Integral
    int maxPtBin = hist->FindFixBin(pTMax);
    double cepMax = hist->Integral( 0, maxPtBin-1); // total size of signal
    cepMax += hist->GetBinContent(maxPtBin) * ( pTMax - hist->GetXaxis()->GetBinLowEdge(maxPtBin) ) 
    / hist->GetXaxis()->GetBinWidth(maxPtBin);
    cepMax -= func->Integral(0, pTMax)/hist->GetBinWidth(1);

    hSig->Scale(cepMax/hSig->Integral());
    SetHistStyle(hSig,1,25);
    //hSig->Draw("same E");

    // Calculted the CEp sample below pTMissing cut and corresponding efficiency of teh cut 
    double integralOfFunc = func->Integral(0, exclusivityCut)/hist->GetBinWidth(1); // background sample in the region below pT missing cut from fit
    double backgroundUncertainity = func->IntegralError(0, exclusivityCut)/hist->GetBinWidth(1); // background sample in the region below pT missing cut from fit
    int pTCutBin = hist->FindFixBin(exclusivityCut);
    double cepSample = hist->Integral( 0, pTCutBin-1); // sample of CEP candidates below pT missing cut 
    cepSample += hist->GetBinContent(pTCutBin) * ( exclusivityCut - hist->GetXaxis()->GetBinLowEdge(pTCutBin) ) 
    / hist->GetXaxis()->GetBinWidth(pTCutBin);


    double mcSample = hSig->Integral( 0, pTCutBin-1); // signal sample in the region below pT missing cut from MC simulation
    mcSample += hSig->GetBinContent(pTCutBin) * ( exclusivityCut - hSig->GetXaxis()->GetBinLowEdge(pTCutBin) ) 
    / hSig->GetXaxis()->GetBinWidth(pTCutBin);


    double effValMC = mcSample/hSig->Integral();
    double effVal = (cepSample-integralOfFunc)/cepMax;

    TH1F *hTot = (TH1F*)hSig->Clone("hTot");
    AddFunctionToHistogram(hTot, func);
    SetHistStyle(hSig,3,24);
    hTot->Draw("same E");

    bool showTextUp = !(part == PROTON || (part == PION && rpCon == ET)); 

    vector<double> vectorOfSignals;
    vectorOfSignals.push_back(cepSignal[0]);
    vectorOfSignals.push_back(cepSignal[1]);
    vectorOfSignals.push_back(cepMax);

    // set pt miss efficiency for later use
    mPtMissEff[part][rpCon] = (effVal+effValMC)*0.5;

    CreateText(0.28, showTextUp ? 0.6 : 0.1,0.35, showTextUp ? 0.8 : 0.3);
    text -> SetTextAlign(12);
    text -> AddText(Form("#varepsilon_{p_{T}^{miss} cut} = %.2f #pm %.2f%s",(effVal+effValMC)*50,abs(effVal-effValMC)*50,"%"));

    pair<double, double> result = mUtil->CalculateMeanAndError(vectorOfSignals);
    double ratio = result.first / (result.first + integralOfFunc);
    double ratio_uncertainty = sqrt( ( pow(integralOfFunc*result.second,2) + pow(result.first*backgroundUncertainity,2) ) /  pow(result.first + integralOfFunc,4) );
    // set background/ signal ratio 
    mBcgFraction[part][rpCon] = integralOfFunc/(result.first + integralOfFunc);
    mBcgFracError[part][rpCon] = ratio_uncertainty;
    //cout<<nameSuffix<<Form(": S/(S+B) = %.4f, sqrt( %.2f^2 * %.2f^2 + %.2f^2 * %.2f^2) / (%.2f+%.2f)^2", ratio_uncertainty, integralOfFunc, result.second, result.first, backgroundUncertainity, result.first, integralOfFunc)<<endl;
    text -> AddText(Form("S/(S+B) = %.3f #pm %.3f",ratio,ratio_uncertainty));
    //text -> AddText(Form("CEP: %i, %i, %.f",cepSignal[0],cepSignal[1],cepMax));
    text -> Draw("same");

    CreateLegend(0.58, showTextUp ? 0.6 : 0.1, 0.8, showTextUp ? 0.8 : 0.3);
    legend->AddEntry(hist, "Data (opposite-sign pairs)","ple");
    legend->AddEntry(hTot, "MC Sig+Bkgd","ple");
    //legend->AddEntry(hSig, "MC Signal","ple");
    legend->AddEntry(func, "Non-excl. bkgd","l");
    legend->Draw("same");

    // Fill histogram for effStudy
    for (int i = 0; i < histToFill[0]->GetNbinsX(); ++i)
    {
        double cut = (histToFill[0]->GetXaxis()->GetBinCenter(1) + 2*i)/1000.0; // From GeV to MeV
        int cutBin = hist->FindFixBin(cut);
        double cep = hist->Integral( 0, cutBin-1); 
        cep += hist->GetBinContent(cutBin) * ( cut - hist->GetXaxis()->GetBinLowEdge(cutBin) ) 
            / hist->GetXaxis()->GetBinWidth(cutBin);

        double mc = hSig->Integral( 0, cutBin-1);
        mc += hSig->GetBinContent(cutBin) * ( cut - hSig->GetXaxis()->GetBinLowEdge(cutBin) ) 
            / hSig->GetXaxis()->GetBinWidth(cutBin);

        double totMc = hTot->Integral( 0, cutBin-1);
        totMc += hTot->GetBinContent(cutBin) * ( cut - hTot->GetXaxis()->GetBinLowEdge(cutBin) ) 
            / hTot->GetXaxis()->GetBinWidth(cutBin);

        double intOfFunc = func->Integral(0, cut)/hist->GetBinWidth(1);

        int iBin = i +1; 
        double cepSigSample = cep-intOfFunc;

        histToFill[0]->SetBinContent(iBin, cepSigSample < cepMax ? cepSigSample : cepMax); 
        histToFill[1]->SetBinContent(iBin,cepMax); 

        histToFill[2]->SetBinContent(iBin, mc < cepMax ? mc: cepMax); 
        histToFill[3]->SetBinContent(iBin,cepMax);

        histToFill[4]->SetBinContent(iBin,cepSigSample < cep ? cepSigSample : cep); 
        histToFill[5]->SetBinContent(iBin,cep); 

        histToFill[6]->SetBinContent(iBin,mc < totMc ? mc : totMc); 
        histToFill[7]->SetBinContent(iBin,totMc); 

        histToFill[8]->SetBinContent(iBin,cepSigSample); 
        histToFill[9]->SetBinContent(iBin,sqrt(cep)); 
        histToFill[9]->SetBinError(iBin,0.5);

        histToFill[10]->SetBinContent(iBin,mc); 
        histToFill[11]->SetBinContent(iBin,sqrt(totMc)); 
        histToFill[11]->SetBinError(iBin,0.5);
    }
    if( DEBUG )
       cerr<<"PlotManager::GenereateMCpTDist() is done..."<<endl;
}//GenereateMCpTDist