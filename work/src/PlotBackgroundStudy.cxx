#include "PlotManager.h"

TH1D* PlotManager::bkgdHistogramComparison(TH2 *hist, int part, int rpCon, bool study) // bin-by-bin method
{
    // Allocate memory for funcVector
    vector<TF1*> funcVector;

    double fitLimitMin = nonExclusiveFitRangeMin[part];
    double fitLimitMax = nonExclusiveFitRangeMax[part];
    TH1D *hBkgdHistogram = (TH1D*)bkgdHistogram(hist, exclusivityCut, fitLimitMin, fitLimitMax, 0, &funcVector);

    if(!study)
        return hBkgdHistogram;

    changeDir(kMAINANA, "BackgroundStudy");
    TString subDirName = mUtil->particleName(part);
    if(rpCon < nRpConfigurations)
        subDirName += mUtil->rpConfigName(rpCon);

    changeSubDir(subDirName);

    TString nameSuffix = hist->GetName();
    nameSuffix = nameSuffix(19, nameSuffix.Length()-19);
    // Create canvas outside the loop
    CreateCanvas(&canvas, "pTMissProjection" + nameSuffix);

    for(int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        TH1D* hMissPtProj = hist->ProjectionY("tmpNameProjectionY_singleBin", i, i);
        SetHistStyle(hMissPtProj);

        hMissPtProj->SetTitle(";p_{T}^{miss} [GeV];d#sigma/dX [nb/GeV]");
        hMissPtProj->GetYaxis()->SetRangeUser(0, 1.2*hMissPtProj->GetMaximum() );  
        hMissPtProj->Draw("hist E");

        if(funcVector[i - 1] != nullptr) // Access funcVector and binsVector with proper indexing
            funcVector[i - 1]->Draw("same");

        CreateText(0.6,0.14,0.95,0.22); 
        text -> AddText(Form("m(" + mUtil->pairLabel(part) + ") #in [%.2f;%.2f] GeV",hBkgdHistogram->GetBinLowEdge(i),hBkgdHistogram->GetBinLowEdge(i)+hBkgdHistogram->GetBinWidth(i)));
        text -> Draw("same");

        DrawMainText(part, rpCon);
        // Write canvas to output file
        WriteCanvas("pTMissProjection"+nameSuffix+Form("%f", hBkgdHistogram->GetBinCenter(i) ));

        // Delete hMissPtProj to avoid memory leaks
        delete hMissPtProj;
    }

    // Close canvas after loop if necessary
    canvas->Close();
    return hBkgdHistogram;
}//bkgdHistogramComparison

Double_t PlotManager::bkgdEvents(const TH1* hPtMiss, Double_t ptMissCut, Double_t fitLimitMin, Double_t fitLimitMax, TF1 *funcReturn) const
{
   static int funcId; ++funcId; 
   TString idStr; idStr.Form("func%d", funcId);
   TF1* func = new TF1(idStr, "[0]*x-[1]*x*x", 0, fitLimitMax);
   func->SetNpx(1e3);
   func->SetParameter(0, 1000000);
   func->SetParameter(1, 1000);
   func->SetParLimits(1,0,10000000);
   double nBkgdEvents = 1e-9;
   if(!hPtMiss){
      std::cerr << "ERROR in PlotManager::bkgdEvents()" << std::endl;
      return nBkgdEvents;
   }

   const int firstBinInFitRange = hPtMiss->GetXaxis()->FindBin(fitLimitMin);
   const int lastBinInFitRange = hPtMiss->GetXaxis()->FindBin(fitLimitMax)-1;
   // check if the fitting range is not empty
   if( hPtMiss->Integral( firstBinInFitRange, lastBinInFitRange) < 1e-5 )
      return nBkgdEvents; 
   // check if there is low statistic / some bins are empty
   for(int i=firstBinInFitRange; i<=lastBinInFitRange; ++i)
      if(hPtMiss->GetBinContent(i) < 1e-5)
         func->FixParameter(1, 0); 

   // Rebin so that there are 4 bins to fit
   //const int nBins = lastBinInFitRange - firstBinInFitRange + 1;
   TH1F* hTemp = new TH1F();
   hPtMiss->Copy( *hTemp );
   hTemp->SetName("VeryTemp");
   //const int nRebin = static_cast<int>( nBins/4 ); 
   //if( nRebin>1 )
      //hTemp->Rebin(nRebin);

   // fit the background by pol2 in fit range outside the signal region
   hTemp->Fit( func, "MQNI", "", fitLimitMin, fitLimitMax);
   if( fabs(func->GetParameter(0)-1) < 0.001 && fabs(func->GetParameter(1)+1) < 0.001 )
   {
      nBkgdEvents = 0;
   } else
   {
      double integralOfFunc = func->Integral(0, ptMissCut);
      if( integralOfFunc < 0 ) // use linear approximation
      {
         func->FixParameter(1, 0);
         hTemp->Fit( func, "MQNI", "", fitLimitMin, fitLimitMax);
         integralOfFunc = func->Integral(0, ptMissCut);
      }
      nBkgdEvents = integralOfFunc / hTemp->GetBinWidth(1);
   }

   delete hTemp;
   
   if(funcReturn)
   {
      /*if( nRebin>1 )
      {
         func->SetParameter(0, func->GetParameter(0)/nRebin );
         func->SetParameter(1, func->GetParameter(1)/nRebin );
      }*/
      *funcReturn = *func;
   }else 
      delete func;
   
   return nBkgdEvents;
}


Double_t PlotManager::integratePtMiss(const TH1* hPtMiss, Double_t ptMissCut) const
{
   int ptMissCutBin = hPtMiss->FindFirstBinAbove(ptMissCut);
   double integral = hPtMiss->Integral(1, ptMissCutBin-1);
   integral += hPtMiss->GetBinContent(ptMissCutBin) * ( ptMissCut - hPtMiss->GetXaxis()->GetBinLowEdge(ptMissCutBin) ) / hPtMiss->GetXaxis()->GetBinWidth(ptMissCutBin);
   return integral;
}



TH1D* PlotManager::bkgdHistogram(const TH2* hMissPtVsX, Double_t ptMissCut, Double_t fitLimitMin, Double_t fitLimitMax, Int_t mode, vector<TF1*> * funcVec) const
{
    TH1D *hBkgd = hMissPtVsX->ProjectionX(TString("Background_")+hMissPtVsX->GetName());
    hBkgd->Reset("ICESM");
    for(int i=1; i<=(hBkgd->GetNbinsX()); ++i)
    {
        TF1 *pTMissExtrapolationFunc = new TF1();
        TH1D* hMissPtProj = hMissPtVsX->ProjectionY("tmpNameProjectionY_singleBin", i, i);
        double nBkgdEvents = bkgdEvents( hMissPtProj, ptMissCut, fitLimitMin, fitLimitMax, pTMissExtrapolationFunc);
        hBkgd->SetBinContent(i, nBkgdEvents);
        //hBkgd->SetBinError(i, sqrt(nBkgdEvents) ); 
        hBkgd->SetBinError(i, 0.0001 * hBkgd->GetBinContent(i) );
        if(funcVec) 
            funcVec->push_back( pTMissExtrapolationFunc );
    }
    return hBkgd;
}


TH1D* PlotManager::backgroundStudy(TH2* hMissPtVsX, int part, int rpCon, TString hName, int bcgStudy)
{
    int leftBin = hMissPtVsX->GetYaxis()->FindBin( nonExclFitRangeMin[part] );
    int rightBin = hMissPtVsX->GetYaxis()->FindBin( nonExclFitRangeMax[part] );
    if(bcgStudy > 0)
    {
        //cout<<Form("backgroundStudy(): %i default range is from %i to %i",bcgStudy, leftBin, rightBin)<<endl;
        leftBin += bcgStudy == 1 ? 1 : -2;
        rightBin += bcgStudy == 1 ? -1 : 3;
        //cout<<Form("backgroundStudy(): %i range is from %i to %i",bcgStudy, leftBin, rightBin)<<endl;
    }
    TH1D *hBkgd = hMissPtVsX->ProjectionX(hName+hMissPtVsX->GetName(), leftBin, rightBin );
    hName += "_SignalPart_";
    TH1D *hSignal = hMissPtVsX->ProjectionX(hName+hMissPtVsX->GetName(), 1, hMissPtVsX->GetYaxis()->FindBin( exclusivityCut ) );

    //cout<<"For hist "<<hName+hMissPtVsX->GetName()<<" signal is "<<hSignal->Integral()<<" and bcg "<<hBkgd->Integral() <<endl;
    double scaleFactor = double(hSignal->Integral())/double(hBkgd->Integral());
    if( mBcgFraction[part][rpCon] > 0)
        scaleFactor*=mBcgFraction[part][rpCon];
    else{
        cout<<"Warrning in PlotManager::backgroundStudy() the total size of non-exclusive sample was not set"<<endl;
        cout<<"    the background is not correctly plotted"<<endl;
        cout<<"    assuming 20% bcg contribution"<<endl;
        scaleFactor*=0.2; //Assuming 20% bcg contribution
    }
    hBkgd->Scale(scaleFactor); 

    return hBkgd;
}


void PlotManager::subtractBackground(TH1* hSignalPlusBkgd, const TH1* hBkgd) const
{
   if( hSignalPlusBkgd->GetNbinsX() != hBkgd->GetNbinsX() )
   {
      std::cerr << "ERROR in Util::subtractBackground(): number of bins of signal ("<< hSignalPlusBkgd->GetName()<<") and background ("<< hBkgd->GetName()<<") histogram is different. Return" << std::endl;
      return;
   }
   hSignalPlusBkgd->Add(hBkgd, -1.);

   for(int i=1; i<=(hSignalPlusBkgd->GetNbinsX()); ++i)
   {
      const double signalPlusBkgdContent = hSignalPlusBkgd->GetBinContent(i);
      if( signalPlusBkgdContent>1e-9 ) 
         continue;
      if( DEBUG )
         std::cerr << "WARNING in subtractBackground(): Background larger than signal" << std::endl;
      hSignalPlusBkgd->SetBinContent(i, 1e-9 );
      hSignalPlusBkgd->SetBinError(i, 1e-9 );
   }
}      
