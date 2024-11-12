#include "PlotManager.h"

void PlotManager::runVertexStudy()
{
   mStudyDir[kVERTEXSTUDY]->cd();

   int minFill = mTree[kVERTEXSTUDY]->GetMinimum("fillNumber");
   int maxFill = mTree[kVERTEXSTUDY]->GetMaximum("fillNumber");
   
   TH1D *hGausFitPar[5];
   TString hName[] = { TString("GausMean"), TString("GausSigma"), TString("EntriesPerFill"), TString("VertexEff") };
   hGausFitPar[0] = new TH1D("hGausMean", "", maxFill-minFill, minFill, maxFill);
   hGausFitPar[1] = new TH1D("hGausSigma", "", maxFill-minFill, minFill, maxFill);
   hGausFitPar[2] = new TH1D("hEntriesPerFill", "", maxFill-minFill, minFill, maxFill);
   hGausFitPar[3] = new TH1D("hVertexEff", "", maxFill-minFill, minFill, maxFill);
   hGausFitPar[4] = new TH1D("hVertexEffData", "", maxFill-minFill, minFill, maxFill);

   double mean, sigma, eff, entries, dataEff;
   double passed, total;
   int fill;
   // loaded in PlotManager.cxx
   // be aware when applying changes 
   TTree *vertexTree = new TTree("vertexTree", "vertexTree");
   vertexTree->Branch("fill", &fill);
   vertexTree->Branch("mean", &mean);
   vertexTree->Branch("sigma", &sigma);   
   vertexTree->Branch("eff", &eff);
   vertexTree->Branch("dataEff", &dataEff);
   vertexTree->Branch("entries", &entries);

   changeDir(kVERTEXSTUDY, "VertexPlots");

   //for (int iFill = minFill; iFill <= minFill+1; ++iFill)
   for (int iFill = minFill; iFill <= maxFill; ++iFill)
   {
      fill = iFill;
      CreateCanvas(&canvas,Form("%i",iFill));
      mTree[kVERTEXSTUDY]->Draw("vertexZInCm>>htemp(44,-220,220)",Form("fillNumber == %i && TMath::Abs(vertexZInCm) < %f",iFill,vertexRange));
      TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
      if( !htemp ) 
         continue;
      if( !htemp->GetEntries() )
         continue;
      passed = htemp->GetEntries();
      SetHistStyle(htemp, mainColor, mainMarker);
      htemp->SetTitle(";z_{vrtx} [cm]; Events / 10 cm");
      htemp->GetYaxis()->SetTitleOffset(2.7);
      htemp->Draw("E");
      DrawSTARInternal();
      CreateText(0.75,0.84,0.88,0.88);
      text -> AddText(Form("Fill %i",iFill));   
      text -> Draw("same");
      WriteCanvas(Form("%i_cut",iFill));
      mTree[kVERTEXSTUDY]->Draw("vertexZInCm>>htemp(44,-220,220)",Form("fillNumber == %i",iFill));
      htemp = (TH1F*)gPad->GetPrimitive("htemp");
      if( !htemp ) 
         continue;
      if( !htemp->GetEntries() )
         continue;
      total = htemp->GetEntries();
      htemp->Fit("gaus","Q","",-zVertexFitRange,zVertexFitRange);
      TF1 *gaus = (TF1*)htemp->GetListOfFunctions()->FindObject("gaus");
      if(!gaus)
         continue;
      gaus->SetLineColor(mainColor);
      TF1 *gausExtended = (TF1*)gaus->Clone();
      gausExtended->SetRange(-200,200);
      gausExtended->SetLineColor(bckgColor);

      SetHistStyle(htemp, mainColor, mainMarker);
      htemp->SetTitle(";z_{vrtx} [cm]; Events / 10 cm");
      htemp->GetYaxis()->SetTitleOffset(2.7);
      htemp->Draw("E");
      gausExtended->Draw("same");
      gaus->Draw("same");
      DrawSTARInternal();
      DrawLine(htemp,vertexRange, true);

      CreateText(0.75,0.84,0.88,0.88);
      text -> AddText(Form("Fill %i",iFill));   
      text -> Draw("same");
      WriteCanvas(Form("%i_main",iFill));
      
      mean = gaus->GetParameter(1);
      sigma = gaus->GetParameter(2);
      eff = 0.5*( TMath::Erf((vertexRange - mean)/(sqrtOfTwo * sigma))  
         - TMath::Erf((-vertexRange - mean)/(sqrtOfTwo * sigma)) );
      double newMean = mean - gaus->GetParError(1);
      double newSigma = sigma - gaus->GetParError(2);
      double effErr = 0.5*( TMath::Erf((vertexRange - newMean)/(sqrtOfTwo * newSigma))  
         - TMath::Erf((-vertexRange - newMean)/(sqrtOfTwo * newSigma)) ); 
      effErr -= eff;
      entries = htemp->GetEntries();

      hGausFitPar[0]->SetBinContent(iFill + 1 - minFill, mean);
      hGausFitPar[0]->SetBinError(iFill + 1 - minFill, gaus->GetParError(1));

       
      hGausFitPar[1]->SetBinContent(iFill + 1 - minFill, sigma);
      hGausFitPar[1]->SetBinError(iFill + 1 - minFill, gaus->GetParError(2));
      
      hGausFitPar[2]->SetBinContent(iFill + 1 - minFill, entries);
      
      hGausFitPar[3]->SetBinContent(iFill + 1 - minFill, eff);
      hGausFitPar[3]->SetBinError(iFill + 1 - minFill, effErr); 

      dataEff = passed / total;
      hGausFitPar[4]->SetBinContent(iFill + 1 - minFill, dataEff);
      hGausFitPar[4]->SetBinError(iFill + 1 - minFill, dataEff*sqrt( (passed+total)/(passed*total) ));

      canvas->Close();

      vertexTree->Fill();
   }

   mStudyDir[kVERTEXSTUDY]->cd();
   TString yName[] = { TString("<z_{vrtx}> [cm]"), TString("#sigma(z_{vrtx}) [cm]"), TString("Entries"), TString("Eff")};
   double yMinRange[] = { -10, 40};
   double yMaxRange[] = { 10, 80};
   for (int iPar = 0; iPar < 3; ++iPar)
   {
      CreateCanvas(&canvas,hName[iPar]);
      SetHistStyle(hGausFitPar[iPar], mainColor, mainMarker);
      hGausFitPar[iPar]->SetTitle(";Fill number;" + yName[iPar]);
      if( iPar < 2)
         hGausFitPar[iPar]->GetYaxis()->SetRangeUser(yMinRange[iPar], yMaxRange[iPar]);
      hGausFitPar[iPar]->Draw("E");

      DrawSTARInternal();
      WriteCanvas(hName[iPar]);
      canvas->Close();
   }

   CreateCanvas(&canvas,hName[3]);
   SetHistStyle(hGausFitPar[3], mainColor, mainMarker);
   hGausFitPar[3]->SetTitle(";Fill number;" + yName[3]);
   SetHistStyle(hGausFitPar[4], bckgColor, mainMarker);
   TH1 *hRatio;
   DrawRatioPlot(hGausFitPar[3], hGausFitPar[4], &hRatio);
   hGausFitPar[3]->GetYaxis()->SetRangeUser(0.66, 1.0);

   hRatio->Divide(hGausFitPar[4]);
   hRatio->Draw("ep");       // Draw the ratio plot

   // Fit the constant function to the histogram
   TF1 *f1 = new TF1("f1", "[0]");
   hRatio->Fit(f1);
   f1->SetLineColor(kRed);
   f1->SetLineStyle(10);
   f1->Draw("same"); 

   CreateText(0.15,0.84,0.25,0.88);
   text -> SetTextColor(kRed);
   text -> AddText(Form("%.3f",f1->GetParameter(0)));   
   text -> Draw("same");

   hRatio->GetYaxis()->SetRangeUser(1.0, 1.045);
   hRatio->GetYaxis()->SetTitle("Fit-data");
   hRatio->GetXaxis()->SetTitle("Fill number");
   canvas->cd();
   TPad* pad = (TPad*) canvas->FindObject(TString::Format("pad1").Data());
   pad->cd();
   CreateLegend(0.65, 0.2, 0.85, 0.35);
   legend->AddEntry(hGausFitPar[3], "Obtained from fit","ep");
   legend->AddEntry(hGausFitPar[4], "Obtained from data","ep");
   legend->Draw("same");
   DrawSTARInternal(0.65, 0.4, 0.85, 0.45);
   //canvas->Update();
   canvas->cd(0);
   WriteCanvas(hName[3]);
   canvas->Close();


   vertexTree->Write();
}