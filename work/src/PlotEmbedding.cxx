#include "PlotManager.h"

using namespace std;


void PlotManager::runEmbeddingStudy()
{  
   changeDir(kEMBEDING, "TpcEff");
   if( DEBUG )
      cerr<<"PlotManager::runEmbeddingStudy() going to runTpcEff"<<endl;
   runTpcEff();
   changeDir(kEMBEDING, "TPCSysStudy");
   if( DEBUG )
      cerr<<"PlotManager::runEmbeddingStudy() going to run runTpcEmbedSysStudy()"<<endl;   
   runTpcEmbedSysStudy();

   changeDir(kEMBEDING, "TofEff");
   if( DEBUG )
      cerr<<"PlotManager::runEmbeddingStudy() going to runTofEff"<<endl;
   runTofEff();

   if( DEBUG )
      cerr<<"PlotManager::runEmbeddingStudy() going to run runEmbedVarDiff()"<<endl;
   runEmbedVarDiff();
}//runEmbeddingStudy

void PlotManager::initEmbeddingStudy( unsigned int particle)
{
   mTree[kEMBEDING] = dynamic_cast<TTree*>( embFile[particle]->Get( nameOfTree[kEMBEDING] ) );
   embedParticle = particle;
   if (!mTree[kEMBEDING])
      cerr<<"Error: cannot open "<<nameOfTree[kEMBEDING]<<endl;
}

void PlotManager::runTpcEff()
{
   // plot pT, eta and phi for the TrueMc and reco
   // divide them to get the eff as a function fo pT, eta and phi
   CreateCanvas(&canvas,"tpcEff");
   canvas->cd();
   canvas->SetLogy();

   TH1F *histReco, *histGenerated;
   TH2F *hTot, *hPassed;

   TString piName[] = {mUtil->particleName(embedParticle) + mUtil->signName(PLUS), mUtil->particleName(embedParticle) + mUtil->signName(MINUS)};
   //TString piTexName[] = { mUtil->particleTag(embedParticle, PLUS), mUtil->particleTag(embedParticle, MINUS)};

   const unsigned int nVal = 4;
   const char* valName[] = {"pTInGev", "phi", "eta", "vertexZInCm"};
   TString valLabel[] = {"p_{T} [GeV]", "#phi", "#eta", "true z_{vrtx}"};
   bool pairFlag[] = {1,1,1,0};

   unsigned int nBins[] = { 80, 40, 40, 20 };
   double min[] = { 0.0, -TMath::Pi(), -1.0,-100 }; 
   double max[] = { 3.0, TMath::Pi(), 1.0,100 }; 

   TString cutLabel[] = { Form("|z_{vrtx}| < %.f cm", vertexRange), Form("|#eta| < %.1f", maxEta), Form("p_{T} > %.1f GeV", minPt[embedParticle])};


   const unsigned int nCuts = 3;
   TString cuts[nCuts];
   TString finalCuts;

   TString valCuts;
   TString valCutLabel;
   bool cutsUsedForVal[ nVal ][ nCuts ] = { { 1,1,0}, { 1,1,1}, { 1,0,1}, { 0,1,1}};


   for (int iPi = 0; iPi < 2; ++iPi)
   {
      changeSubDir(piName[iPi]);

      cuts[0] = Form("TMath::Abs(vertexZInCm_TrueMc) < %f", vertexRange);
      cuts[1] = Form("TMath::Abs(eta%i_TrueMc) < %f && eta%i_TrueMc > %f*vertexZInCm_TrueMc - %f && eta%i_TrueMc < %f*vertexZInCm_TrueMc + %f", iPi, maxEta, iPi, etaVertexSlope, etaVertexShift, iPi, etaVertexSlope, etaVertexShift);
      cuts[2] = Form("pTInGev%i_TrueMc > %f", iPi, minPt[embedParticle]);
      canvas->cd();
      gPad->SetMargin(0.07,0.02,0.105,0.03); // reset margin
      for (unsigned int iVal = 0; iVal < nVal; ++iVal)
      {
         valCuts = "";
         valCutLabel = "";
         for (unsigned int i = 0; i < nCuts; ++i){
            if( !cutsUsedForVal[iVal][i] )
               continue;
            valCuts += cuts[i] + " && ";
            valCutLabel += cutLabel[i] + + ", ";
         }
         valCuts.Resize(valCuts.Length() - 4);
         valCutLabel.Resize(valCutLabel.Length() - 2);

         finalCuts = valCuts;
         if( pairFlag[iVal] )
            mTree[kEMBEDING]->Draw(Form("%s%i_TrueMc>>h2(%i,%f,%f)",valName[iVal], iPi, nBins[iVal], min[iVal], max[iVal]),finalCuts);
         else
            mTree[kEMBEDING]->Draw(Form("%s_TrueMc>>h2(%i,%f,%f)",valName[iVal], nBins[iVal], min[iVal], max[iVal]),finalCuts);
         histGenerated = (TH1F*)gPad->GetPrimitive("h2");
         SetHistStyle(histGenerated);
         histGenerated->SetTitle(";" + valLabel[iVal] + ";counts");
         canvas->SetLogy();
         CreateLegend(0.74, 0.74, 0.97, 0.89);
         legend->AddEntry(histGenerated, "True MC","ple");
         legend->Draw("same");
         WriteCanvas(TString(valName[iVal]) + "_TrueMc");

         if( finalCuts != "")
            finalCuts+=" && ";

         finalCuts += Form("idTruth%i != -1 && ",iPi) + assosiationCut(iPi);
         if( pairFlag[iVal] )
            mTree[kEMBEDING]->Draw(Form("%s%i_TrueMc>>h(%i,%f,%f)",valName[iVal], iPi, nBins[iVal], min[iVal], max[iVal]),finalCuts);
         else
            mTree[kEMBEDING]->Draw(Form("%s_TrueMc>>h(%i,%f,%f)",valName[iVal], nBins[iVal], min[iVal], max[iVal]),finalCuts);
         histReco = (TH1F*)gPad->GetPrimitive("h");
         SetHistStyle(histReco);
         histReco->SetTitle(";" + valLabel[iVal] + ";counts"); 
         CreateLegend(0.74, 0.74, 0.97, 0.89);
         legend->AddEntry(histGenerated, "MC reconstructed","ple");
         legend->Draw("same");
         WriteCanvas(TString(valName[iVal]));

         TEfficiency* pEff = new TEfficiency(*histReco,*histGenerated);
         TGraphAsymmErrors *gr = pEff->CreateGraph();
         SetGraphStyle(gr);

         gr->SetTitle(Form(";%s;Eff", valLabel[iVal].Data())); 
         //gr->SetTitle(Form("%s, %s;%s;Eff", piTexName[iPi].Data(), valCutLabel.Data(), valLabel[iVal].Data())); 
         gr->GetYaxis()->SetRangeUser(0.35, 1.0);
         gr->Draw("ap");
         
         if( iVal == 0){
            CreateLine(minPt[embedParticle],0.35,minPt[embedParticle],0.9);
            line->SetLineWidth(4);
            line->Draw("same");            
         }

         canvas->SetLogy(0);
         WriteCanvas(TString(valName[iVal]) + "_Eff");
      }


      gPad->SetMargin(0.07,0.11,0.105,0.03); // 2D Margin

      finalCuts = cuts[0] + " && " + cuts[1] + " && " + cuts[2];
      mTree[kEMBEDING]->Draw(Form("phi%i_TrueMc:eta%i_TrueMc>>hTot(%i,%f,%f,%i,%f,%f)",iPi, iPi, nBins[2], min[2], max[2],nBins[1], min[1], max[1]),finalCuts,"colz");
      hTot = (TH2F*)gPad->GetPrimitive("hTot");


      finalCuts += Form("&& idTruth%i != -1 && ",iPi) + assosiationCut(iPi);
      mTree[kEMBEDING]->Draw(Form("phi%i_TrueMc:eta%i_TrueMc>>hPassed(%i,%f,%f,%i,%f,%f)",iPi, iPi, nBins[2], min[2], max[2],nBins[1], min[1], max[1]),finalCuts,"colz");
      hPassed = (TH2F*)gPad->GetPrimitive("hPassed");

      hPassed->Divide(hTot);
      SetHistStyle(hPassed);
      hPassed->GetZaxis()->SetTitleOffset(1.2);
      hPassed->GetZaxis()->SetRangeUser(0.35, 1.0);
      hPassed->SetTitle(";" + valLabel[2] + ";" + valLabel[1] + ";Eff"); 
      SetPalletRange(hPassed);
      WriteCanvas("PhiEta_Eff");

      finalCuts = cuts[2];
      mTree[kEMBEDING]->Draw(Form("vertexZInCm_TrueMc:eta%i_TrueMc>>hTot(%i,%f,%f,%i,%f,%f)", iPi,nBins[2], min[2], max[2], nBins[3], min[3], max[3]),finalCuts,"colz");
      hTot = (TH2F*)gPad->GetPrimitive("hTot");

      finalCuts += Form("&& idTruth%i != -1 && ",iPi) + assosiationCut(iPi);
      mTree[kEMBEDING]->Draw(Form("vertexZInCm_TrueMc:eta%i_TrueMc>>hPassed(%i,%f,%f,%i,%f,%f)", iPi, nBins[2], min[2], max[2], nBins[3], min[3], max[3]),finalCuts,"colz");
      hPassed = (TH2F*)gPad->GetPrimitive("hPassed");
      hPassed->Divide(hTot);
      SetHistStyle(hPassed);
      hPassed->GetZaxis()->SetTitleOffset(1.2);
      hPassed->GetZaxis()->SetRangeUser(0.35, 1.0);
      hPassed->SetTitle(";" + valLabel[2] + ";" + valLabel[3] + ";Eff"); 
      drawVertexEtaSpace();
      SetPalletRange(hPassed);
      WriteCanvas("VertexEta_Eff");
      plotAverageEtaVertexEff(hPassed,"VertexEta_Eff");
   }
   canvas->Close();
}//runTpcEff


void PlotManager::runTofEff()
{
   CreateCanvas(&canvas,"tofEff");
   canvas->cd();
   canvas->SetLogy();

   TH1F *histReco, *histGenerated;
   TH2F *hTot, *hPassed;

   TString piName[] = {mUtil->particleName(embedParticle) + mUtil->signName(PLUS), mUtil->particleName(embedParticle) + mUtil->signName(MINUS)};
   //TString piTexName[] = { mUtil->particleTag(embedParticle, PLUS), mUtil->particleTag(embedParticle, MINUS)};

   const unsigned int nVal = 4;
   const char* valName[] = {"pTInGev", "phi", "eta", "vertexZInCm"};
   TString valLabel[] = {"p_{T} [GeV]", "#phi", "#eta", "true z_{vrtx}"};
   bool pairFlag[] = {1,1,1,0};

   unsigned int nBins[] = { 80, 40, 40, 20 };
   double min[] = { 0.0, -TMath::Pi(), -1.0,-100 }; 
   double max[] = { 3.0, TMath::Pi(), 1.0,100 }; 


   TString cutLabel[] = { Form("|z_{vrtx}| < %.f cm", vertexRange), Form("|#eta| < %.1f", maxEta), Form("p_{T} > %.1f GeV", minPt[embedParticle])};

   const unsigned int nCuts = 3;
   TString cuts[nCuts];
   TString finalCuts;

   TString valCuts;
   TString valCutLabel;
   bool cutsUsedForVal[ nVal ][ nCuts ] = { { 1,1,0}, { 1,1,1}, { 1,0,1}, { 0,1,1}};
   TString assosiationCuts;

   for (int iPi = 0; iPi < 2; ++iPi)
   {
      changeSubDir(piName[iPi]);

      cuts[0] = Form("TMath::Abs(vertexZInCm_TrueMc) < %f", vertexRange);
      cuts[1] = Form("TMath::Abs(eta%i_TrueMc) < %f", iPi, maxEta);
      cuts[2] = Form("pTInGev%i_TrueMc > %f", iPi, minPt[embedParticle]);
 
      assosiationCuts = Form("idTruth%i != -1 && ",iPi) + assosiationCut(iPi);
      canvas->cd();
      gPad->SetMargin(0.07,0.02,0.105,0.03); // reset margin
      for (unsigned int iVal = 0; iVal < nVal; ++iVal)
      {
         valCuts = ""; valCutLabel = "";
         for (unsigned int i = 0; i < nCuts; ++i){
            if( !cutsUsedForVal[iVal][i] )
               continue;
            valCuts += cuts[i] + " && ";
            valCutLabel += cutLabel[i] + + ", ";
         }
         valCuts.Resize(valCuts.Length() - 4);
         valCutLabel.Resize(valCutLabel.Length() - 2);

         finalCuts = assosiationCuts + " && " + valCuts;
         if( pairFlag[iVal] )
            mTree[kEMBEDING]->Draw(Form("%s%i>>h2(%i,%f,%f)",valName[iVal], iPi, nBins[iVal], min[iVal], max[iVal]),finalCuts);
         else
            mTree[kEMBEDING]->Draw(Form("%s_TrueMc>>h2(%i,%f,%f)",valName[iVal], nBins[iVal], min[iVal], max[iVal]),finalCuts);
         histGenerated = (TH1F*)gPad->GetPrimitive("h2");
         SetHistStyle(histGenerated);
         histGenerated->SetTitle(";" + valLabel[iVal] + ";counts");
         canvas->SetLogy();
         CreateLegend(0.74, 0.74, 0.97, 0.89);
         legend->AddEntry(histGenerated, "Reco TPC track","ple");
         legend->Draw("same");
         WriteCanvas(TString(valName[iVal]));

         finalCuts += Form(" && tofMatched%i", iPi);
         if( pairFlag[iVal] )
            mTree[kEMBEDING]->Draw(Form("%s%i>>h(%i,%f,%f)",valName[iVal], iPi, nBins[iVal], min[iVal], max[iVal]),finalCuts);
         else
            mTree[kEMBEDING]->Draw(Form("%s_TrueMc>>h(%i,%f,%f)",valName[iVal], nBins[iVal], min[iVal], max[iVal]),finalCuts);
         histReco = (TH1F*)gPad->GetPrimitive("h");
         SetHistStyle(histReco);
         histReco->SetTitle(";" + valLabel[iVal] + ";counts"); 
         CreateLegend(0.74, 0.74, 0.97, 0.89);
         legend->AddEntry(histReco, "Reco TPC-TOF track","ple");
         legend->Draw("same");
         WriteCanvas(TString(valName[iVal]) + "_Tof");

         TEfficiency* pEff = new TEfficiency(*histReco,*histGenerated);
         TGraphAsymmErrors *gr = pEff->CreateGraph();
         SetGraphStyle(gr);

         gr->GetYaxis()->SetRangeUser(0.35, 1.0);
         gr->SetTitle(Form(";%s;Eff", valLabel[iVal].Data())); 
         //gr->SetTitle(Form("%s, %s;%s;Eff", piTexName[iPi].Data(), valCutLabel.Data(), valLabel[iVal].Data())); 
         gr->Draw("ap");

         if( iVal == 0){
            CreateLine(minPt[embedParticle],0.0,minPt[embedParticle],0.8);
            line->SetLineWidth(4);
            line->Draw("same");            
         }

         canvas->SetLogy(0);
         WriteCanvas(TString(valName[iVal]) + "_TofEff");        
      }


      gPad->SetMargin(0.07,0.11,0.105,0.03); // 2D Margin

      finalCuts = assosiationCuts;
      mTree[kEMBEDING]->Draw(Form("phi%i:eta%i>>hTot(%i,%f,%f,%i,%f,%f)",iPi, iPi, nBins[2], min[2], max[2], nBins[1], min[1], max[1]),finalCuts,"colz");
      hTot = (TH2F*)gPad->GetPrimitive("hTot");

      finalCuts += Form(" && tofMatched%i", iPi);
      mTree[kEMBEDING]->Draw(Form("phi%i:eta%i>>hPassed(%i,%f,%f,%i,%f,%f)",iPi, iPi, nBins[2], min[2], max[2], nBins[1], min[1], max[1]),finalCuts,"colz");
      hPassed = (TH2F*)gPad->GetPrimitive("hPassed");
      hPassed->Divide(hTot);
      SetHistStyle(hPassed);
      hPassed->GetZaxis()->SetTitleOffset(1.2);
      hPassed->GetZaxis()->SetRangeUser(0.0, 0.8);
      hPassed->SetTitle(";" + valLabel[2] + ";" + valLabel[1] + ";Eff"); 
      SetPalletRange(hPassed);
      WriteCanvas("PhiEta_TofEff");

      finalCuts = assosiationCuts;
      mTree[kEMBEDING]->Draw(Form("vertexZInCm_TrueMc:eta%i>>hTot(%i,%f,%f,%i,%f,%f)", iPi, nBins[2], min[2], max[2], nBins[3], min[3], max[3]),finalCuts,"colz");
      hTot = (TH2F*)gPad->GetPrimitive("hTot");

      finalCuts += Form(" && tofMatched%i", iPi);
      mTree[kEMBEDING]->Draw(Form("vertexZInCm_TrueMc:eta%i>>hPassed(%i,%f,%f,%i,%f,%f)", iPi, nBins[2], min[2], max[2], nBins[3], min[3], max[3]),finalCuts,"colz");
      hPassed = (TH2F*)gPad->GetPrimitive("hPassed");
      hPassed->Divide(hTot);
      SetHistStyle(hPassed);
      hPassed->GetZaxis()->SetTitleOffset(1.2);
      hPassed->GetZaxis()->SetRangeUser(0.0, 0.8);
      hPassed->SetTitle(";" + valLabel[2] + ";" + valLabel[3] + ";Eff"); 
      drawVertexEtaSpace();
      SetPalletRange(hPassed);
      WriteCanvas("VertexEta_TofEff");

      plotAverageEtaVertexEff(hPassed,"VertexEta_TofEff_" + piName[iPi]);
   }
   canvas->Close();
}//runTofEff


void PlotManager::plotAverageEtaVertexEff(TH2 *hEff, TString histName) 
{
   TH1D *hAverageEff = new TH1D(histName + "_AverageEff",";Eff;",60,0.4,1.0);

   TCanvas *canvasEff;
   CreateCanvas(&canvasEff,histName + "_AverageEff"); 
   SetHistStyle(hAverageEff);
   double eff;

   // Loop through the bins and get the bin content
   for (int i = 1; i <= hEff->GetNbinsX(); ++i) { 
      double eta = hEff->GetXaxis()->GetBinCenter(i);
      if( abs(eta) > maxEta ) continue;
      for (int j = 1; j <= hEff->GetNbinsY(); ++j) {
         double vertex = hEff->GetYaxis()->GetBinCenter(j);
         if( etaVertexSlope*vertex - etaVertexShift > eta || etaVertexSlope*vertex + etaVertexShift < eta )
            continue;
         eff = hEff->GetBinContent(i, j);
         hAverageEff->Fill(eff);
      }
   }

   hAverageEff->SetStats(true);
   hAverageEff->Draw("");

   TFitResultPtr fitResult = hAverageEff->Fit("gaus", "S");

   CreateText(0.1, 0.85, 0.4, 0.95);
   text->AddText( Form("#mu = %0.3f #pm %0.3f",fitResult->Parameter(1),fitResult->Error(1)) );
   text->AddText( Form("#sigma = %0.3f #pm %0.3f",fitResult->Parameter(2),fitResult->Error(2)) );   
   text->Draw("same");

   WriteCanvas(histName + "_AverageEff", canvasEff);
   canvasEff->Close();
}//plotAverageEtaVertexEff


void PlotManager::runEmbedVarDiff()
{
   // plot diff in pT, eta, phi
   changeDir(kEMBEDING, "varDiff");
   CreateCanvas(&canvas,"varDiff");
   canvas->cd();
   canvas->SetLogy();

   TH1F *hist;

   TString piName[] = {mUtil->particleName(embedParticle) + mUtil->signName(PLUS), mUtil->particleName(embedParticle) + mUtil->signName(MINUS)};
   const char* valName[3] = {"pTInGev", "phi", "eta"};
   TString valLabel[3] = {"p_{T} [GeV]", "#phi", "#eta"};
   double range[3] = {1.0,0.1,0.1};

   for (int iPi = 0; iPi < 2; ++iPi)
   {
      changeSubDir(piName[iPi]);
      TString cuts = Form("idTruth%i != -1", iPi); 
      for (int iVal = 0; iVal < 3; ++iVal)
      {
         mTree[kEMBEDING]->Draw(Form("%s%i-%s%i_TrueMc>>h(200,%f,%f)",valName[iVal], iPi,valName[iVal], iPi, -range[iVal], range[iVal]),cuts);
         hist = (TH1F*)gPad->GetPrimitive("h");
         SetHistStyle(hist);
         hist->SetTitle(";#Delta" + valLabel[iVal] + ";counts");
         WriteCanvas(TString(valName[iVal]));
      }
   }
   changeDir(kEMBEDING, "varDiff");
   vector<int> colorSet = { mainColor, bckgColor};
   vector<int> markerStyle = { mainMarker, bckgMarker};

   TString partName[nParticles][nSigns] = { {"#pi^{+}","#pi^{-}"}, {"K^{+}","K^{-}"}, {"p","#bar{p}"}};
   TH1F *hists[nSigns];

   TCanvas *tmpCanvas;
   CreateCanvas(&tmpCanvas,"tmpCanvas");

   CreateLegend(0.84, 0.74, 0.97, 0.89);
   for (int iPi = 0; iPi < nSigns; ++iPi)
   {
      TString cuts = Form("idTruth%i != -1", iPi); 
      tmpCanvas->cd();
      mTree[kEMBEDING]->Draw(Form("TMath::Sqrt( TMath::Power(eta%i_TrueMc-eta%i, 2) + TMath::Power(phi%i_TrueMc-phi%i, 2))>>htemp%i(100,0.0,0.5)",iPi,iPi,iPi,iPi,iPi),cuts);
      
      TH1F* htemp = (TH1F*)gPad->GetPrimitive(Form("htemp%i",iPi));
      hists[iPi] = (TH1F*)htemp->Clone(Form("histPion%i",iPi));

      SetHistStyle(hists[iPi], colorSet[iPi], markerStyle[iPi]);
      hists[iPi]->SetTitle(";#delta(#eta, #phi);counts");
      canvas->cd();
      if( iPi == 0)
         hists[iPi]->Draw("");
      else
         hists[iPi]->Draw("same");
   
      legend->AddEntry(hists[iPi], partName[embedParticle][iPi],"ple");
   }
   legend->Draw("same");
   CreateLine(deltaEtaPhiCut,0,deltaEtaPhiCut,hists[0]->GetMaximum()/2);
   line->SetLineWidth(4);
   line->Draw("same");
   WriteCanvas("deltaEtaPhi"+mUtil->particleName(embedParticle));

   canvas->Close();
   tmpCanvas->Close();
}//runEmbedVarDiff


void PlotManager::runEmbeddingQA()
{
   changeDir(kEMBEDING, "embedQA");
   if( DEBUG )
      cout<<"PlotManager::runEmbeddingQA() is running"<<endl;

   changeSubDir(mUtil->particleName(embedParticle));
   m2FromEmbedding();

   vector<int> tree = { kEMBEDING, kEMBEDING, kMAINANA }; 
   vector<int> colorSet = { 1, bckgColor, mainColor};
   vector<int> markerStyle = { 26, bckgMarker, mainMarker};
   TH1D *hist[nDataSets];

   vector<TString> histList;
   vector<unsigned int> flags; 
   // 0 - TrueMC and Data, vertex scaled; 1 - TrueMC and Data, event info;
   // 2 - all sets, event info; 3 - MCZB and Data, event info;
   // 4 - all sets, tracks info; 5 - MCZB and Data, tracks info

   vector<TString> labelList;
   vector<unsigned int> nBins;
   vector<double> binLow; 
   vector<double> binMax; 

   //AddHist("vertexYInCm","y_{vrtx} [cm]",100,-0.3,0.3,0,histList,labelList,nBins,binLow,binMax,flags);
   //AddHist("vertexXInCm","x_{vrtx} [cm]",100,-0.3,0.3,0,histList,labelList,nBins,binLow,binMax,flags);

   AddHist("vertexZInCm","z_{vrtx} [cm]",100,-100,100,1,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("vertexYInCm","y_{vrtx} [cm]",100,-0.3,0.3,1,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("vertexXInCm","x_{vrtx} [cm]",100,-0.3,0.3,1,histList,labelList,nBins,binLow,binMax,flags);

   int minBin[nParticles] = {-10, -30,-45};
   int maxBin[nParticles] = {25, 25, 25};

   for (int i = 0; i < nSigns; ++i)
      for (int iPart = 0; iPart < nParticles; ++iPart)
         AddHist("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(i),"nSigmaTPC" + mUtil->particleName(iPart) + 
            mUtil->signName(i),100,minBin[iPart],maxBin[iPart],3,histList,labelList,nBins,binLow,binMax,flags);

   AddHist("momentumInGev","p [GeV]",100,0.0,3.0,4,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("pTInGev","p_{T} [GeV]",100,0.0,3.0,4,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("eta","#eta",100,-1.2,1.2,4,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("phi","#phi",100,-3.5,3.5,4,histList,labelList,nBins,binLow,binMax,flags);

   AddHist("dEdxInKevCm","dEdx [keV/cm]",100,1,5,5,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("tofTimeInNs","tofTime [ns]",100,-1000,55000,5,histList,labelList,nBins,binLow,binMax,flags);
   //AddHist("tofLengthInCm","tofLength [cm]",100,-1000,500,5,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("tofLengthInCm","tofLength [cm]",100,150,400,5,histList,labelList,nBins,binLow,binMax,flags);
   //AddHist("dcaXYInCm","DCA_{xy} [cm]",100,0,3.0,5,histList,labelList,nBins,binLow,binMax,flags);
   //AddHist("dcaZInCm","DCA_{z} [cm]",100,0.0,3.0,5,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("nHitsFit","N^{fit}_{hits}",40,10,50,5,histList,labelList,nBins,binLow,binMax,flags);
   AddHist("nHitsDEdx","N^{dE/dx}_{hits}",40,10,50,5,histList,labelList,nBins,binLow,binMax,flags);



   // Draw 1D histos
   TCanvas *tmpCanvas;
   CreateCanvas(&tmpCanvas,"tmpCanvas");

   for (unsigned int iHist = 0; iHist < histList.size(); ++iHist)
   {  
      CreateCanvas(&canvas,histList[iHist]);
      canvas->SetLogy();
      CreateLegend(0.84, 0.74, 0.97, 0.89);

      for (unsigned int set = 0; set < nDataSets; ++set) 
      {
         if(flags[iHist] < 2 && set == MCZB) continue;
         if( (flags[iHist] == 3 || flags[iHist] == 5) && set == MC) continue;

         for (int i = 0; i < ( flags[iHist] > 3 ? nSigns : 1); ++i)
         {            
            TString setTag = set == MC ? "_" + mUtil->dataTagName(TRUEMC) : "";
            TString flagTag = flags[iHist] > 3 ? Form("%i",i) : "";
            TString name = histList[iHist] + flagTag + setTag;

            TString cut = Form("TMath::Abs(eta%i) < %f && pTInGev%i > %f && nHitsDEdx%i >= %i && nHitsFit%i >= %i", i, maxEta, i, minPt[embedParticle], i, minNHitsDEdx[0], i, minNHitsFit[0]);
            if( flags[iHist] <= 3 )
               cut += Form(" && TMath::Abs(eta1) < %f && pTInGev1 > %f && nHitsDEdx1 >= %i && nHitsFit1 >= %i", maxEta, minPt[embedParticle], minNHitsDEdx[0], minNHitsFit[0]);
            
            if( set == DATA){
               cut += Form(" && pairID==%i && TMath::Abs(vertexZInCm) < %f && dcaXYInCm0 < %f && dcaXYInCm1 < %f && TMath::Abs(dcaZInCm0) < %f && TMath::Abs(dcaZInCm1) < %f", embedParticle, vertexRange, maxDcaXY, maxDcaXY, maxDcaZ, maxDcaZ);
            }else if( set == MCZB){
               cut += Form("&& TMath::Abs(vertexZInCm_TrueMc) < %f && TMath::Abs(eta%i_TrueMc) < %f",vertexRange,i,maxEta);
               cut += flags[iHist] > 3 ? Form(" && idTruth%i != -1 && ",i) + assosiationCut(i) : " && idTruth0 != -1 && idTruth1 != -1 && " + assosiationCut(PLUS) + " && " + assosiationCut(MINUS);
            }else{ 
               cut = Form("TMath::Abs(vertexZInCm_TrueMc) < %f && TMath::Abs(eta%i_TrueMc) < %f && ",vertexRange,i,maxEta) + assosiationCut(i);
            }

            // for normalizing vertex dist
            //+= Form(" && TMath::Gaus(%s,0,50,1)/TMath::Gaus(%s,0,96,1) && TMath::Abs(vertexZInCm) < %f", name.Data(), name.Data(), vertexRange);
            tmpCanvas->cd();

            mTree[tree[set]]->Draw(name + Form(">>htemp%i(%i,%.2f,%.2f)",set*100+i,nBins[iHist],binLow[iHist],binMax[iHist]),cut);

            TH1D* htemp = (TH1D*)gPad->GetPrimitive(Form("htemp%i",set*100+i));

            if( i==0 )
               hist[set] = (TH1D*)htemp->Clone("hist" + mUtil->dataSetName(set));
            else
               hist[set]->Add(htemp);
         }

         canvas->cd();
         hist[set]->Scale(1.0/double(hist[set]->Integral()));
         SetHistStyle(hist[set], colorSet[set], markerStyle[set]);
         hist[set]->GetYaxis()->SetTitleOffset(1.4);
         hist[set]->SetTitle(";"+labelList[iHist]+";Normalized counts");

         hist[set]->GetYaxis()->SetRangeUser(1e-4, 10*hist[set]->GetMaximum());

         if( set == 0)
            hist[set]->Draw("hist");
         else
            hist[set]->Draw("same hist");
      
         legend->AddEntry(hist[set], mUtil->dataSetName(set),"ple");
      }
      legend->Draw("same");
      TString canvasName = flags[iHist] > 0 ? histList[iHist] : histList[iHist] + "_norm";
      WriteCanvas(canvasName);
      canvas->Close();
   }

   tmpCanvas->Close();
   //m2FromEmbedding();
}//runEmbeddingQA

void PlotManager::m2FromEmbedding()
{
  if( DEBUG )
      cout<<"PlotManager::m2FromEmbedding() is running"<<endl;

   vector<int> tree = { 0, kEMBEDING, kMAINANA }; 
   vector<int> colorSet = { 1, bckgColor, mainColor};
   vector<int> markerStyle = { 1, bckgMarker, mainMarker};

   TH1D *hist[nDataSets];

   // Draw 1D histos
   TCanvas *tmpCanvas;
   CreateCanvas(&tmpCanvas,"tmpCanvas");

   TString histName = "mSquared";
   CreateCanvas(&canvas,histName);
   canvas->SetLogy();
   CreateLegend();

   for (unsigned int set = MCZB; set < nDataSets; ++set) 
   {
      tmpCanvas->cd();

      TString cut = Form("pairID==%i && TMath::Abs(eta0) < %f && TMath::Abs(eta1) < %f && pTInGev0 > %f && pTInGev1 > %f && nHitsDEdx0 >= %i && nHitsFit0 >= %i && nHitsDEdx1 >= %i && nHitsFit1 >= %i", embedParticle, maxEta, maxEta, minPt[embedParticle], minPt[embedParticle], minNHitsDEdx[0], minNHitsFit[0], minNHitsDEdx[0], minNHitsFit[0]);
      if( set != DATA){
         cut += Form(" && idTruth0 != -1 && idTruth1 != -1 && TMath::Abs(vertexZInCm_TrueMc) < %f && ",vertexRange) + assosiationCut(PLUS) + " && " + assosiationCut(MINUS);
         for (int iPi = 0; iPi < nSigns; ++iPi)
            cut += Form(" && TMath::Abs(eta%i_TrueMc) < %f && eta%i_TrueMc > %f*vertexZInCm_TrueMc - %f && eta%i_TrueMc < %f*vertexZInCm_TrueMc + %f", iPi, maxEta, iPi, etaVertexSlope, etaVertexShift, iPi, etaVertexSlope, etaVertexShift);
      }else{
         cut += Form(" && TMath::Abs(vertexZInCm) < %f && dcaXYInCm0 < %f && dcaXYInCm1 < %f && TMath::Abs(dcaZInCm0) < %f && TMath::Abs(dcaZInCm1) < %f", vertexRange, maxDcaXY, maxDcaXY, maxDcaZ, maxDcaZ);
         for (int iPi = 0; iPi < nSigns; ++iPi)
            cut += Form(" && TMath::Abs(eta%i) < %f && eta%i > %f*vertexZInCm - %f && eta%i < %f*vertexZInCm + %f", iPi, maxEta, iPi, etaVertexSlope, etaVertexShift, iPi, etaVertexSlope, etaVertexShift);
      }

      mTree[tree[set]]->Draw(histName + ">>htempM2(200, -0.5, 1.5)",cut);
      TH1D* htemp = (TH1D*)gPad->GetPrimitive("htempM2");
      hist[set] = (TH1D*)htemp->Clone("hist" + mUtil->dataSetName(set));

      canvas->cd();
      hist[set]->Scale(1.0/double(hist[set]->Integral()));
      SetHistStyle(hist[set], colorSet[set], markerStyle[set]);
      hist[set]->SetTitle(";m^{2}_{TOF} [GeV^{2}];Normalized counts");

      if( set == 0)
         hist[set]->Draw("hist");
      else
         hist[set]->Draw("same hist");
   
      legend->AddEntry(hist[set], mUtil->dataSetName(set),"ple");
   }

   legend->Draw("same");
   WriteCanvas(histName);
   canvas->Close();
   tmpCanvas->Close();
}//m2FromEmbedding

void PlotManager::drawVertexEtaSpace()
{
   double etaIntercept = etaVertexSlope*vertexRange + etaVertexShift;
   double vrtxIntercept = (maxEta - etaVertexShift)/etaVertexSlope;


   CreateLine(-maxEta,vertexRange,etaIntercept,vertexRange);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

   CreateLine(-etaIntercept,-vertexRange,maxEta,-vertexRange);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

   CreateLine(-maxEta,-vrtxIntercept,-maxEta,vertexRange);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

   CreateLine(maxEta,vrtxIntercept,maxEta,-vertexRange);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

   CreateLine(-maxEta,-vrtxIntercept,-etaIntercept,-vertexRange);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

   CreateLine(etaIntercept,vertexRange,maxEta,vrtxIntercept);
   line->SetLineWidth( 2*lineWidth );
   line->Draw("same");

}


void PlotManager::runTpcEmbedSysStudy()
{
   const char* effLabel[nTPCnHitsStudies+1] = {
      "Nominal", "N^{dE/dx}_{hits} loose", "N^{fit}_{hits} loose",
      "N^{dE/dx}_{hits} tight", "N^{fit}_{hits} tight"
   };

   vector<int> colorSet = { 1, 2, 4, 6, 8};
   vector<int> markerStyle = { mainMarker, 24, 25, 26, 28};

   double sigma[nTPCnHitsStudies+1], sigmaErr[nTPCnHitsStudies+1];
   double mean[nTPCnHitsStudies+1], meanErr[nTPCnHitsStudies+1];


   // load TPC plots and calulcate eff hist 
   auto effhist = [&](TString partName, TString histName, TString canvasName)
   {
      CreateCanvas(&canvas,canvasName);
      TH1D *hist[nTPCnHitsStudies+1];
      for (unsigned int iHitStudy = 0; iHitStudy <= nTPCnHitsStudies; ++iHitStudy)
      {
         TString dir = "Embedding/";
         if( iHitStudy > 0 )
            dir += embedSysStudyFiles[iHitStudy - 1] + "/";  
         dir += "TpcEff/" + partName + "/"; 
         TCanvas *cnvs = (TCanvas*)effFile->Get(dir + histName);
         if (!cnvs) {
            std::cerr << "Error in PlotManager::runTpcEmbedSysStudy() while retrieving canvas: "<< histName<<" in study "<<iHitStudy <<std::endl;
            return;
         }
         hist[iHitStudy] = (TH1D*)cnvs->GetPrimitive(histName);
         if (!hist[iHitStudy]) {
            std::cerr << "Error in PlotManager::runTpcEmbedSysStudy() while retrieving histogram for "<< histName<<" in study "<<iHitStudy<< std::endl;
            return;
         }
         SetHistStyle(hist[iHitStudy], colorSet[iHitStudy], markerStyle[iHitStudy]);
         hist[iHitStudy]->SetTitle(";Eff;Counts");
         hist[iHitStudy]->GetXaxis()->SetRangeUser(0.5, 0.95);
         hist[iHitStudy]->GetYaxis()->SetRangeUser(0.0, 120);
         hist[iHitStudy]->SetLineWidth(4);
         mean[iHitStudy] = hist[iHitStudy]->GetMean();
         meanErr[iHitStudy] = hist[iHitStudy]->GetMeanError();
         sigma[iHitStudy] = hist[iHitStudy]->GetStdDev();
         sigmaErr[iHitStudy] = hist[iHitStudy]->GetStdDevError();
         canvas->cd();
         if( iHitStudy == 0)
            hist[iHitStudy]->Draw("hist");
         else
            hist[iHitStudy]->Draw("same hist");
      }
      mCurrDir->cd();
      CreateLegend(0.08,0.6,0.65,0.95);
      for (unsigned int i = 0; i <= nTPCnHitsStudies; ++i)
         legend->AddEntry(hist[i], Form("%s: #mu = %0.3f #pm %0.3f #sigma = %0.3f #pm %0.3f",effLabel[i], mean[i], meanErr[i], sigma[i],sigmaErr[i]), "pl");
      legend->Draw("same");
      WriteCanvas(canvasName);
      canvas->Close();
   };

   TString partName[] = {mUtil->particleName(embedParticle) + mUtil->signName(PLUS), 
      mUtil->particleName(embedParticle) + mUtil->signName(MINUS)};

   effhist(partName[PLUS],"VertexEta_Eff_AverageEff","TPCEff"+partName[PLUS]);
   effhist(partName[MINUS],"VertexEta_Eff_AverageEff","TPCEff"+partName[MINUS]);
}//runTpcEmbedSysStudy




