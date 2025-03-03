#include "PlotManager.h"

void PlotManager::runElasticStudy()
{
   changeDir(kELASTICANA, "vertexStudy");
   if( DEBUG )
      cerr<<"PlotManager::runElasticStudy() going to run vertexStudy()"<<endl;
   vertexStudy();
   changeDir(kELASTICANA, "elasticStudy");
   if( DEBUG )
      cerr<<"PlotManager::runElasticStudy() going to run elasticStudy()"<<endl;
   elasticStudy();
}   


void PlotManager::elasticStudy()
{
   // load histogram from root file, edit and draw them
   auto drawHist = [&](TString histName, TString title, int ID = -1, bool setLogy = false, TString textInput = "")
   {
      TString dir = studyName[kELASTICANA];
      TString histFullName = dir + "/" + histName;

      bool twoDim = ID == -1 ? true : false;

      TH1 *hist = nullptr;
      inFile->GetObject(histFullName, hist);

      if (!hist) {
         cerr<<"Error: cannot loaded histogram "<<histName<<"in PlotElastic::elasticStudy()"<<endl; 
         return;
      }

      CreateCanvas(&canvas, histName);
      if(twoDim)
         gPad->SetMargin(0.10,0.1,0.105,0.03); // 2D Margin
      else if(ID == 8)
         gPad->SetMargin(0.08,0.01,0.105,0.01);
      else
         gPad->SetMargin(0.07,0.04,0.105,0.05);

      if(setLogy)
         canvas->SetLogy();
      SetHistStyle(hist);
      hist->SetTitle(title);

      TString opt = twoDim ? "colz" : "";
      if(twoDim){
         hist->GetZaxis()->SetTitleOffset(1.6);
         hist->GetYaxis()->SetTitleOffset(2.0);
         canvas->SetLogz();
      }
      else
         hist->GetYaxis()->SetRangeUser(1, setLogy ? 100*hist->GetMaximum() : 1.2*hist->GetMaximum());

      if( ID == 8){
         gStyle->SetPaintTextFormat("4.1f");
         hist->GetYaxis()->SetTitleOffset(1.5);
         hist->SetLineColor(1);
         hist->SetMarkerColor(1);
         hist->LabelsOption("d");
         opt = "hist TEXT30";
      }

      hist->Draw(opt);

      DrawSTARInternal();
      DrawSystemDescription();



      if( ID == 1){
         TF1 *fit = new TF1("f1", "gaus(0)");
         fit->SetParameter(0,7.8e+4);
         fit->SetParameter(1,0.0);
         fit->SetParameter(2,1.30e-4);
         hist->Fit(fit,"Q");

         CreateText(0.7,0.65,0.95,0.89); 
         text -> AddText(Form("#sigma = %.2f #pm %.2f [#murad]", fit->GetParameter(2)*pow(10,6), fit->GetParError(2)*pow(10,6)));
         text -> Draw("same");

         CreateLegend(0.7, 0.55, 0.95, 0.65);
         legend->AddEntry(hist, "Data","ple");
         legend->AddEntry(fit, "Gaus","l");
         legend->Draw("same");

      }
      else if( ID == -1){
         // Create and customize the circle
         TEllipse* circle = new TEllipse(0.0, 0.0, nSigma*sigma);
         circle->SetLineColor(kRed);    // Set circle line color
         circle->SetLineWidth(2);       // Set circle line width
         circle->SetFillStyle(0);       // Make the circle transparent
         circle->Draw("same");          // Draw the circle on top of the histogram
      }

      WriteCanvas(histName);
      canvas->Close();

   };

   drawHist("hDeltaTheta", ";#Delta#theta [rad];Number of events",1);
   drawHist("hThetaCorr", ";#Delta#theta_{X} [rad];#Delta#theta_{Y} [rad]");
   drawHist("AnaFlow_"+studyName[kELASTICANA], ";;Number of events",8, true);
}

void PlotManager::vertexStudy()
{

   TH1D *hVertex[2][2]; // [X or Y] [EU-WD or ED-WU]
   hVertex[0][0] = new TH1D("hVertex_X_EU-WD", "hVertex_X_EU-WD", 120, -2, 2);
   hVertex[0][1] = new TH1D("hVertex_X_ED-WU", "hVertex_X_ED-WU", 120, -2, 2);
   hVertex[1][0] = new TH1D("hVertex_Y_EU-WD", "hVertex_Y_EU-WD", 120, -1.5, 2.5);
   hVertex[1][1] = new TH1D("hVertex_Y_ED-WU", "hVertex_Y_ED-WU", 120, -1.5, 2.5);

   TF1 *fit = new TF1("f", "pol1");

   TIter next(inFile->GetDirectory(studyName[kELASTICANA])->GetListOfKeys());
   TKey* key;
   while (( key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());

      if (!cl->InheritsFrom("TH2D")) 
         continue;
      TString name = key->GetName();

      if( !name.BeginsWith("hVertexExtraction") )
         continue;

      TSubString rpArm = name(18,5);
      TSubString coordinate = name(24,1);
      TSubString run = name(26,8);

      int arm = rpArm == "EU-WD" ? 0 : 1;
      int cor = coordinate == "x" ? 0 : 1;

      TH2D *hist = (TH2D*)key->ReadObj();

      hist->Fit(fit,"Q");
      hVertex[cor][arm]->Fill( fit->GetParameter(0) );

      if( run!="18097010")
         continue;

      CreateCanvas(&canvas,name ); 
      gPad->SetMargin(0.07,0.11,0.105,0.03); // 2D Margin
      SetHistStyle(hist);
      TString title = ";Local angle in ";
      title+= cor == X ? "X-Z plane #theta_{X}^{RP} [mrad];X^{RP} [mm]" : "Y-Z plane #theta_{Y}^{RP} [mrad];X^{RP} [mm]";
      hist->SetTitle(title); 
      hist->Draw("colz");

      CreateText(0.1,0.65,0.3,0.89); 
      text -> SetTextAlign(12); 
      text -> AddText("run 18097010");
      text -> AddText(mUtil->armName(arm));
      text -> AddText(Form("%s = %0.3f #pm %0.3f", cor == X ? "<X_{vtx}>" : "<Y_{vtx}>",fit->GetParameter(0), fit->GetParError(0)));
      text -> Draw("same");

      WriteCanvas();
   }


   auto plotHits = [&](TH1D* hVertex, char coordinate, char rpArm) 
   {
      hVertex->Fit("gaus","Q");
      TString hName = hVertex->GetName();
      CreateCanvas(&canvas,hName ); 
      SetHistStyle(hVertex);
      TString title = ";"; // title
      title+= coordinate == X ? "<X_{vtx}>" : "<Y_{vtx}>"; // x-axis title
      title+= ";runs"; // y-axis title
      hVertex->SetTitle(title); 
      hVertex->Draw();

      CreateText(0.7,0.65,0.95,0.89); 
      TF1 *gaus = (TF1*)hVertex->GetListOfFunctions()->FindObject("gaus");
      text -> AddText(Form("#mu = %0.3f #pm %0.3f", gaus->GetParameter(1), gaus->GetParError(1)));
      text -> AddText(Form("#sigma = %0.3f #pm %0.3f", gaus->GetParameter(2), gaus->GetParError(2)));
      text -> AddText(mUtil->armName(rpArm));
      text -> Draw("same");

      WriteCanvas();

   };


   for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
         plotHits(hVertex[i][j], i, j);

}