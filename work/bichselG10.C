//
//root4star [0] .L bichselG10.C 
//root4star [1] .x bichselG10.C("tag") 
//

//________________________________________________________________________________
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include <stdio.h>
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TClassTable.h"
#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "TLegend.h"
#include "TROOT.h"
#else
class Bichsel;
#endif
Bichsel *m_Bichsel = 0;
const Int_t NMasses = 10;
const Double_t Masses[NMasses] = {0.13956995,
               0.493677,
               0.93827231,
               1.87561339,
               0.51099907e-3,
               0.1056584,
               2.80925,
               2.80923, //GEANT3
               3.727417, //GEANT3
               0.13956995,
};
const Int_t   Index[NMasses] = { 2,    3,   4,   5,   0,    1,  6,    7,       8,    -2};
//const Int_t   Index[NMasses] = { 4,    3,   2,   0,   5,    1,  6,    7,       8,    -2};
const Int_t   Colors[NMasses] = { 1,    2,   4,   7,   6,    3,  9,    30,       8,    -2};
const Char_t *Names[NMasses] = {"Pion","Kaon","Proton","Deuteron","Electron","#mu","t","He3","#alpha","2#pi"};
//const Char_t *Names[NMasses] = {"Proton", "Kaon","Pion","Electron", "Deuteron","#mu","t","He3","#alpha","2#pi"};
const Int_t NF = 3;  //         0       1    2     3     4      5   6     7
const Char_t *FNames[8] = {"Girrf","Sirrf","Bz","B70","B60","B70M","dNdx","BzM"};
const Int_t Nlog2dx = 3;
const Double_t log2dx[Nlog2dx] = {0,1,2};
//________________________________________________________________________________
Double_t bichselZ(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZ(TMath::Log10(poverm),dx2)));//TMath::Exp(7.81779499999999961e-01));
   //return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}

void bichselG10(TString input) {  
   if (gClassTable->GetID("StBichsel") < 0 || !m_Bichsel)
   {
      gSystem->Load("libTable");
      gSystem->Load("St_base");
      gSystem->Load("StarClassLibrary");
      gSystem->Load("StBichsel");
      m_Bichsel = Bichsel::Instance();
   }
   //TString inputFileLocation = "/gpfs01/star/pwg/truhlar/Run17_P20ic/" + input + "/28FC38D8C6BCFC2677ACA1B26D01EC8D_985.root";
   TString inputFileLocation = "/gpfs01/star/pwg/truhlar/Run17_P20ic/" + input + "/merged/StRP_production_0000.root";
   TFile* data = TFile::Open(inputFileLocation, "read");
   if (!data)
   {
      cout<<"Error: cannot open "<<inputFileLocation<<endl;
      return;
   }

   TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
   tree->Draw("TMath::Log10(dEdxInKevCm0):TMath::Log10(momentumInGev0)>>hDEdx(200,-1,1,100,0.1,2)","pTMissing < 0.12"); 
   tree->Draw("TMath::Log10(dEdxInKevCm1):TMath::Log10(momentumInGev1)>>+hDEdx","pTMissing < 0.12"); 
   TH2D* hdEdx = (TH2D*)gPad->GetPrimitive("hDEdx");


   TCanvas *cCanvas2D = new TCanvas("hDEdx","hDEdx", 1165.0, 980.0);
   gPad->SetMargin(0.11,0.12,0.125,0.01); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);
   gStyle->SetOptDate(0);
   gStyle->SetLineWidth(4);      //axis line
   gStyle->SetFrameLineWidth(2); //frame line
   gPad->SetTickx();
   gPad->SetTicky(); 
   gStyle->SetOptStat("");
   cCanvas2D->SetGridx(0);
   cCanvas2D->SetGridy(0);
   hdEdx->GetXaxis()->SetTitleFont(43);
   hdEdx->GetYaxis()->SetTitleFont(43);
   hdEdx->GetZaxis()->SetTitleFont(43);
   hdEdx->GetXaxis()->SetLabelFont(43);
   hdEdx->GetYaxis()->SetLabelFont(43);
   hdEdx->GetZaxis()->SetLabelFont(43);
   hdEdx->GetYaxis()->SetRangeUser(0.1,1.7);
   hdEdx->GetXaxis()->SetRangeUser(-0.8,1.0);
   hdEdx->GetXaxis()->SetLabelSize(45);
   hdEdx->GetYaxis()->SetLabelSize(45);
   hdEdx->GetZaxis()->SetLabelSize(45);
   hdEdx->GetXaxis()->SetTitleSize(50);
   hdEdx->GetYaxis()->SetTitleSize(50);
   hdEdx->GetZaxis()->SetTitleSize(50);
   hdEdx->GetXaxis()->SetTitleOffset(1.07);
   hdEdx->GetYaxis()->SetTitleOffset(0.90);
   hdEdx->GetZaxis()->SetTitleOffset(0.50);
   hdEdx->SetStats(0); 
   hdEdx->SetTitle(" ; log_{10} p [GeV] ;log_{10} dE/dx [keV/cm] ");
   hdEdx->Draw("colz");
   cCanvas2D->SetLogz(1);
   
   TPaveText *textPub = new TPaveText(0.15,0.88,0.7,0.95,"brNDC");
   textPub -> SetTextSize(45);
   textPub -> SetFillColor(0);
   textPub -> SetTextFont(43);
   textPub -> SetTextAlign(12);
   textPub -> AddText("p + p #rightarrow p + h^{+}h^{-} + p     #sqrt{s} = 510 GeV");
   textPub -> Draw("same");

   const Char_t *type="Bz";
   TString Type(type);
   TLegend *leg = new TLegend(0.65,0.55,0.8,0.8);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(45);
   leg->SetTextFont(43);
   Double_t xmax = 0.9;
   //  for (int h = 0; h < NMasses; h++) { // Masses
   short lineStyles[] = {1, 2, 5, 6}; 
   for (int h = 0; h < NF; h++) 
   { // Masses
      Int_t f = 1;
      Int_t dx = 1;
      Char_t *FunName = Form("%s%s%i",FNames[f],Names[h],(int)log2dx[dx]);
      cout << "Make " << FunName << endl;
      Double_t xmin = -0.9;
      //    if (h == 0 || h >= 5) xmin = -0.75;
      if (h == 1) xmin = -0.80;
      if (h == 2) xmin = -0.60;
      if (h == 3) xmin = -0.30;
      TF1 *func = new TF1(FunName,bichselZ ,xmin, xmax,2);
      func->SetParameter(0,Masses[h]);
      func->SetParameter(1,1.);

      Int_t color = Colors[h];
      func->SetLineColor(color);
      func->SetMarkerColor(color);
      func->SetLineStyle(lineStyles[h]);
      func->SetLineWidth(4);
      func->Draw("same");
      leg->AddEntry(func,Names[h]);
   }
   leg->Draw("same");

   cCanvas2D->Update();
   cCanvas2D->SaveAs("dEdx.png");
}
