#include "PlotManager.h"

     
void PlotManager::runM2Plots()
{
   vector<TString> stratName = { TString("noM2noPt"), TString("noM2"), TString("StrictDEdx"), TString("FullPID") };
   vector<TString> setName = { TString("BeforeTotCharge"), TString("AfterTotCharge"), TString("Exclusive") };
   vector<TString> partName = { TString("All pairs"), TString("#pi^{+}#pi^{-} (dE/dx)"),
                           TString("K^{+}K^{-} (dE/dx)"), TString("p#bar{p} (dE/dx)") };

   vector<int> colorSet = { 4, 209, 2, 1};
   changeSubDir("m2plots");

   TH1D *hist[ partName.size() ];
   for (unsigned int pid = 0; pid < stratName.size(); ++pid)
   {
      for (unsigned int stat = 0; stat < setName.size(); ++stat)
      {
         CreateCanvas(&canvas, stratName[pid] + "_" + setName[stat]);
         canvas->SetLogy();
         CreateLegend( 0.7, 0.65, 0.97, 0.89);

         // Load and plot hist
         for (unsigned int part = 0; part < partName.size(); ++part)
         {
            hist[part] = (TH1D*)inFile->Get( Form("PID/hMSquared_%i_%i_%i",pid, stat, part) );
            if(!hist[part]){
               cerr<<"Error: cannot loaded in PlotMainAna::runMainAnaPlots()"<<endl; 
               continue;
            } 
            SetHistStyle(hist[part], colorSet[part], 20);
            hist[part]->SetTitle(";m^{2}_{TOF} [GeV^{2}];Number of events");
            if( part == 0)
               hist[part]->Draw("");
            else
               hist[part]->Draw("same");

            legend->AddEntry(hist[part], partName[part],"l");
         }
         legend->Draw("same");

         CreateLine(m2minKaons,0,m2minKaons,hist[0]->GetMaximum()/2);
         line->SetLineStyle(10);
         line->SetLineColor(2);
         line->SetLineWidth(4);
         line->Draw("same");

         CreateLine(m2minProtons,0,m2minProtons,hist[0]->GetMaximum()/2);
         line->SetLineStyle(10);
         line->SetLineColor(1);
         line->SetLineWidth(4);
         line->Draw("same");

         mCurrDir->cd();
         //canvas->Update();
         WriteCanvas(stratName[pid] + "_" + setName[stat]);
         canvas->Close();
      }
   }
}

void PlotManager::drawMissIdProbability() 
{
   const char* histNames[4] = {"hUndefined", "hPions", "hKaons", "hProtons"};
   const char* labels[3][4] = {
     {"#pi^{+}#pi^{-} #rightarrow not identified", "#pi^{+}#pi^{-} #rightarrow #pi^{+}#pi^{-}", "#pi^{+}#pi^{-} #rightarrow K^{+}K^{-}", "#pi^{+}#pi^{-} #rightarrow p#bar{p}"},
     {"K^{+}K^{-} #rightarrow not identified", "K^{+}K^{-} #rightarrow #pi^{+}#pi^{-}", "K^{+}K^{-} #rightarrow K^{+}K^{-}", "K^{+}K^{-} #rightarrow p#bar{p}"},
     {"p#bar{p} #rightarrow not identified", "p#bar{p} #rightarrow #pi^{+}#pi^{-}", "p#bar{p} #rightarrow K^{+}K^{-}", "p#bar{p} #rightarrow p#bar{p}"}
   };
     
   // Number of PADS
   const Int_t Nx = 4; // un/idendefied as
   const Int_t Ny = 3; // true ID
   // Create a canvas and divide it into 3 rows and 4 columns
   outFile->cd();
   CreateCanvas(&canvas, "Probability of miss/identification");
   CanvasPartition(canvas,Nx,Ny,0.03,0.03,0.05,0.01);

   TCanvas* tmpCanvas;
   CreateCanvas(&tmpCanvas, "tmpCanvas"); 

   TPad *pad[Nx][Ny];
   const double textSize = 20; 
   for (Int_t i = 0; i < Ny; i++) {
      tmpCanvas->cd();
      mTree[kEMBEDING] = dynamic_cast<TTree*>( embFile[i]->Get( nameOfTree[kEMBEDING] ) );

      const char* histName = "histTot";
      //const char* cuts = Form("pTInGev0 > %f && pTInGev1 > %f",minPt[i],minPt[i]);
      const char* cuts = Form("treeState == 3 && pTInGev0 > %f && pTInGev1 > %f",minPt[i],minPt[i]);
      const char* binning = "(26,0.2,2.8,26,0.2,2.8)";
      const char* var = "TMath::Min(pTInGev0,pTInGev1):TMath::Max(pTInGev0,pTInGev1)";
      tmpCanvas->cd();
      mTree[kEMBEDING]->Draw(Form("%s>>%s%s", var, histName, binning), cuts, "colz");
      TH2F* histTot = (TH2F*)gPad->GetPrimitive(histName);
      for (Int_t j = 0; j < Nx; j++) {
         histName = histNames[j];

         tmpCanvas->cd();
         mTree[kEMBEDING]->Draw(Form("%s>>%s%s", var, histName, binning), Form("pairID == %i && %s", -1 + j, cuts), "colz");
         TH2F* hist = (TH2F*)gPad->GetPrimitive(histName);
         hist->Divide(histTot);
         SetHistStyle(hist);
         
         if( i == j-1 ){ // true partile == reconstructed one
            hPIDEff[i] = (TH2F*)hist->Clone(Form("pidEff_%s",histNames[j]));
            hPIDEff[i]->SetDirectory(0);
         }
      
         canvas->cd(0);

         // Get the pads previously created.
         pad[j][i] = (TPad*) canvas->FindObject(TString::Format("pad_%d_%d",j,i).Data());
         pad[j][i]->Draw();
         pad[j][i]->SetFillStyle(4000);
         pad[j][i]->SetFrameFillStyle(4000);
         pad[j][i]->cd();

         hist->GetXaxis()->SetLabelSize(textSize);
         hist->GetYaxis()->SetLabelSize(textSize);
         hist->GetXaxis()->SetTitleSize(textSize);
         hist->GetYaxis()->SetTitleSize(textSize);
         hist->GetZaxis()->SetLabelSize(textSize);
         hist->GetZaxis()->SetTitleSize(textSize);
         hist->GetXaxis()->SetTitleOffset(1.2);
         hist->GetYaxis()->SetTitleOffset(1.2);
         hist->SetTitle(";p^{max}_{T} [GeV];p^{min}_{T} [GeV]"); // Remove the default title
         hist->GetZaxis()->SetRangeUser(0.0, 1.0);
         // Size factors
         Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[j][i]->GetAbsWNDC();
         Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[j][i]->GetAbsHNDC();

         hist->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
         hist->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
         if( j < Nx -1){
            hist->GetZaxis()->SetTickLength(0);
            hist->Draw("col");   
         }else
            hist->Draw("colz"); 

         CreateLine(0.2,minPt[i],2.7,minPt[i]);
         line->SetLineColor(1);
         line->SetLineWidth(3);
         line->Draw("same");

         if( i > 0)
         {
            double cut = minPtPair[i]; 
            CreateLine(0.5,cut,2.7,cut);
            line->SetLineColor(1);
            line->SetLineWidth(3);
            line->Draw("same");
         }

         TLatex latex;
         latex.SetTextSize(textSize*1.5); // Make text larger
         latex.SetTextFont(fontStyle); // Set to bold text
         latex.SetNDC(); // Use normalized device coordinates
         latex.DrawLatex(XtoPad(0.1), YtoPad(0.8), labels[i][j]);
      }
   }
   tmpCanvas->Close();
   // Optionally save the canvas to a file        
   changeSubDir("MissID");
   WriteCanvas("missIdProbability");
   canvas->Close();
}

void PlotManager::CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
 
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
 
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
   for (Int_t i=0;i<Nx;i++) {
 
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
 
      for (Int_t j=0;j<Ny;j++) {
 
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
 
         C->cd(0);
 
         auto name = TString::Format("pad_%d_%d",i,j);
         auto pad = (TPad*) C->FindObject(name.Data());
         if (pad) delete pad;
         pad = new TPad(name.Data(),"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
 
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
 
         pad->Draw();
      }
   }
}
 
double PlotManager::XtoPad(double x)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double pw = xu-xl;
   double lm = gPad->GetLeftMargin();
   double rm = gPad->GetRightMargin();
   double fw = pw-pw*lm-pw*rm;
   return (x*fw+pw*lm)/pw;
}
 
double PlotManager::YtoPad(double y)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double ph = yu-yl;
   double tm = gPad->GetTopMargin();
   double bm = gPad->GetBottomMargin();
   double fh = ph-ph*bm-ph*tm;
   return (y*fh+bm*ph)/ph;
}

/*
void bichselG10(TString input) {  
   gSystem->Load("libTable");
   gSystem->Load("St_base");
   gSystem->Load("StarClassLibrary");
   gSystem->Load("StBichsel");
   m_Bichsel = Bichsel::Instance();

   auto bichselZ = [&](Double_t *x,Double_t *par) {
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
      return  TMath::Log10(scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZ(TMath::Log10(poverm),dx2)));
   };

   CreateCanvas(&canvas,"varDiff");
   canvas->cd();
   mTree[kMAINANA]->Draw("TMath::Log10(dEdxInKevCm0):TMath::Log10(momentumInGev0)>>hDEdx(200,-1,1,100,0.1,2)",Form("pTMissing < %f",exclusivityCut),"colz"); 
   mTree[kMAINANA]->Draw("TMath::Log10(dEdxInKevCm1):TMath::Log10(momentumInGev1)>>+hDEdx",Form("pTMissing < %f",exclusivityCut,"colz")); 
   TH2D* hdEdx = (TH2D*)gPad->GetPrimitive("hDEdx");

   SetHistStyle(hdEdx);
   hdEdx->SetTitle(" ; log_{10} p [GeV] ;log_{10} dE/dx [keV/cm] ");
   canvas->SetLogz(1);

   DrawSystemDescription(-1,0.15,0.88,0.7,0.95);
   CreateLegend(0.65,0.55,0.8,0.8);

   TLegend *leg = new TLegend(0.65,0.55,0.8,0.8);

   Double_t xmax = 0.9;
   Double_t xmin = -0.4;
   
   const Char_t *type="Bz";
   TString Type(type);
   const Int_t   Colors[NMasses] = { 1,    2,   4};

   for (int h = 0; h < nParticles; h++) 
   { // Masses
      TF1 *func = new TF1(Form("%s%s%i","Bz",mUtil->particleName(h),1),bichselZ ,xmin, xmax,2);
      func->SetParameter(0,mUtil->mass(h));
      func->SetParameter(1,1.);
      func->SetLineColor(Colors[h]);
      func->SetMarkerColor(Colors[h]);
      func->SetLineStyle(1);
      func->SetLineWidth(4);
      func->Draw("same");
      legend->AddEntry(func,mUtil->particleName(h));
   }
   legend->Draw("same");

   WriteCanvas(TString(valName[iVal]));
}
*/
