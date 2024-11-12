void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);
double XtoPad(double x);
double YtoPad(double x);



void drawMissIdProbability() {
   const double textSize = 20;   
   const int fontStyle = 43;

const TString embedPathPrefix = "/gpfs01/star/pwg/truhlar/Run17_P20ic/";
const TString embedPathSuffix = "/merged/StRP_production_0000.root";
const TString embedFile[] = { TString("piEmbed"), TString("kEmbed"), TString("pEmbed") };
   // Open the ROOT files
   //const char* fileNames[3] = {"piEmbedFinal.root", "kEmbedFinal.root", "pEmbedFinal.root"}; 
   const char* fileNames[3] = {"piEmbed.root", "kEmbed.root", "pEmbed.root"}; 
   // kkEmbedTest StRP_production_0000   kkEmbed
   const char* histNames[4] = {"hUndefined", "hPions", "hKaons", "hProtons"};
   const char* labels[3][4] = {
     {"#pi^{+}#pi^{-} #rightarrow not identified", "#pi^{+}#pi^{-} #rightarrow #pi^{+}#pi^{-}", "#pi^{+}#pi^{-} #rightarrow K^{+}K^{-}", "#pi^{+}#pi^{-} #rightarrow p#bar{p}"},
     {"K^{+}K^{-} #rightarrow not identified", "K^{+}K^{-} #rightarrow #pi^{+}#pi^{-}", "K^{+}K^{-} #rightarrow K^{+}K^{-}", "K^{+}K^{-} #rightarrow p#bar{p}"},
     {"p#bar{p} #rightarrow not identified", "p#bar{p} #rightarrow #pi^{+}#pi^{-}", "p#bar{p} #rightarrow K^{+}K^{-}", "p#bar{p} #rightarrow p#bar{p}"}
   };
     
   // Create a canvas and divide it into 3 rows and 4 columns
   TCanvas* canvas = new TCanvas("canvas", "Probability of miss/identification", 1800, 1200);
   //canvas->Divide(4, 3, 0.01, 0.01);

   // Number of PADS
   const Int_t Nx = 4; // un/idendefied as
   const Int_t Ny = 3; // true ID

   const double minPt[] = {0.25, 0.3, 0.4}; // for pion, kaon and proton in GeV

   // Margins
   Float_t lMargin = 0.03;
   Float_t rMargin = 0.03;
   Float_t bMargin = 0.05;
   Float_t tMargin = 0.01;

   CanvasPartition(canvas,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
   TCanvas* tmpCanvas = new TCanvas("tmpCanvas", "", 1800, 1200);

   TPad *pad[Nx][Ny];

   for (Int_t i = 0; i < Ny; i++) {
      TFile* file = new TFile( embedPathPrefix + embedFile[i] + embedPathSuffix);
      TTree* tree = dynamic_cast<TTree*>(file->Get("embedQATree"));

      tmpCanvas->cd();
      const char* histName = "histTot";
      //const char* cuts = Form("pTInGev0 > %f && pTInGev1 > %f",minPt[i],minPt[i]);
      const char* cuts = Form("treeState == 3 && pTInGev0 > %f && pTInGev1 > %f",minPt[i],minPt[i]);
      const char* binning = "(58,-0.1,2.8,58,-0.1,2.8)";
      const char* var = "TMath::Min(pTInGev0,pTInGev1):TMath::Max(pTInGev0,pTInGev1)";
      tree->Draw(Form("%s>>%s%s", var, histName, binning), cuts, "colz");
      TH1F* histTot = (TH1F*)gPad->GetPrimitive(histName);
      for (Int_t j = 0; j < Nx; j++) {
         histName = histNames[j];

         tmpCanvas->cd();
         tree->Draw(Form("%s>>%s%s", var, histName, binning), Form("pairID == %i && %s", -1 + j, cuts), "colz");
         TH1F* hist = (TH1F*)gPad->GetPrimitive(histName);
         hist->Divide(histTot);

         canvas->cd(0);

         // Get the pads previously created.
         pad[j][i] = (TPad*) canvas->FindObject(TString::Format("pad_%d_%d",j,i).Data());
         pad[j][i]->Draw();
         pad[j][i]->SetFillStyle(4000);
         pad[j][i]->SetFrameFillStyle(4000);
         pad[j][i]->cd();

                  // Size factors
         Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[j][i]->GetAbsWNDC();
         Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[j][i]->GetAbsHNDC();

         hist->GetXaxis()->SetTitleFont(fontStyle);
         hist->GetYaxis()->SetTitleFont(fontStyle);
         hist->GetXaxis()->SetLabelFont(fontStyle);
         hist->GetYaxis()->SetLabelFont(fontStyle);
         hist->GetXaxis()->SetLabelSize(textSize);
         hist->GetYaxis()->SetLabelSize(textSize);
         hist->GetXaxis()->SetTitleSize(textSize);
         hist->GetYaxis()->SetTitleSize(textSize);
         hist->GetXaxis()->SetTitleOffset(1.2);
         hist->GetYaxis()->SetTitleOffset(1.2);
         hist->SetTitle(";p^{max}_{T} [GeV];p^{min}_{T} [GeV]"); // Remove the default title
         hist->GetZaxis()->SetRangeUser(0.0, 1.0);
         hist->SetStats(false);
         hist->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
         hist->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
         if( j < Nx -1)
            hist->GetZaxis()->SetTickLength(0);
         hist->Draw("colz");

         TLine *line = new TLine(0.2,minPt[i],2.7,minPt[i]);
         line->SetLineColor(1);
         line->SetLineWidth(3);
         line->Draw("same");

         if( i > 0)
         {
            double cut = i == 1 ? 0.6 : 1.0; 
            TLine *l = new TLine(0.5,cut,2.7,cut);
            l->SetLineColor(1);
            l->SetLineWidth(3);
            l->Draw("same");
         }

         TLatex latex;
         latex.SetTextSize(textSize); // Make text larger
         latex.SetTextFont(fontStyle); // Set to bold text
         latex.SetNDC(); // Use normalized device coordinates
         latex.DrawLatex(XtoPad(0.1), YtoPad(0.8), labels[i][j]);
      }
   }

    /*
    // Loop over the files
    for (int i = 0; i < 3; ++i) { 
        TFile* file = new TFile(fileNames[i]);
        TTree* tree = dynamic_cast<TTree*>(file->Get("embedQATree"));

        tmpCanvas->cd();
        const char* histName = "histTot";
        tree->Draw(Form("TMath::Min(pTInGev0,pTInGev1):TMath::Max(pTInGev0,pTInGev1)>>%s(29,-0.1,2.8,29,-0.1,2.8)", histName), "", "colz");
        TH1F* histTot = (TH1F*)gPad->GetPrimitive(histName);

        for (int iPart = 0; iPart < 4; ++iPart) { 
            histName = histNames[iPart];

            tmpCanvas->cd();
            tree->Draw(Form("TMath::Min(pTInGev0,pTInGev1):TMath::Max(pTInGev0,pTInGev1)>>%s(29,-0.1,2.8,29,-0.1,2.8)", histName), Form("pairID == %i", -1 + iPart), "colz");
            TH1F* hist = (TH1F*)gPad->GetPrimitive(histName);
            hist->Divide(histTot);
            hist->GetXaxis()->SetTitleFont(fontStyle);
            hist->GetXaxis()->SetTitleFont(fontStyle);
            hist->GetXaxis()->SetLabelFont(fontStyle);
            hist->GetYaxis()->SetLabelFont(fontStyle);
            hist->GetXaxis()->SetLabelSize(labelSize);
            hist->GetYaxis()->SetLabelSize(labelSize);
            hist->GetXaxis()->SetTitleSize(labelSize);
            hist->GetYaxis()->SetTitleSize(labelSize);
            hist->GetXaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitleOffset(0.8);
            hist->SetTitle(";p^{max}_{T} [GeV];p^{min}_{T} [GeV]"); // Remove the default title
            hist->GetZaxis()->SetRangeUser(0.0, 1.0);
            hist->SetStats(false);
            
            canvas->cd(i * 4 + iPart + 1);
            hist->Draw("colz");
            
            // Draw label
            TLatex latex;
            latex.SetTextSize(0.05); // Make text larger
            latex.SetTextFont(42); // Set to bold text
            latex.SetNDC(); // Use normalized device coordinates
            latex.DrawLatex(0.15, 0.85, labels[i][iPart]); // Position at the top left corner
        }
    }
    */
    tmpCanvas->Close();
    // Optionally save the canvas to a file
    canvas->SaveAs("missIdProbability.png");
    //canvas->Close();
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
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
 
double XtoPad(double x)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double pw = xu-xl;
   double lm = gPad->GetLeftMargin();
   double rm = gPad->GetRightMargin();
   double fw = pw-pw*lm-pw*rm;
   return (x*fw+pw*lm)/pw;
}
 
double YtoPad(double y)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double ph = yu-yl;
   double tm = gPad->GetTopMargin();
   double bm = gPad->GetBottomMargin();
   double fh = ph-ph*bm-ph*tm;
   return (y*fh+bm*ph)/ph;
}