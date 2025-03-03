#include "PlotManager.h"

//////// Plotting functions

void PlotManager::CreateCanvas(TCanvas **can, TString canName, double width, double heigth)
{
   *can = new TCanvas(canName,canName,width,heigth); // 1920/2 = 960
   gPad->SetMargin(0.08,0.02,0.11,0.03); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   //gPad->SetLogy(0);

   // setup the drawing style
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetNumberContours(99);
   gStyle->SetLineWidth(lineWidth);      //axis line
   gStyle->SetFrameLineWidth(lineWidth); //frame line
   gStyle->SetPalette(57);
}

void PlotManager::WriteCanvas(TString canName, TCanvas *can)
{
   if( can == nullptr )
      can = canvas;
   if( canName == "")
      canName = can->GetName();
   can->SetName(canName);
   can->Write(canName);
}//WriteCanvas
 

void PlotManager::SetAxisStyle(TAxis *axis)
{
   axis->SetTitleFont(fontStyle);
   axis->SetLabelFont(fontStyle);
   axis->SetLabelSize(txtSize);
   axis->SetTitleSize(tSize);
   axis->SetDecimals(); 
}//SetAxisStyle

void PlotManager::SetAxisStyle(TGaxis *axis)
{
   axis->SetTitleFont(fontStyle);
   axis->SetLabelFont(fontStyle);
   axis->SetLabelSize(txtSize);
   axis->SetTitleSize(tSize);
   axis->SetDecimals(); 
}//SetAxisStyle


void PlotManager::SetHistStyle(TH1* hist, Int_t color, Int_t markStyle)
{
   if (!hist) {
      std::cerr << "Error in PlotManager::SetHistStyle: Invalid histogram pointer!" << std::endl;
      return;
   }
   hist->SetDirectory(0);
   hist->SetStats(false);
   SetAxisStyle(hist->GetXaxis());
   SetAxisStyle(hist->GetYaxis());
   SetAxisStyle(hist->GetZaxis());
   hist->GetXaxis()->SetTitleOffset(1.0);
   hist->GetYaxis()->SetTitleOffset(1.3);
   hist->GetZaxis()->SetTitleOffset(1.0);
   hist->SetLineColor(color);
   hist->SetLineStyle(lineStyle);
   hist->SetLineWidth(lineWidth);  
   hist->SetMarkerSize(markerSize);
   hist->SetMarkerColor(color);
   hist->SetMarkerStyle(markStyle);
   gStyle->SetPalette(57);
}//SetHistStyle

void PlotManager::SetHistBcgStyle(TH1* hist, Int_t color, Int_t markStyle)
{
   if (!hist) {
      std::cerr << "Error in PlotManager::SetHistStyle: Invalid histogram pointer!" << std::endl;
      return;
   }
   hist->SetStats(false);
   hist->SetLineColor(color);
   hist->SetLineStyle(lineStyle);
   hist->SetLineWidth(lineWidth);  
   hist->SetMarkerSize(markerSize);
   hist->SetMarkerColor(color);
   hist->SetMarkerStyle(markStyle);
}//SetHistStyle

void PlotManager::SetGraphStyle(TGraph *graph, Int_t color, Int_t markStyle)
{
   if (!graph) {
      std::cerr << "Error in PlotManager::SetGraphStyle: Invalid Graph pointer!" << std::endl;
      return;
   }
   SetAxisStyle(graph->GetXaxis());
   SetAxisStyle(graph->GetYaxis());
   graph->GetXaxis()->SetTitleOffset(0.9);
   graph->GetYaxis()->SetTitleOffset(1.2);
   graph->SetLineColor(color);
   graph->SetLineStyle(lineStyle);
   graph->SetLineWidth(lineWidth);  
   graph->SetMarkerSize(markerSize);
   graph->SetMarkerColor(color);
   graph->SetMarkerStyle(markStyle);
}//SetGraphStyle


void PlotManager::CreateLegend(double xl, double yl, double xr, double yr)
{
   if(legend)
      delete legend;
   legend = new TLegend(xl, yl, xr, yr);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);
   legend->SetTextSize(txtSize);
   legend->SetTextFont(fontStyle);
   legend->SetMargin(0.2);   
}

void PlotManager::DrawSTARInternal(double xl, double yl, double xr, double yr)
{
   CreateText(xl, yl, xr, yr);
   text -> SetTextFont(capitalStyle);
   text->AddText("STAR");
//   text->AddText("STAR Internal");
   text -> Draw("same");
}

void PlotManager::DrawSTAR(double xl, double yl, double xr, double yr)
{
   CreateText(xl, yl, xr, yr);
   text -> SetTextFont(capitalStyle);
   text->AddText("STAR");
   text -> Draw("same");
}

void PlotManager::CreateText(double xl, double yl, double xr, double yr)
{
   text = new TPaveText(xl, yl, xr, yr,"brNDC");
   text -> SetTextSize(txtSize);
   text -> SetFillColor(0);
   text -> SetFillStyle(0);
   text -> SetBorderSize(0);
   text -> SetTextFont(fontStyle);
   text -> SetTextAlign(32);   
}
 
void PlotManager::SetGPad(double xl, double yl, double xr, double yr)
{
   gPad->SetMargin(xl,yl,xr,yr); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   gPad->SetLogy(0);
   gStyle->SetOptStat("");   
}

void PlotManager::DrawFiducial(int side, bool elasticFV)
{
   if( elasticFV )
   {
      for (int i = 0; i < nRpOrientations; ++i)
      {
         int orientation = (i == 0 ? 1 : -1);
         CreateLine(fpElPxMin,orientation*fpElPyMin,fpElPxMax,orientation*fpElPyMin);
         line->SetLineWidth( 2*lineWidth );
         line->Draw("same");

         CreateLine(fpElPxMin,orientation*fpElPyMax,fpElPxMax,orientation*fpElPyMax);
         line->SetLineWidth( 2*lineWidth );
         line->Draw("same");

         CreateLine(fpElPxMin,orientation*fpElPyMin,fpElPxMin,orientation*fpElPyMax);
         line->SetLineWidth( 2*lineWidth );
         line->Draw("same");  

         CreateLine(fpElPxMax,orientation*fpElPyMin,fpElPxMax,orientation*fpElPyMax);
         line->SetLineWidth( 2*lineWidth );
         line->Draw("same");  
      }
      return;
   }


   for (int i = 0; i < nRpOrientations; ++i)
   {
      int br = nRpOrientations*side + i;
      int orientation = (i == 0 ? 1 : -1);

      double xMin = fpLPRadius[br] < 0.01 ? sqrt(fpPRadius[br] - fpPyMax[br]*fpPyMax[br]) - fpPxCenter[br] : fpLPxMax[br]; // GeV
      double xMax = sqrt(fpPRadius[br] - fpPyMin[br]*fpPyMin[br]) - fpPxCenter[br]; // GeV

      const Int_t n = 100;
      Double_t x[n], y[n];
      for(int i = 0; i < n; ++i){
         x[i] = xMin + ((xMax-xMin)*i)/(n-1);
         y[i] = orientation*sqrt(abs( (x[i] + fpPxCenter[br])*(x[i] + fpPxCenter[br]) - fpPRadius[br] ));
      }

      TGraph* gr = new TGraph(n,x,y);
      gr->SetLineWidth(2*lineWidth);
      gr->Draw("same");

      CreateLine(fpPxMin[br], orientation*fpPyMin[br],xMax, orientation*fpPyMin[br]);
      line->SetLineWidth( 2*lineWidth );
      line->Draw("same");

      CreateLine(fpPxMin[br], orientation*fpPyMax[br],fpPxMin[br],orientation*fpPyMin[br]);
      line->SetLineWidth( 2*lineWidth );
      line->Draw("same");

      if( fpLPRadius[br] < 0.01 ){
         CreateLine(fpPxMin[br], orientation*fpPyMax[br],xMin,orientation*fpPyMax[br]);
         line->SetLineWidth( 2*lineWidth );
         line->Draw("same");
      }else{
         double xMin = fpPxMin[br]; // GeV
         double xMax = fpLPxMax[br]; // GeV
         for(int i = 0; i < n; ++i){
            x[i] = xMin + ((xMax-xMin)*i)/(n-1);
            y[i] = orientation*sqrt(abs( (x[i] + fpLPxCenter[br])*(x[i] + fpLPxCenter[br]) - fpLPRadius[br] ));
         }
         TGraph* gr = new TGraph(n,x,y);
         gr->SetLineWidth(2*lineWidth);
         gr->Draw("same");
      }
   }
}

void PlotManager::CreateLine(double xl, double yl, double xr, double yr)
{
   line = new TLine(xl, yl, xr, yr);
   line->SetLineStyle( lineStyle );
   line->SetLineColor( lineColor );
   line->SetLineWidth( lineWidth );
}//SetLineStyle


void PlotManager::DrawRatioPlot(TH1 *h1, TH1 *h2, TH1 **h3)
{
   TPad *pad1  = new TPad("pad1", "pad1", 0, 0.3, 1.0, 1.0);
   pad1->SetMargin(0.08,0.02,0.0,0.03); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   h1->SetMinimum(1e-3);

   h1->GetYaxis()->SetTitleOffset(1.4);
   h1->Draw();
   h2->Draw("same");
   //gPad->Update();

   canvas->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.02, 1.0, 0.3);
   pad2->SetMargin(0.08,0.02,0.4,0.0); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad
   //canvas->Update();

   (*h3) = (TH1F*)h1->Clone("h3");
   (*h3)->GetYaxis()->SetNdivisions(505);
   SetHistStyle((*h3), mainColor, mainMarker);
   (*h3)->GetYaxis()->SetTitleOffset(1.4);
}//DrawRatioPlot


void PlotManager::DrawMainText(int part, int rpCon) 
{
   DrawSystemDescription(part);
   //DrawForwardProtonKin(rpCon);
   //DrawCentralKin(part);
}//DrawMainText

void PlotManager::DrawSystemDescription(int part, double xl, double yl, double xr, double yr) 
{
   CreateText(xl, yl, xr, yr);
   text -> SetTextAlign(12);
   if( part == -1)
      text -> AddText("         p + p        #sqrt{s} = 510 GeV");
   else
      text -> AddText("         p + p #rightarrow p + " + mUtil->pairLabel(part) + " + p         #sqrt{s} = 510 GeV");
   text -> Draw("same");
}//DrawSystemDescription

void PlotManager::DrawForwardProtonKin(int rpCon, double xl, double yl, double xr, double yr) 
{
   if( rpCon >= nRpConfigurations)
      return;
   CreateText(xl, yl, xr, yr);
   text -> SetTextAlign(12);
   text -> AddText("Forward proton kinematics:   " + mUtil->rpConfigTag(rpCon));
   //text -> AddText(mUtil->rpConfigTag(rpCon));
   //text -> AddText("described in text");
   //text -> AddText(Form("(p_{x} + %.1f)^{2} + p_{y}^{2} < %.2f GeV^{2}",fpPxCenter,fpPRadius));
   //text -> AddText(Form("%.1f GeV < |p_{y}| < %.1f GeV",fpPyMin,fpPyMax));
   //text -> AddText(Form("p_{x} > %.2f GeV",fpPxMin));
   text -> Draw("same");
}//DrawForwardProtonKin

void PlotManager::DrawCentralKin(int part, double xl, double yl, double xr, double yr) 
{
   if( part == PION)
      yl += 0.05;
   CreateText(xl, yl, xr, yr);
   text -> AddText( mUtil->particleLabels(part) + " kinematics:");   
   text -> AddText("described in text");
   //text -> AddText(Form("p_{T} > %.1f GeV",minPt[part]));
   //text -> AddText(Form("|#eta| < %.1f",maxEta));
   if( part != PION)
      text -> AddText(Form("min(p_{T}^{+},p_{T}^{-}) < %.1f",minPtPair[part]));
   text -> Draw("same");
}//DrawCentralKin



///////// Useful functions

void PlotManager::changeDir(int study, TString dirName)
{
   mDir = mStudyDir[study]->GetDirectory(dirName);
   if (mDir) {// Directory exists, change to it
      mDir->cd();
   } else {// Directory doesn't exist, create it
      mStudyDir[study]->mkdir(dirName);
      mDir = mStudyDir[study]->GetDirectory(dirName);
      mDir->cd();
   }
   mCurrDir = mDir;
}

void PlotManager::changeSubDir(TString dirName)
{
   mCurrDir = mDir->GetDirectory(dirName);
   if (mCurrDir) {// Directory exists, change to it
      mCurrDir->cd();
   } else {// Directory doesn't exist, create it
      mDir->mkdir(dirName);
      mCurrDir = mDir->GetDirectory(dirName);
      mCurrDir->cd();
   }
}

// Function to add values of a given function to each bin of a histogram
void PlotManager::AddFunctionToHistogram(TH1* histogram, TF1* function) 
{
   if (!histogram || !function) {
      std::cerr << "PlotManager::AddFunctionToHistogram: Invalid histogram or function pointer." << std::endl;
      return;
   }

   double xmin, xmax;
   function->GetRange(xmin,xmax);

   // Iterate over each bin of the histogram
   for (Int_t iBin = 1; iBin <= histogram->GetNbinsX(); ++iBin) {
      // Get the center of the bin
      Double_t binCenter = histogram->GetBinCenter(iBin);
      if( binCenter > xmax)
         break;
      if( binCenter < xmin)
         continue;

      // Evaluate the function at the bin center
      Double_t functionValue = function->Eval(binCenter);

      // Add the function value to the bin content
      histogram->SetBinContent(iBin, histogram->GetBinContent(iBin) + functionValue);
   }
}


void PlotManager::replaceNegativeBinsWithAbsolute(TH1* histogram) {
   if (!histogram)
      return;

   // Loop over all bins in the histogram
   // If bin content is negative, replace it with its absolute value
   for (Int_t i = 1; i <= histogram->GetNbinsX(); ++i) 
      if (histogram->GetBinContent(i) < 0)
         histogram->SetBinContent(i, TMath::Abs(histogram->GetBinContent(i)));
}

void PlotManager::AddHist(TString name, TString xlabel, unsigned int nBins, double low, double max, unsigned int flag, 
   vector<TString> &histList, vector<TString> &labelList, vector<unsigned int> &nBinsVec, vector<double> &binLow, 
   vector<double> &binMax, vector<unsigned int> &flags)
{
   labelList.push_back(xlabel);
   histList.push_back(name);
   nBinsVec.push_back(nBins);
   binLow.push_back(low);
   binMax.push_back(max);
   flags.push_back(flag);
}//AddHist

void PlotManager::SetPalletRange(TH1* hist, double min, double max)
{
   TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
   if (palette) {
     palette->SetX1NDC(min); // Set the left edge of the palette
     palette->SetX2NDC(max); // Set the right edge of the palette
   }
}

std::vector<float> PlotManager::getBinsVectorF(const TAxis *axis) const
{
   std::vector<float> binsVector;
   for( int i=1; i<=axis->GetNbins(); ++i)
      binsVector.push_back( axis->GetBinLowEdge( i ) );
   binsVector.push_back( axis->GetBinUpEdge( axis->GetNbins() ) );
   return binsVector;
}


std::vector<double> PlotManager::getBinsVectorD(const TAxis *axis) const
{
   std::vector<double> binsVector;
   for( int i=1; i<=axis->GetNbins(); ++i)
      binsVector.push_back( axis->GetBinLowEdge( i ) );
   binsVector.push_back( axis->GetBinUpEdge( axis->GetNbins() ) );
   return binsVector;
}


double* PlotManager::setBinArray(int nBins, double min, double max)
{
    double *dataArray = new double[nBins];
    // Calculate the bin width
    double binWidth = (max - min) / (nBins - 1);
    for (int i = 0; i < nBins; ++i){
        dataArray[i] = min + i * binWidth;
        //cout<< array[i]<<" ";
    }
    //cout<<endl;
    return dataArray;
}

unsigned int PlotManager::getRpCombination()
{
    return getRpDeltaPhi() < 90 ? IT : ET;
}//getRpCombination

double PlotManager::getRpDeltaPhi()
{
    double deltaPhi = TMath::Abs(mCurrentTree->getPhiRp(East) - mCurrentTree->getPhiRp(West))*convertToDegree;
    if(deltaPhi > 180)
        deltaPhi = 360 - deltaPhi;

    return deltaPhi;
}//getRpDeltaPhi

void PlotManager::HighlightBin(TH1 *hist, int bin)
{
    // Create a box to highlight the bin (use TBox for a rectangle)
    TBox *box = new TBox(hist->GetBinCenter(bin) - hist->GetBinWidth(bin)/2, 1, hist->GetBinCenter(bin) + hist->GetBinWidth(bin)/2, hist->GetBinContent(bin));
    box->SetFillColor(4);  // Set the color to red
    box->SetFillStyle(3001);  // Optional: Set a pattern (3004 is cross-hatching)
    
    // Draw the box over the bin
    box->Draw("same");
}//HighlightBin

void PlotManager::DrawLine(TH1 *hist, double x, bool negativeRange)
{
   double yMax = canvas->GetLogy() ? hist->GetMaximum()*0.1 : hist->GetMaximum()*0.7;

   CreateLine(x,1,x,yMax);
   line->Draw("same");
   if(negativeRange){
      CreateLine(-x,1,-x,yMax);
      line->Draw("same");
   }
}//DrawLine
