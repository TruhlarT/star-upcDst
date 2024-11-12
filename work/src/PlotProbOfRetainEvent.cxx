#include "PlotManager.h"
/*
// must be here for unknown reassons
std::map< int, double> mInstLumiPerRun; // in [#mub^{-1}s^{-1}]
std::map< int, int> mTimeOfRun; // in [s]
std::map< int, int> mTotal;
std::map< int, int> mPassed[nBranchesConfigurations]; 
*/
vector<unsigned int> runsToSkip = {};//18084050, 18084051, 18085002, 18085038, 18086022, 18087054, 18091004, 18092025, 18092026, 18093026, 18094005, 18094006, 18094013, 18094016, 18094021, 18094030, 18094036, 18094037, 18094050, 18098007, 18098013, 18098019, 18101035, 18105039, 18106054, 18107030, 18112036, 18120044, 18121008, 18121009, 18121017, 18121023, 18122011, 18122016, 18122017, 18122018, 18122019, 18122020, 18122021, 18122022, 18122023, 18127079, 18129044, 18132042, 18134006, 18134009, 18134047, 18135055, 18136013, 18136042, 18136044, 18136045, 18136046, 18136047, 18136048, 18136050, 18136051, 18137004, 18137005, 18137006, 18147011};
//18058012, 18062053, 18062054, 18063034, 18063092, 18063093, 18063094, 18063095, 18063096, 18076016, 18077014, 18083049, 18083053, 18085016, 18086057, 18092030, 18092031, 18092032, 18092033, 18092035, 18092036, 18092038, 18092056, 18092076, 18092087, 18092093, 18093051, 18094009, 18094012, 18094029, 18094033, 18095013, 18106042, 18118022, 18127054, 18127081, 18136026, 18144008, 18119019, 18132061, 18132064, 18132065, 18133001, 18133002, 18133003, 18133004, 18133005, 18133006, 18133007, 18133008, 18133009, 18133010, 18133012, 18133017, 18133041, 18133042, 18133043, 18133044, 18133045, 18133046, 18133047, 18133048, 18133049, 18133050, 18133051, 18134002, 18134003, 18134004, 18134005, 18134012, 18134013, 18134014, 18134015, 18134017, 18134019, 18134034, 18134035, 18134037, 18134038, 18134040, 18134042, 18134043, 18134044, 18134045, 18134046, 18134048, 18134049, 18134050, 18134051, 18135003, 18135004, 18135005, 18135007, 18135008, 18135009, 18135010, 18135011, 18135012, 18135013, 18135014, 18135026, 18135027, 18135028, 18135031, 18135033, 18135034, 18135035, 18135036, 18135037, 18135038, 18135039, 18135050, 18135052, 18135053, 18135054, 18135056, 18135057, 18135059, 18135060, 18135061, 18135062, 18135063, 18135065, 18136001, 18136002, 18136009, 18136010, 18136011, 18136012, 18136014, 18136015, 18136016, 18136017, 18136018, 18136019, 18136020, 18136021, 18136022, 18136023, 18136024, 18136034, 18136036, 18136037, 18136040, 18136041, 18136049, 18136052, 18137008, 18137009, 18137010, 18137011, 18137012, 18137013, 18137014, 18137015, 18137016, 18137017, 18137018, 18137019, 18137027, 18137028, 18137029, 18142041, 18092005, 18094011, 18094015, 18094019, 18094035, 18094051, 18094052, 18094053, 18094054, 18094057, 18094058, 18094059, 18094062, 18094064, 18094065, 18094066, 18094068, 18094069, 18094070, 18094076, 18094077, 18094078, 18094079, 18094080, 18095001, 18095002, 18095003, 18095004, 18095005, 18095006, 18095007, 18095008, 18095009, 18095010, 18095011, 18095012, 18095013, 18095014, 18095017, 18095018, 18095019, 18095022, 18095023, 18095024, 18085009, 18106006, 18115016, 18115018, 18115019, 18083048, 18118018, 18119018, 18128055, 18129037, 18129038, 18134038, 18134043, 18134044, 18134045, 18134046, 18087013, 18087014, 18087015, 18087016, 18087017, 18087018, 18136012, 18087019, 18087020, 18087021, 18087022, 18087023, 18137028, 18138004, 18138005, 18138006, 18138007, 18138009, 18138010, 18087024, 18138012, 18138013, 18139014, 18140009, 18142041, 18088008, 18088009, 18148054, 18148061, 18148062, 18148063, 18148064, 18088010, 18148067, 18148068, 18148069, 18149002, 18149003, 18149004, 18149005, 18149006, 18149008, 18088011, 18149009, 18149012, 18149019, 18149020, 18149021, 18149022, 18088012, 18149025, 18149026, 18149027, 18149028, 18149029, 18149030, 18088013, 18088014, 18088015, 18089064, 18091018, 18092003, 18092004, 18092005, 18092014, 18092036, 18094015, 18094035, 18094066, 18095007, 18095010, 18095013, 18097044, 18097049, 18097051, 18097052, 18097054, 18098006, 18104021, 18085009, 18083045, 18106006, 18106008, 18106036, 18107005, 18108006, 18057084, 18057087, 18057088, 18057089, 18057090, 18057092, 18057106, 18058007, 18058008, 18058016, 18058018, 18058032, 18059008, 18059009, 18059010, 18059013, 18059016, 18059017, 18059018, 18059021, 18059024, 18059026, 18059066, 18059067, 18059068, 18059069, 18059072, 18059077, 18059078, 18060002, 18060003, 18060007, 18060008, 18060013, 18060015, 18060016, 18060119, 18060120, 18061007, 18061008, 18061027, 18061049, 18061051, 18061052, 18061053, 18061055, 18061056, 18061057, 18061058, 18061073, 18061076, 18061079, 18061091, 18061092, 18061093, 18061094, 18061097, 18061098, 18061099, 18061100, 18061101, 18061102, 18061104, 18062001, 18062002, 18062006, 18062007, 18062008, 18062010, 18062011, 18062012, 18062014, 18062015, 18062016, 18062017, 18062018, 18062021, 18062045, 18062046, 18062048, 18062052, 18062055, 18062056, 18062060, 18062064, 18062065, 18062068, 18062069, 18062070, 18063014, 18063015, 18063017, 18063018, 18063019, 18063099, 18063100, 18063104, 18063105, 18063119, 18064007, 18064008, 18065014, 18065015, 18065033, 18065045, 18065047, 18065060, 18065061, 18065062, 18065071, 18065074, 18065076, 18065077, 18065081, 18065082, 18065083, 18066016, 18066020, 18066045, 18067086, 18067087, 18067088, 18067089, 18067095, 18067096, 18067097, 18068011, 18068012, 18068013, 18068014, 18068016, 18069027, 18069028, 18069029, 18069030, 18069033, 18069034, 18069035, 18069036, 18069037, 18069038, 18069040, 18069041, 18069042, 18069043, 18069061, 18069062, 18069064, 18069065, 18069073, 18070001, 18070002, 18070004, 18070005, 18070006, 18070007, 18070026, 18070027, 18070029, 18070030, 18070032, 18070033, 18070034, 18070040, 18070041, 18070057, 18070058, 18070068, 18070073, 18070074, 18070075, 18070076, 18070077, 18070079, 18070080, 18070081, 18070082, 18071001, 18071002, 18071005, 18071006, 18071007, 18071008, 18071009, 18071020, 18071023, 18071025, 18071029, 18071033, 18071034, 18071035, 18071040, 18071041, 18071042, 18071052, 18071077, 18071079, 18071080, 18071081, 18071084, 18071085, 18071086, 18071087, 18072001, 18072002, 18072003, 18072005, 18072009, 18072010, 18072019, 18072020, 18072022, 18072023, 18072024, 18072025, 18072028, 18072029, 18072030, 18072031, 18072054, 18072063, 18073022, 18073023, 18073024, 18073025, 18073026, 18073028, 18073029, 18073031, 18073032, 18073033, 18073034, 18073043, 18073046, 18073057, 18073058, 18073059, 18074001, 18074002, 18074003, 18074004, 18074005, 18074006, 18074007, 18074012, 18074013, 18074014, 18074015, 18074019, 18074020, 18074021, 18074023, 18074024, 18074026, 18074027, 18074029, 18074030, 18074031, 18074032, 18074051, 18074053, 18074055, 18074059, 18074060, 18075001, 18075002, 18075003, 18075004, 18075009, 18075010, 18075011, 18075012, 18075014, 18075015, 18075016, 18075017, 18075018, 18075019, 18075020, 18081079, 18081080, 18081081, 18082003, 18082004, 18082006, 18082009, 18082013, 18082015, 18082017, 18082019, 18082040, 18082043, 18083001, 18083006, 18083014, 18083018, 18083019, 18085056, 18087008, 18089012, 18089055, 18090027, 18090058, 18090059, 18090070, 18091003, 18091008, 18091025, 18091028, 18091040, 18091051, 18091058, 18091059, 18092001, 18092022, 18092065, 18093012, 18093024, 18093025, 18093038, 18098017, 18098034, 18099004, 18099005, 18100009, 18100016, 18100053, 18101005, 18101011, 18101033, 18101034, 18102009, 18102024, 18103016, 18103018, 18103023, 18104018, 18104046, 18105021, 18105029, 18105043, 18105047, 18106013, 18106033, 18106044, 18106045, 18106068, 18107042, 18108029, 18108036, 18108083, 18108084, 18108091, 18109002, 18110020, 18111021, 18111041, 18111047, 18111051, 18112037, 18116006, 18118003, 18119021, 18119030, 18119032, 18119055, 18120005, 18120037, 18123025, 18125028, 18127080, 18127084, 18128024, 18128025, 18128035, 18128051, 18128059, 18131094, 18138008, 18138011, 18139001, 18140013, 18140026, 18143001, 18143002, 18147004, 18147007, 18149010};



// fit paramters get from itterative method
double p0[nBranchesConfigurations] = { 1.408, 1.429, 1.429, 1.412};
//double p0[nBranchesConfigurations] = { 1.403, 1.424, 1.408, 1.422};
double p0err[nBranchesConfigurations] = { 0.012, 0.012, 0.012, 0.012 };
double p1[nBranchesConfigurations] = { 0.01483, 0.01494, 0.01494, 0.01483};
double p1err[nBranchesConfigurations] = { 0.00007, 0.00007, 0.00007, 0.00007};
//double p1[nBranchesConfigurations] = { 0.01479, 0.01490, 0.01480, 0.01487};

double sigmaGaus[nBranchesConfigurations] = { 99, 99, 99, 99};
double nSgmaRng = 4;

void PlotManager::runProbOfRetainEvent()
{
   // study the probability of retaining the event
   changeDir(kFULLZB, "ProbOfRetainEvent");
   fillMaps();
   runProbStudy(true, "_0");
   runProbStudy(true, "final"); 
   //for (const auto& run : runsToSkip)
   for (unsigned int i = 0; i < runsToSkip.size(); ++i)
      cout<<runsToSkip[i]<<endl;
}//runProbOfRetainEvent

void PlotManager::runProbStudy(bool savePlots, TString plotSuffix)
{
   TH1F *hPassed[nBranchesConfigurations], *hTotal;
   for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
      hPassed[iConf] = new TH1F(Form("hPassed_%i",iConf),"hPassed",mTotal.size(),-0.5,mTotal.size()-0.5);
   hTotal = new TH1F("hTotal","hTotal",mTotal.size(),-0.5,mTotal.size()-0.5);
   
   unsigned int nRuns = 1;
   map<UInt_t, double> binToLumi;
   //for (auto const& [key, val] : mTotal)
   for (std::map<int, int>::const_iterator it = mTotal.begin(); it != mTotal.end(); ++it) 
   {
      int key = it->first;
      int val = it->second;

      if( val == 0)
         continue;

      if( isFarAway(key)){
         runsToSkip.push_back(key);
         continue;
      }

      binToLumi[ nRuns ] = mInstLumiPerRun[key];
      for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
         hPassed[iConf]->SetBinContent(nRuns, mPassed[iConf][key]);
      hTotal->SetBinContent(nRuns, val);
      nRuns++;
   }   

   TEfficiency* pEff;
   TGraphAsymmErrors *gr;

   TF1 *expFit = new TF1("expFit", "[0]*(TMath::Exp(-[1]*x))",50,160);
   for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
   {
      expFit->SetParameters( p0[iConf], p1[iConf]);
      expFit->SetParError( 0, p0err[iConf]);
      expFit->SetParError( 1, p1err[iConf]);
      pEff = new TEfficiency(*hPassed[iConf],*hTotal);
      double axis[6][nRuns]; // x, x_up err, x_low err, y, y_up err, y_low err
      for (unsigned int iBin = 1; iBin < nRuns; ++iBin)
      {
         axis[0][iBin - 1] = binToLumi[iBin];
         axis[1][iBin - 1] = 0;
         axis[2][iBin - 1] = 0;
         axis[3][iBin - 1] = pEff->GetEfficiency(iBin);
         axis[4][iBin - 1] = pEff->GetEfficiencyErrorUp(iBin);
         axis[5][iBin - 1] = pEff->GetEfficiencyErrorLow(iBin);
      }
      gr = new TGraphAsymmErrors( nRuns, axis[0], axis[3], axis[1], axis[2], axis[5], axis[4]);
      if( plotSuffix == "final" )
         gr->Fit(expFit,"EMR");
      p0[iConf] = expFit->GetParameter(0);
      p0err[iConf] = expFit->GetParError(0);
      p1[iConf] = expFit->GetParameter(1);
      p1err[iConf] = expFit->GetParError(1);     
      if(savePlots){
         plotProbOfRetainEvent(gr, expFit, mUtil->branchesConfigurationName(iConf) + plotSuffix );
         plotRetainSystStudy( nRuns, axis[0], axis[3], iConf, plotSuffix);
      }
      
      delete pEff;
   }
   for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
      delete hPassed[iConf];
   delete hTotal;                                           
}//runProbStudy

void PlotManager::fillMaps()
{
   for (std::map<int, double>::const_iterator it = mInstLumiPerRun.begin(); it != mInstLumiPerRun.end(); ++it) 
   {
      int key = it->first;
      for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
         mPassed[iConf][key] = 0;
      mTotal[key] = 0; 
   }
   
   for(Long64_t iev=0; iev<mTree[kFULLZB]->GetEntries(); ++iev)
   { //get the event
      mTree[kFULLZB]->GetEntry(iev); 
      if( skipRun(mRecTree[kFULLZB]->getRunNumber()) )
      {
         continue;
      }
      mTotal[mRecTree[kFULLZB]->getRunNumber()]++;
      for (int iBrConf = 0; iBrConf < nBranchesConfigurations; ++iBrConf)
         if(!mRecTree[kFULLZB]->getBbcDsmBit() && !mRecTree[kFULLZB]->getZdcDsmBit() && mRecTree[kFULLZB]->getNVertecies() == 0 && mRecTree[kFULLZB]->getTofMult() <= 8 && noVetoInRp(iBrConf))
            mPassed[iBrConf][mRecTree[kFULLZB]->getRunNumber()]++;
   }   
}//fillMaps

void PlotManager::plotProbOfRetainEvent(TGraphAsymmErrors *gr, TF1 *expFit, TString name)
{
   CreateCanvas(&canvas,name);        
   gr->SetTitle(" ; L [#mub^{-1}s^{-1}]; Prob. of retaining CEP event");
   SetGraphStyle(gr);
   expFit->SetLineColor(mainColor);
   gr->GetYaxis()->SetRangeUser(0.0, 0.80); 
   gr->GetXaxis()->SetRangeUser(50, 180); 
   gr->Draw("ap");
   expFit->Draw("same");

   DrawSTARInternal();
   CreateText(0.7,0.65,0.95,0.89);
   text -> AddText("#xi = a #times e^{-bL}");   
   text -> AddText(Form("#chi^{2}/ndf = %.0f/%i", expFit->GetChisquare(), expFit->GetNDF()));
   text -> AddText(Form("a = %0.3f #pm %0.3f", expFit->GetParameter(0), expFit->GetParError(0)));
   text -> AddText(Form("b = %0.5f #pm %0.5f", expFit->GetParameter(1), expFit->GetParError(1)));
   text -> AddText("RP: "+name(0,5));
   text -> Draw("same");

   WriteCanvas(name);
   canvas->Close();
}//plotProbOfRetainEvent


bool PlotManager::isFarAway( UInt_t runNumber)
{
   double instLum = mInstLumiPerRun[runNumber];
   for (int iConf = 0; iConf < nBranchesConfigurations; ++iConf)
   {
      double eff = mPassed[iConf][runNumber] / double(mTotal[runNumber]);
      double effFit = p0[iConf]*(TMath::Exp(-p1[iConf]*instLum));
      //double minEff = (p0[iConf]-nSgmaRng*p0err[iConf])*(TMath::Exp(-(p1[iConf]+nSgmaRng*p1err[iConf])*instLum));
      //double maxEff = (p0[iConf]+nSgmaRng*p0err[iConf])*(TMath::Exp(-(p1[iConf]-nSgmaRng*p1err[iConf])*instLum));
      double minEff = effFit - nSgmaRng*sigmaGaus[iConf];
      double maxEff = effFit + nSgmaRng*sigmaGaus[iConf];

      if( eff < minEff || eff > maxEff){
         return true;
      }
   }
   return false;
}//isFarAway

bool PlotManager::skipRun( UInt_t runNumber)
{
   //for (const auto& run : runsToSkip) 
   for (unsigned int i = 0; i < runsToSkip.size(); ++i)
   {
      if(runsToSkip[i] == runNumber)
         return true;
      //if( run > runNumber)
      //   return false;
   }
   return false;
}//skipRun

bool PlotManager::noVetoInRp(unsigned int branchConfig)
{
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      int branch = mUtil->branchPerBranchConfiguration( branchConfig, iSide);
      int vetoBranch = mUtil->oppositeBranch(branch);
      for (int iSt = 0; iSt < nStationPerSide; ++iSt)
         if(mRecTree[kFULLZB]->getRpTrigBits(mUtil->rpPerBranchStationOrder(vetoBranch, iSt)))
            return false;
   }
   return true;
}//noVetoInRp

void PlotManager::plotRetainSystStudy(unsigned int nPoints, double *x, double *y, unsigned int iConf, TString plotSuffix)
{
   TString histName = "Diff::" + mUtil->branchesConfigurationName(iConf) + plotSuffix + "projection";
   double modifier = 1;
   if( plotSuffix == "final")
      modifier = 0.1;

   CreateCanvas(&canvas,histName);
   TH1D *hist = new TH1D(histName, "", 100, -0.4*modifier, 0.6*modifier);
   double yDiff[nPoints];
   for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint)
   {
      double fit = p0[iConf]*(TMath::Exp(-p1[iConf]*x[iPoint]));
      yDiff[iPoint] = fit - y[iPoint];
      hist->Fill(yDiff[iPoint]);
   }


   SetHistStyle(hist, mainColor, mainMarker);
   hist->SetTitle(";;Fit - data");
   gPad->SetLogy();
   hist->Fit("gaus");
   TF1 *gaus = (TF1*)hist->GetListOfFunctions()->FindObject("gaus");
   sigmaGaus[iConf] = gaus->GetParameter(2);
   gaus->SetLineColor(mainColor);
   hist->Draw("E");
   DrawSTARInternal();

   CreateText(0.7,0.65,0.95,0.89); 
   text -> AddText(Form("#chi^{2}/ndf = %.0f/%i", gaus->GetChisquare(), gaus->GetNDF()));
   text -> AddText(Form("A = %.0f #pm %.0f", gaus->GetParameter(0), gaus->GetParError(0)));
   text -> AddText(Form("#mu = %0.4f #pm %0.4f", gaus->GetParameter(1), gaus->GetParError(1)));
   text -> AddText(Form("#sigma = %0.4f #pm %0.4f", gaus->GetParameter(2), gaus->GetParError(2)));
   text -> AddText("RP: "+mUtil->branchesConfigurationName(iConf));
   text -> Draw("same");

   WriteCanvas(histName);
   canvas->Close();  


   histName = "Diff::" + mUtil->branchesConfigurationName(iConf) + plotSuffix;
   CreateCanvas(&canvas,histName); 
   TGraph *g = new TGraph(nPoints,x,yDiff);
   g->SetTitle(" ; L [#mub^{-1}s^{-1}]; Fit - data");
   g->SetMarkerStyle(24);
   g->SetMarkerColor(mainColor);
   g->SetLineColor(mainColor);
   g->Draw("AP");
  
   DrawSTARInternal();

   WriteCanvas(histName);
   canvas->Close();  
}//plotRetainSystStudy
