#include "PlotManager.h"

void PlotManager::runTofQA()
{
   outFile->cd();
   int runNumber;
   double eta;
   TTree *tofQATree = new TTree("tofQATree", "tofQATree");
   tofQATree->Branch("runNumber", &runNumber);
   tofQATree->Branch("eta", &eta);

   inFile->cd("TofQA");
   TIter histkeyList(gDirectory->GetListOfKeys());
   TKey *key;
   while ((key = (TKey*)histkeyList()))
   {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH2")) 
         continue;
      TH2 *h = (TH2*)key->ReadObj();
      TString histName = h->GetName(); 
      if( !histName.Contains("hEtaPhi_") )
         continue;
      runNumber = (histName.Remove(0,8)).Atoi(); // remove "hEtaPhi_" a convert the rest to integer
      if( runNumber == 999)
         continue;

      eta = h->ProjectionX()->GetMean();
      if( eta < -0.051 || eta > 0.096) // +-9 sigma
         cout<<runNumber<<" "<<eta<<endl;
      tofQATree->Fill();


   }
   outFile->cd();
   tofQATree->Write();
}


