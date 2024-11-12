#include "PlotManager.h"

//_____________________________________________________________________________
int main(int argc, char** argv)
{
   cout<<"Starting the PlotManager..."<<endl;

   PlotManager *pm = new PlotManager();

   if( !pm->Init(argc, argv) ){
      cerr<<"PlotManager::main() could not init input files..."<<endl;
      return 1;
   }
   if( runVERTEXSTUDY )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runVertexStudy()"<<endl;
      pm->runVertexStudy();
   }
   if( runMAINANA )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runMainAnaPlots()"<<endl;
      pm->runMainAnaPlots();
   }
   if( runEMBEDDING )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runEmbeddingStudy()"<<endl;
      for (int i = 0; i < nParticles; ++i)
      {
         pm->initEmbeddingStudy(i);
         pm->runEmbeddingStudy();
         if( runMAINANA )
            pm->runEmbeddingQA();
      }
   }

   if( runTOFQA )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runTofQA()"<<endl;
      pm->runTofQA();
   }

   if( runTRIGEFF )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runDsmEffStudy()"<<endl;
      pm->runDsmEffStudy( kTRIGEFF );
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runTofTrigStudy()"<<endl;
      pm->runTofTrigStudy();
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runRpTrigStudy()"<<endl;
      pm->runRpTrigStudy();
   }

   if( runFULLZB ){
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runProbOfRetainEvent()"<<endl;
      pm->runProbOfRetainEvent();
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runDsmEffStudy()"<<endl;
      pm->runDsmEffStudy( kFULLZB );
   }

   if( runRPMCANA )
   {
      if( DEBUG )
         cerr<<"PlotManager::main() going to run runRPMCPlots()"<<endl;
      pm->runRPMCPlots();
   }
   if( DEBUG )
      cerr<<"PlotManager::main() going to delete PlotManager"<<endl;
   delete pm;

   cout<<"PlotManager finished succesfully!"<<endl;
   cout<<"Check FinalPlots.root for the plots!"<<endl;

}// end of main

PlotManager::PlotManager()
{
   mUtil = new Util();
   for (int iPart = 0; iPart < nParticles; ++iPart)
      for (int iRpCon = 0; iRpCon < nRpConfigurations +1; ++iRpCon)
         mBcgFraction[iPart][iRpCon] = -1.0;
}

PlotManager::~PlotManager()
{
   if( inFile )   inFile->Close();
   if( outFile )   outFile->Close();
   if( effFile )   effFile->Close();
   if( graniittiFile )   graniittiFile->Close();
   for (int i = 0; i < nParticles; ++i)
      if( embFile[i] )   
         embFile[i]->Close();
   if(mUtil) delete mUtil;
}

bool PlotManager::ConnectInputFile( TString inputFile, TFile **file) 
{
   cout << "Input from root file: "<< inputFile << endl;
   *file = TFile::Open(inputFile, "read");
   if(!*file)
   {
      cerr<< "PlotManager::ConnectInputFile() Couldn't open "<< inputFile <<endl;
      return false;
   } 

   return true;
}//ConnectInput


//_____________________________________________________________________________
TFile* PlotManager::CreateOutputFile(const string& out) 
{

   TFile *outputFile = TFile::Open(out.c_str(), "recreate");
   if(!outputFile) 
      return 0x0;

   return outputFile;
}//CreateOutputFile

bool PlotManager::Init(int argc, char** argv)
{
   //connect input file
   TString inputAnaFile;
   //inputAnaFile = "/gpfs01/star/pwg/truhlar/Run17_P20ic/ana0304c/merged/StRP_production_0000.root";
   TString inputEmbedFile[nParticles];
   for (int i = 0; i < nParticles; ++i)
      inputEmbedFile[i] = embedPathPrefix + embedFile[i] + embedPathSuffix;

   if(argc > 2)
      return false; // too many arguments


   const string& input = argv[1];
   if( input.find(".root") == string::npos)
   {
      cerr<<"Wrong input argument"<<endl;
      return false;
   }
   inputAnaFile = input.c_str();


   if(!ConnectInputFile(inputAnaFile, &inFile))
   {
      cerr << "PlotManager::Init() Could not open file: "<< inputAnaFile <<endl; 
      return false;
   }

   if(!ConnectInputFile(graniittiFilePath, &graniittiFile))
   {
      cerr << "PlotManager::Init() Could not open file: "<< graniittiFilePath <<endl; 
      return false;
   }

   //open output file
   outFile = CreateOutputFile("FinalPlots.root"); 
   if(!outFile) 
   {
      cerr << "Can not open output file." << endl; 
      return false;
   }

   if(!readLumiFile()) 
      return false;

   if(!loadEfficiencies()) 
      return false;

   // open embedding files
   for (int iPart = 0; iPart < nParticles; ++iPart)
   {
      if(!ConnectInputFile(inputEmbedFile[iPart], &embFile[iPart]))
      {
         cerr << "PlotManager::Init() Could not open file: "<< inputEmbedFile[iPart] <<endl; 
         return false;
      }
   }
      if(runStudy[kEMBEDING])
         mStudyDir[kEMBEDING] = outFile->mkdir(studyName[kEMBEDING]);

   for (int iStudy = 0; iStudy < nStudies; ++iStudy)
   {
      mRecTree[iStudy] = nullptr;
      if( !runStudy[iStudy] || nameOfTree[iStudy]=="" || iStudy == kEMBEDING)
         continue;

      mTree[iStudy] = dynamic_cast<TTree*>( inFile->Get( nameOfTree[iStudy] ) );
      if (!mTree[iStudy])
      {
         cerr<<"Error: cannot open "<<nameOfTree[iStudy]<<endl;
         return false;
      }
      mRecTree[iStudy] = new RecTree(mTree[iStudy], treeBits[iStudy] );
      
      mStudyDir[iStudy] = outFile->mkdir(studyName[iStudy]);
   }   
   cout<<"All studies loaded..."<<endl;

   setTextSize();
   return true;
}


bool PlotManager::readLumiFile()
{
   ifstream lumiFile;
   lumiFile.open( pathToLumiFile );
   if (!lumiFile.is_open() )
   {
      cerr << "\nERROR in PlotManager::readLumiFile(): Problems with opening a file: " + pathToLumiFile << endl;
      return false;
   }   

   int runNumber, time0, time1, fillNumber;
   double lumi, prescale, livetime;
   unsigned int timeDiff;
   string line;
   if( DEBUG )
      cerr<<"PlotManager::readLumiFile() going to readi lumi file: "<<pathToLumiFile<<endl;
   while ( getline(lumiFile, line) )
   {
      stringstream ss(line);
      ss >> runNumber >> time0 >> time1 >> fillNumber >> lumi >> prescale >> livetime;
      timeDiff = time1-time0;
      double instantinousLumi = prescale*lumi*1000000/double(livetime*timeDiff);
      if( instantinousLumi < 60 || instantinousLumi > 160 ){
         //cout<<"Skipping run "<<runNumber<<" luminosity out of range"<<endl;
         //cout<<"prescale: "<<prescale<<" lumi: "<<lumi<<" livetime: "<<livetime<<" timeDiff: "<<timeDiff<<endl;
         continue;
      }
      mInstLumiPerRun[runNumber] = instantinousLumi;
      mTimeOfRun[runNumber] = timeDiff;
   }
   if( DEBUG )
      cerr<<"File "<<pathToLumiFile<<" was read"<<endl;
   lumiFile.close();

   lumiFile.open( pathToIntegLumiFile );
   if (!lumiFile.is_open() )
   {
      cerr << "\nERROR in PlotManager::readLumiFile(): Problems with opening a file: " + pathToIntegLumiFile << endl;
      return false;
   }   

   while ( getline(lumiFile, line) )
   {
      stringstream ss(line);
      ss >> runNumber >> lumi;
      mIntegLumiPerRun[runNumber] = lumi*1e+3;
      mRuns.push_back(runNumber);
   }

   mRuns.push_back( mRuns[mRuns.size()-1] + 1);
   lumiFile.close();

   return true;
}


bool PlotManager::loadEfficiencies()
{
   if(!ConnectInputFile(pathToEffRootFile, &effFile))
   {
      cerr << "PlotManager::loadEfficiencies() Could not open file: "<< pathToEffRootFile <<endl; 
      return false;
   }

   // load vertex cut corrections
   double eff, dataEff;
   //double mean, sigma, entries;
   int fill;

   TTree *vertexTree = dynamic_cast<TTree*>( effFile->Get( studyName[kVERTEXSTUDY] + "/vertexTree" ) ); 
   if(!vertexTree)
   {
      cerr << "PlotManager::loadEfficiencies() Could not open vertexTree..."<<endl; 
      return false;
   }
   vertexTree->SetBranchAddress("fill", &fill);
   //vertexTree->SetBranchAddress("sigma", &sigma);   
   //vertexTree->SetBranchAddress("mean", &mean);
   vertexTree->SetBranchAddress("eff", &eff);
   vertexTree->SetBranchAddress("dataEff", &dataEff);
   //vertexTree->SetBranchAddress("entries", &entries);
 
   for (Long64_t i = 0; i < vertexTree->GetEntries(); i++) {
      vertexTree->GetEntry(i);
      mVertexEffPerFill[fill] = 0.5*(eff + dataEff);
      mVertexEffErrorPerFill[fill] = abs(eff - dataEff);
   }
   // load TPC plots and calulcate eff hist 
   auto effhist = [&](TString dir, TString val, int part)
   {
      // Retrieve the canvas from the file
      TH1D* hist = nullptr;
      
      TCanvas *canvas = (TCanvas*)effFile->Get(dir + "/" + val);
      if (!canvas) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first canvas!" << std::endl;
         return hist;
      }

      TH1D *hPassed = (TH1D*)canvas->GetPrimitive("h");
      if (!hPassed) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first histogram!" << std::endl;
         return hist;
      }

      canvas = (TCanvas*)effFile->Get(dir + "/" + val  + "_TrueMc");
      if (!canvas) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving second canvas!" << std::endl;
         return hist;
      }

      TH1D *hTotal = (TH1D*)canvas->GetPrimitive("h2");
      if (!hTotal) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving second histogram!" << std::endl;
         return hist;
      }

      hPassed->Divide(hTotal);
      return hPassed;
   };


   // load TPC eff 2D plots 
   auto effhist2D = [&](TString dir, TString val, int part)
   {
      // Retrieve the canvas from the file
      TH2D* hist = nullptr;
      TCanvas *canvas = (TCanvas*)effFile->Get(dir + "/" + val  + "_Eff");
      if (!canvas) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first canvas!" << std::endl;
         return hist;
      }

      TH2D *hEff = (TH2D*)canvas->GetPrimitive("hPassed");
      if (!hEff) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first histogram!" << std::endl;
         return hist;
      }

      return hEff;
   };

   // load TPC or TOF average eff
   auto getEff = [&](TString partName, int iStudy, bool TOFEff = false)
   {
      double mean = 0.0;
      TH1D *tmpHist;
      TString dir = "Embedding/";
      if( iStudy > 0 )
         dir += embedSysStudyFiles[iStudy - 1] + "/";  
      if( TOFEff )
         dir += "TofEff/" + partName + "/"; 
      else
         dir += "TpcEff/" + partName + "/"; 

      TString histName = TOFEff ? "VertexEta_TofEff_" + partName + "_AverageEff" : "VertexEta_Eff_AverageEff";

      TCanvas *cnvs = (TCanvas*)effFile->Get(dir + histName);
      if (!cnvs) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving canvas: "<< dir + histName<<" in study "<<iStudy <<std::endl;
         return mean;
      }
      tmpHist = (TH1D*)cnvs->GetPrimitive(histName);
      if (!tmpHist) {
         std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving histogram for "<< histName<<" in study "<<iStudy<< std::endl;
         return mean;
      }
      mean = tmpHist->GetMean();
      return mean;
   };


   for (unsigned int iPart = 0; iPart < nParticles; ++iPart)
   {
      for (unsigned int iCharge = 0; iCharge < nSigns; ++iCharge)
      {
         hTPCPhiEff[iPart][iCharge] = effhist("Embedding/TpcEff/" + mUtil->particleName(iPart) + mUtil->signName(iCharge), "phi", iCharge);
         hTPCEtaZEff[iPart][iCharge] = effhist2D("Embedding/TpcEff/" + mUtil->particleName(iPart) + mUtil->signName(iCharge), "VertexEta", iCharge);
         hTPCPtEff[iPart][iCharge] = effhist("Embedding/TpcEff/" + mUtil->particleName(iPart) + mUtil->signName(iCharge), "pTInGev", iCharge);
         hTPCEtaPhiEff[iPart][iCharge] = effhist2D("Embedding/TpcEff/" + mUtil->particleName(iPart) + mUtil->signName(iCharge), "PhiEta", iCharge);

         for (unsigned int i = 0; i < nTPCnHitsStudies+1; ++i){
            tpcEff[i][2*iPart + iCharge] = getEff( mUtil->particleName(iPart) + mUtil->signName(iCharge), i);
            //cout<<"TPC eff for "<< mUtil->particleName(iPart) + mUtil->signName(iCharge) << " and study "<<i<<" is set for: "<<tpcEff[i][2*iPart + iCharge]<<endl;
         }

         tofEff[2*iPart + iCharge] = getEff( mUtil->particleName(iPart) + mUtil->signName(iCharge), 0, true);
      }
   }
   return true;
}