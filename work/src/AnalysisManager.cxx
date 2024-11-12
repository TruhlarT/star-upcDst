// Run by: ./AnalysisManager file.list
// e.g. ./AnalysisManager /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list
// or you can open just root file
// ./AnalysisManager /star/data01/pwg_tasks/upc02/Part9/18143045/18143045.root
// or you can open n-th root file in file.list
// ./AnalysisManager file.list index 

#include "UpcDstLibreries.h"
#include "RunDef.h"
#include "AnalysisManager.h"


//_____________________________________________________________________________
int main(int argc, char** argv) 
{
   cout<<"Starting the analysis..."<<endl;
   std::cout << "C++ version: " << __cplusplus << std::endl;

   //open output file
   outFile = CreateOutputFile("AnalysisOutput.root"); 
   if(!outFile) 
   {
      cout << "Can not open output file." << endl; 
      return 2;
   }

   runRP = false;
   if( runMAINANA ){
      mAnaVector.push_back(new MainAna(outFile));
      runRP = true;
   }
   if( runVERTEXSTUDY )
      mAnaVector.push_back(new VertexStudy(outFile));
   if( runEMBEDDING )
      mAnaVector.push_back(new Embedding(outFile));
   if( runTOFQA )
      mAnaVector.push_back(new TofQA(outFile));
   if( runTOFEFF ){
      mAnaVector.push_back(new TofEff(outFile));
      runRP = true;
   }
   if( runTRIGEFF )
      mAnaVector.push_back(new TrigEff(outFile));
   if( runFULLZB )
      mAnaVector.push_back(new FullZB(outFile));
   if( runELASTICANA ){
      mAnaVector.push_back(new ElasticAna(outFile));
      runRP = true;
   }
   if( runRPMCANA ){
      mRpMCAna = new RpMCAna(outFile);
      runRP = true;
   }

   //connect input file(s)
   if(!ConnectInput(argc, argv))
   {
      cout << "Wrong input parameters..." << endl; 
      return 1;
   }

   cout<<"Output file created..."<<endl;

   if( runRP )
   {     
      // Load RP off-sets with off-set corrections
      if( !LoadOffsetFile(nameOfOffSetFile, mOffSet) )
         return 3;
      #if !defined ALIGNMENT
         if( !LoadOffsetFile(nameOfOffSetCorrectionFile, mCorrection) )
            return 3;
      #endif
      cout<<"RP offset loaded..."<<endl;
   }
   // Initiate histograms
   Init();

   //ask for number of events
   nEvents = upcTree->GetEntries();
   //nEvents = 1000; // use for debugging and testing
   cout<<"Proccesing "<<nEvents<<" events"<<endl;
   //event loop through RP MC 
   if( runRPMCANA )
      setMcEventsPerZbEvent();

   //event loop through data 
   for(Long64_t iev=0; iev<nEvents; ++iev) 
   { //get the event
      if( DEBUG && iev%1000 == 0)
         cout<<"Analyzing "<<iev<<". event "<<endl;
      upcTree->GetEntry(iev);
      mRunNumber = upcEvt->getRunNumber();
      if( runRP )
         SetRpEvent();

      Make();
      if( runRPMCANA && ZbEvents.size() > 50)
         RunRpMcAna(iev);

      if( runRP )
         ReleaseTheMemoryOfCorrRpEvt();

   }
   if( runALIGNMENT ){
      RunAlignment();
      SaveAlignment(); 
   }


   //close the outputs
   CleanMemory();
   outFile->Write(0, TObject::kOverwrite);
   outFile->Close();

   cout<<"Ending Analysis... GOOD BYE!"<<endl;
   return 0;
}//main

void Make()
{
   // Set RP correction
   TVector3 corr[nRomanPots];
   TVector3 offset[nRomanPots];
   for (int i = 0; i < nRomanPots; ++i){
      corr[i] = mCorrection[i][mRunNumber];
      offset[i] = mOffSet[i][mRunNumber];
   }

   for (unsigned int i = 0; i < mAnaVector.size(); ++i)
   {
      mAnaVector[i]->SetRpPosition(corr, offset); 
      mAnaVector[i]->SetEvent(upcEvt, rpEvt, mcEvt);
      mAnaVector[i]->Make();
   }
}


void Init()
{
   mUtil = new Util();
   for (unsigned int i = 0; i < mAnaVector.size(); ++i){
      mAnaVector[i]->Init();
   }

   if( runRPMCANA )
      mRpMCAna->Init();
}

void RunRpMcAna(Long64_t iev)
{
   if( iZbEv >= ZbEvents.size())
      return;

   if( ZbEvents[iZbEv] != iev )
      return;

   iZbEv++;


   // Set RP correction
   TVector3 corr[nRomanPots];
   TVector3 offset[nRomanPots];
   for (int i = 0; i < nRomanPots; ++i){
      corr[i] = mCorrection[i][mRunNumber];
      offset[i] = mOffSet[i][mRunNumber];
   }
   mRpMCAna->SetRpPosition(corr, offset); 

   unsigned int nMcEventsPerZbEvent;

   nMcEventsPerZbEvent = nMcEvents/ ZbEvents.size(); 

   for (unsigned int iEmbed = 0; iEmbed < nMcEventsPerZbEvent; ++iEmbed)
   {
      if(iMcEvnt == nMcEvents)
         iMcEvnt = 0;

      mcTree->GetEntry(iMcEvnt++);
      mRpMCAna->SetEvent(upcEvt, rpEvt, mcEvt);
      mRpMCAna->SetMCInfo( mc_vrtx, mc_p);
      if(!mRpMCAna->AreTracksInRp())
         continue;

      mRpMCAna->Make();
   }

}

void RunAlignment()
{
/* The offset corrections are set at zero at the beginning. 
   Thus, the afterburner for the first iteration keeps tracks unchanged.
   This function is designed to be run on single run i.e. 1 job for 1 run
*/
   ElasticAna *aligAna =  new ElasticAna(outFile);
   aligAna->Init();
   for (unsigned int iter = 0; iter < nAligIteration; ++iter)
   {
      for(Long64_t iev=0; iev<nEvents; ++iev) 
      { //get the event
         upcTree->GetEntry(iev);
         mRunNumber = upcEvt->getRunNumber();
         for (unsigned int iRP = 0; iRP < nRomanPots; ++iRP){
            mCorrection[iRP][mRunNumber][Z] = 0;
         }
         SetRpEvent(); 
         aligAna->SetEvent(upcEvt, rpEvt, mcEvt);
         aligAna->Make();
         ReleaseTheMemoryOfCorrRpEvt();
      }          
      for (unsigned int iRP = 0; iRP < nRomanPots; ++iRP){
         mCorrection[iRP][mRunNumber] += aligAna->CalculateOffsetCorrForRP(iRP, iter);
      }

   }
   delete aligAna;   
}

void SaveAlignment()
{

   TTree *aligTree = new TTree("aligTree", "aligTree");

   aligTree->Branch("runNumber", &mRunNumber);

   // RP event info
   for (int i = 0; i < nRomanPots; ++i)
   {
      aligTree->Branch("events_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][Z]);
      aligTree->Branch("X_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][X]);
      aligTree->Branch("Y_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][Y]);
   }

   aligTree->Fill();
}


bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots])
{
   string line;
   ifstream file( fileName );
   if (!file.is_open() )
   {
      cout << "\n ERROR in LoadOffsets(): Problems with opening a file: "<< fileName<< endl;
      return false;
   }   
   int RUNNUMBER;
   double var[2*nRomanPots]; // placeholder for offsets or offsets correction in x and y
   while ( getline(file, line) )
   {
      stringstream ss(line);
      ss >> RUNNUMBER;
      for (int i = 0; i < 2*nRomanPots; ++i)
         ss >> var[i]; 

      for (int iRP = 0; iRP < nRomanPots; ++iRP)
         offsets[iRP][RUNNUMBER].SetXYZ(var[2*iRP], var[2*iRP + 1], 0.0); 
   }

   file.close();
   return true;
}

void runAfterburner(StRPEvent *event, StRPEvent *newRpEvent)
{
   for(unsigned int iCl = 0; iCl < event->getNumberOfClusters(); ++iCl){
      StUPCRpsCluster *oldCl = event->getCluster(iCl);
      float position = oldCl->position();
      int RP = oldCl->romanPotId();
      int plane = oldCl->planeId();
      // Set new position
      position -= mCorrection[RP][mRunNumber][mUtil->planeToCoor(plane)];
      // Create new TrackPoint with new position
      StUPCRpsCluster *rpCluster = newRpEvent->addCluster();
      rpCluster->setPosition(position);
      rpCluster->setPositionRMS(oldCl->positionRMS()); 
      rpCluster->setLength(oldCl->length()); 
      rpCluster->setEnergy(oldCl->energy()); 
      rpCluster->setXY(oldCl->xy()); 
      rpCluster->setQuality(oldCl->quality()); 
      rpCluster->setPlaneId(plane);
      rpCluster->setRomanPotId(RP);
   }

   for(unsigned int iTP = 0; iTP < event->getNumberOfTrackPoints(); ++iTP){
      StUPCRpsTrackPoint *oldTP = event->getTrackPoint(iTP);
      TVector3 position = oldTP->positionVec();
      // Set new position
      int RP = oldTP->rpId();
      position -= mCorrection[RP][mRunNumber];
      // Create new TrackPoint with new position
      StUPCRpsTrackPoint *rpTrackPoint = newRpEvent->addTrackPoint();
      rpTrackPoint->setPosition(position);
      rpTrackPoint->setRpId(RP);
      rpTrackPoint->setTime(oldTP->time(0), 0); // pmtId = 0
      rpTrackPoint->setTime(oldTP->time(1), 1); // pmtId = 1
      rpTrackPoint->setQuality(oldTP->quality());    
      for(UInt_t iPlane=0; iPlane < nPlanes; ++iPlane) // 4 planes in a RP
         rpTrackPoint->setClusterId(oldTP->clusterId(iPlane), iPlane);
   }

   for(unsigned int k = 0; k < event->getNumberOfTracks(); ++k)
   {// Get pointer to k-th track in Roman Pot data collection
      StUPCRpsTrack *track = event->getTrack(k);
      track->setEvent(event);

      if( !(track->getTrackPoint(RP1) ? track->getTrackPoint(RP1)->planesUsed()>=nPlanesUsed : false) ||
      !(track->getTrackPoint(RP2) ? track->getTrackPoint(RP2)->planesUsed()>=nPlanesUsed : false))
         continue;

      StUPCRpsTrack *rpTrack = newRpEvent->addTrack();
      rpTrack->setEvent(newRpEvent);
      rpTrack->setFirstTrackPointId(track->getFirstTrackPointId());
      rpTrack->setSecondTrackPointId(track->getSecondTrackPointId());

      StUPCRpsTrackPoint *newTP = track->getTrackPoint(RP1);
      TVector3 position = newTP->positionVec();

      rpTrack->setType(track->type());
      rpTrack->setBranch(track->branch());
      rpTrack->setP( CalculateMomentumVector(position.X(), position.Y(), position.Z(), rpTrack)); // setting the momentum vector
   }
}

TVector3 CalculateMomentumVector(double TPx, double TPy, double TPz, StUPCRpsTrack *track)
{
    /* DX bending angles in the East and West; 0=E, 1=W */
   double beamMomenta[2] = { mUtil->p0() , mUtil->p0()};

   int sign = (TPz < 0 ? -1 : 1 );
   int iSide = (TPz < 0 ? E : W );
   // below calculating momentum vector

   if( track->type() == StUPCRpsTrack::rpsLocal)
   {
      double x_BCS = TPx - mXYZ_IP[X] - sin(mThetaXY_tilt[X])*( TPz - mXYZ_IP[2] ); // x_RP in beam coordinate system
      double y_BCS = TPy - mXYZ_IP[Y] - sin(mThetaXY_tilt[Y])*( TPz - mXYZ_IP[2] ); // y_RP in beam coordinate system
      double localThetaX = x_BCS / abs( TPz);
      double localThetaY = y_BCS / abs( TPz);
      double momentumValue = beamMomenta[iSide];

      TVector3 momentumVector( 0, 0, sign*momentumValue );
      //cout<<"Angle: "<<localThetaY<<endl;
      momentumVector.RotateX( -sign*localThetaY );
      momentumVector.RotateY( sign*localThetaX );
      return momentumVector;
   }

   double localThetaX = track->thetaRp( 0 ) - sign*mThetaXY_tilt[0]; // adding tilt of the beam (set to 0)
   double localThetaY = track->thetaRp( 1 ) - sign*mThetaXY_tilt[1]; 
   double x_BCS =  TPx - mXYZ_IP[0] - sin(mThetaXY_tilt[0])*( TPz - mXYZ_IP[2] ); // x_RP1 in beam coordinate system
   double d2 = abs( TPz ) - mLDX[iSide] - mDistanceFromIPtoDX[iSide]; // distance from DX magnet exit to first RP station
   double thetaX_IP = ( x_BCS - (d2 + 0.5*mLDX[iSide])*localThetaX ) / ( mDistanceFromIPtoDX[iSide]  + 0.5*mLDX[iSide] );
   double xi = 1. / ( 1 + (mBendingAngle[iSide]*(mDistanceFromIPtoDX[iSide] + 0.5*mLDX[iSide])) / ( localThetaX*abs( TPz ) - x_BCS ) );
   double momentumValue = beamMomenta[iSide] * (1.-xi);

   TVector3 momentumVector( 0, 0, sign*momentumValue );
   //cout<<"Angle: "<<localThetaY<<endl;
   momentumVector.RotateX( -sign*localThetaY );
   momentumVector.RotateY( sign*thetaX_IP );
   return momentumVector;
   
}

void SetRpEvent()
{
   correctedRpEvent = new StRPEvent(*origRpEvt); 
   if(AFTERBURNER)
   {
      correctedRpEvent->clearEvent();
      runAfterburner(origRpEvt, correctedRpEvent); 
      rpEvt = correctedRpEvent;
   }else
   {
      rpEvt = origRpEvt;
   }
}

void ReleaseTheMemoryOfCorrRpEvt()
{
   delete correctedRpEvent; 
   correctedRpEvent=0;
}

void CleanMemory()
{
   // Clean up the memory
   for (Ana* ana : mAnaVector)
      delete ana;
   mAnaVector.clear(); // Optional: Remove the pointers from the vector
   if( runRPMCANA )
      delete mRpMCAna;
}

bool ConnectInput(int argc, char** argv) 
{
   int fileId = -1;
   upcTree = 0x0;

   if( argc < 2 || argc > 3)
      return false;

   const string& input = argv[1];
   if(input.find(".root") != string::npos && argc == 2)
   {
      cout << "Input from root file: "<< input << endl;
      inFile = TFile::Open(input.c_str(), "read");
      if(!inFile)
      {
         cout<< "Couldn't open input root file..."<<endl;
         return false;
      } 
      upcTree = dynamic_cast<TTree*>( inFile->Get("mUPCTree") );

      if( runRPMCANA ){
         // open Monte Cralo      
         string mcInput = pathToMC + input.substr(input.length()-13, input.length()); 
         cout << "MC input from root file: "<< mcInput << endl;
         mcFile = TFile::Open(mcInput.c_str(), "read");
         if(!mcFile)
         {
            cout<< "Couldn't open input MC root file..."<<endl;
            return false;
         } 
         mcTree = dynamic_cast<TTree*>( mcFile->Get("picoDST") );    
      }
   }else if(input.find(".list") != string::npos )
   {
      cout << "Input from chain" << endl;
      upcChain = new TChain("mUPCTree");
      mcChain = new TChain("picoDST");
      ifstream instr(input.c_str());
      if (!instr.is_open())
      {
         cout<< "Couldn't open: "<<input.c_str()<<endl;
         return false;
      }

      string line;
      int lineId=0;
      if(argc == 3)
         fileId = atoi(argv[2]);

      while(getline(instr, line)) 
      {
         if(fileId==lineId || fileId== -1)
         {
            upcChain->AddFile(line.c_str());
            inFile = TFile::Open(line.c_str(), "read");
            if(!inFile)
            {
               cout<< "Couldn't open: "<<line.c_str()<<endl;
               return false;
            } 
         }
         if( runRPMCANA ){
            string mcInput = pathToMC + line.substr(line.length()-14, line.length()); 
            mcChain->AddFile(mcInput.c_str());
            inFile = TFile::Open(mcInput.c_str(), "read");
            if(!inFile)
            {
               cout<< "Couldn't open: "<<mcInput.c_str()<<endl;
               return false;
            } 
         }
         lineId++;
      }
      instr.close();
      upcTree = dynamic_cast<TTree*>( upcChain );
      mcTree = dynamic_cast<TTree*>( mcChain );
   }
   
   if(!upcTree) 
      return false;

   upcTree->SetBranchAddress("mUPCEvent", &upcEvt);
   if( runRP )
      upcTree->SetBranchAddress("mRPEvent", &origRpEvt); 
   
   if( runRPMCANA ){
      mcTree->SetBranchAddress("StRPEvent", &mcEvt); 
      mcTree->SetBranchAddress("mX_IP", &mc_vrtx[0]); 
      mcTree->SetBranchAddress("mY_IP", &mc_vrtx[1]); 
      mcTree->SetBranchAddress("mZ_IP", &mc_vrtx[2]); 
      mcTree->SetBranchAddress("mPx_RPE", &mc_p[0][E]); 
      mcTree->SetBranchAddress("mPy_RPE", &mc_p[1][E]); 
      mcTree->SetBranchAddress("mPz_RPE", &mc_p[2][E]); 
      mcTree->SetBranchAddress("mPx_RPW", &mc_p[0][W]); 
      mcTree->SetBranchAddress("mPy_RPW", &mc_p[1][W]); 
      mcTree->SetBranchAddress("mPz_RPW", &mc_p[2][W]); 
   }

   return true;
}//ConnectInput


//_____________________________________________________________________________
TFile *CreateOutputFile(const string& out) {

   TFile *outputFile = TFile::Open(out.c_str(), "recreate");
   if(!outputFile) 
      return 0x0;

   return outputFile;
}//CreateOutputFile


void setMcEventsPerZbEvent()
{
   nMcEvents = mcTree->GetEntries();
   //nMcEvents = 10000; // use for debugging and testing
   iMcEvnt = 0;
   iZbEv = 0;

   if( DEBUG )
      cout<<"Proccesing setMcEventsPerZbEvent:"<<endl;
   for(Long64_t iev=0; iev<nEvents; ++iev)
   {
      if( DEBUG && iev%10000 == 0)
         cout<<"    "<<iev<<". event "<<endl;

      upcTree->GetEntry(iev);
      SetRpEvent();

      if( mRpMCAna->CheckTriggers(&ZBtriggers, upcEvt, nullptr) ) 
         if( mRpMCAna->Veto(upcEvt, rpEvt))
            ZbEvents.push_back(iev);   

      ReleaseTheMemoryOfCorrRpEvt();
   }
   mRunNumber = upcEvt->getRunNumber();
   mRpMCAna->SetRunNumber(mRunNumber);
   cout<<"Proccesing "<<nMcEvents<<" RP MC events and "<<ZbEvents.size()<<" ZB events"<<endl;
}//setMcEventsPerZbEvent
