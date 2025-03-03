#include "Embedding.h"
#include "StPicoPhysicalHelix.h"
#include "RunDef.h"

Embedding::Embedding(TFile *outFile): Ana(outFile){}

void Embedding::Make()
{
   if( DEBUG )
      cout<<"Embedding::Make() called"<<endl;

   PARTICLE PART = PION;
   //PARTICLE PART = KAON;
   //PARTICLE PART = PROTON;

   int gePidPlus;
   int gePidMinus;
   if( PART == PION){
      gePidPlus = 211; // pi+   
   }else if( PART == KAON){
      gePidPlus = 321; // K+ 
   }else{
      gePidPlus = 2212; // p
   }
   gePidMinus = -gePidPlus;    

   hAnaFlow->Fill(kAll);
   if( mUpcEvt->getNumberOfMCParticles() < 2)
      return;
   hAnaFlow->Fill(kTwoMC);

   const TParticle* mc[nSigns];
   int generatedMCParticles = 0;
   for (int i = 0; i < mUpcEvt->getNumberOfMCParticles(); ++i)
   {
      const TParticle* mcPart = mUpcEvt->getMCParticle(i);
      if( abs(mcPart->GetPdgCode()) !=  gePidPlus)
         continue;
      if( mcPart->GetFirstMother() != 1)
         continue;

      mc[ (mcPart->GetPdgCode() > 0) ? PLUS : MINUS] = mcPart;
      generatedMCParticles++;
   }
   if( generatedMCParticles != 2)
   {
      cout<<"Error in Embedding::Make(): wrong number of MC particles = "<<generatedMCParticles<<endl;
      return;
   }

   
   if( mc[PLUS]->GetPdgCode() != gePidPlus || mc[MINUS]->GetPdgCode() != gePidMinus )
   {
      //cout<<"mc[PLUS] = "<<mc[PLUS]->GetPdgCode()<<" mc[MINUS] = "<<mc[MINUS]->GetPdgCode()<<endl;
      return;
   }
   if (abs(mc[PLUS]->Vz() - mc[MINUS]->Vz() ) > 3)
      return;
   hAnaFlow->Fill(kRightMC);

   for (int iMC = 0; iMC < nSigns; ++iMC)
   {
      mRecTree->setPtInGev(   mc[iMC]->Pt(), iMC, TRUEMC);  
      mRecTree->setEta( mc[iMC]->Eta(), iMC, TRUEMC); 
      double phi = mc[iMC]->Phi();
      if( phi > mUtil->pi() )
         phi -= 2*mUtil->pi();
      mRecTree->setPhi( phi, iMC, TRUEMC); 
      mRecTree->setMomentumInGev(   mc[iMC]->P(), iMC, TRUEMC);
      mRecTree->setCharge( iMC == 0 ? 1 : -1, iMC, TRUEMC);
   }

   mRecTree->setVertexXInCm(mc[PLUS]->Vx(), TRUEMC);
   mRecTree->setVertexYInCm(mc[PLUS]->Vy(), TRUEMC);
   mRecTree->setVertexZInCm(mc[PLUS]->Vz(), TRUEMC);
   
   //if( mc[PLUS]->Vz() > vertexRange )
   //   return;
   
   hNRecoVertices->Fill(mUpcEvt->getNumberOfVertices());

   vector<int> vertexIds;
   for (int iVrtx = 0; iVrtx < mUpcEvt->getNumberOfVertices(); ++iVrtx)
   {
      const StUPCVertex* vertex = mUpcEvt->getVertex(iVrtx);
      hVertexDist[X]->Fill(vertex->getPosX() - mc[PLUS]->Vx());
      hVertexDist[Y]->Fill(vertex->getPosY() - mc[PLUS]->Vy());
      hVertexDist[Z]->Fill(vertex->getPosZ() - mc[PLUS]->Vz());
      if( TMath::Abs( vertex->getPosZ() - mc[PLUS]->Vz() ) > 2 )
         continue;
      vertexIds.push_back(iVrtx);
   }

   hNVerticiesSelected->Fill( vertexIds.size() ); 
   if( vertexIds.size() != 0 ){
      hAnaFlow->Fill(kVertexReco);
      mRecTree->setVertexXInCm( mUpcEvt->getVertex( vertexIds[0] )->getPosX() );
      mRecTree->setVertexYInCm( mUpcEvt->getVertex( vertexIds[0] )->getPosY() );
      mRecTree->setVertexZInCm( mUpcEvt->getVertex( vertexIds[0] )->getPosZ() );
   }

   int totalCharge = 0;
   vector<int> hadronId;

   const double partMass = mUtil->mass(PART);//*mUtil->c();
   const double partMassSq = partMass*partMass; 

   for (int i = 0; i < nSigns; ++i)
      idTruth[i] = -1;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);

      if( !trk->getFlag(StUPCTrack::kCEP) ) continue;
      if( !(trk->getQaTruth() > 95) ) continue; // impemented in StUPCFilterMaker.cxx > 90
      int charge = trk->getCharge();
      if( TMath::Abs( charge ) != 1) continue; // impemented in StUPCFilterMaker.cxx
      if( trk->getIdTruth() != 1 && trk->getIdTruth() != 2) continue;
      if( !IsGoodGlobalTrack(trk) ) continue;
      
      int id = trk->getCharge() == 1 ? PLUS : MINUS;

      mRecTree->SaveTrackInfo(trk, id );
      double timeCaclulated = ( mRecTree->getTofLengthInCm(id)/mUtil->c()) * sqrt(1 + ((partMassSq)/ pow(mRecTree->getMomentumInGeV(id),2))) * pow(10.0,7.0);
      tofTimeOrig[id] = mRecTree->getTofTimeInNs(id);
      mRecTree->setTofTimeInNs(timeCaclulated, id);
      idTruth[id] = trk->getIdTruth();
      qaTruth[id] = trk->getQaTruth();
      tofMatched[id] = trk->getFlag(StUPCTrack::kTof);
      hadronId.push_back(trackID);
      totalCharge += static_cast<int>( trk->getCharge() );
   }
   hNCep->Fill(hadronId.size());

   treeState = 0;
   if(hadronId.size() == nSigns){
      hAnaFlow->Fill(kTwoRecoTrack);
      treeState = 1;
      if(totalCharge == 0){
         hAnaFlow->Fill(kOppositeTrack);
         treeState = 2;
         if( idTruth[PLUS] == 1 && idTruth[MINUS] == 2){
            hAnaFlow->Fill(kCorrectID);
            mRecTree->CalculatePID(false, true, PART);
            treeState = 3;
         }
      }
   }

   mRecTree->FillRecTree();

   return;

}

void Embedding::Init()
{
   if( DEBUG )
      cout<<"Embedding::Init() called"<<endl;
   mOutFile->cd();


   mRecTree = new RecTree(nameOfTree[kEMBEDING], treeBits[kEMBEDING], false);  
   mTree = mRecTree->getTTree();
   for (int i = 0; i < nSigns; ++i){
      mTree->Branch(Form("tofTimeOrig%i",i), &tofTimeOrig[i]);
      mTree->Branch(Form("idTruth%i",i), &idTruth[i]);
      mTree->Branch(Form("qaTruth%i",i), &qaTruth[i]);
      mTree->Branch(Form("tofMatched%i",i), &tofMatched[i]);
   }
   mTree->Branch("treeState", &treeState);
   mTree->Branch("vpdTimeDiff", &vpdTimeDiff);

   TString anaCutName[] = { TString("All"), TString("TwoMC"), TString("RightMC"), TString("VertexReco"), TString("2recoTrks"), TString("opossite charge"), TString("Correct ID") };
   hAnaFlow = new TH1D("EmbedAnalysisFlow", "CutsFlow", nAnaCuts-1, 1, nAnaCuts);
   for(int tb=1; tb<nAnaCuts; ++tb) 
      hAnaFlow->GetXaxis()->SetBinLabel(tb, anaCutName[tb-1]);

   hNRecoVertices = new TH1D("hNRecoVertices", ";#reco verticies per event;Counts", 20, 0, 19);
   hNVerticiesSelected = new TH1D("hNVerticiesSelected", ";#selected verticies per event;Counts", 20, 0, 19);
   for (int i = 0; i < nCoordinates; ++i)
      hVertexDist[i] = new TH1D("hVertexDist"+mUtil->coordinateName(i), ";(reco - MC)." + mUtil->coordinateName(i) + ";Counts", 200, -10, 10);

   hNReco = new TH1D("hNReco", ";#reco tracks per event;Counts", 30, 0, 29);
   hNCep = new TH1D("hNCep", ";#cep tracks per event;Counts", 30, 0, 29);
   hNCepTof = new TH1D("hNCepTof", ";#cep tof tracks per event;Counts", 30, 0, 29);

   mOutFile->cd();
   if( DEBUG )
      cout<<"Embedding::Init() finished"<<endl;
}
