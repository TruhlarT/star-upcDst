#include "Ana.h"
#include "RunDef.h"

Ana::Ana(TFile *outFile){
   mOutFile = outFile;
   mUtil = new Util();
}

Ana::~Ana(){
   if(mUtil){ delete mUtil; mUtil = nullptr; } 
}

void Ana::SetEvent(StUPCEvent *upcEvt, StRPEvent *rpEvt, StRPEvent *mcEvt){
   mUpcEvt = upcEvt;
   mRpEvt = rpEvt;
   mMcEvt = mcEvt;
}

//void Ana::SetRpPosition(TVector3 (&corr)[nRomanPots], TVector3 (&offsets)[nRomanPots])
void Ana::SetRpPosition(TVector3* corr, TVector3 *offsets)
{
   mCorrection = corr;
   mOffSet = offsets;
}

// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
//const int Util::triggerID[nTriggers] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
//           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};


bool Ana::RPInFidRange(double x, double y, int br) const{
   if( (br == ED || br == WD) && y > 0 )
      return false;
   if( (br == EU || br == WU) && y < 0 )
      return false;
   return ((fpLPRadius[br] == 0.0 ? abs(y) < fpPyMax[br] : (x + fpLPxCenter[br])*(x + fpLPxCenter[br]) + y*y < fpLPRadius[br]) && abs(y) > fpPyMin[br] && x > fpPxMin[br] && (x + fpPxCenter[br])*(x + fpPxCenter[br]) + y*y < fpPRadius[br]) ? true : false;
}


bool Ana::CheckTriggers(const vector<int> *triggerArray, StUPCEvent *mUpcEvt, TH1D *hTriggerBits) const{
   if( triggerArray->empty() )
      return true;
   bool triggered = false;
   for(int var = 0; var < nTriggers; ++var)
   {
      if(mUpcEvt->isTrigger(triggerID[var]))
      {  //Checked if it is CPT trigger
         if(hTriggerBits)
            hTriggerBits->Fill(var);
         for (unsigned int i = 0; i < triggerArray->size(); ++i)
            if(triggerID[var] == (*triggerArray)[i])
               triggered=true;
      }
   }
   return triggered;
}


void Ana::AnaRpTracks(StRPEvent *event)
{
   for (unsigned int i = 0; i < nSides; ++i)
      mRpTrackIdVec_perSide[i].clear();
   for (unsigned int i = 0; i < nBranches; ++i)
      mRpTrackIdVec_perBranch[i].clear();
   for (unsigned int i = 0; i < nRomanPots; ++i)
      mTrackPointIdVec[i].clear();
   for (unsigned int i = 0; i < nRomanPots*nPlanes; ++i)
      mClusterIdVec[i].clear();

   // Loop over all tracks reconstructed in Roman Pots    
   for(unsigned int k = 0; k < event->getNumberOfTracks(); ++k)
   {
      // Get pointer to k-th track in Roman Pot data collection
      StUPCRpsTrack *trk = event->getTrack(k);
      trk->setEvent(event);
      // Get ID of a side in which this k-th track was reconstructed
      int side = trk->branch()<2 ? E : W;

      if( (trk->getTrackPoint(0) ? trk->getTrackPoint(0)->planesUsed()>=nPlanesUsed : false) &&
         (trk->getTrackPoint(1) ? trk->getTrackPoint(1)->planesUsed()>=nPlanesUsed : false))
      {
         mRpTrackIdVec_perSide[side].push_back( k );
         mRpTrackIdVec_perBranch[trk->branch()].push_back( k );
      }
   }

   // Loop over all tracks points reconstructed in Roman Pots
   for(unsigned int k = 0; k < event->getNumberOfTrackPoints(); ++k)
   {
      StUPCRpsTrackPoint *trkPoint = event->getTrackPoint(k);
      if( trkPoint->planesUsed() < nPlanesUsed) // Skip "bad" track points
         continue;
      mTrackPointIdVec[trkPoint->rpId()].push_back( k );
   }

   // Loop over all clusters reconstructed in Roman Pots
   for(unsigned int k = 0; k < event->getNumberOfClusters(); ++k)
   {
      StUPCRpsCluster *cluster = event->getCluster(k);
      mClusterIdVec[cluster->romanPotId()*nPlanes + cluster->planeId()].push_back( k );
   }
}


/*
void Ana::CreateTofClusters()
{
   vector<TofCluster> tofClusterVec;
   int nTofClusters = 0;

   //CLUSTERING BASED ON ETA & PHI OF TOF ELEMENT
   for(unsigned int k=0; k< mUpcEvt->getNumberOfHits(); ++k)
   {
      const StUPCTofHit* hit = &mUpcEvt->getHit(k);
      int tray_A = hit->getTray();
      int module_A = hit->getModule();
      int cell_A = hit->getCell();  

      // skipping vpd hits
      if( (module_A==0 && tray_A==121) || (module_A==0 && tray_A==122) ) 
         continue;

      double phi_A = mTofPhi[tray_A-1][module_A-1][cell_A-1];
      double eta_A = mTofEta[tray_A-1][module_A-1][cell_A-1];
      double x_A = mTofX[tray_A-1][module_A-1][cell_A-1];
      double y_A = mTofY[tray_A-1][module_A-1][cell_A-1];
      double z_A = mTofZ[tray_A-1][module_A-1][cell_A-1];

      TVector2 vec_A( cos(phi_A), sin(phi_A) );

      bool attachToExistingCluster = false;
      int existingClusterIndex = -1;

      for(unsigned int i=0; i<tofClustersVec.size(); ++i)
      {
         for(unsigned int j=0; j<tofClustersVec[i].mTofHitIndexVec.size(); ++j)
         {
            const StUPCTofHit* hit_B = &mUpcEvt->getHit( tofClustersVec[i].mTofHitIndexVec[j] );
            int tray_B = hit_B->getTray();
            int module_B = hit_B->getModule();
            int cell_B = hit_B->getCell();  

            // skipping vpd hits
            if( (module_B==0 && tray_B==121) || (module_B==0 && tray_B==122) ) 
               continue;

            double phi_B = mTofPhi[tray_B-1][module_B-1][cell_B-1];
            double eta_B = mTofEta[tray_B-1][module_B-1][cell_B-1];

            TVector2 vec_B( cos(phi_B), sin(phi_B) );

            double distanceInPhi = vec_A.DeltaPhi( vec_B );
            double distanceInEta = fabs( eta_A - eta_B );

            if( sqrt(distanceInPhi*distanceInPhi + distanceInEta*distanceInEta) < 0.1)
            {
               attachToExistingCluster = true;
               existingClusterIndex = i;
            }
         }
      }

      if( !attachToExistingCluster )
      {
         tofClustersVec.push_back( TofCluster() );
         existingClusterIndex = nTofClusters;
         ++nTofClusters;
      }

      tofClustersVec[existingClusterIndex].mClusterSize++;
      tofClustersVec[existingClusterIndex].mTofHitIndexVec.push_back( k );
      tofClustersVec[existingClusterIndex].mHitXVec.push_back( x_A );
      tofClustersVec[existingClusterIndex].mHitYVec.push_back( y_A );
      tofClustersVec[existingClusterIndex].mHitZVec.push_back( z_A );
   }
}
*/
