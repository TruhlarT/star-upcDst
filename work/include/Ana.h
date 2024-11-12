#ifndef Ana_h
#define Ana_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"

using namespace std;
using namespace UTIL;

class Ana{
   public:
      Ana(TFile *outFile);
      virtual ~Ana();

      inline bool RPInFidRange(double x, double y) const { return (abs(y) < fpPyMax[EU] && abs(y) > fpPyMin[EU] && x > fpPxMin[EU] && (x + fpPxCenter[EU])*(x + fpPxCenter[EU]) + y*y < fpPRadius[EU]) ? true : false;}
      bool RPInFidRange(double x, double y, int br) const;
      inline bool RPInElFidRange(double x, double y) const { return (abs(y) < fpElPyMax && abs(y) > fpElPyMin && x > fpElPxMin && x < fpElPxMax) ? true : false;}
      inline bool IsInRpRange(double x, double y, int rpId, TVector3 offSet) const { return (abs(y) < abs(offSet[Y]) || abs(y) > abs(offSet[Y]) + 0.03 || x < -0.02 || x > 0.02) ? false : true; };
      inline bool IsGoodGlobalTrack(const StUPCTrack *trk, unsigned int fit = 0, unsigned int dEdx = 0) const { return (trk->getNhitsFit() >= minNHitsFit[fit] && trk->getNhitsDEdx() >= minNHitsDEdx[dEdx]);}
      inline bool IsGoodEtaTrack(const StUPCTrack *trk) const { return (abs(trk->getEta()) < maxEta && trk->getEta() > vertexEtaRange(0) && trk->getEta() < vertexEtaRange(1) );}
      inline bool IsGoodTrack(const StUPCTrack *trk, unsigned int fit = 0, unsigned int dEdx = 0) const { return ( IsGoodGlobalTrack(trk,fit,dEdx) && trk->getDcaXY() < maxDcaXY && trk->getDcaZ() < maxDcaZ && trk->getDcaZ() > minDcaZ && trk->getPt() > minPt[PION]);}
      inline bool IsGoodTofTrack(const StUPCTrack *trk) const {return (trk->getTofTime() > 0 && trk->getTofPathLength() > 0);}
      bool CheckTriggers(const vector<int> *triggerArray, StUPCEvent *mUpcEvt, TH1D *hTriggerBits) const;

      inline double vertexEtaRange(int rSide ) const{ return mRecTree->IsVertexSet() ? etaVertexSlope*mRecTree->getVertexZInCm() + (2*rSide - 1)*etaVertexShift : (2*rSide - 1)*maxEta; } // 0 for min, 1 for max

      void AnaRpTracks(StRPEvent *event);

      virtual void Make(){cout<<"Hi my name is make"<<endl;};
      virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
      void SetEvent(StUPCEvent *upcEvt, StRPEvent *rpEvt, StRPEvent *mcEvt);
      inline void SetRunNumber(UInt_t runN){ mRunNumber = runN; };
      //virtual void SetRpPosition(TVector3 (&corr)[nRomanPots], TVector3 (&offsets)[nRomanPots]);
      virtual void SetRpPosition(TVector3* corr, TVector3 *offsets);

      inline void SetAnaName(TString name){ anaName = name; }
      inline void SetTriggers(const vector<int> *trigg){ trigger = trigg; }
   
   protected:
      TFile *mOutFile;
      RecTree* mRecTree; 

      StUPCEvent *mUpcEvt;
      StRPEvent *mRpEvt, *mMcEvt;

      Util* mUtil;

      UInt_t mRunNumber;

      vector<unsigned int> mRpTrackIdVec_perSide[nSides];
      vector<unsigned int> mRpTrackIdVec_perBranch[nBranches];
      vector<unsigned int> mTrackPointIdVec[nRomanPots];
      vector<unsigned int> mClusterIdVec[nRomanPots*nPlanes];

      // Control plots
      TH1D *hAnalysisFlow; 

      TString anaName;
      const vector<int> *trigger;

      // RP offsets
      TVector3 *mCorrection;
      TVector3 *mOffSet;


};

#endif
