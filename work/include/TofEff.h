#ifndef TofEff_h
#define TofEff_h

#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class TofEff : public Ana{
   public:
      TofEff(TFile *outfile);
      ~TofEff(); 

      void Init() override;
      void Make() override;

   private:
      TTree *tofEffTree;
      
      TVector3 mRPpTBalance;

      double pTMiss, invMass, deltaPhiProtons, pTState, phiState, thetaState;
      double probePt, probeEta, probePhi;
      double tagPt, tagEta, tagPhi;
      bool probeTofHit, probeTof;

      void CalculateTOFEff(unsigned int tagID);
};

#endif
