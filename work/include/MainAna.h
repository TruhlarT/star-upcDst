#ifndef MainAna_h
#define MainAna_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class MainAna : public Ana{
   public:
      MainAna(TFile *outFile);
      ~MainAna(){};
      
      void Make() override;
      void Init() override;

   private:
      RecTree* mHitsTree[nTPCnHitsStudies];
      RecTree* mPidTree[nPidVariation];
      RecTree* mDcaTree[nDcaVariation];
      RecTree* mCurrentTree;
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

      TH1D *hMSquared[4][3][4];

      // RP ADC and TAC
      TH1D *hRpAdc[2*nRomanPots];
      TH1D *hRpAdcInWindow[2*nRomanPots];
      TH1D *hRpTac[2*nRomanPots];

      // Main Ana Plots
      TH2D *hRPFV[nBranches];
      TH1D *hNRPTracksInFV[nBranches];
      TH1D *hNTOFVerticies, *hZVertex, *hNTOFGoodTracks;
      TH1D *hNHitFit, *hNHitDedx, *hDcaXY, *hDCAZ;
      TH2D *hEtaZVertex;


      void runMainAna(RecTree* recTree, VARIATION nHitFitStatus = NOMINAL, VARIATION nHitDEdxStatus = NOMINAL, VARIATION pid = NOMINAL, unsigned int dca = NOMINAL);
      Int_t PIDStartegy(int strategy);
      bool IsPairOf(int type);

      void SaveMissingMomenta(TVector3 missP);
      void FillMSquared();
      
};

#endif 