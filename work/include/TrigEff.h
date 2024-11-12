#ifndef TrigEff_h
#define TrigEff_h

#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class TrigEff : public Ana{
   public:
      TrigEff(TFile *outfile);
      ~TrigEff(); 

      void Init() override;
      void Make() override;

   private:
      void CheckRPTrigEff();
      void FillRPInfo();

      enum ANA_CUT { ALL = 1, TRIG, ONETOFVX, TWOTOFTRKS, nCuts };

      static const TString mCutName[nCuts];

      // RP ADC and TAC
      TH1D *hRpAdc[2*nRomanPots];
      TH1D *hRpAdcInWindow[2*nRomanPots];
      TH1D *hRpTac[2*nRomanPots];
      TH2D *hNTpPerRP, *hNClustersPerPlane;

      // Efficiency plots
      TH1D *hRpTrigEffPtMiss[2][2][nRpConfigurations]; // [total,passed] [no pxpy cut, pxpy cut]
      TH2D *hRpTrigEffPxPy[2][2][nBranches][nRpConfigurations];
      TH2D *hRpTrigEffXY[2][2][nRomanPots][nRpConfigurations];
      TH2D *hRpTrigEffXYOffSub[2][2][nRomanPots][nRpConfigurations];
};

#endif
