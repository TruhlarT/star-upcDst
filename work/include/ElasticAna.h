#ifndef ElasticAna_h
#define ElasticAna_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class ElasticAna : public Ana{
   public:
      ElasticAna(TFile *outfile);
      ~ElasticAna(); 

      void Init() override;
      void Make() override;

      TVector3 CalculateOffsetCorrForRP(unsigned int Rp, unsigned int iter);

      void InitRPMCInfo(); 
      void SetRPMCInfo(double *mc_vrtx, double (*mc_p)[nSides]);

   private:
      // Control plots
      TH1D *hElasticAnaFlow; 
      enum ANA_CUT { ALL = 1, TRIG, TWORPTRKS, INFID, COLLINEAR, nCuts };
      static const TString mCutName[nCuts];
      enum { kAll = 1, kTtrig, kETPattern, kFourPoints, kTRange, kColinear, kFiducial, kDCut, kElastic, kMax};
      static const TString mElAnaCutsName[kMax];

      typedef struct track_t{
        float posX[4];   // point position at RP 
        float posY[4];   // point position at RP 
        float posZ[4];   // point position at RP
        unsigned int nPoints; 
        float thX;
        float thY;
        float phi;
        float t;
        float X0;
        float Y0;  // X,Y track position at z=0
      } TRACK_T;

      track_t trackEW, trackE, trackW;

      TH1D *hTheta[2][2];
      TH1D *hDeltaTheta[2];
      TH2D *hThetaCorr;
      TH1D *hDCutR;
      TH2D *hDCut;

      TH2D *hVertexExtraction[nArms][nCoordinates -1];
      TH1D *hVertexZExtraction[nArms], *hDeltaProjection[nCoordinates-1][nBranches], *hDeltaFitProjection[nRomanPots][nCoordinates-1];
      TH1D *hClusterLength[nRomanPots], *hClusterEnergy[nRomanPots], *hClusterPerPlain[nRomanPots][nPlanes], *hNClusterPerRP[nRomanPots];
      TH1D *hT;
      TH1D *hFindHitDiff;

      // pX, pY and x,y RP accaptance
      TH2D *hAcceptancePxPy[nBranches];
      TH2D *hAcceptanceXY[nRomanPots];

      // RP ADC and TAC
      TH1D *hRpAdc[2*nRomanPots];
      TH1D *hRpAdcInWindow[2*nRomanPots];
      TH1D *hRpTac[2*nRomanPots];

      // efficiency plots
      TH2D *hEffPxPy[2][nBranches];
      TH2D *hEffXY[2][nRomanPots];
      TH2D *hEffXYOffSub[2][nRomanPots];

      // alignment calculation
      vector<double> alignment[nRomanPots][nXYCoordinates];
      TH1D *hAligEvents[nRomanPots][nAligIteration], *hAligXCorr[nRomanPots][nAligIteration], *hAligYCorr[nRomanPots][nAligIteration]; 

      UInt_t mRunNumber;


      bool checkRP(unsigned int iRP);
      void FillEffPlots( unsigned int iRP, bool passed);
      void getFitPoint(double& x, double& y, track_t track, unsigned int rp);

      void FillElasticPlots();
      void FillClusterInfo();
      void FillAccaptancePlots();

      bool IsElasticEventCandidate(int RPStudy = -1);
      bool InFiducial();
      bool IsInDCut();
      bool IsColinear();

      void runMCValidation();
      void runVertexExtraction();
      void closerTest(const vector<TVector3> trackPoints, int arm);
      void runAlignment();
      void runRpDataDrivenEffStudy();

      track_t makeTrack(bitset<8> rpsBits, int copyRP);
      bool IsInGeoWindow( double phi){ return ( abs(abs( phi )- PHI_CENT ) < PHI_WIDTH); }
};

#endif
