#ifndef Embedding_h
#define Embedding_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class Embedding : public Ana{
   public:
      Embedding(TFile *outFile);
      ~Embedding(){};
      
      void Make() override;
      void Init() override;

   private:
      TTree* mTree;
      double tofTimeOrig[nSigns], qaTruth[nSigns];
      int idTruth[nSigns];
      int treeState;
      double vpdTimeDiff;
      bool tofMatched[nSigns];

      enum ANA_CUTS { kAll = 1, kTwoMC, kRightMC, kVertexReco, kTwoRecoTrack, kOppositeTrack, kCorrectID, nAnaCuts };
      TH1D *hAnaFlow, *hNRecoVertices, *hNVerticiesSelected, *hNReco, *hNCep, *hNCepTof;
      TH1D *hVertexDist[nCoordinates];
};

#endif
