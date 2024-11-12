#ifndef VertexStudy_h
#define VertexStudy_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class VertexStudy : public Ana{
   public:
      VertexStudy(TFile *outFile);
      ~VertexStudy(){if(mVertexTree) delete mVertexTree;};
      
      void Make() override;
      void Init() override;

   private:
      RecTree* mVertexTree; 
      void studyVertexRecoEff();
      pair<bool,bool> findPrimaryTracks(const StUPCTrack* &trkPrimPlus, const StUPCTrack* &trkPrimMinus, double vertexZhypo);
      void studyVertexZDistibution();

      TH2D *hDca;
      TH1D *hAnaFlow; 
      enum { kAll = 1, kTwoTracks, kSameSign, kDcaPart, kDcaBeam, kPrimAll, kPrimTwoTracks, kPrimSameSign, kPrimVertex, kMax};

};

#endif
