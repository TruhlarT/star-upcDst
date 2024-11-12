#ifndef TofQA_h
#define TofQA_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class TofQA : public Ana{
   public:
      TofQA(TFile *outFile);
      ~TofQA(){};
      
      void Make() override;
      void Init() override;

   private:
      TH2D *hEtaPhi;
};

#endif
