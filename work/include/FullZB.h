#ifndef FullZB_h
#define FullZB_h

#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class FullZB : public Ana{
   public:
      FullZB(TFile *outfile);
      ~FullZB(); 

      void Init() override;
      void Make() override;
};

#endif
