#ifndef AnalysisManager_h
#define AnalysisManager_h

// include headers
#include "Libreries.h"
#include "RunDef.h"
#include "Util.h"
#include "MainAna.h"
#include "VertexStudy.h"
#include "TofQA.h"
#include "TofEff.h"
#include "TrigEff.h"
#include "FullZB.h"
#include "ElasticAna.h"
#include "RpMCAna.h"
#include "Embedding.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

TFile *inFile, *outFile, *mcFile;
TChain *upcChain, *mcChain;
StUPCEvent *upcEvt;
StRPEvent *rpEvt, *origRpEvt, *mcEvt, *correctedRpEvent;
TTree *upcTree, *mcTree;

Util* mUtil;
RpMCAna* mRpMCAna;
vector<Ana*> mAnaVector;

bool runRP;
double mc_vtx[nCoordinates]; // x,y,z position of vertex from MC
double mc_p[nCoordinates][nSides]; // x,y,z momenta of proton on East/West side from MC

Long64_t nEvents, nMcEvents, iMcEvnt, iZbEv;
vector<Long64_t> ZbEvents;

// event info
UInt_t mRunNumber;

map<unsigned int, TVector3> mCorrection[nRomanPots];
map<unsigned int, TVector3> mOffSet[nRomanPots];

void Make();
void Init();


bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]);
void runAfterburner(StRPEvent *event, StRPEvent *newRpEvent);
TVector3 CalculateMomentumVector(double TPx, double TPy, double TPz, StUPCRpsTrack *track);
void SetRpEvent();
void ReleaseTheMemoryOfCorrRpEvt();
void RunRpMcAna(Long64_t iev);
void CleanMemory();
void RunAlignment();
void SaveAlignment();
void setMcEventsPerZbEvent();

bool ConnectInput(int argc, char** argv);
TFile *CreateOutputFile(const string& out);


#endif
