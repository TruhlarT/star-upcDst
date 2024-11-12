#ifndef Util_h
#define Util_h

//_____________________________________________________________________________
//    Class to handle all useful utils
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "Libreries.h"

using namespace std;
#include "AnaConfig.h"

namespace UTIL {
   class Util; // just tell the compiler to expect a class def

   // Enumerations - very helpful and convenient !
   enum STUDYMAP { kMAINANA = 0, kVERTEXSTUDY, kEMBEDING, kTOFQA, kTRIGEFF, kFULLZB, kELASTICANA, kRPMCANA, kALIGNMENT, nStudies };
   enum SIDE { E=0, East=0, W=1, West=1, nSides };
   enum RPPORIENTATIONS { Up=0, Down=1, nRpOrientations };
   enum RPCONFIGURATION { IT=0, ET=1, nRpConfigurations};
   enum XY_COORDINATE { X = 0, Y, Z, nCoordinates, nXYCoordinates = Z};
   enum TRIGGER_ID { UPC_v1, UPC_v2, SDT, ET_v1, ET_v2, CPT2_v1, CPT2_v2, CPT2noBBCL, Zerobias, JPsiHTTP_v1, 
   JPsiHTTP_v2, JPsiHTTP_v3, SDT_RHICf, ET_RHICf, CPT2_RHICf, CPT2noBBCL_RHICf_v1, CPT2noBBCL_RHICf_v2, nTriggers };
   enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
   enum PLANE_ID {A, B, C, D, nPlanes };
   enum BRANCH_ID { EU, ED, WU, WD, nBranches };
   enum BRANCHES_CONFIGURATION_ID { CONF_EU_WU=0, CONF_ED_WD, CONF_EU_WD, CONF_ED_WU, nBranchesConfigurations };
   enum ARM_ID { EU_WD, ED_WU, nArms };
   enum STATION_ID { E1, E2, W1, W2, nStations };
   enum STATION_ORDER { RP1, RP2, nStationPerSide, nRpPerStation = nStationPerSide};
   enum PARTICLE { PION = 0, KAON, PROTON, nParticles };
   enum TPC_TRACK_TYPE { GLO, PRI, TOF, QUA, nTpcTrkTypes }; // GLO=global(all), PRI=primary, TOF=PRI&TofMatched, QUA=TOF&QualityCuts
   enum BUNCH_CROSSING { CB, AG, nBnchXngsTypes }; // CB=colliding bunches, AG=abort gaps
   enum QSUM_2TRKS { OPPO, SAME, nCharges2Trks };
   enum SIGN { PLUS = 0, MINUS, nSigns };
   enum LIST_OF_EFF_CORRECTIONS { RPACC, TPCRECOEFF, TOFMATCHEFF, nEffCorrections };
   enum ANALYSIS_CUT { ALL = 1, TRIG, TWOINFID, ONETOFVX, ZVERTEX, TWOTOFTRKS, ETA, OPPOSITE, EXCLUSIVE, PIPI, KK, PPBAR, nAnalysisCuts };
   enum RANGE_LIMIT { MIN, MAX };
   enum DATASET { MC = 0, MCZB, DATA, nDataSets };
   enum DATATAG { TRUEMC = 0, RECO, nDataTag };
   enum NHITSVARIATION { nHitsDEdxLoose = 0, nHitsFitLoose, nHitsDEdxTight, nHitsFitTight, nHitsVariation };
}  
using namespace UTIL;

class UTIL::Util{
   public:
      Util();
      ~Util();
      inline TString sideName(UInt_t id) const { if(id<nSides) return mSideName[id]; else{ std::cerr << "ERROR in Util::sideName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString coordinateName(UInt_t id) const { if(id<nCoordinates) return mCoordinateName[id]; else{ std::cerr << "ERROR in Util::coordinateName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString triggerName(UInt_t id) const { if(id<nTriggers) return mTriggerName[id]; else{ std::cerr << "ERROR in Util::triggerName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString armName(UInt_t id) const { if(id<nArms) return mArmName[id]; else{ std::cerr << "ERROR in Util::armName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString branchName(UInt_t id) const { if(id<nBranches) return mBranchName[id]; else{ std::cerr << "ERROR in Util::branchName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString branchesConfigurationName(UInt_t id) const { if(id<nBranchesConfigurations) return mBranchesConfigurationName[id]; else{ std::cerr << "ERROR in Util::branchesConfigurationName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString rpName(UInt_t id) const { if(id<nRomanPots) return mRpName[id]; else{ std::cerr << "ERROR in Util::rpName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString rpConfigName(UInt_t id) const { if(id<nRpConfigurations) return mRpConfigName[id]; else{ std::cerr << "ERROR in Util::rpConfigName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString rpConfigTag(UInt_t id) const { if(id<nRpConfigurations) return mRpConfigTag[id]; else{ std::cerr << "ERROR in Util::rpConfigName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString planeName(UInt_t id) const { if(id<nPlanes) return mPlaneName[id]; else{ std::cerr << "ERROR in Util::planeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString stationName(UInt_t id) const { if(id<nStations) return mStationName[id]; else{ std::cerr << "ERROR in Util::stationName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString particleName(UInt_t id) const { if(id<nParticles) return mParticleName[id]; else{ std::cerr << "ERROR in Util::particleName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString particleTag(UInt_t particle, UInt_t charge) const { if(particle<nParticles && charge<nSigns) return mParticleTag[charge][particle]; else{ std::cerr << "ERROR in Util::particleTag(UInt_t particle, UInt_t charge): parameters out of range" << std::endl; return TString("");} }
      inline TString particleLabels(UInt_t id) const { if(id<nParticles) return mParticleLables[id]; else{ std::cerr << "ERROR in Util::particleLabels(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString pairLabel(UInt_t id) const { if(id<nParticles) return mPairLabel[id]; else{ std::cerr << "ERROR in Util::pairLabel(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString tpcTrackTypeName(UInt_t id) const { if(id<nTpcTrkTypes) return mTpcTrackTypeName[id]; else{ std::cerr << "ERROR in Util::tpcTrackTypeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString bunchXngTypeName(UInt_t id) const { if(id<nBnchXngsTypes) return mBunchCrossingTypeName[id]; else{ std::cerr << "ERROR in Util::bunchXngTypeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString qSum2TrksName(UInt_t id) const { if(id<nCharges2Trks) return mChargeSum2TrksName[id]; else{ std::cerr << "ERROR in Util::qSum2TrksName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString signName(UInt_t id) const { if(id<nSigns) return mSignName[id]; else{ std::cerr << "ERROR in Util::signName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString efficiencyName(UInt_t id) const { if(id<nEffCorrections) return mEfficiencyName[id]; else{ std::cerr << "ERROR in Util::efficiencyName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisCutName(UInt_t id) const { if(id<=nAnalysisCuts) return mCutName[id-1]; else{ std::cerr << "ERROR in Util::analysisCutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString dataSetName(UInt_t id) const { if(id<nDataSets) return mDataSetName[id]; else{ std::cerr << "ERROR in Util::dataSetName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString dataTagName(UInt_t id) const { if(id<nDataTag) return mDataTagName[id]; else{ std::cerr << "ERROR in Util::dataTagName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString nHitVaryName(UInt_t id) const { if(id<nHitsVariation) return mNHitVaryName[id]; else{ std::cerr << "ERROR in Util::nHitVaryName(UInt_t id): id out of range" << std::endl; return TString("");} }

      inline Double_t mass(int name) const{ return mParticleMass[name]; }
      inline Double_t c() const{ return mSpeedOfLight; }
      inline Double_t p0() const{ return mBeamMomentum; }
      inline Double_t pi() const{ return mPi; }
      inline Double_t epsiolon() const{ return mEpsilon; }

      inline Double_t rpZPosition(int rpId) const { return mRpZPosition[rpId]; };
      inline int branchPerRp(int rpId) const { return mBranchPerRp[rpId]; };
      inline int oppositeBranch(int br) const { return mOppositeBranch[br]; };
      inline int sidePerRp(int rpId) const { return mSidePerRp[rpId]; };
      inline int stationOrderPerRp(int rpId) const { return mStationOrderPerRp[rpId]; };
      inline int rpPerBranchStationOrder(int br, int st) const { return mRpPerBranchStationOrder[br][st]; }
      inline int branchPerBranchConfiguration(int brConf, int side) const { return mBranchPerBranchConfiguration[brConf][side]; }
      inline int planeToCoor(int plane) const { return plane%2 ? X : Y;} 
      inline int mapBranchConfiguration(int estOri, int wstOri) const { return mBranchConfigMap[estOri][wstOri];}

      Double_t binomialCoeff(UInt_t, UInt_t) const;
      TVector3 fitLine(const vector<TVector3> trackPoints, double positionOfFit);
      pair<double, double> CalculateMeanAndError(const vector<double>& values);

   private:
      // Labels, names etc. (defined as TString to gain higher functionality than const char*, e.g. defined "+" operator)
      vector<TString> mSideName;
      vector<TString> mCoordinateName;
      vector<TString> mTriggerName;
      vector<TString> mArmName;
      vector<TString> mBranchName;
      vector<TString> mBranchesConfigurationName;
      vector<TString> mRpName;
      vector<TString> mRpConfigName;
      vector<TString> mRpConfigTag;
      vector<TString> mPlaneName;
      vector<TString> mStationName;
      vector<TString> mParticleName;
      vector<TString> mParticleTag[nSigns];
      vector<TString> mParticleLables;
      vector<TString> mPairLabel;
      vector<TString> mTpcTrackTypeName;
      vector<TString> mBunchCrossingTypeName;
      vector<TString> mChargeSum2TrksName;
      vector<TString> mSignName;
      vector<TString> mEfficiencyName;
      vector<TString> mCutName;
      vector<TString> mDataSetName;
      vector<TString> mDataTagName;
      vector<TString> mNHitVaryName;
          
      Double_t mParticleMass[nParticles]; // GeV/c^2
      const Double_t mSpeedOfLight; // m/s
      const Double_t mBeamMomentum; // GeV/c
      const Double_t mPi;
      const Double_t mEpsilon;
      
      Double_t mRpZPosition[nRomanPots];
      BRANCH_ID mBranchPerRp[nRomanPots];
      BRANCH_ID mOppositeBranch[nBranches];
      SIDE mSidePerRp[nRomanPots];
      STATION_ORDER mStationOrderPerRp[nRomanPots];
      RP_ID mRpPerBranchStationOrder[nBranches][nStationPerSide];
      BRANCH_ID mBranchPerBranchConfiguration[nBranchesConfigurations][nSides];   
      BRANCHES_CONFIGURATION_ID mBranchConfigMap[nRpOrientations][nRpOrientations]; 


};

   void SumDistance2(int &, double *, double & sum, double * par, int);
   void line(double z, double *p, double &x, double &y);
   double distance2(double x,double y,double z, double *p);

#endif
