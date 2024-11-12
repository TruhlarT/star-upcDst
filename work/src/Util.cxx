#include "Util.h"

Util::Util(): mSpeedOfLight(299792458), mBeamMomentum(254.867), mPi(3.14159265359), mEpsilon(1e-7)
{
   mSideName.push_back( TString("East"));
   mSideName.push_back( TString("West"));
   
   mCoordinateName.push_back( TString("x"));
   mCoordinateName.push_back( TString("y"));
   mCoordinateName.push_back( TString("z"));
   
   mTriggerName.push_back( TString("UPC") ); // 570702 RP_UPC
   mTriggerName.push_back( TString("UPC") ); // 570712 RP_UPC
   mTriggerName.push_back( TString("SDT") ); // 570703 RP_SDT
   mTriggerName.push_back( TString("ET") ); // 570709 RP_ET
   mTriggerName.push_back( TString("ET") ); // 570719 RP_ET 
   mTriggerName.push_back( TString("CPT2") ); // 570701 RP_CPT2
   mTriggerName.push_back( TString("CPT2") ); // 570711 RP_CPT2
   mTriggerName.push_back( TString("CPT2noBBCL") ); // 570705 RP_CPT2noBBCL 
   mTriggerName.push_back( TString("Zerobias") ); // 570704 RP_Zerobias
   mTriggerName.push_back( TString("JPsi*HTTP" )); // 570209 JPsi*HTTP 
   mTriggerName.push_back( TString("JPsi*HTTP" )); // 570219 JPsi*HTTP
   mTriggerName.push_back( TString("JPsi*HTTP" )); // 570229 JPsi*HTTP
   mTriggerName.push_back( TString("SDT") ); // 590703 RP_SDT
   mTriggerName.push_back( TString("ET") ); // 590709 RP_ET
   mTriggerName.push_back( TString("CPT2") ); // 590701 RP_CPT2
   mTriggerName.push_back( TString("CPT2noBBCL") ); // 590705 RP_CPT2noBBCL
   mTriggerName.push_back( TString("CPT2noBBCL") ); // 590708 RP_CPTnoBBCL

   mArmName.push_back( TString("EU-WD") );
   mArmName.push_back( TString("ED-WU") );

   mBranchName.push_back( TString("EU") );
   mBranchName.push_back( TString("ED") );
   mBranchName.push_back( TString("WU") );
   mBranchName.push_back( TString("WD") );
   
   mBranchesConfigurationName.push_back( TString("EU-WU") );
   mBranchesConfigurationName.push_back( TString("ED-WD") );
   mBranchesConfigurationName.push_back( TString("EU-WD") );
   mBranchesConfigurationName.push_back( TString("ED-WU") );
      
   mRpName.push_back( TString("E1U") );
   mRpName.push_back( TString("E1D") );
   mRpName.push_back( TString("E2U") );
   mRpName.push_back( TString("E2D") );
   mRpName.push_back( TString("W1U") );
   mRpName.push_back( TString("W1D") );
   mRpName.push_back( TString("W2U") );
   mRpName.push_back( TString("W2D") );

   mRpConfigName.push_back( TString("IT") );
   mRpConfigName.push_back( TString("ET") );

   mRpConfigTag.push_back( TString("#Delta#varphi < 90#circ") );
   mRpConfigTag.push_back( TString("#Delta#varphi > 90#circ") );

   mPlaneName.push_back( TString("A") );
   mPlaneName.push_back( TString("B") );
   mPlaneName.push_back( TString("C") );
   mPlaneName.push_back( TString("D") );

   mStationName.push_back( TString("E1") );
   mStationName.push_back( TString("E2") );
   mStationName.push_back( TString("W1") );
   mStationName.push_back( TString("W2") );

   mParticleName.push_back( TString("pion") );
   mParticleName.push_back( TString("kaon") );
   mParticleName.push_back( TString("proton") );


   mParticleTag[PLUS].push_back( TString("#pi^{+}") );
   mParticleTag[MINUS].push_back( TString("#pi^{-}") );
   mParticleTag[PLUS].push_back( TString("K^{+}") );
   mParticleTag[MINUS].push_back( TString("K^{-}") );
   mParticleTag[PLUS].push_back( TString("p") );
   mParticleTag[MINUS].push_back( TString("#bar{p}") );

   for (int i = 0; i < nParticles; ++i){
      mParticleLables.push_back( mParticleTag[PLUS][i] + "," + mParticleTag[MINUS][i] );
      mPairLabel.push_back( mParticleTag[PLUS][i] + mParticleTag[MINUS][i] );
   }

   mTpcTrackTypeName.push_back( TString("global") );
   mTpcTrackTypeName.push_back( TString("primary") );
   mTpcTrackTypeName.push_back( TString("TofMatched") );
   mTpcTrackTypeName.push_back( TString("goodQuality") );
   mBunchCrossingTypeName.push_back( TString("collidingBunches") );
   mBunchCrossingTypeName.push_back( TString("abortGaps") );
   mChargeSum2TrksName.push_back( TString("oppositeSign") );
   mChargeSum2TrksName.push_back( TString("sameSign") );
   mSignName.push_back( TString("Plus") );
   mSignName.push_back( TString("Minus") );
   mEfficiencyName.push_back( TString("RpAcc") );
   mEfficiencyName.push_back( TString("TpcRecoEff") );
   mEfficiencyName.push_back( TString("TofMatchEff") );

   mCutName.push_back( TString("All") );
   mCutName.push_back( TString("CPT") );
   mCutName.push_back( TString("2 RP trks in FV") );
   mCutName.push_back( TString("1 TOF vrtx") );
   mCutName.push_back( TString(Form("|z_{vrtx}| < %.f",vertexRange)) );
   mCutName.push_back( TString("2 good TOF trks") );
   mCutName.push_back( TString("#eta cut") );
   mCutName.push_back( TString("Opposite-sign") );
   mCutName.push_back( TString("Exclusive (p_{T}^{miss} cut)") );
   mCutName.push_back( TString("#pi^{+}#pi^{-}") );
   mCutName.push_back( TString("K^{+}K^{-}") );
   mCutName.push_back( TString("p#bar{p}") );

   mDataSetName.push_back( TString("TrueMC") );
   mDataSetName.push_back( TString("MCZB") );
   mDataSetName.push_back( TString("Data") );
   mDataTagName.push_back( TString("TrueMc") );
   mDataTagName.push_back( TString("Reco") ); 
   

   mNHitVaryName.push_back( TString("nHitsDEdxLoose") );
   mNHitVaryName.push_back( TString("nHitsFitLoose") );
   mNHitVaryName.push_back( TString("nHitsDEdxTight") );
   mNHitVaryName.push_back( TString("nHitsFitTight") );

   mParticleMass[PION] = 0.13957018;
   mParticleMass[KAON] = 0.493667;
   mParticleMass[PROTON] = 0.93827208;
   
   mBranchPerRp[E1U] = EU;
   mBranchPerRp[E2U] = EU;
   mBranchPerRp[E1D] = ED;
   mBranchPerRp[E2D] = ED;
   mBranchPerRp[W1U] = WU;
   mBranchPerRp[W2U] = WU;
   mBranchPerRp[W1D] = WD;
   mBranchPerRp[W2D] = WD;

   mRpZPosition[E1U] = -15.78;
   mRpZPosition[E2U] = -17.58;
   mRpZPosition[E1D] = -15.78;
   mRpZPosition[E2D] = -17.58;
   mRpZPosition[W1U] = 15.78;
   mRpZPosition[W2U] = 17.58;
   mRpZPosition[W1D] = 15.78;
   mRpZPosition[W2D] = 17.58;
   
   mOppositeBranch[EU] = WD;
   mOppositeBranch[ED] = WU;
   mOppositeBranch[WU] = ED;
   mOppositeBranch[WD] = EU;
   
   mSidePerRp[E1U] = E;
   mSidePerRp[E2U] = E;
   mSidePerRp[E1D] = E;
   mSidePerRp[E2D] = E;
   mSidePerRp[W1U] = W;
   mSidePerRp[W2U] = W;
   mSidePerRp[W1D] = W;
   mSidePerRp[W2D] = W;
   
   mStationOrderPerRp[E1U] = RP1;
   mStationOrderPerRp[E2U] = RP2;
   mStationOrderPerRp[E1D] = RP1;
   mStationOrderPerRp[E2D] = RP2;
   mStationOrderPerRp[W1U] = RP1;
   mStationOrderPerRp[W2U] = RP2;
   mStationOrderPerRp[W1D] = RP1;
   mStationOrderPerRp[W2D] = RP2;
   
   mRpPerBranchStationOrder[EU][RP1] = E1U;
   mRpPerBranchStationOrder[EU][RP2] = E2U;
   mRpPerBranchStationOrder[ED][RP1] = E1D;
   mRpPerBranchStationOrder[ED][RP2] = E2D;
   mRpPerBranchStationOrder[WU][RP1] = W1U;
   mRpPerBranchStationOrder[WU][RP2] = W2U;
   mRpPerBranchStationOrder[WD][RP1] = W1D;
   mRpPerBranchStationOrder[WD][RP2] = W2D;

   mBranchPerBranchConfiguration[CONF_EU_WU][E] = EU;
   mBranchPerBranchConfiguration[CONF_EU_WU][W] = WU;
   mBranchPerBranchConfiguration[CONF_ED_WD][E] = ED;
   mBranchPerBranchConfiguration[CONF_ED_WD][W] = WD;
   mBranchPerBranchConfiguration[CONF_EU_WD][E] = EU;
   mBranchPerBranchConfiguration[CONF_EU_WD][W] = WD;
   mBranchPerBranchConfiguration[CONF_ED_WU][E] = ED;
   mBranchPerBranchConfiguration[CONF_ED_WU][W] = WU;  

   mBranchConfigMap[Up][Up] =  CONF_EU_WU;
   mBranchConfigMap[Down][Down] =  CONF_ED_WD;
   mBranchConfigMap[Up][Down] =  CONF_EU_WD;
   mBranchConfigMap[Down][Up] =  CONF_ED_WU;
}

Util::~Util()
{
}

// fit line through n-track points and return position of the fit in z = positionOfFit
TVector3 Util::fitLine(const vector<TVector3> trackPoints, double positionOfFit)
{
   TGraph2D * gr = new TGraph2D();
   //gr->SetPoint(0,x,y,z)
   for (unsigned int iTp = 0; iTp < trackPoints.size(); ++iTp)
      gr->SetPoint(iTp,trackPoints[iTp].X(),trackPoints[iTp].Y(),trackPoints[iTp].Z());
   
   TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
   Double_t arglist[10];
   arglist[0] = -1; // to silence output
   arglist[1] = 0.000001; // tolerance 
   min->ExecuteCommand("SET PRINT",arglist,1); // Silence output
   arglist[0] = 1000; // number of function calls 
   min->SetObjectFit(gr);
   min->SetFCN(SumDistance2);

   double pStart[4] = {1,1,1,1};
   //(Int_t ipar, const char *parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)
   min->SetParameter(0,"x0",pStart[0],0.01,0,0);
   min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
   min->SetParameter(2,"y0",pStart[2],0.01,0,0);
   min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

   min->ExecuteCommand("MIGRAD",arglist,2);

   // get fit parameters
   double parFit[4];
   for (int i = 0; i < 4; ++i) 
      parFit[i] = min->GetParameter(i);

   TVector3 fitResult;
   fitResult.SetZ(positionOfFit);
   line(fitResult[2],parFit,fitResult[0], fitResult[1]);
   return fitResult;
}


Double_t Util::binomialCoeff(UInt_t n, UInt_t k) const
{
   if(k > n)
      { cerr << "ERROR in Util::binomialCoeff(UInt_t n, UInt_t k):  k>n !!!!" << endl; return 1; }
   double result = 1.0;
   for(unsigned int i=1; i<=k; ++i) 
      result *= static_cast<double>(n-i+1)/static_cast<double>(i);
   return result;
}


// function to be minimized 
void SumDistance2(int &, double *, double & sum, double * par, int) 
{ 
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) { 
      double d = distance2(x[i],y[i],z[i],par); 
      sum += d;
   }
}

// define the parameteric line equation 
void line(double z, double *p, double &x, double &y) 
{ 
   x = p[0] + p[1]*z; 
   y = p[2] + p[3]*z;
} 

// calculate distance line-point 
double distance2(double x,double y,double z, double *p) 
{ 
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   ROOT::Math::XYZVector xp(x,y,z); 
   ROOT::Math::XYZVector x0(p[0], p[2], 0. ); 
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
   ROOT::Math::XYZVector u = (x1-x0).Unit(); 
   double d2 = ((xp-x0).Cross(u)) .Mag2(); 
   return d2; 
}


// Function to calculate mean and error on the mean
pair<double, double> Util::CalculateMeanAndError(const vector<double>& values) 
{
   // Check if the vector is empty
   if (values.empty()) {
      cerr << "Error: Input vector is empty." << endl;
      return make_pair(0.0, 0.0);
   }

   // Calculate the mean
   double sum = 0.0;
   for (double value : values) {
      sum += value;
   }
   double mean = sum / values.size();

   // Calculate the standard deviation
   double sumSquaredDiffs = 0.0;
   for (double value : values) {
      sumSquaredDiffs += (value - mean) * (value - mean);
   }
   double variance = sumSquaredDiffs / (values.size() - 1);  // Sample variance
   double stddev = sqrt(variance);

   // Error on the mean (standard deviation divided by sqrt of the number of entries)
   double errorOnMean = stddev / sqrt(values.size());

   // Return a pair containing the mean and the error on the mean
   return make_pair(mean, errorOnMean);
}//CalculateMeanAndError