#include "PlotManager.h"

double PlotManager::getTpcEff(int sysStudy)
{
    if(sysStudy<=0)
        return tpcEff[-sysStudy][ 2*mCurrentTree->getPairID() ]*tpcEff[-sysStudy][ 2*mCurrentTree->getPairID() + 1 ];

    double eff = 1;

    for (int i = 0; i < nSigns; ++i)
    {
        double partEff = 1;
        if(sysStudy == 1)
            partEff = hTPCPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPhi(i) ) ); 
        if(sysStudy == 2)
            partEff = hTPCEtaZEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCEtaZEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getVertexZInCm() ) ); 
        if(sysStudy == 3){
            if( mCurrentTree->getPtInGev(i) < 2.5)
                partEff = hTPCPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPtEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPtInGev(i) ) ); 
            else // some pions pairs very high pT so there is no correction for them and they are overflowing the final histograms 
                partEff = hTPCPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPtEff[mCurrentTree->getPairID()][i]->FindBin( 2.45 ) );
        }
        if(sysStudy == 4)
            partEff = hTPCEtaPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCEtaPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getPhi(i) ) ); 

        if( partEff == 0)
            cout<<"Warning in PlotManager::getTpcEff the efficiency for sysStudy: "<<sysStudy<<" is 0: "<<mCurrentTree->getEta(i)<<" ; "<< mCurrentTree->getVertexZInCm()<<endl;
        eff*=partEff;
    }

    return eff;
}//getTpcEff



double PlotManager::getPIDEff() 
{ 
    //cout<<"MAX: "<< max(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS))<<" MIN: "<<min(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS))<<endl;
    int id = mCurrentTree->getPairID();
    if(hPIDEff[id] == nullptr){
            cerr<<"Error: PlotManager::getPIDEff() the hist for "<<mUtil->particleName(id) <<" is not set..."<<endl;
            return 1.0;
    }
    double x = max(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS));
    double y = min(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS));
    double eff = hPIDEff[id]->GetBinContent( hPIDEff[id]->FindBin(x,y));
    return eff > 0.00001 ? eff : 1.0; // some pions pairs very high pT so there is no correction for them and they are overflowing the final histograms 
}//getPIDEff

/*
double PlotManager::getLuminosity() // done
{
    //return mIntegLumiPerRun[ mCurrentTree->getRunNumber() ];
     return correctedIntegLum; 
    // return (518.075/670.845)*142.006 * 1000; // convert from pb-1 to nb-1
    // From the table: https://www.star.bnl.gov/protected/common/common2017/trigger2017/lumipp500GeV/lum_pertriggerid_pp2017_500GeV.html
    // Get lumi for RP_CPT2noBBCL = 142.006 Lum [pb-1]
    // Get number of events for the trigger = 670.845 M
    // Get number of events analyzed for the trigger (from AnaFlow) = 518.075
}//getLuminosity
*/