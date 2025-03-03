#include "PlotManager.h"


double PlotManager::getTpcEff(UInt_t varApp, UInt_t nHitsFit, UInt_t nHitsDeDx)
{
    if( varApp == NOMINAL)    
        return tpcEff[mUtil->varyNHits(nHitsFit,nHitsDeDx)][ 2*mCurrentTree->getPairID() ]*tpcEff[mUtil->varyNHits(nHitsFit,nHitsDeDx)][ 2*mCurrentTree->getPairID() + 1 ];

    double eff = 1;

    for (unsigned int i = 0; i < nSigns; ++i)
    {
        double partEff = 1;
        if(varApp == 1)
            partEff = hTPCPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPhi(i) ) ); 
        if(varApp == 2)
            partEff = hTPCEtaZEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCEtaZEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getVertexZInCm() ) ); 
        if(varApp == 3)
            partEff = hTPCEtaPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCEtaPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getPhi(i) ) ); 
        /*if(varApp == 3){
            if( mCurrentTree->getPtInGev(i) < 2.5)
                partEff = hTPCPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPtEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPtInGev(i) ) ); 
            else // some pions pairs very high pT so there is no correction for them and they are overflowing the final histograms 
                partEff = hTPCPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTPCPtEff[mCurrentTree->getPairID()][i]->FindBin( 2.45 ) );
        }*/
        if( partEff == 0.0)
            cout<<"Warning in PlotManager::getTpcEff the efficiency for varApp: "<<varApp<<" is 0: "<<mCurrentTree->getEta(i)<<" ; "<< mCurrentTree->getVertexZInCm()<<endl;
        eff*=partEff;
    }

    return eff;
}//getTpcEff


double getTofEff(UInt_t var, UInt_t nHitsFit, UInt_t nHitsDeDx); 


double PlotManager::getTofEff(UInt_t varApp) 
{
    if(varApp == NOMINAL)     
        return tofEff[ 2*mCurrentTree->getPairID() ]*tofEff[ 2*mCurrentTree->getPairID() + 1 ]; 

    double eff = 1;

    for (unsigned int i = 0; i < nSigns; ++i)
    {
        double partEff = 1;

        if(varApp == 1)
            partEff = hTOFPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTOFPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPhi(i) ) ); 
        if(varApp == 2)
            partEff = hTOFEtaZEff[mCurrentTree->getPairID()][i]->GetBinContent( hTOFEtaZEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getVertexZInCm() ) ); 
        if(varApp == 3)
            partEff = hTOFEtaPhiEff[mCurrentTree->getPairID()][i]->GetBinContent( hTOFEtaPhiEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getEta(i), mCurrentTree->getPhi(i) ) ); 

/*        if(varApp == 3){
            if( mCurrentTree->getPtInGev(i) < 2.5)
                partEff = hTOFPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTOFPtEff[mCurrentTree->getPairID()][i]->FindBin( mCurrentTree->getPtInGev(i) ) ); 
            else // some pions pairs very high pT so there is no correction for them and they are overflowing the final histograms 
                partEff = hTOFPtEff[mCurrentTree->getPairID()][i]->GetBinContent( hTOFPtEff[mCurrentTree->getPairID()][i]->FindBin( 2.45 ) );
        }*/

        if( partEff == 0)
            cout<<"Warning in PlotManager::getTofEff the efficiency for varApp: "<<varApp<<" is 0: "<<mCurrentTree->getEta(i)<<" ; "<< mCurrentTree->getVertexZInCm()<<endl;
        eff*=partEff;
    }

    return eff;
};

double PlotManager::getPIDEff(UInt_t var) 
{ 
    //cout<<"MAX: "<< max(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS))<<" MIN: "<<min(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS))<<endl;
    int id = mCurrentTree->getPairID();
    if(hPIDEff[id][var] == nullptr){
            cerr<<"Error: PlotManager::getPIDEff() the hist for "<<mUtil->particleName(id) <<" is not set..."<<endl;
            return 1.0;
    }
    double x = max(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS));
    double y = min(mCurrentTree->getPtInGev(PLUS),mCurrentTree->getPtInGev(MINUS));
    double eff = hPIDEff[id][var]->GetBinContent( hPIDEff[id][var]->FindBin(x,y));
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