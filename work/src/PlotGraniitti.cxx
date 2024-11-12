#include "PlotManager.h"

void PlotManager::initGraniittiPlots()
{
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        for (int iRpCon = 0; iRpCon < nRpConfigurations ; ++iRpCon)
        {
            TString dir = mUtil->particleName(iPart) + "_" + mUtil->rpConfigName(iRpCon);

            hInvMassGran[iPart][iRpCon] = (TH1F*)graniittiFile->Get(dir + "/h_invMass");
            hDeltaPhiGran[iPart][iRpCon] = (TH1F*)graniittiFile->Get(dir + "/h_deltaPhi");

            if( !hInvMassGran[iPart][iRpCon] || !hDeltaPhiGran[iPart][iRpCon] )
               std::cerr << "Error in PlotManager::initGraniittiPlots() while retrieving graniitti plots!" << std::endl;

            graniittiScaleFactor[iPart][iRpCon] = -999.9;
        }
    }
}//initGraniittiPlots


TH1F* PlotManager::GetGraniittiPlot(int part, int rpCon, TString hName, double dataIntegral, bool scale)
{
    TH1F *hist = nullptr;
    TH1F **graniittiPlots;

    if( hName == "DeltaPhi"){
        graniittiPlots = hDeltaPhiGran[part];
        if( rpCon < nRpConfigurations && !scale){
            double tmp = dataIntegral / graniittiPlots[rpCon]->Integral("width");
            tmp = abs(tmp - 1.0) < 0.15 ? 1.0 : tmp;
            double scndFactor = rpCon == IT ? graniittiScaleFactor[part][ET] : graniittiScaleFactor[part][IT];
            graniittiScaleFactor[part][rpCon] = abs(scndFactor - tmp) < 0.1 ? scndFactor : tmp;
            //cout<<"Setting SF for "<<hName<<" part "<< part <<" rpCon "<<rpCon<<" "<<graniittiScaleFactor[part][rpCon]<<" / "<<getGraniittiSF(part, rpCon)<<endl;
        }
    }else if( hName == "InvMass"){
        graniittiPlots = hInvMassGran[part];
    }
    else{
        return hist;
    }

    if( getGraniittiSF(part, rpCon) < 0){
        std::cerr << "Error in PlotManager::GetGraniittiPlot() scale factor is not set!" << std::endl;
        return hist;
    }

    if( rpCon < nRpConfigurations){
        hist = (TH1F*)graniittiPlots[rpCon]->Clone( "hGran" + hName + "_" + mUtil->particleName(part) + mUtil->rpConfigName(rpCon));
        if(scale)
            hist->Scale( getGraniittiSF(part, rpCon) );
    }else{
        hist = (TH1F*)graniittiPlots[IT]->Clone( "hGran" + hName + "_" + mUtil->particleName(part) );
        TH1F* tmp = (TH1F*)graniittiPlots[ET]->Clone( "tmpHist" );
        if(scale){
            hist->Scale( getGraniittiSF(part, IT) );
            tmp->Scale( getGraniittiSF(part, ET) );
            //cout<<"Scaling "<<hName<<" part "<< part <<" ET "<<getGraniittiSF(part, ET)<<" IT "<<getGraniittiSF(part, IT)<<endl;
        }
        hist->Add(tmp);
    }

    SetHistBcgStyle(hist);

    return hist;
}//GetGraniittiPlot

double PlotManager::getGraniittiSF(int part, int rpCon)
{
    return rpCon < nRpConfigurations ? graniittiScaleFactor[part][rpCon] : (graniittiScaleFactor[part][IT] + graniittiScaleFactor[part][ET])/2;
}//scaleFactor