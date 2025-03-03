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

            dir += "_con";

            hInvMassGranCon[iPart][iRpCon] = (TH1F*)graniittiFile->Get(dir + "/h_invMass");
            hDeltaPhiGranCon[iPart][iRpCon] = (TH1F*)graniittiFile->Get(dir + "/h_deltaPhi");

            if( !hInvMassGranCon[iPart][iRpCon] || !hDeltaPhiGranCon[iPart][iRpCon] )
               std::cerr << "Error in PlotManager::initGraniittiPlots() while retrieving graniitti plots!" << std::endl;


            graniittiScaleFactor[iPart][iRpCon] = -999.9;
        }
    }
}//initGraniittiPlots


TH1F* PlotManager::GetGraniittiPlot(int part, int rpCon, TString hName, double dataIntegral, bool scale, bool contiuum)
{
    TH1F *hist = nullptr;
    TH1F **graniittiPlots;

    bool scl = scale && getGraniittiSF(part, rpCon) < 0; 

    if( hName == "DeltaPhi"){
        if(contiuum)
            graniittiPlots = hDeltaPhiGranCon[part];
        else
            graniittiPlots = hDeltaPhiGran[part];
        if( rpCon < nRpConfigurations && scl){
            double tmp = dataIntegral / graniittiPlots[rpCon]->Integral("width");
            tmp = abs(tmp - 1.0) < 0.15 ? 1.0 : tmp;
            //double scndFactor = rpCon == IT ? graniittiScaleFactor[part][ET] : graniittiScaleFactor[part][IT];
            graniittiScaleFactor[part][rpCon] = abs(1.0 - tmp) < 0.1 ? 1.0 : tmp;
            //cout<<"Setting SF for "<<hName<<" part "<< mUtil->particleName(part) <<" rpCon "<<mUtil->rpConfigName(rpCon)<<" scale "<<scale<<" continuum "<< contiuum<< " "<<graniittiScaleFactor[part][rpCon]<<" / "<<getGraniittiSF(part, rpCon)<<endl;
            //cout<<" dataIntegral / graniittiPlots[rpCon]->Integral(width) "<< dataIntegral <<"  / "<< graniittiPlots[rpCon]->Integral("width") << endl;
        }
    }else if( hName == "InvMass"){
        if(contiuum)
            graniittiPlots = hInvMassGranCon[part];
        else
            graniittiPlots = hInvMassGran[part];
    }
    else{
        return hist;
    }

    if( getGraniittiSF(part, rpCon) < 0 && scale){
        std::cerr << "Error in PlotManager::GetGraniittiPlot() scale factor is not set!" << std::endl;
        return hist;
    }

    if( rpCon < nRpConfigurations){
        hist = (TH1F*)graniittiPlots[rpCon]->Clone( "hGran" + hName + "_" + mUtil->particleName(part) + mUtil->rpConfigName(rpCon));
        if( scale )
            hist->Scale( getGraniittiSF(part, rpCon) );
    }else{
        hist = (TH1F*)graniittiPlots[IT]->Clone( "hGran" + hName + "_" + mUtil->particleName(part) );
        TH1F* tmp = (TH1F*)graniittiPlots[ET]->Clone( "tmpHist" );
        if( scale ){
            hist->Scale( getGraniittiSF(part, IT) );
            tmp->Scale( getGraniittiSF(part, ET) );
        }
        //cout<<"Scaling "<<hName<<" part "<< mUtil->particleName(part) <<" ET "<<getGraniittiSF(part, ET)<<" IT "<<getGraniittiSF(part, IT)<<endl;
        hist->Add(tmp);
    }

    short lineStyles[] = {1, 2, 5, 6}; 
    short lnColor[] = { bckgColor, 1, 9, 6};
    short lnStyle = hName == "InvMass" ? 0 : (rpCon == ET ? 2 : 0);
    if(contiuum)
        lnStyle++;

    SetHistBcgStyle(hist, lnColor[lnStyle]);
    hist->SetLineStyle(lineStyles[lnStyle]);
    hist->SetLineWidth(2*lineWidth); 

    return hist;
}//GetGraniittiPlot

double PlotManager::getGraniittiSF(int part, int rpCon)
{
    return rpCon < nRpConfigurations ? graniittiScaleFactor[part][rpCon] : (graniittiScaleFactor[part][IT] + graniittiScaleFactor[part][ET])/2;
}//scaleFactor