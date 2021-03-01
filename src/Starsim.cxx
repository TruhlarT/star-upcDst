// Run by: ./Starsim file.list
// e.g. ./Starsim /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list
// or you can open just root file
// ./Starsim /star/data01/pwg_tasks/upc02/Part9/18143045/18143045.root
// or you can open n-th root file in file.list
// ./Starsim file.list index 


// c++ headers
#include <iostream>
#include <string>     // std::string, std::to_string
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <TH3.h>
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TGraphAsymmErrors.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>

// picoDst headers
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"

using namespace std;

enum {kAll = 1, kTOF2t, kSameVrtx, kTotCH0, kMissPt, kMaxCount};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};

const double particleMass[nParticles] = { 0.13957, 0.49367, 0.93827}; // GeV /c^2
const double speedOfLight = 299792458; // m/s
const double beamMomentum = 254.9; // GeV
const double PI = TMath::Pi();


TH1F* hPt[3];
TH1F* hZ[3];
TH1F* hEta[3];
TH1F* hPhi[3];
TH1F* hTestPhiOfMatched;

TH1F* hDeltaPhi[2];
TH1F* hDeltaEta[2];
TH2F* hDeltaCorr[2];
TH1F* hDelta[2];
TH1F* hDeltaMatched;

TH1F* hDcaXY[2];
TH1F* hDcaZ[2];
TH1F* hNhitsFit[2];
TH1F* hNhitsDEdx[2];
TH1F* hNSigma[2];

TH1F* hPolarTheta;
TH2F* hThetaCom;

TH3F* hEffPlot[3]; // 0 = after reconstruction, 1 = MC, 3 = data/MC and 4 different binning
TH3F* hEffPlotMy[3][4]; // 0 = after reconstruction, 1 = MC, 3 = data/MC 


TH1F* hNumberTrack[3]; // 1 = MC, 2 = matched with MC
TH1F* hNumberMatchTrack;
TH1F* hNumberTOFmatchTrack;

TH1F* hMyEff[3][4];
TH1F* hRafalEff[3];
TH1F* hPhiMCTest[6];

TFile *infile;
TFile *outfile;
TChain *upcChain;
StUPCEvent *upcEvt;

TTree *upcTree;

TTree *recTree, *bcgTree;
Int_t nTracks, totalCharge, nTofTrks, nPrimtracks, matchedTracks; 
UInt_t runNumber;

Int_t xMax, yMax, zMax;
Double_t phiMC[6];
Double_t etaMC[6];
Double_t zMC[6];
Double_t ptMC[6];


void Init();
void Make();
TFile *CreateOutputTree(const string& out);
bool ConnectInput(int argc, char** argv);
bool IsGoodTrack(const StUPCTrack* trk);
int IsMatchedToMC(const StUPCTrack* trk);
void FillTrackHist(const StUPCTrack* trk, int i);

void TestMCMatching();

//_____________________________________________________________________________
int main(int argc, char** argv) {
    //connect input file(s)
    if(!ConnectInput(argc, argv))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }
    //open output file
    outfile = CreateOutputTree("AnalysisOutput.root"); 
    if(!outfile) 
    {
        cout << "Can not open output file." << endl; 
        return -1;
    }

    xMax = 20;
    yMax = 20;
    zMax = 24;
    Init();

    //ask for number of events
    Long64_t nev = upcTree->GetEntries();
    cout<<"Going to process "<<nev<<" events"<<endl;
    //return 0;
    //event loop
    //nev = 1000;
    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        if(iev%10000 == 0)
            cout<<"Processing "<<iev<<" event"<<endl;
        upcTree->GetEntry(iev); 
        Make();
    } 

    hEffPlot[2] = (TH3F*)hEffPlot[0]->Clone("effRafal");
    hEffPlot[2]->Divide(hEffPlot[1]);

    for (int iBins = 0; iBins < 4; ++iBins)
    {
        hEffPlotMy[2][iBins] = (TH3F*)hEffPlotMy[0][iBins]->Clone(Form("effMy_%i",iBins));
        hEffPlotMy[2][iBins]->Divide(hEffPlotMy[1][iBins]); 
    }

    Int_t x, y, z;

    x = 1;
    y = 1;
    z = 1;
    for (int i = 1; i < xMax*yMax*zMax +1; ++i)
    {
        hRafalEff[0]->SetBinContent(i, hEffPlot[2]->GetBinContent( x, y, z));
        hRafalEff[1]->SetBinContent(i, hEffPlot[0]->GetBinContent( x, y, z));
        hRafalEff[2]->SetBinContent(i, hEffPlot[1]->GetBinContent( x, y, z));

        x++;
        if(x == xMax+1)
        {
            x = 1;
            y++;
            if(y == yMax+1)
            {
                y = 1;
                z++;
            }
        }    
    }

    Int_t xMaxNew, yMaxNew, zMaxNew;
    for (int iBins = 0; iBins < 4; ++iBins)
    {
        x = 1;
        y = 1;
        z = 1;
        xMaxNew = xMax - iBins*4;
        yMaxNew = yMax - iBins*4;
        zMaxNew = zMax - iBins*4; 
        for (int i = 1; i < xMaxNew*yMaxNew*zMaxNew +1; ++i)
        {

            hMyEff[0][iBins]->SetBinContent(i, hEffPlotMy[2][iBins]->GetBinContent( x, y, z));
            hMyEff[1][iBins]->SetBinContent(i, hEffPlotMy[0][iBins]->GetBinContent( x, y, z));
            hMyEff[2][iBins]->SetBinContent(i, hEffPlotMy[1][iBins]->GetBinContent( x, y, z));

            x++;
            if(x == xMaxNew+1)
            {
                x = 1;
                y++;
                if(y == yMaxNew+1)
                {
                    y = 1;
                    z++;
                }
            }    
        }
    }
    outfile->Write();
    outfile->Close();
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

void Init(){

    TString partLabel[2] = { TString("_allTracks"), TString("_goodTracks")};

    const double zvrtxBins[21] = {-200,-100,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,100,200};
    const double ptBins[21] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.6,2,3};
    const double etaBins[25] = {-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2};

    for (int i = 0; i < 2; ++i)
    {
        hPt[i] =  new TH1F("transMomentum" + partLabel[i], "transMomentum" + partLabel[i], 100, 0, 5);
        hEta[i] = new TH1F("eta" + partLabel[i], "eta" + partLabel[i], 120, -6, 6);
        hPhi[i] =  new TH1F("phi" + partLabel[i], "phi" + partLabel[i], 70, 0, 7);
        hZ[i] =  new TH1F("zVrtx" + partLabel[i], "zVrtx" + partLabel[i], 100, -200, 200);
        hDcaXY[i] = new TH1F("DcaXY" + partLabel[i], "hDcaXY" + partLabel[i], 100, -0.5, 5);
        hDcaZ[i] = new TH1F("DcaZ" + partLabel[i], "hDcaZ" + partLabel[i], 100, -5, 5);
        hNhitsFit[i] = new TH1F("NhitsFit" + partLabel[i], "hNhitsFit" + partLabel[i], 60, -0.5, 59.5);
        hNhitsDEdx[i] = new TH1F("NhitsDEdx" + partLabel[i], "hNhitsDEdx" + partLabel[i], 60, -0.5, 59.5);
        hNSigma[i] = new TH1F("NSigma" + partLabel[i], "hNSigma" + partLabel[i], 100, -7, 7);
    }

    hPt[2] =  new TH1F("transMomentum_MC", "transMomentum_MC", 100, 0, 5);
    hEta[2] = new TH1F("eta_MC", "eta_MC", 120, -6, 6);
    hPhi[2] =  new TH1F("phi_MC", "phi_MC", 70, 0, 7);
    hZ[2] =  new TH1F("zVrtx_MC", "zVrtx_MC", 100, -200, 200);



    hEffPlot[0] = new TH3F("effDataSMC", "number of particles from smeared MC", 20, zvrtxBins, 20, ptBins, 24, etaBins); // x = z_vertex, y = p_T, z = eta
    hEffPlot[1] = new TH3F("effDataMC", "number of particles from MC", 20, zvrtxBins, 20, ptBins, 24, etaBins); // x = z_vertex, y = p_T, z = eta
    hEffPlot[2] = new TH3F("effRafal", "efficiencies smeared/MC", 20, zvrtxBins, 20, ptBins, 24, etaBins); // x = z_vertex, y = p_T, z = eta

    for (int i = 0; i < 4; ++i)
    {
        hEffPlotMy[0][i] = new TH3F(Form("effDataSMCmy_%i", i), "number of particles from smeared MC", xMax - i*4, 0, 6.3, yMax - i*4, 0, 3, zMax - i*4, -1.2, 1.2); // x = phi, y = p_T, z = eta
        hEffPlotMy[1][i] = new TH3F(Form("effDataMCmy_%i", i), "number of particles from MC", xMax - i*4, 0, 6.3, yMax - i*4, 0, 3, zMax - i*4, -1.2, 1.2); // x = phi, y = p_T, z = eta
        hEffPlotMy[2][i] = new TH3F(Form("effMy_%i", i), "efficiencies smeared/MC", xMax - i*4, 0, 6.3, yMax - i*4, 0, 3, zMax - i*4, -1.2, 1.2); // x = phi, y = p_T, z = eta
        
        hMyEff[0][i] = new TH1F(Form("myEff_%i", i), "my efficiencies phi:pT:eta", (xMax-i*4)*(yMax-i*4)*(zMax-i*4), -0.5, (xMax-i*4)*(yMax-i*4)*(zMax-i*4) -0.5);
        hMyEff[1][i] = new TH1F(Form("mySMC_%i", i), "number of SMC particles phi:pT:eta", (xMax-i*4)*(yMax-i*4)*(zMax-i*4), -0.5, (xMax-i*4)*(yMax-i*4)*(zMax-i*4) -0.5);
        hMyEff[2][i] = new TH1F(Form("myMC_%i", i), "number of MC particles phi:pT:eta", (xMax-i*4)*(yMax-i*4)*(zMax-i*4), -0.5, (xMax-i*4)*(yMax-i*4)*(zMax-i*4) -0.5);
    }


    hRafalEff[0] = new TH1F("rafalEff", "my efficiencies z:pT:eta", xMax*yMax*zMax, -0.5, xMax*yMax*zMax -0.5);
    hRafalEff[1] = new TH1F("rafalSMC", "number of SMC particles z:pT:eta", xMax*yMax*zMax, -0.5, xMax*yMax*zMax -0.5);
    hRafalEff[2] = new TH1F("rafalMC", "number of MC particles z:pT:eta", xMax*yMax*zMax, -0.5, xMax*yMax*zMax -0.5);

    hNumberTrack[0] = new TH1F("NumberTracks", "Number of primary tracks", 40, -0.5, 39.5);
    hNumberTrack[1] = new TH1F("NumberTracksMC", "Number of Tracks in MC", 40, -0.5, 39.5);
    hNumberTrack[2] = new TH1F("NumberMtchTracks", "Number of matched Tracks w MC", 40, -0.5, 39.5);
    hNumberTOFmatchTrack = new TH1F("hNumberTOFmatchTrack", "Number of TOF matched tracks", 40, -0.5, 39.5);
    hNumberMatchTrack = new TH1F("hNumberMatchTrack", "Number of matched tracks MC for each reco track", 10, -0.5, 9.5);


    hDeltaEta[0] = new TH1F("deltaEta", "(eta_track - eta_MC)^2", 1000, 0.0, 2.5);
    hDeltaPhi[0] = new TH1F("deltaPhi", "(Phi_track - Phi_MC)^2", 1000, 0.0, 2.5);
    hDelta[0] = new TH1F("delta", "deltaPhiSq + deltaEtaSq", 10000, 0.0, 5);
    hDeltaMatched = new TH1F("deltaMatched", "deltaPhiSq + deltaEtaSq for matched tracks", 1000, 0.0, 0.03);
    hDeltaCorr[0] = new TH2F("deltaCorr", "deltaEtaSq vs deltaPhiSq", 1000, 0.0, 1.5, 1000, 0.0, 1.5);

    hDeltaEta[1] = new TH1F("deltaEta_wide", "(eta_track - eta_MC)^2", 1000, 0.0, 25);
    hDeltaPhi[1] = new TH1F("deltaPhi_wide", "(Phi_track - Phi_MC)^2", 1000, 0.0, 50);
    hDelta[1] = new TH1F("delta_wide", "deltaPhiSq + deltaEtaSq", 1000, 0.0, 50);
    hDeltaCorr[1] = new TH2F("deltaCorr_wide", "deltaEtaSq vs deltaPhiSq", 1000, 0.0, 15, 1000, 0.0, 15);

    hPolarTheta = new TH1F("PolarTheta", "PolarTheta", 100, -5, 5);
    hThetaCom = new TH2F("thetaCom","theta Comparison", 100, -6, 6, 100, -6, 6);

    for (int i = 0; i < 6; ++i)
    {
        hPhiMCTest[i] = new TH1F(Form("phiMCTest_%i",i), Form("phiMCTest_%i",i), 100, 0, 6.5); 
    }
    hTestPhiOfMatched = new TH1F("hTestPhiOfMatched" , "hTestPhiOfMatched", 100, 0, 6.5);
    outfile->cd();
}

void Make(){
 
    Int_t mcTracks = upcEvt->getNumberOfMCParticles();
    hNumberTrack[1]->Fill(mcTracks);
    if(mcTracks < 6)
        return;
    // loop over TX MC tracks
    for (int i = 0; i < 6; ++i)
    {

        TParticle *part = upcEvt->getMCParticle(i); 
        phiMC[i] = part->Phi();
        hPhiMCTest[i]->Fill(phiMC[i]);
        if(phiMC[i] < (PI/3)*i || phiMC[i] > (PI/3)*(i+1))
            return;

        etaMC[i] = part->Eta();
        ptMC[i] = part->Pt();
        zMC[i] = part->Vz();           
        hPt[2]->Fill(part->Pt());
        hEta[2]->Fill(part->Eta());
        hPhi[2]->Fill(part->Phi());
        hZ[2]->Fill(part->Vz());

        hPolarTheta->Fill(part->Theta());
        double theta = part->Theta();
        double zDistance = 210;
        if( theta > PI/2)
        {
            theta = PI - theta;
            zDistance = -zDistance;
        }
        hThetaCom->Fill(theta, TMath::ATan(220/(zDistance - part->Vz())));
        //if(part->Pt() < 0.15 || theta < TMath::Abs(TMath::ATan(220/(zDistance - part->Vz()))))
        //    continue;
        hEffPlot[1]->Fill(part->Vz(), part->Pt(), part->Eta());
        for (int iBins = 0; iBins < 4; ++iBins)
            hEffPlotMy[1][iBins]->Fill(part->Phi(), part->Pt(), part->Eta());
    }
//    TestMCMatching();

    nPrimtracks = 0;
    nTofTrks = 0;   
    int indexMC, nMatched;
    nMatched = 0; 
    for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
    {
    // get TPC track object
        const StUPCTrack* trk = upcEvt->getTrack(j);
        FillTrackHist(trk,0);

        if(!IsGoodTrack(trk))
            continue;

        const StUPCVertex* vertex = trk->getVertex();
        if(vertex->getPosZ() > 80 || vertex->getPosZ() < -80)
            continue; 

        matchedTracks = 0;
        indexMC = IsMatchedToMC(trk);
        hNumberMatchTrack->Fill(matchedTracks);
        if(indexMC < 0)
            continue;
        nMatched++;
        FillTrackHist(trk,1);

        Double_t phiToMC = trk->getPhi();
        if( phiToMC < 0)
            phiToMC = 2*3.14159265359 + phiToMC;

        hEffPlot[0]->Fill(zMC[indexMC], ptMC[indexMC], etaMC[indexMC]);
        for (int iBins = 0; iBins < 4; ++iBins)
            hEffPlotMy[0][iBins]->Fill(phiMC[indexMC], ptMC[indexMC], etaMC[indexMC]);


    }
    hNumberTrack[0]->Fill(nPrimtracks);
    hNumberTrack[2]->Fill(nMatched);
    hNumberTOFmatchTrack->Fill(nTofTrks);

    //recTree->Fill();

}



//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    recTree = new TTree("recTree", "recTree");

    return outputFile;

}//CreateOutputTree

bool ConnectInput(int argc, char** argv) {
    int fileId = -1;
    upcTree = 0x0;

    if( argc == 2 || argc == 3)
    {
        const string& input = argv[1];
        if(input.find(".root") != string::npos && argc == 2)
        {
            cout << "Input from root file: "<< input << endl;
            infile = TFile::Open(input.c_str(), "read");
            if(!infile)
            {
                cout<< "Couldn't open input root file..."<<endl;
                return false;
            } 
            upcTree = dynamic_cast<TTree*>( infile->Get("mUPCTree") );
        }else if(input.find(".list") != string::npos )
        {
            cout << "Input from chain" << endl;
            upcChain = new TChain("mUPCTree");
            ifstream instr(input.c_str());
            if (!instr.is_open())
            {
                cout<< "Couldn't open: "<<input.c_str()<<endl;
                return false;
            }

            string line;
            int lineId=0;
            if(argc == 3)
                fileId = atoi(argv[2]);

            while(getline(instr, line)) 
            {
                if(fileId==lineId || fileId== -1)
                {
                    upcChain->AddFile(line.c_str());
                    infile = TFile::Open(line.c_str(), "read");
                    if(!infile)
                    {
                        cout<< "Couldn't open: "<<line.c_str()<<endl;
                        return false;
                    } 
                }
                lineId++;
            }
            instr.close();
            upcTree = dynamic_cast<TTree*>( upcChain );
        }
    }
    
    if(!upcTree) 
        return false;

    upcTree->SetBranchAddress("mUPCEvent", &upcEvt); 

    return true;
}//ConnectInput



bool IsGoodTrack(const StUPCTrack* trk)
{
    if( !trk->getFlag(StUPCTrack::kPrimary)) 
        return false;

    nPrimtracks++;
    if( !trk->getFlag(StUPCTrack::kTof)) 
        return false;

    nTofTrks++;
    if( trk->getPt() < 0.4 ||
        trk->getDcaXY() > 1.5 ||
        (trk->getDcaZ() > 1 || trk->getDcaZ() < -1) ||
        (trk->getEta() > 1.0 || trk->getEta() < -1.0) ||
        trk->getNhitsFit() < 25 ||
        trk->getNhitsDEdx() < 15 ||
        (trk->getNSigmasTPCProton() < -3 || trk->getNSigmasTPCProton() > 3))
        return false;
 

    return true;
}


int IsMatchedToMC(const StUPCTrack* trk)
{
    double deltaPhi, deltaEta, phiToMC, delta;
    double deltaPhiSq;
    int matchedIndex = -1;
    for (int i = 0; i < 6; ++i)
    {
        deltaEta = (trk->getEta() - etaMC[i])*(trk->getEta() - etaMC[i]);
        phiToMC = trk->getPhi();
        if( phiToMC < 0)
            phiToMC = 2*PI + phiToMC;
        deltaPhi = phiToMC - phiMC[i];
        if(TMath::Abs(deltaPhi) > PI)
            deltaPhi = 2*PI - TMath::Abs(deltaPhi);
        deltaPhiSq = deltaPhi*deltaPhi; 
        delta = deltaPhiSq + deltaEta;
        
        hDeltaPhi[0]->Fill(deltaPhiSq);
        hDeltaEta[0]->Fill(deltaEta);
        hDelta[0]->Fill(delta);
        hDeltaCorr[0]->Fill(deltaEta,deltaPhiSq);
        hDeltaPhi[1]->Fill(deltaPhiSq);
        hDeltaEta[1]->Fill(deltaEta);
        hDelta[1]->Fill(delta);
        hDeltaCorr[1]->Fill(deltaEta,deltaPhiSq);
        if(delta < 0.0225)
        {
            matchedIndex = i;
            matchedTracks++;
            if(matchedTracks == 2)
            {
                hTestPhiOfMatched->Fill(phiToMC);
                hDeltaMatched->Fill(delta);
                /*cout<<"Reco: "<< phiToMC <<" : "<< trk->getEta() <<endl;
                cout<<"MC1: "<< phiMC[0] <<" : "<< etaMC[0] <<endl;
                cout<<"MC2: "<< phiMC[1] <<" : "<< etaMC[1] <<endl;
                cout<<"MC3: "<< phiMC[2] <<" : "<< etaMC[2] <<endl;
                cout<<"MC4: "<< phiMC[3] <<" : "<< etaMC[3] <<endl;
                cout<<"MC5: "<< phiMC[4] <<" : "<< etaMC[4] <<endl;
                cout<<"MC6: "<< phiMC[5] <<" : "<< etaMC[5] <<endl;*/
            }
        }
    }
    return matchedIndex;
}

void FillTrackHist(const StUPCTrack* trk, int i)
{
    double phiToMC = trk->getPhi();
    if( phiToMC < 0)
        phiToMC = 2*PI + phiToMC;
    const StUPCVertex* vertex = trk->getVertex();

    hDcaXY[i]->Fill(trk->getDcaXY());
    hDcaZ[i]->Fill(trk->getDcaZ());
    hNhitsFit[i]->Fill(trk->getNhitsFit());
    hNhitsDEdx[i]->Fill(trk->getNhitsDEdx());
    hNSigma[i]->Fill(trk->getNSigmasTPCPion());
    hPt[i]->Fill(trk->getPt());
    hEta[i]->Fill(trk->getEta());
    hPhi[i]->Fill(phiToMC);
    hZ[i]->Fill(vertex->getPosZ());
}


void TestMCMatching()
{

    double deltaPhi, deltaEta, phiToMC, delta;
    double deltaPhiSq;
    int mtch = 0;
    for (int j = 0; j < 6; ++j)
    {

        for (int i = 0; i < 6; ++i)
        {
            deltaEta = (etaMC[j] - etaMC[i])*(etaMC[j] - etaMC[i]);
            deltaPhi = phiMC[j] - phiMC[i];
            deltaPhiSq = deltaPhi*deltaPhi; 
            delta = deltaPhiSq + deltaEta;
            
            if(delta < 0.001)
            {
                mtch++;
            }
        }
    }
}