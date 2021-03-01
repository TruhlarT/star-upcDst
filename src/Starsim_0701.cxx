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


TString summaryLabels[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("2 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString summaryLabels2[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("4 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString systemState[4] = { TString("TPC2t"), TString("TOF2trk"), TString("Q0"), TString("Excl")};
TString systemLabels[2] = { TString("2part"), TString("4part")};
TString sideLabel[nSides] = { TString("East"), TString("West")};
TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};


TH1F* hPt[2];
TH1F* hZ[2];
TH1F* hEta[2];
TH1F* hPhi[2];

TH3F* hEffPlot[3]; // 0 = after reconstruction, 1 = MC, 3 = data/MC and 4 different binning
TH3F* hEffPlotMy[3][4]; // 0 = after reconstruction, 1 = MC, 3 = data/MC 


TH1F* hNumberTrack[2]; // 1 = MC
TH1F* hNumberTOFmatchTrack;

TH1F* hMyEff[3][4];
TH1F* hRafalEff[3];


TFile *infile;
TFile *outfile;
TChain *upcChain;
StUPCEvent *upcEvt;

TTree *upcTree;

TTree *recTree, *bcgTree;
Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;

Int_t xMax, yMax, zMax;


void Init();
void Make();
TFile *CreateOutputTree(const string& out);
bool ConnectInput(int argc, char** argv);
void Clear();

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

    TString partLabel[2] = { TString(""), TString("_MC")};

    const double zvrtxBins[21] = {-200,-100,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,100,200};
    const double ptBins[21] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.6,2,3};
    const double etaBins[25] = {-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2};

    for (int i = 0; i < 2; ++i)
    {
        hPt[i] =  new TH1F("transMomentum" + partLabel[i], "transMomentum" + partLabel[i], 100, 0, 5);
        hEta[i] = new TH1F("eta" + partLabel[i], "eta" + partLabel[i], 120, -6, 6);
        hPhi[i] =  new TH1F("phi" + partLabel[i], "phi" + partLabel[i], 70, 0, 7);
        hZ[i] =  new TH1F("zVrtx" + partLabel[i], "zVrtx" + partLabel[i], 100, -200, 200);
    }





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

    hNumberTrack[0] = new TH1F("NumberTracks", "Number of Tracks", 40, -0.5, 39.5);
    hNumberTrack[1] = new TH1F("NumberTracksMC", "Number of Tracks in MC", 40, -0.5, 39.5);
    hNumberTOFmatchTrack = new TH1F("hNumberTOFmatchTrack", "Number of TOF matched tracks", 40, -0.5, 39.5);



    outfile->cd();
}

void Make(){
 

    Clear();

    Int_t mcTracks = upcEvt->getNumberOfMCParticles();
    hNumberTrack[1]->Fill(mcTracks);
    if(mcTracks < 6)
        return;
    // loop over TX MC tracks
    for (int i = 0; i < 6; ++i)
    {

        TParticle *part = upcEvt->getMCParticle(i);            
        hPt[1]->Fill(part->Pt());
        hEta[1]->Fill(part->Eta());
        hPhi[1]->Fill(part->Phi());
        hZ[1]->Fill(part->Vz());

        hEffPlot[1]->Fill(part->Vz(), part->Pt(), part->Eta());
        for (int iBins = 0; iBins < 4; ++iBins)
            hEffPlotMy[1][iBins]->Fill(part->Phi(), part->Pt(), part->Eta());
    }


    int nPrimtracks = 0;
    nTofTrks = 0;    
    for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
    {
    // get TPC track object
        const StUPCTrack* trk = upcEvt->getTrack(j);
        if( !trk->getFlag(StUPCTrack::kPrimary)) 
            continue;
        nPrimtracks++;

        //if(!trk->getFlag(StUPCTrack::kTof) )
        //    continue;

        nTofTrks++;
        if(nTofTrks > 6)
            continue;

        Double_t phiToMC = trk->getPhi();
        if( phiToMC < 0)
            phiToMC = 2*3.14159265359 + phiToMC;

        const StUPCVertex* vertex = trk->getVertex();

        hEffPlot[0]->Fill(vertex->getPosZ(), trk->getPt(), trk->getEta());
        for (int iBins = 0; iBins < 4; ++iBins)
            hEffPlotMy[0][iBins]->Fill(phiToMC, trk->getPt(), trk->getEta());

        hPt[0]->Fill(trk->getPt());
        hEta[0]->Fill(trk->getEta());
        hPhi[0]->Fill(phiToMC);
        hZ[0]->Fill(vertex->getPosZ());

    }
    hNumberTrack[0]->Fill(nPrimtracks);
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



void Clear()
{    
    nTofTrks = 0;

}



