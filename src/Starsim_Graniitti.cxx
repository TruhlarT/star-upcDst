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

TH2F *hPtCor, *hEtaCor, *hPhiCor;
TH1F* hPt[2];
TH1F* hEta[2];
TH1F* hPhi[2];
TH1F* hMCPDG;

TH1I* hAnalysisFlow; // control plot
TH1I* hAnalysisFlow2; // control plot for 4 pion state
TH1F* hNumberTrack;
TH1F* hNumberTOFmatchTrack;

TH1D* hInvMassPions, *hInvMassPionsMC;

 
TH2D* hdEdxVsqP[4][2];
TH2D* hdEdx[4][2];

TFile *infile;
TFile *outfile;
TChain *upcChain;
StUPCEvent *upcEvt;

TTree *upcTree;

TTree *recTree, *bcgTree;
Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;
Double_t VPDSumEast, VPDSumWest, VPDTimeDiff;


Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t deltaTOF, mSquared, pairRapidity;
Double_t deltaDeltaTOF[nParticles];
Double_t deltaTOFExpected[nParticles];

/////////////////////////////////

Bool_t fourPiState;
UInt_t BBCSmall[nSides]; // BBC truncated sum, small tiles
UInt_t BBCLarge[nSides]; // BBC truncated sum, large tiles 

Double_t dEdx[4];
Double_t momentum[4];
Double_t transMomentum[4];
Double_t TOFtime[4];
Double_t TOFlength[4];
Double_t charge[4];
Double_t nSigmaTPC[nParticles][4];
Double_t vertexId[4];
Double_t vertexesZ[4];
Double_t DcaXY[4];
Double_t DcaZ[4];
Double_t NhitsFit[4];
Double_t NhitsDEdx[4];
Double_t Eta[4];
Double_t Phi[4];
Double_t Chi2[4];


void Init();
void Make();
void FillPlots(int state, int configuration, int isFourPiState );
TFile *CreateOutputTree(const string& out);
bool ConnectInput(int argc, char** argv);
void Clear();
void CalculatePID();

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

    Init();

    //ask for number of events
    Long64_t nev = upcTree->GetEntries();
    cout<<"Proccesing "<<nev<<" events"<<endl;
    //return 0;
    //event loop
    //nev = 1000;
    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        upcTree->GetEntry(iev); 
        Make();
    } 

    TH1D* hEff = new TH1D("efficiencies", "efficiencies (accepted/MC)", 64, 0.3, 3.5);
    hEff = (TH1D*)hInvMassPions->Clone("hist");
    hEff->Divide(hInvMassPionsMC);
   /* TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

    TGraphAsymmErrors err;
    err.Divide(hInvMassPions, hInvMassPionsMC,"pois");
    err.SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; Accepted / Graniitti");
    //err.GetYaxis()->SetRangeUser(0,3);
    err.SetMarkerColor(4);
    err.SetMarkerSize(1);
    err.SetMarkerStyle(20);
    err.SetLineColor(4);
    err.SetLineStyle(1);
    err.SetLineWidth(1);
    err.Draw("AP");
    newCanvas->Update();
    newCanvas->Write("eff_Pions");

    //close the outputs
    newCanvas->Close();
    */outfile->Write();
    outfile->Close();
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

void Init(){

    TString partLabel[2] = { TString(""), TString("_MC")};

    for (int i = 0; i < 2; ++i)
    {
        hPt[i] =  new TH1F("transMomentum" + partLabel[i], "transMomentum" + partLabel[i], 100, 0, 5);
        hEta[i] = new TH1F("eta" + partLabel[i], "eta" + partLabel[i], 120, -6, 6);
        hPhi[i] =  new TH1F("phi" + partLabel[i], "phi" + partLabel[i], 70, 0, 7);
    }

    hPtCor =  new TH2F("transMomentumCor", "transMomentumCor (x-axis is MC)", 100, 0, 5, 100, 0, 5);
    hEtaCor = new TH2F("etaCor", "etaCor (x-axis is MC)", 120, -6, 6, 120, -6, 6);
    hPhiCor =  new TH2F("phiCor", "phiCor (x-axis is MC)", 70, 0, 7, 70, 0, 7);
    hMCPDG =  new TH1F("MCPDG", "MCPDG (x-axis is MC)", 1000, -500, 500);

    hInvMassPionsMC = new TH1D("invMassPionsMC", "inv. mass pions MC", 64, 0.3, 3.5);
    hInvMassPions = new TH1D("invMassPion", "inv. mass pions", 64, 0.3, 3.5);


    hAnalysisFlow = new TH1I("AnalysisFlow", "CutsFlow", kMaxCount-1, 1, kMaxCount);
    hAnalysisFlow2 = new TH1I("AnalysisFlow2", "CutsFlow", kMaxCount-1, 1, kMaxCount); 
    for(int tb=1; tb<kMaxCount; ++tb) 
        hAnalysisFlow->GetXaxis()->SetBinLabel(tb, summaryLabels[tb-1]);

    for(int tb=1; tb<kMaxCount; ++tb) 
        hAnalysisFlow2->GetXaxis()->SetBinLabel(tb, summaryLabels2[tb-1]);


    hNumberTrack = new TH1F("NumberTracks", "Number of Tracks in RP", 40, -0.5, 39.5);
    hNumberTOFmatchTrack = new TH1F("hNumberTOFmatchTrack", "Number of TOF matched tracks", 40, -0.5, 39.5);

    for(int i=0; i<4;++i){

        for( int state = 0; state < 2; ++state)
        {
            hdEdxVsqP[i][state] = new TH2D("dEdxVsqP_"+systemState[i%4]+"_"+"Combi" + systemLabels[state],"dE/dx Vs #frac{q}{e} P",200,-2,2,100,0,20);
            hdEdx[i][state] = new TH2D("dEdx_"+systemState[i%4]+"_"+"Combi" + systemLabels[state],"#log_{10} dE/dx [keV/cm] vs #log_{10} p [GeV/c]",200,-1,1,100,0.1,2);
        }
    }


    outfile->cd();
}

void Make(){
    hAnalysisFlow->Fill(kAll);
    hAnalysisFlow2->Fill(kAll);

    Clear();

    Double_t pTMC[2];
    Double_t etaMC[2];
    Double_t phiMC[2];

    // loop over TX MC tracks
    for (int i = 0; i < 2; ++i)
    {
        TParticle *part = upcEvt->getMCParticle(i);            
        hPt[1]->Fill(part->Pt());
        hEta[1]->Fill(part->Eta());
        hPhi[1]->Fill(part->Phi());
        hMCPDG->Fill(part->GetPdgCode());

        int iCharge = 0;
        if(part->GetPdgCode() < 0)
            iCharge = 1;

        pTMC[iCharge] = part->Pt();
        etaMC[iCharge] = part->Eta();
        phiMC[iCharge] = part->Phi();
       
    }

    if(etaMC[0] < 0.7 && etaMC[0] > -0.7 && pTMC[0] > 0.2 && etaMC[1] < 0.7 && etaMC[1] > -0.7 && pTMC[1] > 0.2)
    {
        TLorentzVector trkLVector, centralVector;
        trkLVector.SetPtEtaPhiM( pTMC[0], etaMC[0], phiMC[0], particleMass[Pion]);
        centralVector = trkLVector;
        trkLVector.SetPtEtaPhiM( pTMC[1], etaMC[1], phiMC[1], particleMass[Pion]);
        centralVector += trkLVector;
        hInvMassPionsMC->Fill(centralVector.M());
    }


    int nPrimtracks = 0;
    for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
    {
    // get TPC track object
        const StUPCTrack* trk = upcEvt->getTrack(j);
        if( !trk->getFlag(StUPCTrack::kPrimary)) 
            continue;
        nPrimtracks++;
        if(nPrimtracks == 3)
            break;

        Double_t phiToMC = trk->getPhi();
        if( phiToMC < 0)
            phiToMC = 2*3.14159265359 + phiToMC;

        hPt[0]->Fill(trk->getPt());
        hEta[0]->Fill(trk->getEta());
        hPhi[0]->Fill(phiToMC);

        int iCharge = 0;
        if(trk->getCharge() < 0)
            iCharge = 1;

        hPtCor->Fill(pTMC[iCharge], trk->getPt());
        hEtaCor->Fill(etaMC[iCharge], trk->getEta());
        hPhiCor->Fill(phiMC[iCharge], phiToMC);


    }

    TLorentzVector centralTracksTotalFourMomentum[nParticles];
    //TLorentzVector tmpTLorentzVector;
    // loop over all TPC tracks
    for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
    {
    // get TPC track object
        const StUPCTrack* trk = upcEvt->getTrack(j);

        TLorentzVector trkLVector[nParticles];
        for (int iPart = 0; iPart < nParticles; ++iPart)
            trk->getLorentzVector(trkLVector[iPart], particleMass[iPart]);

        if( !trk->getFlag(StUPCTrack::kPrimary)) 
            continue;


        if(!trk->getFlag(StUPCTrack::kTof) )
        {
            continue;
        } 

        if(nTofTrks>3)
        {
            nTofTrks++;
            continue;
        }

        //Calculate deltaPhi for 2 hadron state
        //if(nTofTrks == 0) tmpTLorentzVector = trkLVector[Pion];
        //if(nTofTrks == 1) deltaPhi = trkLVector[Pion].Angle(tmpTLorentzVector.Vect());

        // read basic information about the track     
        dEdx[nTofTrks] = trk->getDEdxSignal()*1000000;
        momentum[nTofTrks] = trkLVector[0].P();
        transMomentum[nTofTrks] = trk->getPt();
        charge[nTofTrks] = trk->getCharge();
        TOFtime[nTofTrks] = trk->getTofTime();
        TOFlength[nTofTrks] = trk->getTofPathLength();
        DcaXY[nTofTrks] = trk->getDcaXY();
        DcaZ[nTofTrks] = trk->getDcaZ();
        NhitsFit[nTofTrks] = trk->getNhitsFit();
        NhitsDEdx[nTofTrks] = trk->getNhitsDEdx();
        Eta[nTofTrks] = trk->getEta();
        Phi[nTofTrks] = trk->getPhi();
        Chi2[nTofTrks] = trk->getChi2();
        vertexId[nTofTrks] = trk->getVertexId();
        vertexesZ[nTofTrks] = trk->getVertex()->getPosZ();
        totalCharge += static_cast<int>( trk->getCharge() );
        nSigmaTPC[Pion][nTofTrks] = trk->getNSigmasTPCPion();
        nSigmaTPC[Kaon][nTofTrks] = trk->getNSigmasTPCKaon();
        nSigmaTPC[Proton][nTofTrks] = trk->getNSigmasTPCProton();

        for (int iPart = 0; iPart < nParticles; ++iPart)
            centralTracksTotalFourMomentum[iPart] += trkLVector[iPart];
        
        ++nTofTrks;
    }



    hNumberTOFmatchTrack->Fill(nTofTrks);
    pairRapidity = centralTracksTotalFourMomentum[0].Rapidity();
    for (int iPart = 0; iPart < nParticles; ++iPart)
        invMass[iPart] = centralTracksTotalFourMomentum[iPart].M();

    FillPlots(0,0,0);
    FillPlots(0,0,1);
    fourPiState = false;
    if( nTofTrks==4)
    {
        fourPiState = true;
        hAnalysisFlow2->Fill(kTOF2t);
    }else if (nTofTrks!=2)
    {
        return;
    }


    if(!fourPiState)
        hAnalysisFlow->Fill(kTOF2t);

    if( vertexId[0] != vertexId[1])
        return;

    if(fourPiState)
    {
        if(vertexId[2] != vertexId[3] || vertexId[2] != vertexId[1])
            return;
        hAnalysisFlow2->Fill(kSameVrtx);
        FillPlots(1,0,1);
    }else
    {
        hAnalysisFlow->Fill(kSameVrtx);
        FillPlots(1,0,0);
    }


    for(int iPart = 0; iPart < nParticles; ++iPart)
        chiPair[iPart] = nSigmaTPC[iPart][0]*nSigmaTPC[iPart][0] + nSigmaTPC[iPart][1]*nSigmaTPC[iPart][1];



    runNumber = upcEvt->getRunNumber();
    BBCLarge[East] = upcEvt->getBBCLargeEast();
    BBCSmall[East] = upcEvt->getBBCSmallEast();
    BBCLarge[West] = upcEvt->getBBCLargeWest();
    BBCSmall[West] = upcEvt->getBBCSmallWest();
    VPDTimeDiff = upcEvt->getVPDTimeDiff();
    VPDSumWest = upcEvt->getVPDSumWest();
    VPDSumEast = upcEvt->getVPDSumEast();

    if(!fourPiState)
       CalculatePID();

    if(totalCharge != 0){
        bcgTree->Fill();
        return;
    }

    if(fourPiState)
    {
        hAnalysisFlow2->Fill(kTotCH0);
        FillPlots(2,0,1);
    } else
    {
        hAnalysisFlow->Fill(kTotCH0);
        FillPlots(2,0,0);      
    }

    if(fourPiState)
    {
        hAnalysisFlow2->Fill(kMissPt);
        FillPlots(3,0,1);
    } else
    {
        hAnalysisFlow->Fill(kMissPt);
        FillPlots(3,0,0);      
    }

    recTree->Fill();



    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && NhitsFit[0] >=25 && 
                NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 && 
                DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && 
                DcaXY[0] < 1.5 && DcaXY[1] < 1.5 && Eta[0] > -0.7 && Eta[0] < 0.7 && 
                Eta[1] > -0.7 && Eta[1] < 0.7 && transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
    {
        //PID condition
        if(chiPair[Pion] > 9 && chiPair[Kaon] > 9 && chiPair[Proton] < 9 && mSquared > 0.6) // it is... proton!
        {

        }
        else if(chiPair[Pion] > 9 && chiPair[Kaon] < 9 && chiPair[Proton] > 9 && mSquared > 0.15) // it is... kaon!
        {

        }
        else if( chiPair[Pion] < 12) // it is... pion!
        {
            hInvMassPions->Fill(invMass[Pion]);
        }
    }
}



//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    recTree = new TTree("recTree", "recTree");

// PID and some quality event info
    recTree->Branch("deltaTOF", &deltaTOF, "deltaTOF/D");
    recTree->Branch("mSquared", &mSquared, "mSquared/D");
    recTree->Branch("pairRapidity", &pairRapidity, "pairRapidity/D"); 
    recTree->Branch("nSigTrk1Pion", &nSigmaTPC[Pion][0], "nSigTrk1Pion/D");
    recTree->Branch("nSigTrk2Pion", &nSigmaTPC[Pion][1], "nSigTrk2Pion/D");
    recTree->Branch("nSigTrk3Pion", &nSigmaTPC[Pion][2], "nSigTrk3Pion/D");
    recTree->Branch("nSigTrk4Pion", &nSigmaTPC[Pion][3], "nSigTrk4Pion/D");
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->Branch("invMass" + particleLables[iPart], &invMass[iPart], "invMass"  + particleLables[iPart] + "/D");
        recTree->Branch("chiPair" + particleLables[iPart], &chiPair[iPart], "chiPair"  + particleLables[iPart] + "/D");
        recTree->Branch("deltaTOFExpected" + particleLables[iPart], &deltaTOFExpected[iPart], "deltaTOFExpected"  + particleLables[iPart] + "/D");
        recTree->Branch("deltaDeltaTOF" + particleLables[iPart], &deltaDeltaTOF[iPart], "deltaDeltaTOF"  + particleLables[iPart] + "/D");  
    }


// Vertex info
    recTree->Branch("vertexZ", &vertexesZ[0], "vertexZ/D");

// Central track info
    for (int i = 0; i < 4; ++i)
    {
        recTree->Branch(Form("dEdx%i",i), &dEdx[i], Form("dEdx%i/D",i));
        recTree->Branch(Form("momentum%i",i), &momentum[i], Form("momentum%i/D",i));
        recTree->Branch(Form("transMomentum%i",i), &transMomentum[i], Form("transMomentum%i/D",i));
        recTree->Branch(Form("charge%i",i), &charge[i], Form("charge%i/D",i));
        recTree->Branch(Form("TOFtime%i",i), &TOFtime[i], Form("TOFtime%i/D",i));
        recTree->Branch(Form("TOFlength%i",i), &TOFlength[i], Form("TOFlength%i/D",i));
        recTree->Branch(Form("DcaXY%i",i), &DcaXY[i], Form("DcaXY%i/D",i));
        recTree->Branch(Form("DcaZ%i",i), &DcaZ[i], Form("DcaZ%i/D",i));
        recTree->Branch(Form("NhitsFit%i",i), &NhitsFit[i], Form("NhitsFit%i/D",i));
        recTree->Branch(Form("NhitsDEdx%i",i), &NhitsDEdx[i], Form("NhitsDEdx%i/D",i));
        recTree->Branch(Form("Eta%i",i), &Eta[i], Form("Eta%i/D",i));
        recTree->Branch(Form("Phi%i",i), &Phi[i], Form("Phi%i/D",i));
        recTree->Branch(Form("Chi2%i",i), &Chi2[i], Form("Chi2%i/D",i));
    }


// event info
    recTree->Branch("fourPiState", &fourPiState, "fourPiState/O");
    recTree->Branch("runNumber", &runNumber, "runNumber/i");
    recTree->Branch("VPDTimeDiff", &VPDTimeDiff, "VPDTimeDiff/D");
    recTree->Branch("VPDSumWest", &VPDSumWest, "VPDSumWest/D");
    recTree->Branch("VPDSumEast", &VPDSumEast, "VPDSumEast/D");
    for (int i = 0; i < nSides; ++i)
    {
        recTree->Branch("BBCSmall" + sideLabel[i], &BBCSmall[i], "BBCSmall" + sideLabel[i] + "/i");
        recTree->Branch("BBCLarge" + sideLabel[i], &BBCLarge[i], "BBCLarge" + sideLabel[i] + "/i"); 
    }

// Setting background Tree
    bcgTree = recTree->CloneTree(0);
    bcgTree->SetName("Background");

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

void FillPlots(int state, int configuration, int isFourPiState){

    for(int trackId = 0; trackId < nTofTrks; ++trackId)
    {
        hdEdxVsqP[state][isFourPiState]->Fill(momentum[trackId]*charge[trackId],dEdx[trackId]);
        hdEdx[state][isFourPiState]->Fill(log10(momentum[trackId]),log10(dEdx[trackId]));
    }

}

void Clear(){
    totalCharge = 0;
    nTofTrks = 0;
    runNumber = 0;
    pairRapidity = -999;


    for (int i = 0; i < 4; ++i)
    {
        dEdx[i] = -999;
        momentum[i] = -999;
        transMomentum[i] = -999;
        TOFtime[i] = -999;
        TOFlength[i] = -999;
        charge[i] = -999;
        nSigmaTPC[Pion][i] = -999;
        nSigmaTPC[Kaon][i] = -999;
        nSigmaTPC[Proton][i] = -999;
        vertexId[i] = -999;
        vertexesZ[i] = -999;
        DcaXY[i] = -999;
        DcaZ[i] = -999;
        NhitsFit[i] = -999;
        NhitsDEdx[i] = -999;
        Eta[i] = -999;
        Phi[i] = -999;
        Chi2[i] = -999;    
    }
}



void CalculatePID(){

    if(TOFtime[0] < 0 || TOFtime[1] < 0 || TOFlength[0] < 0 || TOFlength[1] < 0)
    {
        mSquared = -999.0;
        deltaTOF = -999;   
        for (int iPart = 0; iPart < nParticles; ++iPart)
        {
            deltaTOFExpected[iPart] = -999;
            deltaDeltaTOF[iPart] = -999;
        }

        return;
    }
    deltaTOF = TOFtime[1] - TOFtime[0];

    double speedOfLight2 = speedOfLight*speedOfLight;
    double speedOfLight4 = speedOfLight2*speedOfLight2;
    double length1Squared = TOFlength[0]*TOFlength[0]/(100*100); // convert TOFlength from cm to m
    double length2Squared = TOFlength[1]*TOFlength[1]/(100*100); // convert TOFlength from cm to m
    double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
    double deltaTime4 = deltaTime2*deltaTime2;
    double oneOverMomentum1sq = 1/(momentum[0]*momentum[0]);
    double oneOverMomentum2sq = 1/(momentum[1]*momentum[1]);
    double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
    double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
    double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
    mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
    
//////////////////////
      
/////////////// Calculate deltaDeltaTOFPion
    double expectedTime1, expectedTime2;

    for (int i = 0; i < nParticles; ++i)
    {
        expectedTime1 = (TOFlength[0] / speedOfLight) * sqrt(1 + pow(particleMass[i]/ momentum[0], 2)) * pow(10,7); // in ns
        expectedTime2 = (TOFlength[1] / speedOfLight) * sqrt(1 + pow(particleMass[i]/ momentum[1],2)) * pow(10,7); // in ns             
        deltaTOFExpected[i] = expectedTime2 - expectedTime1;        
        deltaDeltaTOF[i] = deltaTOF - deltaTOFExpected[i];
    }

}
        

