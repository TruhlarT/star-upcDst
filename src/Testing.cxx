// Run by: ./Testing file.list
// e.g. ./Testing /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list
// or you can open just root file
// ./Testing /star/data01/pwg_tasks/upc02/Part9/18143045/18143045.root
// or you can open n-th root file in file.list
// ./Testing file.list index 


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

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"

#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCTofHit.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"

using namespace std;

// Configuration setting
const bool IsRP = false;
const bool IsTrigger = false;
const bool DEBUG = false;


enum { kAllEvents = 1, kTrigger,  kRP, kRPinFiducial, kTPCTOF, kSameVertex, kZeroCharge,
       kZVertex, kHitsFit, kHitsDEdx, kDCAz, kDCAxy, kEta, kMandelT, kExclusive, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum XY_COORDINATE {X, Y, nCoordinates};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum ARM_ID { EU_WU, ED_WD, nArms };
enum STATION_ID { E1, E2, W1, W2, nStations };
enum RP_CONFIGURATION {EUD = 0, EDU = 1, IUU = 2, IDD = 3, nConfiguration}; 

const double particleMass[nParticles] = { 0.13957, 0.49367, 0.93827}; // GeV /c^2
const double speedOfLight = 299792458; // m/s
const double beamMomentum = 254.9; // GeV
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP

Bool_t trigger[nTriggers];


TString branchName[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD") };
TString summaryLabels[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("2 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString summaryLabels2[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("4 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString systemID[3] = { TString("Combi"), TString("InEl"), TString("El")};
TString systemState[4] = { TString("TPC2t"), TString("TOF2trk"), TString("Q0"), TString("Excl")};
TString configLabels[4] = { TString("EUD"), TString("EDU"), TString("IUU"), TString("IDD")};
TString systemLabels[2] = { TString("2part"), TString("4part")};
TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString sideLabel[nSides] = { TString("East"), TString("West")};
TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};


TH1I* hAnaFlow; // control plot
TH1F* hTriggerBits; // control plot
TH1F* hConfiguration;
TH1F* hNumberTrackPerBranch[nBranches]; // number of track per each brunch (west up, west down, east up, east down)
TH1F* hNumberTrack;
TH1F* hTrackPlanesUsed;
TH1F* hTrackPointPlanesUsed;
TH1F* hNumberTOFmatchTrack;
TH1F* hInvalidRPtrack; 

static const int nVar = 20;
static const int startTrackIndex = 1;
TString varName[nVar] = {TString("TOFmultiplicity"), TString("VertexZ"), TString("VertexId"), TString("NSigmaTPCPion"), 
                    TString("NSigmaTPCKaon"), TString("NSigmaTPCProton"), TString("dEdx"), TString("Momentum"), 
                    TString("Pt"), TString("Charge"), TString("TofTime"), TString("TofPathLength"), 
                    TString("Nhits"), TString("NhitsFit"), TString("NhitsDEdx"), TString("DCAZ"), 
                    TString("DCAXY"), TString("Eta"), TString("Phi"), TString("Chi2") };
Int_t nBin[] = { 50, 200, 101, 201, 201, 201, 100, 100, 100, 11, 400, 200, 61, 61, 61, 100, 100, 100, 100, 100};  // Binning for histograms
Double_t minBin[] = { -0.5, -200, -0.5, -100.5, -100, -100.5, 0, 0, 0, -5.5, -2000, 100, -0.5, -0.5, -0.5, -1.5, 0, -2, -3.5, 0};
Double_t maxBin[] = { 49.5, 200, 100, 100.5, 100.5, 100.5, 100, 20, 20, 5.5, 60000, 500, 60.5, 60.5, 60.5, 1.5, 3.5, 3.5, 3.5, 10};
TH1F* hIneterestingVar[nVar][3]; // 3 states: TPC primary, TPC matched with TOF, Golden  
Double_t ineterestingVar[nVar];
vector<double> trackVar[nVar];
Double_t varToTree[nVar][4];

TH1F* hEtaTest[2];
TH1F* hInvMass[3][nParticles];
TH2D* hEtaPhi[3];

TH1F* hMissingPt[12];
TH2D* hdEdxVsqP[12];
TH2D* hdEdx[12];

TH1F* hTofHitTime;

TFile *infile;
TFile *outfile;
TChain *upcChain;
TChain *rpChain;
StUPCEvent *upcEvt;
StRPEvent *rpEvt;

TTree *upcTree;
TTree *rpTree;

TTree *recTree, *bcgTree;
Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;
Double_t VPDSumEast, VPDSumWest, VPDTimeDiff;


Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t missingPt, deltaTOF, mSquared, pairRapidity;

/////////////////////////////////

Bool_t elastic, fourPiState;

static const int nVarRP = 15;
TString varNameRP[nVarRP] = {TString("BBCSmall"), TString("BBCLarge"), TString("PxRP"), TString("PyRP"), 
                    TString("ThetaRP"), TString("PhiRP"), TString("TimeRP"), TString("momentumRP"), 
                    TString("PtRP"), TString("EtaRP"), TString("XRP"), TString("YRP"), 
                    TString("ZRP"), TString("TransMomSqRP"), TString("XiRP")};
Int_t nBinRP[] = { 400, 400, 100, 100, 100, 100, 100, 200, 100, 100, 100, 100, 100, 100, 100};  // Binning for histograms
Double_t minBinRP[] = { 0, 0, -1, -1, 0, -4, -0.01, 100, 0, -10, -0.02, -0.08, -16, -2.0, -0.5};
Double_t maxBinRP[] = { 4000, 4000, 1, 1, 0.03, 4, 0.01, 400, 2, 10, 0.06, 0.08, 16, 0, 1};
TH1F* hRPVar[nVarRP][3]; // 3 states: TPC primary, TPC matched with TOF, Golden  
Double_t rpVariable[nVarRP][nSides];
 
Int_t RPConfiguration;
Double_t ADC[nRomanPots][2]; // RP_ID   0,    1,    2,   3,   4,   5,   6, 7
Double_t TAC[nRomanPots][2]; // RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D


void Init();
void Make();
void FillPlots(int state, int configuration);
TFile *CreateOutputTree(const string& out);
bool ConnectInput(int argc, char** argv);
void Clear();
void CalculatePID();
int FindRPConfig(vector<int> rpTrackIdVec_perBranch[nBranches], vector<int> rpTrackIdVec_perSide[nSides]);
void FillTracksQualityPlots(int state);
void SaveVariable();

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

    //close the outputs
    outfile->Write();
    outfile->Close();
    cout<<"Ending Testing... GOOD BYE!"<<endl;
    return 0;
}//main

void Init(){

    if(DEBUG)
        cout<<"Starting Init()"<<endl;

    hAnaFlow = new TH1I("AnaFlow", "CutsFlow", kMax-1, 1, kMax);

    hEtaPhi[0] = new TH2D("EtaPhiAll","EtaPhiAll",200,-2,2,100,-3.5,3.5);
    hEtaPhi[1] = new TH2D("EtaPhiMatched","EtaPhiMatched",200,-2,2,100,-3.5,3.5);
    hEtaPhi[2] = new TH2D("EtaPhiNotMtch","EtaPhiNotMtch",200,-2,2,100,-3.5,3.5);

    for (int i = 0; i < nParticles; ++i)
    {
        hInvMass[0][i] = new TH1F("invMassNonExc"+particleLables[i], "invMassNonExc El + Inel "+particleLables[i], 200, 0.0, 5.0);
        hInvMass[1][i] = new TH1F("invMassNonExcEl"+particleLables[i], "invMassNonExc El "+particleLables[i], 200, 0.0, 5.0);
        hInvMass[2][i] = new TH1F("invMassNonExcInel"+particleLables[i], "invMassNonExc Inel "+particleLables[i], 200, 0.0, 5.0);
    }


    hInvalidRPtrack = new TH1F("invalidRPtrack","invalidRPtrack", 6,0,4);

    if(IsTrigger)
    {
        hTriggerBits = new TH1F("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
        for(int tb=0; tb<nTriggers; ++tb)
        {
            TString label; label.Form("%d",triggerID[tb]);
            hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
        }
    }

    hConfiguration = new TH1F("Configuration", "Track Configuration", 4, -0.5, 3.5);
    for(int tb=0; tb<4; ++tb) 
        hConfiguration->GetXaxis()->SetBinLabel(tb+1, configLabels[tb]);

    for(int i=0; i<nBranches; ++i)
        hNumberTrackPerBranch[i] = new TH1F("NumberTracksPerBranch,"+branchName[i],"Number of tracks in branch "+branchName[i], 8, -0.5, 7.5);

    hNumberTrack = new TH1F("NumberTracks", "Number of Tracks in RP", 40, -0.5, 39.5);
    hNumberTOFmatchTrack = new TH1F("hNumberTOFmatchTrack", "Number of TOF matched tracks", 40, -0.5, 39.5);
    hTrackPlanesUsed = new TH1F("hTrackPlanesUsed", "Number of planes used for track reconstruction", 40, -0.5, 39.5);
    hTrackPointPlanesUsed = new TH1F("hTrackPointPlanesUsed", "Number of planes used for track point reconstruction", 40, -0.5, 39.5); 


    hTofHitTime = new TH1F("TofHitTime", "TOF hit time", 400, -2000, 60000);

    for (int i = 0; i < 3; ++i)
    {
        TString name = "";
        if(i==0){
            outfile->mkdir("TPCprimary")->cd();
        }else if(i==1){
            outfile->cd();
            outfile->mkdir("TOFmatched")->cd();
            name = "_TOFmatched";
        }else if(i==2){
            outfile->cd();
            outfile->mkdir("Golden")->cd();
            name = "_Golden";
        }

        for (int j = 0; j < nVar; ++j)
            hIneterestingVar[j][i] = new TH1F("h" + varName[j] + name,"h" + varName[j] + name, nBin[j], minBin[j], maxBin[j]);
    }

    if(IsRP)
    {
        TDirectory* RPDir = outfile->mkdir("RomanPot");
        RPDir->cd();
        for (int i = 0; i < 3; ++i)
        {
            TString name = "";
            if(i==0){
                RPDir->mkdir("TPCprimary")->cd();
            }else if(i==1){
                RPDir->mkdir("TOFmatched")->cd();
                name = "_TOFmatched";
            }else if(i==2){
                RPDir->mkdir("Golden")->cd();
                name = "_Golden";
            }

            for (int j = 0; j < nVarRP; ++j)
                hRPVar[j][i] = new TH1F("h" + varNameRP[j] + name,"h" + varNameRP[j] + name, nBinRP[j], minBinRP[j], maxBinRP[j]);
        }
    }

    outfile->cd();
    for(int i=0; i<12;++i){
        if(i==0){
            outfile->mkdir("All")->cd();
        }else if(i==4){
            outfile->cd();
            outfile->mkdir("Inelastic")->cd();
        }else if(i==8){
            outfile->cd();
            outfile->mkdir("Elastic")->cd();
        }

        hdEdxVsqP[i] = new TH2D("dEdxVsqP_"+systemState[i%4]+"_"+systemID[i/4],"dE/dx Vs #frac{q}{e} P",200,-2,2,100,0,20);
        hdEdx[i] = new TH2D("dEdx_"+systemState[i%4]+"_"+systemID[i/4],"#log_{10} dE/dx [keV/cm] vs #log_{10} p [GeV/c]",200,-1,1,100,0.1,2);
        hMissingPt[i] = new TH1F( "MissingPt_"+systemState[i%4]+"_"+systemID[i/4], "p^{miss}_{T} [GeV/c]", 200, 0, 2 );
        
    }


    outfile->cd();


    if(DEBUG)
        cout<<"Ending Init()"<<endl;

}

void Make(){
    hAnaFlow->Fill(kAllEvents);

    Clear();

// Trigger part starts here...
    if(IsTrigger) 
    {
        bool CPTtrigger = false;
        for(int var = 0; var < nTriggers; ++var)
        {
            if(upcEvt->isTrigger(triggerID[var]))
            {
                hTriggerBits->Fill(var);
                trigger[var] = true;
                if(var==3 || var==7 || var==9 || var==12 || var==14 || var==15)
                    CPTtrigger=true;
                    // RP_CPT2, RP_CPT2noBBCL, RP_CPT2, RP_CPT2, RP_CPT2noBBCL, RP_CPTnoBBCL
                //if(var==7 || var==14 )
                    // RP_CPT2noBBCL, RP_CPT2noBBCL
            }
        }

        if(!CPTtrigger)
            return;
        hAnaFlow->Fill(kTrigger);
    }

// TOF hit calib - testing code

    for(int TOFhit = 0; TOFhit < upcEvt->getNumberOfHits(); ++TOFhit)
    {
    // get TPC track object
        const StUPCTofHit* hit = upcEvt->getHit(TOFhit);
    
    }    


// Trigger part ends here...
// Roman Pots part start here...
    // Vector below will be filled with indices of good-quality tracks
    vector<int> rpTrackIdVec_perBranch[nBranches];
    vector<int> rpTrackIdVec_perSide[nSides];
    if(IsRP) 
    {
        nTracks = (Int_t) rpEvt->getNumberOfTracks();

        int numberOfTracks = 0;
        int numberOfTracksPerBranch[nBranches] = {0, 0, 0, 0};

        // Loop over all tracks reconstructed in Roman Pots  
        for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
        {
            // Get pointer to k-th track in Roman Pot data collection
            StUPCRpsTrack *trk = rpEvt->getTrack(k);
            trk->setEvent(rpEvt);
            // Get ID of a branch in which this k-th track was reconstructed
            int j = trk->branch();
            int side = j < 2 ? E : W;

            ++numberOfTracks;
            ++numberOfTracksPerBranch[j];
            hTrackPlanesUsed->Fill(trk->planesUsed());
            if(trk->getTrackPoint(0))
                hTrackPointPlanesUsed->Fill(trk->getTrackPoint(0)->planesUsed());

            if(trk->getTrackPoint(1))
                hTrackPointPlanesUsed->Fill(trk->getTrackPoint(1)->planesUsed());

            // If track is global (made of track-points in 1st and 2nd station)
            // and all 8 planes were used in reconstruction - store its ID
            if( trk->type()==StUPCRpsTrack::rpsGlobal && trk->planesUsed()==8)
            {
                rpTrackIdVec_perBranch[j].push_back( k );
                rpTrackIdVec_perSide[side].push_back( k );
            } 
            // a bit looser selection
            /* if( (trk->getTrackPoint(0) ? trk->getTrackPoint(0)->planesUsed()>=3 : false) &&
                (trk->getTrackPoint(1) ? trk->getTrackPoint(1)->planesUsed()>=3 : false) && 
                trk->type()==StUPCRpsTrack::rpsGlobal && trk->planesUsed()<=8)
            {
                rpTrackIdVec_perBranch[j].push_back( k );
                rpTrackIdVec_perSide[side].push_back( k );
            } */
        }

        hNumberTrack->Fill(numberOfTracks);
        for(int i=0; i<nBranches; ++i) 
            hNumberTrackPerBranch[i]->Fill(numberOfTracksPerBranch[i]);

        RPConfiguration = FindRPConfig(rpTrackIdVec_perBranch, rpTrackIdVec_perSide);
        if( RPConfiguration == -1)
            return;
        
        elastic = false;
        if(RPConfiguration > 1)
            elastic = true;
            
        hAnaFlow->Fill(kRP);
        rpVariable[0][0] = upcEvt->getBBCSmallEast();
        rpVariable[1][0] = upcEvt->getBBCLargeEast();
        rpVariable[0][1] = upcEvt->getBBCSmallWest();
        rpVariable[1][1] = upcEvt->getBBCLargeWest();
        double x,y;
        bool protonInRange[nSides] = {false, false};
        for (int iSide = 0; iSide < nSides; ++iSide)
        {
            StUPCRpsTrack* trackRP = rpEvt->getTrack(rpTrackIdVec_perSide[iSide][0]);
            if(trackRP)
            {   
                StUPCRpsTrackPoint* trackPoint = trackRP->getTrackPoint(0);

                if(trackPoint)
                {
                    rpVariable[10][iSide] = trackPoint->x(); 
                    rpVariable[11][iSide] = trackPoint->y(); 
                    rpVariable[12][iSide] = trackPoint->z();
                }else
                {
                    hInvalidRPtrack->Fill(3);
                    return;
                }
                rpVariable[2][iSide] = trackRP->pVec().X();
                rpVariable[3][iSide] = trackRP->pVec().Y();
                rpVariable[4][iSide] = trackRP->thetaRp();
                rpVariable[5][iSide] = trackRP->phiRp(); 
                rpVariable[6][iSide] = trackRP->time();
                rpVariable[7][iSide] = trackRP->p();
                rpVariable[8][iSide] = trackRP->pt();
                rpVariable[9][iSide] = trackRP->eta();
                rpVariable[13][iSide] = trackRP->t(beamMomentum);
                rpVariable[14][iSide] = trackRP->xi(beamMomentum);

                x = rpVariable[2][iSide];
                y = rpVariable[3][iSide];
                if( abs(y) < 0.8 && abs(y) > 0.4 && x > -0.27 && (x + 0.6)*(x + 0.6) + y*y < 1.25 )
                    protonInRange[iSide] = true;
            }else
            {
                hInvalidRPtrack->Fill(1);
                return;
            }
        }

        if(!protonInRange[0] || !protonInRange[1])
            return;
        
        hAnaFlow->Fill(kRPinFiducial);

        hConfiguration->Fill(RPConfiguration);
    }
// Roman Pots part ends here...

// Central part starts here....
    TLorentzVector centralTracksTotalFourMomentum[nParticles];



    ineterestingVar[0] = upcEvt->getNumberOfHits();
    // loop over all TPC tracks
    for(int TPCTrack = 0; TPCTrack < upcEvt->getNumberOfTracks(); ++TPCTrack)
    {
    // get TPC track object
        const StUPCTrack* trk = upcEvt->getTrack(TPCTrack);

        TLorentzVector trkLVector[nParticles];
        for (int iPart = 0; iPart < nParticles; ++iPart)
            trk->getLorentzVector(trkLVector[iPart], particleMass[iPart]);

        if( !trk->getFlag(StUPCTrack::kPrimary)) 
            continue;

        hEtaPhi[0]->Fill(trk->getEta(),trk->getPhi());

        ineterestingVar[1] = trk->getVertex()->getPosZ();
        ineterestingVar[2] = trk->getVertexId();
        ineterestingVar[3] = trk->getNSigmasTPCPion();
        ineterestingVar[4] = trk->getNSigmasTPCKaon();
        ineterestingVar[5] = trk->getNSigmasTPCProton();
        ineterestingVar[6] = trk->getDEdxSignal()*1000000;
        ineterestingVar[7] = trkLVector[0].P();
        ineterestingVar[8] = trk->getPt();
        ineterestingVar[9] = trk->getCharge();
        ineterestingVar[10] = trk->getTofTime();
        ineterestingVar[11] = trk->getTofPathLength();
        ineterestingVar[12] = trk->getNhits();
        ineterestingVar[13] = trk->getNhitsFit();
        ineterestingVar[14] = trk->getNhitsDEdx();
        ineterestingVar[15] = trk->getDcaZ();
        ineterestingVar[16] = trk->getDcaXY();
        ineterestingVar[17] = trk->getEta();
        ineterestingVar[18] = trk->getPhi();
        ineterestingVar[19] = trk->getChi2();

        FillTracksQualityPlots(0);

        if(!trk->getFlag(StUPCTrack::kTof) )
        {
            hEtaPhi[2]->Fill(trk->getEta(),trk->getPhi());
            continue;
        } 

        hEtaPhi[1]->Fill(trk->getEta(),trk->getPhi());

        // read basic information about the track   
        for (int iVar = startTrackIndex; iVar < nVar; ++iVar)
            trackVar[iVar].push_back(ineterestingVar[iVar]);

        totalCharge += static_cast<int>( trk->getCharge() );

        for (int iPart = 0; iPart < nParticles; ++iPart)
            centralTracksTotalFourMomentum[iPart] += trkLVector[iPart];
        
        ++nTofTrks;
    }
// Central part ends here...
    if(nTofTrks == 0)
        return;
    if(DEBUG)
        cout<<"Exiting TPC track loop with "<<nTofTrks<<" NTOFtrks.."<<endl;
// Some other stuff...
    FillTracksQualityPlots(1);
    if(DEBUG)
        cout<<"Filled"<<endl;
    hNumberTOFmatchTrack->Fill(nTofTrks);
    if(IsRP)
        missingPt = (centralTracksTotalFourMomentum[0].Vect() + rpEvt->getTrack( rpTrackIdVec_perSide[E][0] )->pVec() + rpEvt->getTrack( rpTrackIdVec_perSide[W][0] )->pVec() ).Pt();
    pairRapidity = centralTracksTotalFourMomentum[0].Rapidity();
    for (int iPart = 0; iPart < nParticles; ++iPart)
        invMass[iPart] = centralTracksTotalFourMomentum[iPart].M();

    FillPlots(0,RPConfiguration);
    if(DEBUG)
        cout<<"Plots 0"<<endl;

    runNumber = upcEvt->getRunNumber();
    VPDTimeDiff = upcEvt->getVPDTimeDiff();
    VPDSumWest = upcEvt->getVPDSumWest();
    VPDSumEast = upcEvt->getVPDSumEast();

    if( nTofTrks == 2)
    {
        if(trackVar[2][0] == trackVar[2][1])
        {
            hAnaFlow->Fill(kSameVertex);
            FillPlots(1,RPConfiguration);
            if(DEBUG)
                cout<<"Plots 1"<<endl;

            for(int iPart = 0; iPart < nParticles; ++iPart)
                chiPair[iPart] = trackVar[3 + iPart][0]*trackVar[3 + iPart][0] + trackVar[3 + iPart][1]*trackVar[3 + iPart][1];

            if(IsRP)
            {
                for (int iRP = 0; iRP < nRomanPots; ++iRP)
                {
                    ADC[iRP][0] = rpEvt->adc(iRP,0);
                    ADC[iRP][1] = rpEvt->adc(iRP,1);
                    TAC[iRP][0] = rpEvt->tac(iRP,0);
                    TAC[iRP][1] = rpEvt->tac(iRP,1);
                }
            }


            CalculatePID();

            if(DEBUG)
                cout<<"PIDed"<<endl;

            if(totalCharge != 0)
            {
                if( missingPt < 0.1 )
                    //bcgTree->Fill();
                return;
            }

            hAnaFlow->Fill(kZeroCharge);
            FillPlots(2,RPConfiguration);                  
            if(DEBUG)
                cout<<"Plots 2"<<endl;


            if(trackVar[1][0] < 80 && trackVar[1][0] > -80)
            {
                hAnaFlow->Fill(kZVertex);
                if(trackVar[13][0] >=25 && trackVar[13][1] >= 25)
                {
                    hAnaFlow->Fill(kHitsFit);
                    if(trackVar[14][0] >= 15 && trackVar[14][1] >= 15)
                    {  
                        hAnaFlow->Fill(kHitsDEdx);
                        if(trackVar[15][0] < 1 && trackVar[15][0] > -1 && trackVar[15][1] < 1 && trackVar[15][1] > -1)
                        {
                            hAnaFlow->Fill(kDCAz);
                            if(trackVar[16][0] < 1.5 && trackVar[16][1] < 1.5)
                            {
                                hAnaFlow->Fill(kDCAxy);
                                if(trackVar[17][0] > -0.7 && trackVar[17][0] < 0.7 && trackVar[17][1] > -0.7 && trackVar[17][1] < 0.7)
                                {
                                    hAnaFlow->Fill(kEta);
                                    if(IsRP)
                                    {
                                        if(rpVariable[13][0] < -0.12 && rpVariable[13][1] < -0.12 && rpVariable[13][0] > -1.0  && rpVariable[13][1] > -1.0)
                                        {
                                            hAnaFlow->Fill(kMandelT);               
                                            if(missingPt <= 0.1)
                                            {
                                                hAnaFlow->Fill(kExclusive);
                                                FillTracksQualityPlots(2);
                                            }
                                        }
                                    }else
                                    {
                                        FillTracksQualityPlots(2);
                                    }
                                    if(DEBUG)
                                        cout<<"Filled Golden"<<endl;
                                }
                            }
                        }
                    }
                }
            }
        } 
    }

    SaveVariable();
    if(DEBUG)
        cout<<"Filling tree"<<endl;
    recTree->Fill();
    if(DEBUG)
        cout<<"Tree filled"<<endl;

}


//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    if(DEBUG)
        cout<<"Creating output file"<<endl;
    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    recTree = new TTree("recTree", "recTree");

// PID and some quality event info
    recTree->Branch("missingPt", &missingPt, "missingPt/D");
    recTree->Branch("mSquared", &mSquared, "mSquared/D");
    recTree->Branch("pairRapidity", &pairRapidity, "pairRapidity/D"); 
    
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->Branch("invMass" + particleLables[iPart], &invMass[iPart], "invMass"  + particleLables[iPart] + "/D");
        recTree->Branch("chiPair" + particleLables[iPart], &chiPair[iPart], "chiPair"  + particleLables[iPart] + "/D");
    }




// Central track info
    for (int i = 0; i < 4; ++i)
        for (int iVar = startTrackIndex; iVar < nVar; ++iVar)
            recTree->Branch(varName[iVar] + Form("%i",i), &varToTree[iVar][i], varName[iVar] + Form("%i/D",i));


// RP track info  
    if(IsRP)
    {
        for (int i = 0; i < nSides; ++i)
            for (int iVar = 0; iVar < nVarRP; ++iVar)
                recTree->Branch(varNameRP[iVar] + sideLabel[i], &rpVariable[iVar], varNameRP[iVar] + sideLabel[i] + "/D");

        // RP event info
        recTree->Branch("elastic", &elastic, "elastic/O");
        for (int i = 0; i < nRomanPots; ++i)
        {
            recTree->Branch("ADC_" + rpNames[i] + "V", &ADC[i][0], "ADC_" + rpNames[i] + "V/D");
            recTree->Branch("ADC_" + rpNames[i] + "H", &ADC[i][1], "ADC_" + rpNames[i] + "H/D");
            recTree->Branch("TAC_" + rpNames[i] + "V", &TAC[i][0], "TAC_" + rpNames[i] + "V/D");
            recTree->Branch("TAC_" + rpNames[i] + "H", &TAC[i][1], "TAC_" + rpNames[i] + "H/D");
        }
    }


// event info
    recTree->Branch("runNumber", &runNumber, "runNumber/i");
    recTree->Branch("VPDTimeDiff", &VPDTimeDiff, "VPDTimeDiff/D");
    recTree->Branch("VPDSumWest", &VPDSumWest, "VPDSumWest/D");
    recTree->Branch("VPDSumEast", &VPDSumEast, "VPDSumEast/D");
    recTree->Branch("RP_CPT2_570701", &trigger[3], "RP_CPT2_570701/O");
    recTree->Branch("RP_CPT2noBBCL_570705", &trigger[7], "RP_CPT2noBBCL_570705/O");
    recTree->Branch("RP_CPT2_570711", &trigger[9], "RP_CPT2_570711/O");
    recTree->Branch("RP_CPT2_590701", &trigger[12], "RP_CPT2_590701/O");
    recTree->Branch("RP_CPT2noBBCL_590705", &trigger[14], "RP_CPT2noBBCL_590705/O");
    recTree->Branch("RP_CPTnoBBCL_590708", &trigger[15], "RP_CPTnoBBCL_590708/O");

// Setting background Tree
    bcgTree = recTree->CloneTree(0);
    bcgTree->SetName("Background");

    if(DEBUG)
        cout<<"Output file created"<<endl;

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
    if(IsRP)
        upcTree->SetBranchAddress("mRPEvent", &rpEvt); 

    return true;
}//ConnectInput

void FillPlots(int state, int configuration){
    for(int tmp = 0; tmp < 3; ++tmp)
    {
        if(tmp == 1 && (configuration!=IUU && configuration!=IDD))
            continue;

        if(tmp == 2 && (configuration!=EUD && configuration!=EDU))
            continue;

        for(int trackId = 0; trackId < nTofTrks; ++trackId)
        {
            hdEdxVsqP[4*tmp+state]->Fill(trackVar[7][trackId]*trackVar[9][trackId],trackVar[6][trackId]);
            hdEdx[4*tmp+state]->Fill(log10(trackVar[7][trackId]),log10(trackVar[6][trackId]));
        }

        hMissingPt[4*tmp+state]->Fill(missingPt); 
    }

}

void Clear(){
    if(DEBUG)
        cout<<"Clearing"<<endl;
    totalCharge = 0;
    nTofTrks = 0;

    for (int i = startTrackIndex; i < nVar; ++i)
        trackVar[i].clear();

    for (int i = 0; i < nTriggers; ++i)
        trigger[i] = false;
    
    if(DEBUG)
        cout<<"Cleared"<<endl;
}



void CalculatePID(){

    if(trackVar[10][0] < 0 || trackVar[10][1] < 0 || trackVar[11][0] < 0 || trackVar[11][1] < 0)
    {
        mSquared = -999.0;
        return;
    }
    
    double speedOfLight2 = speedOfLight*speedOfLight;
    double speedOfLight4 = speedOfLight2*speedOfLight2;
    double length1Squared = trackVar[11][0]*trackVar[11][0]/(100*100); // convert TOFlength from cm to m
    double length2Squared = trackVar[11][1]*trackVar[11][1]/(100*100); // convert TOFlength from cm to m
    double deltaTOF = trackVar[10][1] - trackVar[10][0];
    double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
    double deltaTime4 = deltaTime2*deltaTime2;
    double oneOverMomentum1sq = 1/(trackVar[7][0]*trackVar[7][0]);
    double oneOverMomentum2sq = 1/(trackVar[7][1]*trackVar[7][1]);
    double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
    double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
    double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
    mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);      
}
        

int FindRPConfig(vector<int> rpTrackIdVec_perBranch[nBranches], vector<int> rpTrackIdVec_perSide[nSides]){
    int config = -1;

    //If exactly one good-quality track in each branch in the arm and there is no track in the other RP do some staff
    if(rpTrackIdVec_perSide[E].size()!=1 || rpTrackIdVec_perSide[W].size()!=1)
        return config;

    for(int i=0; i<nConfiguration; ++i)
    { 
    // Define IDs of branches based on known ID of the arm
        int branch[nSides];
        switch(i) {
            case EUD :
            {
                branch[E] = EU;
                branch[W] = WD;
            }
            break;
            case EDU :
            {
                branch[E] = ED;
                branch[W] = WU;
            }
            break; 
            case IUU :
            {
                branch[E] = EU;
                branch[W] = WU;
            }
            break; 
            case IDD :
            {
                branch[E] = ED;
                branch[W] = WD;
            }
            break;        
        }

        if( rpTrackIdVec_perBranch[ branch[E] ].size()==1 
        && rpTrackIdVec_perBranch[ branch[W] ].size()==1)
            return i;
    }

    return config;
}

void FillTracksQualityPlots(int state)
{
    if(state == 0)
    {
        for (int i = 0; i < nVar; ++i)
            hIneterestingVar[i][state]->Fill(ineterestingVar[i]);

        if(IsRP)
            for (int i = 0; i < nVarRP; ++i)
                for (int j = 0; j < nSides; ++j)
                    hRPVar[i][state]->Fill(rpVariable[i][j]);
    }else
    {
        if(IsRP)
            for (int i = 0; i < nVarRP; ++i)
                for (int j = 0; j < nSides; ++j)
                    hRPVar[i][state]->Fill(rpVariable[i][j]);
        
        for (int i = 0; i < startTrackIndex; ++i)
            hIneterestingVar[i][state]->Fill(ineterestingVar[i]);

        for (int iTrck = 0; iTrck < nTofTrks; ++iTrck)
            for (int i = startTrackIndex; i < nVar; ++i)
                hIneterestingVar[i][state]->Fill(trackVar[i][iTrck]);
    }
}

void SaveVariable()
{
    for (int ij = 0; ij < 4; ++ij)
        for (int i = 0; i < nVar; ++i)
            trackVar[i].push_back(-999);

    for (int ij = 0; ij < 4; ++ij)
        for (int i = 0; i < nVar; ++i)
            varToTree[i][ij] = trackVar[i][ij];

}
