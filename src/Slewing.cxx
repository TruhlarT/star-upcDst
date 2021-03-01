// Run by: ./Analysis input index
// e.g. ./Analysis /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list -1
// or you can open just root file
// ./Analysis /gpfs01/star/pwg/jaroslav/test/star-upcDst/trees/RP_central_elastic_17_v0/merge_files/StUPCRP_central_elastic_17_v0_0000.root


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
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"

using namespace std;

enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum ARM_ID { EU_WU, ED_WD, nArms };
enum STATION_ID { E1, E2, W1, W2, nStations };
enum RP_CONFIGURATION {EUD, EDU, IUU, IDD, nConfiguration}; 

const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP

TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};
TString branchNames[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD")};


TH1F* hInvalidRPtrack; 
TH1I* hNprimVrtMatchTOF;

TFile *infile;
TFile *outfile;
TChain *upcChain;
TChain *rpChain;
StUPCEvent *upcEvt;
StRPEvent *rpEvt;

TH1D* hPtCharged[2];
TH1D* hRatio;

TTree *upcTree;
TTree *rpTree;

TTree *recTree[nRomanPots];

Double_t ADC[2][2]; 
Double_t TAC[2][2];
Double_t rpX[2], rpY[2], rpZ[2]; 

int configurationTable[4][4] = { {0, 2, 5, 7}, {1, 3, 4, 6}, {0, 2, 4, 6}, {1, 3, 5, 7}};


void Init();
void Make();
TFile *CreateOutputTree(const string& out);
bool ConnectInput(const string& in, int fileId);

//_____________________________________________________________________________
int main(int argc, char** argv) {
    //open output file
    outfile = CreateOutputTree("AnalysisOutput.root"); 
    if(!outfile) 
    {
        cout << "Can not open output file." << endl; 
        return -1;
    }

    Init();
    int fileId = 0;
    const string& input = argv[1];
    if( input.find(".root") == string::npos )
        fileId = atoi(argv[2]);

    if(!ConnectInput(argv[1], fileId))
    {
        cout << "No input." << endl; 
        return 1;
    }
    //ask for number of events
    Long64_t nev = upcTree->GetEntries();
    cout<<"Proccesing "<<nev<<" events"<<endl;
    //event loop
    //nev = 1000;
    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        upcTree->GetEntry(iev); 
        Make();
    } 

    for (int i = 1; i < 100; ++i)
    {
        if(hPtCharged[0]->GetBinContent(i) != 0)    
            hRatio->SetBinContent(i,hPtCharged[1]->GetBinContent(i)/hPtCharged[0]->GetBinContent(i));
        else
            hRatio->SetBinContent(i,0);
    }
    
    //close the outputs
    outfile->Write();
    outfile->Close();
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

void Init(){
    hInvalidRPtrack = new TH1F("invalidRPtrack","invalidRPtrack", 6,0,4);
    hNprimVrtMatchTOF = new TH1I("nPrimVrtxMatchTOF", "Primary vertices with at least 1 matched TOF tracks", 10, 0 ,10);

    hPtCharged[0] = new TH1D("ptPositove", "Pt of negative primary tracks", 100, 0 ,10);
    hPtCharged[1] = new TH1D("ptNegative", "Pt of positive primary tracks", 100, 0 ,10);

    hRatio = new TH1D("ratio", "Ratio of positive and negative primary tracks", 100, 0 ,10);
}

void Make(){

    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var)
    {
        if(upcEvt->isTrigger(triggerID[var]))
        {
            //if(var==3 || var==7 || var==9 || var==12 || var==14 || var==15)
            if(var == 8 || var == 11 || var == 16) // Elastic triggers
                // RP_CPT2, RP_CPT2noBBCL, RP_CPT2, RP_CPT2, RP_CPT2noBBCL, RP_CPTnoBBCL
                CPTtrigger=true;
        }
    }

    if(!CPTtrigger)
        return;


    for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
    {
        const StUPCTrack* trk = upcEvt->getTrack(j);
        if( !trk->getFlag(StUPCTrack::kPrimary)) 
            continue;
        
        if(trk->getCharge() > 0)
            hPtCharged[1]->Fill(trk->getPt());
        else
            hPtCharged[0]->Fill(trk->getPt());
    }


    // Vector below will be filled with indexes of good-quality tracks
    vector<int> rpTrackIdVec_perBranch[nBranches];
    vector<int> rpTrackIdVec_perSide[nSides];

    // Loop over all tracks reconstructed in Roman Pots  
    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
    {
        // Get pointer to k-th track in Roman Pot data collection
        StUPCRpsTrack *trk = rpEvt->getTrack(k);
        trk->setEvent(rpEvt);
        // Get ID of a branch in which this k-th track was reconstructed
        int j = trk->branch();
        int side = j<2 ? E : W;

        if( trk->type()==StUPCRpsTrack::rpsGlobal && trk->planesUsed()==8)
        {
        	rpTrackIdVec_perBranch[j].push_back( k );
        	rpTrackIdVec_perSide[side].push_back( k );
        } 

    }

    vector<int> vertexID;
    vector<int> nTOFmTracks;
    int vrtxId;
    bool newVertex;


// Loop over arms - check if have good-quality tracks, selecting branch combination
    for(int i=0; i<2; ++i)
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
        //If exactly one good-quality track in each branch in the arm and there is no track in the other RP do some staff
        if( rpTrackIdVec_perBranch[ branch[E] ].size()==1 
        && rpTrackIdVec_perBranch[ branch[W] ].size()==1
        && rpTrackIdVec_perSide[E].size()==1 && rpTrackIdVec_perSide[W].size()==1)
        {

            vertexID.clear();
            nTOFmTracks.clear();
            for(int j = 0; j < upcEvt->getNumberOfTracks(); ++j)
            {
            // get TPC track object
                const StUPCTrack* trk = upcEvt->getTrack(j);

                if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof) ) 
                    continue;

                vrtxId = trk->getVertexId();
                newVertex = true;
                for(std::vector<int>::size_type i = 0; i != vertexID.size(); ++i) 
                {
                    if( vrtxId == vertexID[i])
                    {
                        newVertex = false;
                        nTOFmTracks[i]++;
                        break;
                    }
                }

                if(newVertex)
                {
                    vertexID.push_back(vrtxId);
                    nTOFmTracks.push_back(1);
                }

            }

            hNprimVrtMatchTOF->Fill(vertexID.size());

            for (int iSide = 0; iSide < 2; ++iSide)
            {
                for(int rpInBranch = 0; rpInBranch < 2; ++rpInBranch)
                {
                    int rpIndex = configurationTable[i][iSide*2 + rpInBranch];
                    ADC[rpInBranch][0] = rpEvt->adc(rpIndex,0);
                    ADC[rpInBranch][1] = rpEvt->adc(rpIndex,1);
                    TAC[rpInBranch][0] = rpEvt->tac(rpIndex,0);
                    TAC[rpInBranch][1] = rpEvt->tac(rpIndex,1);

                    StUPCRpsTrack* trackRP = rpEvt->getTrack(rpTrackIdVec_perSide[iSide][0]); 
                    if(trackRP)
                    {
                        StUPCRpsTrackPoint* trackPoint = trackRP->getTrackPoint(rpInBranch);
                        if(trackPoint)
                        {
                            rpX[rpInBranch] = trackPoint->x(); 
                            rpY[rpInBranch] = trackPoint->y(); 
                            rpZ[rpInBranch] = trackPoint->z();
                        }else
                        {
                            hInvalidRPtrack->Fill(3);
                            continue;
                        }

                    }else
                    {
                        hInvalidRPtrack->Fill(1);
                        continue;
                    }

                }
                recTree[branch[iSide]]->Fill();
            } 
        }
    } // end of loop over arms
}


//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    for (int i = 0; i < nBranches; ++i)
    {
        recTree[i] = new TTree(branchNames[i], branchNames[i]);
        for (int j = 0; j < 2; ++j)
        {
            recTree[i]->Branch(Form("ADC_V_%i",j+1), &ADC[j][0], Form("ADC_V_%i/D",j+1));
            recTree[i]->Branch(Form("ADC_H_%i",j+1), &ADC[j][1], Form("ADC_H_%i/D",j+1));
            recTree[i]->Branch(Form("TAC_V_%i",j+1), &TAC[j][0], Form("TAC_V_%i/D",j+1));
            recTree[i]->Branch(Form("TAC_H_%i",j+1), &TAC[j][1], Form("TAC_H_%i/D",j+1));
            recTree[i]->Branch(Form("rpX_%i",j+1), &rpX[j], Form("rpX_%i/D",j+1));
            recTree[i]->Branch(Form("rpY_%i",j+1), &rpY[j], Form("rpY_%i/D",j+1));
            recTree[i]->Branch(Form("rpZ_%i",j+1), &rpZ[j], Form("rpZ_%i/D",j+1));
        }
    }


 	return outputFile;

}//CreateOutputTree

bool ConnectInput(const string& in, int fileId) {
    //input from file or chain
    upcTree = 0x0;
    if( in.find(".root") != string::npos ) 
    {
        cout << "Input from root file" << endl;
        infile = TFile::Open(in.c_str(), "read");
        if(!infile) 
            return false;
        upcTree = dynamic_cast<TTree*>( infile->Get("mUPCTree") );
    } else 
    {
        cout << "Input from chain" << endl;
        upcChain = new TChain("mUPCTree");
        ifstream instr(in.c_str());
        string line;
        int lineId=0;
        while(getline(instr, line)) 
        {
            if(fileId==lineId || fileId== -1)
                upcChain->AddFile(line.c_str());
            lineId++;
        }
        instr.close();
        upcTree = dynamic_cast<TTree*>( upcChain );
    }

    if(!upcTree) 
        return false;

    upcTree->SetBranchAddress("mUPCEvent", &upcEvt);
    upcTree->SetBranchAddress("mRPEvent", &rpEvt); 

    return true;

}//ConnectInput