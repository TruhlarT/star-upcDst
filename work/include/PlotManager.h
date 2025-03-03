#ifndef PlotManager_h
#define PlotManager_h

// include headers
#include "Libreries.h"
#include "RecTree.h"
#include "Util.h"
#include "RunDef.h"


class PlotManager{
    public:
        PlotManager();
        ~PlotManager();

        bool Init(int argc, char** argv);
        void runVertexStudy();
        void runMainAnaPlots();
        void initEmbeddingStudy( unsigned int particle);
        void runEmbeddingStudy();
        void runEmbeddingQA();
        void runTofQA();
        void runDsmEffStudy(int study);
        void runElasticStudy();
        void runTofTrigStudy();
        void runRpTrigStudy();
        void runProbOfRetainEvent();
        void runRPMCPlots();

    private:
        TFile *embFile[nParticles]; 
        TFile *outFile, *inFile, *effFile, *graniittiFile;
        TTree *mTree[nStudies];
        RecTree *mRecTree[nStudies];
        RecTree *mCurrentTree;
        Util *mUtil;

        TDirectory *mStudyDir[nStudies];
        TDirectory *mDir, *mCurrDir;

        TCanvas *canvas;
        TLegend *legend;
        TPaveText *text;
        TLine *line;

        // holder for runs
        std::vector<double> mRuns;

        // Luminosity info
        std::map< int, double> mInstLumiPerRun; // in [#mub^{-1}s^{-1}]
        std::map< int, double> mIntegLumiPerRun; // in [#pb^{-1}]
        std::map< int, int> mTimeOfRun; // in [s]

        // efficiency info
        std::map< int, double> mVertexEffPerFill; 
        std::map< int, double> mVertexEffErrorPerFill; 

        // TPC and TOF eff
        double tpcEff[nTPCnHitsStudies+1][nParticles*nSigns]; // nominal, nDedxLoose, nFitLoose, nDedxTight, nFitTight
        double tofEff[nParticles*nSigns]; 

        // PID efficiency
        TH2F* hPIDEff[nParticles][nPidVariation+1];
        // size of non-exclusive background
        double mBcgFraction[nParticles][nRpConfigurations+1];
        double mBcgFracError[nParticles][nRpConfigurations+1];
        double mPtMissEff[nParticles][nRpConfigurations+1];

        double txtSize;
        double tSize;
        inline void setTextSize(double s1 = textSize, double s2 = titleSize ){ txtSize = s1; tSize = s2;};

        //init functions
        bool ConnectInputFile( TString inputFile, TFile **file); 
        TFile *CreateOutputFile(const string& out);
        bool readLumiFile();
        bool loadEfficiencies();
        
        // useful functions
        void changeDir(int study, TString dirName);
        void changeSubDir(TString dirName);
        void AddFunctionToHistogram(TH1* histogram, TF1* function);
        void replaceNegativeBinsWithAbsolute(TH1* histogram);
        void AddHist(TString name, TString xlabel, unsigned int nBins, double low, double max, unsigned int flag, 
   vector<TString> &histList, vector<TString> &labelList, vector<unsigned int> &nBinsVec, vector<double> &binLow, 
   vector<double> &binMax, vector<unsigned int> &flags);

///////// MainAnaPlots /////////
        TGraphAsymmErrors *systUncertainty;
        //histograms
        static constexpr unsigned int nTotalTPCTOFSysStudies = nTPCAppSysStudies + nTOFAppSysStudies + nTPCnHitsStudies + nPidVariation + nDcaVariation;
        TH2F *hMissingPtVsInvMass[nParticles][nRpConfigurations+1][nTotalTPCTOFSysStudies+1];
        TH1F *hInvMassCorr[nParticles][nRpConfigurations+1][nTotalTPCTOFSysStudies+1];
        TH1F *hBcgSysStudy[nParticles][nRpConfigurations+1];

        TH1F *hDeltaPhi[nParticles][nRpConfigurations+1][nTotalTPCTOFSysStudies+1]; 
        TH2F *hMissingPtVsDeltaPhi[nParticles][nRpConfigurations+1][nTotalTPCTOFSysStudies+1];      

        struct hGroup {
            TH1F** hMain;
            TH2F** hVsPt;
            TH1F** hMainBcgSubtracted;
            TH1F* hBcgSysStudy;

            int part;
            int rpCon;
            TString hName;
            TString xLabel;
            TString yLabel;

            hGroup() : hMain(nullptr), hVsPt(nullptr), hBcgSysStudy(nullptr) {}
        };

        vector<hGroup> mainAnaHists;

        // eff
        TH1D *hTPCPhiEff[nParticles][nSigns];
        TH2D *hTPCEtaZEff[nParticles][nSigns];
        TH1D *hTPCPtEff[nParticles][nSigns];
        TH2D *hTPCEtaPhiEff[nParticles][nSigns];

        TH1D *hTOFPhiEff[nParticles][nSigns];
        TH2D *hTOFEtaZEff[nParticles][nSigns];
        TH1D *hTOFPtEff[nParticles][nSigns];
        TH2D *hTOFEtaPhiEff[nParticles][nSigns];

///////// PlotMainAna /////////
        void runAnaPlots();
        void initMainAnaPlots();
        void PlotMainAnaPlots();
        void loopThroughTree(UInt_t tpcAppMethod, UInt_t tofAppMethod, UInt_t nHitsFit, UInt_t nHitsDeDx, UInt_t pidMethod, UInt_t dcaMethod, int id);
        void correctMainPlotsForBinWidth();
        void subtractBackgroundFromMainPlots();
        void integrateCrossSection(hGroup mainHists, TString hName);
///////// PlotPID /////////
        void runM2Plots();
        void drawMissIdProbability(); 
        void CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin, Float_t rMargin, Float_t bMargin, Float_t tMargin);
        double XtoPad(double x);
        double YtoPad(double y);

///////// PlotBackgroundStudy /////////
        TH1D* backgroundStudy(TH2 *hist, int part, int rpCon, TString hName = "NonExclScaledBkgd_", int bcgStudy = -1);
        TH1D* bkgdHistogram(const TH2*, Double_t=exclusivityCut, Double_t=0.2, Double_t=0.28, Int_t=0, vector<TF1*> * = nullptr) const;
        TH1D* bkgdHistogramComparison(TH2* hMissPtVsX, int part, int rpCon, bool study);
        void subtractBackground(TH1*, const TH1*) const;
        Double_t bkgdEvents(const TH1*, Double_t=exclusivityCut, Double_t=0.2, Double_t=0.28, TF1* =nullptr) const;
        Double_t integratePtMiss(const TH1*, Double_t=exclusivityCut) const;

///////// PlotGetEfficiency /////////
        inline double getRetainCorrection() { return 1.39 * exp(-0.01468*mInstLumiPerRun[ mCurrentTree->getRunNumber() ]); }; // not use
        inline double getRpEff() { return rpEff[ mCurrentTree->getBranch(E) ]*rpEff[ mCurrentTree->getBranch(W) ]; };
        inline double getVertexCutEff() { return mVertexEffPerFill[mCurrentTree->getFillNumber()]; };
        inline double getVertexReconstructionEff(UInt_t dca = NOMINAL) { return vertexRecoEff[dca]; };
        inline double getPtMissEff() { return mPtMissEff[mCurrentTree->getPairID()][getRpCombination()]; };
        inline double getLuminosity() { return correctedIntegLum; };

        double getPIDEff(UInt_t var);
        double getTpcEff(UInt_t varApp, UInt_t nHitsFit, UInt_t nHitsDeDx);
        double getTofEff(UInt_t varApp); 
///////// PlotSysStudy /////////
        void PlotSysStudy();
        void DrawSysUncertainty(hGroup mainHists, bool drawCanvas = false);
        void PlotVariation(TH1* hist, const char* outputFileName); 
        void runSysStudy();

///////// PlotGraniitti /////////
        double graniittiScaleFactor[nParticles][nRpConfigurations];

        TH1F *hInvMassGran[nParticles][nRpConfigurations];
        TH1F *hDeltaPhiGran[nParticles][nRpConfigurations];
        TH1F *hInvMassGranCon[nParticles][nRpConfigurations];
        TH1F *hDeltaPhiGranCon[nParticles][nRpConfigurations]; 
        void initGraniittiPlots();
        TH1F* GetGraniittiPlot(int part, int rpCon, TString hName, double dataIntegral, bool scale, bool continuum = false);
        double getGraniittiSF(int part, int rpCon);
///////// PlotProbOfRetainEvent /////////
        std::map< int, int> mTotal;
        std::map< int, int> mPassed[nBranchesConfigurations]; 

        void runProbStudy(bool savePlots, TString plotSuffix);
        void fillMaps();
        void plotProbOfRetainEvent(TGraphAsymmErrors *gr, TF1 *expFit, TString name);
        bool isFarAway( UInt_t runNumber);
        bool skipRun( UInt_t runNumber);
        bool noVetoInRp(unsigned int branchConfig);
        void plotRetainSystStudy(unsigned int nPoints, double *x, double *y, unsigned int iConf, TString plotSuffix);

///////// PlotMissingPtPlots /////////
        void runMissingPt();
        void PlotMissingPt(int part, int rpCon);
        void DrawEffStudy(int part, int rpCon, TH1F* hEff[]);
        void DrawExclusiveLine(double xl, double yl, double xr, double yr);
        std::tuple<double, double, double> FitpXYMissing( TString name, int part, TH1 *hist);
        void GenereateMCpTDist(TH1 *hist, double* fitParam, int part, int rpCon, int* cepSignal, TH1F* histToFill[]);

///////// Central Embedding /////////
        unsigned int embedParticle;

        void runTpcEff();
        void runTofEff();
        void runEmbedVarDiff();
        void m2FromEmbedding();
        void drawVertexEtaSpace();
        void plotAverageEtaVertexEff(TH2 *hEff, TString histName);
        inline TString assosiationCut(int iPart) { return Form("TMath::Sqrt( TMath::Power(eta%i_TrueMc-eta%i, 2) + TMath::Power(phi%i_TrueMc-phi%i, 2)) < %f",iPart,iPart,iPart,iPart,deltaEtaPhiCut); };
        void runTpcEmbedSysStudy();
///////// RP plots /////////
        // MC and DATA comparison for elastic events
        void runRPMCComparison();
        // RP efficiency study
        void runRPEff();
        void plotGradientAndFindStableRegions(TH2D *hist, TString hName, int side);
        void runRPEffStudy();
///////// Elastic ana //////////////
        void vertexStudy();
        void elasticStudy();
///////// Trigger efficiency /////////

        // Efficiency plots
        TH1D *hDSMEff[2];

        enum DSMEFF { kRPET = 1, kRPIT, kRPETVETO, kRPITVETO, kTOFA, kTOFB, kBBCE, kBBCW, kBBCLE, kBBCLW, kZDCE, kZDCW, nDSMEFF };

        void CheckDSMEff(RecTree* recTree);
        void DrawDataSetInfo(TString dataSetName);
        void CalculateEfficiency(TString effString, unsigned int nPassed, unsigned int nTotal);

///////// Usefull plotting functions /////////
        void CreateCanvas(TCanvas **can, TString canName, double width = 1920.0, double heigth = 1080.0);
        void WriteCanvas(TString canName = "", TCanvas *can = nullptr);
        void SetAxisStyle(TAxis *axis);
        void SetAxisStyle(TGaxis *axis);
        void SetHistStyle(TH1* hist, Int_t color = mainColor, Int_t markStyle = mainMarker);
        void SetHistBcgStyle(TH1* hist, Int_t color = bckgColor, Int_t markStyle = bckgMarker);
        void SetGraphStyle(TGraph *graph, Int_t color = mainColor, Int_t markStyle = mainMarker);
        void CreateLegend(double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
        void CreateText(double xl, double yl, double xr, double yr);
        void SetGPad(double xl = 0.09, double yl = 0.16, double xr = 0.09, double yr = 0.015); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
        void DrawFiducial(int side, bool elasticFV = false);
        void CreateLine(double xl, double yl, double xr, double yr);
        void DrawRatioPlot(TH1 *h1, TH1 *h2, TH1 **h3);

        void DrawSTARInternal(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
        void DrawSTAR(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
        void DrawMainText(int part, int rpCon);            
        void DrawSystemDescription(int part = -1, double xl = 0.35, double yl = 0.88, double xr = 0.88, double yr = 0.95);
        //void DrawForwardProtonKin(double xl = 0.64, double yl = 0.63, double xr = 0.9, double yr = 0.84);
        void DrawForwardProtonKin(int rpCon = nRpConfigurations, double xl = 0.6, double yl = 0.68, double xr = 0.9, double yr = 0.84);
        void DrawCentralKin(int part, double xl = 0.47, double yl = 0.63, double xr = 0.62, double yr = 0.84);
        void SetPalletRange(TH1* hist, double min = 0.88, double max = 0.9);

        vector<float> getBinsVectorF(const TAxis *) const;
        vector<double> getBinsVectorD(const TAxis *) const;
        double* setBinArray(int nBins, double min, double max);

        unsigned int getRpCombination();
        double getRpDeltaPhi();
        void HighlightBin(TH1 *hist, int bin);
        void DrawLine(TH1 *hist, double x, bool negativeRange = false);
};

#endif
