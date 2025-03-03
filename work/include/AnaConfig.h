#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

const unsigned int nAligIteration = 4;
const bool DEBUG = false;

const int nPlanesUsed = 3;
const bool AFTERBURNER = true;
 
const TString studyName[] = { TString("MainAna"), TString("VertexStudy"), TString("Embedding"), 
            TString("TofQA"), TString("TofEff"), TString("TrigEff"), TString("FullZB"), TString("ElasticAna"), 
            TString("RpMcAna"), TString("Alignment") };
const TString nameOfTree[] = { TString("recTree"), TString("vertexTree"), TString("embedQATree"), TString(""), TString("tofEffTree"),
            TString("TrigEffTree"), TString("FullZBTree"), TString("ElasticAnaTree"), TString("RpMcAna"), TString("")}; // name MCZBTree set in Util
const bitset<8> treeBits[] = { bitset<8>(string("00011111")), bitset<8>(string("10001011")), bitset<8>(string("10001110")), 
            bitset<8>(string("00000000")), bitset<8>(string("000010011")), bitset<8>(string("01100011")), bitset<8>(string("01100001")), 
            bitset<8>(string("01110001")), bitset<8>(string("01110011")), bitset<8>(string("00000000")) };
/*
    MainAna    TString("recTree")  bitset<8>(string("00011111"))
    VertexStudy     TString("vertexTree")   bitset<8>(string("10001011"))
    Embedding   TString("embedQATree")  bitset<8>(string("10001110"))
    TofQA       TString("") bitset<8>(string("00000000"))
    TofEff      TString("tofEffTree")   bitset<8>(string("000010011"))
    TrigEff     TString("TrigEffTree")  bitset<8>(string("01100011"))
    FullZB      TString("FullZBTree")   bitset<8>(string("01100001"))
    ElasticAna  TString("ElasticAnaTree")   bitset<8>(string("01110001"))
    RpMcAna     TString("RPMCAnaTree")  bitset<8>(string("01110011"))
    Alignment   TString("") bitset<8>(string("00000000"))
*/



const TString inputDir = "/gpfs01/star/pwg_tasks/upc02/CEP_input_files/";
const TString pathToLumiFile = inputDir + "lists/luminosityForZB.list";
const TString pathToIntegLumiFile = inputDir + "lists/lumPerRun.list";
const TString pathToEffRootFile = inputDir + "eff.root";
const TString pathToOffSetFile = inputDir + "OffSetsRun17.list";
const TString pathToOffSetCorrectionFile = inputDir + "OffSetsCorrectionsRun17.list";


const TString embedPathPrefix = "/gpfs01/star/pwg/truhlar/Run17_P20ic/";
const TString embedPathSuffix = "/merged/StRP_production_0000.root";

const TString graniittiFilePath = "/star/data01/pwg_tasks/upc02/Graniitti/graniittiPlots.root";

// To plot PID eff and study TPC and TOF eff
const TString embedFile[] = { TString("piEmbed"), TString("kEmbed"), TString("pEmbed") };
const TString embedSysStudyFiles[] = { TString("nFitLoose"), TString("nFitTight"),
                            TString("nDedxLoose"), TString("nDedxTight")};
//const TString embedSysStudyFiles[] = { TString("nHitsFitLoose"), TString("nHitsFitTight"),
//                            TString("nHitsDedxLoose"), TString("nHitsDedxTight")};

// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 
// 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL 
// 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const std::vector<int> triggerID = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const std::vector<int> CEPtriggers = { 570705, 590705 }; 
const std::vector<int> ZBtriggers = { 570704 };
const std::vector<int> ElasticTriggers = { 570709, 570719, 590709 };
const std::vector<int> noTriggers = {}; // for MC analysis

// RP MC 
const string pathToMC = "/star/data01/pwg_tasks/upc02/run17MC/";
const double maxMcToRecoDist = 0.08; // max distance between MC generated and reconstructed is 80 MeV in px, py space
const double maxDiffBetweenRpHitMCAndReco = 0.01; // max distance between MC generated and reconstructed hit in RP in meters

// RP Eff
const double maxDistBetweenFitAndHit = 0.012; // max distance between found hit in RP under study and the fit

// Elastic event conditions:
const unsigned int nSigma = 3;
const double sigma = 0.00013; // rad, for Main period
const double PHI_CENT  = 90.0*0.0174533; // geometrical window used in elastic paper 510 GeV
const double PHI_WIDTH = 12.0*0.0174533; // geometrical window used in elastic paper 510 GeV
const double DCUT = 3*0.0007; // in meters

// TPC good track quality cuts
const int minNHitsFit[] = { 20, 17, 22}; //{ 20, 17, 22}; // 17 for loose, 20 for standard, 22 for tight cut
const int minNHitsDEdx[] = { 15, 12, 17}; // 12 for loose, 15 for standard, 17 for tight cut
const double minPt[] = {0.25, 0.3, 0.3}; // for pion, kaon and proton in GeV
const double minPtPair[] = {-999.0, 0.7, 1.1}; // // for kaon and proton in GeV
const double maxDcaXY[] = {1.5, 2.0, 1.25}; // cm
const double minDcaZ[] = {-1.0, -1.5, -0.8}; // cm
const double maxDcaZ[] = {1.0, 1.5, 0.8}; // cm
const double maxEta = 0.9;

// eta-zVertex phase-space limitation eta = - 1/250 * zVrtx -+ 0.9 = a*zVrtx + b
const double etaVertexSlope = -1/250.0;
const double etaVertexShift = 0.9;

// vertex range in z-coordinate
const double vertexRange = 100.0; // cm
const double zVertexFitRange = 120; // cm

// Exclusivity cut
const bool keepNonExclusive = true;
const double exclusivityCut = 0.12; // GeV

// pT missing study
const Long64_t nMCGenerratedExclusiveEvents = 200000;
// for pTMissing fit in pT missing fit study and in background subtraction via bin-by-bin method
const double nonExclusiveFitRangeMin[] = { 0.14, 0.14, 0.14};
const double nonExclusiveFitRangeMax[] = {0.3,0.35,0.35};
// for background subtraction study
const double nonExclFitRangeMin[] = { 0.16, 0.16, 0.16}; // pion, kaon, proton
const double nonExclFitRangeMax[] = {0.22,0.24,0.3};

// Central hadrons PID cuts
const double m2minKaons[] = {0.15, 0.1, 0.2}; // GeV^2
const double m2minProtons[] = {0.6, 0.55, 0.7}; // GeV^2
const double minChiSq[] = {9, 12, 7};
const double minChiSqPions[] = {12, 15, 9};

// Central embedding
const double deltaEtaPhiCut = .15;

//Forward proton kinematics
const double fpPxMin[] = {-0.23, -0.25, -0.21, -0.19}; // GeV
const double fpPyMin[] = {0.42, 0.48, 0.46, 0.46}; // GeV
const double fpPyMax[] = {0.86, 0.84, 0.84, 0.88};// GeV
const double fpPxCenter[] = {0.64, 0.7, 0.6, 0.7}; // GeV
const double fpPRadius[] = {1.36, 1.5, 1.3, 1.5}; // GeV^2
const double fpLPxCenter[] = {0.0, -0.25, -0.28, 0.0}; // GeV
const double fpLPxMax[] = {0.0, 0.06, 0.03, 0.0}; // GeV for better plotting
const double fpLPRadius[] = {0.0, 0.959, 0.946, 0.0}; // GeV^2

//Forward proton kinematics for the limited elastic fiducial volume
const double fpElPxMin = -0.17; // GeV
const double fpElPxMax = 0.21; // GeV
const double fpElPyMin = 0.48; // GeV
const double fpElPyMax = 0.74; // GeV

// Plot setting
const double titleSize = 45;//50;
const double textSize = 40;//45;
const int fontStyle = 43;
const int capitalStyle = 63;

const int mainColor = 4;
const int bckgColor = 2;

const int mainMarker = 20;
const int bckgMarker = 25;
const double markerSize = 2.0;

const int lineStyle = 1;
const double lineWidth = 2.0;
const int lineColor = 1;

// Usefull constants
const double sqrtOfTwo = 1.41421356237;
const double convertToDegree = 57.2957795;

// main histograms' settings
const int massNBins[] = {84, 27, 10};
//const int massNBins[] = {43, 27, 10};

//const double massBins[3][100] = {{0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6},
const double massBins[3][100] = {{0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.1, 3.2, 3.3, 3.4, 3.5},
                            {0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.1, 2.2, 2.3, 2.4},    
                            {1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5}}; 


                            
                       
const unsigned int minRun = 83026; // 180083026
const unsigned int maxRun = 148048; // 180148048


// Collider settings: default values for pp run17 0
//Forward proton reconstruction
const double vertex_X = -0.0010; // meters \pm 0.0004 // -0.00105 the old value
const double vertex_Y = -0.0014; // meters \pm 0.0002 // 0.0011 the old value
const double vertex_Z = .0; // meters
const double beamTiltAngle_X = -0.000070; // rad
const double beamTiltAngle_Y = -0.000070; // rad

const double mXYZ_IP[] = { vertex_X, vertex_Y, vertex_Z}; /* collision coordinates at the IP; 0=X, 1=Y, 2=Z */
const double mThetaXY_tilt[] = { beamTiltAngle_X, beamTiltAngle_Y}; /* tilt angles of the beam at collision; 0=X, 1=Y */
const double mDistanceFromIPtoDX[] = { 9.8, 9.8}; /* distance from the IP to the DX magnet in the East and West; 0=E, 1=W */
const double mLDX[] = { 3.7, 3.7};     /* length of DX in the East and West; 0=E, 1=W */
const double mBendingAngle[] = { 0.018832292, 0.018826657}; 


// different TPC, TOF, PID corr. and applictaion method
constexpr unsigned int nTPCAppSysStudies = 3; // 4 different ways how to apply the TPC eff
constexpr unsigned int nTOFAppSysStudies = 3; // 4 different ways how to apply the TPC eff
constexpr unsigned int nTPCnHitsStudies = 4; // 4 TPC nHits variations
constexpr unsigned int nPidVariation = 2; // Loose and Tight
constexpr unsigned int nDcaVariation = 4; // Loose and Tight

// Corrections and efficiencies
const double rpEff[] = {0.898, 0.870, 0.896, 0.892};
const double vertexRecoEff[] = {0.8964, 0.9153, 0.8693, 0.8967, 0.8784};

// sys. uncertainities of corrections for each particle type without luminosity
const double relativEffSysLow[] = {0.050, 0.061, 0.064};
const double relativEffSysHigh[] = {0.085, 0.092, 0.094};

// Corrected integrated luminosity
const double correctedIntegLum = 26.2541673925 * 1000; // convert from pb-1 to nb-1
const double relLumUncertainty = 0.064;

// pT smearing parameters
const double smearParam[6][5] = {{0.01116, -0.00154, 0.00035, 0.00921, 0.00051}, // pionPlus
                            {0.01250, -0.00205, 0.00039, 0.00704, 0.00107}, // pionMinus
                            {-0.00532, 0.00460, 0.00016, 0.02284, -0.00275}, // kaonPlus
                            {0.01622, -0.00512, 0.00169, 0.00333, 0.00268}, // kaonMinus
                            {-0.03793, 0.02444, -0.00220, 0.03921, -0.00501}, // protonPlus
                            {-0.01580, 0.01403, -0.00057, 0.02112, -0.00032}}; // protonMinus

#endif //ifndef CONFIG
