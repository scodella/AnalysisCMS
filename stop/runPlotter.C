#include "../test/HistogramReader.h"

// Constants
//------------------------------------------------------------------------------
const Bool_t datadriven = false;
const Bool_t allplots  = false;
const Bool_t dosystematics = true;
const Bool_t postfitplots = true;
const Bool_t paperstyle = true;
const Bool_t regionlegend = true;
const Bool_t relativeratio = true;

//const TString inputdir  = "../rootfiles/nominalFullStatus/";
//const TString inputdir  = "../minitrees/rootfiles/Zpeak_kinematic/";
//const TString inputdir  = "../minitrees/rootfilesOct17/nominal/";
//const TString inputdir  = "../minitrees/rootfilesFakePM/nominal/";
//const TString inputdir  = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/nominal/";
//const TString inputdir  = "../minitrees/rootfiles/WZtoWWptcut20_veto/";
//const TString inputdir  = "../minitrees/rootfiles/WZtoWWvetoMetCorr2/";
//const TString inputdir  = "../rootfiles/fake/";
//const TString inputdir  = "../rootfiles/SS/";
//const TString inputdir  = "../rootfiles/SSp/";
//const TString inputdir  = "../rootfiles/SSm/";
//const TString inputdir  = "../minitrees/rootfiles/ttZalexis/";
//const TString inputdir  = "../minitrees/rootfilesOct17/WZ3L/";
//const TString inputdir  = "../minitrees/rootfiles/ttZWZ3LNoMetCut/";
//const TString inputdir  = "../minitrees/rootfiles/ZZ/";
//const TString inputdir  = "../minitrees/rootfiles3R/ZpeakkfM/";
//const TString inputdir  = "../minitrees/rootfiles3R/ZZkfM/";
//const TString inputdir  = "../minitrees/rootfiles3R/ZZkfPt/";
//const TString inputdir  = "../minitrees/rootfiles3R/ZZkfdPhi/";
//const TString inputdir  = "../minitrees/rootfilesOct17/Zpeak_ptll/";
//const TString inputdir  = "../minitrees/rootfilesOct17/ZpeakMetCorr/";
//const TString inputdir  = "../minitrees/rootfiles/ZpeakDYcorrections/";
////const TString inputdir  = "../minitrees/rootfilesOct17/WZtoWWveto/";
//const TString inputdir  = "../minitrees/rootfilesFakePM/nominal/";
//const TString inputdir  = "../minitrees/rootfilesNoTopPt/nominal/";
//const TString inputdir  = "../minitrees/rootfilesCWR/nominal/";
//const TString inputdir  = "../minitrees/rootfilesOct17/CWRkinematic/";
//const TString inputdir  = "../minitrees/rootfilesFakePM/kinematiclowMet/";
//const TString inputdir  = "../minitrees/rootfiles/kinematic/";
//const TString inputdir  = "../../PlotsConfigurations/Configurations/T2tt/DatacardsTChiWWint/MassPointChiWW_mChi1_Xm288_Xm1/Postfit/";
//const TString inputdir  = "../../PlotsConfigurations/Configurations/T2tt/DatacardsPostfitPaperV2/MassPoint2tt_mStop-350to400_Sm350_Xm225/Postfit/"; const TString signal = "T2tt";
//const TString inputdir  = "/eos/cms/store/user/scodella/Stop/rootfiles/PostfitPaperV2/MassPoint2tt_mStop-350to400_Sm350_Xm225/Postfit/"; const TString signal = "T2tt";
const TString inputdir  = "/eos/cms/store/user/scodella/Stop/rootfiles/PostfitPaperV2/MassPointChiSlep_Xm500_Xm200/Postfit/"; const TString signal = "TChi";
//const TString inputdir  = "../../PlotsConfigurations/Configurations/T2tt/DatacardsTChiWWPostfitZeroStat2/MassPointChiWW_Xm100_Xm1/Postfit/"; const TString signal = "TChiWWlow";

const TString outputdir = "figures/";

//const TString signal = "";

const TString sl  = "#font[12]{l}";
const TString sll = "#font[12]{ll}";
const TString sm  = "#font[50]{p}_{T}^{miss}";
const TString pt  = "#font[50]{p}_{T}";
const TString mt2 = "#font[50]{M}_{T2}";

enum {linY, logY};

TString LSP = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
TString CHR = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
TString STP = "#tilde{t}"; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// runPlotter
//
// option = "hist"         --> all distributions normalized to the luminosity
// option = "nostack,hist" --> signal and top distributions normalized to one
//
//   Draw(TString  hname,                  Name of the histogram.
//        TString  xtitle       = "",      Title of the x-axis.
//        Int_t    ngroup       = -1,      Number of bins to be merged into one bin.
//        Int_t    precision    = 0,       Number of decimal digits.
//        TString  units        = "NULL",  Units of the histogram.
//        Bool_t   setlogy      = false,   Set it to true (false) for logarithmic (linear) scale.
//        Bool_t   moveoverflow = true,    Set it to true to plot the events out of range.
//        Float_t  xmin         = -999,
//        Float_t  xmax         = -999,
//        Float_t  ymin         = -999,
//        Float_t  ymax         = -999);
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void runPlotter(TString level,
		TString option = "hist")
{
  gInterpreter->ExecuteMacro("../test/PaperStyle.C");

  TString tok;

  Ssiz_t from = 0;

  TString analysis = (level.Tokenize(tok, from, "/")) ? tok : "NONE";

  if (analysis.EqualTo("NONE")) return;

  float lumi = lumi_fb_Full2016;
  //float lumi = lumi_fb_2016_susy;

  Bool_t scale = logY;

  int firstchannel = ee;
  int lastchannel  = ll;

  HistogramReader plotter(inputdir + analysis, outputdir);

  plotter.SetStackOption(option);
  plotter.SetPublicStyle(false);
  if (paperstyle) {
    plotter.SetSavePdf(true); 
    plotter.SetIsPreliminary(false);
    //plotter.SetDynamicRatioAxis(true);
  }
  if (regionlegend) plotter.SetDrawRegionLegend(true);

  if (option.Contains("nostack"))
    {
      plotter.SetDrawRatio(false);
    }
  else
    {
      plotter.SetLuminosity(lumi, postfitplots);
      plotter.SetDrawRatio (true);
      plotter.SetDrawRatioRel (relativeratio);
    }
  
  float SF_ttZ = 1., SF_ZMet = 1., SF_DY = 1., SF_WZ = 1.;
  if (!postfitplots && !inputdir.Contains("fake") && !inputdir.Contains("SS") && !inputdir.Contains("ttZ") && !inputdir.Contains("WZ3L")) {
    if (level.Contains("_SR") || level.Contains("_VRggg")) {// && inputdir.Contains("DYcorr"))) {
      SF_ttZ = 1.44, SF_WZ = 0.97; 
      if (!postfitplots && !inputdir.Contains("Zpeakk") && !inputdir.Contains("ZZ")) {
	if (level.Contains("NoJet")) {
	  SF_ZMet = 1.;//0.74; // +/- 0.19
	} else if (level.Contains("NoTag")) {
	  SF_ZMet = 1.21;// +/- 0.17 
	} else { 
	  SF_ZMet = 1.06;// +/- 0.12
	} 
      }
    }
  }

  // Get the data
  //----------------------------------------------------------------------------
  plotter.AddData("01_Data", "Data", color_Data);
  //plotter.SetDrawRatio(false);

  TString DYCorr = "_DYcorr";
  if (!level.Contains("_SR")) DYCorr = "";
  if (inputdir.Contains("/WZ") || inputdir.Contains("fake") || inputdir.Contains("SS") || 
      inputdir.Contains("ttZ") || inputdir.Contains("/ZZ") || inputdir.Contains("/Zpeak") ||
      inputdir.Contains("Full") || inputdir.Contains("kinematic")) DYCorr = "";
  if (inputdir.Contains("CWR")) DYCorr = "_DYcorr";

  if (level.Contains("_DYnorm")) {
    SF_DY = 0.683211799801126896;
    level.ReplaceAll("_DYnorm", "");
  } else if (level.Contains("_DYcorrVR1norm")) {
    DYCorr = "_DYcorrVR1norm";
    level.ReplaceAll("_DYcorrVR1norm", "");
  }

  TString TopPt = "_NoTopPt";
  if (postfitplots) TopPt = "";

  
  std::cout << "DYCorr = " << DYCorr << " TopPt = " << TopPt << std::endl;

  // Add processes
  //----------------------------------------------------------------------------
  if (paperstyle) {
    plotter.AddProcess("98_Others",    "Minor bkg.",      color_VVV);
    plotter.AddProcess("07_ZJetsHT" + DYCorr + "_bis",     "Drell-Yan",   color_ZJets,  roc_background, SF_DY);
    plotter.AddProcess("03_ZZ",        "ZZ (#rightarrow 2l2#nu)",   color_VZ, roc_background, SF_ZMet);
    plotter.AddProcess("10_TTZ",       "t#bar{t}Z",      color_TTZ,  roc_background, SF_ttZ);
    plotter.AddProcess("02_WZTo3LNu",  "WZ (#rightarrow 3l)",  color_WZTo3LNu,  roc_background, SF_WZ);
    plotter.AddProcess("06_WW",        "WW",       color_WW, roc_background);
    plotter.AddProcess("05_ST",        "tW",       color_ST);
    plotter.AddProcess("04_TTTo2L2Nu" + TopPt, "t#bar{t}",       color_TTTo2L2Nu, roc_background);
  } else {
    //plotter.AddProcess("16_tZq",      "tZq",      49);
    //plotter.AddProcess("16_STs",      "ST s-ch.",      49);
    //plotter.AddProcess("17_STt",      "ST t-ch.",      46);
    //plotter.AddProcess("TTToSemiLepton", "t#bar{t} Semilep.",  49);
    ////plotter.AddProcess("14_HZ",        "HZ",       color_HZ);
    if (inputdir.Contains("SS/") || inputdir.Contains("ttZ") || inputdir.Contains("WZ3L") || 
	inputdir.Contains("/ZZ")) 
      plotter.AddProcess("13_VVV",      "VVV",      color_VVV);
    else 
      plotter.AddProcess("15_VZ3V",      "VVV + VZ",      color_VVV);
    //plotter.AddProcess("15_VZ",      "VZ",      color_VVV);
    if (inputdir.Contains("/ZZ") || inputdir.Contains("/WZ") || inputdir.Contains("/ttZ")) { 
      plotter.AddProcess("14_ZZTo4L",        "ZZ (#rightarrow 4l)",       color_ZZ4L,  roc_background);
      ////plotter.AddProcess("14a_ZZTo4L",        "qqZZ (#rightarrow 4l)",       49,  roc_background);// 1.256/1.212);
      ////plotter.AddProcess("14b_ZZTo4L",        "ggZZ (#rightarrow 4l)",       48,  roc_background);
      ////plotter.AddProcess("14c_ZZTo4L",        "H#rightarrow ZZ",       47,  roc_background);
    }
    if (inputdir.Contains("SS/")) 
      plotter.AddProcess("03_VZ",        "VZ (#rightarrow 2l)",       color_VZ,  roc_background);
    else if (!inputdir.Contains("/ZZ")) plotter.AddProcess("03_ZZ",        "ZZ (#rightarrow 2l2#nu)",   color_VZ, roc_background, SF_ZMet);
    ////plotter.AddProcess("03a_ZZ",        "qqZZ (#rightarrow 2l2#nu)",       48,  roc_background);
    ////plotter.AddProcess("03b_ZZ",        "ggZZ (#rightarrow 2l2#nu)",       47,  roc_background);
    //plotter.AddProcess("15_VZ",        "VZ (#rightarrow 2l2q)",       color_VZ2L2Q,  roc_background);
    //plotter.AddProcess("11_Wg",        "W#gamma",  color_Wg);
    //plotter.AddProcess("15_WgStar",    "W#gamma*", color_WgStar);
    plotter.AddProcess("09_TTW",       "t#bar{t}W",      color_TTV);
    plotter.AddProcess("10_TTZ",       "t#bar{t}Z",      color_TTZ,  roc_background, SF_ttZ);
    plotter.AddProcess("11_HWW",       "HWW",      color_HWW);
    if (!inputdir.Contains("/WZ")) plotter.AddProcess("02_WZTo3LNu",  "WZ (#rightarrow 3l)",  color_WZTo3LNu,  roc_background, SF_WZ);
    if (!inputdir.Contains("/WZ")) plotter.AddProcess("06_WW",        "WW",       color_WW, roc_background);
    plotter.AddProcess("05_ST",        "tW",       color_ST);
    plotter.AddProcess("07_ZJetsHT" + DYCorr,     "Z+jets",   color_ZJets,  roc_background, SF_DY);
    if (inputdir.Contains("SS")) plotter.AddProcess("TTToSemiLepton", "t#bar{t} Semilep.",  41, roc_background);
    plotter.AddProcess("04_TTTo2L2Nu", "t#bar{t}",       color_TTTo2L2Nu, roc_background);
    //else plotter.AddProcess("TTJets", "#bar{t}t",       color_TTTo2L2Nu);
    //if (inputdir.Contains("SS/")) plotter.AddProcess("WJetsToLNu", "WJets",      color_WJets);
    if (inputdir.Contains("/WZ")) plotter.AddProcess("02_WZTo3LNu",  "WZ (#rightarrow 3l)",  color_WZTo3LNu,  roc_background, SF_WZ);
    if (inputdir.Contains("/WZ")) plotter.AddProcess("06_WW",        "WW",       color_WW);
  }

  if (postfitplots) {
    //plotter.SetDynamicRatioAxis(true);
    plotter.AddPrefit("99_TotalBackground", "Pre-fit", 9);//color_Prefit);
    plotter.AddPostfit("00_Total", "Post-fit", 2);//kRed+2);
    plotter.AddPostfitSM("99_TotalBackground", "Post-fit", 2);//kRed+4);
  }
  
  if (signal=="T2tt") {
    
    if (paperstyle) {
      plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm225", STP+" "+STP+",   "+STP+"#rightarrow t"+LSP+",   (m_{"+STP+"} = 350 GeV,  m_{"+LSP+"} = 225 GeV)",kViolet);
    } else {
      plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm225", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,225)",kViolet);
      plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm175", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,175)",kRed);
      plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm263", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,263)",kOrange);
      plotter.AddSignal("T2tt_mStop-150to250_Sm250_Xm125", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (250,125)",kBlue);
      plotter.AddSignal("T2tt_mStop-400to1200_Sm450_Xm325","#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (450,325)",kCyan);
    }
    
  } else if (signal=="TChi") {

    if (paperstyle) {
      if (inputdir.Contains("Xm500_Xm200"))
	plotter.AddSignal("TChiSlep_Xm500_Xm200", CHR+CHR+",   "+CHR+"#rightarrow l#nu"+LSP+",   (m_{"+CHR+"} = 500 GeV,  m_{"+LSP+"} = 200 GeV)", kViolet);
      else if (inputdir.Contains("Xm100_Xm1"))
	plotter.AddSignal("TChiWW_Xm100_Xm1", CHR+CHR+",   "+CHR+"#rightarrow W"+LSP+",   (m_{"+CHR+"} = 100 GeV,  m_{"+LSP+"} = 1 GeV)", kViolet);
      else if (inputdir.Contains("Xm125_Xm1"))
	plotter.AddSignal("TChiWW_Xm125_Xm1", CHR+CHR+",   "+CHR+"#rightarrow W"+LSP+",   (m_{"+CHR+"} = 125 GeV,  m_{"+LSP+"} = 1 GeV)", kViolet);
      else if (inputdir.Contains("Xm150_Xm1"))
	plotter.AddSignal("TChiWW_Xm150_Xm1", CHR+CHR+",   "+CHR+"#rightarrow W"+LSP+",   (m_{"+CHR+"} = 150 GeV,  m_{"+LSP+"} = 1 GeV)", kViolet);
      else if (inputdir.Contains("Xm175_Xm1"))
	plotter.AddSignal("TChiWW_Xm175_Xm1", CHR+CHR+",   "+CHR+"#rightarrow W"+LSP+",   (m_{"+CHR+"} = 175 GeV,  m_{"+LSP+"} = 1 GeV)", kViolet);
      else 
	plotter.AddSignal("TChiSlep_Xm500_Xm200", CHR+CHR+",   "+CHR+"#rightarrow l#nu"+LSP+",   (m_{"+CHR+"} = 500 GeV,  m_{"+LSP+"} = 200 GeV)", kViolet);
    } else {
      plotter.AddSignal("TChiSlep_Xm200_Xm1",   "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (200,  1)", kRed);
      plotter.AddSignal("TChiSlep_Xm500_Xm200", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (500,200)", kViolet);
      plotter.AddSignal("TChiSlep_Xm800_Xm400", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (800,400)", kOrange);
      plotter.AddSignal("TChiWW_Xm200_Xm25",  "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (200, 25)", kBlue);
      plotter.AddSignal("TChiWW_Xm300_Xm100", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (300,100)", kCyan);
    }
    
  } else if (signal=="Mix") {

    plotter.AddSignal("T2bW_Sm250_Xm1", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (250,1)",kRed);
    plotter.AddSignal("T2bW_Sm250_Xm50", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (250,50)",kRed);
    plotter.AddSignal("T2bW_Sm200_Xm25", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (200,25)",kRed);
    //plotter.AddSignal("T2tt_mStop-150to250_Sm250_Xm125", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (250,125)",kRed);
    //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm225", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,225)",kRed);
    //plotter.AddSignal("T2tt_mStop-250to350_Sm300_Xm125", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (300,125)",kRed);
    //plotter.AddSignal("TChiSlep_Xm200_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (200,1)", kBlue);
    //plotter.AddSignal("TChiSlep_Xm500_Xm200", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (500,200)", kCyan);

  } else if (signal=="TChiWWlow") {

    plotter.AddSignal("TChiWW_all_Xm100_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (100, 1)", kRed);
    plotter.AddSignal("TChiWW_mChi1bis_Xm113_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (113, 1)", kCyan);
    plotter.AddSignal("TChiWW_Xm125_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (125, 1)", kViolet);
    plotter.AddSignal("TChiWW_mChi1_Xm138_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (138, 1)", kGreen+3);
    //plotter.AddSignal("TChiWW_Xm150_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (150, 1)", kGreen);
    //plotter.AddSignal("TChiWW_Xm175_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (175, 1)", kBlue);
    //plotter.AddSignal("TChiWW_Xm200_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (200, 1)", kCyan);
   
  } else if (signal=="TChiWWint") {

    plotter.AddSignal("TChiWW_mChi1_Xm238_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (238, 1)", kRed);
    plotter.AddSignal("TChiWW_Xm250_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (250, 1)", kViolet);
    plotter.AddSignal("TChiWW_mChi1_Xm263_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (263, 1)", kGreen+3);
    plotter.AddSignal("TChiWW_Xm275_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (275, 1)", kGreen);
    plotter.AddSignal("TChiWW_mChi1_Xm288_Xm1", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (288, 1)", kBlue);
   
  }else if (signal=="TChiSleplow") {

    plotter.AddSignal("TChiSlep_Xm200_Xm100",   "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (200,100)", kRed);
    plotter.AddSignal("TChiSlep_Xm200_Xm125", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (200,125)", kViolet);
    plotter.AddSignal("TChiSlep_Xm250_Xm125", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow l#nu#tilde{#chi}^{0}_{1} (250,125)", kOrange);
    plotter.AddSignal("TChiSlep_Xm250_Xm150",  "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (250, 150)", kBlue);
    plotter.AddSignal("TChiSlep_Xm300_Xm175", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (300,175)", kCyan);
   
  }
  
  if (inputdir.Contains("rootfiles/nominalX")) {

    // Draw events by cut
    //----------------------------------------------------------------------------
    
    plotter.SetDrawYield(false);
  
    gSystem->mkdir(outputdir + level, kTRUE);
    
    for (int i=firstchannel; i<=lastchannel; i++)
      {
	plotter.LoopEventsByCut(analysis, "h_counterLum_" + schannel[i]);
	
	TString title = (i < lastchannel) ? lchannel[i] : "inclusive";
	
	plotter.SetTitle(title);
	
	plotter.Draw(analysis + "/h_counterLum_" + schannel[i] + "_evolution", "", -1, 0, "NULL", logY, false);
      }
    
    
    // Draw events by channel
    //----------------------------------------------------------------------------
    plotter.SetDrawYield(false);
    
    for (int j=0; j<=njetbin; j++)
      {
	if (j != njetbin) continue;
	
	TString jetbin = (j < njetbin) ? Form("/%djet", j) : "";
	
	gSystem->mkdir(outputdir + level + jetbin, kTRUE);
	
	plotter.LoopEventsByChannel(level + jetbin);
	
	plotter.Draw(level + jetbin + "/h_counterLum_evolution", "", -1, 0, "NULL", scale, false);
      }
    
  }
  
  if (dosystematics) {
    if (postfitplots) {
      plotter.AddSystematic("Stop", "Postfit");
    } else {
      plotter.AddSystematic("Stop", "Statistics");
      plotter.AddSystematic("Stop", "Luminosity");
      plotter.AddSystematic("Stop", "Trigger");
      plotter.AddSystematic("Stop", "MT2llTop");
      plotter.AddSystematic("Stop", "MT2llWW");
      plotter.AddSystematic("Stop", "Fake");
      plotter.AddSystematic("Stop", "Idiso");
      plotter.AddSystematic("Stop", "JES");
      plotter.AddSystematic("Stop", "MET");
      plotter.AddSystematic("Stop", "PDF");
      plotter.AddSystematic("Stop", "Q2");
      plotter.AddSystematic("Stop", "Reco");
      plotter.AddSystematic("Stop", "Toppt");
      plotter.AddSystematic("Stop", "Isrnjet");
      plotter.AddSystematic("Stop", "Metfastsim");
      plotter.AddSystematic("Stop", "Pileup");
      plotter.AddSystematic("Stop", "Fastsim");
      plotter.AddSystematic("Stop", "BtagFS");
      plotter.AddSystematic("Stop", "Btag");
      plotter.AddSystematic("Stop", "Btaglight");
      plotter.AddSystematic("Stop", "ttZSF");
      plotter.AddSystematic("Stop", "WZSF");
      plotter.AddSystematic("Stop", "ZZSF");
      if (!inputdir.Contains("rootfilesCWR")) plotter.AddSystematic("Stop", "ZZshape");
      ////plotter.AddSystematic("Stop", "DYSF");
      plotter.AddSystematic("Stop", "DYshape");
      if (level.Contains("_SR")) plotter.AddSystematic("Stop", "DYnojet");
      ////plotter.AddSystematic("Stop", "normWZ");
      if (inputdir.Contains("rootfilesCWR")) plotter.AddSystematic("Stop", "normWW");
      if (inputdir.Contains("rootfilesCWR")) plotter.AddSystematic("Stop", "normTtbar");
      if (inputdir.Contains("rootfilesCWR")) plotter.AddSystematic("Stop", "normTW");
      plotter.AddSystematic("Stop", "normDY");
      ///plotter.AddSystematic("Stop", "normTTW");
      ///plotter.AddSystematic("Stop", "normHWW");
      ///plotter.AddSystematic("Stop", "normVVV");
      ////plotter.AddSystematic("Stop", "Muscale");
    }
  }
  
  // Draw distributions
  //----------------------------------------------------------------------------
  if (!option.Contains("nostack")) plotter.SetDrawYield(true);
  if (paperstyle) plotter.SetDrawYield(false);
  //plotter.SetDrawYield(false);

  float m2l_xmin   = (level.Contains("WZ")) ?  60 :   0;  // [GeV]
  float m2l_xmax   = (level.Contains("WZ")) ? 120 : 300;  // [GeV]
  int   m2l_ngroup = (level.Contains("WZ")) ?   2 :   5;

  for (int j=0; j<=njetbin; j++)
    {
      if (j != njetbin) continue;   
         
      TString jetbin = (j < njetbin) ? Form("/%djet", j) : "";

      gSystem->mkdir(outputdir + level + jetbin, kTRUE);

      TString prefix = level + jetbin + "/h_"; 

      for (int i=firstchannel; i<=lastchannel+1; i++)
	{
	  TString suffix = (i <= lastchannel) ? "_" + schannel[i] : "_sf";

	  TString title = (i < lastchannel) ? lchannel[i] : "inclusive";

	  if (suffix=="_sf") title = "ee+#mu#mu";

	  if (postfitplots) {
	    if (suffix=="_ee" || suffix=="_mm") continue;
	  } 

	  if (inputdir.Contains("/Zpeak")) 
	    if (suffix=="_sf" || suffix=="_em") continue;

	  if (inputdir.Contains("CWR")) 
	    if (suffix=="_ee" || suffix=="_mm") continue;

	  //if (suffix!="_ll") continue;

	  plotter.SetTitle(title);
	  
	  if (postfitplots) {
	    if (level.Contains("_SR3") && !inputdir.Contains("MassPointChi")) {
	      plotter.Draw(prefix + "MT2llisr_prefit"     + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	      plotter.Draw(prefix + "MT2llisr_fit_b"      + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	      plotter.Draw(prefix + "MT2llisr_fit_s"      + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    } else {
	      plotter.Draw(prefix + "MT2ll_prefit"     + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	      plotter.Draw(prefix + "MT2ll_fit_b"      + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	      plotter.Draw(prefix + "MT2ll_fit_s"      + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    }
	    continue;
	  }
	  
	  if (inputdir.Contains("DYcorrections")) {
	    plotter.Draw(prefix + "dphiLLbin1"     + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    plotter.Draw(prefix + "dphiLLbin2"     + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    plotter.Draw(prefix + "dphiLLbin3"     + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    plotter.Draw(prefix + "dphiLLbin4"     + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    plotter.Draw(prefix + "dphiLLbin5"     + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    continue;
	  }
	  
	  // Common histograms
	  //--------------------------------------------------------------------
	  if (inputdir.Contains("/ttZ/")) {
	    plotter.Draw(prefix + "Counter"        + suffix, "m_{ll} (" + sll + ")",               1, 0, "GeV",  linY, false, 80, 100);
	    continue;
	  }
	  if (inputdir.Contains("noZcut")) {
	    plotter.Draw(prefix + "M1ll"     + suffix, "m_{ll}",                                  3, 0, "GeV",  linY, false, 60,  120);
	    plotter.Draw(prefix + "M2ll"     + suffix, "m_{ll}",                                  4, 0, "GeV",  linY, false,  0,  180);
	    continue;
	  } 
	  /*plotter.Draw(prefix + "m2l" + suffix, "m_{" + sll + "}", m2l_ngroup, 0, "GeV", linY, true, m2l_xmin, m2l_xmax);
	  plotter.Draw(prefix + "lep1pt"         + suffix, "leading lepton p_{T}",              10, 0, "GeV",  scale, true, 0,  250);
	  plotter.Draw(prefix + "lep2pt"         + suffix, "trailing lepton p_{T}",              5, 0, "GeV",  scale, true, 0,  150);
	  plotter.Draw(prefix + "lep1eta"        + suffix, "leading lepton #eta",               -1, 1, "NULL", scale);
	  plotter.Draw(prefix + "lep2eta"        + suffix, "trailing lepton #eta",              -1, 1, "NULL", scale);
	  plotter.Draw(prefix + "lep1phi"        + suffix, "leading lepton #phi",                5, 2, "rad",  scale);
	  plotter.Draw(prefix + "lep2phi"        + suffix, "trailing lepton #phi",               5, 2, "rad",  scale);
	  plotter.Draw(prefix + "nvtx"           + suffix, "number of vertices",                -1, 0, "NULL", scale, true, 0,   30);*/
	  plotter.Draw(prefix + "MT2ll"        + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	  //plotter.Draw(prefix + "MT2ll"        + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  linY, false, 0, 140);
	  if (level.Contains("_SR3") && signal=="T2tt")  plotter.Draw(prefix + "MT2llisr"        + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	  //if (level.Contains("_SR3")) plotter.Draw(prefix + "MT2llisr"        + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	  
	  if (inputdir.Contains("Oct17/nominal/")) continue;
	  if (inputdir.Contains("NoTopPt/nominal/")) continue; 
	  if (dosystematics && !inputdir.Contains("rootfilesCWR")) continue;
	  int METrebin = 1;
	  if (level.Contains("01_")) METrebin = 1;
	  else if (inputdir.Contains("/ZZ") || inputdir.Contains("/WZ")) METrebin = 4;
	  float minMET = (level.Contains("_SR")) ? 140. : 0.;
	  plotter.Draw(prefix + "MET"     + suffix, sm, METrebin, 0, "GeV",  scale, true, minMET,  400);
	  //plotter.Draw(prefix + "njet20"         + suffix, "number of jets",                    -1, 0, "NULL", linY, true, 0., 2.);
	  //plotter.Draw(prefix + "mt2ll"   + suffix, mt2 + "(" + sll + ")",         10, 0, "GeV",  scale, true, 0, 400);continue;
	  if (inputdir.Contains("gkinematic")) {
	    plotter.Draw(prefix + "dphiisrmet"         + suffix, "#Delta#phi(ISR jet," + sm + ")",             8, 0, "NULL", linY, false);
	  }
	  if (inputdir.Contains("../rootfilesX/")) {
	    plotter.Draw(prefix + "lep1pt"         + suffix, "leading lepton p_{T}",              10, 0, "GeV",  scale, true, 0,  250); // LS
	    plotter.Draw(prefix + "lep2pt"         + suffix, "trailing lepton p_{T}",              5, 0, "GeV",  scale, true, 0,  150);
	    plotter.Draw(prefix + "lep1eta"        + suffix, "leading lepton #eta",               -1, 1, "units", scale);
	    plotter.Draw(prefix + "lep2eta"        + suffix, "trailing lepton #eta",              -1, 1, "units", scale);
	    plotter.Draw(prefix + "lep1phi"        + suffix, "leading lepton #phi",                5, 2, "rad",  scale);
	    plotter.Draw(prefix + "lep2phi"        + suffix, "trailing lepton #phi",               5, 2, "rad",  scale);
	    plotter.Draw(prefix + "nbjet30csvv2m"  + suffix, "number of 30 GeV csvv2m b-jets",    -1, 0, "units", scale);
	    plotter.Draw(prefix + "njet"           + suffix, "number of jets",                    -1, 0, "units", scale, true);
	  } else if (inputdir.Contains("kinematic") || inputdir.Contains("rootfilesCWR")) {
	    plotter.Draw(prefix + "Lep1Pt"         + suffix, "leading lepton " + pt,              10, 0, "GeV",  scale, true, 0,  250); // LS
	    plotter.Draw(prefix + "Lep2Pt"         + suffix, "trailing lepton " + pt,              5, 0, "GeV",  scale, true, 0,  150);
	    plotter.Draw(prefix + "JetPt"          + suffix, "jet " + pt,                          5, 0, "GeV",  scale, true, 0,  250);
	    plotter.Draw(prefix + "nbjet"          + suffix, "number of b-tagged jets",           -1, 0, "units", scale, true, -0.5, 3.5);
	    plotter.Draw(prefix + "njet20"         + suffix, "number of jets",                    -1, 0, "units", scale, true, -0.5, 6.5);
	    plotter.Draw(prefix + "nbjet"          + suffix, "number of b-tagged jets",           -1, 0, "units", linY, true, -0.5, 3.5);
	    plotter.Draw(prefix + "njet20"         + suffix, "number of jets",                    -1, 0, "units", linY, true, -0.5, 6.5);
	    plotter.Draw(prefix + "njetISR"        + suffix, "number of ISR jets",                -1, 0, "units", scale, true, -0.5, 1.5);
	    plotter.Draw(prefix + "njetISR"        + suffix, "number of ISR jets",                -1, 0, "units", linY, true, -0.5, 1.5);
	    int rebinISR = (level.Contains("_SR3")) ? 4 : 2;
	    plotter.Draw(prefix + "ptjetISR"       + suffix, "ISR jet " + pt,               rebinISR, 0, "GeV",  scale, true, 0, 600.);
	    plotter.Draw(prefix + "dphijetISR"     + suffix, "#Delta#phi(ISR jet, " + sm +")",    10, 2, "rad", scale);
	    plotter.Draw(prefix + "mt2LL"          + suffix, mt2 + "(" + sll + ")",               20, 0, "GeV",  scale, true, 0, 200);
	    if (inputdir.Contains("kinematic")) {
	      plotter.Draw(prefix + "Jet1Pt"         + suffix, "leading jet " + pt,                 10, 0, "GeV",  scale, true, 0,  250); 
	      plotter.Draw(prefix + "Jet2Pt"         + suffix, "trailing jet " + pt,                 5, 0, "GeV",  scale, true, 0,  150);
	      plotter.Draw(prefix + "Lep1Eta"        + suffix, "leading lepton #eta",               -1, 1, "NULL", scale);
	      plotter.Draw(prefix + "Lep2Eta"        + suffix, "trailing lepton #eta",              -1, 1, "NULL", scale);
	      plotter.Draw(prefix + "Lep1Phi"        + suffix, "leading lepton #phi",                5, 2, "rad",  scale);
	      plotter.Draw(prefix + "Lep2Phi"        + suffix, "trailing lepton #phi",               5, 2, "rad",  scale);
	    }
	    continue;
	  } else if (inputdir.Contains("Full")) {
	    plotter.Draw(prefix + "m2l" + suffix, "m_{" + sll + "}", -1, 0, "GeV", scale, true, 81., 101.);
	    continue;
	  }
	  //continue;
	  //plotter.Draw(prefix + "Counter"        + suffix, "m_{ll} (" + sll + ")",               1, 0, "GeV",  linY, false, 80, 100);
	  //plotter.Draw(prefix + "MT2_Met"      + suffix, mt2 + "-Met",                         1, 0, "GeV",  scale, false);
	  //continue;
	  
	  if (inputdir.Contains("../rootfiles/nominal")) {
	    plotter.Draw(prefix + "mt2ll"   + suffix, mt2 + "(" + sll + ")",         10, 0, "GeV",  scale, true, 0, 400);
	    continue;
	  }
	  if (inputdir.Contains("ZZ")) {
	    plotter.Draw(prefix + "M1ll"     + suffix, "m_{ll}",                                  1, 0, "GeV",  linY, false, 60,  120);
	    plotter.Draw(prefix + "M2ll"     + suffix, "m_{ll}",                                  1, 0, "GeV",  linY, false, 60,  120);
	  } else if (inputdir.Contains("Zpeak")) {
	    //plotter.Draw(prefix + "HT"     + suffix, "#sum_{jet} p_{T}",                                  2, 0, "GeV",  scale, true, 0,  800);
	    //plotter.Draw(prefix + "HT"     + suffix, "#sum_{jet} p_{T}",                                  2, 0, "GeV",  linY, true, 0,  800);
	    //plotter.Draw(prefix + "MT2ll_HTm200" + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    //plotter.Draw(prefix + "MT2ll_HTp200" + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    //plotter.Draw(prefix + "MT2ll_HTm150" + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    //plotter.Draw(prefix + "MT2ll_HTp150" + suffix, mt2 + "(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	    //plotter.Draw(prefix + "njet20dphilmet"           + suffix, "number of 20 GeV jets",             -1, 0, "NULL", scale);
	    //plotter.Draw(prefix + "njet20dphilmet"           + suffix, "number of 20 GeV jets",             -1, 0, "NULL", linY);
	    plotter.Draw(prefix + "ptLL"              + suffix, "#font[50]{p}_{T}^{ll}",             1, 0, "GeV", scale);
	    plotter.Draw(prefix + "ptLLbins"              + suffix, "#font[50]{p}_{T}^{ll}",             1, 0, "GeV", scale);
	    plotter.Draw(prefix + "dphillMET"           + suffix, "#Delta#phi(" +sll + "," + sm + ")",             10, 0, "NULL", scale);
	    //plotter.Draw(prefix + "dphillMET"           + suffix, "#Delta#phi(" +sll + "," + sm + ")",             10, 0, "NULL", linY);
	    plotter.Draw(prefix + "dphiLL"              + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", scale);
	    plotter.Draw(prefix + "dphiLL"              + suffix, "#Delta#phi(lep1,lep2)",             10, 0, "NULL", linY);
	    continue;
	  }
	  //plotter.Draw(prefix + "lep1pt"      + suffix, "lep1pt",                         1, 0, "GeV",  scale, true, 0, 200);
	  //plotter.Draw(prefix + "maxjetpt"      + suffix, "leading jet pt",                  2, 0, "GeV",  scale, false, 0, 600);
	  //plotter.Draw(prefix + "dphiminlepmet"    + suffix, "#Delta#phi(lep,E_{T}^{miss})",      10, 2, "rad",  linY, false);
	  //plotter.Draw(prefix + "dphillMET"    + suffix, "#Delta#phi(ll,E_{T}^{miss})",      10, 2, "rad",  linY, false);
	  //plotter.Draw(prefix + "njet20"           + suffix, "number of jets",             -1, 0, "NULL", scale, true, 0, 7);
	  //plotter.Draw(prefix + "njet30"           + suffix, "number of 30 GeV jets",             -1, 0, "NULL", scale);
	  //plotter.Draw(prefix + "njet20"           + suffix, "number of jets",             -1, 0, "NULL", linY, true, 0, 7);
	  //plotter.Draw(prefix + "njet30"           + suffix, "number of 30 GeV jets",             -1, 0, "NULL", linY);
	  if (inputdir.Contains("../rootfiles/nominal")) {
	    //plotter.Draw(prefix + "jet1pt"      + suffix, "jep1pt",                         10, 0, "GeV",  scale, true, 0, 200);
	    //plotter.Draw(prefix + "m2l" + suffix, "m_{" + sll + "}", m2l_ngroup, 0, "GeV", logY, true, m2l_xmin, m2l_xmax);
	    plotter.Draw(prefix + "dphillmet"      + suffix, "#Delta#phi(" +sll + "," + sm + ")",  10, 2, "rad",  scale); 
	    plotter.Draw(prefix + "dphillmet"      + suffix, "#Delta#phi(" +sll + "," + sm + ")",  10, 2, "rad",  linY); 
	    plotter.Draw(prefix + "njet"           + suffix, "number of jets",             -1, 0, "NULL", scale);
	    plotter.Draw(prefix + "njet"           + suffix, "number of jets",             -1, 0, "NULL", linY);
	    //plotter.Draw(prefix + "HTvisible_Met"+ suffix, "HTvisible-Met",                      1, 0, "GeV",  scale, false);
	    //plotter.Draw(prefix + "ht"           + suffix, "H_{T}",                             20, 0, "GeV",  scale, true, 0, 1500);
	    //plotter.Draw(prefix + "htvisible"    + suffix, "#sum_{jet,lepton} p_{T}",           20, 0, "GeV",  scale, true, 0, 1500);
	    //plotter.Draw(prefix + "htjets"       + suffix, "#sum_{jet} p_{T}",                  20, 0, "GeV",  scale, true, 0, 1500);
	    //plotter.Draw(prefix + "htnojets"     + suffix, "p_{T}^{lep1} + p_{T}^{lep2} + MET", 20, 0, "GeV",  scale, true, 0, 1500);
	    //plotter.Draw(prefix + "metPfType1"     + suffix, sm,                                  10, 0, "GeV",  scale, true, 0,  400);
	    //plotter.Draw(prefix + "Counter"        + suffix, "m_{ll} (" + sll + ")",               1, 0, "GeV",  scale, false, 80, 100);
	    //plotter.Draw(prefix + "mt2ll"        + suffix, mt2 + "(" + sll + ")",               10, 0, "GeV",  scale, false, 0, 400);
	    //plotter.Draw(prefix + "dphilmet1"    + suffix, "#Delta#phi(lep1,E_{T}^{miss})",      10, 2, "rad",  scale, false);
	    //plotter.Draw(prefix + "dphilmet2"    + suffix, "#Delta#phi(lep2,E_{T}^{miss})",      10, 2, "rad",  scale, false);
	  }
	  /*plotter.Draw(prefix + "nbjet20cmvav2l" + suffix, "number of 20 GeV cmvav2l b-jets",   -1, 0, "NULL", scale);
	  plotter.Draw(prefix + "metPfType1Phi"  + suffix, sm + " #phi",                         5, 2, "rad",  scale);
	  plotter.Draw(prefix + "jet1eta"        + suffix, "leading jet #eta",                  -1, 1, "NULL", scale, false);
	  plotter.Draw(prefix + "jet2eta"        + suffix, "trailing jet #eta",                 -1, 1, "NULL", scale, false);
	  plotter.Draw(prefix + "jet1phi"        + suffix, "leading jet #phi",                   5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "jet2phi"        + suffix, "trailing jet #phi",                  5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "jet1pt"         + suffix, "leading jet p_{T}",                  5, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "jet2pt"         + suffix, "trailing jet p_{T}",                 5, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "pt2l"         + suffix, "p_{T}^{#font[12]{ll}}",             10, 0, "GeV",  scale, true, 0,  300);
	  plotter.Draw(prefix + "dphill"         + suffix, "#Delta#phi(lep1,lep2)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "detall"         + suffix, "#Delta#eta(lep1,lep2)",              5, 2, "rad",  scale, true, 0, 5);*/
																   
	  if (!allplots) continue;
 
	  plotter.Draw(prefix + "dyll"         + suffix, "lepton #Delta#eta",                 -1, 3, "NULL", scale);
	  plotter.Draw(prefix + "dphimetjet"   + suffix, "min #Delta#phi(jet," + sm + ")",     5, 2, "rad",  scale);
	  plotter.Draw(prefix + "dphimetbbll"  + suffix, "#Delta#phi(bbll," + sm + ")",        5, 2, "rad",  scale);
	  plotter.Draw(prefix + "mllbb"        + suffix, "m_{" + sll + "bb}",                 10, 0, "GeV",  scale, false, 0, 600);
	  plotter.Draw(prefix + "meff"         + suffix, "m_{eff}",                           10, 0, "GeV",  scale, false, 0, 600);
	  plotter.Draw(prefix + "ptbll"        + suffix, "p_{T}^{llmet}",                     10, 0, "GeV",  scale, false, 0, 600);
	  plotter.Draw(prefix + "mt2bb"        + suffix, mt2 + "(bb)" ,                       10, 0, "GeV",  scale, false, 0, 800);
	  plotter.Draw(prefix + "mt2lblb"      + suffix, mt2 + "(" + sl + "b" + sl + "b)",    10, 0, "GeV",  scale, false, 0, 800);
	  plotter.Draw(prefix + "MetMeff_Met"  + suffix, "MetMEff-Met",                        1, 0, "GeV",  scale, false);
	  plotter.Draw(prefix + "dphimetptbll" + suffix, "#Delta#phi(llmet," + sm + ")",       5, 2, "rad",  scale);
	  plotter.Draw(prefix + "dphijet1met"  + suffix, "#Delta#phi(jet1,E_{T}^{miss})",      5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphijet2met"  + suffix, "#Delta#phi(jet2,E_{T}^{miss})",      5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphijj"       + suffix, "#Delta#phi(jet1,jet2)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphijjmet"    + suffix, "#Delta#phi(jj,E_{T}^{miss})",        5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphilep1jet1" + suffix, "#Delta#phi(lep1,jet1)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphilep1jet2" + suffix, "#Delta#phi(lep1,jet2)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphilep2jet1" + suffix, "#Delta#phi(lep2,jet1)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphilep2jet2" + suffix, "#Delta#phi(lep2,jet2)",              5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "dphillstar"   + suffix, "#Delta#phi*(lep1,lep2)",             5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "metTtrkPhi"   + suffix, "track E_{T}^{miss} #phi",            5, 2, "rad",  scale, false);
	  plotter.Draw(prefix + "drll"         + suffix, "#DeltaR(lep1,lep2)",                 5, 1, "NULL", scale, false);
	  plotter.Draw(prefix + "jet1mass"     + suffix, "leading jet mass",                  -1, 0, "GeV",  scale, true, 0,   50);
	  plotter.Draw(prefix + "jet2mass"     + suffix, "trailing jet mass",                 -1, 0, "GeV",  scale, true, 0,   50);
	  plotter.Draw(prefix + "mc"           + suffix, "m_{c}",                             10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "metTtrk"      + suffix, "track E_{T}^{miss}",                10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "mpmet"        + suffix, "min projected E_{T}^{miss}",        10, 0, "GeV",  logY,  true, 0,  200);
	  plotter.Draw(prefix + "metPuppi"     + suffix, "PUPPI E_{T}^{miss}",                10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "mth"          + suffix, "m_{T}^{H}",                         10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "metmeff"      + suffix, "MET/Meff",                           5, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "mtw1"         + suffix, "m_{T}^{W,1}",                       10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "mtw2"         + suffix, "m_{T}^{W,2}",                       10, 0, "GeV",  scale, true, 0,  400);
	  plotter.Draw(prefix + "ptww"         + suffix, "p_{T}^{WW}",                        10, 0, "GeV",  scale, true, 0,  600);
	  plotter.Draw(prefix + "sumjpt12"     + suffix, "p_{T}^{jet1} + p_{T}^{jet2}",       10, 0, "GeV",  scale, true, 0,  600);
	  plotter.Draw(prefix + "sumpt12"      + suffix, "p_{T}^{lep1} + p_{T}^{lep2}",       10, 0, "GeV",  scale, true, 0,  600);
	  plotter.Draw(prefix + "MR"           + suffix, "MR",                                 5, 0, "GeV",  scale);
	  plotter.Draw(prefix + "R2"           + suffix, "R2",                                 5, 0, "NULL", scale, false);
	  plotter.Draw(prefix + "Rpt"          + suffix, "Rpt",                                5, 0, "NULL", scale, false);
	  plotter.Draw(prefix + "invGamma"     + suffix, "invGamma",                           5, 0, "NULL", scale, false);
	  plotter.Draw(prefix + "Mdr"          + suffix, "Mdr",                                5, 0, "GeV",  scale);
	  plotter.Draw(prefix + "DeltaPhiRll"  + suffix, "DeltaPhiRll",                        5, 2, "rad",  scale, false);
	  
	} 
      
    }
  
  // Copy index.php in every directory
  //----------------------------------------------------------------------------
  gSystem->Exec("for dir in $(find ./ -type d); do cp -n ../index.php $dir/; done");
  gSystem->Exec("rm -f index.php");
}


//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
# ifndef __CINT__
int main(int argc, char ** argv)
{
  if (argc < 2) {
    
    printf("\n rm -rf %s\n\n", outputdir.Data());
    
    for (int i=0; i<ncut; i++) printf(" ./runPlotter %s\n", scut[i].Data());
    
    printf("\n");
    
    exit(0);
  }
  
  if (argc == 2)
    runPlotter(argv[1]);
  else
    runPlotter(argv[1], argv[2]);
  
  return 0;
}
# endif
