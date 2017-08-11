#include "../test/HistogramReader.h"

// Constants
//------------------------------------------------------------------------------
const Bool_t datadriven = false;
const Bool_t allplots  = false;
const Bool_t dosystematics = true;

//const TString inputdir  = "../rootfiles/ZpeakJetCR/";
const TString inputdir  = "../minitrees/rootfiles3R/nominal/";
const TString outputdir = "figures/";

const TString sl  = "#font[12]{l}";
const TString sll = "#font[12]{ll}";
const TString sm  = "E_{T}^{miss}";

enum {linY, logY};


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
  plotter.SetSavePdf    (false);

  if (option.Contains("nostack"))
    {
      plotter.SetDrawRatio(false);
    }
  else
    {
      plotter.SetLuminosity(lumi);
      plotter.SetDrawRatio (true);
    }

  float SF_ttZ = 1., SF_ZMet = 1.;
  if (level.Contains("_SR")) { SF_ttZ = 1.56; SF_ZMet = 1.39; }

  // Get the data
  //----------------------------------------------------------------------------
  //plotter.AddData("01_Data", "data", color_Data);


  // Add processes
  //----------------------------------------------------------------------------
  //plotter.AddProcess("14_HZ",        "HZ",       color_HZ);
  plotter.AddProcess("13_VVV",      "VVV",      color_VVV);
  plotter.AddProcess("03_VZ",        "VZ (#rightarrow 2l)",       color_VZ,  roc_background, SF_ZMet);
  //plotter.AddProcess("11_Wg",        "W#gamma",  color_Wg);
  //plotter.AddProcess("15_WgStar",    "W#gamma*", color_WgStar);
  plotter.AddProcess("09_TTW",       "t#bar{t}W",      color_TTV);
  plotter.AddProcess("10_TTZ",       "t#bar{t}Z",      color_TTZ,  roc_background, SF_ttZ);
  plotter.AddProcess("11_HWW",       "HWW",      color_HWW);
  plotter.AddProcess("02_WZTo3LNu",  "WZ (#rightarrow 3l)",       color_WZTo3LNu,  roc_background);
  plotter.AddProcess("06_WW",        "WW",       color_WW);
  plotter.AddProcess("05_ST",        "tW",       color_ST);
  plotter.AddProcess("07_ZJetsHT",     "Z+jets",   color_ZJets,  roc_background, SF_ZMet);
  if (inputdir.Contains("SS")) plotter.AddProcess("TTToSemiLepton", "t#bar{t} Semilep.",  41);
  if (!inputdir.Contains("SS")) plotter.AddProcess("04_TTTo2L2Nu", "t#bar{t}",       color_TTTo2L2Nu);
  else plotter.AddProcess("04_TTTo2L2Nu", "t#bar{t}",       color_TTTo2L2Nu);
  //else plotter.AddProcess("TTJets", "#bar{t}t",       color_TTTo2L2Nu);
  //if (inputdir.Contains("SS")) plotter.AddProcess("WJetsToLNu", "WJets",      color_WJets);

  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm175", "#tilde{t}#rightarrow t#tilde{#chi}^{0}_{1} (m_{#tilde{t}}=350, m_{#tilde{#chi}^{0}_{1}}=175)",2);
  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm225", "#tilde{t}#rightarrow t#tilde{#chi}^{0}_{1} (m_{#tilde{t}}=350, m_{#tilde{#chi}^{0}_{1}}=225)",3);
  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm263", "#tilde{t}#rightarrow t#tilde{#chi}^{0}_{1} (m_{#tilde{t}}=350, m_{#tilde{#chi}^{0}_{1}}=263)",6);
  //plotter.AddSignal("T2tt_mStop-150to250_Sm250_Xm125", "#tilde{t}#rightarrow t#tilde{#chi}^{0}_{1} (m_{#tilde{t}}=250, m_{#tilde{#chi}^{0}_{1}}=125)",4);
  //plotter.AddSignal("T2tt_mStop-400to1200_Sm450_Xm325","#tilde{t}#rightarrow t#tilde{#chi}^{0}_{1} (m_{#tilde{t}}=450, m_{#tilde{#chi}^{0}_{1}}=325)",7);
  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm225", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,225)",kViolet);
  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm175", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,175)",kRed);
  //plotter.AddSignal("T2tt_mStop-350to400_Sm350_Xm263", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (350,263)",kOrange);
  //plotter.AddSignal("T2tt_mStop-150to250_Sm250_Xm125", "#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (250,125)",kBlue);
  //plotter.AddSignal("T2tt_mStop-400to1200_Sm450_Xm325","#tilde{t} #tilde{t}, #tilde{t} #rightarrow t#tilde{#chi}^{0}_{1} (450,325)",kCyan);

  //plotter.AddSignal("TChiSlep_Xm500_Xm400",  "#tilde{#chi}^{#pm}#rightarrow #tilde{l}#tilde{#nu} (500,400)",  2);
  //plotter.AddSignal("TChiSlep_Xm300_Xm200",  "#tilde{#chi}^{#pm}#rightarrow #tilde{l}#tilde{#nu} (300,200)",  3);
  //plotter.AddSignal("TChiSlep_Xm100_Xm1",  "#tilde{#chi}^{#pm}#rightarrow #tilde{l}#tilde{#nu} (100,1)",  4);
  plotter.AddSignal("TChiSlep_Xm200_Xm1",   "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{l}#tilde{#nu} (200,  1)", kRed);
  plotter.AddSignal("TChiSlep_Xm500_Xm200", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{l}#tilde{#nu} (500,200)", kViolet);
  plotter.AddSignal("TChiSlep_Xm800_Xm400", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{l}#tilde{#nu} (800,400)", kOrange);
  plotter.AddSignal("TChiWW_Xm200_Xm25",  "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (200, 25)", kBlue);
  plotter.AddSignal("TChiWW_Xm300_Xm100", "#tilde{#chi}^{#pm}#tilde{#chi}^{#pm}, #tilde{#chi}^{#pm}#rightarrow#tilde{#chi}^{0}_{1}W (300,100)", kCyan);

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

  if (dosystematics) {
    plotter.AddSystematic("Stop", "Luminosity");
    plotter.AddSystematic("Stop", "Trigger");
    plotter.AddSystematic("Stop", "MT2llTop");
    plotter.AddSystematic("Stop", "MT2llWW");
    plotter.AddSystematic("Stop", "ttZSF");
    plotter.AddSystematic("Stop", "ZMetSF");
    plotter.AddSystematic("Stop", "Fake");
    plotter.AddSystematic("Stop", "Idiso");
    plotter.AddSystematic("Stop", "JES");
    plotter.AddSystematic("Stop", "MET");
    plotter.AddSystematic("Stop", "PDF");
    plotter.AddSystematic("Stop", "Q2");
    plotter.AddSystematic("Stop", "Reco");
    plotter.AddSystematic("Stop", "Toppt");
    plotter.AddSystematic("Stop", "Metfastsim");
    plotter.AddSystematic("Stop", "Pileup");
    plotter.AddSystematic("Stop", "Fastsim");
    plotter.AddSystematic("Stop", "BtagFS");
    plotter.AddSystematic("Stop", "Btag");
  }
  
  // Draw distributions
  //----------------------------------------------------------------------------
  if (!option.Contains("nostack")) plotter.SetDrawYield(true);

  float m2l_xmin   = (level.Contains("WZ")) ?  60 :   0;  // [GeV]
  float m2l_xmax   = (level.Contains("WZ")) ? 120 : 300;  // [GeV]
  int   m2l_ngroup = (level.Contains("WZ")) ?   2 :   5;
  
  for (int j=0; j<=njetbin; j++)
    {
      if (j != njetbin) continue;   
         
      TString jetbin = (j < njetbin) ? Form("/%djet", j) : "";

      gSystem->mkdir(outputdir + level + jetbin, kTRUE);

      TString prefix = level + jetbin + "/h_";

      for (int i=firstchannel; i<=lastchannel; i++)
	{
	  TString suffix = "_" + schannel[i];
	  
	  TString title = (i < lastchannel) ? lchannel[i] : "inclusive";

	  plotter.SetTitle(title);
	  
	  // Common histograms
	  //--------------------------------------------------------------------
	  /*plotter.Draw(prefix + "m2l" + suffix, "m_{" + sll + "}", m2l_ngroup, 0, "GeV", linY, true, m2l_xmin, m2l_xmax);
	  plotter.Draw(prefix + "lep1pt"         + suffix, "leading lepton p_{T}",              10, 0, "GeV",  scale, true, 0,  250);
	  plotter.Draw(prefix + "lep2pt"         + suffix, "trailing lepton p_{T}",              5, 0, "GeV",  scale, true, 0,  150);
	  plotter.Draw(prefix + "lep1eta"        + suffix, "leading lepton #eta",               -1, 1, "NULL", scale);
	  plotter.Draw(prefix + "lep2eta"        + suffix, "trailing lepton #eta",              -1, 1, "NULL", scale);
	  plotter.Draw(prefix + "lep1phi"        + suffix, "leading lepton #phi",                5, 2, "rad",  scale);
	  plotter.Draw(prefix + "lep2phi"        + suffix, "trailing lepton #phi",               5, 2, "rad",  scale);
	  plotter.Draw(prefix + "nvtx"           + suffix, "number of vertices",                -1, 0, "NULL", scale, true, 0,   30);*/
	  plotter.Draw(prefix + "MT2ll"        + suffix, "M_{T2}(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	  //if (level.Contains("_SR3"))  plotter.Draw(prefix + "MT2llisr"        + suffix, "M_{T2}(" + sll + ")",               1, 0, "GeV",  scale, false, 0, 140);
	  //plotter.Draw(prefix + "Counter"        + suffix, "m_{ll} (" + sll + ")",               1, 0, "GeV",  linY, false, 80, 100);
	  //plotter.Draw(prefix + "MT2_Met"      + suffix, "M_{T2}-Met",                         1, 0, "GeV",  scale, false);
	  //plotter.Draw(prefix + "MET"     + suffix, sm,                                  1, 0, "GeV",  scale, false, 0,  140);
	  //plotter.Draw(prefix + "lep1pt"      + suffix, "lep1pt",                         1, 0, "GeV",  scale, true, 0, 200);
	  //plotter.Draw(prefix + "maxjetpt"      + suffix, "leading jet pt",                  2, 0, "GeV",  scale, false, 0, 600);
	  //plotter.Draw(prefix + "dphiminlepmet"    + suffix, "#Delta#phi(lep,E_{T}^{miss})",      10, 2, "rad",  linY, false);
	  //plotter.Draw(prefix + "njet20"           + suffix, "number of 20 GeV jets",             -1, 0, "NULL", scale);
	  //plotter.Draw(prefix + "njet30"           + suffix, "number of 30 GeV jets",             -1, 0, "NULL", scale);
	  //plotter.Draw(prefix + "njet20"           + suffix, "number of 20 GeV jets",             -1, 0, "NULL", linY);
	  //plotter.Draw(prefix + "njet30"           + suffix, "number of 30 GeV jets",             -1, 0, "NULL", linY);
	  if (inputdir.Contains("rootfiles/nominal")) {
	    //plotter.Draw(prefix + "jet1pt"      + suffix, "jep1pt",                         10, 0, "GeV",  scale, true, 0, 200);
	    //plotter.Draw(prefix + "m2l" + suffix, "m_{" + sll + "}", m2l_ngroup, 0, "GeV", logY, true, m2l_xmin, m2l_xmax);
	    //plotter.Draw(prefix + "dphillmet"      + suffix, "#Delta#phi(" +sll + "," + sm + ")",  10, 2, "rad",  scale); 
	    //plotter.Draw(prefix + "njet"           + suffix, "number of 20 GeV jets",             -1, 0, "NULL", scale);
	    //plotter.Draw(prefix + "metPfType1"     + suffix, sm,                                  10, 0, "GeV",  scale, true, 0,  140);
	    //plotter.Draw(prefix + "Counter"        + suffix, "m_{ll} (" + sll + ")",               1, 0, "GeV",  scale, false, 80, 100);
	    //plotter.Draw(prefix + "mt2ll"        + suffix, "M_{T2}(" + sll + ")",               10, 0, "GeV",  scale, false, 0, 400);
	    //plotter.Draw(prefix + "dphilmet1"    + suffix, "#Delta#phi(lep1,E_{T}^{miss})",      10, 2, "rad",  scale, false);
	    //plotter.Draw(prefix + "dphilmet2"    + suffix, "#Delta#phi(lep2,E_{T}^{miss})",      10, 2, "rad",  scale, false);
	  }
	  /*plotter.Draw(prefix + "nbjet20cmvav2l" + suffix, "number of 20 GeV cmvav2l b-jets",   -1, 0, "NULL", scale);
	  plotter.Draw(prefix + "nbjet30csvv2m"  + suffix, "number of 30 GeV csvv2m b-jets",    -1, 0, "NULL", scale);
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
	  plotter.Draw(prefix + "mt2bb"        + suffix, "M_{T2}(bb)" ,                       10, 0, "GeV",  scale, false, 0, 800);
	  plotter.Draw(prefix + "mt2lblb"      + suffix, "M_{T2}(" + sl + "b" + sl + "b)",    10, 0, "GeV",  scale, false, 0, 800);
	  plotter.Draw(prefix + "MetMeff_Met"  + suffix, "MetMEff-Met",                        1, 0, "GeV",  scale, false);
	  plotter.Draw(prefix + "HTvisible_Met"+ suffix, "HTvisible-Met",                      1, 0, "GeV",  scale, false);
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
	  plotter.Draw(prefix + "ht"           + suffix, "H_{T}",                             20, 0, "GeV",  scale, true, 0, 1500);
	  plotter.Draw(prefix + "htvisible"    + suffix, "#sum_{jet,lepton} p_{T}",           20, 0, "GeV",  scale, true, 0, 1500);
	  plotter.Draw(prefix + "htjets"       + suffix, "#sum_{jet} p_{T}",                  20, 0, "GeV",  scale, true, 0, 1500);
	  plotter.Draw(prefix + "htnojets"     + suffix, "p_{T}^{lep1} + p_{T}^{lep2} + MET", 20, 0, "GeV",  scale, true, 0, 1500);
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
