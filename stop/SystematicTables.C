#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TH1.h"
#include "THStack.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <fstream>
#include <iostream>

#include "../include/CutsStop.h"

int const nSRs = 12;
bool UseSR[nSRs]; 

void SetUseSR(TString Analysis) {

  if (Analysis=="_Stop") {
    UseSR[0] = 0; UseSR[1] = 0; UseSR [2] = 1; UseSR [3] = 1;
    UseSR[4] = 0; UseSR[5] = 0; UseSR [6] = 1; UseSR [7] = 1;
    UseSR[8] = 0; UseSR[9] = 0; UseSR[10] = 1; UseSR[11] = 1; 
  } else if (Analysis=="_Chargino") {
    UseSR[0] = 1; UseSR[1] = 1; UseSR [2] = 0; UseSR [3] = 0;
    UseSR[4] = 1; UseSR[5] = 1; UseSR [6] = 0; UseSR [7] = 0;
    UseSR[8] = 1; UseSR[9] = 0; UseSR[10] = 0; UseSR[11] = 0; 
  }

}

TString ReadSystematic(TString SystName, TString SRname, TString Channel, TString Type, TString Aname) {

  TString isrstring = SRname.Contains("_isr") ? "isr" : "";
  SRname.ReplaceAll("_isr", "");

  TString TableName;
  if       (Aname=="_Chargino" && Type!="_Background")    TableName = "./Tables/Syst_02_" + SRname + "_MT2ll" + isrstring + Channel + "_TChiSlep_Xm500_Xm200";
  else if  (Aname=="_Stop" && Type!="_Background")  { TableName = "./Tables/Syst_02_" + SRname + "_MT2ll" + isrstring + Channel + "_T2tt_mStop-350to400_Sm350_Xm225";}
  else  TableName = "./TablesSystPaper/Syst_02_" + SRname + "_MT2ll" + isrstring + Channel + Type;
  //TString TableName = "./TablesSystPaper/Syst_02_" + SRname + "_MT2ll" + isrstring + Channel + Type;
  /*if (SRname.Contains("Tag") || (SRname.Contains("Veto") && (isrstring!="" || !SRname.Contains("SR3")))) 
    TableName.ReplaceAll("Syst/", "T2ttNoTopPt/");
  else
  TableName.ReplaceAll("Syst/", "TChiNoTopPt/");*/
  std::ifstream outFile(TableName + ".txt",std::ios::out);

  float maxsyst = -1., minsyst = 999.;
  float imaxsyst = -1., iminsyst = 999.;

  float bin[7], ibin[7];
  
  int firstbin = 4;//(Type.Contains("TChi")) ? 4 : 0;

  while (outFile) {

    TString systname; 
    outFile >> systname >> bin[0] >> bin[1] >> bin[2] >> bin[3] >> bin[4] >> bin[5] >> bin[6];

    if (systname!=SystName) continue;

    for (int ib = firstbin; ib<7; ib++) 
      if (bin[ib]>=1.0) ibin[ib] = ceil(bin[ib]-0.5);
      else ibin[ib] = 0;

    for (int ib = firstbin; ib<7; ib++) {
      if (ibin[ib]>imaxsyst)
	imaxsyst = ibin[ib];
      if (ibin[ib]<iminsyst)
	iminsyst = ibin[ib];
    }
    for (int ib = firstbin; ib<7; ib++) {
      if (bin[ib]>maxsyst)
	maxsyst = bin[ib];
      if (bin[ib]<minsyst)
	minsyst = bin[ib];
    }
    
  }
  
  outFile.close();
  std::cout<<maxsyst<<endl; 
  if (maxsyst==-1.) return "missing";

  TString thissyst = "$";
  cout << SystName << " " << minsyst << " " << maxsyst << " " << iminsyst << " " << imaxsyst << " " << endl;//fabs(maxsyst-minsyst)/minsyst<< endl;
  if ((imaxsyst>=1. && iminsyst<1.) || (iminsyst>0. && fabs(imaxsyst-iminsyst)/iminsyst>=0.125)) {
    thissyst += int(iminsyst); thissyst += "$-$"; thissyst += int(imaxsyst); thissyst += "\\%$"; 
  } /*else if (iminsyst<1. && imaxsyst<1.) {
    int x1 = ceil(10.*minsyst);
    int x2 = ceil(10.*maxsyst);
    if (x2-x1>=2) {
      thissyst += "0."; thissyst += x1; thissyst += "$-$"; 
      if (x2<10) { thissyst += "0."; thissyst += x2;  }
      else thissyst += "1";
      thissyst += "\\%$";
    } else {
      thissyst += "0."; thissyst += x2; thissyst += "\\%$";
    }
  } else {
    if (maxsyst<1.) {
      thissyst += thissyst += "0."; thissyst += int(10.*maxsyst); thissyst += "\\%$";
    } else if (maxsyst<10.) {
      thissyst +=  maxsyst; thissyst.Remove(4); thissyst += "\\%$";
    } else {
      thissyst +=  int(maxsyst); thissyst += "\\%$"; 
    }
    }*/
  else if (imaxsyst<1.) {
    thissyst += "<1\\%$";
  } else {
    thissyst +=  maxsyst; thissyst.Remove(4); thissyst += "\\%$";
  }

  return thissyst;
    
}

void SystematicTables(TString Type, TString Analysis, TString Channel) {

  TString SRname[nSRs] = {"SR1_NoJet", "SR1_NoTag", "SR1_Veto", "SR1_Tag", 
			  "SR2_NoJet", "SR2_NoTag", "SR2_Veto", "SR2_Tag", 
			  "SR3_Veto",  "SR3_Tag",   "SR3_Veto_isr", "SR3_Tag_isr"};
  TString SRtitle[nSRs] = {"SR1^{\\text{0Jet}}_{\\text{0Tag}}", "SR1^{\\text{Jets}}_{\\text{0Tag}}", "SR1_{\\text{0Tag}}", "SR1_{\\text{Tags}}", 
			   "SR2^{\\text{0Jet}}_{\\text{0Tag}}", "SR2^{\\text{Jets}}_{\\text{0Tag}}", "SR2_{\\text{0Tag}}", "SR2_{\\text{Tags}}", 
			   "SR3_{\\text{0Tag}}", "SR3_{\\text{Tags}}", "SR3^{\\text{ISR}}_{\\text{0Tag}}", "SR3^{\\text{ISR}}_{\\text{Tags}}"};

  int const nSyst = 26;
  TString SystName[nSyst] = {"Statistics", "Luminosity", "Trigger", "MT2llTop", "MT2llWW",
			     "Fake", "Idiso", "JES", "MET", "PDF", "Q2", "Reco", "Toppt",
			     "Isrnjet", "Metfastsim", "Pileup", "Fastsim", "BtagFS", 
			     "Btag", "Btaglight", "ttZSF", "WZSF", "ZZSF", "ZZshape",
			     "DYshape", "normDY"};
  int IsSignal[nSyst]     = {1, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0,
			     0, 1, 0, 0, 0,  
			     0, 0, 0, 0, 0, 0, 
			     0, 0};
//  int IsSignal[nSyst]     = {0, 1, 1, 0, 0,
//			     0, 1, 1, 1, 0, 1, 1, 0,
//			     1, 1, 1, 1, 1,  
//			     1, 1, 0, 0, 0, 0, 
//			     0, 0};
  int IsBackground[nSyst] = {1, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0,
			     0, 1, 0, 0, 0,  
			     0, 0, 0, 0, 0, 0, 
			     0, 0};
//  int IsBackground[nSyst] = {0, 1, 1, 1, 1,
//			     1, 1, 1, 1, 1, 1, 1, 1,
//			     0, 0, 0, 0, 0,  
//			     1, 1, 1, 1, 1, 1, 
//			     1, 1};
  TString SystDescr[nSyst] = {"MC statistics", "Luminosity", "Trigger", "\\mtll shape (Top)", "\\mtll shape (WW)",
			      "Nonprompt leptons", "Lepton ID and isolation", "JES", "Unclustered energy", "PDFs", "QCD scale", "Track reconsruction", "\\ttbar \\pt reweighting",
			      "ISR reweighting", "\\ptmiss FastSim", "Pileup", "Leptons Fastsim", "b-tagging FastSim", 
			      "b-tagging", "b-tagging (light jets)", "ttZ normalization", "WZ normalization", "ZZ normalization", "ZZ k-factors",
			      "\\mtll shape (Drell-Yan)", "Drell-Yan normalization"};
  
  TString AnalysisName[2] = {"_Stop", "_Chargino"};

  
  TString TableName = "./TablesSyst/Systematic" + Type + Analysis + Channel;
  std::ofstream inFile(TableName + ".tex",std::ios::out);
  
  //inFile << "\\begin{table}[htb]" << endl;
  //inFile << "\\tiny" << endl;
  inFile << "\\begin{center}" << endl;
  if (Analysis=="_Stop")
    inFile << "\\begin{tabular}{|l|cccccc|}" << endl;
  else if (Analysis=="_Chargino") 
    inFile << "\\begin{tabular}{|l|ccccc|}" << endl;
  else if (Analysis=="_All") 
    inFile << "\\begin{tabular}{|l|cccccc|ccccc|}" << endl;
  inFile << "\\hline\\hline" << endl;
  inFile << " Systematic ";
  for (int na = 0; na<2; na++) {
    if (Analysis=="_All" || Analysis==AnalysisName[na]) {
      SetUseSR(AnalysisName[na]);
      for (int sr = 0; sr<nSRs; sr++) 
	if (UseSR[sr]) 
      inFile << " & " << SRtitle[sr];
    }
  }
  inFile << " \\\\ " << endl;
  inFile << "\\hline" << endl;

  for (int is = 0; is<nSyst; is++) {
    if ( (Type=="_Background" && IsBackground[is]==1) ||
	 (Type!="_Background" && IsSignal[is]==1) ) {

      //if (SystName[is]=="DYnojet" && !SRname[sr].Contains("NoJet")) continue;
      //if (SystName[is]=="normDY" && SRname[sr].Contains("NoJet")) continue;
      if (SystName[is]=="DYnojet" && Analysis=="_Stop") continue;

      inFile << SystDescr[is];

      for (int na = 0; na<2; na++) {
	if (Analysis=="_All" || Analysis==AnalysisName[na]) {
	  SetUseSR(AnalysisName[na]);
	  for (int sr = 0; sr<nSRs; sr++) 
	    if (UseSR[sr]) {
	
	      TString systvalue; TString ThisSyst = SystName[is];
	      if (SystName[is]=="normDY" && SRname[sr].Contains("NoJet")) ThisSyst = "DYnojet";
	      systvalue = ReadSystematic(ThisSyst, SRname[sr], Channel, Type, Analysis);
	      inFile << " & " << systvalue;
	      
	    }
	}
      }


      inFile << " \\\\ " << endl;

    }
  }

  inFile << " \\hline\\hline" << endl;
  inFile << "\\end{tabular}" << endl;
  inFile << "\\end{center}" << endl;
  //inFile << "\\end{table}" << endl;
     
  inFile.close();

}
