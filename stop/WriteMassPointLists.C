#include <map>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
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
#include "TLorentzVector.h"

#include <fstream>
#include <iostream>

TString Model = "T2bW";
TString MinitreesDirectory = "/afs/cern.ch/work/s/scodella/Stop/CodeDevelopment/CMSSW_8_0_26_patch1/src/AnalysisCMS/minitrees/nominal/Stop/";
TString SystematicsDirectory = "/eos/cms/store/caf/user/scodella/BTV/Feb2017_Summer16_stop_ghent_others_isr/MCl2stop__SFWeights__";
const int nSystematicsToList = 2;
TString SystematicToList[nSystematicsToList] = {"JES", "MET"}; 
TString Variation[2] = {"up", "do"};

void WriteMassPointLists() {

  // Open output files
  ofstream MassPointList; MassPointList.open("./MassPointLists/MassPointList_" + Model + ".txt");
  ofstream MinitreesList; MinitreesList.open("./MassPointLists/samples_" + Model + "_Ghent_Summer16_minitrees.txt");
  ofstream SystematicsList[2*nSystematicsToList];
  for (int is = 0; is<nSystematicsToList; is++) 
    for (int vr = 0; vr<2; vr++) 
      SystematicsList[vr+2*is].open("./MassPointLists/samples_" + Model + "_Ghent_Summer16_" + SystematicToList[is] + Variation[vr] + ".txt"); 
  
  // Get mass plan histograms
  TFile *HistoFile = TFile::Open("./MassPlanHistograms/MassPlan_" + Model + ".root");
  TH2D *MassPlan = (TH2D*) HistoFile->Get("MassPlan");
  cout << MassPlan->GetEntries() << endl;

  // Loop over mass plan and fill lists
  for (int iXMass = 0; iXMass<2000; iXMass++) {

    bool NewXMass = true;

    for (int iNeutralinoMass = 0; iNeutralinoMass<2000; iNeutralinoMass++) {
       
	float XMass = 1.*iXMass;
	float NeutralinoMass  = 1.*iNeutralinoMass;
	
	int nEntries = MassPlan->GetBinContent(MassPlan->FindBin(XMass, NeutralinoMass));
	
	if (nEntries>0) {

	  if (NewXMass) {
	    
	    MassPointList << "#" << endl;
	    MinitreesList << "#" << endl;
	    for (int is = 0; is<nSystematicsToList; is++) 
	      for (int vr = 0; vr<2; vr++) 
		SystematicsList[vr+2*is] << "#" << endl;

	    NewXMass = false;

	  }
	  
	  TString ModelString = Model; ModelString.ReplaceAll("T", "");
	  TString TreeName = Model;
	  if (Model=="T2tt") {

	    TString AddString = "_mStop-";
	    if (iXMass<=250) AddString += "150to250";
	    else if (iXMass<350) AddString += "250to350";
	    else if (iXMass<400) AddString += "350to400";
	    else AddString += "400to1200";
	    ModelString += AddString;
	    TreeName += AddString;

	  }  

	  TreeName += ".root";

	  if (Model.Contains("T2")) ModelString += "_Sm";
	  else ModelString += "_Xm";

	  ModelString += iXMass;
	  ModelString += "_Xm";
	  ModelString += iNeutralinoMass;

	  MassPointList << ModelString << endl;

	  MinitreesList << MinitreesDirectory << TreeName <<  " 3 " << iXMass << " " << iNeutralinoMass << endl;
	  
	  for (int is = 0; is<nSystematicsToList; is++) 
	    for (int vr = 0; vr<2; vr++) 
	      SystematicsList[vr+2*is] << SystematicsDirectory << SystematicToList[is] << Variation[vr] << "__hadd/latino_" << TreeName << " " << iXMass << " " << iNeutralinoMass << endl;
	  
	}
	
    }
    
  }
	    
  MassPointList << "#" << endl;
  MinitreesList << "#" << endl;
  for (int is = 0; is<nSystematicsToList; is++) 
    for (int vr = 0; vr<2; vr++) 
      SystematicsList[vr+2*is] << "#" << endl;
  
}
