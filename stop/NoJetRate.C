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

bool DoBackgroungProcesses = false;
bool DoSignalMassPoints = true;

TString RootFilesDirectory = "../minitrees/rootfiles3R/";
TString MassPointFileName = "../../PlotsConfigurations/Configurations/T2tt/MassPointList_TChiWW.txt";

int const nObservables = 4;
TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen"};

int const nChannels = 4;
TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};

void NoJetRateByProcess(TString ProcessName, TString Variation) {
  
  gSystem->mkdir(RootFilesDirectory + "Nojetrateup/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "Nojetratedo/Stop/", kTRUE);

  TFile *InputRootFile = TFile::Open(RootFilesDirectory + "nominal/Stop/" + ProcessName + ".root");
      
  int VariationSign = (Variation=="up") ?    1 :   -1;
    
  TFile *OutputFile = new TFile(RootFilesDirectory + "Nojetrate" + Variation + "/Stop/" + ProcessName + ".root", "recreate");
  
  for (int sr = 0; sr<ncut; sr++) {
    
    if (!scut[sr].Contains("_SR")) continue;
    if (!scut[sr].Contains("_NoJet")) continue;

    float RateUncertainty = (scut[sr].Contains("_SR1")) ? 0.20 : 0.35;
    
    TString NoTagRegion = scut[sr]; NoTagRegion.ReplaceAll("_NoJet", "_NoTag");
    
    OutputFile->cd();
    
    gDirectory->mkdir(scut[sr]);
    gDirectory->mkdir(NoTagRegion);
    
    for (int ob = 0; ob<nObservables; ob++) 
      for (int ch = 0; ch<nChannels; ch++) {
	
	TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];
	
	OutputFile->cd(scut[sr]);
	
	TH1F *OutputHistoNoJet = (TH1F*) InputRootFile->Get(scut[sr] + "/" + HistoName);
	
	float NominalNoJetBinContent[100];
	for (int ib = 1; ib<=OutputHistoNoJet->GetNbinsX(); ib++)
	  NominalNoJetBinContent[ib] = OutputHistoNoJet->GetBinContent(ib);
	
	OutputHistoNoJet->Scale(1.+VariationSign*RateUncertainty);
	
	OutputHistoNoJet->Write();

	float NoJetToNoTagMigration[100];
	for (int ib = 1; ib<=OutputHistoNoJet->GetNbinsX(); ib++) {
	  float ScaledNoJetBinContent = OutputHistoNoJet->GetBinContent(ib);
	  NoJetToNoTagMigration[ib] = -1.*(ScaledNoJetBinContent - NominalNoJetBinContent[ib]);
	}

	OutputFile->cd(NoTagRegion);
	
	TH1F *OutputHistoNoTag = (TH1F*) InputRootFile->Get(NoTagRegion + "/" + HistoName);
	
	for (int ib = 1; ib<=OutputHistoNoJet->GetNbinsX(); ib++) {
	  float NominalNoTagBinContent = OutputHistoNoTag->GetBinContent(ib);
	  OutputHistoNoTag->SetBinContent(ib, NominalNoTagBinContent+NoJetToNoTagMigration[ib]);
	}
	
	OutputHistoNoTag->Write();
	
      }
   
  }
 
  OutputFile->Close();
  
  InputRootFile->Close();
  
}

void NoJetRate() {
  
  int const nBackgroundProcesses = 11;
  TString BackgroundProcessName[nBackgroundProcesses] = {"02_WZTo3LNu", "03_VZ",  "04_TTTo2L2Nu",
							 "05_ST",       "06_WW",  "07_ZJets",
							 "07_ZJetsHT",  "09_TTW", "10_TTZ",
							 "11_HWW",      "13_VVV"};
  
  if (DoBackgroungProcesses) {

    for (int pr = 0; pr<nBackgroundProcesses; pr++) {

      NoJetRateByProcess(BackgroundProcessName[pr], "up"); 
      NoJetRateByProcess(BackgroundProcessName[pr], "do"); 

    }

  }

  if (DoSignalMassPoints) {

    ifstream MassPointFile; MassPointFile.open(MassPointFileName);

    while (MassPointFile) {

      TString MassPoint;
      MassPointFile >> MassPoint;

      if (MassPoint!="" && !MassPoint.Contains("#")) {

	NoJetRateByProcess("T" + MassPoint, "up"); 
	NoJetRateByProcess("T" + MassPoint, "do"); 
    
      }

    }

  }

}
