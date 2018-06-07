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

TFile *myfile0;

TH1F *GetHisto(TString HistoName) {
  
  TH1F *H1;
  if (HistoName.Contains("Veto") && HistoName.Contains("_sf")) {
    HistoName.ReplaceAll("Veto", "NoJet");
    HistoName.ReplaceAll("_sf",  "_ee");
    H1 = (TH1F*) myfile0->Get(HistoName);
    H1->SetDirectory(0);
    HistoName.ReplaceAll("NoJet", "NoTag");
    TH1F *H2 = (TH1F*) myfile0->Get(HistoName);
    H2->SetDirectory(0);
    H1->Add(H2);
    HistoName.ReplaceAll("_ee", "_mm");
    TH1F *H3 = (TH1F*) myfile0->Get(HistoName);
    H3->SetDirectory(0);
    H1->Add(H3);
    HistoName.ReplaceAll("NoTag", "NoJet");
    TH1F *H4 = (TH1F*) myfile0->Get(HistoName);
    H4->SetDirectory(0);
    H1->Add(H4);
  } else if (HistoName.Contains("Veto") && !HistoName.Contains("_sf")) {
    HistoName.ReplaceAll("Veto", "NoJet");
    H1 = (TH1F*) myfile0->Get(HistoName);
    H1->SetDirectory(0);
    HistoName.ReplaceAll("NoJet", "NoTag");
    TH1F *H2 = (TH1F*) myfile0->Get(HistoName);
    H2->SetDirectory(0);
    H1->Add(H2);
  } else if (!HistoName.Contains("Veto") && HistoName.Contains("_sf")) {
    HistoName.ReplaceAll("_sf", "_ee");
    H1 = (TH1F*) myfile0->Get(HistoName);
    H1->SetDirectory(0);
    HistoName.ReplaceAll("_ee", "_mm");
    TH1F *H2 = (TH1F*) myfile0->Get(HistoName);
    H2->SetDirectory(0);
    H1->Add(H2);
  } else {
    H1 = (TH1F*) myfile0->Get(HistoName);
    H1->SetDirectory(0);
  }

  return H1;

}

void EfficiencyMaps(TString Signal, TString Analysis, TString Channel, bool entries = false) {

  int const nSRs = 12;
  TString SRname[nSRs] = {"SR1_NoJet", "SR1_NoTag", "SR1_Veto", "SR1_Tag", 
			  "SR2_NoJet", "SR2_NoTag", "SR2_Veto", "SR2_Tag", 
			  "SR3_Veto",  "SR3_Tag",   "SR3_Veto_isr", "SR3_Tag_isr"};
  TString SRtitle[nSRs] = {"SR1^{\\text{0Jet}}_{\\text{0Tag}}", "SR1_^{\\text{Jets}}_{\\text{0Tag}}", "SR1_{\\text{0Tag}}", "SR1_{\\text{Tags}}", 
			   "SR2^{\\text{0Jet}}_{\\text{0Tag}}", "SR2_^{\\text{Jets}}_{\\text{0Tag}}", "SR2_{\\text{0Tag}}", "SR2_{\\text{Tags}}", 
			   "SR3_{\\text{0Tag}}", "SR3_{\\text{Tags}}", "SR3^{\\text{ISR}}_{\\text{0Tag}}", "SR3^{\\text{ISR}}_{\\text{Tags}}"};
  bool UseSR[nSRs]; 
  if (Analysis=="_Stop") {
    UseSR[0] = 0; UseSR[1] = 0; UseSR [2] = 1; UseSR [3] = 1;
    UseSR[4] = 0; UseSR[5] = 0; UseSR [6] = 1; UseSR [7] = 1;
    UseSR[8] = 0; UseSR[9] = 0; UseSR[10] = 1; UseSR[11] = 1; 
  } else if (Analysis=="_Chargino") {
    UseSR[0] = 1; UseSR[1] = 1; UseSR [2] = 0; UseSR [3] = 0;
    UseSR[4] = 1; UseSR[5] = 1; UseSR [6] = 0; UseSR [7] = 0;
    UseSR[8] = 1; UseSR[9] = 0; UseSR[10] = 0; UseSR[11] = 0; 
  }

  int nSRsignal = 0;
  for (int sr = 0; sr<nSRs; sr++)
    if (UseSR[sr]) 
      nSRsignal++;

  TH2F *EffMap  = new TH2F("EffMap",  "", 7, 0., 140., nSRsignal, 0., nSRsignal);
  TH2F *StatMap = new TH2F("StatMap", "", 7, 0., 140., nSRsignal, 0., nSRsignal);

  EffMap->SetXTitle("M_{T2}(ll) [GeV]");
  StatMap->SetXTitle("M_{T2}(ll) [GeV]");
  
  TString Syst = entries ? "noweight" : "nominal"; 
  myfile0 = new TFile("../minitrees/rootfilesOct17/" + Syst + "/Stop/" + Signal + ".root", "read");

  int isr = 1;

  for (int sr = 0; sr<nSRs; sr++)
    if (UseSR[sr]) {

      EffMap->GetYaxis()->SetBinLabel(isr, SRtitle[sr]);
      StatMap->GetYaxis()->SetBinLabel(isr, SRtitle[sr]);

      TString HistoName = "Stop/02_" + SRname[sr] + "/h_MT2ll" + Channel; 
      TString HistoNameGen = HistoName; 
      HistoNameGen.ReplaceAll("SR1", "SR1gen");
      HistoNameGen.ReplaceAll("SR2", "SR2gen");
      HistoNameGen.ReplaceAll("SR3", "SR3gen");
      HistoNameGen.ReplaceAll("MT2ll", "MT2llgen");
      if (SRname[sr].Contains("_isr")) {
	HistoName.ReplaceAll("_isr", "");
	HistoName.ReplaceAll("MT2ll", "MT2llisr");
	HistoNameGen.ReplaceAll("_isr", "");
	HistoNameGen.ReplaceAll("MT2ll", "MT2llisr");
      }

      TH1F *ThisHisto = GetHisto(HistoName);
      TH1F *ThisHistoGen = GetHisto(HistoNameGen);

      for (int ib = 1; ib<=7; ib++) {
	
	//float XXX = ThisHisto->GetBinContent(ib);
	//float EEE = ThisHistoGen->GetBinContent(ib);
	float XXX = 35.896*(ThisHisto->GetBinContent(ib)+ThisHistoGen->GetBinContent(ib))/2.;
        float EEE = 35.896*(ThisHisto->GetBinError(ib)+ThisHistoGen->GetBinError(ib))/2.;
	float value = EEE/XXX;
	if (entries) value = ThisHisto->GetBinContent(ib);
	EffMap->SetBinContent(ib, isr, value);
	//EffMap->SetBinError(ib, isr, EEE);
	
      }	

      isr++;
      
    }
  
  gStyle->SetOptStat(0);
  //gStyle->SetPaintTextFormat("5.2f");

  TCanvas *CC = new TCanvas("CC", "", 600, 400);
  CC->Divide(1, 1); 

  EffMap->Draw("COLZtext");
  //EffMap->Draw("COLZtexte");
  
  TString channelFlag = (Channel=="_em") ? "_emu" : "_eemumu";
  TString PlotName = entries ? "SignalEntries" : "SignalStatError";
  CC->Print("../Plots/MapMySignal/" + PlotName + "Map_" + Signal + channelFlag + ".pdf");
  CC->Print("../Plots/MapMySignal/" + PlotName + "Map_" + Signal + channelFlag + ".png");

}
