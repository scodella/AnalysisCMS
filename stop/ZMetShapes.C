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

float Lumi = 35.867;

TString RootFilesDirectory = "../minitrees/rootfiles3R/";

int const nObservables = 2;
TString ObservableName[nObservables] = {"MT2ll", "MT2llisr"};

int const nChannels = 4;
TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};

float StatZero = 1.84102;

// DY NoJet
float SF_DY_0     = 4.06;
float SF_DY_0_err = 2.38;

TH1F *GetDYNoJet(TFile* InputFile, TString HistoName, bool doStat, int doSyst) {

  TH1F *ThisHisto = (TH1F*) InputFile->Get(HistoName);
  ThisHisto->SetDirectory(0);

  float NativeAverageWeight = ThisHisto->Integral()/ThisHisto->GetEntries();

  float ThisSF = SF_DY_0 + doSyst*SF_DY_0_err;

  ThisHisto->Scale(ThisSF);

  if (doStat) {
    int foundZeros = 0, ib = 1;
    while (foundZeros<2 && ib<=ThisHisto->GetNbinsX()) {
      if (ThisHisto->GetBinContent(ib)==0.) {
	if (foundZeros==0) ThisHisto->SetBinError(ib, StatZero*ThisSF*NativeAverageWeight);
	if (foundZeros==1) ThisHisto->SetBinError(ib, StatZero*NativeAverageWeight);
	foundZeros++;
      }
      ib++;
    }
  }

  for (int ib = 1; ib<=ThisHisto->GetNbinsX(); ib++) 
    if (ThisHisto->GetBinContent(ib)==0.) 
      ThisHisto->SetBinContent(ib, 0.001/Lumi);
      
  return ThisHisto;

}

float SF_DY_1     = 1.58;
float SF_DY_1_err = 0.50;

TH1F *GetDYJet(TFile* InputFile, TString HistoName, bool doStat, int doSyst) {

  TH1F *ThisHisto = (TH1F*) InputFile->Get(HistoName);
  ThisHisto->SetDirectory(0);
  
  float NativeAverageWeight = ThisHisto->Integral()/ThisHisto->GetEntries();
  
  float ThisSF = SF_DY_1 + doSyst*SF_DY_1_err;

  for (int ib = 6; ib<=7; ib++) {
    ThisHisto->SetBinContent(ib, ThisSF*ThisHisto->GetBinContent(ib));
    ThisHisto->SetBinError  (ib, ThisSF*ThisHisto->GetBinError(ib));
  }

  if (doStat) {
    int foundZeros = 0, ib = 1;
    while (foundZeros<1 && ib<=4) {
      if (ThisHisto->GetBinContent(ib)==0.) {
	if (foundZeros==0) ThisHisto->SetBinError(ib, StatZero*NativeAverageWeight);
	foundZeros++;
      }
      ib++;
    }
    foundZeros = 0, ib = 7;
    while (foundZeros<1 && ib>=4) {
      if (ThisHisto->GetBinContent(ib)==0.) {
	if (foundZeros==0) ThisHisto->SetBinError(ib, StatZero*NativeAverageWeight);
	foundZeros++;
      }
      ib--;
    }
  }

  for (int ib = 1; ib<=ThisHisto->GetNbinsX(); ib++) 
    if (ThisHisto->GetBinContent(ib)==0.) 
      ThisHisto->SetBinContent(ib, 0.001/Lumi);
      
  return ThisHisto;
  
}

float SF_ZZ_0     = 1.00;
float SF_ZZ_0_err = 0.10;

TH1F *GetZZNoJet(TFile* InputFile, TString HistoName, bool doStat, int doSyst) {

  TH1F *ThisHisto = (TH1F*) InputFile->Get(HistoName);
  ThisHisto->SetDirectory(0);
  
  float NativeAverageWeight = ThisHisto->Integral()/ThisHisto->GetEntries();
  
  float ThisSF = SF_ZZ_0 + doSyst*SF_ZZ_0_err;

  ThisHisto->Scale(ThisSF);

  if (doStat) {
    for (int ib = 0; ib<=ThisHisto->GetNbinsX(); ib++) 
      if (ThisHisto->GetBinContent(ib)==0.)
	ThisHisto->SetBinError(ib, StatZero*NativeAverageWeight);
  }

  for (int ib = 1; ib<=ThisHisto->GetNbinsX(); ib++) 
    if (ThisHisto->GetBinContent(ib)==0.) 
      ThisHisto->SetBinContent(ib, 0.001/Lumi);
      
  return ThisHisto;

}

float SF_ZZ_1     = 1.51;
float SF_ZZ_1_err  = 0.21;

TH1F *GetZZJet(TFile* InputFile, TString HistoName, bool doStat, int doSyst) {

  TH1F *ThisHisto = (TH1F*) InputFile->Get(HistoName);
  ThisHisto->SetDirectory(0);
  
  float NativeAverageWeight = ThisHisto->Integral()/ThisHisto->GetEntries();
  
  float ThisSF = SF_ZZ_1 + doSyst*SF_ZZ_1_err;

  ThisHisto->Scale(ThisSF);

  if (doStat) {
    for (int ib = 0; ib<=ThisHisto->GetNbinsX(); ib++) 
      if (ThisHisto->GetBinContent(ib)==0.)
	ThisHisto->SetBinError(ib, StatZero*NativeAverageWeight);
  }

  for (int ib = 1; ib<=ThisHisto->GetNbinsX(); ib++) 
    if (ThisHisto->GetBinContent(ib)==0.) 
      ThisHisto->SetBinContent(ib, 0.001/Lumi);
      
  return ThisHisto;

}

void ZMetShapeSystematics(TString SystematicName) {
  
  for (int vr = 0; vr<2; vr++) { 

    TString Variation = vr==0 ? "up" : "do";

    TFile *ZZRootFile = TFile::Open(RootFilesDirectory + SystematicName + Variation + "/Stop/03_VZ.root");
    TFile *DYRootFile = TFile::Open(RootFilesDirectory + SystematicName + Variation + "/Stop/07_ZJetsHT.root");
          
    TFile *OutputFileZZ = new TFile(RootFilesDirectory + SystematicName + Variation + "/Stop/03_VZ_scaled.root", "recreate");
    TFile *OutputFileDY = new TFile(RootFilesDirectory + SystematicName + Variation + "/Stop/07_ZJetsHT_scaled.root", "recreate");
  
    for (int sr = 0; sr<ncut; sr++) {
      
      if (!scut[sr].Contains("_SR")) continue;
      if (scut[sr].Contains("_SRs")) continue;
      if (!scut[sr].Contains("_NoJet")) continue;
        
      TString NoTagRegion = scut[sr]; NoTagRegion.ReplaceAll("_NoJet", "_NoTag");
      TString TagRegion   = scut[sr]; TagRegion.ReplaceAll("_NoJet", "_Tag");
        
      OutputFileZZ->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
      
      OutputFileDY->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
    
      for (int ob = 0; ob<nObservables; ob++) 
	for (int ch = 0; ch<nChannels; ch++) {
	  
	  TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];
	  
	  OutputFileZZ->cd(scut[sr]);
	  TH1F *OutputHistoZZNoJet = GetZZNoJet(ZZRootFile, scut[sr] + "/" + HistoName, false, 0);
	  OutputHistoZZNoJet->Write();
	  
	  OutputFileDY->cd(scut[sr]);
	  TH1F *OutputHistoDYNoJet = GetDYNoJet(DYRootFile, scut[sr] + "/" + HistoName, false, 0);
	  OutputHistoDYNoJet->Write();
	  
	  OutputFileZZ->cd(NoTagRegion);
	  TH1F *OutputHistoZZJet = GetZZJet(ZZRootFile, NoTagRegion + "/" + HistoName, false, 0);
	  OutputHistoZZJet->Write();
	  
	  OutputFileDY->cd(NoTagRegion);
	  TH1F *OutputHistoDYJet = GetDYJet(DYRootFile, NoTagRegion + "/" + HistoName, false, 0);
	  OutputHistoDYJet->Write();
	  
	  OutputFileZZ->cd(TagRegion);
	  TH1F *OutputHistoZZTag = GetZZJet(ZZRootFile, TagRegion + "/" + HistoName, false, 0);
	  OutputHistoZZTag->Write();
	  
	  OutputFileDY->cd(TagRegion);
	  TH1F *OutputHistoDYTag = GetDYJet(DYRootFile, TagRegion + "/" + HistoName, false, 0);
	  OutputHistoDYTag->Write();

	}
	  	
    }
   
    ZZRootFile->Close();
    DYRootFile->Close();

    OutputFileZZ->Close();
    OutputFileDY->Close();

  }
   
}

void ZMetSpecificSystematics() {

  gSystem->mkdir(RootFilesDirectory + "ZZnojetup/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "DYnojetup/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "ZMETjetup/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "ZZnojetdo/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "DYnojetdo/Stop/", kTRUE);
  gSystem->mkdir(RootFilesDirectory + "ZMETjetdo/Stop/", kTRUE);

  TFile *ZZRootFile = TFile::Open(RootFilesDirectory + "nominal/Stop/03_VZ.root");
  TFile *DYRootFile = TFile::Open(RootFilesDirectory + "nominal/Stop/07_ZJetsHT.root");

  for (int vr = 0; vr<2; vr++) { 

    TString Variation = vr==0 ? "up" : "do";
          
    TFile *OutputFileZZ_ZZnojet = new TFile(RootFilesDirectory + "ZZnojet" + Variation + "/Stop/03_VZ_scaled.root", "recreate");
    TFile *OutputFileDY_DYnojet = new TFile(RootFilesDirectory + "DYnojet" + Variation + "/Stop/07_ZJetsHT_scaled.root", "recreate");
    TFile *OutputFileZZ_ZMETjet = new TFile(RootFilesDirectory + "ZMETjet" + Variation + "/Stop/03_VZ_scaled.root", "recreate");
    TFile *OutputFileDY_ZMETjet = new TFile(RootFilesDirectory + "ZMETjet" + Variation + "/Stop/07_ZJetsHT_scaled.root", "recreate");
  
    for (int sr = 0; sr<ncut; sr++) {
      
      if (!scut[sr].Contains("_SR")) continue;
      if (scut[sr].Contains("_SRs")) continue;
      if (!scut[sr].Contains("_NoJet")) continue;
        
      TString NoTagRegion = scut[sr]; NoTagRegion.ReplaceAll("_NoJet", "_NoTag");
      TString TagRegion   = scut[sr]; TagRegion.ReplaceAll("_NoJet", "_Tag");
        
      OutputFileZZ_ZZnojet->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);

      OutputFileZZ_ZMETjet->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
      
      OutputFileDY_DYnojet->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);

      OutputFileDY_ZMETjet->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
    
      for (int ob = 0; ob<nObservables; ob++) 
	for (int ch = 0; ch<nChannels; ch++) {
	  
	  TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];

	  // ZZnojet	  
	  OutputFileZZ_ZZnojet->cd(scut[sr]);
	  TH1F *OutputHistoZZ_NoJet_ZZnojet = GetZZNoJet(ZZRootFile, scut[sr] + "/" + HistoName, false, 1-2*vr);
	  OutputHistoZZ_NoJet_ZZnojet->Write();

	  OutputFileZZ_ZZnojet->cd(NoTagRegion);
	  TH1F *OutputHistoZZ_NoTag_ZZnojet = GetZZJet(ZZRootFile, NoTagRegion + "/" + HistoName, false, 0);
	  OutputHistoZZ_NoTag_ZZnojet->Write();

	  OutputFileZZ_ZZnojet->cd(TagRegion);
	  TH1F *OutputHistoZZ_Tag_ZZnojet = GetZZJet(ZZRootFile, TagRegion + "/" + HistoName, false, 0);
	  OutputHistoZZ_Tag_ZZnojet->Write();
	  
	  // DYnojet
	  OutputFileDY_DYnojet->cd(scut[sr]);
	  TH1F *OutputHistoDY_NoJet_DYnojet = GetDYNoJet(DYRootFile, scut[sr] + "/" + HistoName, false, 1-2*vr);
	  OutputHistoDY_NoJet_DYnojet->Write();

	  OutputFileDY_DYnojet->cd(NoTagRegion);
	  TH1F *OutputHistoDY_NoTag_DYnojet = GetDYJet(DYRootFile, NoTagRegion + "/" + HistoName, false, 0);
	  OutputHistoDY_NoTag_DYnojet->Write();

	  OutputFileDY_DYnojet->cd(TagRegion);
	  TH1F *OutputHistoDY_Tag_DYnojet = GetDYJet(DYRootFile, TagRegion + "/" + HistoName, false, 0);
	  OutputHistoDY_Tag_DYnojet->Write();

	  // ZMETjet ZZ
	  OutputFileZZ_ZMETjet->cd(scut[sr]);
	  TH1F *OutputHistoZZ_NoJet_ZMETjet = GetZZNoJet(ZZRootFile, scut[sr] + "/" + HistoName, false, 0);
	  OutputHistoZZ_NoJet_ZMETjet->Write();

	  OutputFileZZ_ZMETjet->cd(NoTagRegion);
	  TH1F *OutputHistoZZ_NoTag_ZMETjet = GetZZJet(ZZRootFile, NoTagRegion + "/" + HistoName, false, 1-2*vr);
	  OutputHistoZZ_NoTag_ZMETjet->Write();

	  OutputFileZZ_ZMETjet->cd(TagRegion);
	  TH1F *OutputHistoZZ_Tag_ZMETjet = GetZZJet(ZZRootFile, TagRegion + "/" + HistoName, false, 1-2*vr);
	  OutputHistoZZ_Tag_ZMETjet->Write();

	  // ZMETjet DY
	  OutputFileDY_ZMETjet->cd(scut[sr]);
	  TH1F *OutputHistoDY_NoJet_ZMETjet = GetDYNoJet(DYRootFile, scut[sr] + "/" + HistoName, false, 0);
	  OutputHistoDY_NoJet_ZMETjet->Write();

	  OutputFileDY_ZMETjet->cd(NoTagRegion);
	  TH1F *OutputHistoDY_NoTag_ZMETjet = GetDYJet(DYRootFile, NoTagRegion + "/" + HistoName, false, 2*vr-1);
	  OutputHistoDY_NoTag_ZMETjet->Write();

	  OutputFileDY_ZMETjet->cd(TagRegion);
	  TH1F *OutputHistoDY_Tag_ZMETjet = GetDYJet(DYRootFile, TagRegion + "/" + HistoName, false, 2*vr-1);
	  OutputHistoDY_Tag_ZMETjet->Write();

	}
	  	
    }
  
    OutputFileZZ_ZZnojet->Close();
    OutputFileDY_DYnojet->Close();
    OutputFileZZ_ZMETjet->Close();
    OutputFileDY_ZMETjet->Close();
   
  }

  ZZRootFile->Close();
  DYRootFile->Close();
   
}

void ZMetShapesNominal() {
  
  TFile *ZZRootFile = TFile::Open(RootFilesDirectory + "nominal/Stop/03_VZ.root");
  TFile *DYRootFile = TFile::Open(RootFilesDirectory + "nominal/Stop/07_ZJetsHT.root");
  
  TFile *OutputFileZZ = new TFile(RootFilesDirectory + "nominal/Stop/03_VZ_scaled.root", "recreate");
  TFile *OutputFileDY = new TFile(RootFilesDirectory + "nominal/Stop/07_ZJetsHT_scaled.root", "recreate");
  
    for (int sr = 0; sr<ncut; sr++) {
      
      if (!scut[sr].Contains("_SR")) continue;
      if (scut[sr].Contains("_SRs")) continue;
      if (!scut[sr].Contains("_NoJet")) continue;
        
      TString NoTagRegion = scut[sr]; NoTagRegion.ReplaceAll("_NoJet", "_NoTag");
      TString TagRegion   = scut[sr]; TagRegion.ReplaceAll("_NoJet", "_Tag");
        
      OutputFileZZ->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
      
      OutputFileDY->cd();
      gDirectory->mkdir(scut[sr]);
      gDirectory->mkdir(NoTagRegion);
      gDirectory->mkdir(TagRegion);
    
      for (int ob = 0; ob<nObservables; ob++) 
	for (int ch = 0; ch<nChannels; ch++) {
	  
	  TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];
	  
	  OutputFileZZ->cd(scut[sr]);
	  TH1F *OutputHistoZZNoJet = GetZZNoJet(ZZRootFile, scut[sr] + "/" + HistoName, true, 0);
	  OutputHistoZZNoJet->Write();
	  
	  OutputFileDY->cd(scut[sr]);
	  TH1F *OutputHistoDYNoJet = GetDYNoJet(DYRootFile, scut[sr] + "/" + HistoName, true, 0);
	  OutputHistoDYNoJet->Write();
	  
	  OutputFileZZ->cd(NoTagRegion);
	  TH1F *OutputHistoZZJet = GetZZJet(ZZRootFile, NoTagRegion + "/" + HistoName, true, 0);
	  OutputHistoZZJet->Write();
	  
	  OutputFileDY->cd(NoTagRegion);
	  TH1F *OutputHistoDYJet = GetDYJet(DYRootFile, NoTagRegion + "/" + HistoName, true, 0);
	  OutputHistoDYJet->Write();
	  
	  OutputFileZZ->cd(TagRegion);
	  TH1F *OutputHistoZZTag = GetZZJet(ZZRootFile, TagRegion + "/" + HistoName, true, 0);
	  OutputHistoZZTag->Write();
	  
	  OutputFileDY->cd(TagRegion);
	  TH1F *OutputHistoDYTag = GetDYJet(DYRootFile, TagRegion + "/" + HistoName, true, 0);
	  OutputHistoDYTag->Write();

	}
      
    }
    
    ZZRootFile->Close();
    DYRootFile->Close();
    
    OutputFileZZ->Close();
    OutputFileDY->Close();
    
}

void ZMetShapes() {

  ZMetShapeSystematics("Btag");
  ZMetShapeSystematics("Idiso");
  ZMetShapeSystematics("JES");
  ZMetShapeSystematics("MET");
  ZMetShapeSystematics("PDF");
  ZMetShapeSystematics("Q2");
  ZMetShapeSystematics("Reco");

  ZMetSpecificSystematics();

  ZMetShapesNominal();

}
