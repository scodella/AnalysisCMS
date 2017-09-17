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

void JetPlots() {

  TFile *StopW = TFile::Open("../rootfiles/nominal/Stop/T2tt_mStop-350to400_Sm350_Xm263.root");
  TFile *StopC = TFile::Open("../rootfiles/nominal/Stop/T2tt_mStop-350to400_Sm350_Xm225.root");
  TFile *StopT = TFile::Open("../rootfiles/nominal/Stop/T2tt_mStop-350to400_Sm350_Xm175.root");
  TFile *TTbar = TFile::Open("../rootfiles/nominal/Stop/04_TTTo2L2Nu.root");
  TFile *WW    = TFile::Open("../rootfiles/nominal/Stop/06_WW.root");

  TH1F *TagStopW = (TH1F*) StopW->Get("Stop/00_Zveto/h_nbjet30csvv2m");
  TH1F *TagStopC = (TH1F*) StopC->Get("Stop/00_Zveto/h_nbjet30csvv2m");
  TH1F *TagStopT = (TH1F*) StopT->Get("Stop/00_Zveto/h_nbjet30csvv2m");
  TH1F *TagTTbar = (TH1F*) TTbar->Get("Stop/00_Zveto/h_nbjet30csvv2m");
  TH1F *TagWW    = (TH1F*) WW   ->Get("Stop/00_Zveto/h_nbjet30csvv2m");

  TH1F *JetStopW = (TH1F*) StopW->Get("Stop/00_Zveto/h_njet");
  TH1F *JetStopC = (TH1F*) StopC->Get("Stop/00_Zveto/h_njet");
  TH1F *JetStopT = (TH1F*) StopT->Get("Stop/00_Zveto/h_njet");
  TH1F *JetTTbar = (TH1F*) TTbar->Get("Stop/00_Zveto/h_njet");
  TH1F *JetWW    = (TH1F*) WW   ->Get("Stop/00_Zveto/h_njet");

  TH1F *JetVetoStopW = (TH1F*) StopW->Get("Stop/01_NoTag/h_njet");
  TH1F *JetVetoStopC = (TH1F*) StopC->Get("Stop/01_NoTag/h_njet");
  TH1F *JetVetoStopT = (TH1F*) StopT->Get("Stop/01_NoTag/h_njet");
  TH1F *JetVetoTTbar = (TH1F*) TTbar->Get("Stop/01_NoTag/h_njet");
  TH1F *JetVetoWW    = (TH1F*) WW   ->Get("Stop/01_NoTag/h_njet");

  TCanvas *CC = new TCanvas("CC", "", 600, 400);
  CC->cd(); CC->SetBorderMode(0); CC->SetFillColor(0);
  
  TagStopW->SetMarkerStyle(20);
  TagStopW->SetMarkerColor(kOrange);
  TagStopC->SetMarkerStyle(21);
  TagStopC->SetMarkerColor(kViolet);
  TagStopT->SetMarkerStyle(22);
  TagStopT->SetMarkerColor(kRed);
  TagTTbar->SetMarkerStyle(23);
  TagTTbar->SetMarkerColor(kYellow);
  TagWW   ->SetMarkerStyle(24);
  TagWW   ->SetMarkerColor(kAzur-9);

  TagStop->SetXTitle("");

}
