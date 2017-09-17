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

void MT2GenReco(TString Sample, TString Type, TString TXm = "") {

  const int nChannels = 4, nVariables = 2, nMetCuts = 5, nBtagCuts = 2;
  TString Channel[nChannels] = {"_ee", "_mm", "_em", "_ll"};
  TString Variable[nVariables] = {"MT2ll", "MET"};
  TString MetCut[nMetCuts] = {"01", "02_VR1", "02_SR1", "02_SR2", "02_SR3"};
  TString BtagCut[nBtagCuts] = {"_NoTag", "_Tag"};
  
  TString RootFileName = "";
  if (Sample=="Signal") RootFileName = "../minitrees/rootfiles/nominal/Stop/T2tt_isr_mStop-" + Type + ".root";
  if (Sample=="ttbarpowheg") RootFileName = "../minitrees/rootfiles/nominal/Stop/04_TTTo2L2Nu.root";

  if (RootFileName=="") return;

  TFile *RootFile = TFile::Open(RootFileName);

  TCanvas *CL = new TCanvas("CL", "", 900, 300); 
  CL->Divide(3, 1);
  
  TPad *PadL1 = (TPad*)CL->GetPad(1);
  TPad *PadL2 = (TPad*)CL->GetPad(2);
  TPad *PadL3 = (TPad*)CL->GetPad(3);

  PadL1->SetLogy();
  PadL2->SetLogy();
  PadL3->SetLogy();

  // Run2015B Setting
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  CL->Range(0,0,1,1);
  CL->SetFillColor(10);
  CL->SetBorderMode(0);
  CL->SetBorderSize(2);
  CL->SetTickx(1);
  CL->SetTicky(1);
  CL->SetLeftMargin(0.16);
  CL->SetRightMargin(0.02);
  CL->SetTopMargin(0.05);
  CL->SetBottomMargin(0.13);
  CL->SetFrameFillColor(0);
  CL->SetFrameFillStyle(0);
  CL->SetFrameBorderMode(0);
  
  PadL1->SetFillColor(0);
  PadL1->SetBorderMode(0);
  PadL1->SetBorderSize(2);
  //PadL1->SetGridy();
  //PadL1->SetLogx();
  PadL1->SetTickx(1);
  PadL1->SetTicky(1);
  PadL1->SetLeftMargin(0.16);
  PadL1->SetRightMargin(0.02);
  //PadL1->SetTopMargin(0.05);
  //PadL1->SetBottomMargin(0.31);
  PadL1->SetTopMargin(0.065);
  PadL1->SetBottomMargin(0.13);
  PadL1->SetFrameFillStyle(0);
  PadL1->SetFrameBorderMode(0);
  PadL1->SetFrameFillStyle(0);
  PadL1->SetFrameBorderMode(0);
  PadL1->Draw();
  
  PadL2->SetFillColor(0);
  PadL2->SetBorderMode(0);
  PadL2->SetBorderSize(2);
  //PadL2->SetGridy();
  //PadL2->SetLogx();
  PadL2->SetTickx(1);
  PadL2->SetTicky(1);
  PadL2->SetLeftMargin(0.16);
  PadL2->SetRightMargin(0.02);
  //PadL2->SetTopMargin(0.05);
  //PadL2->SetBottomMargin(0.31);
  PadL2->SetTopMargin(0.065);
  PadL2->SetBottomMargin(0.13);
  PadL2->SetFrameFillStyle(0);
  PadL2->SetFrameBorderMode(0);
  PadL2->SetFrameFillStyle(0);
  PadL2->SetFrameBorderMode(0);
  PadL2->Draw();
  
  PadL3->SetFillColor(0);
  PadL3->SetBorderMode(0);
  PadL3->SetBorderSize(2);
  //PadL3->SetGridy();
  //PadL3->SetLogx();
  PadL3->SetTickx(1);
  PadL3->SetTicky(1);
  PadL3->SetLeftMargin(0.16);
  PadL3->SetRightMargin(0.02);
  //PadL3->SetTopMargin(0.05);
  //PadL3->SetBottomMargin(0.31);
  PadL3->SetTopMargin(0.065);
  PadL3->SetBottomMargin(0.13);
  PadL3->SetFrameFillStyle(0);
  PadL3->SetFrameBorderMode(0);
  PadL3->SetFrameFillStyle(0);
  PadL3->SetFrameBorderMode(0);
  PadL3->Draw();
  // End Run2015B Setting 

  float lumi = 9.983;
  TString Lumi = "9.983";// Lumi += lumi;      
  
  TLatex *text1 = new TLatex(0.96,0.95125, Lumi + " fb^{-1} (13 Tev)"); 
  text1->SetNDC();                                              
  text1->SetTextAlign(31);                          
  text1->SetTextFont(42);    
  text1->SetTextSize(0.04);   
  text1->SetLineWidth(2);    
  text1->Draw();  

  for (int sr = 1; sr<=3; sr++) {
  
    TString HistoDir = "Stop/02_SR"; HistoDir += sr;
    TH1F *SRTR = (TH1F*) RootFile->Get(HistoDir + "_Tag/h_MT2llgen_ll");
    TH1F *SRVR = (TH1F*) RootFile->Get(HistoDir + "_NoTag/h_MT2llgen_ll");    
    TH1F *SRTG = (TH1F*) RootFile->Get(HistoDir + "gen_Tag/h_MT2llgen_ll");
    TH1F *SRVG = (TH1F*) RootFile->Get(HistoDir + "gen_NoTag/h_MT2llgen_ll");

    SRTR->Add(SRVR); SRTG->Add(SRVG); 

    if (sr==1) PadL1->cd();
    if (sr==2) PadL2->cd();
    if (sr==3) PadL3->cd();

    SRTR->SetLineWidth(2.); 
    SRTG->SetLineWidth(2.); 
    SRTG->SetLineColor(2); 

    SRTR->Scale(lumi);
    SRTG->Scale(lumi);

    SRTR->SetXTitle("M_{T2}(ll) [GeV]");
    SRTR->SetYTitle("Events / 20 GeV");
    
    SRTR->GetXaxis()->SetLabelFont(42);
    SRTR->GetYaxis()->SetLabelFont(42);
    SRTR->GetXaxis()->SetTitleFont(42);
    SRTR->GetYaxis()->SetTitleFont(42);
    
    SRTR->GetYaxis()->SetTitleSize(0.06);
    SRTR->GetYaxis()->SetLabelSize(0.05);
    SRTR->GetXaxis()->SetTitleSize(0.06);
    SRTR->GetXaxis()->SetLabelSize(0.05);
    SRTR->GetXaxis()->SetTitleOffset(0.95);
    SRTR->GetYaxis()->SetTitleOffset(1.25);
    
    cout << SRTR->Integral() << " " << SRTG->Integral() << " " << (SRTR->Integral()-SRTG->Integral())/SRTR->Integral() << endl;

    SRTR->DrawCopy("histo");
    SRTG->DrawCopy("histosame");

    TLegend  *leg = new TLegend(0.23, 0.2, 0.39, 0.3);
    leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
    leg->SetTextColor(1); leg->SetTextSize(0.035);
    leg->SetTextFont(62);
    TString LegTitle;
    if (Sample=="Signal") LegTitle = "T2tt " + TXm;
    else if (Sample=="ttbarpowheg") LegTitle = "t#bar{t} powheg";
    leg->SetHeader(LegTitle);
    leg->AddEntry(SRTR, "pf-MET",  "l");
    leg->AddEntry(SRTG, "gen-MET", "l");

    leg->Draw();

  }

  TString PlotName = RootFileName;
  PlotName.ReplaceAll("../minitrees/rootfiles/nominal/Stop/", "");
  PlotName.ReplaceAll(".root", "");
  //CL->Print("../Plots/MT2GenReco/" + PlotName + ".png");

}
