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

int const maxHistos = 10;
TH1F *Histo[maxHistos];
int HistoColor[maxHistos], nHistos, Rebinning;
int HistoMarker[maxHistos];
float MaxXValue, MinYValue, MaxYValue;
float LegX1, LegY1, LegX2, LegY2;
TString XTitle, YTitle, HistoLegend[maxHistos];
bool LogScale, Normalize, UseMarker, SumOverflow;

void MakePlot(TString PlotTitle) {

  TCanvas *CC = new TCanvas("CC", "", 600, 400);
  CC->Divide(1, 1); 
  TPad *PD = (TPad*) CC->GetPad(1);
  if (LogScale)  PD->SetLogy();
  PD->cd();

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  CC->Range(0,0,1,1);
  CC->SetFillColor(10);
  CC->SetBorderMode(0);
  CC->SetBorderSize(2);
  CC->SetTickx(1);
  CC->SetTicky(1);
  CC->SetLeftMargin(0.16);
  CC->SetRightMargin(0.02);
  CC->SetTopMargin(0.05);
  CC->SetBottomMargin(0.13);
  CC->SetFrameFillColor(0);
  CC->SetFrameFillStyle(0);
  CC->SetFrameBorderMode(0);
      
  TLegend  *leg = new TLegend(LegX1, LegY1, LegX2, LegY2);
  leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
  leg->SetTextColor(1); leg->SetTextSize(0.035);
  leg->SetTextFont(62); 
  //leg->SetHeader(TagLegend[tagflag]); //leg->SetMargin(0.2); 

  int nBins = Histo[0]->GetNbinsX();
  cout <<"S"<<endl;
  for (int hs = 0; hs<nHistos; hs++) {

    Histo[hs]->Rebin(Rebinning);
cout <<"A"<<endl;
    if (Normalize) 
      Histo[hs]->Scale(1./Histo[hs]->Integral(0, nBins+1));
    
    if (SumOverflow) { 
      int LastBin = nBins;
      if (MaxXValue>-999.) LastBin = Histo[hs]->GetBin(MaxXValue);
      float Overflow = Histo[hs]->Integral(LastBin, nBins+1);
      float OverflowError = 0.;
      for (int ib = LastBin; ib<=nBins+1; ib++) 
      OverflowError = Histo[hs]->GetBinError(ib)*Histo[hs]->GetBinError(ib);
	OverflowError = sqrt(OverflowError);
      Histo[hs]->SetBinContent(nBins, Overflow);
      Histo[hs]->SetBinError(nBins, OverflowError);
    }
    cout <<"A"<<endl;
    if (MaxXValue>-999.) {
      float MinXValue = Histo[hs]->GetBinLowEdge(1);
      Histo[hs]->GetXaxis()->SetRangeUser(MinXValue, MaxXValue);
    }

    if (MinYValue>-999.) Histo[hs]->SetMinimum(MinYValue);
    if (MaxYValue>-999.) Histo[hs]->SetMaximum(MaxYValue);

    Histo[hs]->SetXTitle(XTitle);
    float BW = Histo[hs]->GetBinWidth(1);
    TString BinWidth;
    if (BW>1.) BinWidth = int(BW); // to fix
    else { BW *= 100; BinWidth = "0."; BinWidth += int(BW); }
    YTitle.ReplaceAll("BinWidth", BinWidth);
    Histo[hs]->SetYTitle(YTitle);
    Histo[hs]->GetYaxis()->CenterTitle();
    Histo[hs]->SetLineColor(HistoColor[hs]);
    Histo[hs]->SetMarkerColor(HistoColor[hs]);
    Histo[hs]->SetMarkerStyle(HistoMarker[hs]);

    TString LegendOption = UseMarker ? "lpe" : "f";
    leg->AddEntry(Histo[hs], HistoLegend[hs], LegendOption);
    
    TString HistoOption = hs==0 ? "" : "same";
    if (UseMarker) HistoOption += "pe";
    else HistoOption += "histo";
cout <<"A"<<endl;
    Histo[hs]->Draw(HistoOption);

  }

  leg->Draw();

  CC->Print("../Plots/ShapeComparisons/" + PlotTitle + ".png");

}

TH1F* GetHistogram(TFile *inputFile, TString Cut, TString Observable) {

  TH1F *ThisHisto;

  if (!Observable.Contains("_sf")) 
    ThisHisto = (TH1F*) inputFile->Get("Stop/" + Cut + "/h_" + Observable);
  else {
    Observable.ReplaceAll("_sf", "_ee");
    ThisHisto = (TH1F*) inputFile->Get("Stop/" + Cut + "/h_" + Observable);
    Observable.ReplaceAll("_ee", "_mm");
    TH1F* ThisHisto2 = (TH1F*) inputFile->Get("Stop/" + Cut + "/h_" + Observable);
    ThisHisto->Add(ThisHisto2);
  }

  return ThisHisto;

}

TH1F* GetHistogramCut(TFile *inputFile, TString Region, TString Cut, TString Observable) {

  TH1F *NoJet = GetHistogram(inputFile, Region + "NoJet", Observable);
  TH1F *NoTag = GetHistogram(inputFile, Region + "NoTag", Observable);

  if (Cut=="NoJet") return NoJet;
  if (Cut=="NoTag") return NoTag;

  NoJet->Add(NoTag);
  return NoJet;

}

TH1F* GetSRs(TFile *inputFile, TString Cut, TString Observable) {

  TH1F *SR1 = GetHistogramCut(inputFile, "02_SR1_", Cut, Observable);
  TH1F *SR2 = GetHistogramCut(inputFile, "02_SR2_", Cut, Observable);
  TH1F *SR3 = GetHistogramCut(inputFile, "02_SR3_", Cut, Observable);

  SR1->Add(SR2); SR1->Add(SR3);
  return SR1;

}

void WZShapes() {

  TFile *WZTo2TV = TFile::Open("../minitrees/rootfiles/nominalK/Stop/02_WZTo3LNu.root");
 
  Histo[0] = GetHistogram(WZTo2TV, "01_NoTag", "MET_ll");
  HistoColor[0] = 1;
  HistoLegend[0] = "WZ: 2 tight and 1 lost leptons";
  HistoMarker[0] = 20;

  TFile *WZTo3T = TFile::Open("../minitrees/rootfiles/invertveto/Stop/WZTo3LNu.root");

  Histo[1] = GetHistogram(WZTo3T, "01_NoTag", "MET_ll");
  HistoColor[1] = 2;
  HistoLegend[1] = "WZ: 3 tight leptons";
  HistoMarker[1] = 21;

  MaxXValue = 400.; MinYValue = 1e-04; MaxYValue = -999.;
  nHistos = 2; Rebinning = 1; SumOverflow = false;
  XTitle = "MET [GeV]"; YTitle = "events / BinWidth GeV";
  LogScale = true; Normalize = true; UseMarker = true;
  LegX1 = 0.5; LegY1 = 0.5; LegX2 = 0.8; LegY2 = 0.8;

  MakePlot("WZ_3T_vs_2T1L_MET_ll");

  // 
  Histo[0] = GetSRs(WZTo2TV, "Veto", "MT2ll_ll");
  HistoColor[0] = 1;
  HistoLegend[0] = "WZ: 2 tight and 1 lost leptons";
  HistoMarker[0] = 20;

  Histo[1] = GetSRs(WZTo3T, "Veto", "MT2ll_ll");
  HistoColor[1] = 2;
  HistoLegend[1] = "WZ: 3 tight leptons";
  HistoMarker[1] = 21;

  MaxXValue = -999.; MinYValue = 1e-02; MaxYValue = -999.;
  nHistos = 2; Rebinning = 1; SumOverflow = false;
  XTitle = "M_{T2}(ll) [GeV]"; YTitle = "events / BinWidth GeV";
  LogScale = true; Normalize = true; UseMarker = true;
  LegX1 = 0.5; LegY1 = 0.65; LegX2 = 0.8; LegY2 = 0.8;

  MakePlot("WZ_3T_vs_2T1L_MT2ll_ll");
 
}

void DYShapes() {

  TFile *DYZ = TFile::Open("../minitrees/rootfiles3R/ZpeakK/Stop/07_ZJetsHT.root");

  Histo[0] = GetSRs(DYZ, "NoJet", "dphiLL_ll");
  HistoColor[0] = 1;
  HistoLegend[0] = "DY: |M_{ll} - M_{Z}| < 15 GeV";
  HistoMarker[0] = 20;

  TFile *DYN = TFile::Open("../minitrees/rootfiles3R/nominalK/Stop/07_ZJetsHT.root");

  Histo[1] = GetSRs(DYN, "NoJet", "dphiLL_ll");
  HistoColor[1] = 2;
  HistoLegend[1] = "DY: |M_{ll} - M_{Z}| > 15 GeV";
  HistoMarker[1] = 21;

  MaxXValue = -999.; MinYValue = 1e-04; MaxYValue = -999.;
  nHistos = 2; Rebinning = 10; SumOverflow = false;
  XTitle = "#Delta#phi(ll)"; YTitle = "events / BinWidth ";
  LogScale = false; Normalize = true; UseMarker = true;
  LegX1 = 0.5; LegY1 = 0.5; LegX2 = 0.8; LegY2 = 0.8;

  MakePlot("DY_NoJet_DPhi");

}
