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

#include "../include/Constants.h"

bool PlotISR = true;

const int  nbackgrounds = 10;
  
const TString backname[nbackgrounds] = {
  "04_TTTo2L2Nu",
  "05_ST",
  "06_WW",
  "07_ZJetsHT",
  "10_TTZ",
  "09_TTW",
  "03_VZ",
  "02_WZTo3LNu",
  "13_VVV",
  "11_HWW"
};

const int nmasspoints = 3;

const TString masspointname[nmasspoints] = {
  "T2tt_isr_mStop-350to400_Sm350_Xm175",
  "T2tt_isr_mStop-350to400_Sm350_Xm225",
  "T2tt_isr_mStop-350to400_Sm350_Xm263"
  //"T2bW_isr_Sm350_Xm175",
  //"T2bW_isr_Sm350_Xm150"
};

const TString masspointlegend[nmasspoints] = {
  "m_{stop}=350, m_{X}=175",
  "m_{stop}=350, m_{X}=225",
  "m_{stop}=350, m_{X}=263"
};

/*const int nmasspoints = 2;

const TString masspointname[nmasspoints] = {
  "T2tt_isr_mStop-350to400_Sm350_Xm225",
  "T2tt_isrnorew_mStop-350to400_Sm350_Xm225"
};

const TString masspointlegend[nmasspoints] = {
  "m_{stop}=350, m_{X}=225 ISR rew.",
  "m_{stop}=350, m_{X}=225"
  };*/

TString HistogramFileDirectory = "../minitrees/rootfiles/nominal/Stop/";

void SignalBackground(TString FigureOfMerit = "Ratio", TString Variable = "MT2ll", TString Channel = "_ll", float Lumi = -1.) {

  if (Lumi==-1.) Lumi = lumi_fb_Full2016;

  int nSignalRegions = 0;
  TString SignalRegionName[100];
  
  for (int ct = 0; ct<ncut; ct++) {
    if (scut[ct].Contains("_SR") && scut[ct].Contains("_NoTag")) {

	SignalRegionName[nSignalRegions] = scut[ct];
	nSignalRegions++;
	
    }
  }
 
  if (nSignalRegions==0) {

    cout << "Error: no sigal regions found! Exiting" << endl;
    return;

  }

  TFile *TemplateFile = TFile::Open(HistogramFileDirectory + masspointname[0] + ".root", "read");
  TH1F *TemplateHisto = (TH1F*) TemplateFile->Get(SignalRegionName[0] + "/h_" + Variable + Channel);
  
  int nVariableBins = TemplateHisto->GetNbinsX();
  float LowerEdge = TemplateHisto->GetBinLowEdge(1);
  float HigherEdge = TemplateHisto->GetBinLowEdge(nVariableBins+1);

  TemplateFile->Close();

  TH1F *ResultHisto[nmasspoints][2], *ResultHistoISR[nmasspoints][2];

  float SignalYield[nmasspoints][2][100];
  float SignalError[nmasspoints][2][100];
  float SignalYieldISR[nmasspoints][2][100];
  float SignalErrorISR[nmasspoints][2][100];

  for (int mp = 0; mp<nmasspoints; mp++) {
    
    TString HistoName = "MassPoint"; HistoName += mp;  
    ResultHisto[mp][0] = new TH1F(HistoName + "_Tagged",   "", nSignalRegions*nVariableBins, LowerEdge, nSignalRegions*HigherEdge);
    ResultHisto[mp][1] = new TH1F(HistoName + "_NoTagged", "", nSignalRegions*nVariableBins, LowerEdge, nSignalRegions*HigherEdge);
    ResultHistoISR[mp][0] = new TH1F(HistoName + "_Tagged_ISRcut",   "", nSignalRegions*nVariableBins, LowerEdge, nSignalRegions*HigherEdge);
    ResultHistoISR[mp][1] = new TH1F(HistoName + "_NoTagged_ISRcut", "", nSignalRegions*nVariableBins, LowerEdge, nSignalRegions*HigherEdge);
    
    TFile *HistoFile = TFile::Open(HistogramFileDirectory + masspointname[mp] + ".root", "read");
  
    int longbin[2] = {1, 1};

    for (int ct = 0; ct<ncut; ct++) {
      if (scut[ct].Contains("_SR")) {
	
	int tagflag = 0;
	if (scut[ct].Contains("_NoTag")) tagflag = 1;

	TH1F *ThisHisto = (TH1F*) HistoFile->Get(scut[ct] + "/h_" + Variable + Channel);
	TH1F *ThisHistoISR = (TH1F*) HistoFile->Get(scut[ct] + "/h_" + Variable + "isr" + Channel);
	
	for (int ib = 1; ib<=nVariableBins; ib++) {

	  SignalYield[mp][tagflag][longbin[tagflag]] = Lumi*ThisHisto->GetBinContent(ib);
	  SignalError[mp][tagflag][longbin[tagflag]] = Lumi*ThisHisto->GetBinError(ib);
	  SignalYieldISR[mp][tagflag][longbin[tagflag]] = Lumi*ThisHistoISR->GetBinContent(ib);
	  SignalErrorISR[mp][tagflag][longbin[tagflag]] = Lumi*ThisHistoISR->GetBinError(ib);
	  longbin[tagflag]++;

	}

      }
    }

  }

  float BackgroundYield[2][100];
  float BackgroundYieldISR[2][100];

  for (int tagflag = 0; tagflag<2; tagflag++) 
    for (int ib = 1; ib<=nSignalRegions*nVariableBins; ib++) 
      BackgroundYield[tagflag][ib] = BackgroundYieldISR[tagflag][ib] = 0.;

  for (int mp = 0; mp<nbackgrounds; mp++) {
    
    TFile *HistoFile = TFile::Open(HistogramFileDirectory + backname[mp] + ".root", "read");
  
    int longbin[2] = {1, 1};

    for (int ct = 0; ct<ncut; ct++) {
      if (scut[ct].Contains("_SR")) {
	
	int tagflag = 0;
	if (scut[ct].Contains("_NoTag")) tagflag = 1;

	TH1F *ThisHisto = (TH1F*) HistoFile->Get(scut[ct] + "/h_" + Variable + Channel);
	TH1F *ThisHistoISR = (TH1F*) HistoFile->Get(scut[ct] + "/h_" + Variable + "isr" + Channel);
	
	for (int ib = 1; ib<=nVariableBins; ib++) {
	  
	  float ThisContent = ThisHisto->GetBinContent(ib);
	  if (ThisContent<0.) ThisContent = 0.001;
	  BackgroundYield[tagflag][longbin[tagflag]] += Lumi*ThisContent;
	  ThisContent = ThisHistoISR->GetBinContent(ib);
	  if (ThisContent<0.) ThisContent = 0.001;
	  BackgroundYieldISR[tagflag][longbin[tagflag]] += Lumi*ThisContent;
	  longbin[tagflag]++;

	}

      }
    }

  }

  for (int mp = 0; mp<nmasspoints; mp++) {
    for (int tagflag = 0; tagflag<2; tagflag++) {
      for (int ib = 1; ib<=nSignalRegions*nVariableBins; ib++) {
	
	float ThisResult = 0., ThisResultError = 0.;
	float ThisResultISR = 0., ThisResultErrorISR = 0.;

	if (FigureOfMerit=="Ratio") {
	  ThisResult = SignalYield[mp][tagflag][ib]/BackgroundYield[tagflag][ib];
	  ThisResultError = SignalError[mp][tagflag][ib]/BackgroundYield[tagflag][ib];
	  ThisResultISR = SignalYieldISR[mp][tagflag][ib]/BackgroundYieldISR[tagflag][ib];
	  ThisResultErrorISR = SignalErrorISR[mp][tagflag][ib]/BackgroundYieldISR[tagflag][ib];
	} else if (FigureOfMerit=="Significance") {
	  ThisResult = SignalYield[mp][tagflag][ib]/sqrt(SignalYield[mp][tagflag][ib]+BackgroundYield[tagflag][ib]);
	  ThisResultISR = SignalYieldISR[mp][tagflag][ib]/sqrt(SignalYieldISR[mp][tagflag][ib]+BackgroundYieldISR[tagflag][ib]);
	} else if (FigureOfMerit=="SignificanceSyst") {
	  ThisResult = SignalYield[mp][tagflag][ib]/sqrt(SignalYield[mp][tagflag][ib]+BackgroundYield[tagflag][ib]+pow(0.1*BackgroundYield[tagflag][ib],2));
	  ThisResultISR = SignalYieldISR[mp][tagflag][ib]/sqrt(SignalYieldISR[mp][tagflag][ib]+BackgroundYieldISR[tagflag][ib]+pow(0.1*BackgroundYieldISR[tagflag][ib],2));
	} else if (FigureOfMerit=="Signal") {
	  ThisResult = SignalYield[mp][tagflag][ib];
	  ThisResultISR = SignalYieldISR[mp][tagflag][ib];
	} else if (FigureOfMerit=="Background") {
	  ThisResult = BackgroundYield[tagflag][ib];
	  ThisResultISR = BackgroundYieldISR[tagflag][ib];
	}

	ResultHisto[mp][tagflag]->SetBinContent(ib, ThisResult );
	ResultHisto[mp][tagflag]->SetBinError(ib, ThisResultError );
	ResultHistoISR[mp][tagflag]->SetBinContent(ib, ThisResultISR );
	ResultHistoISR[mp][tagflag]->SetBinError(ib, ThisResultErrorISR );
	
      }
    }
  }

  TCanvas *CC = new TCanvas("CC", "", 600, 400);
  CC->Divide(1, 1); 
  TPad *PD = (TPad*) CC->GetPad(1);
  if (FigureOfMerit=="Signal" || FigureOfMerit=="Background") PD->SetLogy();
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

  TString TagFlag[2] = {"_Tag", "_NoTag"};
  TString TagLegend[2] = {"Events with >=1 b-tag", "Events with no b-tags"};
  TString LumiFlag = "_Lumi"; LumiFlag += int(Lumi+0.5);
  
  int MassPointColor[3] = {1, 2, 4};
  int MassPointMarker[3] = {20, 25, 22};

    TString TitleY = FigureOfMerit;
    if (FigureOfMerit=="Ratio") TitleY = "Signal / Background";

  for (int tagflag = 0; tagflag<2; tagflag++) {

    float Maximum = 1.;
    if (FigureOfMerit=="Ratio") 
      Maximum = 1.2*(1-tagflag) + 2.5*tagflag;
    else if (FigureOfMerit=="Significance")
      Maximum = 2.*(1-tagflag) + 3.*tagflag;
    else if (FigureOfMerit=="SignificanceSyst")
      Maximum = 2.*(1-tagflag) + 3.*tagflag;
    else if (FigureOfMerit=="Signal")
      Maximum = 100.;
    else if (FigureOfMerit=="Background")
      Maximum = 10000.;

    TString Option = "histo";
      
    TLegend  *leg = new TLegend(0.13, 0.6, 0.29, 0.85);
    leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
    leg->SetTextColor(1); leg->SetTextSize(0.035);
    leg->SetTextFont(62); 
    leg->SetHeader(TagLegend[tagflag]); //leg->SetMargin(0.2); 

    for (int mp = 0; mp<nmasspoints; mp++) {

      ResultHisto[mp][tagflag]->SetMaximum(Maximum);
      ResultHisto[mp][tagflag]->SetLineWidth(2);
      ResultHisto[mp][tagflag]->SetLineStyle(1);
      ResultHisto[mp][tagflag]->SetMarkerStyle(MassPointMarker[mp]);
      ResultHisto[mp][tagflag]->SetMarkerColor(MassPointColor[mp]);
      ResultHisto[mp][tagflag]->SetLineColor(MassPointColor[mp]);

      ResultHistoISR[mp][tagflag]->SetMaximum(Maximum);
      ResultHistoISR[mp][tagflag]->SetLineWidth(2);
      ResultHistoISR[mp][tagflag]->SetLineStyle(2);
      ResultHistoISR[mp][tagflag]->SetMarkerStyle(MassPointMarker[mp]);
      ResultHistoISR[mp][tagflag]->SetMarkerColor(MassPointColor[mp]);
      ResultHistoISR[mp][tagflag]->SetLineColor(MassPointColor[mp]);
      
      ResultHisto[mp][tagflag]->SetYTitle(TitleY);
      TString Space = "                   ";
      ResultHisto[mp][tagflag]->SetXTitle(Space + "SR1" + Space + Space + "SR2" + Space + Space + "SR3      MT2 [GeV]");
      for (int ib = 1; ib<=nVariableBins*nSignalRegions; ib++) {
	int xbin = 20*(((ib-1)%nVariableBins)+1);
	TString xBin; xBin += xbin; 
	ResultHisto[mp][tagflag]->GetXaxis()->SetBinLabel(ib, "     " + xBin);
      }

      ResultHisto[mp][tagflag]->Draw(Option);

      leg->AddEntry(ResultHisto[mp][tagflag], masspointlegend[mp], "f");

      Option = "histosame";

      if (PlotISR) {

	ResultHistoISR[mp][tagflag]->Draw(Option);
	leg->AddEntry(ResultHistoISR[mp][tagflag], masspointlegend[mp] + "with ISR jet", "f");
    
      }

    }

    leg->Draw();

    TLine *SR1Line = new TLine(140., 0., 140., Maximum);
    SR1Line->SetLineColor(38);
    SR1Line->SetLineStyle(2);
    SR1Line->Draw("same");

    TLine *SR2Line = new TLine(280., 0., 280., Maximum);
    SR2Line->SetLineColor(38);
    SR2Line->SetLineStyle(2);
    SR2Line->Draw("same");

    TString PlotName = FigureOfMerit + "_" + Variable + Channel + TagFlag[tagflag] + LumiFlag + ".png";
    CC->Print("SignalBackground/" + PlotName);

  }

}
