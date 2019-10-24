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

void TheoreticalVariations() {

  //TString RootFilesDirectory = "../minitrees/rootfiles3R/";
  TString RootFilesDirectory = "../minitrees/rootfiles/";

 int const nUncertainties = 2;
 TString UncertaintyName[nUncertainties] = {"Q2", "PDF"};
 
 /*int const nProcesses = 11;
 TString ProcessName[nProcesses] = {"02_WZTo3LNu", "03_VZ",  "04_TTTo2L2Nu",
				    "05_ST",       "06_WW",  "07_ZJets",
				    "07_ZJetsHT",  "09_TTW", "10_TTZ",
				    "11_HWW",      "13_VVV"};
 */
 int const nProcesses = 1;
 //TString ProcessName[nProcesses] = {"07_ZJetsHT_DYcorr"};
 TString ProcessName[nProcesses] = {"05_ST"};

 int const nObservables = 4;
 TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen"};

 int const nChannels = 4;
 TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};
 
 for (int un = 0; un<nUncertainties; un++) {
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "up/Stop/", kTRUE);
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "do/Stop/", kTRUE);
 }

 for (int pr = 0; pr<nProcesses; pr++) {

   TFile *InputRootFile = TFile::Open(RootFilesDirectory + "theory/Stop/" + ProcessName[pr] + ".root");

   for (int un = 0; un<nUncertainties; un++) {

     for (int vr = 0; vr<2; vr++) {

       TString Variation = (vr==0) ? "up" : "do";

       TFile *OutputFile = new TFile(RootFilesDirectory + UncertaintyName[un] + Variation + "/Stop/" + ProcessName[pr] + ".root", "recreate");

       for (int sr = 0; sr<ncut; sr++) {

	 if (!scut[sr].Contains("SR") && !scut[sr].Contains("Zveto")) continue;
	 if (!scut[sr].Contains("SR")) continue;
	 if (scut[sr].Contains("ZZ")) continue;

	 //cout << ProcessName[pr] << " " << UncertaintyName[un] << " " << scut[sr] << endl;

	 OutputFile->cd();
	 
	 gDirectory->mkdir(scut[sr]);
	 
	 for (int ob = 0; ob<nObservables; ob++) 
	   for (int ch = 0; ch<nChannels; ch++) {
	     
	     OutputFile->cd(scut[sr]);

	     TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];
	     TH1F *OutputHisto = new TH1F(HistoName, "", 7, 0., 140.);
	     
	     float BinContent[7];
	     for (int ib = 0; ib<7; ib++)
	       if (Variation=="up") BinContent[ib] = -1.;
	       else if (Variation=="do") BinContent[ib] = 99999999.;
	     
	     for (int tv = 0; tv<111; tv++) {

	       if (UncertaintyName[un]=="Q2") {
		 if (tv>=9 || tv==5 || tv==7) continue;
	       } else if (UncertaintyName[un]=="PDF" && tv>0 && tv<9) continue;

	       TString InputHistoName = HistoName + "_TV"; InputHistoName += tv;
	       TH1F *InputHisto = (TH1F*) InputRootFile->Get(scut[sr] + "/" + InputHistoName);

	       for (int ib = 0; ib<7; ib++) {
		 
		 float ThisBinContent = InputHisto->GetBinContent(ib+1);
		 if (Variation=="up" && ThisBinContent>BinContent[ib]) BinContent[ib] = ThisBinContent;
		 else  if (Variation=="do" && ThisBinContent<BinContent[ib]) BinContent[ib] = ThisBinContent;

	       }

	     }

	     for (int ib = 0; ib<7; ib++)
	       OutputHisto->SetBinContent(ib+1, BinContent[ib]);

	     OutputFile->cd(scut[sr]);

	     OutputHisto->Write();

	   }
	 
       }
       
       OutputFile->Close();
       
     }

   }
   
   InputRootFile->Close();

 }

}  
