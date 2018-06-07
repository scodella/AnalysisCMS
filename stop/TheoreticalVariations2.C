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

void TheoreticalVariations2() {

  TString RootFilesDirectory = "../minitrees/rootfiles/";

 int const nUncertainties = 2;
 TString UncertaintyName[nUncertainties] = {"Q2", "PDF"};
 /*
 int const nProcesses = 12;
 TString ProcessName[nProcesses] = {"02_WZTo3LNu", "03_ZZa", "04_TTTo2L2Nu_NoTopPt",
				    "05_ST",       "06_WW", //"07_ZJetsHT",
				    "07_ZJetsHT_DYcorr", "09_TTW", "10_TTZ",
				    "11_HWW",      "13_VVV", "15_VZ", "15_VZ3V"};
 */
 int const nProcesses = 1;
 //TString ProcessName[nProcesses] = {"04_TTTo2L2Nu"};
 //TString ProcessName[nProcesses] = {"07_ZJetsHT"};
 //TString ProcessName[nProcesses] = {"15_VZ"};
 TString ProcessName[nProcesses] = {"03_ZZa"};

 //int const nObservables = 4;
 //TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen"};
 int const nObservables = 14;
 TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen", 
					 "MET",   "mt2LL",    "njet20",   "nbjet",       "njetISR", 
					 "dphijetISR", "ptjetISR", "JetPt", "Lep1Pt",    "Lep2Pt"};

 int const nChannels = 4;
 TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};
 
 for (int un = 0; un<nUncertainties; un++) {
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "up/Stop/", kTRUE);
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "do/Stop/", kTRUE);
 }

 for (int pr = 0; pr<nProcesses; pr++) {

   cout << "Processing " << ProcessName[pr] << endl;

   TFile *InputRootFile = TFile::Open(RootFilesDirectory + "theory/Stop/" + ProcessName[pr] + ".root");

   for (int un = 0; un<nUncertainties; un++) {

     for (int vr = 0; vr<2; vr++) {

       TString Variation = (vr==0) ? "up" : "do";

       cout << "         " << UncertaintyName[un] << Variation << endl;

       TFile *OutputFile = new TFile(RootFilesDirectory + UncertaintyName[un] + Variation + "/Stop/" + ProcessName[pr] + ".root", "recreate");

       for (int sr = 0; sr<ncut; sr++) {

	 if (scut[sr].Contains("Zveto")) continue;
	 if (!scut[sr].Contains("SR") && !scut[sr].Contains("VR")) continue;
	 if (scut[sr].Contains("ZZ")) continue;

	 //cout << ProcessName[pr] << " " << UncertaintyName[un] << " " << scut[sr] << endl;

	 OutputFile->cd();
	 
	 gDirectory->mkdir(scut[sr]);
	 
	 for (int ob = 0; ob<nObservables; ob++) 
	   for (int ch = 0; ch<nChannels; ch++) {

	     //cout << "         " << UncertaintyName[un] << Variation << " " << ObservableName[ob] << ChannelName[ch] << endl;
	     
	     OutputFile->cd(scut[sr]);

	     int nBins = 7; float minX = 0., maxX = 140.;
	     /*
	h_MET_theoreticalvariation               [i][j][tv] = new TH1D("h_MET"                + suffix, "",   80,    0,  800);
	h_njet20_theoreticalvariation            [i][j][tv] = new TH1D("h_njet20"             + suffix, "",    7, -0.5,  6.5);
	h_nbjet_theoreticalvariation             [i][j][tv] = new TH1D("h_nbjet"              + suffix, "",    7, -0.5,  6.5);
	h_njetISR_theoreticalvariation           [i][j][tv] = new TH1D("h_njetISR"            + suffix, "",    7, -0.5,  6.5);
 	h_dphijetISR_theoreticalvariation        [i][j][tv] = new TH1D("h_dphijetISR"         + suffix, "",  100,    0,  3.2);
	h_ptjetISR_theoreticalvariation          [i][j][tv] = new TH1D("h_ptjetISR"           + suffix, "",   80,     0,  800);	
	h_Lep1Pt_theoreticalvariation            [i][j][tv] = new TH1D("h_Lep1Pt"             + suffix, "",  2000,    0, 2000);
	h_Lep2Pt_theoreticalvariation            [i][j][tv] = new TH1D("h_Lep2Pt"             + suffix, "",  2000,    0, 2000);
	h_JetPt_theoreticalvariation             [i][j][tv] = new TH1D("h_JetPt"              + suffix, "",  2000,    0, 2000);
	h_mt2LL_theoreticalvariation             [i][j][tv] = new TH1F("h_mt2LL"              + suffix, "",  2000,    0, 2000);*/
	     if (!ObservableName[ob].Contains("MT2ll")) {
	       if (ObservableName[ob].Contains("MET"))          { nBins =   80; minX =    0.; maxX =  800.; }
	       if (ObservableName[ob].Contains("njet20"))       { nBins =    7; minX =  -0.5; maxX =  6.5; }
	       if (ObservableName[ob].Contains("nbjet"))        { nBins =    7; minX =  -0.5; maxX =  6.5; }
	       if (ObservableName[ob].Contains("njetISR"))      { nBins =    7; minX =  -0.5; maxX =  6.5; }
	       if (ObservableName[ob].Contains("ptjetISR"))     { nBins =   80; minX =    0.; maxX =  800.; }
	       if (ObservableName[ob].Contains("dphijetISR"))   { nBins =  100; minX =     0; maxX =  3.2; }
	       if (ObservableName[ob].Contains("Lep1Pt"))       { nBins = 2000; minX =    0.; maxX =  2000.; }
	       if (ObservableName[ob].Contains("Lep2Pt"))       { nBins = 2000; minX =    0.; maxX =  2000.; }
	       if (ObservableName[ob].Contains("JetPt"))        { nBins = 2000; minX =    0.; maxX =  2000.; }
	       if (ObservableName[ob].Contains("mt2LL"))        { nBins = 2000; minX =    0.; maxX =  2000.; }
	     }

	     TString HistoName = "h_" + ObservableName[ob] + ChannelName[ch];
	     TH1D *OutputHisto = new TH1D(HistoName, "", nBins, minX, maxX);
	     
	     float BinContent[2200];
	     for (int ib = 0; ib<nBins+2; ib++)
	       if (Variation=="up") BinContent[ib] = -1.;
	       else if (Variation=="do") BinContent[ib] = 99999999.;
	     
	     for (int tv = 0; tv<111; tv++) {

	       if (UncertaintyName[un]=="Q2") {
		 if (tv>=9 || tv==5 || tv==7) continue;
	       } else if (UncertaintyName[un]=="PDF" && tv>0 && tv<9) continue;

	       TString InputHistoName = HistoName + "_TV"; InputHistoName += tv;
	       TH1F *InputHisto = (TH1F*) InputRootFile->Get(scut[sr] + "/" + InputHistoName);

	       for (int ib = 0; ib<nBins+2; ib++) {
		 
		 float ThisBinContent = InputHisto->GetBinContent(ib);
		 if (Variation=="up" && ThisBinContent>BinContent[ib]) BinContent[ib] = ThisBinContent;
		 else  if (Variation=="do" && ThisBinContent<BinContent[ib]) BinContent[ib] = ThisBinContent;

	       }

	     }

	     for (int ib = 0; ib<nBins+2; ib++)
	       OutputHisto->SetBinContent(ib, BinContent[ib]);

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


void TheoreticalVariationsSignal() {

  TString RootFilesDirectory = "../minitrees/rootfiles/";

 int const nUncertainties = 2;
 TString UncertaintyName[nUncertainties] = {"Q2", "PDF"};

 int const nObservables = 4;
 TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen"};

 int const nChannels = 4;
 TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};
 
 for (int un = 0; un<nUncertainties; un++) {
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "up/Stop/", kTRUE);
   gSystem->mkdir(RootFilesDirectory + UncertaintyName[un] + "do/Stop/", kTRUE);
 }

 //TString MassPointFileName = "../../PlotsConfigurations/Configurations/T2tt/MassPointList_T2bW.txt";
 //TString MassPointFileName = "../../PlotsConfigurations/Configurations/T2tt/MassPointList_TChiSlep.txt";
 //TString MassPointFileName = "../../PlotsConfigurations/Configurations/T2tt/MassPointList_TChiWW_Ext1.txt";
 TString MassPointFileName = "../../PlotsConfigurations/Configurations/T2tt/MassPointList_TChiWW_mChi1ter.txt";
 ifstream MassPointFile; MassPointFile.open(MassPointFileName);

 int kkk = 0;	
 
 while (MassPointFile) {	

   TString MassPoint;
   MassPointFile >> MassPoint;	

   if (MassPoint=="" || MassPoint.Contains("#")) continue;

   kkk++;

   //if (kkk<30) continue;

   TFile *InputRootFile = TFile::Open(RootFilesDirectory + "theory/Stop/T" + MassPoint + ".root");

   for (int un = 0; un<nUncertainties; un++) {

     for (int vr = 0; vr<2; vr++) {

       TString Variation = (vr==0) ? "up" : "do";

       TFile *OutputFile = new TFile(RootFilesDirectory + UncertaintyName[un] + Variation + "/Stop/T" + MassPoint + ".root", "recreate");

       for (int sr = 0; sr<ncut; sr++) {

	 if (scut[sr].Contains("Zveto")) continue;
	 if (!scut[sr].Contains("SR") && !scut[sr].Contains("VR")) continue;
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
