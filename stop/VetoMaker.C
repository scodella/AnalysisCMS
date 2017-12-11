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
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TDirectory.h"
#include <fstream>
#include <iostream>


void VetoMaker(TString RootFileName = "../minitrees/rootfiles/nominal/Stop/04_TTTo2L2Nu_Test.root",
	       TString metCut       = "02_VR1",
	       TString Variable     = "MT2ll",
	       TString channel      = "_ll",
	       TString kind         = "_fake")
{ 
  TFile* RootFile = new TFile(RootFileName, "update");

  TH1F* histo  = (TH1F*)RootFile->Get("Stop/" + metCut + "_NoTag/h_" + Variable + kind + channel);
  TH1F* histo1 = (TH1F*)RootFile->Get("Stop/" + metCut + "_NoJet/h_" + Variable + kind + channel); 

  RootFile->cd();

  gDirectory->mkdir("Stop/"+ metCut + "_" + "Veto/");

  RootFile->cd("Stop/"+ metCut + "_" + "Veto/");

  TH1F* newHisto = (TH1F*)histo->Clone("newHisto");

  newHisto->Add(histo1);

  newHisto->Write("newHisto", TObject::kOverwrite);
    
  RootFile->Close();
}
