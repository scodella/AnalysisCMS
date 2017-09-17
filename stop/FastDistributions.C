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

void FastDistributions(TString Variable, TString Cut) {

  TString DirName = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/nominal/Stop/";

  TChain DY("latino");
  //DY.Add(DirName + "DYJetsToLL_M-10to50-LO.root");
  //DY.Add(DirName + "DYJetsToLL_M-5to50_HT-*root");
  //DY.Add(DirName + "DYJetsToLL_M-50-LO-ext1__part*.root");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-70to100__part*");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-100to200_ext1__part*");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-200to400_ext1__part*");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-400to600.root");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-600to800__part*");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-800to1200.root");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-1200to2500.root");
  //DY.Add(DirName + "DYJetsToLL_M-50_HT-2500toInf.root"); 
  DY.Add(DirName + "DYJetsToLL_M-50__part*.root"); 

  TChain ZZ("latino");
  //ZZ.Add(DirName + "02_WZ*.root");
  //ZZ.Add(DirName + "03_VZ.root");
  ZZ.Add(DirName + "TTTo2L2Nu_*.root");

  TCanvas *CC = new TCanvas("CC","",800,400);
  CC->Divide(2, 1);
  
  CC->cd(1);
  DY.Draw(Variable, Cut);
  
  CC->cd(2);
  ZZ.Draw(Variable, Cut);


}
