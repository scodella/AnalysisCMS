#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
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
#include "../include/Constants.h"
#include "../include/CutsStop.h"
#include <fstream>
#include <iostream>

//----------------------------------------------
// Global Variables
// ---------------------------------------------

TString FileAddressm = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/WZtoWW/Stop/";
TString FileAddress  = "histo_Minitree_PhilMetCheck/";
//TString FileAddress = "Feb2017_rootfiles_Ghent_ReMiniAODSpring16/nominal/Stop/";
const int nsample = 1;
TString SampleName [nsample] = {"DoubleEG_Run2016B-03Feb2017_ver2-v2.root"};
//TString  _sample   [nsample] = {"T2tt", "TTbar", "ST", "DYJets", "WW", "TTV", "VZ"};



/* Functiion to read a simple Histo from AnalysisCMS output
 *  
 * sample    = only the .root output
 * histoname = h_variable 
 * scut = see include/CutsStop.h
 * channel = ee,mm,em
 * jetBin = " ", 0jet, 1jet, 2jet
 * lumi = 5.8 for runB; 35.867 for all 2016
 *
 *
 * */

 // Just to plot faster 
 void Plot ()

 {
  TFile* WWmimic = new TFile(FileAddressm + "02_WZTo3LNu.root" , "read"); 
  TFile* WW = new TFile(FileAddress + "06_WW.root" , "read"); 
  // variable: MT2ll, dphillMET
  TH1F* h_WWm_var   = (TH1F*) WWmimic -> Get("Stop/02_SRs_NoTag/h_MT2ll_ll");
  TH1F* h_WWm_var_1 = (TH1F*) WWmimic -> Get("Stop/02_SR1_NoJet/h_MT2ll_ll");
  TH1F* h_WWm_var_2 = (TH1F*) WWmimic -> Get("Stop/02_SR2_NoJet/h_MT2ll_ll");
  TH1F* h_WWm_var_3 = (TH1F*) WWmimic -> Get("Stop/02_SR3_NoJet/h_MT2ll_ll");

  h_WWm_var-> Add(h_WWm_var_1);
  h_WWm_var-> Add(h_WWm_var_2);
  h_WWm_var-> Add(h_WWm_var_3);
  float integralm = h_WWm_var->Integral(-1,-1);
  h_WWm_var->Scale(1/integralm);

  TH1F* h_WW_var = (TH1F*) WW -> Get("mt2ll");
//  TH1F* h_WW_dphillmet_1 = (TH1F*) WW -> Get("Stop/02_SR1_NoJet/h_dphillMET_ll");
//  TH1F* h_WW_dphillmet_2 = (TH1F*) WW -> Get("Stop/02_SR2_NoJet/h_dphillMET_ll");
//  TH1F* h_WW_dphillmet_3 = (TH1F*) WW -> Get("Stop/02_SR3_NoJet/h_dphillMET_ll");
//
//  h_WW_dphillmet-> Add(h_WW_dphillmet_1);
//  h_WW_dphillmet-> Add(h_WW_dphillmet_2);
//  h_WW_dphillmet-> Add(h_WW_dphillmet_3);
  float integral = h_WW_var->Integral(-1,-1);
  h_WW_var->Scale(1/integral);

  TCanvas* c_var = new TCanvas (); 
  TLegend* leg= new TLegend(0.75, 0.9, 0.9, 0.7);
  gStyle->SetOptStat(0000);
  h_WW_var  -> SetLineColor(kAzure-9);
  h_WW_var  -> SetLineWidth(2);
  h_WWm_var -> SetLineWidth(2);
  h_WWm_var -> SetLineColor(kOrange-2);
  leg -> AddEntry(h_WW_var, "WW");
  leg -> AddEntry(h_WWm_var, "WZ");
  h_WW_var -> Draw("E");
  h_WWm_var -> Draw("Esame");
  leg -> Draw();
} 

 TH1F* pathReader()
 {

  TString _inputdir ="/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/nominal" ;
  std::vector<TString> _mcfilename;
  std::vector<TString> _systematics;
  
  _mcfilename.push_back("06_WW"); _mcfilename.push_back("05_ST");
  _systematics.push_back("Btag"); _systematics.push_back("Reco");
  
  TFile* myfile1 = new TFile(_inputdir + "/../" + _systematics.at(0) + "do/Stop/" + _mcfilename.at(0) + ".root", "read"); 
  TH1F* myhisto = (TH1F*) myfile1 -> Get("Stop/00_Zveto/h_MT2ll_ll");

  return myhisto; 
 }




/* This function reads a .root from AnalysisCMS and return a selected histogram */

TH1F *reader0 (TString sample, float lumi, TString histoname, TString scut, TString channel, TString jetBin="")
{

  TFile* file0  = new TFile (sample, "read"); 
  TH1F*  histo0 = (TH1F*) file0->Get(scut + "/" + jetBin + histoname + "_" + channel);
  std::cout << "Histo Integral     " <<  histo0 -> Integral(-1,-1) << std::endl; 
  std::cout << "Scale to lumi = " << lumi << std::endl; 
  histo0 -> Scale(lumi);
  std::cout << "Histo Scaled Integral     " <<  histo0 -> Integral(-1,-1) << std::endl;  
  for (Int_t ibin=1; ibin<=histo0->GetNbinsX(); ibin++) {
     std::cout << " Bin number  =  " << ibin << std::endl; 
     std::cout << " Bin content =  " << histo0-> GetBinContent(ibin) << std::endl;
     std::cout << " Bin error   =  " << histo0-> GetBinError (ibin)  << std::endl; 
  }
  return histo0;
}



/* This function reads a TTree in a filename and return the number of entries*/

float reader1 (TString filename)
{ 
  
   TFile* file = TFile::Open(filename);

   TTree* latino = (TTree*)file->Get("latino");

   float N = latino -> GetEntries(); 

   return N; 
}
  
/* This function complete for n iteration the file name and uses recursively reader1 giving the sum */ 

void user (TString filename, int maxIter)
{
  float entries = 0; 
  //you have several parts ex: latino_TTTo2L2Nu_ext1_0001__part8.root 9,10, ......34
  for (int i = 0; i <= maxIter; i++)
   {
    TString file = filename;
    file += std::to_string(i);
    file += ".root";
    std::cout<<file<<std::endl;
    float n = reader1(file);
    entries += n;      
   }
  
  std::cout << entries << std::endl; 
}
