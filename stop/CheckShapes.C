#include"TString.h"
#include"TH1F.h"
#include"TFile.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPad.h"

bool verbose = true;
bool info    = true;

void CheckShapes(TString ProcessName, TString UncertaintyName, TString ShapeName, TString channel, TString HistoDirectory = "../minitrees/rootfiles2R/", bool sameFlavour = false) {


  if (info){
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "ProcessName       " << ProcessName << std::endl; 
   std::cout << "UncertaintyName   " << UncertaintyName << std::endl; 
   std::cout << "ShapeName + channel  " << ShapeName+"_"+channel << std::endl; 
   std::cout << "HistoDirectory  "   << HistoDirectory << std::endl; 
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << HistoDirectory + UncertaintyName + "up/Stop/" + ProcessName + ".root" << std::endl;
   std::cout << HistoDirectory +  "nominal/Stop/" + ProcessName + ".root" << std::endl;
   std::cout << HistoDirectory + UncertaintyName + "do/Stop/" + ProcessName + ".root" << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;
  }
  //TString HistoDirectory = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/";

  TH1F *HUP; TH1F *HCN; TH1F *HDO;

  TFile *_fileu = TFile::Open(HistoDirectory + UncertaintyName + "up/Stop/" + ProcessName + ".root");
  TFile *_filec = TFile::Open(HistoDirectory + "nominal/Stop/"              + ProcessName + ".root");
  TFile *_filed = TFile::Open(HistoDirectory + UncertaintyName + "do/Stop/" + ProcessName + ".root");

  if (sameFlavour){
 
    std::cout << "sameFlavour" << std::endl;
    
    HUP = (TH1F*) _fileu->Get(ShapeName+"_ee"); TH1F *HUP_mm = (TH1F*) _fileu->Get(ShapeName+"_mm");
    HCN = (TH1F*) _filec->Get(ShapeName+"_ee"); TH1F *HCN_mm = (TH1F*) _filec->Get(ShapeName+"_mm");
    HDO = (TH1F*) _filed->Get(ShapeName+"_ee"); TH1F *HDO_mm = (TH1F*) _filed->Get(ShapeName+"_mm");

    HUP -> Add(HUP_mm); 
    HCN -> Add(HCN_mm);
    HDO -> Add(HDO_mm);

    channel = "_sameFlavour";  
  }
  else{

    std::cout << "sameFlavour = false" << std::endl;

    HUP = (TH1F*) _fileu->Get(ShapeName+"_"+channel);
    HCN = (TH1F*) _filec->Get(ShapeName+"_"+channel);
    HDO = (TH1F*) _filed->Get(ShapeName+"_"+channel);
  }

  gSystem->mkdir("Pileup_Uncert/", kTRUE); 
  TCanvas* c = new TCanvas ( "PU_uncertainty", "", 1200,1000); 
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
  //pad1->cd();
  //ad1->SetLogy(1);

  //TPad *pad1 = new TPad();
  //pad1->cd();
  //pad1->SetLogy(); 

  gStyle->SetOptStat(""); 
  TString title = ProcessName.ReplaceAll(".root", "");  

  HUP->SetLineColor(2);
  HDO->SetLineColor(4);

  float min = HDO->GetMinimum();
  HCN->SetMinimum(0.9*min);
  float max = HCN->GetMaximum();
  HCN->SetMaximum(2*max);
  //HCN->SetMaximum(1.2*max);
 
  HCN -> SetTitle(title + "_" + channel);
  HCN -> SetXTitle ("MT2ll  (GeV)");

  TLegend* leg;
  //leg = new TLegend(0.1,0.7,0.3,0.9);
  leg = new TLegend(0.75, 0.9, 0.9, 0.7);
  leg -> AddEntry(HCN, "nominal" , "lp");
  leg -> AddEntry(HUP, "Pileup Up" , "lp");
  leg -> AddEntry(HDO, "Pileup Down" , "lp");
 
  gPad->SetLogy();
  
  HCN->Draw();
  HUP->Draw("same");
  HDO->Draw("same");
  leg -> Draw();

  if (verbose) {
    cout << "bin  HUP/HCN       %      HDO/HCN       %" <<std::endl;          
    for (int ibin = 1; ibin<=HCN->GetNbinsX(); ibin++)
     {
      cout << ibin << "   " << HUP->GetBinContent(ibin)/HCN->GetBinContent(ibin) << "   " << 100*(HCN->GetBinContent(ibin)-HUP->GetBinContent(ibin))/HUP->GetBinContent(ibin) << "   "<< HDO->GetBinContent(ibin)/HCN->GetBinContent(ibin) << " " << 100*(HCN->GetBinContent(ibin)-HDO->GetBinContent(ibin))/HDO->GetBinContent(ibin)<< endl;
     }
    cout << "Total   " << HUP->Integral()/HCN->Integral() << "   " << 100*(HCN->Integral()-HUP->Integral())/HCN->Integral() << "   " << HDO->Integral()/HCN->Integral() << "   " << 100*(HCN->Integral()-HDO->Integral())/HCN->Integral()<< endl;
  }
  c -> Print("Pileup_Uncert/PU_uncertainty"+ title + "_" + channel +".png");
}

 
