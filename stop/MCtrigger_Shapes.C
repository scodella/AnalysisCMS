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

void MCtrigger_Shapes(TString ProcessName, TString channel, TString cut = "Stop/02_SRs/", TString HistoDirectoryi1 = "../freezing_rootfiles/nominal/Stop/", TString HistoDirectoryi2 = "../MCtrigger_rootfiles/nominal/Stop/", TString HistoDirectoryi3 = "../minitrees/rootfiles/nominal/Stop/",  bool sameFlavour = false) {


  if (info){
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "HistoDirectoryi1  "   << HistoDirectoryi1 << std::endl; 
   std::cout << "HistoDirectoryi2  "   << HistoDirectoryi2 << std::endl; 
   std::cout << "HistoDirectoryi3  "   << HistoDirectoryi3 << std::endl; 
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "channel  "<<  channel << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;
  }
  //TString HistoDirectory = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/";

  TH1F *h_1; TH1F *h_2; TH1F *h_3;

  TFile *_file1 = TFile::Open(HistoDirectoryi2 + ProcessName + ".root"); // MC trigger selection
  TFile *_file2 = TFile::Open(HistoDirectoryi1 + ProcessName + ".root"); // Efficiency from Data
  TFile *_file3 = TFile::Open(HistoDirectoryi3 + ProcessName + ".root"); // No trigger weight

  if (sameFlavour){
 
    std::cout << "sameFlavour" << std::endl;
    
    h_1  = (TH1F*) _file1->Get(cut+"h_MT2ll_ee"); TH1F *h_1_mm = (TH1F*) _file1->Get(cut+"h_MT2ll_mm");
    h_2  = (TH1F*) _file2->Get(cut+"h_MT2ll_ee"); TH1F *h_2_mm = (TH1F*) _file2->Get(cut+"h_MT2ll_mm");
    h_3  = (TH1F*) _file3->Get(cut+"h_MT2ll_ee"); TH1F *h_3_mm = (TH1F*) _file3->Get(cut+"h_MT2ll_mm");

    h_1 -> Add(h_1_mm); 
    h_2 -> Add(h_2_mm);
    h_3 -> Add(h_3_mm);

    channel = "_sameFlavour";  
  }
  else{

    std::cout << "sameFlavour = false" << std::endl;

    h_1 = (TH1F*) _file1->Get(cut+"h_MT2ll_"+channel);
    h_2 = (TH1F*) _file2->Get(cut+"h_MT2ll_"+channel);
    h_3 = (TH1F*) _file3->Get(cut+"h_MT2ll_"+channel);
  }

  gSystem->mkdir("MCTrigger_Eff/", kTRUE); 
  TCanvas* c = new TCanvas ( "MCtrigger", "", 1200,1000); 
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
  //pad1->cd();
  //ad1->SetLogy(1);

  //TPad *pad1 = new TPad();
  //pad1->cd();
  //pad1->SetLogy(); 

  gStyle->SetOptStat(""); 
  TString title = ProcessName.ReplaceAll(".root", "");  

  h_1->SetLineColor(2);
  h_2->SetLineColor(4);

  float min = h_1->GetMinimum();
  h_3->SetMinimum(0.9*min);
  float max = h_2->GetMaximum();
  h_3->SetMaximum(2*max);
  //HCN->SetMaximum(1.2*max);
 
  h_3 -> SetTitle(title + "_" + channel);
  h_3 -> SetXTitle ("MT2ll  (GeV)");

  TLegend* leg;
  //leg = new TLegend(0.1,0.7,0.3,0.9);
  leg = new TLegend(0.75, 0.9, 0.9, 0.7);
  leg -> AddEntry(h_3, "No trigger weight" , "lp");
  leg -> AddEntry(h_1, "MC trigger select." , "lp");
  leg -> AddEntry(h_2, "trigger eff from Data" , "lp");
 
  gPad->SetLogy();
  
  h_3->Draw();
  h_1->Draw("same");
 // h_2->Draw("same");
  leg -> Draw();

  if (verbose) {
    cout << "bin  MC trigger select/trigger eff from Data" <<std::endl;          
    for (int ibin = 1; ibin<=h_1->GetNbinsX(); ibin++)
     {
      cout << ibin << "   " << h_1->GetBinContent(ibin)/h_2->GetBinContent(ibin) << endl;
     }
    cout << "Total MC trigger select." << h_1->Integral()<< endl;
    cout << "Total No trigger weight"  << h_3->Integral()<< endl;
    cout << "Total trigger eff from Data" << h_2->Integral()<< endl;
  }
  c -> Print("MCTrigger_Eff/MCTrigger_"+ title + "_" + channel +".png");
}

 
