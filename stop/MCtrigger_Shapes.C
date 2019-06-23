#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMath.h"
#include "TFrame.h"
#include "TAttMarker.h"
#include "TLine.h"

bool verbose = true;
bool info    = true;


void MCtrigger_Shapes_Builder( TString channel,
                               TString ProcessName, 
                               //TString ProcessName = "TTTo2L2Nu__part0", 
                               TString cut = "Stop/02_SRs/", 
                               //TString HistoDirectoryi1 = "../MCtrigger_rootfiles/nominal/Stop/", 
                               TString HistoDirectoryi1 = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfilesTrgEmu/nominal/Stop/",  
                               //TString HistoDirectoryi2 = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfilesTrg/nominal/Stop/",  
                               TString HistoDirectoryi2 = "../freezing_rootfiles/nominal/Stop/",  
                               //TString HistoDirectoryi3 = "../minitrees/rootfiles/nominal/Stop/",
                               bool sameFlavour = false) 

{
  if (info){
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "HistoDirectoryi1  "   << HistoDirectoryi1 << std::endl; 
   std::cout << "HistoDirectoryi2  "   << HistoDirectoryi2 << std::endl; 
//   std::cout << "HistoDirectoryi3  "   << HistoDirectoryi3 << std::endl; 
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "channel  "<<  channel << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;
  }

  TH1F *h_1; TH1F *h_2; TH1F *h_3;

  TFile *_file1 = TFile::Open(HistoDirectoryi1 + ProcessName + ".root"); // MC trigger selection
  TFile *_file2 = TFile::Open(HistoDirectoryi2 + ProcessName + ".root"); // Efficiency from Data
//  TFile *_file3 = TFile::Open(HistoDirectoryi3 + ProcessName + ".root"); // No trigger weight

  if (sameFlavour){
 
    std::cout << "sameFlavour" << std::endl;
    
    h_1  = (TH1F*) _file1->Get(cut+"h_MT2ll_ee"); TH1F *h_1_mm = (TH1F*) _file1->Get(cut+"h_MT2ll_mm");
    h_2  = (TH1F*) _file2->Get(cut+"h_MT2ll_ee"); TH1F *h_2_mm = (TH1F*) _file2->Get(cut+"h_MT2ll_mm");
//  h_3  = (TH1F*) _file3->Get(cut+"h_MT2ll_ee"); TH1F *h_3_mm = (TH1F*) _file3->Get(cut+"h_MT2ll_mm");

    h_1 -> Add(h_1_mm); 
    h_2 -> Add(h_2_mm);
//    h_3 -> Add(h_3_mm);

    channel = "_sameFlavour";  
  }
  else{

    std::cout << "sameFlavour = false" << std::endl;

    h_1 = (TH1F*) _file1->Get(cut+"h_MT2ll_"+channel);
    h_2 = (TH1F*) _file2->Get(cut+"h_MT2ll_"+channel);
//    h_3 = (TH1F*) _file3->Get(cut+"h_MT2ll_"+channel);
  }

  // Normalization to unity
  h_1->Scale(1/h_1->Integral());
  h_2->Scale(1/h_2->Integral());

  gSystem->cd("/afs/cern.ch/user/b/bchazinq/www/public/");
  gSystem->mkdir("MCTrigger_Eff_Shape/", kTRUE);
  //gSystem->CopyFile("index.php", "MCTrigger_Eff/.", kTRUE);
  gSystem->cd("-");

  
  TCanvas* canvas = NULL;  
  TPad* pad1 = NULL;
  TPad* pad2 = NULL; 
  
  //TCanvas* c = new TCanvas ( "MCtrigger_"+ ProcessName, "", 1200,1000); 
  TCanvas* c = new TCanvas ( "MCtrigger_"+ ProcessName, "", 550, 600);
  
  pad1 = new TPad();
  pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  
  pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  

  TString title = ProcessName.ReplaceAll(".root", "");  

  
  pad1->cd();
  gStyle->SetOptStat(""); 
  pad1->SetLogy(); 
  //gPad->SetLogy();


  float min = h_1->GetMinimum();
  h_1->SetMinimum(0.9*min);
  float max = h_2->GetMaximum();
  h_1->SetMaximum(2*max);
  //HCN->SetMaximum(1.2*max);

  //axis title
  //----------- 
  TAxis* xaxis = (TAxis*)h_1->GetXaxis();
  TAxis* yaxis = (TAxis*)h_1->GetYaxis();
  xaxis->SetTitleSize(.10);
  xaxis->SetTitle("MT2ll  (GeV)");
  xaxis->SetLabelSize(.08);
  yaxis->SetTitle("events/20 GeV");
  yaxis->SetTitleSize(.05);
  yaxis->SetTitleOffset(1);
  yaxis->CenterTitle(); 

  // Legend
  //----------------------------------------------------------------------------
  TLegend* leg;
  //leg = new TLegend(0.1,0.7,0.3,0.9);
  leg = new TLegend(0.65, 0.9, 0.9, 0.7);
  leg -> AddEntry(h_1, "MC trigger emulation" , "lp");
  leg -> AddEntry(h_2, "trigger eff from Data" , "lp");
  //leg -> AddEntry(h_3, "No trigger weight" , "lp");
  leg -> SetTextSize (0.03);
  TString processIN[]  = {"02_WZTo3LNu","03_ZZ","04_TTTo2L2Nu","05_ST","06_WW","07_ZJetsHT","07_ZJetsHT_DYcorr","09_TTW","10_TTZ","11_HWW","13_VVV","14_ZZTo4l","15_VZ", "15_VZ3V"};
  TString processOUT[] = {"WZTo3LNu",   "ZZ",       "ttbar",   "tW",   "WW",     "ZJets",   "ZJetsHT_DYcorr","ttbarW","ttbarZ",   "HWW",   "VVV",   "ZZTo4l",   "VZ", "VVV+VZ"};
  for (int i=0; i<14;i++)
    { 
      if (ProcessName == processIN[i]) {leg -> SetHeader(processOUT[i] + "  channel " + channel); break;}
      if (i==13 && ProcessName != processIN[i]){std::cout<<"WARNING this sample:  "<< ProcessName << "is not in the rootfile"<<std::endl;}
    } 
    
  h_1->SetMarkerStyle(kPlus);
  h_2->SetMarkerStyle(kPlus);
  h_1->SetMarkerSize(4);
  h_2->SetMarkerSize(4);
  h_1->SetMarkerColor(2);
  h_2->SetMarkerColor(4);
  
  // Draw
  //-----------------------------------------------------------------------------
  h_1->Draw("hist p0");
  h_2->Draw("hist p0, same");
  //h_3->Draw("same");
  leg -> Draw();
  
  //------------------------------------------------------------------------------- 
  // Ratio
  //-------------------------------------------------------------------------------
  pad2->cd();
  TH1F* ratio = (TH1F*)h_1->Clone("ratio");
  float ymin =0.95, ymax =1.05;
  for (Int_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {
    float ratioVal = -999; //float ratioErr = 999; 
    ratioVal = h_1->GetBinContent(ibin)/h_2->GetBinContent(ibin);
    //ratioErr = (ratioVal * ratioVal)*( TMath::Power(h_1->GetBinError(ibin)/ h_1->GetBinContent(ibin),2) + TMath::Power(h_2->GetBinError(ibin)/ h_2->GetBinContent(ibin),2));
    ratio -> SetBinContent(ibin, ratioVal);
    //ratio->SetBinError (ibin, ratioErr); 
    //if (ratioVal+ratioErr>ymax) ymax = (ratioVal+ratioErr) + (ratioVal+ratioErr)*0.01;
    //if (ratioVal>0. && ratioVal-ratioErr<ymin) ymin = (ratioVal-ratioErr) - (ratioVal+ratioErr)*0.01 ;  
    //float min = ratioVal - 0.02*(h_2->GetBinContent(ibin)); if (abs(min) < abs(ymin)) ymin=min;
    //float max = ratioVal + 0.02*(h_2->GetBinContent(ibin)); if (abs(max) > abs(ymax)) ymin=max;
  }
  ratio->SetMarkerColor(1);
  ratio->SetMarkerSize(1);
  ratio->SetMarkerStyle(kFullCircle);
  ratio->GetYaxis()->SetRangeUser(ymin, ymax);
 
  //axix title
  //------------------------------------------------------------------------------
  TAxis* yaxis2 = (TAxis*)ratio->GetYaxis();
  yaxis2->SetLabelSize(.08);
  yaxis2->SetTitle("McTrigg / EffData");
  yaxis2->SetTitleSize(.055);
  yaxis2->SetTitleOffset(1);
  yaxis->CenterTitle();
  //draw
  //-----------------------------------------------------------------------------
  ratio->Draw("hist p0");
  
  c->Update();
  pad2->Update();
  TLine *line=new TLine(0,1.0,140,1.);
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(kDotted);
  line->Draw("same");
  c->Modified();
  c->Update();
  
  if (verbose) {
    cout << "bin  MC trigger select/trigger eff from Data    %    " <<std::endl;          
    for (int ibin = 1; ibin<=h_1->GetNbinsX(); ibin++)
     {
      cout << ibin << "   " << h_1->GetBinContent(ibin)/h_2->GetBinContent(ibin) <<"                            " << Form("%.2f", (h_1->GetBinContent(ibin)/h_2->GetBinContent(ibin))*100) <<endl;
    //  cout << "   MCtrigger_rootfiles   " << h_1->GetBinContent(ibin) << "   freezing_rootfiles  " << h_2->GetBinContent(ibin) <<endl; 
     }
    cout << "                                           "<< endl; 
    cout << "                                           "<< endl; 
    cout << "Total MC trigger select." << h_1->Integral()<< endl;
    cout << "Total trigger eff from Data" << h_2->Integral()<< endl;
    cout << "                                           "<< endl;
    cout << "Total_MC_trigger_select./Total_trigger_eff_from_Data    %"<< endl;
    cout << h_1->Integral()/h_2->Integral() << "                     " << Form("%.2f", (h_1->Integral()/h_2->Integral())*100) << endl;    
//    cout << "Total No trigger weight"  << h_3->Integral()<< endl;
  }
  //c ->SaveAs("/afs/cern.ch/user/b/bchazinq/www/MCTrigger_Eff/MCTrigger_"+ title + "_" + channel +".png");
  c ->Print("/afs/cern.ch/user/b/bchazinq/www/public/MCTrigger_Eff_Shape/MCTrigger_"+ title + "_" + channel +".png");
  gSystem->cd("/afs/cern.ch/user/b/bchazinq/work/CMSSW_8_0_5/src/AnalysisCMS/stop/");
  //std::cout<<  gSystem->pwd() <<std::endl; 
}

 
void MCtrigger_Shapes()
{
  TString sampleIN[]  = {"02_WZTo3LNu", "03_ZZ", "04_TTTo2L2Nu","05_ST","06_WW","07_ZJetsHT_DYcorr","10_TTZ"};
  //"09_TTW","11_HWW","13_VVV","15_VZ","07_ZJetsHT", "14_ZZTo4l","15_VZ3V"
  for (int i=0; i<7;i++)
   {
     MCtrigger_Shapes_Builder("ll",sampleIN[i]);
   }
} 
