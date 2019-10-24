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

bool verbose = false;
bool info    = true;


void Shapes_Builder( TString channel,
                     TString ProcessName1, 
                     TString ProcessName2, 
                     //TString ProcessName = "TTTo2L2Nu__part0", 
                     TString cut = "Stop/02_SRs/",
                     TString cutName = "BaseSelection", 
                     //TString cut = "Stop/02_SRs_NoTag/",
                     //TString cutName = "BaseSelection_NoTag", 
                     //TString HistoDirectoryi1 = "../MCtrigger_rootfiles/nominal/Stop/", 
                     //TString HistoDirectoryi1 = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfilesTrgEmu/nominal/Stop/",  
                     //TString HistoDirectoryi2 = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfilesTrg/nominal/Stop/",  
 
                     TString HistoDirectoryi1 = "../freezing_rootfiles/nominal/Stop/",  //3rd Lepton veto 
                     //TString HistoDirectoryi2 = "../freezing_rootfiles/nominal/Stop/",  //3rd Lepton veto 
                     
                     //TString HistoDirectoryi1 = "../rootfiles_No3LepVeto_Use/nominal/Stop/", //3rd Lepton no veto 
                     TString HistoDirectoryi2 = "../rootfiles_No3LepVeto_Use/nominal/Stop/", //3rd Lepton no veto 

                     //TString HistoDirectoryi3 = "../minitrees/rootfiles/nominal/Stop/",
                     bool sameFlavour = false) 

{
  if (info){
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "ProcessName1   " << ProcessName1 << std::endl;
   std::cout << "ProcessName2   " << ProcessName2 << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "HistoDirectoryi1  "   << HistoDirectoryi1 << std::endl; 
   std::cout << "HistoDirectoryi2  "   << HistoDirectoryi2 << std::endl; 
//   std::cout << "HistoDirectoryi3  "   << HistoDirectoryi3 << std::endl; 
   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "channel  "<<  channel << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;
  }


  //======================================================== 
  // User Configuartion
  //======================================================== 

    //Output folder
    //----------------------------------------------
    TString outputFolder = "3rdVeto_NoVeto/";
    TString outputSuffix;

    //Output suffix
    //----------------------------------------------
    if (ProcessName2 == "") outputSuffix="";
    else outputSuffix= "Yes3rdLepVeto";
   
    //Ratio axis limits
    //----------------------------------------------
    float r_ymin =0.7, r_ymax =1.3;
    //float ymin =0.0, ymax =0.1;

    //Set normalized distributions
    //----------------------------------------------
    bool normalized = false; 

    //Legend and ratio y axis title
    //---------------------------------------------- 
    // For ProcessName1 & ProcessName2 
     // TString ratio_Ytitle = " T2tt / ttbar";
     // TString h1_legendEntry = "T2tt";  
     // TString h2_legendEntry = "ttbar";
      TString legHeader      = "3rd lepton No Veto";
    //----------------------------------------------
    // For ProcessName1 Only  
    TString ratio_Ytitle = " veto / No veto";
    TString h1_legendEntry = "3rd lepton veto";  
    TString h2_legendEntry = "No 3rd lepton veto";

  //======================================================== 
  //======================================================== 

  TH1F *h_1; TH1F *h_2; TH1F *h_3;

  TFile *_file1 = TFile::Open(HistoDirectoryi1 + ProcessName1 + ".root"); // MC trigger selection
  TFile *_file2;
  if (ProcessName2 == "") _file2 = TFile::Open(HistoDirectoryi2 + ProcessName1 + ".root"); // Efficiency from Data
  else _file2 = TFile::Open(HistoDirectoryi2 + ProcessName2 + ".root");
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
  //--------------------------------------------------
  if (normalized) 
   {
     h_1->Scale(1/h_1->Integral());
     h_2->Scale(1/h_2->Integral());
   }
  //--------------------------------------------------
  
  gSystem->cd("/afs/cern.ch/user/b/bchazinq/www/STOP/");
  gSystem->mkdir(outputFolder, kTRUE);
  //gSystem->CopyFile("index.php", "MCTrigger_Eff/.", kTRUE);
  gSystem->cd("-");

  
  TCanvas* canvas = NULL;  
  TPad* pad1 = NULL;
  TPad* pad2 = NULL; 
  
  //TCanvas* c = new TCanvas ( "MCtrigger_"+ ProcessName, "", 1200,1000); 
  TCanvas* c = new TCanvas ( cutName + ProcessName1+ "_"+ ProcessName2, "", 550, 600);
  
  pad1 = new TPad();
  pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  
  pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
 
  TString title;  
  if (ProcessName2 == "") {title = ProcessName1.ReplaceAll(".root", "");}  
  else { title = ProcessName1.ReplaceAll(".root", "") + "Vs" + ProcessName2.ReplaceAll(".root", "");}
  
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
  leg -> AddEntry(h_1, h1_legendEntry, "lp");
  leg -> AddEntry(h_2, h2_legendEntry, "lp");
  //leg -> AddEntry(h_3, "No trigger weight" , "lp");
  leg -> SetTextSize (0.03);
  TString processIN[]  = {"02_WZTo3LNu","03_ZZ","04_TTTo2L2Nu","05_ST","06_WW","07_ZJetsHT","07_ZJetsHT_DYcorr","09_TTW","10_TTZ","11_HWW","13_VVV","14_ZZTo4l","15_VZ", "15_VZ3V","T2tt_mStop-350to400_Sm350_Xm225", "T2tt_mStop-350to400_Sm350_Xm263", "T2tt_mStop-350to400_Sm350_Xm175", "TChiSlep_Xm500_Xm200"};
  TString processOUT[] = {"WZTo3LNu",   "ZZ",       "ttbar",   "tW",   "WW",     "ZJets",   "ZJetsHT_DYcorr","ttbarW","ttbarZ",   "HWW",   "VVV",   "ZZTo4l",   "VZ", "VVV+VZ", "T2tt(350,225)", "T2tt(350,263)", "T2tt(350,175)", "TChiSlep(500,200)"};

  if (ProcessName2 == "")
   {
    for (int i=0; i<18;i++)
      { 
        if (ProcessName1 == processIN[i]) {leg -> SetHeader(processOUT[i] + "  channel " + channel); break;}
        if (i==17 && ProcessName1 != processIN[i]){std::cout<<"WARNING this sample:  "<< ProcessName1 << "is not in the rootfile"<<std::endl;}
      } 
    }
   else
    {
     leg -> SetHeader( legHeader + "  channel " + channel);
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
  float ymin = r_ymin, ymax = r_ymax;
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
  yaxis2->SetTitle(ratio_Ytitle);
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
    cout << "bin h1 / h2    %    " <<std::endl;          
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
  c ->Print("/afs/cern.ch/user/b/bchazinq/www/STOP/" + outputFolder+ "comparation_"+ title + "_" + channel + "_" +cutName + "_" +outputSuffix + ".png");
  gSystem->cd("/afs/cern.ch/user/b/bchazinq/work/CMSSW_8_0_5/src/AnalysisCMS/stop/");
  //std::cout<<  gSystem->pwd() <<std::endl; 
}

 
void Compare_Histos()
{
  TString sampleIN[]  = {"T2tt_mStop-350to400_Sm350_Xm175" };
  //TString sampleIN[]  = {"T2tt_mStop-350to400_Sm350_Xm225","T2tt_mStop-350to400_Sm350_Xm263", "TChiSlep_Xm500_Xm200", "04_TTTo2L2Nu"  };
  //TString sampleIN[]  = {"02_WZTo3LNu", "03_ZZ", "04_TTTo2L2Nu","05_ST","06_WW","07_ZJetsHT_DYcorr","10_TTZ"};
  //"09_TTW","11_HWW","13_VVV","15_VZ","07_ZJetsHT", "14_ZZTo4l","15_VZ3V"
  for (int i=0; i<1;i++)
   {
     Shapes_Builder("ll",sampleIN[i],"");
   }
  //Shapes_Builder("ll","T2tt_mStop-350to400_Sm350_Xm263","04_TTTo2L2Nu");
} 
