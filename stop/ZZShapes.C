#include "../include/CutsStop.h"

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
#include "TLorentzVector.h"
#include "TAxis.h"
#include <fstream>
#include <iostream>


bool verbose = true;

float BB[4] = {1.27037, 1.26317, 1.6111, 1.11252};

TH1D *GetHisto(TFile *filein, TString hname) {

  TH1D *H0 = new TH1D("H0", "", 7, 0., 140.); 

  TH1D *H1 = (TH1D*)filein->Get("Stop/02_SR1_Tag/h_MT2ll_ll");
  TH1D *H2 = (TH1D*)filein->Get("Stop/02_SR2_Tag/h_MT2ll_ll");
  TH1D *H3 = (TH1D*)filein->Get("Stop/02_SR3_Tag/h_MT2ll_ll");
  TH1D *H4 = (TH1D*)filein->Get("Stop/02_SR1_NoTag/h_MT2ll_ll");
  TH1D *H5 = (TH1D*)filein->Get("Stop/02_SR2_NoTag/h_MT2ll_ll");
  TH1D *H6 = (TH1D*)filein->Get("Stop/02_SR3_NoTag/h_MT2ll_ll");
  TH1D *H7 = (TH1D*)filein->Get("Stop/02_SR1_NoJet/h_MT2ll_ll");
  TH1D *H8 = (TH1D*)filein->Get("Stop/02_SR2_NoJet/h_MT2ll_ll");
  TH1D *H9 = (TH1D*)filein->Get("Stop/02_SR3_NoJet/h_MT2ll_ll");

  if (hname.Contains("Tag") || hname.Contains("All")) {
    if (hname.Contains("SR1") || hname.Contains("SRs")) H0->Add(H1);
    if (hname.Contains("SR2") || hname.Contains("SRs")) H0->Add(H2);
    if (hname.Contains("SR3") || hname.Contains("SRs")) H0->Add(H3);
  } 
  if (hname.Contains("NoTag") || hname.Contains("All") || hname.Contains("Veto")) {
    if (hname.Contains("SR1") || hname.Contains("SRs")) H0->Add(H4);
    if (hname.Contains("SR2") || hname.Contains("SRs")) H0->Add(H5);
    if (hname.Contains("SR3") || hname.Contains("SRs")) H0->Add(H6);
  } 
  if (hname.Contains("NoJet") || hname.Contains("All") || hname.Contains("Veto")) {
    if (hname.Contains("SR1") || hname.Contains("SRs")) H0->Add(H7);
    if (hname.Contains("SR2") || hname.Contains("SRs")) H0->Add(H8);
    if (hname.Contains("SR3") || hname.Contains("SRs")) H0->Add(H9);
  }
  //cout << H0->Integral() << endl;
  return H0;

}
  
//TString HistoDirectory = "../minitrees/rootfiles/";
TString HistoDirectory = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/";

void PlotZZShapes(TString ShapeName) {

  TH1::SetDefaultSumw2();

  TFile *_fileM = TFile::Open(HistoDirectory + "nominal/Stop/03_ZZ.root");
//  TFile *_fileD = TFile::Open(HistoDirectory + "kfdPhi/Stop/03_ZZ.root");
//  TFile *_fileP = TFile::Open(HistoDirectory + "kfPt/Stop/03_ZZ.root");
//  TFile *_fileN = TFile::Open(HistoDirectory + "kfNo/Stop/03_ZZ.root");

  TH1D *HM = GetHisto(_fileM, ShapeName);
//  TH1D *HD = GetHisto(_fileD, ShapeName);
//  TH1D *HP = GetHisto(_fileP, ShapeName);
//  TH1D *HN = GetHisto(_fileN, ShapeName);

//  HM->SetMarkerStyle(20);
//  HD->SetLineColor(2);
//  HP->SetLineColor(3);
//  HN->SetLineColor(4);

  TCanvas *CC = new TCanvas("CC", "", 600, 600);

/*  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.33, 1, 1.00);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.33);

  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();

  CC->SetFillColor(10);
  CC->SetBorderMode(0);
  CC->SetBorderSize(2);
  pad1->SetFillColor(10);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad2->SetFillColor(10);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  gStyle->SetOptStat(0);

  pad1->cd();

//  float Integ = HM->Integral(); HM->Scale(1./Integ);
//  Integ = HD->Integral(); HD->Scale(1./Integ);
//  Integ = HP->Integral(); HP->Scale(1./Integ);
//  Integ = HN->Integral(); HN->Scale(1./Integ);
  //HM->Scale(BB[0]/BB[0]);
  //HD->Scale(BB[0]/BB[1]);
  //HP->Scale(BB[0]/BB[2]);
  //HN->Scale(BB[0]/BB[3]);

  TLegend  *leg = new TLegend(0.2, 0.5, 0.5, 0.8);
  leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
  leg->SetTextColor(1); leg->SetTextSize(0.035);
  leg->SetTextFont(62); 
  
  leg->AddEntry(HM, "no k-factors", "f");
//  leg->AddEntry(HN, "ZZ-mass", "lpe");
//  leg->AddEntry(HP, "ZZ-p_{T}", "f");
//  leg->AddEntry(HD, "ZZ-#Delta#phi", "f");

  HM->SetXTitle("M_{T2}^{ll} [GeV]");
  HM->SetYTitle("Fraction of events / 20 GeV");

  HM->GetXaxis()->SetLabelOffset(999.);
  HM->GetXaxis()->SetTitleOffset(999.);
  HM->GetXaxis()->SetLabelSize(0.07);
  HM->GetXaxis()->SetTitleSize(0.07);
  HM->GetYaxis()->SetLabelSize(0.045);
  HM->GetYaxis()->SetTitleSize(0.05);
  HM->SetLineColor(1);
*/
  HM->Draw();
//  HM->DrawCopy("p");
//  HD->DrawCopy("samehisto");
//  HP->DrawCopy("samehisto");
//  HN->DrawCopy("samehisto");

//  leg->Draw();

/*  pad2->cd();

  HD->Divide(HM);
  HP->Divide(HM);
  HN->Divide(HM);

  HD->SetXTitle("M_{T2}^{ll} [GeV]");
  HD->SetYTitle("ZZ-x/ZZ-mass");

  HD->GetXaxis()->SetLabelSize(0.09);
  HD->GetXaxis()->SetTitleSize(0.1);
  HD->GetYaxis()->SetLabelSize(0.09);
  HD->GetYaxis()->SetTitleSize(0.1);

  HD->SetMinimum(0.7);
  HD->SetMaximum(1.3);

  HD->DrawCopy("histoe");
  HP->DrawCopy("samehistoe");
  HN->DrawCopy("samehistoe");
*/
}

void ZZShapes() {

  int const nObservables = 4;
  TString ObservableName[nObservables] = {"MT2ll", "MT2llgen", "MT2llisr", "MT2llisrgen"};
  
  int const nChannels = 4;
  TString ChannelName[nChannels] = {"_ee", "_em", "_mm", "_ll"};
  
  gSystem->mkdir(HistoDirectory + "ZZshapeup/Stop/", kTRUE);
  gSystem->mkdir(HistoDirectory + "ZZshapedo/Stop/", kTRUE);

  TFile *_fileM = TFile::Open(HistoDirectory + "nominal/Stop/03_ZZ.root");
  TFile *_fileD = TFile::Open(HistoDirectory + "kfdPhi/Stop/03_ZZ.root");
  TFile *_fileP = TFile::Open(HistoDirectory + "kfPt/Stop/03_ZZ.root");
  TFile *_fileN = TFile::Open(HistoDirectory + "kfNo/Stop/03_ZZ.root");
  
  for (int vr = 0; vr<2; vr++) {
    
    TString Variation = (vr==0) ? "up" : "do";

    TFile *OutputFile = new TFile(HistoDirectory + "ZZshape" + Variation + "/Stop/03_ZZ.root", "recreate");
    
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
	  TH1F *OutputHisto;
	  if (vr==0) {
	    OutputHisto = (TH1F*) _fileP->Get(scut[sr] + "/" + HistoName);
	    OutputHisto->Scale(BB[0]/BB[2]);
	  } else {
	    OutputHisto = (TH1F*) _fileN->Get(scut[sr] + "/" + HistoName);
	    OutputHisto->Scale(BB[0]/BB[3]);
	  } 

	  OutputHisto->Write();
	  
	}
         
    }
       
    OutputFile->Close();
       
  }

}


//------------------------------------------------------------------------------
//  kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
//------------------------------------------------------------------------------
float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {        
        k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
        k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
        k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
        k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
        k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
        k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
        k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
        k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
        k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
        k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
        k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
        k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
        k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
        k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
        k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
        k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
        k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
        k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
        k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
        k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
        k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
        k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
        k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
        k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
        k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
        k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
        k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
        k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
        k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
        k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);
    }

    if (finalState==2) {
       k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
       k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
       k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
       k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
       k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
       k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
       k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
       k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
       k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
       k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
       k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
       k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
       k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
       k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
       k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
       k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
       k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
       k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
       k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
       k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
       k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
       k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
       k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
       k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
       k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
       k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
       k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
       k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
       k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
       k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

//------------------------------------------------------------------------------
//  kfactor_qqZZ_qcd_M   (float GENmassZZ, int finalState);
//------------------------------------------------------------------------------
float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.13223956942*(abs(GENmassZZ)>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}


//------------------------------------------------------------------------------
//  kfactor_qqZZ_qcd_Pt  (float GENpTZZ,   int finalState);
//------------------------------------------------------------------------------
float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.57222952415*(abs(GENpTZZ)>100.0);
    }

    if (finalState==2) {
        k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.58485762205*(abs(GENpTZZ)>100.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}



TH1F *ApplyKFactor(TH1F *INHISTO, int Type, int State) {

  int nBins = INHISTO->GetNbinsX();
  float minX = INHISTO->GetBinLowEdge(1);
  float maxX = INHISTO->GetBinLowEdge(nBins+1);
  TH1F *OUTHISTO = new TH1F("out", "", nBins, minX, maxX);

  for (int ib = 1; ib<=nBins; ib++) {

    float xx = OUTHISTO->GetBinCenter(ib);
    float yy = INHISTO->GetBinContent(ib);

    if (Type==0) yy *= kfactor_qqZZ_qcd_dPhi(xx, State);
    if (Type==1) yy *= kfactor_qqZZ_qcd_M   (xx, State);
    if (Type==2) yy *= kfactor_qqZZ_qcd_Pt  (xx, State);

    OUTHISTO->SetBinContent(ib, yy);

  }

  return OUTHISTO;

}

void TestKFactors() {

  TFile *_file0 = TFile::Open("../minitrees/multilepton/Stop/ZZTo4L__part0.root");
  //TFile *_file0 = TFile::Open("../minitrees/nominalTrg/Stop/ZZTo4L__part0.root");
  TTree *latino = (TTree*)_file0->Get("latino");
  
  TH1F *M = new TH1F("M", "", 40, 0., 1000.);
  latino->Project("M", "ZZmass", "nlepton>=4 && ntightlepton>=3 && lep1pt>20. && lep2pt>20. && lep3pt>20. && lep4pt>20. && (fabs(lep1id)!=fabs(lep2id) || fabs(lep1id)!=fabs(lep3id) || fabs(lep1id)!=fabs(lep4id))");
  //latino->Project("M", "ZZmass");
  TH1F *M2 = ApplyKFactor(M, 1, 2);
  cout << M2->Integral()/M->Integral() << endl;
  TString M2Legend = "ZZ-mass";
  M2->SetXTitle("M(ZZ) [GeV]");
  M->SetYTitle("Arbitrary units / 25 GeV");
  /*
  TH1F *M = new TH1F("M", "", 80, 0., 400.);
  latino->Project("M", "ZZpt", "nlepton>=4 && ntightlepton>=3 && lep1pt>20. && lep2pt>20. && lep3pt>20. && lep4pt>20. && (fabs(lep1id)!=fabs(lep2id) || fabs(lep1id)!=fabs(lep3id) || fabs(lep1id)!=fabs(lep4id))");
  //latino->Project("M", "ZZpt");
  TH1F *M2 = ApplyKFactor(M, 2, 2);
  cout << M2->Integral()/M->Integral() << endl;
  TString M2Legend = "ZZ-p_{T}";
  M2->SetXTitle("p_{T}(ZZ) [GeV]");
  M->SetYTitle("Arbitrary units / 5 GeV");
  
  TH1F *M = new TH1F("M", "", 32, 0., 3.2);
  latino->Project("M", "fabs(ZZdphi)", "nlepton>=4 && ntightlepton>=3 && lep1pt>20. && lep2pt>20. && lep3pt>20. && lep4pt>20. && (fabs(lep1id)!=fabs(lep2id) || fabs(lep1id)!=fabs(lep3id) || fabs(lep1id)!=fabs(lep4id))");
  //latino->Project("M", "ZZdphi");
  TH1F *M2 = ApplyKFactor(M, 0, 2);
  cout << M2->Integral()/M->Integral() << endl;
  TString M2Legend = "ZZ-#Delta#phi";
  M2->SetXTitle("#Delta#phi(ZZ) [GeV]");
  M->SetYTitle("Arbitrary units / 0.1");
  */
  TCanvas *CC = new TCanvas("CC", "", 600, 600);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.33, 1, 1.00);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.33);

  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();

  CC->SetFillColor(10);
  CC->SetBorderMode(0);
  CC->SetBorderSize(2);
  pad1->SetFillColor(10);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad2->SetFillColor(10);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  gStyle->SetOptStat(0);

  pad1->cd();

  M->GetXaxis()->SetLabelOffset(999.);
  M->GetXaxis()->SetTitleOffset(999.);
  M->GetXaxis()->SetTitleSize(0.07);
  M->GetYaxis()->SetLabelSize(0.045);
  M->GetYaxis()->SetTitleSize(0.05);
  M2->SetLineColor(2);

  float maxY = M2->GetMaximum();
  M->SetMaximum(maxY*1.1);
  M->DrawCopy("histo");
  M2->DrawCopy("histosame");
  
   TLegend  *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg->SetFillColor(kWhite); leg->SetBorderSize(0.);
  leg->SetTextColor(1); leg->SetTextSize(0.035);
  leg->SetTextFont(62); 
  
  leg->AddEntry(M, "no k-factors", "f");
  leg->AddEntry(M2, M2Legend, "f");
  
  leg->Draw();
 
  pad2->cd();

  M2->Divide(M);
  M2->SetYTitle("k-factors");

  //M2->GetXaxis()->SetLabelOffset(0.9);
  //M2->GetXaxis()->SetTitleOffset(0.9);
  M2->GetXaxis()->SetLabelSize(0.09);
  M2->GetXaxis()->SetTitleSize(0.1);
  M2->GetYaxis()->SetLabelSize(0.09);
  M2->GetYaxis()->SetTitleSize(0.1);
  M2->GetYaxis()->SetTitleOffset(0.5);

  M2->SetMinimum(0.6);
  M2->SetMaximum(1.6);

  M2->DrawCopy("histo");

}
