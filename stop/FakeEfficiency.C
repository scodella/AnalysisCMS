#include "TInterpreter.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>

float         trig_DbleEle;
float         trig_SnglEle;
float         trig_DbleMu;
float         trig_SnglMu;
float         trig_EleMu;

vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_flavour;
vector<float>   *std_vector_lepton_isTightLepton;
vector<float>   *std_vector_lepton_isLooseLepton;

void FakeEfficiency() {

  TString RootSystem = "root://eoscms.cern.ch/";
  TString TreeDirectory = "/eos/cms/store/user/scodella/Stop/LatinoSkims/Feb2017_Run2016E_ReMiniAOD/l2stop__hadd__TrigMakerData/";

  const int nSamples = 5;
  TString SampleName[nSamples] = {"DoubleMuon", "SingleMuon", "DoubleEG", "SingleElectron", "MuonEG"};
  
  const int nPeriods = 8;
  TString PeriodName[nPeriods] = {"Run2016B-03Feb2017_ver2-v2", "Run2016C-03Feb2017-v1", "Run2016D-03Feb2017-v1",
				  "Run2016E-03Feb2017-v1",      "Run2016F-03Feb2017-v1", "Run2016G-03Feb2017-v1", 
				  "Run2016H-03Feb2017_ver2-v1", "Run2016H-03Feb2017_ver3-v1"};

  for (int is = 0; is<nSamples; is++) 
    if (SampleName[is].Contains("Single")) 
      for (int ip = 0; ip<nPeriods; ip++) { 
	if (ip<6 && SampleName[is].Contains("Muon")) continue;
	const int nPtEdges = 7;
	float PtEdge[nPtEdges] = {20., 25., 30., 40., 50., 80., 150.};
	TH2F *TriggerEvents[2][2];
	TriggerEvents[0][0] = new TH2F("SingleEventsOS", "", nPtEdges-1, PtEdge, nPtEdges-1, PtEdge);
	TriggerEvents[1][0] = new TH2F("SingleEventsSS", "", nPtEdges-1, PtEdge, nPtEdges-1, PtEdge);
	TriggerEvents[0][1] = new TH2F("DoubleEventsOS", "", nPtEdges-1, PtEdge, nPtEdges-1, PtEdge);
	TriggerEvents[1][1] = new TH2F("DoubleEventsSS", "", nPtEdges-1, PtEdge, nPtEdges-1, PtEdge);
	
	TString TreeName = RootSystem + TreeDirectory + "latino_" + SampleName[is] + "_" + PeriodName[ip];
	if (TreeName.Contains("Run2016B")) TreeName.ReplaceAll("Run2016E", "Run2016B");
	if (TreeName.Contains("Run2016C")) TreeName.ReplaceAll("Run2016E", "Run2016C");
	if (TreeName.Contains("Run2016D")) TreeName.ReplaceAll("Run2016E", "Run2016D");
	if (TreeName.Contains("Run2016F")) TreeName.ReplaceAll("Run2016E", "Run2016F");
	if (TreeName.Contains("Run2016G")) TreeName.ReplaceAll("Run2016E", "Run2016G");
	if (TreeName.Contains("Run2016H")) TreeName.ReplaceAll("Run2016E", "Run2016H");
	
	TFile *TreeFile = TFile::Open(TreeName + ".root");
	
	TTree *ThisTree = (TTree*) TreeFile->Get("latino");  
	ThisTree->SetBranchAddress("trig_SnglMu",                     &trig_SnglMu);
	ThisTree->SetBranchAddress("trig_SnglEle",                    &trig_SnglEle);
	ThisTree->SetBranchAddress("trig_DbleMu",                     &trig_DbleMu);
	ThisTree->SetBranchAddress("trig_DbleEle",                    &trig_DbleEle);
	ThisTree->SetBranchAddress("trig_EleMu",                      &trig_EleMu);
	ThisTree->SetBranchAddress("std_vector_lepton_pt",            &std_vector_lepton_pt);
	ThisTree->SetBranchAddress("std_vector_lepton_isTightLepton", &std_vector_lepton_isTightLepton);
	ThisTree->SetBranchAddress("std_vector_lepton_isLooseLepton", &std_vector_lepton_isLooseLepton);
	ThisTree->SetBranchAddress("std_vector_lepton_flavour",       &std_vector_lepton_flavour);
	
	Int_t nentries = ThisTree->GetEntries();
	
	cout << "Looping on " << SampleName[is] + "_" + PeriodName[ip] << endl;
	          
	for (Int_t i = 0; i<nentries; i++) {
        
	  ThisTree->GetEntry(i);

	  if (i!=0 && i%(nentries/10)==0) cout << "  Done " << 10*i/(nentries/10) << "%" << endl;

	  bool PassSingleTrigger = false;
	  if (SampleName[is].Contains("SingleMuon")     && trig_SnglMu)  PassSingleTrigger = true;
	  if (SampleName[is].Contains("SingleElectron") && !trig_SnglMu && trig_SnglEle) PassSingleTrigger = true;

	  if (PassSingleTrigger) {

	    float Lep1Pt = -1., Lep2Pt = -1.; int Sign = -1;

	    if (std_vector_lepton_pt->size()>=2)
	      if (std_vector_lepton_isTightLepton->at(0) && std_vector_lepton_isTightLepton->at(1))
		if (std_vector_lepton_pt->at(0)>25. && std_vector_lepton_pt->at(1)>20.) {
		  
		  bool PassThirdLeptonVeto = true;
		  if (std_vector_lepton_pt->size()>2)
		    if (std_vector_lepton_isLooseLepton->at(2)==1)
		       PassThirdLeptonVeto = false;

		  if (PassThirdLeptonVeto) {

		    float MaxPt = PtEdge[nPtEdges-1];
		    Lep1Pt = (std_vector_lepton_pt->at(0)<MaxPt) ? std_vector_lepton_pt->at(0) : MaxPt-1.;
		    Lep2Pt = (std_vector_lepton_pt->at(1)<MaxPt) ? std_vector_lepton_pt->at(1) : MaxPt-1.;
		    Sign = (std_vector_lepton_flavour->at(0)*std_vector_lepton_flavour->at(1)<0) ? 0 : 1;
		    
		  }

		}
	    
	    if (Sign>=0) {

	      TriggerEvents[Sign][0]->Fill(Lep1Pt, Lep2Pt);
		
	      bool PassDoubleLeptonTrigger = false;
	      if (SampleName[is].Contains("SingleMuon")     && (trig_DbleMu  || trig_EleMu)) 
		PassDoubleLeptonTrigger = true;
	      if (SampleName[is].Contains("SingleElectron") && (trig_DbleEle || trig_EleMu)) 
		PassDoubleLeptonTrigger = true;
	      
	      if (PassDoubleLeptonTrigger)
		TriggerEvents[Sign][1]->Fill(Lep1Pt, Lep2Pt);
	      
	    }
	    
	  }
	  
	}

	TString HistoFileName = "./FakeEfficiency/" + SampleName[is] + "_" + PeriodName[ip];
	TFile *HistoFile = new TFile(HistoFileName + ".root", "recreate");

	TriggerEvents[0][0]->Write();
	TriggerEvents[1][0]->Write();
	TriggerEvents[0][1]->Write();
	TriggerEvents[1][1]->Write();

	HistoFile->Close();

      }

}

void CompareFakeEfficiency() {

  TFile *HistoFile = TFile::Open("./FakeEfficiency/SingleLepton.root");
  
  TH2D *SingleEventsOS = (TH2D*) HistoFile->Get("SingleEventsOS");
  TH2D *SingleEventsSS = (TH2D*) HistoFile->Get("SingleEventsSS");
  TH2D *DoubleEventsOS = (TH2D*) HistoFile->Get("DoubleEventsOS");
  TH2D *DoubleEventsSS = (TH2D*) HistoFile->Get("DoubleEventsSS");

  TCanvas *CC = new TCanvas("CL", "", 900, 450); 
  CC->Divide(2, 1);
  
  TPad *PadC1 = (TPad*)CC->GetPad(1);
  TPad *PadC2 = (TPad*)CC->GetPad(2);

  // Run2015B Setting
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  CC->Range(0,0,1,1);
  CC->SetFillColor(10);
  CC->SetBorderMode(0);
  CC->SetBorderSize(2);
  CC->SetTickx(1);
  CC->SetTicky(1);
  CC->SetLeftMargin(0.16);
  CC->SetRightMargin(0.02);
  CC->SetTopMargin(0.05);
  CC->SetBottomMargin(0.13);
  CC->SetFrameFillColor(0);
  CC->SetFrameFillStyle(0);
  CC->SetFrameBorderMode(0);
  
  PadC1->SetFillColor(0);
  PadC1->SetBorderMode(0);
  PadC1->SetBorderSize(2);
  //PadC1->SetGridy();
  //PadC1->SetLogx();
  PadC1->SetTickx(1);
  PadC1->SetTicky(1);
  PadC1->SetLeftMargin(0.16);
  PadC1->SetRightMargin(0.02);
  //PadC1->SetTopMargin(0.05);
  //PadC1->SetBottomMargin(0.31);
  PadC1->SetTopMargin(0.065);
  PadC1->SetBottomMargin(0.13);
  PadC1->SetFrameFillStyle(0);
  PadC1->SetFrameBorderMode(0);
  PadC1->SetFrameFillStyle(0);
  PadC1->SetFrameBorderMode(0);
  PadC1->Draw();
  
  PadC2->SetFillColor(0);
  PadC2->SetBorderMode(0);
  PadC2->SetBorderSize(2);
  //PadC2->SetGridy();
  //PadC2->SetLogx();
  PadC2->SetTickx(1);
  PadC2->SetTicky(1);
  PadC2->SetLeftMargin(0.16);
  PadC2->SetRightMargin(0.02);
  //PadC2->SetTopMargin(0.05);
  //PadC2->SetBottomMargin(0.31);
  PadC2->SetTopMargin(0.065);
  PadC2->SetBottomMargin(0.13);
  PadC2->SetFrameFillStyle(0);
  PadC2->SetFrameBorderMode(0);
  PadC2->SetFrameFillStyle(0);
  PadC2->SetFrameBorderMode(0);
  PadC2->Draw();

  PadC1->cd();

  //DoubleEventsSS->DrawCopy("TEXTCOLZ");

  DoubleEventsOS->Divide(SingleEventsOS);
  DoubleEventsSS->Divide(SingleEventsSS);
  //DoubleEventsSS->Divide(DoubleEventsOS);

  PadC1->cd();

  DoubleEventsOS->GetZaxis()->SetRangeUser(0.8, 1.0);
  DoubleEventsOS->GetXaxis()->SetLabelFont(42);
  DoubleEventsOS->GetYaxis()->SetLabelFont(42);
  DoubleEventsOS->GetXaxis()->SetTitleFont(42);
  DoubleEventsOS->GetYaxis()->SetTitleFont(42);
  DoubleEventsOS->GetYaxis()->SetTitleSize(0.06);
  DoubleEventsOS->GetYaxis()->SetLabelSize(0.05);
  DoubleEventsOS->GetXaxis()->SetTitleSize(0.06);
  DoubleEventsOS->GetXaxis()->SetLabelSize(0.05);
  DoubleEventsOS->GetXaxis()->SetTitleOffset(0.95);
  DoubleEventsOS->GetYaxis()->SetTitleOffset(1.25);
  DoubleEventsOS->SetXTitle("Leading lepton p_{T} [GeV]");
  DoubleEventsOS->SetYTitle("Trailing lepton p_{T} [GeV]");
  DoubleEventsOS->Draw("TEXTCOLZ");

  PadC2->cd();

  DoubleEventsSS->GetZaxis()->SetRangeUser(0.8, 1.0);
  DoubleEventsSS->GetXaxis()->SetLabelFont(42);
  DoubleEventsSS->GetYaxis()->SetLabelFont(42);
  DoubleEventsSS->GetXaxis()->SetTitleFont(42);
  DoubleEventsSS->GetYaxis()->SetTitleFont(42);
  DoubleEventsSS->GetYaxis()->SetTitleSize(0.06);
  DoubleEventsSS->GetYaxis()->SetLabelSize(0.05);
  DoubleEventsSS->GetXaxis()->SetTitleSize(0.06);
  DoubleEventsSS->GetXaxis()->SetLabelSize(0.05);
  DoubleEventsSS->GetXaxis()->SetTitleOffset(0.95);
  DoubleEventsSS->GetYaxis()->SetTitleOffset(1.25);
  DoubleEventsSS->SetXTitle("Leading lepton p_{T} [GeV]");
  DoubleEventsSS->SetYTitle("Trailing lepton p_{T} [GeV]");
  DoubleEventsSS->Draw("TEXTCOLZ");  

}
