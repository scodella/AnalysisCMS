#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

TH1D *RebinThis(TH1D *inHisto) {

  int const nbins = 11;
  float ptLLbins[nbins+1] = {0., 40., 80., 120., 160., 200., 250., 300., 400., 500., 600., 800.}; 
  //int const nbins = 23;
  //float ptLLbins[nbins+1] = {0, 10., 20., 30., 40., 50., 60, 70., 80., 90., 100., 110., 120., 130., 140., 160., 180., 200., 250., 300., 400., 500., 600., 800.}; 
 
  TH1D *outHisto = new TH1D("outH", "", nbins, ptLLbins);

  int ih = 1;
  float xx = 0., ex2 = 0.;

  for (int ib = 1; ib<=inHisto->GetNbinsX()+1; ib++) {

    if ( outHisto->GetBinLowEdge(ih)<inHisto->GetBinCenter(ib) && 
	 inHisto->GetBinCenter(ib)<outHisto->GetBinLowEdge(ih+1) ) {

      xx += inHisto->GetBinContent(ib);
      ex2 += inHisto->GetBinError(ib)*inHisto->GetBinError(ib);
  
    } else {

      outHisto->SetBinContent(ih, xx);
      outHisto->SetBinError(ih, sqrt(ex2));

      xx = inHisto->GetBinContent(ib);
      ex2 = inHisto->GetBinError(ib)*inHisto->GetBinError(ib);

      ih++;

    }  

  }

  return outHisto;
 
}

void DYcorr(TString hname) {

  TFile *FF = TFile::Open(hname + ".root");

  TH1D* inMC   = (TH1D*) FF->Get("MC");
  TH1D* inDY   = (TH1D*) FF->Get("Z+jets");
  TH1D* inWZ   = (TH1D*) FF->Get("WZ (#rightarrow 3l)");
  TH1D* inZZ   = (TH1D*) FF->Get("ZZ (#rightarrow 2l2#nu)");
  TH1D* inTT   = (TH1D*) FF->Get("t#bar{t}");
  TH1D* inData = (TH1D*) FF->Get("data");
  
  TH1D *MC   = RebinThis(inMC);
  TH1D *DY   = RebinThis(inDY);
  TH1D *WZ   = RebinThis(inWZ);
  TH1D *ZZ   = RebinThis(inZZ);
  TH1D *TT   = RebinThis(inTT);
  TH1D *Data = RebinThis(inData); 
  //Data->Draw();
    
  if (hname.Contains("DPhi")) DY->Rebin(10);

  int lastBin = MC->GetNbinsX();
  float lastbin = MC->GetBinContent(lastBin) + MC->GetBinContent(lastBin+1); MC->SetBinContent(lastBin, lastbin);
  lastbin = DY->GetBinContent(lastBin) + DY->GetBinContent(lastBin+1); DY->SetBinContent(lastBin, lastbin);
  lastbin = WZ->GetBinContent(lastBin) + WZ->GetBinContent(lastBin+1); WZ->SetBinContent(lastBin, lastbin);
  lastbin = ZZ->GetBinContent(lastBin) + ZZ->GetBinContent(lastBin+1); ZZ->SetBinContent(lastBin, lastbin);
  lastbin = TT->GetBinContent(lastBin) + TT->GetBinContent(lastBin+1); TT->SetBinContent(lastBin, lastbin);
  lastbin = Data->GetBinContent(lastBin) + Data->GetBinContent(lastBin+1); Data->SetBinContent(lastBin, lastbin);
  
  //MC->Add(ZZ, -1); ZZ->Scale(1.17); MC->Add(ZZ, +1); 

  MC->Add(DY, -1);
  Data->Add(MC, -1);
  Data->Divide(DY);
  Data->SetXTitle("p_{T}(ll) [GeV]");
  Data->SetYTitle("");
  
  TF1* F1;
  //if (hname.Contains("DPhi")) F1 = new TF1("F1", "[0]*(1.+[1]*x)/(1.+[2]*x)", 0., 800.);
  if (hname.Contains("DPhi")) F1 = F1 = new TF1("F1", "[0]+[1]*TMath::Erf((x+[2])/[3]) + [4]*TMath::Erf((x+[5])/[6])", 0., 800.);
  if (hname.Contains("Pt")) {
    F1 = new TF1("F1", "[0]+[1]*TMath::Erf((x+[2])/[3]) + [4]*TMath::Erf((x+[5])/[6])", 0., 800.);
    //F1->SetParameters(1.02852, -0.0949640,-19.0422, 10.4487, 0.0758834, -56.1146, 41.1653);
  }
  //[0]+[1]*TMath::Erf((_gen_ptll+[2])/[3]) + [4]*TMath::Erf((_gen_ptll+[5])/[6]
  Data->Fit(F1);
  Data->SetLineColor(1);
  Data->SetMarkerColor(1);
  Data->SetMarkerStyle(20);
  Data->Draw("pe");
  cout << F1->GetExpFormula("p") << endl;
  
}


void CompareCorrections() {

TF1 *CCC = new TF1("CCC", "1.97776+0.413621*TMath::Erf((x+-404.316)/143.214)+0.413537*TMath::Erf((x+-404.316)/143.214)", 0., 800.);
TF1 *ZZp = new TF1("ZZp", "2.0805+0.67786*TMath::Erf((x+-479.08)/288.418)+0.67786*TMath::Erf((x+-479.08)/288.418)", 0., 800.); 

 CCC->SetLineColor(1);
 ZZp->SetLineColor(kOrange+3);

 CCC->Draw();
 ZZp->Draw("same");

}
