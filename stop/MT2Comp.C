void MT2Comp(TString Selection, TString Region, bool NoNorm = false) {

  TFile *F1 = TFile::Open(Selection + "_" + Region + ".root");
  TFile *F2 = TFile::Open(Selection + "_" + Region + "_MetCorr.root");

  TH1D* MC1   = (TH1D*) F1->Get("MC");
  TH1D* MC2   = (TH1D*) F2->Get("MC");
  
  TCanvas *canvas = new TCanvas("cname", "", 550, 720);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  
  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  
  pad1->cd();
  pad1->SetLogy();

  if (!NoNorm) {
    float MC1integral = MC1->Integral(); MC1->Scale(1./MC1integral);
    float MC2integral = MC2->Integral(); MC2->Scale(1./MC2integral);
  }

  MC1-> SetLineColor(1); MC2-> SetLineColor(2);
  MC1-> SetLineStyle(1); MC2-> SetLineStyle(1); 
  MC1-> SetLineWidth(2); MC2-> SetLineWidth(2); 
  MC1-> SetFillStyle(0); MC2-> SetFillStyle(0); 
  MC1-> SetFillColor(0); MC2-> SetFillStyle(0); 

  Float_t padw = gPad->XtoPixel(gPad->GetX2());
  Float_t padh = gPad->YtoPixel(gPad->GetY1());
  Float_t size = (padw < padh) ? padw : padh;
  size = 20. / size;  // Like this label size is always 20 pixels

  float hmax = MC1->GetMaximum();
  hmax *= 3.5;

  MC1->GetYaxis()->SetTitleSize(size);
  MC1->SetYTitle(" events / 20 GeV ");
  MC1->GetYaxis()->SetLabelSize(size);
  MC2->GetYaxis()->SetTitleOffset(0.5);

  MC1->DrawCopy("hist");
  MC2->DrawCopy("histsame");

  TLegend* legend = new TLegend(0.3, 0.74, 0.7, 0.84);
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  ( size);

  TString Header = Selection + Region;
  Header.ReplaceAll("Zpeak", "|M_{ee,#mu#mu}-M_{Z}|<15 GeV  " );
  Header.ReplaceAll("_Tag", "  Tag");
  Header.ReplaceAll("_Veto", "  Veto");
  Header.ReplaceAll("SRs", "p_{T}^{miss}>140 GeV"); 
  Header.ReplaceAll("SR1", "140<p_{T}^{miss}<200 GeV"); 
  Header.ReplaceAll("SR2", "200<p_{T}^{miss}<300 GeV"); 
  Header.ReplaceAll("SR3", "p_{T}^{miss}>300 GeV"); 

  TLatex* tl = new TLatex(0.15, 0.87, Header);
  tl->SetNDC();
  tl->SetTextSize(size);
  tl->SetTextFont(42);
  tl->Draw();

  legend->AddEntry(MC1, "MC", "l");
  legend->AddEntry(MC2, "MC corrected", "l");

  legend->Draw();

  pad2->cd();

  padw = gPad->XtoPixel(gPad->GetX2());
  padh = gPad->YtoPixel(gPad->GetY1());
  size = (padw < padh) ? padw : padh;
  size = 20. / size;  // Like this label size is always 20 pixels
  
  MC2->Divide(MC1);

  float ydiff = 0.;
  for (int ib = 1; ib<=MC2->GetNbinsX(); ib++) 
    if (fabs(MC2->GetBinContent(ib)-1.)>ydiff) ydiff = fabs(MC2->GetBinContent(ib)-1.);
  ydiff *= 20.; ydiff += 0.5; ydiff *= 2.; int iydiff = int(ydiff); ydiff = float(iydiff)/20.;

  MC2->SetMinimum(1. - ydiff);
  MC2->SetMaximum(1. + ydiff);
  MC2->SetMarkerStyle(20);
  MC2->SetMarkerSize(1.);
  MC2->SetMarkerColor(1);

  MC2->GetXaxis()->SetTitleSize(size);
  MC2->GetYaxis()->SetTitleSize(size);
  MC2->GetXaxis()->SetLabelSize(size);
  MC2->GetYaxis()->SetLabelSize(size);
  MC2->GetXaxis()->SetTitleOffset(1.);
  MC2->GetYaxis()->SetTitleOffset(0.5);
  MC2->SetYTitle(" MC corr. / MC");
  MC2->SetXTitle(" M_{T2}(ll) [GeV]");

  MC2->Draw("phist");

  TString NormFlag = NoNorm ? "_NoNorm" : "";
  canvas->Print("../Plots/MET/" + Selection + "_" + Region + NormFlag + "_MT2Comp.png");

}
