bool verbose = true;

void CheckShapes(TString ProcessName, TString UncertaintyName, TString ShapeName, TString HistoDirectory = "../minitrees/rootfiles2R/") {
  
  //TString HistoDirectory = "/eos/cms/store/user/scodella/Stop/MiniTrees/minitrees_36fb/rootfiles/";

  TFile *_fileu = TFile::Open(HistoDirectory + UncertaintyName + "up/Stop/" + ProcessName + ".root");
  TFile *_filec = TFile::Open(HistoDirectory + "nominal/Stop/"              + ProcessName + ".root");
  TFile *_filed = TFile::Open(HistoDirectory + UncertaintyName + "do/Stop/" + ProcessName + ".root");

  TH1F *HUP = (TH1F*) _fileu->Get(ShapeName);
  TH1F *HCN = (TH1F*) _filec->Get(ShapeName);
  TH1F *HDO = (TH1F*) _filed->Get(ShapeName);

  HUP->SetLineColor(2);
  HDO->SetLineColor(4);

  float min = HDO->GetMinimum();
  HCN->SetMinimum(0.9*min);
  float max = HCN->GetMaximum();
  HCN->SetMaximum(1.2*max);

  HCN->Draw();
  HUP->Draw("same");
  HDO->Draw("same");

  if (verbose) {
    for (int ibin = 1; ibin<=HCN->GetNbinsX(); ibin++) 
      cout << ibin << " " << HUP->GetBinContent(ibin)/HCN->GetBinContent(ibin) << " " << HDO->GetBinContent(ibin)/HCN->GetBinContent(ibin) << endl;
    cout << "Total " << HUP->Integral()/HCN->Integral() << " " << HDO->Integral()/HCN->Integral() << endl;
  }
}
