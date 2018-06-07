bool verbose = true;

void CheckShapes(TString ProcessName, TString UncertaintyName, TString ShapeName, TString HistoDirectory = "../minitrees/rootfiles/") {
  
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
  HCN->Draw("same");

  if (verbose) {
    for (int ibin = 1; ibin<=HCN->GetNbinsX(); ibin++) 
      cout << ibin << " " << HUP->GetBinContent(ibin)/HCN->GetBinContent(ibin) << " " << HDO->GetBinContent(ibin)/HCN->GetBinContent(ibin) << endl;
    cout << "Total " << HUP->Integral()/HCN->Integral() << " " << HDO->Integral()/HCN->Integral() << endl;
  }

}

void CheckShapes(TString ProcessName, TString UncertaintyName, TString ShapeName, TString HistoDir1, TString HistoDir2) {
  
  TFile *_file1 = TFile::Open(HistoDir1 + UncertaintyName + "/Stop/" + ProcessName + ".root");
  TFile *_file2 = TFile::Open(HistoDir2 + UncertaintyName + "/Stop/" + ProcessName + ".root");
  
  TH1F *H1 = (TH1F*) _file1->Get(ShapeName);
  TH1F *H2 = (TH1F*) _file2->Get(ShapeName);
  
  H1->SetLineColor(2);
  H2->SetLineColor(4);

  float min = H1->GetMinimum();
  H1->SetMinimum(0.9*min);
  float max = H1->GetMaximum();
  H1->SetMaximum(1.2*max);

  H1->Draw();
  H2->Draw("same");

  if (verbose) {
    cout << " " << H1->GetEntries() << " " << H2->GetEntries() << endl;
    cout << " " << H1->Integral() << " " << H2->Integral() << endl;
    for (int ibin = 1; ibin<=H1->GetNbinsX(); ibin++) 
      cout << ibin << " " << H2->GetBinContent(ibin)/H1->GetBinContent(ibin) << endl;
    cout << "Total " << H2->Integral()/H1->Integral() << endl;
  }

}
