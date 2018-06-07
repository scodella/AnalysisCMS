void MT2Gen(TString FileName, TString channel = "_ll") {

  TFile *UU = new TFile("../minitrees/rootfilesOct17/nominal/Stop/" + FileName + ".root", "read");

  //TH1F *C1V = (TH1F*) UU->Get("Stop/01_NoTag/h_MT2llgen" + channel);
  TH1F *D1T = (TH1F*) UU->Get("Stop/02_SR1gen_Tag/h_MT2llgen" + channel);
  TH1F *D1J = (TH1F*) UU->Get("Stop/02_SR1gen_NoTag/h_MT2llgen" + channel);
  TH1F *D10 = (TH1F*) UU->Get("Stop/02_SR1gen_NoJet/h_MT2llgen" + channel);
  TH1F *D2T = (TH1F*) UU->Get("Stop/02_SR2gen_Tag/h_MT2llgen" + channel);
  TH1F *D2J = (TH1F*) UU->Get("Stop/02_SR2gen_NoTag/h_MT2llgen" + channel);
  TH1F *D20 = (TH1F*) UU->Get("Stop/02_SR2gen_NoJet/h_MT2llgen" + channel);
  TH1F *D3T = (TH1F*) UU->Get("Stop/02_SR3gen_Tag/h_MT2llgen" + channel);
  TH1F *D3J = (TH1F*) UU->Get("Stop/02_SR3gen_NoTag/h_MT2llgen" + channel);
  TH1F *D30 = (TH1F*) UU->Get("Stop/02_SR3gen_NoJet/h_MT2llgen" + channel);

  D1T->Add(D1J); D1T->Add(D10); D1T->Add(D2T); D1T->Add(D2J); D1T->Add(D20); D1T->Add(D3T); D1T->Add(D3J); D1T->Add(D30); 

  TH1F *R1T = (TH1F*) UU->Get("Stop/02_SR1_Tag/h_MT2ll" + channel);
  TH1F *R1J = (TH1F*) UU->Get("Stop/02_SR1_NoTag/h_MT2ll" + channel);
  TH1F *R10 = (TH1F*) UU->Get("Stop/02_SR1_NoJet/h_MT2ll" + channel);
  TH1F *R2T = (TH1F*) UU->Get("Stop/02_SR2_Tag/h_MT2ll" + channel);
  TH1F *R2J = (TH1F*) UU->Get("Stop/02_SR2_NoTag/h_MT2ll" + channel);
  TH1F *R20 = (TH1F*) UU->Get("Stop/02_SR2_NoJet/h_MT2ll" + channel);
  TH1F *R3T = (TH1F*) UU->Get("Stop/02_SR3_Tag/h_MT2ll" + channel);
  TH1F *R3J = (TH1F*) UU->Get("Stop/02_SR3_NoTag/h_MT2ll" + channel);
  TH1F *R30 = (TH1F*) UU->Get("Stop/02_SR3_NoJet/h_MT2ll" + channel);

  R1T->Add(R1J); R1T->Add(R10); R1T->Add(R2T); R1T->Add(R2J); R1T->Add(R20); R1T->Add(R3T); R1T->Add(R3J); R1T->Add(R30); 

  D1T->SetLineColor(2);
  D1T->SetMarkerStyle(20);
  D1T->SetMarkerColor(2);
  D1T->Draw();
  R1T->SetLineColor(1);
  R1T->Draw("same");

  cout << D1T->GetEntries() << " " << D1T->Integral() << " " << R1T->GetEntries() << " " << R1T->Integral() << " " << D1T->Integral()/R1T->Integral() << endl;
  cout << D1T->Integral()/D1T->GetEntries() << "        " << R1T->Integral()/R1T->GetEntries() << endl;

}
