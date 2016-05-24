#define AnalysisControl_cxx
#include "../include/AnalysisControl.h"


//------------------------------------------------------------------------------
// AnalysisControl
//------------------------------------------------------------------------------
AnalysisControl::AnalysisControl(TTree* tree, TString systematic) : AnalysisCMS(tree, systematic)
{
  SetSaveMinitree(false);
}


//------------------------------------------------------------------------------
// Loop
//------------------------------------------------------------------------------
void AnalysisControl::Loop(TString analysis, TString filename, float luminosity)
{
  if (fChain == 0) return;

  Setup(analysis, filename, luminosity);


  // Define histograms
  //----------------------------------------------------------------------------
  root_output->cd();

  for (int j=0; j<ncut; j++) {

    for (int k=0; k<=njetbin; k++) {

      TString sbin = (k < njetbin) ? Form("/%djet", k) : "";

      TString directory = scut[j] + sbin;

      root_output->cd();

      if (k < njetbin) gDirectory->mkdir(directory);

      root_output->cd(directory);

      for (int i=ee; i<=ll; i++) {

	TString suffix = "_" + schannel[i];

	DefineHistograms(i, j, k, suffix);
      }
    }
  }

  root_output->cd();


  // Loop over events
  //----------------------------------------------------------------------------
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    PrintProgress(jentry, _nentries);

    EventSetup();


    // Analysis
    //--------------------------------------------------------------------------
    if (Lepton1.v.Pt() < 10.) continue;
    if (Lepton2.v.Pt() < 10.) continue;

    _nelectron = 0;

    if (abs(Lepton1.flavour) == ELECTRON_FLAVOUR) _nelectron++;
    if (abs(Lepton2.flavour) == ELECTRON_FLAVOUR) _nelectron++;

    if      (_nelectron == 2) _channel = ee;
    else if (_nelectron == 1) _channel = em;
    else if (_nelectron == 0) _channel = mm;
    
    _m2l  = mll;   // Needs l2Sel
    _pt2l = ptll;  // Needs l2Sel

    bool opposite_sign = (Lepton1.flavour * Lepton2.flavour < 0);

    bool pass;


    // Z+jets
    //--------------------------------------------------------------------------
    pass = true;

    pass &= opposite_sign;
    pass &= (Lepton1.v.Pt() > 20.);
    pass &= (Lepton2.v.Pt() > 20.);
    pass &= (std_vector_lepton_pt->at(2) < 10.);
    pass &= (_m2l > 12.);
    pass &= (_nbjet20cmvav2l == 0);

    FillLevelHistograms(Control_00_ZJets, pass);


    // Top
    //--------------------------------------------------------------------------
    pass = true;

    pass &= opposite_sign;
    pass &= (Lepton1.v.Pt() > 20.);
    pass &= (Lepton2.v.Pt() > 20.);
    pass &= (std_vector_lepton_pt->at(2) < 10.);
    pass &= (_m2l > 12.);
    pass &= (_nbjet20cmvav2l > 0);
    pass &= (_njet > 1);
    pass &= (_channel == em || fabs(_m2l - Z_MASS) > 15.);
    pass &= (_channel == em || MET.Et() > 40.);

    FillLevelHistograms(Control_01_Top, pass);


    // AN-15-325, latinos
    // WW cross section measurement at sqrt(s) = 13 TeV
    //--------------------------------------------------------------------------
    pass = true;

    pass &= opposite_sign;
    pass &= (Lepton1.v.Pt() > 20.);
    pass &= (Lepton2.v.Pt() > 20.);
    pass &= (std_vector_lepton_pt->at(2) < 10.);
    pass &= (_m2l > 12.);
    pass &= (_nbjet20cmvav2l == 0);
    pass &= (MET.Et() > 20.);
    pass &= (mpmet > 20.);
    pass &= (_pt2l > 30.);
    pass &= (_channel == em || fabs(_m2l - Z_MASS) > 15.);
    pass &= (_channel == em || MET.Et() > 40.);
    pass &= (_channel == em || mpmet > 40.);
    pass &= (_channel == em || _pt2l > 45.);

    FillLevelHistograms(Control_02_WW0j, pass && _njet == 0);
    FillLevelHistograms(Control_03_WW1j, pass && _njet == 1);

    if (pass && _njet == 0 && _channel == em) GetRecoWeightsLHE(list_vectors_weights_0jet);
    if (pass && _njet == 1 && _channel == em) GetRecoWeightsLHE(list_vectors_weights_1jet);
  }


  EndJob();
}


//------------------------------------------------------------------------------
// FillLevelHistograms
//------------------------------------------------------------------------------
void AnalysisControl::FillLevelHistograms(int  icut,
					  bool pass)
{
  if (!pass) return;

  FillHistograms(_channel, icut, _jetbin);
  FillHistograms(_channel, icut, njetbin);
}