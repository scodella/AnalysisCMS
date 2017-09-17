#ifndef AnalysisStop_h
#define AnalysisStop_h
 
#include "AnalysisCMS.h"
#include "../../BTagSFUtil/BTagSFUtil.C"
#include <map>


class AnalysisStop : public AnalysisCMS
{
 public :

  AnalysisStop(TTree* tree, TString systematic);
  AnalysisStop(TFile* file, TString systematic,	int FillAllHistograms = 1);

  void BookAnalysisHistograms(); 

  void BookSystematicHistograms();

  void BookTheoreticalVariationsHistograms(); 

  void GetAnalysisVariables(); 

  void FillAnalysisHistograms(int     ichannel,
			      int     icut,
			      int     ijet); 

  void FillSystematicHistograms(int     ichannel,
				int     icut); 

  void FillTheoreticalVariationsHistograms(int     ichannel,
					   int     icut);

  void FillLevelHistograms   (int     icut,
			      bool    pass);

  void SaveSystematicHistograms();

  void SaveTheoreticalVariationsHistograms();

  void Loop                  (TString analysis,
			      TString sample,
			      float   luminosity,
			      float   StopRefMass = -1.,
			      float   NeutralinoRefMass = -1.);

  void CorrectEventWeight    ();

  void SetSUSYProductionMap  ();

  void GetMiniTree           (TFile *MiniTreeFile, TString systematic);

  bool PassFastsimJetsCleanup();

  TString FastSimDataset;
  BTagSFUtil *BTagSF, *BTagSF_Upb, *BTagSF_Dob, *BTagSF_UpFSb, *BTagSF_DoFSb;

  TString SUSYProductionProcess;

  typedef pair<int, int> MassPoint;
  typedef pair<float, float> StopCrossSection;
  typedef pair<StopCrossSection, int> MassPointParameters;

  typedef map<MassPoint, MassPointParameters> MassPointMap;
  MassPointMap StopNeutralinoMap;
    
  typedef map<MassPoint, float> MassPointISRMap;
  MassPointISRMap StopNeutralinoISRMap;

  // Analysis histograms
  //----------------------------------------------------------------------------
  TH1D*                  h_mt2lblbcomb      [nchannel][ncut][njetbin+1];
  TH1D*                  h_mt2bbtrue        [nchannel][ncut][njetbin+1];
  TH1D*                  h_mt2lblbtrue      [nchannel][ncut][njetbin+1];
  TH1D*                  h_mt2lblbmatch     [nchannel][ncut][njetbin+1];
  TH1D*                  h_mlb1comb         [nchannel][ncut][njetbin+1];
  TH1D*                  h_mlb2comb         [nchannel][ncut][njetbin+1];
  TH1D*                  h_mlb1true         [nchannel][ncut][njetbin+1];
  TH1D*                  h_mlb2true         [nchannel][ncut][njetbin+1];
  TH2D*                  h_mt2lblbvsmlbtrue [nchannel][ncut][njetbin+1];
  TH1D*                  h_nisrjet          [nchannel][ncut][njetbin+1];
  TH1D*                  h_maxjetpt         [nchannel][ncut][njetbin+1];
  TH1D*                  h_njet20           [nchannel][ncut][njetbin+1];
  TH1D*                  h_njet20dphilmet   [nchannel][ncut][njetbin+1];
  TH1D*                  h_njet30           [nchannel][ncut][njetbin+1];
  TH1D*                  h_MET              [nchannel][ncut][njetbin+1];
  TH1D*                  h_Counter          [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiminlepmet    [nchannel][ncut][njetbin+1];
  TH1D*                  h_HT               [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTm150     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTp150     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTm200     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTp200     [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphillMET        [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLL           [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin1       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin2       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin3       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin4       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin5       [nchannel][ncut][njetbin+1];
  TH1D*                  h_M1ll             [nchannel][ncut][njetbin+1];
  TH1D*                  h_M2ll             [nchannel][ncut][njetbin+1];

  int _SaveHistograms, _DoTheoreticalVariations;

  float _metmeff, _MT2ll, _MT2llgen, _MT2llfake;
  TH1D*                  h_metmeff            [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2ll              [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llgen           [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisr           [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisrgen        [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2ll_nvtxup       [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2ll_nvtxdo       [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llgen_nvtxup    [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llgen_nvtxdo    [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisr_nvtxup    [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisr_nvtxdo    [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisrgen_nvtxup [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2llisrgen_nvtxdo [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2ll_fake         [nchannel][ncut][njetbin+1];
  TH1F*                  h_MT2ll_truth        [nchannel][ncut][njetbin+1];
  TH1F*                  h_MET_fake           [nchannel][ncut][njetbin+1];
  TH1F*                  h_MET_truth          [nchannel][ncut][njetbin+1];

  bool  _hasisrjet;

  float _dphiminlmet, _m2lZ2;

  int   _nLeptonsMatched, _njet30;

  float _MT2_Met; int NbinsMT2 = 7; int NbinsMet = 5;
  float vMinMT2 = 0., vMinMet = 0., vMaxMT2 = 140., vMaxMet = 500.;
  TH1D*                  h_MT2_Met          [nchannel][ncut][njetbin+1];

  float _HTvisible_Met; int NbinsHTvisible = 5;
  float vMinHTvisible = 0., vMaxHTvisible = 400.;
  TH1D*                  h_HTvisible_Met    [nchannel][ncut][njetbin+1];

  float _MetMeff_Met; int NbinsMetMeff = 8;
  float vMinMetMeff = 0., vMaxMetMeff = 0.8;
  TH1D*                  h_MetMeff_Met      [nchannel][ncut][njetbin+1];

  float _R2_Met; int NbinsR2 = 7;
  float vMinR2 = 0., vMaxR2 = 1.4;
  TH1D*                  h_R2_Met           [nchannel][ncut][njetbin+1];

  // Systematic output
  TFile*                 root_output_systematic[nsystematic];
  TH1F*                  h_MT2ll_systematic       [nchannel][ncut][nsystematic];
  TH1F*                  h_MT2llgen_systematic    [nchannel][ncut][nsystematic];
  TH1F*                  h_MT2llisr_systematic    [nchannel][ncut][nsystematic];
  TH1F*                  h_MT2llisrgen_systematic [nchannel][ncut][nsystematic];

  // Theoretical variations output
  const int nTheoreticalVariations = 111;
  float _TheoreticalVariationRenormalization[111];
  TFile*                 root_output_theoreticalvariations;
  TH1F*                  h_MT2ll_theoreticalvariation        [nchannel][ncut][111];
  TH1F*                  h_MT2llgen_theoreticalvariation     [nchannel][ncut][111];
  TH1F*                  h_MT2llisr_theoreticalvariation     [nchannel][ncut][111];
  TH1F*                  h_MT2llisrgen_theoreticalvariation  [nchannel][ncut][111];

  //
  bool _applyDYcorrections;

  // Tools for ISR jet reweighting 
  // https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf
  bool _applyisrreweighting = true;
  float _event_weight_Isrnjetup, _event_weight_Isrnjetdo;
  float _event_weight_Fakeup,    _event_weight_Fakedo;

  void ApplyISRReweighting();

  const int nISRMultiplicityBins = 7;
  float ISRBinWeight[7] = {1., 0.920, 0.821, 0.715, 0.662, 0.561, 0.511};

  const int nISRPtBins = 8;
  float ISRPtBinWeight[8] = {1., 1.052, 1.179, 1.150, 1.057, 1.000, 0.912, 0.783};
  float ISRPtBinEdge[8] = {0., 50., 100., 150., 200., 300., 400., 600.};

  bool        ShapeWZtoWW         ();

  bool        ShapeZZ             ();

  bool        ShapeZWtoZ          ();

};

#endif
