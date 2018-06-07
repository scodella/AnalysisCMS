#ifndef AnalysisStop_h
#define AnalysisStop_h
 
#include "AnalysisCMS.h"
#include "../../BTagSFUtil/BTagSFUtil.C"
#include <map>
#include "TEfficiency.h"

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

  float PileupSyst(float runperiod, float putrue, int variation);

  TString FastSimDataset;
  BTagSFUtil *BTagSF, *BTagSF_Upb, *BTagSF_Dob, *BTagSF_Upl, *BTagSF_Dol, *BTagSF_UpFSb, *BTagSF_DoFSb;

  TEfficiency *TrgEff_ee, *TrgEff_mm, *TrgEff_em;
  void ApplyStopTriggerEfficiency();

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
  TH1D*                  h_njet20dphilmet   [nchannel][ncut][njetbin+1];
  TH1D*                  h_njet30           [nchannel][ncut][njetbin+1];
  TH1D*                  h_Counter          [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiminlepmet    [nchannel][ncut][njetbin+1];
  TH1D*                  h_HT               [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTm150     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTp150     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTm200     [nchannel][ncut][njetbin+1];
  TH1D*                  h_MT2ll_HTp200     [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphillMET        [nchannel][ncut][njetbin+1];
  TH2D*                  h_ptdphiLL         [nchannel][ncut][njetbin+1];
  TH1D*                  h_ptLL             [nchannel][ncut][njetbin+1];
  TH1D*                  h_m2L              [nchannel][ncut][njetbin+1];
  TH1D*                  h_ptLLbins         [nchannel][ncut][njetbin+1];
  TH1D*                  h_genptLL          [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLL           [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin1       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin2       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin3       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin4       [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphiLLbin5       [nchannel][ncut][njetbin+1];
  TH1D*                  h_M1ll             [nchannel][ncut][njetbin+1];
  TH1D*                  h_M2ll             [nchannel][ncut][njetbin+1];
  TH1D*                  h_njetISR          [nchannel][ncut][njetbin+1];
  TH1D*                  h_ptjetISR         [nchannel][ncut][njetbin+1];
  TH1F*                  h_mt2LL            [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphijetISR       [nchannel][ncut][njetbin+1];
  TH1D*                  h_MET              [nchannel][ncut][njetbin+1];
  TH1D*                  h_njet20           [nchannel][ncut][njetbin+1];
  TH1D*                  h_nbjet            [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep1Pt           [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep2Pt           [nchannel][ncut][njetbin+1];
  TH1D*                  h_JetPt            [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep1Eta          [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep2Eta          [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep1Phi          [nchannel][ncut][njetbin+1];
  TH1D*                  h_Lep2Phi          [nchannel][ncut][njetbin+1];
  TH1D*                  h_Jet1Pt           [nchannel][ncut][njetbin+1];
  TH1D*                  h_Jet2Pt           [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphil1MET        [nchannel][ncut][njetbin+1];
  TH1D*                  h_dphil2MET        [nchannel][ncut][njetbin+1];
  TH1D*                  h_METphi           [nchannel][ncut][njetbin+1];

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
  TH2F*                  h_MT2ll_MET          [nchannel][ncut][njetbin+1];
  TH1F*                  h_dphiisrmet         [nchannel][ncut][njetbin+1];

  bool  _hasisrjet; int _njetISR; float _ptjetISR, _dphijetISR;

  float _dphiminlmet, _m2lZ2, _dphiisrmet;

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

  TH1D*                  h_njetISR_systematic          [nchannel][ncut][nsystematic];
  TH1D*                  h_ptjetISR_systematic         [nchannel][ncut][nsystematic];
  TH1F*                  h_mt2LL_systematic            [nchannel][ncut][nsystematic];
  TH1D*                  h_dphijetISR_systematic       [nchannel][ncut][nsystematic];
  TH1D*                  h_MET_systematic              [nchannel][ncut][nsystematic];
  TH1D*                  h_njet20_systematic           [nchannel][ncut][nsystematic];
  TH1D*                  h_nbjet_systematic            [nchannel][ncut][nsystematic];
  TH1D*                  h_Lep1Pt_systematic           [nchannel][ncut][nsystematic];
  TH1D*                  h_Lep2Pt_systematic           [nchannel][ncut][nsystematic];
  TH1D*                  h_JetPt_systematic            [nchannel][ncut][nsystematic];

  // Theoretical variations output
  const int nTheoreticalVariations = 111;
  float _TheoreticalVariationRenormalization[111];
  TFile*                 root_output_theoreticalvariations;
  TH1F*                  h_MT2ll_theoreticalvariation        [nchannel][ncut][111];
  TH1F*                  h_MT2llgen_theoreticalvariation     [nchannel][ncut][111];
  TH1F*                  h_MT2llisr_theoreticalvariation     [nchannel][ncut][111];
  TH1F*                  h_MT2llisrgen_theoreticalvariation  [nchannel][ncut][111];

  TH1D*                  h_njetISR_theoreticalvariation          [nchannel][ncut][111];
  TH1D*                  h_ptjetISR_theoreticalvariation         [nchannel][ncut][111];
  TH1F*                  h_mt2LL_theoreticalvariation            [nchannel][ncut][111];
  TH1D*                  h_dphijetISR_theoreticalvariation       [nchannel][ncut][111];
  TH1D*                  h_MET_theoreticalvariation              [nchannel][ncut][111];
  TH1D*                  h_njet20_theoreticalvariation           [nchannel][ncut][111];
  TH1D*                  h_nbjet_theoreticalvariation            [nchannel][ncut][111];
  TH1D*                  h_Lep1Pt_theoreticalvariation           [nchannel][ncut][111];
  TH1D*                  h_Lep2Pt_theoreticalvariation           [nchannel][ncut][111];
  TH1D*                  h_JetPt_theoreticalvariation            [nchannel][ncut][111];

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

  float _massZcandidate;

  bool        ShapeZZ             ();

  bool        ShapeZWtoZ          ();

  vector<int> _lheEvent;
  TTree *lheTree; 
  bool hasLHE;

  int const nEtaEdges = 6; int const nRuns = 152;
  float EtaEdge[6] = {0., 1., 1.4442, 1.566, 2., 2.5};
  float _elescale[5][2][152];
  int Run1[152] = {1, 3, 273158, 273446, 273494, 273725, 273728, 274094, 274172, 274241, 274244, 274251, 274316, 274317, 274336, 274344, 274388, 274422, 274440, 274442, 274968, 274969, 274970, 274998, 274999, 275000, 275068, 275073, 275124, 275282, 275310, 275311, 275337, 275344, 275375, 275376, 275657, 275776, 275782, 275832, 275835, 275847, 275886, 275911, 275913, 275920, 276242, 276244, 276282, 276283, 276315, 276361, 276363, 276384, 276437, 276454, 276501, 276525, 276527, 276543, 276582, 276585, 276587, 276655, 276659, 276776, 276794, 276808, 276811, 276831, 276870, 276935, 276950, 277069, 277076, 277087, 277096, 277112, 277148, 277168, 277194, 277305, 277981, 278018, 278167, 278175, 278240, 278308, 278310, 278345, 278349, 278406, 278509, 278769, 278801, 278808, 278820, 278822, 278873, 278923, 278962, 278969, 278975, 279024, 279115, 279479, 279654, 279667, 279694, 279715, 279716, 279760, 279766, 279794, 279841, 279844, 279931, 281613, 281693, 281707, 281726, 281797, 281976, 282037, 282092, 282708, 282735, 282800, 282814, 282842, 283042, 283052, 283270, 283283, 283305, 283308, 283353, 283358, 283408, 283416, 283478, 283548, 283820, 283830, 283865, 283876, 283877, 283885, 283934, 283946, 283964, 284025};
  int Run2[152] = {2, 273157, 273445, 273493, 273724, 273727, 274093, 274171, 274240, 274243, 274250, 274315, 274316, 274335, 274343, 274387, 274421, 274439, 274441, 274967, 274968, 274969, 274997, 274998, 274999, 275067, 275072, 275123, 275281, 275309, 275310, 275336, 275343, 275374, 275375, 275656, 275775, 275781, 275831, 275834, 275846, 275885, 275910, 275912, 275919, 276241, 276243, 276281, 276282, 276314, 276360, 276362, 276383, 276436, 276453, 276500, 276524, 276526, 276542, 276581, 276584, 276586, 276654, 276658, 276775, 276793, 276807, 276810, 276830, 276869, 276934, 276949, 277068, 277075, 277086, 277095, 277111, 277147, 277167, 277193, 277304, 277980, 278017, 278166, 278174, 278239, 278307, 278309, 278344, 278348, 278405, 278508, 278768, 278800, 278807, 278819, 278821, 278872, 278922, 278961, 278968, 278974, 279023, 279114, 279478, 279653, 279666, 279693, 279714, 279715, 279759, 279765, 279793, 279840, 279843, 279930, 281612, 281692, 281706, 281725, 281796, 281975, 282036, 282091, 282707, 282734, 282799, 282813, 282841, 283041, 283051, 283269, 283282, 283304, 283307, 283352, 283357, 283407, 283415, 283477, 283547, 283819, 283829, 283864, 283875, 283876, 283884, 283933, 283945, 283963, 284024, 284044};
  void  FillElescale();
  float GetElescale(float eleeta, float eler9, int elerun);

};

#endif
