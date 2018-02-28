#define AnalysisStop_cxx
#include "../include/AnalysisStop.h"
#include <fstream>
#include <iostream>

//------------------------------------------------------------------------------
// AnalysisStop
//------------------------------------------------------------------------------
AnalysisStop::AnalysisStop(TTree* tree, TString systematic) : AnalysisCMS(tree, systematic)
{
  _applyDYcorrections = systematic.Contains("DYcorr") ? true : false;
  systematic.ReplaceAll("DYcorr", "");
  _applytopptreweighting = systematic.Contains("NoTopPt") ? false : true;
  systematic.ReplaceAll("NoTopPt", "");
  if (systematic.Contains("nominal")) {
    SetSaveMinitree(true);
    _SaveHistograms = 0;
  } else if (systematic.Contains("multilepton") || _systematic.Contains("MetTrigger")) {
    SetSaveMinitree(true);
    _SaveHistograms = -1;
  } else {
    SetSaveMinitree(false);
    if (systematic.Contains("Minitrees")) SetSaveMinitree(true);
    _SaveHistograms = 1;
    if (systematic.Contains("JES") || systematic.Contains("MET") || systematic.Contains("DYshape") || systematic.Contains("Pileup") || systematic.Contains("SS")) _SaveHistograms = 2;
  }
  _DoTheoreticalVariations = 0;
  if (systematic=="theory") { _DoTheoreticalVariations = 1; _SaveHistograms = -1; }
}

AnalysisStop::AnalysisStop(TFile* MiniTreeFile, TString systematic, int SaveHistograms) 
{
  _applyDYcorrections = systematic.Contains("DYcorr") ? true : false;
  _applytopptreweighting = true;
  systematic.ReplaceAll("DYcorr", "");
  SetSaveMinitree(false);
  GetMiniTree(MiniTreeFile, systematic);
  _systematic = systematic;
  _SaveHistograms = SaveHistograms;
  if (systematic!="nominal") _SaveHistograms = 1;
  if (systematic.Contains("JES") || systematic.Contains("MET") || systematic.Contains("DYshape") || systematic.Contains("Pileup")) _SaveHistograms = 2;
  _DoTheoreticalVariations = (_SaveHistograms>3) ? 1 : 0;
  if (systematic=="theory") { _DoTheoreticalVariations = 1; _SaveHistograms = -1; }
}

//------------------------------------------------------------------------------
// Loop
//------------------------------------------------------------------------------
void AnalysisStop::Loop(TString analysis, TString filename, float luminosity, float StopRefMass, float NeutralinoRefMass)
{
  if (fChain == 0) return;

  TString MassPointFlag = "";

  SUSYProductionProcess = "";
  if (filename.Contains("T2tt")) SUSYProductionProcess = "T2tt";
  if (filename.Contains("T2bW")) SUSYProductionProcess = "T2bW";
  if (filename.Contains("T2bt")) SUSYProductionProcess = "T2bt";
  if (filename.Contains("_dM10to80")) SUSYProductionProcess += "_dM10to80";
  if (filename.Contains("TChiWW")) SUSYProductionProcess = "TChiWW";
  if (filename.Contains("TChiSlep")) SUSYProductionProcess = "TChiSlep";
  if (filename.Contains("TChiStau")) SUSYProductionProcess = "TChiStau";

  if (SUSYProductionProcess!="") {

    SetSUSYProductionMap();

    if (StopRefMass==-1.) {
  
      if (filename.Contains("150to250"))       { StopRefMass = 150.; NeutralinoRefMass =  25.; }
      else if (filename.Contains("250to350"))  { StopRefMass = 275.; NeutralinoRefMass = 150.; }
      else if (filename.Contains("350to400"))  { StopRefMass = 350.; NeutralinoRefMass = 225.; }
      else if (filename.Contains("400to1200")) { StopRefMass = 450.; NeutralinoRefMass = 325.; }
      else                                     { StopRefMass = 350.; NeutralinoRefMass = 225.; }

    }

    int iStopRefMass = StopRefMass, iNeutralinoRefMass = NeutralinoRefMass;
    MassPointFlag = "_Sm"; MassPointFlag += iStopRefMass; MassPointFlag += "_Xm"; MassPointFlag += iNeutralinoRefMass;
    if (SUSYProductionProcess.Contains("TChi")) MassPointFlag.ReplaceAll("_Sm", "_Xm");

  }

  //_applytopptreweighting = false;//true;
  if (!_applytopptreweighting && filename.Contains("TTTo2L2Nu")) {
    MassPointFlag = "_NoTopPt";
  }

  int applyZZkfactors = -1;
  if (filename.Contains("ZZTo") && !filename.Contains("qqHToZZ") && !filename.Contains("ZZTo2L2Q")) {
    if (filename.Contains("gg")) applyZZkfactors = 0;
    else if (filename.Contains("ZZTo4")) applyZZkfactors = 1;
    else applyZZkfactors = 2;
    ReadGluGluHZZkfactors();
  }

  if (_applyDYcorrections) {
    if (!filename.Contains("DY") && !filename.Contains("ZJets")) _applyDYcorrections = false;
    else MassPointFlag = "_DYcorr";
  }
  if (_systematic.Contains("DYshape") && (filename.Contains("DY") || filename.Contains("ZJets"))) _applyDYcorrections = true;

  Setup(analysis, filename, luminosity, MassPointFlag);
 
  // Define histograms
  //----------------------------------------------------------------------------
  if (_SaveHistograms>=0) BookAnalysisHistograms();
  if (_SaveHistograms>=3) BookSystematicHistograms();
  if (_DoTheoreticalVariations==1) BookTheoreticalVariationsHistograms();
  
  root_output->cd();

  FastSimDataset = _isfastsim ? "_T2" : "";
  BTagSF       = new BTagSFUtil("mujets", "CSVv2", "Medium",   0, FastSimDataset);
  BTagSF_Upb   = new BTagSFUtil("mujets", "CSVv2", "Medium",  +1, FastSimDataset);
  BTagSF_Dob   = new BTagSFUtil("mujets", "CSVv2", "Medium",  -1, FastSimDataset);
  BTagSF_Upl   = new BTagSFUtil("mujets", "CSVv2", "Medium",  +3, FastSimDataset);
  BTagSF_Dol   = new BTagSFUtil("mujets", "CSVv2", "Medium",  -3, FastSimDataset);
  BTagSF_UpFSb = new BTagSFUtil("mujets", "CSVv2", "Medium", +11, FastSimDataset);
  BTagSF_DoFSb = new BTagSFUtil("mujets", "CSVv2", "Medium", -11, FastSimDataset);
  /*BTagSF       = new BTagSFUtil("mujets", "DeepCSV", "Medium",   0, FastSimDataset);
  BTagSF_Upb   = new BTagSFUtil("mujets", "DeepCSV", "Medium",  +1, FastSimDataset);
  BTagSF_Dob   = new BTagSFUtil("mujets", "DeepCSV", "Medium",  -1, FastSimDataset);
  BTagSF_UpFSb = new BTagSFUtil("mujets", "DeepCSV", "Medium", +11, FastSimDataset);
  BTagSF_DoFSb = new BTagSFUtil("mujets", "DeepCSV", "Medium", -11, FastSimDataset);*/

  // Stop trigger efficiencies
  TFile TriggerEfficiencyFile("/afs/cern.ch/work/s/scodella/Stop/CodeDevelopment/CMSSW_8_0_26_patch1/src/AnalysisCMS/stop/trigger/CopyOfFrameWork/output_stop_pt.root");
  TrgEff_ee = (TEfficiency *) TriggerEfficiencyFile.Get("DATAdenominatoree2D_clone");
  TrgEff_mm = (TEfficiency *) TriggerEfficiencyFile.Get("DATAdenominatormm2D_clone");
  TrgEff_em = (TEfficiency *) TriggerEfficiencyFile.Get("DATAdenominatorem2D_clone");

  // Loop over events
  //----------------------------------------------------------------------------
  int nnn = 0, lll = 0;
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {
    
    if (!_isminitree) {

      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;

    }
    
    fChain->GetEntry(jentry);
    
    PrintProgress(jentry, _nentries);

    bool pass_masspoint = true;
    if (SUSYProductionProcess!="") {
      float susyMx = (SUSYProductionProcess.Contains("T2")) ? susyMstop : susyMChargino;
      pass_masspoint = (susyMx==StopRefMass && susyMLSP==NeutralinoRefMass) ? true : false;
    }
    if (!pass_masspoint && !_saveminitree) continue;

    if (!_isminitree) {

      // I found an event with metPfType1 = nan, this cause MT2 calculation to stuck
      if (metPfType1!=metPfType1) continue;
      EventSetup(2.4, 20.);

      if (_systematic.Contains("MetTrigger")) {
	int METtriggers[] = {63, 64, 65, 66, 68, 79, 80, 81, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92};
	for(int a = 0; a < 18; a++) {
	  if(std_vector_trigger->at(METtriggers[a]) != 0) {
	    minitree->Fill();  
	    break;
	  }
	}
	continue;
      }

      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17#Cleaning_up_of_fastsim_jets_from
      if (!PassFastsimJetsCleanup()) continue;
      
      if (_event_weight==0.) continue;
 
    }

    // Get analysis variables
    //--------------------------------------------------------------------------
    
    GetAnalysisVariables();

    // Analysis
    //--------------------------------------------------------------------------
    
    if (!_isminitree) {

      // To merge DY_LO and DY_HTbins( always ht >70 GeV) samples it is needed require ht > 70 in DY_LO sample; 
      // The analysis default sample is DY_NNLO  
      if (filename.Contains("DY") && filename.Contains("LO") && !filename.Contains("HT")) 
	if (_htgen>70.) continue;

    }

    // 3L sidebands
    if (_systematic.Contains("WZ3L") || _systematic.Contains("ttZ")) {
      
      if (_nlepton<3) continue;
      if (MET.Et()<140. && !_systematic.Contains("NoMetCut")) continue; 
      if (_ntightlepton<3) continue; 

      //if (fabs(_lep3mid)==24 || fabs(_lep3mid)==23 || fabs(_lep3mid)==15) _nLeptonsMatched++;   

      //if (_systematic.Contains("WZ3L") && AnalysisLeptons[0].v.Pt()<30.) continue;
      //if (_systematic.Contains("ttZ") && AnalysisLeptons[0].v.Pt()<25.) continue;
      if (AnalysisLeptons[0].v.Pt()<25.) continue;
      if (AnalysisLeptons[1].v.Pt()<20.) continue;
      if (AnalysisLeptons[2].v.Pt()<20.) continue;
      
      float DMZ = 999.; int Zlep1 = -1, Zlep2 = -1, Wlep = -1;
      for (int l1 = 0; l1<3; l1++) {
	for (int l2 = l1+1; l2<3; l2++) {
	  if (AnalysisLeptons[l1].flavour*AnalysisLeptons[l2].flavour<0. && 
	      fabs(AnalysisLeptons[l1].flavour)==fabs(AnalysisLeptons[l2].flavour)) {
	    
	    TLorentzVector Zcand = AnalysisLeptons[l1].v + AnalysisLeptons[l2].v;
	    float DMZcand = fabs( Zcand.M() - Z_MASS );
	    if (DMZcand<DMZ) {
	      
	      Zlep1 = l1;
	      Zlep2 = l2;
	      
	      if (l1==0 && l2==1) Wlep = 2;
	      if (l1==0 && l2==2) Wlep = 1;
	      if (l1==1 && l2==2) Wlep = 0;
	      
	      DMZ = DMZcand;
	      
	    }
	    
	  }
	}
      }

      if (DMZ>=10. || Zlep1<0 || Zlep2<0 || Wlep<0) continue; 
      else {

	if (abs(AnalysisLeptons[Wlep].flavour)==ELECTRON_FLAVOUR) _channel = ee;
	else _channel = mm;

	if (_systematic.Contains("WZ3L") && AnalysisLeptons[0].v.Pt()<30.) {
	  FillLevelHistograms(Stop_00_VZ, _leadingPtCSVv2M<20.);
	} else if (_systematic.Contains("ttZ")) {
	  FillLevelHistograms(Stop_00_ttZ_Tag,  _leadingPtCSVv2M>20.);// && _njet>=2); 
	  FillLevelHistograms(Stop_00_ttZ_2Tag, _leadingPtCSVv2M>30.);//_trailingPtCSVv2M>20.); 
	}
	
      }

      continue;
      
    } // end 3L sidebands
    
    if (!_isminitree) {

      if (!_systematic.Contains("fake") && !_systematic.Contains("invertveto") && 
	  !_systematic.Contains("WZtoWW") && !_systematic.Contains("ZZ") && !_systematic.Contains("MetTrigger") && 
	  !_systematic.Contains("ZWtoZ") && !_systematic.Contains("multilepton") && _nlepton != 2) continue; // 2 and only 2 leptons
      //if (_systematic.Contains("invertveto") && _nlepton < 3) continue;
      //if (_systematic.Contains("invertvetotight") && _ntightlepton<3) continue;
      if (_systematic.Contains("fake") && (_nlepton <= 2 || _ntightlepton>2)) continue;      

      //if (!_systematic.Contains("SS")  && !_systematic.Contains("multilepton") && Lepton1.flavour * Lepton2.flavour > 0) continue; 
      if (_systematic.Contains("SS") && Lepton1.flavour * Lepton2.flavour < 0) continue;
      if (_systematic.Contains("SSp") && (Lepton1.flavour <0 || Lepton2.flavour < 0)) continue;
      if (_systematic.Contains("SSm") && (Lepton1.flavour >0 || Lepton2.flavour > 0)) continue;

      if (_ismc) CorrectEventWeight(); 

      float LeadingLeptonPtCut = (_systematic.Contains("multilepton")) ? 20. : 25.;
      if (Lepton1.v.Pt() < LeadingLeptonPtCut) continue;
      if (Lepton2.v.Pt() < 20.) continue;
      
      _nelectron = 0;

      if (abs(Lepton1.flavour) == ELECTRON_FLAVOUR) _nelectron++;
      if (abs(Lepton2.flavour) == ELECTRON_FLAVOUR) _nelectron++;
    
      if      (_nelectron == 2) _channel = ee;
      else if (_nelectron == 1) _channel = em;
      else if (_nelectron == 0) _channel = mm; 
    
      _m2l  = mll;
      _pt2l = ptll;
      
      if (_systematic.Contains("multilepton")) {
	minitree->Fill();      
	continue;
      }
      
    }
      
    _njet30 = 0;
    for (int ij = 0; ij<_jet_pt.size(); ij++)
      if (_jet_pt.at(ij)>=30.) _njet30++;

    bool passShapeAnalysis = true;
    if (_systematic.Contains("WZtoWW")) passShapeAnalysis = ShapeWZtoWW();
    if (_systematic.Contains("ZWtoZ"))  passShapeAnalysis = ShapeZWtoZ();
    if (_systematic.Contains("ZZ"))     passShapeAnalysis = ShapeZZ();
    if (!passShapeAnalysis) continue;
    
    if (_systematic.Contains("invertvetotight") && _ntightlepton<3) continue;
    if (!_systematic.Contains("SS") && Lepton1.flavour * Lepton2.flavour > 0) continue;
    if (Lepton1.v.Pt() < 25.) continue;

    if (_systematic.Contains("lowMet") && MET.Et()>100.) continue;

    //ApplyStopTriggerEfficiency();

    // ISR reweighting
    //if (_applyisrreweighting && _isfastsim && (_isminitree || _systematic!="nominal")) ApplyISRReweighting();
    if (_applyisrreweighting && _isfastsim && (_isminitree || !_saveminitree)) ApplyISRReweighting();
 
    // DY corrections
    if (_applyDYcorrections && (_isminitree || !_saveminitree)) {
      float DYweight = 1.;
      // Latinos
      //float DYweight = (1.02852 - 0.0949640*TMath::Erf((_gen_ptll-19.0422)/10.4487) + 0.0758834*TMath::Erf((_gen_ptll-56.1146)/41.1653));
      // AN v4 	
      /*
      if (_njet>=1) {
	float x = fabs(dphill);
	if (MET.Et()<140.) DYweight = 1.03909 +0.414243*x -0.485282*x*x + 0.10607*x*x*x;
	else { 
	  if (_systematic.Contains("DYshapedo")) DYweight = 1.;
	  else if (_systematic.Contains("DYshapeup")) DYweight = 1.23218 +1.4524*x -1.3196*x*x + 0.253482*x*x*x;
	  else DYweight = 1.16803 +1.10048*x -1.13608*x*x + 0.23263*x*x*x;
	} 
      }
      // AN v6 
      if (MET.Et()>140.) {
	//float x = fabs(dphill);
	//float DYfit = 1.38577 + 3.36723*x -2.78402*x*x + 0.536200*x*x*x; 
	////float DYfit = 1.37517 + 2.17002*x -1.83012*x*x + 0.352897*x*x*x; 
	float x = _pt2l;
	float DYfit = 1.58986+0.368033*TMath::Erf((x+-318.317)/256.58)+0.367954*TMath::Erf((x+-318.316)/256.58);
	////float DYfit = 1.08702*(1.+0.00178185*x)/(1.-0.000300219*x); 
	if (_systematic.Contains("DYshapedo")) DYweight = 1.;
	else if (_systematic.Contains("DYshapeup")) DYweight = DYfit;
	else DYweight = (DYfit+1.)/2.;
      }
      */
      if (MET.Et()>100.) {
	// VR1 
	/*if (_MT2ll<20.) DYweight = 0.638916;
	else if (_MT2ll<40.) DYweight = 0.625556;
	else if (_MT2ll<60.) DYweight = 0.690579;
	else if (_MT2ll<80.) DYweight = 0.775666;
	else if (_MT2ll<100.) DYweight = 1.15934;
	else if (_MT2ll<120.) DYweight = 1.4206;
	else DYweight = 1.15749; */ 
	// VR1 norm
	if (_MT2ll<20.) DYweight = 0.904089;
	else if (_MT2ll<40.) DYweight = 0.885184;
	else if (_MT2ll<60.) DYweight = 0.977194;
	else if (_MT2ll<80.) DYweight = 1.09759;
	else if (_MT2ll<100.) DYweight = 1.64051;
	else if (_MT2ll<120.) DYweight = 2.0102;
	else DYweight = 1.63789;
	// SRs
	/*if (_MT2ll<20.) DYweight = 1.44905;
	else if (_MT2ll<40.) DYweight = 1.34401;
	else if (_MT2ll<60.) DYweight = 0.695682;
	else if (_MT2ll<80.) DYweight = 0.680000;
	else if (_MT2ll<100.) DYweight = 0.670264;
	else if (_MT2ll<120.) DYweight = 2.24221;
	else DYweight = 1.77837;*/
	if (_systematic.Contains("DYshapedo")) {
	  if (MassPointFlag=="_DYcorr") DYweight = 1.;
	  else DYweight = 1. - (DYweight-1.);
	} else if (_systematic.Contains("DYshapeup")) { 
	  if (MassPointFlag=="_DYcorr") DYweight = 1. + 2*(DYweight-1.);
	}
	if (MET.Et()<140. && MassPointFlag=="_DYcorr") DYweight *= 0.70;
      }
      _event_weight *= DYweight;
      _event_weight_Btagup     *= DYweight; 
      _event_weight_Btagdo     *= DYweight;
      _event_weight_Btaglightup*= DYweight; 
      _event_weight_Btaglightdo*= DYweight;
      _event_weight_BtagFSup   *= DYweight; 
      _event_weight_BtagFSdo   *= DYweight; 
      _event_weight_Idisoup    *= DYweight; 
      _event_weight_Idisodo    *= DYweight;  
      _event_weight_Idisoeleup *= DYweight; 
      _event_weight_Idisoeledo *= DYweight;  
      _event_weight_Idisomuup  *= DYweight; 
      _event_weight_Idisomudo  *= DYweight; 
      _event_weight_Triggerup  *= DYweight; 
      _event_weight_Triggerdo  *= DYweight; 
      _event_weight_Recoup     *= DYweight; 
      _event_weight_Recodo     *= DYweight; 
      _event_weight_Fastsimup  *= DYweight; 
      _event_weight_Fastsimdo  *= DYweight;
      _event_weight_Toppt      *= DYweight;
      _event_weight_Isrnjetup  *= DYweight;
      _event_weight_Isrnjetdo  *= DYweight;
    } 

    // For ZZ k-factors
    if ((_isminitree || !_saveminitree) && applyZZkfactors>-1) {
      float kFactor = 1.;
      if (applyZZkfactors==0) {
	TString Variation = "Nominal";
	if (_systematic=="PDFup") Variation = "PDFScaleDn";
	if (_systematic=="PDFdo") Variation = "PDFScaleUp";
	if (_systematic=="Q2up")  Variation = "QCDScaleDn";
	if (_systematic=="Q2do")  Variation = "QCDScaleUp";
	kFactor = kfactor_ggHZZ_qcd(_ZZmass, Variation);
      } else {
	int finalState = 2;
	if (applyZZkfactors==1) {
	  if (fabs(_lep1id)==fabs(_lep2id) && fabs(_lep1id)==fabs(_lep3id) && 
	      (_nlepton<4 || fabs(_lep1id)==fabs(_lep4id)) ) finalState = 1;
	}
	if (_systematic.Contains("kfPt"))   kFactor = kfactor_qqZZ_qcd_Pt  (_ZZpt,   finalState); 
	else if (_systematic.Contains("kfdPhi")) kFactor = kfactor_qqZZ_qcd_dPhi(_ZZdphi, finalState); 
	else if (!_systematic.Contains("kfNo"))    kFactor = kfactor_qqZZ_qcd_M   (_ZZmass, finalState); 
      }
      if (filename.Contains("ZZTo4L") && !filename.Contains("HToZZ")) kFactor *= 1.256/1.212;
      _event_weight *= kFactor;
      _event_weight_Btagup     *= kFactor; 
      _event_weight_Btagdo     *= kFactor;
      _event_weight_Btaglightup*= kFactor; 
      _event_weight_Btaglightdo*= kFactor;
      _event_weight_BtagFSup   *= kFactor; 
      _event_weight_BtagFSdo   *= kFactor; 
      _event_weight_Idisoup    *= kFactor; 
      _event_weight_Idisodo    *= kFactor;  
      _event_weight_Idisoeleup *= kFactor; 
      _event_weight_Idisoeledo *= kFactor;  
      _event_weight_Idisomuup  *= kFactor; 
      _event_weight_Idisomudo  *= kFactor; 
      _event_weight_Triggerup  *= kFactor; 
      _event_weight_Triggerdo  *= kFactor; 
      _event_weight_Recoup     *= kFactor; 
      _event_weight_Recodo     *= kFactor; 
      _event_weight_Fastsimup  *= kFactor; 
      _event_weight_Fastsimdo  *= kFactor;
      _event_weight_Toppt      *= kFactor;
      _event_weight_Isrnjetup  *= kFactor;
      _event_weight_Isrnjetdo  *= kFactor;
    }
    
    if (_systematic=="Topptdo")
      _event_weight = _event_weight_Toppt*_event_weight_Toppt/_event_weight;

    if (_ismc && !_isfastsim && _systematic.Contains("Pileup")) {
      if (_systematic=="Pileupup") _event_weight *= PileupSyst(iRunPeriod, trpu, +1);
      if (_systematic=="Pileupdo") _event_weight *= PileupSyst(iRunPeriod, trpu, -1);
    }

    // Fake uncertainty
    _event_weight_Fakeup = _event_weight;
    _event_weight_Fakedo = _event_weight;
    bool _applyFakeRateCorrection = _ismc && !_isfastsim;
    if (_systematic.Contains("SS") || _systematic.Contains("fake") || _systematic.Contains("kinematic")) _applyFakeRateCorrection = false; 
    if (_nLeptonsMatched<2) {
      float fakeRateUp = 1.5;
      float fakeRateDo = 0.5;
      if (_applyFakeRateCorrection && !_saveminitree) {
	cout << " Applying fake rescaling " << endl;
	float fakeRateCorrection = 1.135;
	fakeRateUp = 1.32;
	fakeRateDo = 0.95;
	_event_weight *= fakeRateCorrection;
	_event_weight_Btagup     *= fakeRateCorrection; 
	_event_weight_Btagdo     *= fakeRateCorrection;
	_event_weight_Btaglightup*= fakeRateCorrection; 
	_event_weight_Btaglightdo*= fakeRateCorrection;
	_event_weight_BtagFSup   *= fakeRateCorrection; 
	_event_weight_BtagFSdo   *= fakeRateCorrection; 
	_event_weight_Idisoup    *= fakeRateCorrection; 
	_event_weight_Idisodo    *= fakeRateCorrection;  
	_event_weight_Idisoeleup *= fakeRateCorrection; 
	_event_weight_Idisoeledo *= fakeRateCorrection;  
	_event_weight_Idisomuup  *= fakeRateCorrection; 
	_event_weight_Idisomudo  *= fakeRateCorrection; 
	_event_weight_Triggerup  *= fakeRateCorrection; 
	_event_weight_Triggerdo  *= fakeRateCorrection; 
	_event_weight_Recoup     *= fakeRateCorrection; 
	_event_weight_Recodo     *= fakeRateCorrection;
	_event_weight_Fastsimup  *= fakeRateCorrection; 
	_event_weight_Fastsimdo  *= fakeRateCorrection;
	_event_weight_Toppt      *= fakeRateCorrection;
	_event_weight_Isrnjetup  *= fakeRateCorrection;
	_event_weight_Isrnjetdo  *= fakeRateCorrection;
      }
      _event_weight_Fakeup *= fakeRateUp;
      _event_weight_Fakedo *= fakeRateDo;
    }
    if (_systematic=="Fakeup") _event_weight = _event_weight_Fakeup;
    if (_systematic=="Fakedo") _event_weight = _event_weight_Fakedo;
    
    if (_systematic.Contains("MetCorr") && _ismc) {
      float xmet = (MET.Et()<400.) ? MET.Et() : 400.;
      if (_leadingPtCSVv2M>=20.) 
	_event_weight *= (0.925552+0.000384281*xmet-1.58635e-06*xmet*xmet);
      else if (_systematic.Contains("Zpeak"))
	_event_weight *= (1.19644-0.00115806*xmet+1.95527e-06*xmet*xmet);
      else if (_systematic.Contains("WZtoWW")) 
	//_event_weight *= (1.04613-0.000528796*xmet-2.12198e-07*xmet*xmet);
	_event_weight *= (1.08594-0.000319475*xmet-1.39091e-06*xmet*xmet);
    } 

    if (_systematic.Contains("noweight")) _event_weight = 1.;
      
    if (_isfastsim && _DoTheoreticalVariations==1) {
      hasLHE = false;
      for (int il = 0; il<_lheEvent.size(); il++) {
	if (_lheEvent.at(il)==event) {
	  hasLHE = true;
	  lheTree->GetEntry(il); // This fills std_vector_LHE_weight for FastSim
	  break;
	}
      }
    }
      
    // Fill histograms
    //--------------------------------------------------------------------------
    bool pass = true;
    bool pass_blind = true; 

    // Blinding policy: Just CR with 10  /fb.  For us blinded () = Met < 140; 
    if (filename.Contains("Data") || filename.Contains("PromptReco") || filename.Contains("23Sep2016") || _filename.Contains("03Feb2017")) 
      {
	
	if (_systematic.Contains("Blind")) pass_blind = false;
	//if (MET.Et()<140.) pass_blind = true;
	if (run < 276502) pass_blind = true;
      }

    if (_systematic.Contains("BlindMT2") && _mt2ll>40.) pass_blind = false;

    // Fill histograms
    // -----------------------------------------------------------------------------
    
    // Basics Stop
    //-------------------------------------------------------------------------
    FillLevelHistograms(Stop_00_Has2Leptons, !_isminitree && _m2l>20. && pass && pass_blind && pass_masspoint);

    pass &= _m2l>20.;

    //FillLevelHistograms(Stop_00_mll20, pass && !_isminitree && pass_blind && pass_masspoint); 
   
    if (_systematic.Contains("ZZlepton")) {
      if (_leadingPtCSVv2M >= 20.) {
	FillLevelHistograms(Stop_ZZ_Tag_loose, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_Tag, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_Tag_tight, pass && pass_blind && pass_masspoint);
	if (MET.Et()>=140.) {
	  FillLevelHistograms(Stop_ZZ_SRs_Tag_loose, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_SRs_Tag, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_SRs_Tag_tight, pass && pass_blind && pass_masspoint);
	}
      } else if (_njet<1) {
	FillLevelHistograms(Stop_ZZ_NoJet_loose, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_NoJet, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_NoJet_tight, pass && pass_blind && pass_masspoint);
	if (MET.Et()>=140.) {
	  FillLevelHistograms(Stop_ZZ_SRs_NoJet_loose, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_SRs_NoJet, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_SRs_NoJet_tight, pass && pass_blind && pass_masspoint);
	}
      } else {
	FillLevelHistograms(Stop_ZZ_NoTag_loose, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_NoTag, pass && pass_blind && pass_masspoint);
	if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_NoTag_tight, pass && pass_blind && pass_masspoint);
	if (MET.Et()>=140.) {
	  FillLevelHistograms(Stop_ZZ_SRs_NoTag_loose, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=3) FillLevelHistograms(Stop_ZZ_SRs_NoTag, pass && pass_blind && pass_masspoint);
	  if (_ntightlepton>=4) FillLevelHistograms(Stop_ZZ_SRs_NoTag_tight, pass && pass_blind && pass_masspoint);
	}
      }
      continue;
    }
    
    // Look in the Z-peak
    bool Zpeak = pass && _channel!=em && fabs(_m2l - Z_MASS)<15. ;  
    if (_systematic.Contains("Zpeak") && !Zpeak) continue;
    
    /*if (MET.Et()>140. && _njet==0 && dphill>3.2-0.32) {
      cout << "SSSLLL " << _nLeptonsMatched << " M2l " << _m2l << "  dphill " << 180./3.14159*dphill << "  dphilmet " <<  180./3.14159*_dphiminlmet << "  dphillmet " << 180./3.14159*_dphillmet << " mt2 " << _MT2ll << endl;
      cout << " event = " << run << " " << event << "  " << ZDecayType << " " << _htgen << endl;
      cout << " Lep1 " << Lepton1.v.Pt() << " " << Lepton1.v.Eta() << " " << 180./3.14159*Lepton1.v.Phi() << " " << Lepton1.flavour << endl;
      cout << " Lep2 " << Lepton2.v.Pt() << " " << Lepton2.v.Eta() << " " << 180./3.14159*Lepton2.v.Phi() << " " << Lepton2.flavour << endl;
      cout << " MET " << MET.Et() << " " << 180./3.14159*MET.Phi() << endl;
      cout << " LepPt " << (Lepton1.v+Lepton2.v).Pt() << endl;
      TVector3 JJJ; JJJ.SetXYZ(0.,0.,0.);
      cout << " jet" << endl;
      if (_njet==0) {
	for (int i=0; i<std_vector_jet_pt->size(); i++) 
	  if (std_vector_jet_pt ->at(i)>=0.) {
	    TVector3 JJT; JJT.SetPtEtaPhi(std_vector_jet_pt->at(i), std_vector_jet_eta->at(i), std_vector_jet_phi->at(i));
	    JJJ += JJT;
	    cout << i << "    " << std_vector_jet_pt ->at(i) << " " << std_vector_jet_eta ->at(i) << " " << 180./3.14159*std_vector_jet_phi ->at(i) << "    " << 180./3.14159*(MET.Phi()-std_vector_jet_phi ->at(i)) << endl;
	  }
      }
      cout << " JetPT " << JJJ.Pt() << " " << JJJ.Phi() << endl;
      cout << " lep " << endl;
      for (int i=0; i<std_vector_lepton_pt->size(); i++) {
	if (std_vector_lepton_pt->at(i)>=0.)
	  cout << "    " << std_vector_lepton_pt->at(i) << " " << std_vector_lepton_eta->at(i) << " " << std_vector_lepton_phi->at(i) << " " << std_vector_lepton_flavour->at(i) << " " << std_vector_lepton_isLooseLepton->at(i) << " " << std_vector_lepton_isTightLepton->at(i) << endl;
	  }
      } continue;*/

    //bool Zveto = ( _channel == em || fabs(_m2l - Z_MASS) > 15. ); 
    //if (_systematic.Contains("Zveto") && !Zveto) continue;

    if (!_isminitree) { 

      // Zpeak level
      FillLevelHistograms(Stop_05_Zpeak,      Zpeak && pass_blind && pass_masspoint); // 2 OS Leptons, mll > 20, Z peak
      //FillLevelHistograms(Stop_05_NoTagZpeak, Zpeak && (_leadingPtCSVv2M <  20.) && pass_masspoint); // 2 OS Leptons, mll > 20, Z peak + 0 b Tag (VET0) 
      //FillLevelHistograms(Stop_05_TagZpeak,   Zpeak && (_leadingPtCSVv2M >= 20.) && pass_masspoint); // 2 OS Leptons, mll > 20, Z peak + 1 b Tag         
      
      // DY estimation -> Routin Level
      // ---------------------------------------------------------------
      //FillLevelHistograms(Stop_04_Routin,      pass && pass_masspoint);  // 2 OS Leptons, mll > 20 
      //FillLevelHistograms(Stop_04_Jet2Routin,  pass &&  jetpt2 >= 30. && pass_masspoint);  // 2 OS Leptons, mll > 20,  2 Jets
      //FillLevelHistograms(Stop_04_TagRoutin,   pass && (_leadingPtCSVv2M >= 20.) && pass_masspoint);  // 2 OS Leptons, mll > 20, 1 b Tag
      //FillLevelHistograms(Stop_04_NoTagRoutin, pass && (_leadingPtCSVv2M <  20.) && pass_masspoint);  // 2 OS Leptons, mll > 20, 0 b Tag (VET0)
      
      bool Zveto = ( _channel == em || fabs(_m2l - Z_MASS) > 15. ); 
      bool _isZSelection = _systematic.Contains("WZ")    || _systematic.Contains("ZWtoZ") || 
	_systematic.Contains("Zpeak") || _systematic.Contains("ZZ");
      pass &= Zveto || _isZSelection;
    
      // Leave this line at the end of this if or the results on latino trees and minitrees will be inconsistent
      if (pass && _saveminitree && (MET.Et()>=100. || !_systematic.Contains("Zpeak"))) minitree->Fill();

    }

    FillLevelHistograms(Stop_00_Zveto, pass && pass_blind && pass_masspoint);

    if (_systematic.Contains("SS")) {
      if (pass && pass_blind && MET.Et()>=50.) {
	FillLevelHistograms(Stop_00_SS_Tag,  _leadingPtCSVv2M >= 30.);
	FillLevelHistograms(Stop_00_SS_2Jet, jetpt2>=30.);
	FillLevelHistograms(Stop_00_SS,      _leadingPtCSVv2M >= 30. && jetpt2>=30.);
      } 
    }

    //~~~~~~~~~~~~~~~~~~~~~~ MC normalization study ~~~~~~~~~

    bool WW = pass && _njet == 0;
    FillLevelHistograms(Stop_00_WWsel, WW && pass_blind && pass_masspoint);
    FillLevelHistograms(Stop_00_WWselMET, WW && MET.Et()>= 50 && pass_blind && pass_masspoint);    
    bool TTbar = pass && _njet >1 && _leadingPtCSVv2M >= 20.; 
    FillLevelHistograms(Stop_00_TTsel, TTbar && pass_blind && pass_masspoint);
    FillLevelHistograms(Stop_00_TTselMET, TTbar && MET.Et()>= 50 && pass_blind && pass_masspoint);
      
   //~~~~~~~~~~~~~~~~~~~~~~ Analysis Cuts ~~~~~~~~~~~~~~~~~~ 

   // Tag SELECTION -> Bin0Tag & Bin1Tag;  used in minitrees and latino trees

    FillLevelHistograms(Stop_01_Tag,       pass && (_leadingPtCSVv2M >= 20.) && pass_blind && pass_masspoint);
    FillLevelHistograms(Stop_01_NoTag,     pass && (_leadingPtCSVv2M <  20.) && pass_blind && pass_masspoint);
    
    if (_leadingPtCSVv2M >= 20.) {
      FillLevelHistograms(Stop_02_VR1_Tag,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_blind && pass_masspoint);
      if (_leadingPtCSVv2T >= 30. && njet>1)
	FillLevelHistograms(Stop_02_VR1_Tag2Jet,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1_Tag,   pass && (MET.Et()>=140. && MET.Et()<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2_Tag,   pass && (MET.Et()>=200. && MET.Et()<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3_Tag,   pass && (MET.Et()>=300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1gen_Tag, pass && (metGenpt>=140. && metGenpt<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2gen_Tag, pass && (metGenpt>=200. && metGenpt<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3gen_Tag, pass && (metGenpt>=300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SRs_Tag,   pass && (MET.Et()>=140.) && pass_blind && pass_masspoint);
    }
    
    if (_leadingPtCSVv2M <  20.) {
      FillLevelHistograms(Stop_02_VR1_NoTag,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_blind && pass_masspoint);
      if (njet<1) FillLevelHistograms(Stop_02_VR1_NoJet,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_masspoint);
    }

    if (_njet<1) {
      //FillLevelHistograms(Stop_02_VR1_NoJet,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1_NoJet,   pass && (MET.Et()>=140. && MET.Et()<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2_NoJet,   pass && (MET.Et()>=200. && MET.Et()<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3_NoJet,   pass && (MET.Et()>=300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1gen_NoJet, pass && (metGenpt>=140. && metGenpt<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2gen_NoJet, pass && (metGenpt>=200. && metGenpt<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3gen_NoJet, pass && (metGenpt>=300.) && pass_blind && pass_masspoint);
      //} 
    } else if (_leadingPtCSVv2M <  20.) {
      //FillLevelHistograms(Stop_02_VR1_NoTag,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_blind && pass_masspoint);
      //if (_njet<1)
      //FillLevelHistograms(Stop_02_VR1_NoJet,   pass && (MET.Et()>=100. && MET.Et()<140.) && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1_NoTag,   pass && (MET.Et()>=140. && MET.Et()<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2_NoTag,   pass && (MET.Et()>=200. && MET.Et()<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3_NoTag,   pass && (MET.Et()>=300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR1gen_NoTag, pass && (metGenpt>=140. && metGenpt<200.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR2gen_NoTag, pass && (metGenpt>=200. && metGenpt<300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SR3gen_NoTag, pass && (metGenpt>=300.) && pass_blind && pass_masspoint);
      FillLevelHistograms(Stop_02_SRs_NoTag,   pass && (MET.Et()>=140.) && pass_blind && pass_masspoint);
    }
    
    FillLevelHistograms(Stop_02_SRs, pass && (MET.Et()>=140.) && pass_blind && pass_masspoint);
    
  }

  if (_SaveHistograms>=3) SaveSystematicHistograms();
  if (_DoTheoreticalVariations==1) SaveTheoreticalVariationsHistograms();
 
  EndJob();
}

//------------------------------------------------------------------------------
// BookAnalysisHistograms
//------------------------------------------------------------------------------
void AnalysisStop::BookAnalysisHistograms()
{

  TH1::SetDefaultSumw2();

  for (int j=0; j<ncut; j++) {
  
    if (_isminitree && j<Stop_00_SS_Tag && j!=Stop_00_VZ && j!=Stop_00_ttZ_Tag && j!=Stop_00_ttZ_2Tag) continue;

    if (!_systematic.Contains("SS") && j>=Stop_00_SS_Tag && j<=Stop_00_SS) continue;
      
    if (_systematic.Contains("WZ3L") && j!=Stop_00_VZ && !_systematic.Contains("ttZ")) continue; 
    if (_systematic.Contains("ttZ") && j!=Stop_00_ttZ_Tag && j!=Stop_00_ttZ_2Tag && !_systematic.Contains("WZ3L")) continue; 
    if ((_systematic.Contains("WZ3L") || _systematic.Contains("ttZ")) && 
	j!=Stop_00_VZ && j!=Stop_00_ttZ_Tag && j!=Stop_00_ttZ_2Tag) continue;

    if (_systematic.Contains("ZZlepton") && !scut[j].Contains("/ZZ_")) continue;
    if (!_systematic.Contains("ZZlepton") && scut[j].Contains("/ZZ_")) continue;

    for (int k=0; k<=njetbin; k++) {

      if ((_SaveHistograms>=2 || _systematic!="nominal") && k<njetbin) continue; 

      TString sbin = (k < njetbin) ? Form("/%djet", k) : "";

      TString directory = scut[j] + sbin;

      root_output->cd();

      if (_SaveHistograms>=2 || _systematic!="nominal" || k<njetbin) gDirectory->mkdir(directory);

      root_output->cd(directory);

      for (int i=ee; i<=ll; i++) {

	TString suffix = "_" + schannel[i];
	
	if (_SaveHistograms==0) DefineHistograms(i, j, k, suffix);

	h_MT2ll             [i][j][k] = new TH1F("h_MT2ll"              + suffix, "",    7,    0,  140);
	h_MT2llgen          [i][j][k] = new TH1F("h_MT2llgen"           + suffix, "",    7,    0,  140);
	h_MT2llisr          [i][j][k] = new TH1F("h_MT2llisr"           + suffix, "",    7,    0,  140);
	h_MT2llisrgen       [i][j][k] = new TH1F("h_MT2llisrgen"        + suffix, "",    7,    0,  140);

	if (_SaveHistograms>=2 && _systematic=="nominal") {

	  h_MT2ll_nvtxup      [i][j][k] = new TH1F("h_MT2ll_nvtxup"       + suffix, "",    7,    0,  140);
	  h_MT2ll_nvtxdo      [i][j][k] = new TH1F("h_MT2ll_nvtxdo"       + suffix, "",    7,    0,  140);
	  h_MT2llgen_nvtxup   [i][j][k] = new TH1F("h_MT2llgen_nvtxup"    + suffix, "",    7,    0,  140);
	  h_MT2llgen_nvtxdo   [i][j][k] = new TH1F("h_MT2llgen_nvtxdo"    + suffix, "",    7,    0,  140);
	  h_MT2llisr_nvtxup   [i][j][k] = new TH1F("h_MT2llisr_nvtxup"    + suffix, "",    7,    0,  140);
	  h_MT2llisr_nvtxdo   [i][j][k] = new TH1F("h_MT2llisr_nvtxdo"    + suffix, "",    7,    0,  140);
	  h_MT2llisrgen_nvtxup[i][j][k] = new TH1F("h_MT2llisrgen_nvtxup" + suffix, "",    7,    0,  140);
	  h_MT2llisrgen_nvtxdo[i][j][k] = new TH1F("h_MT2llisrgen_nvtxdo" + suffix, "",    7,    0,  140);
	  
	}

	if (_SaveHistograms>=2) continue;

	if (_systematic.Contains("Zpeak")) {
	  h_HT                [i][j][k] = new TH1D("h_HT"                + suffix, "",   80,    0,  800);
	  h_MT2ll_HTm150      [i][j][k] = new TH1D("h_MT2ll_HTm150"      + suffix, "",    7,    0,  140);
	  h_MT2ll_HTp150      [i][j][k] = new TH1D("h_MT2ll_HTp150"      + suffix, "",    7,    0,  140);
	  h_MT2ll_HTm200      [i][j][k] = new TH1D("h_MT2ll_HTm200"      + suffix, "",    7,    0,  140);
	  h_MT2ll_HTp200      [i][j][k] = new TH1D("h_MT2ll_HTp200"      + suffix, "",    7,    0,  140);
	  h_njet20dphilmet    [i][j][k] = new TH1D("h_njet20dphilmet"    + suffix, "",   10,    0,   10);
	} else if (_systematic.Contains("ZZ") || _systematic.Contains("WZ")) {
	  h_M1ll              [i][j][k] = new TH1D("h_M1ll"              + suffix, "",   60,  60., 120.);
	  h_M2ll              [i][j][k] = new TH1D("h_M2ll"              + suffix, "",  180,   0., 180.);
	} 

	h_dphillMET         [i][j][k] = new TH1D("h_dphillMET"          + suffix, "",  100,    0,  3.2);
	h_ptdphiLL          [i][j][k] = new TH2D("h_ptdphiLL"           + suffix, "",   80,    0,  800, 10,  0,  3.2);
	h_dphiLL            [i][j][k] = new TH1D("h_dphiLL"             + suffix, "",  100,    0,  3.2);
	float ptLLbins[24] = {0, 10., 20., 30., 40., 50., 60, 70., 80., 90., 100., 110., 120., 130., 140., // 15
			    160., 180., 200., 250., 300., 400., 500., 600., 800.}; // 9
	h_ptLL              [i][j][k] = new TH1D("h_ptLL"               + suffix, "",   80,    0,  800);
	h_ptLLbins          [i][j][k] = new TH1D("h_ptLLbins"           + suffix, "",   23, ptLLbins);
	h_genptLL           [i][j][k] = new TH1D("h_genptLL"            + suffix, "",   80,    0,  800);
	h_MET               [i][j][k] = new TH1D("h_MET"                + suffix, "",   80,    0,  800);
	h_Counter           [i][j][k] = new TH1D("h_Counter"            + suffix, "",    1,   80,  100);
	h_njet20            [i][j][k] = new TH1D("h_njet20"             + suffix, "",    7, -0.5,  6.5);
	h_nbjet             [i][j][k] = new TH1D("h_nbjet"              + suffix, "",    7, -0.5,  6.5);
	h_njet30            [i][j][k] = new TH1D("h_njet30"             + suffix, "",    7, -0.5,  6.5);	
	h_Lep1Pt            [i][j][k] = new TH1D("h_Lep1Pt"             + suffix, "",  2000,    0, 2000);
	h_Lep2Pt            [i][j][k] = new TH1D("h_Lep2Pt"             + suffix, "",  2000,    0, 2000);
	h_Lep1Phi           [i][j][k] = new TH1D("h_Lep2Phi"            + suffix, "",  200, -3.2,  3.2);
	h_Lep2Eta           [i][j][k] = new TH1D("h_Lep1Eta"            + suffix, "",   60,  -3.,   3.);
	h_Lep1Eta           [i][j][k] = new TH1D("h_Lep2Eta"            + suffix, "",   60,  -3.,   3.);
	h_Lep2Phi           [i][j][k] = new TH1D("h_Lep1Phi"            + suffix, "",  200, -3.2,  3.2);
 	h_dphil1MET         [i][j][k] = new TH1D("h_dphil1MET"          + suffix, "",  100,    0,  3.2);
 	h_dphil2MET         [i][j][k] = new TH1D("h_dphil2MET"          + suffix, "",  100,    0,  3.2);
        h_METphi            [i][j][k] = new TH1D("h_METphi"             + suffix, "",  200, -3.2,  3.2);

	if (_systematic.Contains("DYcorrections")) {
	  h_dphiLLbin1            [i][j][k] = new TH1D("h_dphiLLbin1"   + suffix, "",   80,    0, 3.14159);
	  h_dphiLLbin2            [i][j][k] = new TH1D("h_dphiLLbin2"   + suffix, "",   80,    0, 3.14159);
	  h_dphiLLbin3            [i][j][k] = new TH1D("h_dphiLLbin3"   + suffix, "",   80,    0, 3.14159);
	  h_dphiLLbin4            [i][j][k] = new TH1D("h_dphiLLbin4"   + suffix, "",   80,    0, 3.14159);
	  h_dphiLLbin5            [i][j][k] = new TH1D("h_dphiLLbin5"   + suffix, "",   80,    0, 3.14159);
	  continue;
	}

	if (_systematic.Contains("ZZ")) continue;
	if (_systematic.Contains("WZ")) continue;
	if (_systematic.Contains("ZWtoZ")) continue;

	h_dphiminlepmet     [i][j][k] = new TH1D("h_dphiminlepmet"      + suffix, "",  100,    0,  3.2);

	if (_systematic.Contains("Zpeak")) continue;

	h_MET_fake          [i][j][k] = new TH1F("h_MET_fake"           + suffix, "",  100,    0, 1000);
	h_MET_truth         [i][j][k] = new TH1F("h_MET_truth"          + suffix, "",  100,    0, 1000);
	h_MT2ll_fake        [i][j][k] = new TH1F("h_MT2ll_fake"         + suffix, "",    7,    0,  140);
	h_MT2ll_truth       [i][j][k] = new TH1F("h_MT2ll_truth"        + suffix, "",    7,    0,  140);

	if (!_systematic.Contains("kinematic")) continue;	

	h_Jet1Pt            [i][j][k] = new TH1D("h_Jet1Pt"             + suffix, "",  2000,    0, 2000);
	h_Jet2Pt            [i][j][k] = new TH1D("h_Jet2Pt"             + suffix, "",  2000,    0, 2000);
	h_JetPt             [i][j][k] = new TH1D("h_JetPt"              + suffix, "",  2000,    0, 2000);

	h_MT2ll_MET         [i][j][k] = new TH2F("h_MT2ll_MET"          + suffix, "",   80,    0,  800,    7,    0,  140);

	h_dphiisrmet        [i][j][k] = new TH1F("h_dphiisrmet"         + suffix, "",  128,    0,  3.2);

	h_mt2lblbcomb       [i][j][k] = new TH1D("h_mt2lblbcomb"        + suffix, "", 3000,    0, 3000);
	h_mt2bbtrue         [i][j][k] = new TH1D("h_mt2bbtrue"          + suffix, "", 3000,    0, 3000);
	h_mt2lblbtrue       [i][j][k] = new TH1D("h_mt2lblbtrue"        + suffix, "", 3000,    0, 3000);
	h_mt2lblbmatch      [i][j][k] = new TH1D("h_mt2lblbmatch"       + suffix, "", 3000,    0, 3000);
	h_mlb1comb          [i][j][k] = new TH1D("h_mlb1comb"           + suffix, "", 3000,    0, 3000);
	h_mlb2comb          [i][j][k] = new TH1D("h_mlb2comb"           + suffix, "", 3000,    0, 3000);
	h_mlb1true          [i][j][k] = new TH1D("h_mlb1true"           + suffix, "", 3000,    0, 3000);
	h_mlb2true          [i][j][k] = new TH1D("h_mlb2true"           + suffix, "", 3000,    0, 3000);
	h_mt2lblbvsmlbtrue  [i][j][k] = new TH2D("h_mt2lblbvsmlbtrue"   + suffix, "",  100,    0, 1000,  100,    0, 1000);
	h_nisrjet           [i][j][k] = new TH1D("h_nisrjet"            + suffix, "",   10,    0,   10);

	h_metmeff           [i][j][k] = new TH1D("h_metmeff"            + suffix, "",  500,    0,    5);

	
	h_MT2_Met           [i][j][k] = new TH1D("h_MT2_Met" + suffix, "", NbinsMT2*NbinsMet, vMinMT2, vMinMT2 + NbinsMet*(vMaxMT2-vMinMT2));
	h_HTvisible_Met     [i][j][k] = new TH1D("h_HTvisible_Met" + suffix, "", NbinsHTvisible*NbinsMet, vMinHTvisible, vMinHTvisible + 
NbinsMet*(vMaxHTvisible-vMinHTvisible));
	h_MetMeff_Met       [i][j][k] = new TH1D("h_MetMeff_Met" + suffix, "", NbinsMetMeff*NbinsMet, vMinMetMeff, vMinMetMeff + NbinsMet*(vMaxMetMeff-vMinMetMeff));
	h_R2_Met            [i][j][k] = new TH1D("h_R2_Met" + suffix, "", NbinsR2*NbinsMet, vMinR2, vMinR2 + NbinsMet*(vMaxR2-vMinR2));

      }
    }
  }

}


//------------------------------------------------------------------------------
// BookSystematicHistograms
//------------------------------------------------------------------------------
void AnalysisStop::BookSystematicHistograms()
{

  TH1::SetDefaultSumw2();

  for (int is = 0; is<nsystematic; is++) {
 
    if (ssystematic[is]=="nominal" || !systematicfromweight[is]) continue;
    if (ssystematic[is].Contains("Isrnjet") && !_isfastsim) continue;
    
    TString prefix = (_isminitree) ? "minitrees/" : "";
    gSystem->mkdir(prefix + "rootfiles/" + ssystematic[is] + "/" + _analysis, kTRUE);
    TString _systematiclongname = _longname;
    _systematiclongname.ReplaceAll("nominal/", ssystematic[is] + "/");
    root_output_systematic[is] = new TFile(prefix + "rootfiles/" + _systematiclongname + ".root", "recreate");

    root_output_systematic[is]->cd();

    for (int j=0; j<ncut; j++) {
  
      if (_isminitree && j<Stop_00_Zveto) continue; 
      
      if (_systematic.Contains("WZ3L") && j!=Stop_00_VZ && !_systematic.Contains("ttZ")) continue; 
      if (_systematic.Contains("ttZ") && j!=Stop_00_ttZ_Tag && j!=Stop_00_ttZ_2Tag && !_systematic.Contains("WZ3L")) continue; 
      if ((_systematic.Contains("WZ3L") || _systematic.Contains("ttZ")) && 
	j!=Stop_00_VZ && j!=Stop_00_ttZ_Tag && j!=Stop_00_ttZ_2Tag) continue;

      if (_systematic.Contains("ZZ") && !scut[j].Contains("/ZZ_")) continue;
      if (!_systematic.Contains("ZZ") && scut[j].Contains("/ZZ_")) continue;

      TString directory = scut[j];

      root_output_systematic[is]->cd();

      gDirectory->mkdir(directory);

      root_output_systematic[is]->cd(directory);
      
      for (int i=ee; i<=ll; i++) {

	TString suffix = "_" + schannel[i];

	h_MT2ll_systematic        [i][j][is] = new TH1F("h_MT2ll"            + suffix, "",    7,    0,  140);
	h_MT2llgen_systematic     [i][j][is] = new TH1F("h_MT2llgen"         + suffix, "",    7,    0,  140);
	h_MT2llisr_systematic     [i][j][is] = new TH1F("h_MT2llisr"         + suffix, "",    7,    0,  140);
	h_MT2llisrgen_systematic  [i][j][is] = new TH1F("h_MT2llisrgen"      + suffix, "",    7,    0,  140);
	
      }

    }

  }

}

//------------------------------------------------------------------------------
// BookTheoreticalVariationsHistograms
//------------------------------------------------------------------------------
void AnalysisStop::BookTheoreticalVariationsHistograms()
{
  
  // Theoretical variation renormalizations
  TH1F *ListVectorWeights;
  
  if (!_isfastsim) {

    TFile* file = TFile::Open(_filename);
    ListVectorWeights = (TH1F*)file->Get("list_vectors_weights");

  } else if (_isminitree) {
    
    _lheEvent.clear();

    TString _lhefilename;

    TString _longnameroot = _longname + ".root";
    TString tok;
    Ssiz_t from = 0;
    const char* delim = "/";
    
    while (_longnameroot.Tokenize(tok, from, delim)) {
      if (tok.Contains(".root")) {
	_lhefilename = tok.ReplaceAll(".root", "");
      }
    }

    TFile* file = TFile::Open("/eos/user/s/scodella/Stop/LatinoTrees/latino_" + _lhefilename + ".root");

    if (file) {

      ListVectorWeights = (TH1F*)file->Get("SumGenWeights");

      lheTree = (TTree*) file->Get("latino");

      UInt_t _lheevent;
      lheTree->SetBranchAddress("event", &_lheevent);
      lheTree->SetBranchAddress("std_vector_GEN_weight", &std_vector_LHE_weight);

      int nlheentries = lheTree->GetEntries();
      for (int il = 0; il<nlheentries; il++) {

	lheTree->GetEntry(il);
	_lheEvent.push_back(_lheevent);
	
      }

    }

  }

  if (!ListVectorWeights) { _DoTheoreticalVariations = 0; return; }

  int _lheoffset = _isfastsim ? 1 : 0;
  for (int tv = 0; tv<nTheoreticalVariations; tv++) 
    _TheoreticalVariationRenormalization[tv] = ListVectorWeights->GetBinContent(tv+1+_lheoffset)/ListVectorWeights->GetBinContent(1+_lheoffset);

  TH1::SetDefaultSumw2();
    
  TString prefix = (_isminitree) ? "minitrees/" : "";
  gSystem->mkdir(prefix + "rootfiles/theory/" + _analysis, kTRUE);
  TString _systematiclongname = _longname;
  _systematiclongname.ReplaceAll("nominal/", "theory/");
  root_output_theoreticalvariations = new TFile(prefix + "rootfiles/" + _systematiclongname + ".root", "recreate");

  root_output_theoreticalvariations->cd();

  for (int j=0; j<ncut; j++) {
  
    if (!scut[j].Contains("SR") && !scut[j].Contains("VR")) continue;
    if (scut[j].Contains("/ZZ_")) continue;

    TString directory = scut[j];
      
    root_output_theoreticalvariations->cd();
    
    gDirectory->mkdir(directory);
    
    root_output_theoreticalvariations->cd(directory);
    
    for (int i=ee; i<=ll; i++) {
      
      for (int tv = 0; tv<nTheoreticalVariations; tv++) {
	
	TString suffix = "_" + schannel[i] + "_TV"; suffix += tv;
	
	h_MT2ll_theoreticalvariation        [i][j][tv] = new TH1F("h_MT2ll"            + suffix, "",    7,    0,  140);
	h_MT2llgen_theoreticalvariation     [i][j][tv] = new TH1F("h_MT2llgen"         + suffix, "",    7,    0,  140);
	h_MT2llisr_theoreticalvariation     [i][j][tv] = new TH1F("h_MT2llisr"         + suffix, "",    7,    0,  140);
	h_MT2llisrgen_theoreticalvariation  [i][j][tv] = new TH1F("h_MT2llisrgen"      + suffix, "",    7,    0,  140);
	
      }
      
    }

  }

}

//------------------------------------------------------------------------------
// GetAnalysisVariables
//------------------------------------------------------------------------------
void AnalysisStop::GetAnalysisVariables()
{

  if (_isminitree) { 

    // Met
    MET.SetPtEtaPhiM(metPfType1, 0.0, metPfType1Phi, 0.0); 
    _jetbin = (_njet < njetbin) ? _njet : njetbin - 1;

    // Leptons
    AnalysisLeptons.clear();
    
    Lepton lep; TLorentzVector tlv;
    lep.type      = 1;
    lep.flavour   = _lep1id;
    tlv.SetPtEtaPhiM(_lep1pt, _lep1eta, _lep1phi, _lep1mass);
    lep.v = tlv;
    AnalysisLeptons.push_back(lep);
    lep.type      = 1;
    lep.flavour   = _lep2id;
    tlv.SetPtEtaPhiM(_lep2pt, _lep2eta, _lep2phi, _lep2mass);
    lep.v = tlv;
    AnalysisLeptons.push_back(lep);
    if (_nlepton>=3) {
      lep.type      = _lep3type;
      lep.flavour   = _lep3id;
      lep.idisoW    = _lep3idisoW;
      tlv.SetPtEtaPhiM(_lep3pt, _lep3eta, _lep3phi, _lep3mass);
      lep.v = tlv;
      AnalysisLeptons.push_back(lep);
      if (_nlepton>=4) {
	lep.type      = _lep4type;
	lep.flavour   = _lep4id;
	lep.idisoW    = _lep4idisoW;
	tlv.SetPtEtaPhiM(_lep4pt, _lep4eta, _lep4phi, _lep4mass);
	lep.v = tlv;
	AnalysisLeptons.push_back(lep);
      }
    }
    Lepton1 = AnalysisLeptons[0];
    Lepton2 = AnalysisLeptons[1];
    
    _pt2l = (Lepton1.v+Lepton2.v).Pt();

  }

  // Jets
  _jet_pt.clear();
  _jet_eta.clear();
  _jet_phi.clear();
  
  for (int i=0; i<std_vector_jet_pt->size(); i++) {

    float pt  = std_vector_jet_pt ->at(i);
    if (pt < 20.) continue;

    float eta = std_vector_jet_eta->at(i);
    if (abs(eta) > 2.4) continue;

    float phi = std_vector_jet_phi->at(i);

    _jet_eta.push_back(eta); 
    _jet_phi.push_back(phi);
    _jet_pt .push_back(pt); 
    
  }
  
  _metmeff = MET.Et()/_meff;
  _MT2ll = (_mt2ll<140.) ? _mt2ll : 139.;
  _MT2llgen = (_mt2llgen<140.) ? _mt2llgen : 139.;
  _MT2llfake = (_mt2ll<160.) ? _mt2ll : 159.;

  // NbinsMet;
  //float tempMet = MET.Et(); if (tempMet<vMinMet) tempMet = vMinMet; if (tempMet>=vMaxMet) tempMet = vMaxMet - 0.1;
  //float widthMetBin = (vMaxMet - vMinMet)/NbinsMet;
  int binMet;// = (tempMet - vMinMet)/widthMetBin;
  if (MET.Et()<100.) binMet = 0;
  else if (MET.Et()<140.) binMet = 1;
  else if (MET.Et()<200.) binMet = 2;
  else if (MET.Et()<300.) binMet = 3;
  else binMet = 4;

  // _MT2_Met
  float tempMT2 = _mt2ll; if (tempMT2<vMinMT2) tempMT2 = vMinMT2; if (tempMT2>=vMaxMT2) tempMT2 = vMaxMT2 - 0.1;
  _MT2_Met = tempMT2 + (vMaxMT2 - vMinMT2)*binMet;

  // _HTvisible_Met
  float tempHTvisible = _htvisible; if (tempHTvisible<vMinHTvisible) tempHTvisible = vMinHTvisible; if (tempHTvisible>=vMaxHTvisible) tempHTvisible = vMaxHTvisible - 
0.1;
  _HTvisible_Met = tempHTvisible + (vMaxHTvisible - vMinHTvisible)*binMet;

  // _MetMeff_Met
  float tempMetMeff = _metmeff; if (tempMetMeff<vMinMetMeff) tempMetMeff = vMinMetMeff; if (tempMetMeff>=vMaxMetMeff) tempMetMeff = vMaxMetMeff - 0.01;
  _MetMeff_Met = tempMetMeff + (vMaxMetMeff - vMinMetMeff)*binMet;

  // _R2_Met
  float tempR2 = _metmeff; if (tempR2<vMinR2) tempR2 = vMinR2; if (tempR2>=vMaxR2) tempR2 = vMaxR2 - 0.01;
  _R2_Met = tempR2 + (vMaxR2 - vMinR2)*binMet;

  
  // Fake
  _nLeptonsMatched = 0;
  if (fabs(_lep1mid)==24 || fabs(_lep1mid)==23 || fabs(_lep1mid)==15) _nLeptonsMatched++; 
  if (fabs(_lep2mid)==24 || fabs(_lep2mid)==23 || fabs(_lep2mid)==15) _nLeptonsMatched++;
  if (_lep1mid==_lep2mid && (fabs(_lep1mid)==21 || fabs(_lep1mid)<=6)) _nLeptonsMatched = 2;
  //if (fabs(_lep1isfake)<0.1) _nLeptonsMatched++; 
  //if (fabs(_lep2isfake)<0.1) _nLeptonsMatched++;
  
  // ISR jet
  _dphiisrmet = (jetpt1>150. && jetpt1!=_leadingPtCSVv2M) ? acos(cos(metPfType1Phi-jetphi1)) : -1.;
  _hasisrjet = (jetpt1>150. && jetpt1!=_leadingPtCSVv2M && acos(cos(metPfType1Phi-jetphi1))>2.5) ? true : false;

  // 
  TLorentzVector Lep1; Lep1.SetPtEtaPhiM(_lep1pt, _lep1eta, _lep1phi, 0.1); 
  TLorentzVector Lep2; Lep2.SetPtEtaPhiM(_lep2pt, _lep2eta, _lep2phi, 0.1); 
  float dphil1met = fabs(Lep1.DeltaPhi(MET));
  float dphil2met = fabs(Lep2.DeltaPhi(MET));
  _dphiminlmet  = min(dphil1met, dphil2met);

}

//------------------------------------------------------------------------------
// FillAnalysisHistograms
//------------------------------------------------------------------------------
void AnalysisStop::FillAnalysisHistograms(int ichannel,
					  int icut,
					  int ijet)
{
  if (ichannel != ll) FillAnalysisHistograms(ll, icut, ijet);
 
  if ((_SaveHistograms>=2 || _systematic!="nominal") && ijet<njetbin) return;

  h_MT2ll            [ichannel][icut][ijet]->Fill(_MT2ll,          _event_weight);
  h_MT2llgen         [ichannel][icut][ijet]->Fill(_MT2llgen,       _event_weight);

  if (_hasisrjet) {
    h_MT2llisr         [ichannel][icut][ijet]->Fill(_MT2ll,          _event_weight);
    h_MT2llisrgen      [ichannel][icut][ijet]->Fill(_MT2llgen,       _event_weight);
  }

  if (_SaveHistograms>=2 && _systematic=="nominal") {

    if (nvtx>=20) {
      h_MT2ll_nvtxup   [ichannel][icut][ijet]->Fill(_MT2ll,       _event_weight);
      h_MT2llgen_nvtxup[ichannel][icut][ijet]->Fill(_MT2llgen,    _event_weight);
    } else {
      h_MT2ll_nvtxdo   [ichannel][icut][ijet]->Fill(_MT2ll,       _event_weight);
      h_MT2llgen_nvtxdo[ichannel][icut][ijet]->Fill(_MT2llgen,    _event_weight);
    }

    if (_hasisrjet) {
      if (nvtx>=20) {
	h_MT2llisr_nvtxup   [ichannel][icut][ijet]->Fill(_MT2ll,       _event_weight);
	h_MT2llisrgen_nvtxup[ichannel][icut][ijet]->Fill(_MT2llgen,    _event_weight);
      } else {
	h_MT2llisr_nvtxdo   [ichannel][icut][ijet]->Fill(_MT2ll,       _event_weight);
	h_MT2llisrgen_nvtxdo[ichannel][icut][ijet]->Fill(_MT2llgen,    _event_weight);
      }
    }
    
    return;

  }

  if (_SaveHistograms>=2) return; 

  if (_systematic.Contains("Zpeak")) {
    h_HT             [ichannel][icut][ijet]->Fill(_htjets,           _event_weight);
    if (_htjets<200.) h_MT2ll_HTm200 [ichannel][icut][ijet]->Fill(_MT2ll,           _event_weight);
    else              h_MT2ll_HTp200 [ichannel][icut][ijet]->Fill(_MT2ll,           _event_weight);
    if (_htjets<150.) h_MT2ll_HTm150 [ichannel][icut][ijet]->Fill(_MT2ll,           _event_weight);
    else              h_MT2ll_HTp150 [ichannel][icut][ijet]->Fill(_MT2ll,           _event_weight);
    if (_dphiminlmet<2.)
      h_njet20dphilmet [ichannel][icut][ijet]->Fill(_njet,           _event_weight);
  } else if (_systematic.Contains("ZZ") || _systematic.Contains("WZ")) {
    h_M1ll            [ichannel][icut][ijet]->Fill(_m2l,              _event_weight);
    if (_systematic.Contains("ZZ")) h_M2ll            [ichannel][icut][ijet]->Fill(_m2lZ2,            _event_weight);
    else if (_systematic.Contains("WZ")) h_M2ll            [ichannel][icut][ijet]->Fill(_massZcandidate,            _event_weight);
  }

  h_dphillMET        [ichannel][icut][ijet]->Fill(_dphillmet,      _event_weight);
  h_dphiLL           [ichannel][icut][ijet]->Fill(dphill,          _event_weight);
  h_ptLL             [ichannel][icut][ijet]->Fill(_pt2l,           _event_weight);
  h_ptLLbins         [ichannel][icut][ijet]->Fill(_pt2l,           _event_weight);
  h_ptdphiLL         [ichannel][icut][ijet]->Fill(_pt2l, dphill,   _event_weight);
  h_genptLL          [ichannel][icut][ijet]->Fill(_gen_ptll,       _event_weight);
  h_MET              [ichannel][icut][ijet]->Fill(MET.Et(),        _event_weight);
  h_Counter          [ichannel][icut][ijet]->Fill(90.,             _event_weight);
  h_njet20           [ichannel][icut][ijet]->Fill(_njet,           _event_weight);
  h_nbjet            [ichannel][icut][ijet]->Fill(_nbjet30csvv2m,  _event_weight);
  h_njet30           [ichannel][icut][ijet]->Fill(_njet30,         _event_weight);
  h_Lep1Pt           [ichannel][icut][ijet]->Fill(_lep1pt,         _event_weight);
  h_Lep2Pt           [ichannel][icut][ijet]->Fill(_lep2pt,         _event_weight);
  h_Lep1Phi          [ichannel][icut][ijet]->Fill(_lep1phi,        _event_weight);
  h_Lep2Phi          [ichannel][icut][ijet]->Fill(_lep2phi,        _event_weight);
  h_Lep1Eta          [ichannel][icut][ijet]->Fill(_lep1eta,        _event_weight);
  h_Lep2Eta          [ichannel][icut][ijet]->Fill(_lep2eta,        _event_weight);
  h_dphil1MET        [ichannel][icut][ijet]->Fill(dphilmet1,       _event_weight);
  h_dphil2MET        [ichannel][icut][ijet]->Fill(dphilmet2,       _event_weight);
  h_METphi           [ichannel][icut][ijet]->Fill(MET.Phi(),       _event_weight);

  if (_systematic.Contains("DYcorrections")) {
    float pi = acos(-1.);
    if (_dphillmet<pi/8.)         h_dphiLLbin1 [ichannel][icut][ijet]->Fill(dphill, _event_weight);
    else if (_dphillmet<3.*pi/8.) h_dphiLLbin2 [ichannel][icut][ijet]->Fill(dphill, _event_weight);
    else if (_dphillmet<5.*pi/8.) h_dphiLLbin3 [ichannel][icut][ijet]->Fill(dphill, _event_weight);
    else if (_dphillmet<7.*pi/8.) h_dphiLLbin4 [ichannel][icut][ijet]->Fill(dphill, _event_weight);
    else                          h_dphiLLbin5 [ichannel][icut][ijet]->Fill(dphill, _event_weight);
    return;
  }
  
  if (_systematic.Contains("ZZ")) return;
  if (_systematic.Contains("WZ")) return;
  if (_systematic.Contains("ZWtoZ")) return;

  h_dphiminlepmet    [ichannel][icut][ijet]->Fill(_dphiminlmet,    _event_weight);
  
  if (_systematic.Contains("Zpeak")) return;

  if (_nLeptonsMatched==2) {
    h_MT2ll_truth    [ichannel][icut][ijet]->Fill(_MT2ll,      _event_weight);
    h_MET_truth      [ichannel][icut][ijet]->Fill(MET.Et(),        _event_weight);
  } else { 
    h_MT2ll_fake     [ichannel][icut][ijet]->Fill(_MT2ll,      _event_weight);
    h_MET_fake       [ichannel][icut][ijet]->Fill(MET.Et(),        _event_weight);
  }

  if (!_systematic.Contains("kinematic")) return;
  
  h_Jet1Pt           [ichannel][icut][ijet]->Fill(jetpt1,         _event_weight);
  h_Jet2Pt           [ichannel][icut][ijet]->Fill(jetpt2,         _event_weight);
  for (int ij = 0; ij<_jet_pt.size(); ij++) h_JetPt[ichannel][icut][ijet]->Fill(_jet_pt.at(ij), _event_weight);

  h_MT2ll_MET        [ichannel][icut][ijet]->Fill(MET.Et(),  _MT2ll,      _event_weight);

  h_dphiisrmet       [ichannel][icut][ijet]->Fill(_dphiisrmet,     _event_weight);

  h_mt2lblbcomb      [ichannel][icut][ijet]->Fill(_mt2lblbcomb,    _event_weight);
  h_mt2bbtrue        [ichannel][icut][ijet]->Fill(_mt2bbtrue,      _event_weight);
  h_mt2lblbtrue      [ichannel][icut][ijet]->Fill(_mt2lblbtrue,    _event_weight);
  h_mt2lblbmatch     [ichannel][icut][ijet]->Fill(_mt2lblbmatch,   _event_weight);
  h_mlb1comb         [ichannel][icut][ijet]->Fill(_mlb1comb,       _event_weight);
  h_mlb2comb         [ichannel][icut][ijet]->Fill(_mlb2comb,       _event_weight);
  h_mlb1true         [ichannel][icut][ijet]->Fill(_mlb1true,       _event_weight);
  h_mlb2true         [ichannel][icut][ijet]->Fill(_mlb2true,       _event_weight);
  h_mt2lblbvsmlbtrue [ichannel][icut][ijet]->Fill(_mlb1true, _mt2lblbtrue,       _event_weight);
  h_mt2lblbvsmlbtrue [ichannel][icut][ijet]->Fill(_mlb2true, _mt2lblbtrue,       _event_weight);
  h_nisrjet          [ichannel][icut][ijet]->Fill(_nisrjet,        _event_weight);
  h_metmeff          [ichannel][icut][ijet]->Fill(_metmeff,        _event_weight);

  h_MT2_Met          [ichannel][icut][ijet]->Fill(_MT2_Met,        _event_weight);
  h_HTvisible_Met    [ichannel][icut][ijet]->Fill(_HTvisible_Met,  _event_weight);
  h_MetMeff_Met      [ichannel][icut][ijet]->Fill(_MetMeff_Met,    _event_weight);
  h_R2_Met           [ichannel][icut][ijet]->Fill(_R2_Met,         _event_weight);
  

}


//------------------------------------------------------------------------------
// FillSystematicHistograms
//------------------------------------------------------------------------------
void AnalysisStop::FillSystematicHistograms(int ichannel,
					    int icut)
{
  if (ichannel != ll) FillSystematicHistograms(ll, icut);

  for (int is = 0; is<nsystematic; is++) {
  
    if (ssystematic[is]=="nominal" || !systematicfromweight[is]) continue;
    if (ssystematic[is].Contains("Isrnjet") && !_isfastsim) continue;
    
    root_output_systematic[is]->cd();

    float _event_weight_systematic = -999.;
    if (ssystematic[is]=="Triggerup")   _event_weight_systematic = _event_weight_Triggerup;
    if (ssystematic[is]=="Triggerdo")   _event_weight_systematic = _event_weight_Triggerdo;
    if (ssystematic[is]=="Recoup")      _event_weight_systematic = _event_weight_Recoup;
    if (ssystematic[is]=="Recodo")      _event_weight_systematic = _event_weight_Recodo;
    if (ssystematic[is]=="Idisoup")     _event_weight_systematic = _event_weight_Idisoup;
    if (ssystematic[is]=="Idisodo")     _event_weight_systematic = _event_weight_Idisodo;
    if (ssystematic[is]=="Fastsimup")   _event_weight_systematic = _event_weight_Fastsimup;
    if (ssystematic[is]=="Fastsimdo")   _event_weight_systematic = _event_weight_Fastsimdo;
    if (ssystematic[is]=="Btagup")      _event_weight_systematic = _event_weight_Btagup;
    if (ssystematic[is]=="Btagdo")      _event_weight_systematic = _event_weight_Btagdo;
    if (ssystematic[is]=="Btaglightup") _event_weight_systematic = _event_weight_Btaglightup;
    if (ssystematic[is]=="Btaglightdo") _event_weight_systematic = _event_weight_Btaglightdo;
    if (ssystematic[is]=="BtagFSup")    _event_weight_systematic = _event_weight_BtagFSup;
    if (ssystematic[is]=="BtagFSdo")    _event_weight_systematic = _event_weight_BtagFSdo;
    if (ssystematic[is]=="Topptup")     _event_weight_systematic = _event_weight_Toppt;
    if (ssystematic[is]=="Topptdo")     _event_weight_systematic = _event_weight*_event_weight/_event_weight_Toppt;
    if (ssystematic[is]=="Isrnjetup")   _event_weight_systematic = _event_weight_Isrnjetup;
    if (ssystematic[is]=="Isrnjetdo")   _event_weight_systematic = _event_weight_Isrnjetdo;
    if (ssystematic[is]=="Fakeup")      _event_weight_systematic = _event_weight_Fakeup;
    if (ssystematic[is]=="Fakedo")      _event_weight_systematic = _event_weight_Fakedo;

    if (_event_weight_systematic==-999.) 
      printf("\n\n Bad name for systematics, please check!\n");

    h_MT2ll_systematic      [ichannel][icut][is]->Fill(_MT2ll,    _event_weight_systematic);
    h_MT2llgen_systematic   [ichannel][icut][is]->Fill(_MT2llgen, _event_weight_systematic);
    if (_hasisrjet) {
      h_MT2llisr_systematic      [ichannel][icut][is]->Fill(_MT2ll,    _event_weight_systematic);
      h_MT2llisrgen_systematic   [ichannel][icut][is]->Fill(_MT2llgen, _event_weight_systematic);
    }

  }

}


//------------------------------------------------------------------------------
// FillTheoreticalVariationsHistograms
//------------------------------------------------------------------------------
void AnalysisStop::FillTheoreticalVariationsHistograms(int ichannel,
						       int icut)
{
  if (!scut[icut].Contains("SR") && !scut[icut].Contains("VR")) return;
  if (scut[icut].Contains("/ZZ_")) return;

  if (ichannel != ll) FillTheoreticalVariationsHistograms(ll, icut);
    
  root_output_theoreticalvariations->cd();

  bool PrintedError = false;

  int _lheoffset = _isfastsim ? 1 : 0;
  
  for (int tv = 0; tv<nTheoreticalVariations; tv++) {

    int tvo = tv + _lheoffset;

    float _TheoreticalVariationRatio = std_vector_LHE_weight->at(tvo)/std_vector_LHE_weight->at(_lheoffset);
    if (std_vector_LHE_weight->at(tvo)==-999 || std_vector_LHE_weight->at(_lheoffset)==-999) {
      _TheoreticalVariationRatio = _TheoreticalVariationRenormalization[tv];
      if (ichannel==ll && !PrintedError) {
	cout << "  Bad LHE weight " << event << " " << std_vector_LHE_weight->at(tvo) << " " << std_vector_LHE_weight->at(_lheoffset) << endl;
	PrintedError = true;
      }
    }
    
    float _event_weight_theoreticalvariation = _event_weight;
    if (!_isfastsim || hasLHE) {
      _event_weight_theoreticalvariation *= _TheoreticalVariationRatio;
      _event_weight_theoreticalvariation /= _TheoreticalVariationRenormalization[tv];
    }
    
    h_MT2ll_theoreticalvariation      [ichannel][icut][tv]->Fill(_MT2ll,    _event_weight_theoreticalvariation);
    h_MT2llgen_theoreticalvariation   [ichannel][icut][tv]->Fill(_MT2llgen, _event_weight_theoreticalvariation);
    if (_hasisrjet) {
      h_MT2llisr_theoreticalvariation      [ichannel][icut][tv]->Fill(_MT2ll,    _event_weight_theoreticalvariation);
      h_MT2llisrgen_theoreticalvariation   [ichannel][icut][tv]->Fill(_MT2llgen, _event_weight_theoreticalvariation);
    }

  }

}

//------------------------------------------------------------------------------
// SaveSystematicHistograms
//------------------------------------------------------------------------------
void AnalysisStop::SaveSystematicHistograms()
{

  printf("\n\n Writing histograms for systematics.\n");

  for (int is = 0; is<nsystematic; is++) {
  
    if (ssystematic[is]=="nominal" || !systematicfromweight[is]) continue;
    if (ssystematic[is].Contains("Isrnjet") && !_isfastsim) continue;

    root_output_systematic[is]->cd();
    root_output_systematic[is]->Write("", TObject::kOverwrite);
    root_output_systematic[is]->Close();

  }

}

//------------------------------------------------------------------------------
// SaveTheoreticalVariationsHistograms
//------------------------------------------------------------------------------
void AnalysisStop::SaveTheoreticalVariationsHistograms()
{

  printf("\n\n Writing histograms for theoretical variations.\n");

  root_output_theoreticalvariations->cd();
  root_output_theoreticalvariations->Write("", TObject::kOverwrite);
  root_output_theoreticalvariations->Close();
  
}

//------------------------------------------------------------------------------
// FillLevelHistograms
//------------------------------------------------------------------------------
void AnalysisStop::FillLevelHistograms(int  icut,
				       bool pass)
{
  if (!pass) return;
  
  if (_SaveHistograms==0) {
    
    FillHistograms(_channel, icut, _jetbin);
    FillHistograms(_channel, icut, njetbin);
  }

  if (_SaveHistograms>=0) {
    
    FillAnalysisHistograms(_channel, icut, _jetbin);
    FillAnalysisHistograms(_channel, icut, njetbin);

  }

  if (_SaveHistograms>=3) FillSystematicHistograms(_channel, icut);
  if (_DoTheoreticalVariations==1) FillTheoreticalVariationsHistograms(_channel, icut);
 
}

float PileUpSystP1M[90] = { 1.66904, 1.2121, 1.14547, 1.17193, 1.12362, 1.14139, 1.20212, 1.71672, 1.5997, 1.5458, 1.47727, 1.38479, 1.32886, 1.29153, 1.25883, 1.20155, 1.14502, 1.10661, 1.07659, 1.04871, 1.02428, 1.00312, 0.982422, 0.958406, 0.927708, 0.891305, 0.853997, 0.819795, 0.788817, 0.758488, 0.726557, 0.692518, 0.656706, 0.619119, 0.579366, 0.537338, 0.493498, 0.448753, 0.404212, 0.360973, 0.319954, 0.281835, 0.247086, 0.216087, 0.18936, 0.168051, 0.154977, 0.156909, 0.18881, 0.277265, 0.445063, 0.657962, 0.82459, 0.905039, 0.926348, 0.920609, 0.904571, 0.884892, 0.863938, 0.842522, 0.820929, 0.799275, 0.777618, 0.755998, 0.73445, 0.713001, 0.691682, 0.670516, 0.64953, 0.628747, 0.608188, 0.587877, 0.567832, 0.548072, 0.528616, 0.509479, 0.490678, 0.472225, 0.454135, 0.436419, 0.419087, 0.402149, 0.385613, 0.369485, 0.353773, 0.33848, 0.32361, 0.309166, 0.29515, 0.281562  };
float PileUpSystP1P[90] = { 0.596435, 0.83675, 0.865677, 0.860908, 0.903698, 0.874877, 0.882876, 0.656778, 0.610035, 0.657271, 0.669408, 0.710769, 0.746008, 0.771425, 0.78863, 0.810422, 0.847939, 0.884053, 0.910585, 0.933347, 0.955776, 0.976453, 0.995153, 1.0139, 1.03586, 1.065, 1.10268, 1.14602, 1.19032, 1.23333, 1.27689, 1.32485, 1.3801, 1.44392, 1.51751, 1.60349, 1.70588, 1.82928, 1.97837, 2.15799, 2.37316, 2.62912, 2.93164, 3.28683, 3.70079, 4.17768, 4.71363, 5.27955, 5.77728, 5.96274, 5.45787, 4.18824, 2.78028, 1.8282, 1.35224, 1.1549, 1.08731, 1.07443, 1.0837, 1.10203, 1.12445, 1.14916, 1.17554, 1.20342, 1.23278, 1.26366, 1.29613, 1.33025, 1.36612, 1.40383, 1.44347, 1.48515, 1.52899, 1.5751, 1.62361, 1.67465, 1.72838, 1.78495, 1.84452, 1.90728, 1.97341, 2.04312, 2.11661, 2.19413, 2.27593, 2.36226, 2.45341, 2.54968, 2.65141, 2.75894 };
float PileUpSystP2M[90] = { 1.03051, 1.36949, 0.941357, 1.14162, 1.11684, 1.06049, 1.14511, 1.33389, 1.03644, 2.14613, 2.36898, 1.71482, 1.42921, 1.2979, 1.24562, 1.19742, 1.13932, 1.09099, 1.0655, 1.06859, 1.08475, 1.09223, 1.0888, 1.07508, 1.05079, 1.02403, 1.00059, 0.976649, 0.946524, 0.909366, 0.868148, 0.826521, 0.787116, 0.751226, 0.718591, 0.687299, 0.654443, 0.617476, 0.575379, 0.52887, 0.47969, 0.429795, 0.380882, 0.334214, 0.290651, 0.250734, 0.214756, 0.182827, 0.154921, 0.130912, 0.110602, 0.0937429, 0.0800511, 0.0692168, 0.0609123, 0.0548037, 0.0505736, 0.0479621, 0.0468207, 0.047149, 0.0490551, 0.0525922, 0.0574874, 0.0629216, 0.0676229, 0.070378, 0.0706089, 0.0685185, 0.064773, 0.0600994, 0.0550633, 0.0500267, 0.0451911, 0.0406542, 0.0364547, 0.0325986, 0.0290777, 0.025875, 0.022973, 0.02035, 0.0179875, 0.0158573, 0.013966, 0.0122029, 0.0106342, 0.00929965, 0.0076075, 0.0103809, 0.000314326, 0  };
float PileUpSystP2P[90] = { 0.978101, 0.725415, 1.0339, 0.870258, 0.910537, 0.927313, 0.932128, 0.746435, 0.913067, 0.739255, 0.399735, 0.508298, 0.649039, 0.740407, 0.79206, 0.817669, 0.849986, 0.889286, 0.923568, 0.942041, 0.938948, 0.926031, 0.919579, 0.921326, 0.930358, 0.948269, 0.970873, 0.992159, 1.01294, 1.03875, 1.07328, 1.1165, 1.16616, 1.21945, 1.2739, 1.32794, 1.38196, 1.43928, 1.50582, 1.58865, 1.69445, 1.82878, 1.9966, 2.20308, 2.45433, 2.75803, 3.12381, 3.56351, 4.09132, 4.7237, 5.47881, 6.3751, 7.42842, 8.64707, 10.0245, 11.5303, 13.1033, 14.651, 16.0609, 17.2203, 18.0317, 18.4161, 18.3127, 17.7033, 16.6642, 15.3929, 14.1529, 13.171, 12.573, 12.3873, 12.5843, 13.1141, 13.929, 14.9936, 16.2867, 17.8012, 19.5414, 21.5223, 23.7668, 26.3066, 29.1796, 32.4324, 36.1193, 40.3045, 45.0674, 50.4456, 56.715, 63.5715, 72.0067, 84.6335  };

float AnalysisStop::PileupSyst(float runperiod, float putrue, int variation) {

  if (putrue<0.) putrue = 0.;
  if (putrue>=90.) putrue = 89.;
  int ipu = putrue;
  int period = runperiod;
  if (period<=4 && variation==+1) return PileUpSystP1P[ipu];
  else if (period<=4 && variation==-1) return PileUpSystP1M[ipu];
  else if (period>=5 && variation==+1) return PileUpSystP2P[ipu];
  else if (period>=5 && variation==-1) return PileUpSystP2M[ipu];
  return 1.;

}

void AnalysisStop::SetSUSYProductionMap() {

  if (SUSYProductionProcess=="T2tt") {

  MassPoint mp0 (150, 1);
  StopCrossSection cs0 (249.409, 36.7821);
  MassPointParameters mpp0 (cs0, 808201);
  StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

  MassPoint mp1 (150, 25);
  StopCrossSection cs1 (249.409, 36.7821);
  MassPointParameters mpp1 (cs1, 811719);
  StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

  MassPoint mp2 (150, 50);
  StopCrossSection cs2 (249.409, 36.7821);
  MassPointParameters mpp2 (cs2, 797743);
  StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

  MassPoint mp3 (150, 63);
  StopCrossSection cs3 (249.409, 36.7821);
  MassPointParameters mpp3 (cs3, 815018);
  StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

  MassPoint mp4 (167, 1);
  StopCrossSection cs4 (151.469, 22.263);
  MassPointParameters mpp4 (cs4, 783879);
  StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

  MassPoint mp5 (175, 1);
  StopCrossSection cs5 (121.416, 17.7681);
  MassPointParameters mpp5 (cs5, 757245);
  StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

  MassPoint mp6 (175, 25);
  StopCrossSection cs6 (121.416, 17.7681);
  MassPointParameters mpp6 (cs6, 733926);
  StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

  MassPoint mp7 (175, 50);
  StopCrossSection cs7 (121.416, 17.7681);
  MassPointParameters mpp7 (cs7, 745981);
  StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

  MassPoint mp8 (175, 75);
  StopCrossSection cs8 (121.416, 17.7681);
  MassPointParameters mpp8 (cs8, 763353);
  StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

  MassPoint mp9 (175, 88);
  StopCrossSection cs9 (121.416, 17.7681);
  MassPointParameters mpp9 (cs9, 759045);
  StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

  MassPoint mp10 (183, 1);
  StopCrossSection cs10 (98.4784, 14.1473);
  MassPointParameters mpp10 (cs10, 730357);
  StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

  MassPoint mp11 (192, 25);
  StopCrossSection cs11 (78.4483, 11.3431);
  MassPointParameters mpp11 (cs11, 711773);
  StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

  MassPoint mp12 (200, 1);
  StopCrossSection cs12 (64.5085, 9.29555);
  MassPointParameters mpp12 (cs12, 960691);
  StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

  MassPoint mp13 (200, 25);
  StopCrossSection cs13 (64.5085, 9.29555);
  MassPointParameters mpp13 (cs13, 958340);
  StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

  MassPoint mp14 (200, 50);
  StopCrossSection cs14 (64.5085, 9.29555);
  MassPointParameters mpp14 (cs14, 959423);
  StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

  MassPoint mp15 (200, 75);
  StopCrossSection cs15 (64.5085, 9.29555);
  MassPointParameters mpp15 (cs15, 958283);
  StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

  MassPoint mp16 (200, 100);
  StopCrossSection cs16 (64.5085, 9.29555);
  MassPointParameters mpp16 (cs16, 985915);
  StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

  MassPoint mp17 (200, 113);
  StopCrossSection cs17 (64.5085, 9.29555);
  MassPointParameters mpp17 (cs17, 988029);
  StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

  MassPoint mp18 (208, 25);
  StopCrossSection cs18 (53.4447, 7.65327);
  MassPointParameters mpp18 (cs18, 966052);
  StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

  MassPoint mp19 (217, 50);
  StopCrossSection cs19 (43.4633, 6.22129);
  MassPointParameters mpp19 (cs19, 944392);
  StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

  MassPoint mp20 (225, 1);
  StopCrossSection cs20 (36.3818, 5.17309);
  MassPointParameters mpp20 (cs20, 912791);
  StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

  MassPoint mp21 (225, 25);
  StopCrossSection cs21 (36.3818, 5.17309);
  MassPointParameters mpp21 (cs21, 935799);
  StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

  MassPoint mp22 (225, 50);
  StopCrossSection cs22 (36.3818, 5.17309);
  MassPointParameters mpp22 (cs22, 922708);
  StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

  MassPoint mp23 (225, 75);
  StopCrossSection cs23 (36.3818, 5.17309);
  MassPointParameters mpp23 (cs23, 956935);
  StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

  MassPoint mp24 (225, 100);
  StopCrossSection cs24 (36.3818, 5.17309);
  MassPointParameters mpp24 (cs24, 957650);
  StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

  MassPoint mp25 (225, 125);
  StopCrossSection cs25 (36.3818, 5.17309);
  MassPointParameters mpp25 (cs25, 932550);
  StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

  MassPoint mp26 (225, 138);
  StopCrossSection cs26 (36.3818, 5.17309);
  MassPointParameters mpp26 (cs26, 924977);
  StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

  MassPoint mp27 (233, 50);
  StopCrossSection cs27 (30.6565, 4.35198);
  MassPointParameters mpp27 (cs27, 909789);
  StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

  MassPoint mp28 (242, 75);
  StopCrossSection cs28 (25.4398, 3.58399);
  MassPointParameters mpp28 (cs28, 867073);
  StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

  MassPoint mp29 (250, 1);
  StopCrossSection cs29 (21.5949, 3.03613);
  MassPointParameters mpp29 (cs29, 878608);
  StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

  MassPoint mp30 (250, 25);
  StopCrossSection cs30 (21.5949, 3.03613);
  MassPointParameters mpp30 (cs30, 889926);
  StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

  MassPoint mp31 (250, 50);
  StopCrossSection cs31 (21.5949, 3.03613);
  MassPointParameters mpp31 (cs31, 874381);
  StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

  MassPoint mp32 (250, 75);
  StopCrossSection cs32 (21.5949, 3.03613);
  MassPointParameters mpp32 (cs32, 859765);
  StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

  MassPoint mp33 (250, 100);
  StopCrossSection cs33 (21.5949, 3.03613);
  MassPointParameters mpp33 (cs33, 872729);
  StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

  MassPoint mp34 (250, 125);
  StopCrossSection cs34 (21.5949, 3.03613);
  MassPointParameters mpp34 (cs34, 877004);
  StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

  MassPoint mp35 (250, 150);
  StopCrossSection cs35 (21.5949, 3.03613);
  MassPointParameters mpp35 (cs35, 893386);
  StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

  MassPoint mp36 (250, 163);
  StopCrossSection cs36 (21.5949, 3.03613);
  MassPointParameters mpp36 (cs36, 878259);
  StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

  MassPoint mp37 (258, 75);
  StopCrossSection cs37 (18.4347, 2.56587);
  MassPointParameters mpp37 (cs37, 853986);
  StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

  MassPoint mp38 (267, 100);
  StopCrossSection cs38 (15.5256, 2.16481);
  MassPointParameters mpp38 (cs38, 845095);
  StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

  MassPoint mp39 (275, 1);
  StopCrossSection cs39 (13.3231, 1.89919);
  MassPointParameters mpp39 (cs39, 861956);
  StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

  MassPoint mp40 (275, 25);
  StopCrossSection cs40 (13.3231, 1.89919);
  MassPointParameters mpp40 (cs40, 824537);
  StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

  MassPoint mp41 (275, 50);
  StopCrossSection cs41 (13.3231, 1.89919);
  MassPointParameters mpp41 (cs41, 851997);
  StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

  MassPoint mp42 (275, 75);
  StopCrossSection cs42 (13.3231, 1.89919);
  MassPointParameters mpp42 (cs42, 844693);
  StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

  MassPoint mp43 (275, 100);
  StopCrossSection cs43 (13.3231, 1.89919);
  MassPointParameters mpp43 (cs43, 838466);
  StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

  MassPoint mp44 (275, 125);
  StopCrossSection cs44 (13.3231, 1.89919);
  MassPointParameters mpp44 (cs44, 829187);
  StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

  MassPoint mp45 (275, 150);
  StopCrossSection cs45 (13.3231, 1.89919);
  MassPointParameters mpp45 (cs45, 840658);
  StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

  MassPoint mp46 (275, 175);
  StopCrossSection cs46 (13.3231, 1.89919);
  MassPointParameters mpp46 (cs46, 846634);
  StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

  MassPoint mp47 (275, 188);
  StopCrossSection cs47 (13.3231, 1.89919);
  MassPointParameters mpp47 (cs47, 845385);
  StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

  MassPoint mp48 (283, 100);
  StopCrossSection cs48 (11.5185, 1.62631);
  MassPointParameters mpp48 (cs48, 851914);
  StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

  MassPoint mp49 (292, 125);
  StopCrossSection cs49 (9.79779, 1.36209);
  MassPointParameters mpp49 (cs49, 783430);
  StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

  MassPoint mp50 (300, 1);
  StopCrossSection cs50 (8.51615, 1.18564);
  MassPointParameters mpp50 (cs50, 975907);
  StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

  MassPoint mp51 (300, 25);
  StopCrossSection cs51 (8.51615, 1.18564);
  MassPointParameters mpp51 (cs51, 934485);
  StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

  MassPoint mp52 (300, 50);
  StopCrossSection cs52 (8.51615, 1.18564);
  MassPointParameters mpp52 (cs52, 983303);
  StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

  MassPoint mp53 (300, 75);
  StopCrossSection cs53 (8.51615, 1.18564);
  MassPointParameters mpp53 (cs53, 978530);
  StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

  MassPoint mp54 (300, 100);
  StopCrossSection cs54 (8.51615, 1.18564);
  MassPointParameters mpp54 (cs54, 975157);
  StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

  MassPoint mp55 (300, 125);
  StopCrossSection cs55 (8.51615, 1.18564);
  MassPointParameters mpp55 (cs55, 966457);
  StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

  MassPoint mp56 (300, 150);
  StopCrossSection cs56 (8.51615, 1.18564);
  MassPointParameters mpp56 (cs56, 970794);
  StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

  MassPoint mp57 (300, 175);
  StopCrossSection cs57 (8.51615, 1.18564);
  MassPointParameters mpp57 (cs57, 980640);
  StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

  MassPoint mp58 (300, 200);
  StopCrossSection cs58 (8.51615, 1.18564);
  MassPointParameters mpp58 (cs58, 962116);
  StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

  MassPoint mp59 (300, 213);
  StopCrossSection cs59 (8.51615, 1.18564);
  MassPointParameters mpp59 (cs59, 960632);
  StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

  MassPoint mp60 (308, 125);
  StopCrossSection cs60 (7.43297, 1.03471);
  MassPointParameters mpp60 (cs60, 954600);
  StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

  MassPoint mp61 (317, 150);
  StopCrossSection cs61 (6.39537, 0.887432);
  MassPointParameters mpp61 (cs61, 962749);
  StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

  MassPoint mp62 (325, 25);
  StopCrossSection cs62 (5.60471, 0.774257);
  MassPointParameters mpp62 (cs62, 952721);
  StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

  MassPoint mp63 (325, 50);
  StopCrossSection cs63 (5.60471, 0.774257);
  MassPointParameters mpp63 (cs63, 962633);
  StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

  MassPoint mp64 (325, 75);
  StopCrossSection cs64 (5.60471, 0.774257);
  MassPointParameters mpp64 (cs64, 948917);
  StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

  MassPoint mp65 (325, 100);
  StopCrossSection cs65 (5.60471, 0.774257);
  MassPointParameters mpp65 (cs65, 960304);
  StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

  MassPoint mp66 (325, 125);
  StopCrossSection cs66 (5.60471, 0.774257);
  MassPointParameters mpp66 (cs66, 930616);
  StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

  MassPoint mp67 (325, 150);
  StopCrossSection cs67 (5.60471, 0.774257);
  MassPointParameters mpp67 (cs67, 956900);
  StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

  MassPoint mp68 (325, 175);
  StopCrossSection cs68 (5.60471, 0.774257);
  MassPointParameters mpp68 (cs68, 950666);
  StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

  MassPoint mp69 (325, 200);
  StopCrossSection cs69 (5.60471, 0.774257);
  MassPointParameters mpp69 (cs69, 934111);
  StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

  MassPoint mp70 (325, 225);
  StopCrossSection cs70 (5.60471, 0.774257);
  MassPointParameters mpp70 (cs70, 958895);
  StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

  MassPoint mp71 (325, 238);
  StopCrossSection cs71 (5.60471, 0.774257);
  MassPointParameters mpp71 (cs71, 957518);
  StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

  MassPoint mp72 (333, 150);
  StopCrossSection cs72 (4.93598, 0.677722);
  MassPointParameters mpp72 (cs72, 908690);
  StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

  MassPoint mp73 (342, 175);
  StopCrossSection cs73 (4.2853, 0.589713);
  MassPointParameters mpp73 (cs73, 919345);
  StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

  MassPoint mp74 (350, 1);
  StopCrossSection cs74 (3.78661, 0.5183);
  MassPointParameters mpp74 (cs74, 907345);
  StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

  MassPoint mp75 (350, 50);
  StopCrossSection cs75 (3.78661, 0.5183);
  MassPointParameters mpp75 (cs75, 898779);
  StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

  MassPoint mp76 (350, 75);
  StopCrossSection cs76 (3.78661, 0.5183);
  MassPointParameters mpp76 (cs76, 910911);
  StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

  MassPoint mp77 (350, 100);
  StopCrossSection cs77 (3.78661, 0.5183);
  MassPointParameters mpp77 (cs77, 900543);
  StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

  MassPoint mp78 (350, 125);
  StopCrossSection cs78 (3.78661, 0.5183);
  MassPointParameters mpp78 (cs78, 901150);
  StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

  MassPoint mp79 (350, 150);
  StopCrossSection cs79 (3.78661, 0.5183);
  MassPointParameters mpp79 (cs79, 905247);
  StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

  MassPoint mp80 (350, 175);
  StopCrossSection cs80 (3.78661, 0.5183);
  MassPointParameters mpp80 (cs80, 883578);
  StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

  MassPoint mp81 (350, 200);
  StopCrossSection cs81 (3.78661, 0.5183);
  MassPointParameters mpp81 (cs81, 894275);
  StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

  MassPoint mp82 (350, 225);
  StopCrossSection cs82 (3.78661, 0.5183);
  MassPointParameters mpp82 (cs82, 894541);
  StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

  MassPoint mp83 (350, 250);
  StopCrossSection cs83 (3.78661, 0.5183);
  MassPointParameters mpp83 (cs83, 896650);
  StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

  MassPoint mp84 (350, 263);
  StopCrossSection cs84 (3.78661, 0.5183);
  MassPointParameters mpp84 (cs84, 907384);
  StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

  MassPoint mp85 (358, 175);
  StopCrossSection cs85 (3.35736, 0.463444);
  MassPointParameters mpp85 (cs85, 882289);
  StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

  MassPoint mp86 (367, 200);
  StopCrossSection cs86 (2.93791, 0.403858);
  MassPointParameters mpp86 (cs86, 885117);
  StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

  MassPoint mp87 (375, 75);
  StopCrossSection cs87 (2.61162, 0.361649);
  MassPointParameters mpp87 (cs87, 856669);
  StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

  MassPoint mp88 (375, 100);
  StopCrossSection cs88 (2.61162, 0.361649);
  MassPointParameters mpp88 (cs88, 882061);
  StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

  MassPoint mp89 (375, 125);
  StopCrossSection cs89 (2.61162, 0.361649);
  MassPointParameters mpp89 (cs89, 890334);
  StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

  MassPoint mp90 (375, 150);
  StopCrossSection cs90 (2.61162, 0.361649);
  MassPointParameters mpp90 (cs90, 891883);
  StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

  MassPoint mp91 (375, 175);
  StopCrossSection cs91 (2.61162, 0.361649);
  MassPointParameters mpp91 (cs91, 865630);
  StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

  MassPoint mp92 (375, 200);
  StopCrossSection cs92 (2.61162, 0.361649);
  MassPointParameters mpp92 (cs92, 913123);
  StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

  MassPoint mp93 (375, 225);
  StopCrossSection cs93 (2.61162, 0.361649);
  MassPointParameters mpp93 (cs93, 890175);
  StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

  MassPoint mp94 (375, 250);
  StopCrossSection cs94 (2.61162, 0.361649);
  MassPointParameters mpp94 (cs94, 883829);
  StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

  MassPoint mp95 (375, 275);
  StopCrossSection cs95 (2.61162, 0.361649);
  MassPointParameters mpp95 (cs95, 889172);
  StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

  MassPoint mp96 (375, 288);
  StopCrossSection cs96 (2.61162, 0.361649);
  MassPointParameters mpp96 (cs96, 877750);
  StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

  MassPoint mp97 (383, 200);
  StopCrossSection cs97 (2.33031, 0.319632);
  MassPointParameters mpp97 (cs97, 816009);
  StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

  MassPoint mp98 (392, 225);
  StopCrossSection cs98 (2.05132, 0.279655);
  MassPointParameters mpp98 (cs98, 712720);
  StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

  MassPoint mp99 (400, 1);
  StopCrossSection cs99 (1.83537, 0.251418);
  MassPointParameters mpp99 (cs99, 698496);
  StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

  MassPoint mp100 (400, 50);
  StopCrossSection cs100 (1.83537, 0.251418);
  MassPointParameters mpp100 (cs100, 701993);
  StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

  MassPoint mp101 (400, 100);
  StopCrossSection cs101 (1.83537, 0.251418);
  MassPointParameters mpp101 (cs101, 685608);
  StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

  MassPoint mp102 (400, 125);
  StopCrossSection cs102 (1.83537, 0.251418);
  MassPointParameters mpp102 (cs102, 696052);
  StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

  MassPoint mp103 (400, 150);
  StopCrossSection cs103 (1.83537, 0.251418);
  MassPointParameters mpp103 (cs103, 703616);
  StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

  MassPoint mp104 (400, 175);
  StopCrossSection cs104 (1.83537, 0.251418);
  MassPointParameters mpp104 (cs104, 697935);
  StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

  MassPoint mp105 (400, 200);
  StopCrossSection cs105 (1.83537, 0.251418);
  MassPointParameters mpp105 (cs105, 701844);
  StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

  MassPoint mp106 (400, 225);
  StopCrossSection cs106 (1.83537, 0.251418);
  MassPointParameters mpp106 (cs106, 701488);
  StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

  MassPoint mp107 (400, 250);
  StopCrossSection cs107 (1.83537, 0.251418);
  MassPointParameters mpp107 (cs107, 697102);
  StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

  MassPoint mp108 (400, 275);
  StopCrossSection cs108 (1.83537, 0.251418);
  MassPointParameters mpp108 (cs108, 686788);
  StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

  MassPoint mp109 (400, 300);
  StopCrossSection cs109 (1.83537, 0.251418);
  MassPointParameters mpp109 (cs109, 704446);
  StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

  MassPoint mp110 (400, 313);
  StopCrossSection cs110 (1.83537, 0.251418);
  MassPointParameters mpp110 (cs110, 699713);
  StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

  MassPoint mp111 (408, 225);
  StopCrossSection cs111 (1.64598, 0.224102);
  MassPointParameters mpp111 (cs111, 629945);
  StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

  MassPoint mp112 (417, 250);
  StopCrossSection cs112 (1.45754, 0.197237);
  MassPointParameters mpp112 (cs112, 554977);
  StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

  MassPoint mp113 (425, 125);
  StopCrossSection cs113 (1.31169, 0.177095);
  MassPointParameters mpp113 (cs113, 480381);
  StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

  MassPoint mp114 (425, 150);
  StopCrossSection cs114 (1.31169, 0.177095);
  MassPointParameters mpp114 (cs114, 482995);
  StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

  MassPoint mp115 (425, 175);
  StopCrossSection cs115 (1.31169, 0.177095);
  MassPointParameters mpp115 (cs115, 481395);
  StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

  MassPoint mp116 (425, 200);
  StopCrossSection cs116 (1.31169, 0.177095);
  MassPointParameters mpp116 (cs116, 503771);
  StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

  MassPoint mp117 (425, 225);
  StopCrossSection cs117 (1.31169, 0.177095);
  MassPointParameters mpp117 (cs117, 486428);
  StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

  MassPoint mp118 (425, 250);
  StopCrossSection cs118 (1.31169, 0.177095);
  MassPointParameters mpp118 (cs118, 494006);
  StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

  MassPoint mp119 (425, 275);
  StopCrossSection cs119 (1.31169, 0.177095);
  MassPointParameters mpp119 (cs119, 481194);
  StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

  MassPoint mp120 (425, 300);
  StopCrossSection cs120 (1.31169, 0.177095);
  MassPointParameters mpp120 (cs120, 487595);
  StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

  MassPoint mp121 (425, 325);
  StopCrossSection cs121 (1.31169, 0.177095);
  MassPointParameters mpp121 (cs121, 482190);
  StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

  MassPoint mp122 (425, 338);
  StopCrossSection cs122 (1.31169, 0.177095);
  MassPointParameters mpp122 (cs122, 471920);
  StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

  MassPoint mp123 (433, 250);
  StopCrossSection cs123 (1.17767, 0.15845);
  MassPointParameters mpp123 (cs123, 441590);
  StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

  MassPoint mp124 (442, 275);
  StopCrossSection cs124 (1.04898, 0.142727);
  MassPointParameters mpp124 (cs124, 372555);
  StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

  MassPoint mp125 (450, 1);
  StopCrossSection cs125 (0.948333, 0.127607);
  MassPointParameters mpp125 (cs125, 345199);
  StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

  MassPoint mp126 (450, 50);
  StopCrossSection cs126 (0.948333, 0.127607);
  MassPointParameters mpp126 (cs126, 360986);
  StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

  MassPoint mp127 (450, 100);
  StopCrossSection cs127 (0.948333, 0.127607);
  MassPointParameters mpp127 (cs127, 340979);
  StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

  MassPoint mp128 (450, 150);
  StopCrossSection cs128 (0.948333, 0.127607);
  MassPointParameters mpp128 (cs128, 348945);
  StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

  MassPoint mp129 (450, 175);
  StopCrossSection cs129 (0.948333, 0.127607);
  MassPointParameters mpp129 (cs129, 355531);
  StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

  MassPoint mp130 (450, 200);
  StopCrossSection cs130 (0.948333, 0.127607);
  MassPointParameters mpp130 (cs130, 338247);
  StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

  MassPoint mp131 (450, 225);
  StopCrossSection cs131 (0.948333, 0.127607);
  MassPointParameters mpp131 (cs131, 353388);
  StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

  MassPoint mp132 (450, 250);
  StopCrossSection cs132 (0.948333, 0.127607);
  MassPointParameters mpp132 (cs132, 342538);
  StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

  MassPoint mp133 (450, 275);
  StopCrossSection cs133 (0.948333, 0.127607);
  MassPointParameters mpp133 (cs133, 341142);
  StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

  MassPoint mp134 (450, 300);
  StopCrossSection cs134 (0.948333, 0.127607);
  MassPointParameters mpp134 (cs134, 347312);
  StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

  MassPoint mp135 (450, 325);
  StopCrossSection cs135 (0.948333, 0.127607);
  MassPointParameters mpp135 (cs135, 331272);
  StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

  MassPoint mp136 (450, 350);
  StopCrossSection cs136 (0.948333, 0.127607);
  MassPointParameters mpp136 (cs136, 345538);
  StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

  MassPoint mp137 (450, 363);
  StopCrossSection cs137 (0.948333, 0.127607);
  MassPointParameters mpp137 (cs137, 341472);
  StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

  MassPoint mp138 (458, 275);
  StopCrossSection cs138 (0.858396, 0.115469);
  MassPointParameters mpp138 (cs138, 320588);
  StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

  MassPoint mp139 (467, 300);
  StopCrossSection cs139 (0.768552, 0.103094);
  MassPointParameters mpp139 (cs139, 277504);
  StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

  MassPoint mp140 (475, 175);
  StopCrossSection cs140 (0.697075, 0.0933565);
  MassPointParameters mpp140 (cs140, 249037);
  StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

  MassPoint mp141 (475, 200);
  StopCrossSection cs141 (0.697075, 0.0933565);
  MassPointParameters mpp141 (cs141, 245396);
  StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

  MassPoint mp142 (475, 225);
  StopCrossSection cs142 (0.697075, 0.0933565);
  MassPointParameters mpp142 (cs142, 252206);
  StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

  MassPoint mp143 (475, 250);
  StopCrossSection cs143 (0.697075, 0.0933565);
  MassPointParameters mpp143 (cs143, 250149);
  StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

  MassPoint mp144 (475, 275);
  StopCrossSection cs144 (0.697075, 0.0933565);
  MassPointParameters mpp144 (cs144, 259063);
  StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

  MassPoint mp145 (475, 300);
  StopCrossSection cs145 (0.697075, 0.0933565);
  MassPointParameters mpp145 (cs145, 238918);
  StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

  MassPoint mp146 (475, 325);
  StopCrossSection cs146 (0.697075, 0.0933565);
  MassPointParameters mpp146 (cs146, 255027);
  StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

  MassPoint mp147 (475, 350);
  StopCrossSection cs147 (0.697075, 0.0933565);
  MassPointParameters mpp147 (cs147, 258151);
  StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

  MassPoint mp148 (475, 375);
  StopCrossSection cs148 (0.697075, 0.0933565);
  MassPointParameters mpp148 (cs148, 246829);
  StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

  MassPoint mp149 (475, 388);
  StopCrossSection cs149 (0.697075, 0.0933565);
  MassPointParameters mpp149 (cs149, 243991);
  StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

  MassPoint mp150 (483, 300);
  StopCrossSection cs150 (0.633519, 0.0848849);
  MassPointParameters mpp150 (cs150, 236700);
  StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

  MassPoint mp151 (492, 325);
  StopCrossSection cs151 (0.56929, 0.0761869);
  MassPointParameters mpp151 (cs151, 200396);
  StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

  MassPoint mp152 (500, 1);
  StopCrossSection cs152 (0.51848, 0.0693711);
  MassPointParameters mpp152 (cs152, 201790);
  StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

  MassPoint mp153 (500, 50);
  StopCrossSection cs153 (0.51848, 0.0693711);
  MassPointParameters mpp153 (cs153, 204320);
  StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

  MassPoint mp154 (500, 100);
  StopCrossSection cs154 (0.51848, 0.0693711);
  MassPointParameters mpp154 (cs154, 193841);
  StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

  MassPoint mp155 (500, 150);
  StopCrossSection cs155 (0.51848, 0.0693711);
  MassPointParameters mpp155 (cs155, 196526);
  StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

  MassPoint mp156 (500, 200);
  StopCrossSection cs156 (0.51848, 0.0693711);
  MassPointParameters mpp156 (cs156, 198566);
  StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

  MassPoint mp157 (500, 225);
  StopCrossSection cs157 (0.51848, 0.0693711);
  MassPointParameters mpp157 (cs157, 202707);
  StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

  MassPoint mp158 (500, 250);
  StopCrossSection cs158 (0.51848, 0.0693711);
  MassPointParameters mpp158 (cs158, 193604);
  StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

  MassPoint mp159 (500, 275);
  StopCrossSection cs159 (0.51848, 0.0693711);
  MassPointParameters mpp159 (cs159, 201321);
  StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

  MassPoint mp160 (500, 300);
  StopCrossSection cs160 (0.51848, 0.0693711);
  MassPointParameters mpp160 (cs160, 191044);
  StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

  MassPoint mp161 (500, 325);
  StopCrossSection cs161 (0.51848, 0.0693711);
  MassPointParameters mpp161 (cs161, 190896);
  StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

  MassPoint mp162 (500, 350);
  StopCrossSection cs162 (0.51848, 0.0693711);
  MassPointParameters mpp162 (cs162, 190200);
  StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

  MassPoint mp163 (500, 375);
  StopCrossSection cs163 (0.51848, 0.0693711);
  MassPointParameters mpp163 (cs163, 188442);
  StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

  MassPoint mp164 (500, 400);
  StopCrossSection cs164 (0.51848, 0.0693711);
  MassPointParameters mpp164 (cs164, 196993);
  StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

  MassPoint mp165 (500, 413);
  StopCrossSection cs165 (0.51848, 0.0693711);
  MassPointParameters mpp165 (cs165, 189000);
  StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

  MassPoint mp166 (508, 325);
  StopCrossSection cs166 (0.473193, 0.0630664);
  MassPointParameters mpp166 (cs166, 184242);
  StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

  MassPoint mp167 (517, 350);
  StopCrossSection cs167 (0.42723, 0.0569597);
  MassPointParameters mpp167 (cs167, 163267);
  StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

  MassPoint mp168 (525, 225);
  StopCrossSection cs168 (0.390303, 0.0520832);
  MassPointParameters mpp168 (cs168, 149626);
  StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

  MassPoint mp169 (525, 250);
  StopCrossSection cs169 (0.390303, 0.0520832);
  MassPointParameters mpp169 (cs169, 146799);
  StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

  MassPoint mp170 (525, 275);
  StopCrossSection cs170 (0.390303, 0.0520832);
  MassPointParameters mpp170 (cs170, 153427);
  StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

  MassPoint mp171 (525, 300);
  StopCrossSection cs171 (0.390303, 0.0520832);
  MassPointParameters mpp171 (cs171, 154624);
  StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

  MassPoint mp172 (525, 325);
  StopCrossSection cs172 (0.390303, 0.0520832);
  MassPointParameters mpp172 (cs172, 138642);
  StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

  MassPoint mp173 (525, 350);
  StopCrossSection cs173 (0.390303, 0.0520832);
  MassPointParameters mpp173 (cs173, 149976);
  StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

  MassPoint mp174 (525, 375);
  StopCrossSection cs174 (0.390303, 0.0520832);
  MassPointParameters mpp174 (cs174, 148360);
  StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

  MassPoint mp175 (525, 400);
  StopCrossSection cs175 (0.390303, 0.0520832);
  MassPointParameters mpp175 (cs175, 142655);
  StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

  MassPoint mp176 (525, 425);
  StopCrossSection cs176 (0.390303, 0.0520832);
  MassPointParameters mpp176 (cs176, 138180);
  StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

  MassPoint mp177 (525, 438);
  StopCrossSection cs177 (0.390303, 0.0520832);
  MassPointParameters mpp177 (cs177, 134106);
  StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

  MassPoint mp178 (533, 350);
  StopCrossSection cs178 (0.356725, 0.0474963);
  MassPointParameters mpp178 (cs178, 148390);
  StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

  MassPoint mp179 (542, 375);
  StopCrossSection cs179 (0.323163, 0.0429481);
  MassPointParameters mpp179 (cs179, 120506);
  StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

  MassPoint mp180 (550, 1);
  StopCrossSection cs180 (0.296128, 0.0392923);
  MassPointParameters mpp180 (cs180, 114315);
  StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

  MassPoint mp181 (550, 50);
  StopCrossSection cs181 (0.296128, 0.0392923);
  MassPointParameters mpp181 (cs181, 109953);
  StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

  MassPoint mp182 (550, 100);
  StopCrossSection cs182 (0.296128, 0.0392923);
  MassPointParameters mpp182 (cs182, 106444);
  StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

  MassPoint mp183 (550, 150);
  StopCrossSection cs183 (0.296128, 0.0392923);
  MassPointParameters mpp183 (cs183, 114667);
  StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

  MassPoint mp184 (550, 200);
  StopCrossSection cs184 (0.296128, 0.0392923);
  MassPointParameters mpp184 (cs184, 107117);
  StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

  MassPoint mp185 (550, 250);
  StopCrossSection cs185 (0.296128, 0.0392923);
  MassPointParameters mpp185 (cs185, 106871);
  StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

  MassPoint mp186 (550, 275);
  StopCrossSection cs186 (0.296128, 0.0392923);
  MassPointParameters mpp186 (cs186, 114239);
  StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

  MassPoint mp187 (550, 300);
  StopCrossSection cs187 (0.296128, 0.0392923);
  MassPointParameters mpp187 (cs187, 105867);
  StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

  MassPoint mp188 (550, 325);
  StopCrossSection cs188 (0.296128, 0.0392923);
  MassPointParameters mpp188 (cs188, 109650);
  StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

  MassPoint mp189 (550, 350);
  StopCrossSection cs189 (0.296128, 0.0392923);
  MassPointParameters mpp189 (cs189, 121720);
  StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

  MassPoint mp190 (550, 375);
  StopCrossSection cs190 (0.296128, 0.0392923);
  MassPointParameters mpp190 (cs190, 111063);
  StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

  MassPoint mp191 (550, 400);
  StopCrossSection cs191 (0.296128, 0.0392923);
  MassPointParameters mpp191 (cs191, 120961);
  StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

  MassPoint mp192 (550, 425);
  StopCrossSection cs192 (0.296128, 0.0392923);
  MassPointParameters mpp192 (cs192, 110705);
  StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

  MassPoint mp193 (550, 450);
  StopCrossSection cs193 (0.296128, 0.0392923);
  MassPointParameters mpp193 (cs193, 106443);
  StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

  MassPoint mp194 (550, 463);
  StopCrossSection cs194 (0.296128, 0.0392923);
  MassPointParameters mpp194 (cs194, 114809);
  StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

  MassPoint mp195 (558, 375);
  StopCrossSection cs195 (0.271976, 0.0359305);
  MassPointParameters mpp195 (cs195, 100950);
  StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

  MassPoint mp196 (567, 400);
  StopCrossSection cs196 (0.246349, 0.0326119);
  MassPointParameters mpp196 (cs196, 89804);
  StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

  MassPoint mp197 (575, 275);
  StopCrossSection cs197 (0.226118, 0.0300151);
  MassPointParameters mpp197 (cs197, 85820);
  StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

  MassPoint mp198 (575, 300);
  StopCrossSection cs198 (0.226118, 0.0300151);
  MassPointParameters mpp198 (cs198, 91177);
  StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

  MassPoint mp199 (575, 325);
  StopCrossSection cs199 (0.226118, 0.0300151);
  MassPointParameters mpp199 (cs199, 84924);
  StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

  MassPoint mp200 (575, 350);
  StopCrossSection cs200 (0.226118, 0.0300151);
  MassPointParameters mpp200 (cs200, 87200);
  StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

  MassPoint mp201 (575, 375);
  StopCrossSection cs201 (0.226118, 0.0300151);
  MassPointParameters mpp201 (cs201, 80105);
  StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

  MassPoint mp202 (575, 400);
  StopCrossSection cs202 (0.226118, 0.0300151);
  MassPointParameters mpp202 (cs202, 90596);
  StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

  MassPoint mp203 (575, 425);
  StopCrossSection cs203 (0.226118, 0.0300151);
  MassPointParameters mpp203 (cs203, 79896);
  StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

  MassPoint mp204 (575, 450);
  StopCrossSection cs204 (0.226118, 0.0300151);
  MassPointParameters mpp204 (cs204, 84527);
  StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

  MassPoint mp205 (575, 475);
  StopCrossSection cs205 (0.226118, 0.0300151);
  MassPointParameters mpp205 (cs205, 80622);
  StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

  MassPoint mp206 (575, 488);
  StopCrossSection cs206 (0.226118, 0.0300151);
  MassPointParameters mpp206 (cs206, 76205);
  StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

  MassPoint mp207 (583, 400);
  StopCrossSection cs207 (0.207962, 0.0275786);
  MassPointParameters mpp207 (cs207, 83316);
  StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

  MassPoint mp208 (592, 425);
  StopCrossSection cs208 (0.189289, 0.024915);
  MassPointParameters mpp208 (cs208, 64819);
  StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

  MassPoint mp209 (600, 1);
  StopCrossSection cs209 (0.174599, 0.02306);
  MassPointParameters mpp209 (cs209, 65444);
  StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

  MassPoint mp210 (600, 50);
  StopCrossSection cs210 (0.174599, 0.02306);
  MassPointParameters mpp210 (cs210, 67689);
  StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

  MassPoint mp211 (600, 100);
  StopCrossSection cs211 (0.174599, 0.02306);
  MassPointParameters mpp211 (cs211, 77545);
  StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

  MassPoint mp212 (600, 150);
  StopCrossSection cs212 (0.174599, 0.02306);
  MassPointParameters mpp212 (cs212, 70701);
  StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

  MassPoint mp213 (600, 200);
  StopCrossSection cs213 (0.174599, 0.02306);
  MassPointParameters mpp213 (cs213, 61618);
  StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

  MassPoint mp214 (600, 250);
  StopCrossSection cs214 (0.174599, 0.02306);
  MassPointParameters mpp214 (cs214, 72780);
  StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

  MassPoint mp215 (600, 300);
  StopCrossSection cs215 (0.174599, 0.02306);
  MassPointParameters mpp215 (cs215, 72083);
  StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

  MassPoint mp216 (600, 325);
  StopCrossSection cs216 (0.174599, 0.02306);
  MassPointParameters mpp216 (cs216, 67829);
  StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

  MassPoint mp217 (600, 350);
  StopCrossSection cs217 (0.174599, 0.02306);
  MassPointParameters mpp217 (cs217, 70760);
  StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

  MassPoint mp218 (600, 375);
  StopCrossSection cs218 (0.174599, 0.02306);
  MassPointParameters mpp218 (cs218, 64281);
  StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

  MassPoint mp219 (600, 400);
  StopCrossSection cs219 (0.174599, 0.02306);
  MassPointParameters mpp219 (cs219, 68997);
  StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

  MassPoint mp220 (600, 425);
  StopCrossSection cs220 (0.174599, 0.02306);
  MassPointParameters mpp220 (cs220, 74713);
  StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

  MassPoint mp221 (600, 450);
  StopCrossSection cs221 (0.174599, 0.02306);
  MassPointParameters mpp221 (cs221, 72238);
  StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

  MassPoint mp222 (600, 475);
  StopCrossSection cs222 (0.174599, 0.02306);
  MassPointParameters mpp222 (cs222, 73545);
  StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

  MassPoint mp223 (600, 500);
  StopCrossSection cs223 (0.174599, 0.02306);
  MassPointParameters mpp223 (cs223, 62982);
  StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

  MassPoint mp224 (600, 513);
  StopCrossSection cs224 (0.174599, 0.02306);
  MassPointParameters mpp224 (cs224, 70385);
  StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

  MassPoint mp225 (608, 425);
  StopCrossSection cs225 (0.161398, 0.0211267);
  MassPointParameters mpp225 (cs225, 57022);
  StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

  MassPoint mp226 (617, 450);
  StopCrossSection cs226 (0.14728, 0.01944);
  MassPointParameters mpp226 (cs226, 61131);
  StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

  MassPoint mp227 (625, 325);
  StopCrossSection cs227 (0.136372, 0.0174504);
  MassPointParameters mpp227 (cs227, 56053);
  StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

  MassPoint mp228 (625, 350);
  StopCrossSection cs228 (0.136372, 0.0174504);
  MassPointParameters mpp228 (cs228, 55670);
  StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

  MassPoint mp229 (625, 375);
  StopCrossSection cs229 (0.136372, 0.0174504);
  MassPointParameters mpp229 (cs229, 51913);
  StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

  MassPoint mp230 (625, 400);
  StopCrossSection cs230 (0.136372, 0.0174504);
  MassPointParameters mpp230 (cs230, 52024);
  StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

  MassPoint mp231 (625, 425);
  StopCrossSection cs231 (0.136372, 0.0174504);
  MassPointParameters mpp231 (cs231, 58821);
  StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

  MassPoint mp232 (625, 450);
  StopCrossSection cs232 (0.136372, 0.0174504);
  MassPointParameters mpp232 (cs232, 57771);
  StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

  MassPoint mp233 (625, 475);
  StopCrossSection cs233 (0.136372, 0.0174504);
  MassPointParameters mpp233 (cs233, 53659);
  StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

  MassPoint mp234 (625, 500);
  StopCrossSection cs234 (0.136372, 0.0174504);
  MassPointParameters mpp234 (cs234, 51350);
  StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

  MassPoint mp235 (625, 525);
  StopCrossSection cs235 (0.136372, 0.0174504);
  MassPointParameters mpp235 (cs235, 55581);
  StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

  MassPoint mp236 (625, 538);
  StopCrossSection cs236 (0.136372, 0.0174504);
  MassPointParameters mpp236 (cs236, 54603);
  StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

  MassPoint mp237 (633, 450);
  StopCrossSection cs237 (0.125996, 0.0165449);
  MassPointParameters mpp237 (cs237, 48772);
  StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

  MassPoint mp238 (642, 475);
  StopCrossSection cs238 (0.115573, 0.0147355);
  MassPointParameters mpp238 (cs238, 52247);
  StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

  MassPoint mp239 (650, 1);
  StopCrossSection cs239 (0.107045, 0.0138336);
  MassPointParameters mpp239 (cs239, 39561);
  StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

  MassPoint mp240 (650, 50);
  StopCrossSection cs240 (0.107045, 0.0138336);
  MassPointParameters mpp240 (cs240, 41932);
  StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

  MassPoint mp241 (650, 100);
  StopCrossSection cs241 (0.107045, 0.0138336);
  MassPointParameters mpp241 (cs241, 39789);
  StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

  MassPoint mp242 (650, 150);
  StopCrossSection cs242 (0.107045, 0.0138336);
  MassPointParameters mpp242 (cs242, 45693);
  StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

  MassPoint mp243 (650, 200);
  StopCrossSection cs243 (0.107045, 0.0138336);
  MassPointParameters mpp243 (cs243, 40678);
  StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

  MassPoint mp244 (650, 250);
  StopCrossSection cs244 (0.107045, 0.0138336);
  MassPointParameters mpp244 (cs244, 45334);
  StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

  MassPoint mp245 (650, 300);
  StopCrossSection cs245 (0.107045, 0.0138336);
  MassPointParameters mpp245 (cs245, 39559);
  StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

  MassPoint mp246 (650, 350);
  StopCrossSection cs246 (0.107045, 0.0138336);
  MassPointParameters mpp246 (cs246, 42018);
  StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

  MassPoint mp247 (650, 375);
  StopCrossSection cs247 (0.107045, 0.0138336);
  MassPointParameters mpp247 (cs247, 39344);
  StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

  MassPoint mp248 (650, 400);
  StopCrossSection cs248 (0.107045, 0.0138336);
  MassPointParameters mpp248 (cs248, 42883);
  StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

  MassPoint mp249 (650, 425);
  StopCrossSection cs249 (0.107045, 0.0138336);
  MassPointParameters mpp249 (cs249, 40883);
  StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

  MassPoint mp250 (650, 450);
  StopCrossSection cs250 (0.107045, 0.0138336);
  MassPointParameters mpp250 (cs250, 41270);
  StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

  MassPoint mp251 (650, 475);
  StopCrossSection cs251 (0.107045, 0.0138336);
  MassPointParameters mpp251 (cs251, 42292);
  StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

  MassPoint mp252 (650, 500);
  StopCrossSection cs252 (0.107045, 0.0138336);
  MassPointParameters mpp252 (cs252, 39517);
  StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

  MassPoint mp253 (650, 525);
  StopCrossSection cs253 (0.107045, 0.0138336);
  MassPointParameters mpp253 (cs253, 42652);
  StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

  MassPoint mp254 (650, 550);
  StopCrossSection cs254 (0.107045, 0.0138336);
  MassPointParameters mpp254 (cs254, 41964);
  StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

  MassPoint mp255 (650, 563);
  StopCrossSection cs255 (0.107045, 0.0138336);
  MassPointParameters mpp255 (cs255, 40391);
  StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

  MassPoint mp256 (658, 475);
  StopCrossSection cs256 (0.0991824, 0.0128381);
  MassPointParameters mpp256 (cs256, 39792);
  StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

  MassPoint mp257 (667, 500);
  StopCrossSection cs257 (0.0910543, 0.0118196);
  MassPointParameters mpp257 (cs257, 35953);
  StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

  MassPoint mp258 (675, 375);
  StopCrossSection cs258 (0.0844877, 0.0110428);
  MassPointParameters mpp258 (cs258, 32183);
  StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

  MassPoint mp259 (675, 400);
  StopCrossSection cs259 (0.0844877, 0.0110428);
  MassPointParameters mpp259 (cs259, 32162);
  StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

  MassPoint mp260 (675, 425);
  StopCrossSection cs260 (0.0844877, 0.0110428);
  MassPointParameters mpp260 (cs260, 33487);
  StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

  MassPoint mp261 (675, 450);
  StopCrossSection cs261 (0.0844877, 0.0110428);
  MassPointParameters mpp261 (cs261, 33464);
  StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

  MassPoint mp262 (675, 475);
  StopCrossSection cs262 (0.0844877, 0.0110428);
  MassPointParameters mpp262 (cs262, 33683);
  StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

  MassPoint mp263 (675, 500);
  StopCrossSection cs263 (0.0844877, 0.0110428);
  MassPointParameters mpp263 (cs263, 34330);
  StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

  MassPoint mp264 (675, 525);
  StopCrossSection cs264 (0.0844877, 0.0110428);
  MassPointParameters mpp264 (cs264, 37133);
  StopNeutralinoMap.insert(std::make_pair(mp264, mpp264));

  MassPoint mp265 (675, 550);
  StopCrossSection cs265 (0.0844877, 0.0110428);
  MassPointParameters mpp265 (cs265, 30405);
  StopNeutralinoMap.insert(std::make_pair(mp265, mpp265));

  MassPoint mp266 (675, 575);
  StopCrossSection cs266 (0.0844877, 0.0110428);
  MassPointParameters mpp266 (cs266, 32896);
  StopNeutralinoMap.insert(std::make_pair(mp266, mpp266));

  MassPoint mp267 (675, 588);
  StopCrossSection cs267 (0.0844877, 0.0110428);
  MassPointParameters mpp267 (cs267, 33017);
  StopNeutralinoMap.insert(std::make_pair(mp267, mpp267));

  MassPoint mp268 (683, 500);
  StopCrossSection cs268 (0.0783936, 0.0102976);
  MassPointParameters mpp268 (cs268, 25776);
  StopNeutralinoMap.insert(std::make_pair(mp268, mpp268));

  MassPoint mp269 (692, 525);
  StopCrossSection cs269 (0.0721663, 0.00956121);
  MassPointParameters mpp269 (cs269, 25630);
  StopNeutralinoMap.insert(std::make_pair(mp269, mpp269));

  MassPoint mp270 (700, 1);
  StopCrossSection cs270 (0.0670476, 0.00894609);
  MassPointParameters mpp270 (cs270, 26111);
  StopNeutralinoMap.insert(std::make_pair(mp270, mpp270));

  MassPoint mp271 (700, 50);
  StopCrossSection cs271 (0.0670476, 0.00894609);
  MassPointParameters mpp271 (cs271, 27465);
  StopNeutralinoMap.insert(std::make_pair(mp271, mpp271));

  MassPoint mp272 (700, 100);
  StopCrossSection cs272 (0.0670476, 0.00894609);
  MassPointParameters mpp272 (cs272, 25649);
  StopNeutralinoMap.insert(std::make_pair(mp272, mpp272));

  MassPoint mp273 (700, 150);
  StopCrossSection cs273 (0.0670476, 0.00894609);
  MassPointParameters mpp273 (cs273, 20171);
  StopNeutralinoMap.insert(std::make_pair(mp273, mpp273));

  MassPoint mp274 (700, 200);
  StopCrossSection cs274 (0.0670476, 0.00894609);
  MassPointParameters mpp274 (cs274, 22547);
  StopNeutralinoMap.insert(std::make_pair(mp274, mpp274));

  MassPoint mp275 (700, 250);
  StopCrossSection cs275 (0.0670476, 0.00894609);
  MassPointParameters mpp275 (cs275, 26741);
  StopNeutralinoMap.insert(std::make_pair(mp275, mpp275));

  MassPoint mp276 (700, 300);
  StopCrossSection cs276 (0.0670476, 0.00894609);
  MassPointParameters mpp276 (cs276, 28242);
  StopNeutralinoMap.insert(std::make_pair(mp276, mpp276));

  MassPoint mp277 (700, 350);
  StopCrossSection cs277 (0.0670476, 0.00894609);
  MassPointParameters mpp277 (cs277, 23129);
  StopNeutralinoMap.insert(std::make_pair(mp277, mpp277));

  MassPoint mp278 (700, 400);
  StopCrossSection cs278 (0.0670476, 0.00894609);
  MassPointParameters mpp278 (cs278, 23085);
  StopNeutralinoMap.insert(std::make_pair(mp278, mpp278));

  MassPoint mp279 (700, 425);
  StopCrossSection cs279 (0.0670476, 0.00894609);
  MassPointParameters mpp279 (cs279, 25370);
  StopNeutralinoMap.insert(std::make_pair(mp279, mpp279));

  MassPoint mp280 (700, 450);
  StopCrossSection cs280 (0.0670476, 0.00894609);
  MassPointParameters mpp280 (cs280, 24735);
  StopNeutralinoMap.insert(std::make_pair(mp280, mpp280));

  MassPoint mp281 (700, 475);
  StopCrossSection cs281 (0.0670476, 0.00894609);
  MassPointParameters mpp281 (cs281, 26015);
  StopNeutralinoMap.insert(std::make_pair(mp281, mpp281));

  MassPoint mp282 (700, 500);
  StopCrossSection cs282 (0.0670476, 0.00894609);
  MassPointParameters mpp282 (cs282, 26271);
  StopNeutralinoMap.insert(std::make_pair(mp282, mpp282));

  MassPoint mp283 (700, 525);
  StopCrossSection cs283 (0.0670476, 0.00894609);
  MassPointParameters mpp283 (cs283, 26332);
  StopNeutralinoMap.insert(std::make_pair(mp283, mpp283));

  MassPoint mp284 (700, 550);
  StopCrossSection cs284 (0.0670476, 0.00894609);
  MassPointParameters mpp284 (cs284, 20934);
  StopNeutralinoMap.insert(std::make_pair(mp284, mpp284));

  MassPoint mp285 (700, 575);
  StopCrossSection cs285 (0.0670476, 0.00894609);
  MassPointParameters mpp285 (cs285, 27086);
  StopNeutralinoMap.insert(std::make_pair(mp285, mpp285));

  MassPoint mp286 (700, 600);
  StopCrossSection cs286 (0.0670476, 0.00894609);
  MassPointParameters mpp286 (cs286, 24494);
  StopNeutralinoMap.insert(std::make_pair(mp286, mpp286));

  MassPoint mp287 (700, 613);
  StopCrossSection cs287 (0.0670476, 0.00894609);
  MassPointParameters mpp287 (cs287, 26506);
  StopNeutralinoMap.insert(std::make_pair(mp287, mpp287));

  MassPoint mp288 (708, 525);
  StopCrossSection cs288 (0.0624336, 0.00835443);
  MassPointParameters mpp288 (cs288, 27069);
  StopNeutralinoMap.insert(std::make_pair(mp288, mpp288));

  MassPoint mp289 (717, 550);
  StopCrossSection cs289 (0.0575708, 0.00775986);
  MassPointParameters mpp289 (cs289, 23062);
  StopNeutralinoMap.insert(std::make_pair(mp289, mpp289));

  MassPoint mp290 (725, 425);
  StopCrossSection cs290 (0.0536438, 0.00728504);
  MassPointParameters mpp290 (cs290, 18160);
  StopNeutralinoMap.insert(std::make_pair(mp290, mpp290));

  MassPoint mp291 (725, 450);
  StopCrossSection cs291 (0.0536438, 0.00728504);
  MassPointParameters mpp291 (cs291, 22343);
  StopNeutralinoMap.insert(std::make_pair(mp291, mpp291));

  MassPoint mp292 (725, 475);
  StopCrossSection cs292 (0.0536438, 0.00728504);
  MassPointParameters mpp292 (cs292, 20603);
  StopNeutralinoMap.insert(std::make_pair(mp292, mpp292));

  MassPoint mp293 (725, 500);
  StopCrossSection cs293 (0.0536438, 0.00728504);
  MassPointParameters mpp293 (cs293, 22427);
  StopNeutralinoMap.insert(std::make_pair(mp293, mpp293));

  MassPoint mp294 (725, 525);
  StopCrossSection cs294 (0.0536438, 0.00728504);
  MassPointParameters mpp294 (cs294, 24821);
  StopNeutralinoMap.insert(std::make_pair(mp294, mpp294));

  MassPoint mp295 (725, 550);
  StopCrossSection cs295 (0.0536438, 0.00728504);
  MassPointParameters mpp295 (cs295, 22408);
  StopNeutralinoMap.insert(std::make_pair(mp295, mpp295));

  MassPoint mp296 (725, 575);
  StopCrossSection cs296 (0.0536438, 0.00728504);
  MassPointParameters mpp296 (cs296, 19038);
  StopNeutralinoMap.insert(std::make_pair(mp296, mpp296));

  MassPoint mp297 (725, 600);
  StopCrossSection cs297 (0.0536438, 0.00728504);
  MassPointParameters mpp297 (cs297, 22627);
  StopNeutralinoMap.insert(std::make_pair(mp297, mpp297));

  MassPoint mp298 (725, 625);
  StopCrossSection cs298 (0.0536438, 0.00728504);
  MassPointParameters mpp298 (cs298, 22703);
  StopNeutralinoMap.insert(std::make_pair(mp298, mpp298));

  MassPoint mp299 (725, 638);
  StopCrossSection cs299 (0.0536438, 0.00728504);
  MassPointParameters mpp299 (cs299, 22341);
  StopNeutralinoMap.insert(std::make_pair(mp299, mpp299));

  MassPoint mp300 (733, 550);
  StopCrossSection cs300 (0.0499888, 0.00679985);
  MassPointParameters mpp300 (cs300, 17080);
  StopNeutralinoMap.insert(std::make_pair(mp300, mpp300));

  MassPoint mp301 (742, 575);
  StopCrossSection cs301 (0.0462725, 0.00633304);
  MassPointParameters mpp301 (cs301, 21237);
  StopNeutralinoMap.insert(std::make_pair(mp301, mpp301));

  MassPoint mp302 (750, 1);
  StopCrossSection cs302 (0.0431418, 0.00593006);
  MassPointParameters mpp302 (cs302, 18315);
  StopNeutralinoMap.insert(std::make_pair(mp302, mpp302));

  MassPoint mp303 (750, 50);
  StopCrossSection cs303 (0.0431418, 0.00593006);
  MassPointParameters mpp303 (cs303, 18103);
  StopNeutralinoMap.insert(std::make_pair(mp303, mpp303));

  MassPoint mp304 (750, 100);
  StopCrossSection cs304 (0.0431418, 0.00593006);
  MassPointParameters mpp304 (cs304, 22648);
  StopNeutralinoMap.insert(std::make_pair(mp304, mpp304));

  MassPoint mp305 (750, 150);
  StopCrossSection cs305 (0.0431418, 0.00593006);
  MassPointParameters mpp305 (cs305, 18839);
  StopNeutralinoMap.insert(std::make_pair(mp305, mpp305));

  MassPoint mp306 (750, 200);
  StopCrossSection cs306 (0.0431418, 0.00593006);
  MassPointParameters mpp306 (cs306, 18964);
  StopNeutralinoMap.insert(std::make_pair(mp306, mpp306));

  MassPoint mp307 (750, 250);
  StopCrossSection cs307 (0.0431418, 0.00593006);
  MassPointParameters mpp307 (cs307, 17522);
  StopNeutralinoMap.insert(std::make_pair(mp307, mpp307));

  MassPoint mp308 (750, 300);
  StopCrossSection cs308 (0.0431418, 0.00593006);
  MassPointParameters mpp308 (cs308, 21378);
  StopNeutralinoMap.insert(std::make_pair(mp308, mpp308));

  MassPoint mp309 (750, 350);
  StopCrossSection cs309 (0.0431418, 0.00593006);
  MassPointParameters mpp309 (cs309, 17361);
  StopNeutralinoMap.insert(std::make_pair(mp309, mpp309));

  MassPoint mp310 (750, 400);
  StopCrossSection cs310 (0.0431418, 0.00593006);
  MassPointParameters mpp310 (cs310, 18457);
  StopNeutralinoMap.insert(std::make_pair(mp310, mpp310));

  MassPoint mp311 (750, 450);
  StopCrossSection cs311 (0.0431418, 0.00593006);
  MassPointParameters mpp311 (cs311, 19891);
  StopNeutralinoMap.insert(std::make_pair(mp311, mpp311));

  MassPoint mp312 (750, 475);
  StopCrossSection cs312 (0.0431418, 0.00593006);
  MassPointParameters mpp312 (cs312, 18658);
  StopNeutralinoMap.insert(std::make_pair(mp312, mpp312));

  MassPoint mp313 (750, 500);
  StopCrossSection cs313 (0.0431418, 0.00593006);
  MassPointParameters mpp313 (cs313, 19258);
  StopNeutralinoMap.insert(std::make_pair(mp313, mpp313));

  MassPoint mp314 (750, 525);
  StopCrossSection cs314 (0.0431418, 0.00593006);
  MassPointParameters mpp314 (cs314, 15547);
  StopNeutralinoMap.insert(std::make_pair(mp314, mpp314));

  MassPoint mp315 (750, 550);
  StopCrossSection cs315 (0.0431418, 0.00593006);
  MassPointParameters mpp315 (cs315, 20206);
  StopNeutralinoMap.insert(std::make_pair(mp315, mpp315));

  MassPoint mp316 (750, 575);
  StopCrossSection cs316 (0.0431418, 0.00593006);
  MassPointParameters mpp316 (cs316, 20945);
  StopNeutralinoMap.insert(std::make_pair(mp316, mpp316));

  MassPoint mp317 (750, 600);
  StopCrossSection cs317 (0.0431418, 0.00593006);
  MassPointParameters mpp317 (cs317, 20420);
  StopNeutralinoMap.insert(std::make_pair(mp317, mpp317));

  MassPoint mp318 (750, 625);
  StopCrossSection cs318 (0.0431418, 0.00593006);
  MassPointParameters mpp318 (cs318, 23907);
  StopNeutralinoMap.insert(std::make_pair(mp318, mpp318));

  MassPoint mp319 (750, 650);
  StopCrossSection cs319 (0.0431418, 0.00593006);
  MassPointParameters mpp319 (cs319, 16680);
  StopNeutralinoMap.insert(std::make_pair(mp319, mpp319));

  MassPoint mp320 (758, 575);
  StopCrossSection cs320 (0.0403137, 0.00557285);
  MassPointParameters mpp320 (cs320, 22200);
  StopNeutralinoMap.insert(std::make_pair(mp320, mpp320));

  MassPoint mp321 (767, 600);
  StopCrossSection cs321 (0.0372964, 0.00517853);
  MassPointParameters mpp321 (cs321, 19945);
  StopNeutralinoMap.insert(std::make_pair(mp321, mpp321));

  MassPoint mp322 (775, 475);
  StopCrossSection cs322 (0.0348796, 0.00486909);
  MassPointParameters mpp322 (cs322, 18999);
  StopNeutralinoMap.insert(std::make_pair(mp322, mpp322));

  MassPoint mp323 (775, 500);
  StopCrossSection cs323 (0.0348796, 0.00486909);
  MassPointParameters mpp323 (cs323, 17163);
  StopNeutralinoMap.insert(std::make_pair(mp323, mpp323));

  MassPoint mp324 (775, 525);
  StopCrossSection cs324 (0.0348796, 0.00486909);
  MassPointParameters mpp324 (cs324, 17133);
  StopNeutralinoMap.insert(std::make_pair(mp324, mpp324));

  MassPoint mp325 (775, 550);
  StopCrossSection cs325 (0.0348796, 0.00486909);
  MassPointParameters mpp325 (cs325, 15054);
  StopNeutralinoMap.insert(std::make_pair(mp325, mpp325));

  MassPoint mp326 (775, 575);
  StopCrossSection cs326 (0.0348796, 0.00486909);
  MassPointParameters mpp326 (cs326, 19633);
  StopNeutralinoMap.insert(std::make_pair(mp326, mpp326));

  MassPoint mp327 (775, 600);
  StopCrossSection cs327 (0.0348796, 0.00486909);
  MassPointParameters mpp327 (cs327, 17283);
  StopNeutralinoMap.insert(std::make_pair(mp327, mpp327));

  MassPoint mp328 (775, 625);
  StopCrossSection cs328 (0.0348796, 0.00486909);
  MassPointParameters mpp328 (cs328, 20492);
  StopNeutralinoMap.insert(std::make_pair(mp328, mpp328));

  MassPoint mp329 (775, 650);
  StopCrossSection cs329 (0.0348796, 0.00486909);
  MassPointParameters mpp329 (cs329, 18582);
  StopNeutralinoMap.insert(std::make_pair(mp329, mpp329));

  MassPoint mp330 (783, 600);
  StopCrossSection cs330 (0.0326196, 0.00457813);
  MassPointParameters mpp330 (cs330, 18739);
  StopNeutralinoMap.insert(std::make_pair(mp330, mpp330));

  MassPoint mp331 (792, 625);
  StopCrossSection cs331 (0.0302563, 0.00427359);
  MassPointParameters mpp331 (cs331, 16790);
  StopNeutralinoMap.insert(std::make_pair(mp331, mpp331));

  MassPoint mp332 (800, 1);
  StopCrossSection cs332 (0.0283338, 0.00401518);
  MassPointParameters mpp332 (cs332, 18533);
  StopNeutralinoMap.insert(std::make_pair(mp332, mpp332));

  MassPoint mp333 (800, 50);
  StopCrossSection cs333 (0.0283338, 0.00401518);
  MassPointParameters mpp333 (cs333, 18635);
  StopNeutralinoMap.insert(std::make_pair(mp333, mpp333));

  MassPoint mp334 (800, 100);
  StopCrossSection cs334 (0.0283338, 0.00401518);
  MassPointParameters mpp334 (cs334, 23126);
  StopNeutralinoMap.insert(std::make_pair(mp334, mpp334));

  MassPoint mp335 (800, 150);
  StopCrossSection cs335 (0.0283338, 0.00401518);
  MassPointParameters mpp335 (cs335, 19048);
  StopNeutralinoMap.insert(std::make_pair(mp335, mpp335));

  MassPoint mp336 (800, 200);
  StopCrossSection cs336 (0.0283338, 0.00401518);
  MassPointParameters mpp336 (cs336, 18703);
  StopNeutralinoMap.insert(std::make_pair(mp336, mpp336));

  MassPoint mp337 (800, 250);
  StopCrossSection cs337 (0.0283338, 0.00401518);
  MassPointParameters mpp337 (cs337, 22584);
  StopNeutralinoMap.insert(std::make_pair(mp337, mpp337));

  MassPoint mp338 (800, 300);
  StopCrossSection cs338 (0.0283338, 0.00401518);
  MassPointParameters mpp338 (cs338, 16831);
  StopNeutralinoMap.insert(std::make_pair(mp338, mpp338));

  MassPoint mp339 (800, 350);
  StopCrossSection cs339 (0.0283338, 0.00401518);
  MassPointParameters mpp339 (cs339, 18809);
  StopNeutralinoMap.insert(std::make_pair(mp339, mpp339));

  MassPoint mp340 (800, 400);
  StopCrossSection cs340 (0.0283338, 0.00401518);
  MassPointParameters mpp340 (cs340, 19428);
  StopNeutralinoMap.insert(std::make_pair(mp340, mpp340));

  MassPoint mp341 (800, 450);
  StopCrossSection cs341 (0.0283338, 0.00401518);
  MassPointParameters mpp341 (cs341, 19534);
  StopNeutralinoMap.insert(std::make_pair(mp341, mpp341));

  MassPoint mp342 (800, 500);
  StopCrossSection cs342 (0.0283338, 0.00401518);
  MassPointParameters mpp342 (cs342, 21312);
  StopNeutralinoMap.insert(std::make_pair(mp342, mpp342));

  MassPoint mp343 (800, 525);
  StopCrossSection cs343 (0.0283338, 0.00401518);
  MassPointParameters mpp343 (cs343, 17702);
  StopNeutralinoMap.insert(std::make_pair(mp343, mpp343));

  MassPoint mp344 (800, 550);
  StopCrossSection cs344 (0.0283338, 0.00401518);
  MassPointParameters mpp344 (cs344, 18987);
  StopNeutralinoMap.insert(std::make_pair(mp344, mpp344));

  MassPoint mp345 (800, 575);
  StopCrossSection cs345 (0.0283338, 0.00401518);
  MassPointParameters mpp345 (cs345, 15937);
  StopNeutralinoMap.insert(std::make_pair(mp345, mpp345));

  MassPoint mp346 (800, 600);
  StopCrossSection cs346 (0.0283338, 0.00401518);
  MassPointParameters mpp346 (cs346, 20128);
  StopNeutralinoMap.insert(std::make_pair(mp346, mpp346));

  MassPoint mp347 (800, 625);
  StopCrossSection cs347 (0.0283338, 0.00401518);
  MassPointParameters mpp347 (cs347, 14755);
  StopNeutralinoMap.insert(std::make_pair(mp347, mpp347));

  MassPoint mp348 (800, 650);
  StopCrossSection cs348 (0.0283338, 0.00401518);
  MassPointParameters mpp348 (cs348, 19720);
  StopNeutralinoMap.insert(std::make_pair(mp348, mpp348));

  MassPoint mp349 (808, 625);
  StopCrossSection cs349 (0.0265622, 0.00379026);
  MassPointParameters mpp349 (cs349, 21277);
  StopNeutralinoMap.insert(std::make_pair(mp349, mpp349));

  MassPoint mp350 (817, 650);
  StopCrossSection cs350 (0.0247104, 0.00355087);
  MassPointParameters mpp350 (cs350, 19143);
  StopNeutralinoMap.insert(std::make_pair(mp350, mpp350));

  MassPoint mp351 (825, 525);
  StopCrossSection cs351 (0.0230866, 0.00333435);
  MassPointParameters mpp351 (cs351, 21129);
  StopNeutralinoMap.insert(std::make_pair(mp351, mpp351));

  MassPoint mp352 (825, 550);
  StopCrossSection cs352 (0.0230866, 0.00333435);
  MassPointParameters mpp352 (cs352, 19159);
  StopNeutralinoMap.insert(std::make_pair(mp352, mpp352));

  MassPoint mp353 (825, 575);
  StopCrossSection cs353 (0.0230866, 0.00333435);
  MassPointParameters mpp353 (cs353, 18458);
  StopNeutralinoMap.insert(std::make_pair(mp353, mpp353));

  MassPoint mp354 (825, 600);
  StopCrossSection cs354 (0.0230866, 0.00333435);
  MassPointParameters mpp354 (cs354, 20323);
  StopNeutralinoMap.insert(std::make_pair(mp354, mpp354));

  MassPoint mp355 (825, 625);
  StopCrossSection cs355 (0.0230866, 0.00333435);
  MassPointParameters mpp355 (cs355, 19352);
  StopNeutralinoMap.insert(std::make_pair(mp355, mpp355));

  MassPoint mp356 (825, 650);
  StopCrossSection cs356 (0.0230866, 0.00333435);
  MassPointParameters mpp356 (cs356, 19509);
  StopNeutralinoMap.insert(std::make_pair(mp356, mpp356));

  MassPoint mp357 (833, 650);
  StopCrossSection cs357 (0.0216993, 0.0031511);
  MassPointParameters mpp357 (cs357, 21277);
  StopNeutralinoMap.insert(std::make_pair(mp357, mpp357));

  MassPoint mp358 (850, 1);
  StopCrossSection cs358 (0.0189612, 0.00278768);
  MassPointParameters mpp358 (cs358, 18872);
  StopNeutralinoMap.insert(std::make_pair(mp358, mpp358));

  MassPoint mp359 (850, 50);
  StopCrossSection cs359 (0.0189612, 0.00278768);
  MassPointParameters mpp359 (cs359, 21196);
  StopNeutralinoMap.insert(std::make_pair(mp359, mpp359));

  MassPoint mp360 (850, 100);
  StopCrossSection cs360 (0.0189612, 0.00278768);
  MassPointParameters mpp360 (cs360, 18684);
  StopNeutralinoMap.insert(std::make_pair(mp360, mpp360));

  MassPoint mp361 (850, 150);
  StopCrossSection cs361 (0.0189612, 0.00278768);
  MassPointParameters mpp361 (cs361, 19904);
  StopNeutralinoMap.insert(std::make_pair(mp361, mpp361));

  MassPoint mp362 (850, 200);
  StopCrossSection cs362 (0.0189612, 0.00278768);
  MassPointParameters mpp362 (cs362, 17242);
  StopNeutralinoMap.insert(std::make_pair(mp362, mpp362));

  MassPoint mp363 (850, 250);
  StopCrossSection cs363 (0.0189612, 0.00278768);
  MassPointParameters mpp363 (cs363, 17776);
  StopNeutralinoMap.insert(std::make_pair(mp363, mpp363));

  MassPoint mp364 (850, 300);
  StopCrossSection cs364 (0.0189612, 0.00278768);
  MassPointParameters mpp364 (cs364, 18397);
  StopNeutralinoMap.insert(std::make_pair(mp364, mpp364));

  MassPoint mp365 (850, 350);
  StopCrossSection cs365 (0.0189612, 0.00278768);
  MassPointParameters mpp365 (cs365, 19063);
  StopNeutralinoMap.insert(std::make_pair(mp365, mpp365));

  MassPoint mp366 (850, 400);
  StopCrossSection cs366 (0.0189612, 0.00278768);
  MassPointParameters mpp366 (cs366, 16887);
  StopNeutralinoMap.insert(std::make_pair(mp366, mpp366));

  MassPoint mp367 (850, 450);
  StopCrossSection cs367 (0.0189612, 0.00278768);
  MassPointParameters mpp367 (cs367, 21780);
  StopNeutralinoMap.insert(std::make_pair(mp367, mpp367));

  MassPoint mp368 (850, 500);
  StopCrossSection cs368 (0.0189612, 0.00278768);
  MassPointParameters mpp368 (cs368, 21169);
  StopNeutralinoMap.insert(std::make_pair(mp368, mpp368));

  MassPoint mp369 (850, 550);
  StopCrossSection cs369 (0.0189612, 0.00278768);
  MassPointParameters mpp369 (cs369, 18999);
  StopNeutralinoMap.insert(std::make_pair(mp369, mpp369));

  MassPoint mp370 (850, 575);
  StopCrossSection cs370 (0.0189612, 0.00278768);
  MassPointParameters mpp370 (cs370, 16400);
  StopNeutralinoMap.insert(std::make_pair(mp370, mpp370));

  MassPoint mp371 (850, 600);
  StopCrossSection cs371 (0.0189612, 0.00278768);
  MassPointParameters mpp371 (cs371, 14517);
  StopNeutralinoMap.insert(std::make_pair(mp371, mpp371));

  MassPoint mp372 (850, 625);
  StopCrossSection cs372 (0.0189612, 0.00278768);
  MassPointParameters mpp372 (cs372, 18129);
  StopNeutralinoMap.insert(std::make_pair(mp372, mpp372));

  MassPoint mp373 (850, 650);
  StopCrossSection cs373 (0.0189612, 0.00278768);
  MassPointParameters mpp373 (cs373, 19614);
  StopNeutralinoMap.insert(std::make_pair(mp373, mpp373));

  MassPoint mp374 (875, 575);
  StopCrossSection cs374 (0.015625, 0.00233698);
  MassPointParameters mpp374 (cs374, 17851);
  StopNeutralinoMap.insert(std::make_pair(mp374, mpp374));

  MassPoint mp375 (875, 600);
  StopCrossSection cs375 (0.015625, 0.00233698);
  MassPointParameters mpp375 (cs375, 20189);
  StopNeutralinoMap.insert(std::make_pair(mp375, mpp375));

  MassPoint mp376 (875, 625);
  StopCrossSection cs376 (0.015625, 0.00233698);
  MassPointParameters mpp376 (cs376, 18632);
  StopNeutralinoMap.insert(std::make_pair(mp376, mpp376));

  MassPoint mp377 (875, 650);
  StopCrossSection cs377 (0.015625, 0.00233698);
  MassPointParameters mpp377 (cs377, 18532);
  StopNeutralinoMap.insert(std::make_pair(mp377, mpp377));

  MassPoint mp378 (900, 1);
  StopCrossSection cs378 (0.0128895, 0.00195954);
  MassPointParameters mpp378 (cs378, 19667);
  StopNeutralinoMap.insert(std::make_pair(mp378, mpp378));

  MassPoint mp379 (900, 50);
  StopCrossSection cs379 (0.0128895, 0.00195954);
  MassPointParameters mpp379 (cs379, 21651);
  StopNeutralinoMap.insert(std::make_pair(mp379, mpp379));

  MassPoint mp380 (900, 100);
  StopCrossSection cs380 (0.0128895, 0.00195954);
  MassPointParameters mpp380 (cs380, 18211);
  StopNeutralinoMap.insert(std::make_pair(mp380, mpp380));

  MassPoint mp381 (900, 150);
  StopCrossSection cs381 (0.0128895, 0.00195954);
  MassPointParameters mpp381 (cs381, 17432);
  StopNeutralinoMap.insert(std::make_pair(mp381, mpp381));

  MassPoint mp382 (900, 200);
  StopCrossSection cs382 (0.0128895, 0.00195954);
  MassPointParameters mpp382 (cs382, 19494);
  StopNeutralinoMap.insert(std::make_pair(mp382, mpp382));

  MassPoint mp383 (900, 250);
  StopCrossSection cs383 (0.0128895, 0.00195954);
  MassPointParameters mpp383 (cs383, 16297);
  StopNeutralinoMap.insert(std::make_pair(mp383, mpp383));

  MassPoint mp384 (900, 300);
  StopCrossSection cs384 (0.0128895, 0.00195954);
  MassPointParameters mpp384 (cs384, 19115);
  StopNeutralinoMap.insert(std::make_pair(mp384, mpp384));

  MassPoint mp385 (900, 350);
  StopCrossSection cs385 (0.0128895, 0.00195954);
  MassPointParameters mpp385 (cs385, 19783);
  StopNeutralinoMap.insert(std::make_pair(mp385, mpp385));

  MassPoint mp386 (900, 400);
  StopCrossSection cs386 (0.0128895, 0.00195954);
  MassPointParameters mpp386 (cs386, 19405);
  StopNeutralinoMap.insert(std::make_pair(mp386, mpp386));

  MassPoint mp387 (900, 450);
  StopCrossSection cs387 (0.0128895, 0.00195954);
  MassPointParameters mpp387 (cs387, 19605);
  StopNeutralinoMap.insert(std::make_pair(mp387, mpp387));

  MassPoint mp388 (900, 500);
  StopCrossSection cs388 (0.0128895, 0.00195954);
  MassPointParameters mpp388 (cs388, 17254);
  StopNeutralinoMap.insert(std::make_pair(mp388, mpp388));

  MassPoint mp389 (900, 550);
  StopCrossSection cs389 (0.0128895, 0.00195954);
  MassPointParameters mpp389 (cs389, 16796);
  StopNeutralinoMap.insert(std::make_pair(mp389, mpp389));

  MassPoint mp390 (900, 600);
  StopCrossSection cs390 (0.0128895, 0.00195954);
  MassPointParameters mpp390 (cs390, 18222);
  StopNeutralinoMap.insert(std::make_pair(mp390, mpp390));

  MassPoint mp391 (900, 625);
  StopCrossSection cs391 (0.0128895, 0.00195954);
  MassPointParameters mpp391 (cs391, 18390);
  StopNeutralinoMap.insert(std::make_pair(mp391, mpp391));

  MassPoint mp392 (900, 650);
  StopCrossSection cs392 (0.0128895, 0.00195954);
  MassPointParameters mpp392 (cs392, 17897);
  StopNeutralinoMap.insert(std::make_pair(mp392, mpp392));

  MassPoint mp393 (925, 625);
  StopCrossSection cs393 (0.0106631, 0.00165071);
  MassPointParameters mpp393 (cs393, 21218);
  StopNeutralinoMap.insert(std::make_pair(mp393, mpp393));

  MassPoint mp394 (925, 650);
  StopCrossSection cs394 (0.0106631, 0.00165071);
  MassPointParameters mpp394 (cs394, 20376);
  StopNeutralinoMap.insert(std::make_pair(mp394, mpp394));

  MassPoint mp395 (950, 1);
  StopCrossSection cs395 (0.00883465, 0.0013886);
  MassPointParameters mpp395 (cs395, 23521);
  StopNeutralinoMap.insert(std::make_pair(mp395, mpp395));

  MassPoint mp396 (950, 50);
  StopCrossSection cs396 (0.00883465, 0.0013886);
  MassPointParameters mpp396 (cs396, 20772);
  StopNeutralinoMap.insert(std::make_pair(mp396, mpp396));

  MassPoint mp397 (950, 100);
  StopCrossSection cs397 (0.00883465, 0.0013886);
  MassPointParameters mpp397 (cs397, 20889);
  StopNeutralinoMap.insert(std::make_pair(mp397, mpp397));

  MassPoint mp398 (950, 150);
  StopCrossSection cs398 (0.00883465, 0.0013886);
  MassPointParameters mpp398 (cs398, 17523);
  StopNeutralinoMap.insert(std::make_pair(mp398, mpp398));

  MassPoint mp399 (950, 200);
  StopCrossSection cs399 (0.00883465, 0.0013886);
  MassPointParameters mpp399 (cs399, 20573);
  StopNeutralinoMap.insert(std::make_pair(mp399, mpp399));

  MassPoint mp400 (950, 250);
  StopCrossSection cs400 (0.00883465, 0.0013886);
  MassPointParameters mpp400 (cs400, 21709);
  StopNeutralinoMap.insert(std::make_pair(mp400, mpp400));

  MassPoint mp401 (950, 300);
  StopCrossSection cs401 (0.00883465, 0.0013886);
  MassPointParameters mpp401 (cs401, 17749);
  StopNeutralinoMap.insert(std::make_pair(mp401, mpp401));

  MassPoint mp402 (950, 350);
  StopCrossSection cs402 (0.00883465, 0.0013886);
  MassPointParameters mpp402 (cs402, 17035);
  StopNeutralinoMap.insert(std::make_pair(mp402, mpp402));

  MassPoint mp403 (950, 400);
  StopCrossSection cs403 (0.00883465, 0.0013886);
  MassPointParameters mpp403 (cs403, 15832);
  StopNeutralinoMap.insert(std::make_pair(mp403, mpp403));

  MassPoint mp404 (950, 450);
  StopCrossSection cs404 (0.00883465, 0.0013886);
  MassPointParameters mpp404 (cs404, 15418);
  StopNeutralinoMap.insert(std::make_pair(mp404, mpp404));

  MassPoint mp405 (950, 500);
  StopCrossSection cs405 (0.00883465, 0.0013886);
  MassPointParameters mpp405 (cs405, 18750);
  StopNeutralinoMap.insert(std::make_pair(mp405, mpp405));

  MassPoint mp406 (950, 550);
  StopCrossSection cs406 (0.00883465, 0.0013886);
  MassPointParameters mpp406 (cs406, 18086);
  StopNeutralinoMap.insert(std::make_pair(mp406, mpp406));

  MassPoint mp407 (950, 600);
  StopCrossSection cs407 (0.00883465, 0.0013886);
  MassPointParameters mpp407 (cs407, 18997);
  StopNeutralinoMap.insert(std::make_pair(mp407, mpp407));

  MassPoint mp408 (950, 650);
  StopCrossSection cs408 (0.00883465, 0.0013886);
  MassPointParameters mpp408 (cs408, 19072);
  StopNeutralinoMap.insert(std::make_pair(mp408, mpp408));

  MassPoint mp409 (1000, 1);
  StopCrossSection cs409 (0.00615134, 0.00100238);
  MassPointParameters mpp409 (cs409, 19536);
  StopNeutralinoMap.insert(std::make_pair(mp409, mpp409));

  MassPoint mp410 (1000, 50);
  StopCrossSection cs410 (0.00615134, 0.00100238);
  MassPointParameters mpp410 (cs410, 17885);
  StopNeutralinoMap.insert(std::make_pair(mp410, mpp410));

  MassPoint mp411 (1000, 100);
  StopCrossSection cs411 (0.00615134, 0.00100238);
  MassPointParameters mpp411 (cs411, 19056);
  StopNeutralinoMap.insert(std::make_pair(mp411, mpp411));

  MassPoint mp412 (1000, 150);
  StopCrossSection cs412 (0.00615134, 0.00100238);
  MassPointParameters mpp412 (cs412, 23375);
  StopNeutralinoMap.insert(std::make_pair(mp412, mpp412));

  MassPoint mp413 (1000, 200);
  StopCrossSection cs413 (0.00615134, 0.00100238);
  MassPointParameters mpp413 (cs413, 15383);
  StopNeutralinoMap.insert(std::make_pair(mp413, mpp413));

  MassPoint mp414 (1000, 250);
  StopCrossSection cs414 (0.00615134, 0.00100238);
  MassPointParameters mpp414 (cs414, 19039);
  StopNeutralinoMap.insert(std::make_pair(mp414, mpp414));

  MassPoint mp415 (1000, 300);
  StopCrossSection cs415 (0.00615134, 0.00100238);
  MassPointParameters mpp415 (cs415, 14698);
  StopNeutralinoMap.insert(std::make_pair(mp415, mpp415));

  MassPoint mp416 (1000, 350);
  StopCrossSection cs416 (0.00615134, 0.00100238);
  MassPointParameters mpp416 (cs416, 17436);
  StopNeutralinoMap.insert(std::make_pair(mp416, mpp416));

  MassPoint mp417 (1000, 400);
  StopCrossSection cs417 (0.00615134, 0.00100238);
  MassPointParameters mpp417 (cs417, 16450);
  StopNeutralinoMap.insert(std::make_pair(mp417, mpp417));

  MassPoint mp418 (1000, 450);
  StopCrossSection cs418 (0.00615134, 0.00100238);
  MassPointParameters mpp418 (cs418, 15813);
  StopNeutralinoMap.insert(std::make_pair(mp418, mpp418));

  MassPoint mp419 (1000, 500);
  StopCrossSection cs419 (0.00615134, 0.00100238);
  MassPointParameters mpp419 (cs419, 17173);
  StopNeutralinoMap.insert(std::make_pair(mp419, mpp419));

  MassPoint mp420 (1000, 550);
  StopCrossSection cs420 (0.00615134, 0.00100238);
  MassPointParameters mpp420 (cs420, 22344);
  StopNeutralinoMap.insert(std::make_pair(mp420, mpp420));

  MassPoint mp421 (1000, 600);
  StopCrossSection cs421 (0.00615134, 0.00100238);
  MassPointParameters mpp421 (cs421, 20898);
  StopNeutralinoMap.insert(std::make_pair(mp421, mpp421));

  MassPoint mp422 (1000, 650);
  StopCrossSection cs422 (0.00615134, 0.00100238);
  MassPointParameters mpp422 (cs422, 25491);
  StopNeutralinoMap.insert(std::make_pair(mp422, mpp422));

  MassPoint mp423 (1050, 1);
  StopCrossSection cs423 (0.00432261, 0.000725589);
  MassPointParameters mpp423 (cs423, 18683);
  StopNeutralinoMap.insert(std::make_pair(mp423, mpp423));

  MassPoint mp424 (1050, 50);
  StopCrossSection cs424 (0.00432261, 0.000725589);
  MassPointParameters mpp424 (cs424, 15579);
  StopNeutralinoMap.insert(std::make_pair(mp424, mpp424));

  MassPoint mp425 (1050, 100);
  StopCrossSection cs425 (0.00432261, 0.000725589);
  MassPointParameters mpp425 (cs425, 14860);
  StopNeutralinoMap.insert(std::make_pair(mp425, mpp425));

  MassPoint mp426 (1050, 150);
  StopCrossSection cs426 (0.00432261, 0.000725589);
  MassPointParameters mpp426 (cs426, 20257);
  StopNeutralinoMap.insert(std::make_pair(mp426, mpp426));

  MassPoint mp427 (1050, 200);
  StopCrossSection cs427 (0.00432261, 0.000725589);
  MassPointParameters mpp427 (cs427, 17442);
  StopNeutralinoMap.insert(std::make_pair(mp427, mpp427));

  MassPoint mp428 (1050, 250);
  StopCrossSection cs428 (0.00432261, 0.000725589);
  MassPointParameters mpp428 (cs428, 19038);
  StopNeutralinoMap.insert(std::make_pair(mp428, mpp428));

  MassPoint mp429 (1050, 300);
  StopCrossSection cs429 (0.00432261, 0.000725589);
  MassPointParameters mpp429 (cs429, 19647);
  StopNeutralinoMap.insert(std::make_pair(mp429, mpp429));

  MassPoint mp430 (1050, 350);
  StopCrossSection cs430 (0.00432261, 0.000725589);
  MassPointParameters mpp430 (cs430, 18029);
  StopNeutralinoMap.insert(std::make_pair(mp430, mpp430));

  MassPoint mp431 (1050, 400);
  StopCrossSection cs431 (0.00432261, 0.000725589);
  MassPointParameters mpp431 (cs431, 18420);
  StopNeutralinoMap.insert(std::make_pair(mp431, mpp431));

  MassPoint mp432 (1050, 450);
  StopCrossSection cs432 (0.00432261, 0.000725589);
  MassPointParameters mpp432 (cs432, 18297);
  StopNeutralinoMap.insert(std::make_pair(mp432, mpp432));

  MassPoint mp433 (1050, 500);
  StopCrossSection cs433 (0.00432261, 0.000725589);
  MassPointParameters mpp433 (cs433, 16860);
  StopNeutralinoMap.insert(std::make_pair(mp433, mpp433));

  MassPoint mp434 (1050, 550);
  StopCrossSection cs434 (0.00432261, 0.000725589);
  MassPointParameters mpp434 (cs434, 18895);
  StopNeutralinoMap.insert(std::make_pair(mp434, mpp434));

  MassPoint mp435 (1050, 600);
  StopCrossSection cs435 (0.00432261, 0.000725589);
  MassPointParameters mpp435 (cs435, 16704);
  StopNeutralinoMap.insert(std::make_pair(mp435, mpp435));

  MassPoint mp436 (1050, 650);
  StopCrossSection cs436 (0.00432261, 0.000725589);
  MassPointParameters mpp436 (cs436, 18113);
  StopNeutralinoMap.insert(std::make_pair(mp436, mpp436));

  MassPoint mp437 (1100, 1);
  StopCrossSection cs437 (0.00307413, 0.000532983);
  MassPointParameters mpp437 (cs437, 17516);
  StopNeutralinoMap.insert(std::make_pair(mp437, mpp437));

  MassPoint mp438 (1100, 50);
  StopCrossSection cs438 (0.00307413, 0.000532983);
  MassPointParameters mpp438 (cs438, 19165);
  StopNeutralinoMap.insert(std::make_pair(mp438, mpp438));

  MassPoint mp439 (1100, 100);
  StopCrossSection cs439 (0.00307413, 0.000532983);
  MassPointParameters mpp439 (cs439, 18053);
  StopNeutralinoMap.insert(std::make_pair(mp439, mpp439));

  MassPoint mp440 (1100, 150);
  StopCrossSection cs440 (0.00307413, 0.000532983);
  MassPointParameters mpp440 (cs440, 21117);
  StopNeutralinoMap.insert(std::make_pair(mp440, mpp440));

  MassPoint mp441 (1100, 200);
  StopCrossSection cs441 (0.00307413, 0.000532983);
  MassPointParameters mpp441 (cs441, 20849);
  StopNeutralinoMap.insert(std::make_pair(mp441, mpp441));

  MassPoint mp442 (1100, 250);
  StopCrossSection cs442 (0.00307413, 0.000532983);
  MassPointParameters mpp442 (cs442, 21487);
  StopNeutralinoMap.insert(std::make_pair(mp442, mpp442));

  MassPoint mp443 (1100, 300);
  StopCrossSection cs443 (0.00307413, 0.000532983);
  MassPointParameters mpp443 (cs443, 18002);
  StopNeutralinoMap.insert(std::make_pair(mp443, mpp443));

  MassPoint mp444 (1100, 350);
  StopCrossSection cs444 (0.00307413, 0.000532983);
  MassPointParameters mpp444 (cs444, 17150);
  StopNeutralinoMap.insert(std::make_pair(mp444, mpp444));

  MassPoint mp445 (1100, 400);
  StopCrossSection cs445 (0.00307413, 0.000532983);
  MassPointParameters mpp445 (cs445, 21811);
  StopNeutralinoMap.insert(std::make_pair(mp445, mpp445));

  MassPoint mp446 (1100, 450);
  StopCrossSection cs446 (0.00307413, 0.000532983);
  MassPointParameters mpp446 (cs446, 18963);
  StopNeutralinoMap.insert(std::make_pair(mp446, mpp446));

  MassPoint mp447 (1100, 500);
  StopCrossSection cs447 (0.00307413, 0.000532983);
  MassPointParameters mpp447 (cs447, 20169);
  StopNeutralinoMap.insert(std::make_pair(mp447, mpp447));

  MassPoint mp448 (1100, 550);
  StopCrossSection cs448 (0.00307413, 0.000532983);
  MassPointParameters mpp448 (cs448, 15726);
  StopNeutralinoMap.insert(std::make_pair(mp448, mpp448));

  MassPoint mp449 (1100, 600);
  StopCrossSection cs449 (0.00307413, 0.000532983);
  MassPointParameters mpp449 (cs449, 18043);
  StopNeutralinoMap.insert(std::make_pair(mp449, mpp449));

  MassPoint mp450 (1100, 650);
  StopCrossSection cs450 (0.00307413, 0.000532983);
  MassPointParameters mpp450 (cs450, 17519);
  StopNeutralinoMap.insert(std::make_pair(mp450, mpp450));

  MassPoint mp451 (1150, 1);
  StopCrossSection cs451 (0.00221047, 0.000396247);
  MassPointParameters mpp451 (cs451, 20192);
  StopNeutralinoMap.insert(std::make_pair(mp451, mpp451));

  MassPoint mp452 (1150, 50);
  StopCrossSection cs452 (0.00221047, 0.000396247);
  MassPointParameters mpp452 (cs452, 20634);
  StopNeutralinoMap.insert(std::make_pair(mp452, mpp452));

  MassPoint mp453 (1150, 100);
  StopCrossSection cs453 (0.00221047, 0.000396247);
  MassPointParameters mpp453 (cs453, 18944);
  StopNeutralinoMap.insert(std::make_pair(mp453, mpp453));

  MassPoint mp454 (1150, 150);
  StopCrossSection cs454 (0.00221047, 0.000396247);
  MassPointParameters mpp454 (cs454, 18300);
  StopNeutralinoMap.insert(std::make_pair(mp454, mpp454));

  MassPoint mp455 (1150, 200);
  StopCrossSection cs455 (0.00221047, 0.000396247);
  MassPointParameters mpp455 (cs455, 17542);
  StopNeutralinoMap.insert(std::make_pair(mp455, mpp455));

  MassPoint mp456 (1150, 250);
  StopCrossSection cs456 (0.00221047, 0.000396247);
  MassPointParameters mpp456 (cs456, 17251);
  StopNeutralinoMap.insert(std::make_pair(mp456, mpp456));

  MassPoint mp457 (1150, 300);
  StopCrossSection cs457 (0.00221047, 0.000396247);
  MassPointParameters mpp457 (cs457, 19473);
  StopNeutralinoMap.insert(std::make_pair(mp457, mpp457));

  MassPoint mp458 (1150, 350);
  StopCrossSection cs458 (0.00221047, 0.000396247);
  MassPointParameters mpp458 (cs458, 16007);
  StopNeutralinoMap.insert(std::make_pair(mp458, mpp458));

  MassPoint mp459 (1150, 400);
  StopCrossSection cs459 (0.00221047, 0.000396247);
  MassPointParameters mpp459 (cs459, 20796);
  StopNeutralinoMap.insert(std::make_pair(mp459, mpp459));

  MassPoint mp460 (1150, 450);
  StopCrossSection cs460 (0.00221047, 0.000396247);
  MassPointParameters mpp460 (cs460, 18274);
  StopNeutralinoMap.insert(std::make_pair(mp460, mpp460));

  MassPoint mp461 (1150, 500);
  StopCrossSection cs461 (0.00221047, 0.000396247);
  MassPointParameters mpp461 (cs461, 21031);
  StopNeutralinoMap.insert(std::make_pair(mp461, mpp461));

  MassPoint mp462 (1150, 550);
  StopCrossSection cs462 (0.00221047, 0.000396247);
  MassPointParameters mpp462 (cs462, 18437);
  StopNeutralinoMap.insert(std::make_pair(mp462, mpp462));

  MassPoint mp463 (1150, 600);
  StopCrossSection cs463 (0.00221047, 0.000396247);
  MassPointParameters mpp463 (cs463, 18597);
  StopNeutralinoMap.insert(std::make_pair(mp463, mpp463));

  MassPoint mp464 (1150, 650);
  StopCrossSection cs464 (0.00221047, 0.000396247);
  MassPointParameters mpp464 (cs464, 19505);
  StopNeutralinoMap.insert(std::make_pair(mp464, mpp464));

  MassPoint mp465 (1200, 1);
  StopCrossSection cs465 (0.00159844, 0.000296045);
  MassPointParameters mpp465 (cs465, 16938);
  StopNeutralinoMap.insert(std::make_pair(mp465, mpp465));

  MassPoint mp466 (1200, 50);
  StopCrossSection cs466 (0.00159844, 0.000296045);
  MassPointParameters mpp466 (cs466, 17111);
  StopNeutralinoMap.insert(std::make_pair(mp466, mpp466));

  MassPoint mp467 (1200, 100);
  StopCrossSection cs467 (0.00159844, 0.000296045);
  MassPointParameters mpp467 (cs467, 21532);
  StopNeutralinoMap.insert(std::make_pair(mp467, mpp467));

  MassPoint mp468 (1200, 150);
  StopCrossSection cs468 (0.00159844, 0.000296045);
  MassPointParameters mpp468 (cs468, 20130);
  StopNeutralinoMap.insert(std::make_pair(mp468, mpp468));

  MassPoint mp469 (1200, 200);
  StopCrossSection cs469 (0.00159844, 0.000296045);
  MassPointParameters mpp469 (cs469, 21315);
  StopNeutralinoMap.insert(std::make_pair(mp469, mpp469));

  MassPoint mp470 (1200, 250);
  StopCrossSection cs470 (0.00159844, 0.000296045);
  MassPointParameters mpp470 (cs470, 18727);
  StopNeutralinoMap.insert(std::make_pair(mp470, mpp470));

  MassPoint mp471 (1200, 300);
  StopCrossSection cs471 (0.00159844, 0.000296045);
  MassPointParameters mpp471 (cs471, 18112);
  StopNeutralinoMap.insert(std::make_pair(mp471, mpp471));

  MassPoint mp472 (1200, 350);
  StopCrossSection cs472 (0.00159844, 0.000296045);
  MassPointParameters mpp472 (cs472, 18028);
  StopNeutralinoMap.insert(std::make_pair(mp472, mpp472));

  MassPoint mp473 (1200, 400);
  StopCrossSection cs473 (0.00159844, 0.000296045);
  MassPointParameters mpp473 (cs473, 20219);
  StopNeutralinoMap.insert(std::make_pair(mp473, mpp473));

  MassPoint mp474 (1200, 450);
  StopCrossSection cs474 (0.00159844, 0.000296045);
  MassPointParameters mpp474 (cs474, 23181);
  StopNeutralinoMap.insert(std::make_pair(mp474, mpp474));

  MassPoint mp475 (1200, 500);
  StopCrossSection cs475 (0.00159844, 0.000296045);
  MassPointParameters mpp475 (cs475, 18878);
  StopNeutralinoMap.insert(std::make_pair(mp475, mpp475));

  MassPoint mp476 (1200, 550);
  StopCrossSection cs476 (0.00159844, 0.000296045);
  MassPointParameters mpp476 (cs476, 21553);
  StopNeutralinoMap.insert(std::make_pair(mp476, mpp476));

  MassPoint mp477 (1200, 600);
  StopCrossSection cs477 (0.00159844, 0.000296045);
  MassPointParameters mpp477 (cs477, 19744);
  StopNeutralinoMap.insert(std::make_pair(mp477, mpp477));

  MassPoint mp478 (1200, 650);
  StopCrossSection cs478 (0.00159844, 0.000296045);
  MassPointParameters mpp478 (cs478, 21272);
  StopNeutralinoMap.insert(std::make_pair(mp478, mpp478));
  
  MassPoint mpisr0 (258, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 1.09939));

  MassPoint mpisr1 (267, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 1.09652));

  MassPoint mpisr2 (275, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 1.10285));

  MassPoint mpisr3 (275, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 1.10219));

  MassPoint mpisr4 (275, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 1.10108));

  MassPoint mpisr5 (275, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 1.09992));

  MassPoint mpisr6 (275, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 1.09982));

  MassPoint mpisr7 (275, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 1.09519));

  MassPoint mpisr8 (275, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 1.09265));

  MassPoint mpisr9 (275, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 1.09216));

  MassPoint mpisr10 (275, 188);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 1.09456));

  MassPoint mpisr11 (283, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 1.10093));

  MassPoint mpisr12 (292, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 1.10276));

  MassPoint mpisr13 (300, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 1.10687));

  MassPoint mpisr14 (300, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 1.10644));

  MassPoint mpisr15 (300, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 1.10644));

  MassPoint mpisr16 (300, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 1.10535));

  MassPoint mpisr17 (300, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 1.10603));

  MassPoint mpisr18 (300, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 1.1038));

  MassPoint mpisr19 (300, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 1.09642));

  MassPoint mpisr20 (300, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 1.09563));

  MassPoint mpisr21 (300, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 1.09772));

  MassPoint mpisr22 (300, 213);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 1.0967));

  MassPoint mpisr23 (308, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 1.1043));

  MassPoint mpisr24 (317, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 1.1044));

  MassPoint mpisr25 (325, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 1.1088));

  MassPoint mpisr26 (325, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 1.11033));

  MassPoint mpisr27 (325, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 1.11037));

  MassPoint mpisr28 (325, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 1.10851));

  MassPoint mpisr29 (325, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 1.10818));

  MassPoint mpisr30 (325, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 1.10746));

  MassPoint mpisr31 (325, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 1.10279));

  MassPoint mpisr32 (325, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 1.10115));

  MassPoint mpisr33 (325, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 1.10065));

  MassPoint mpisr34 (325, 238);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 1.10038));

  MassPoint mpisr35 (333, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 1.10925));

  MassPoint mpisr36 (342, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 1.10865));

  MassPoint mpisr37 (350, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 1.1139));

  MassPoint mpisr38 (350, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 1.11552));

  MassPoint mpisr39 (350, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 1.11271));

  MassPoint mpisr40 (350, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 1.11288));

  MassPoint mpisr41 (350, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 1.11276));

  MassPoint mpisr42 (350, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 1.11163));

  MassPoint mpisr43 (350, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 1.11045));

  MassPoint mpisr44 (350, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 1.10544));

  MassPoint mpisr45 (350, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 1.10282));

  MassPoint mpisr46 (350, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 1.10173));

  MassPoint mpisr47 (350, 263);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 1.10243));

  MassPoint mpisr48 (358, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 1.1118));

  MassPoint mpisr49 (367, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 1.11048));

  MassPoint mpisr50 (375, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 1.11803));

  MassPoint mpisr51 (375, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 1.11762));

  MassPoint mpisr52 (375, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 1.11712));

  MassPoint mpisr53 (375, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 1.11556));

  MassPoint mpisr54 (375, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 1.11525));

  MassPoint mpisr55 (375, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 1.11384));

  MassPoint mpisr56 (375, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 1.10938));

  MassPoint mpisr57 (375, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 1.10722));

  MassPoint mpisr58 (375, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 1.10616));

  MassPoint mpisr59 (375, 288);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 1.10725));

  MassPoint mpisr60 (383, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 1.11511));

  MassPoint mpisr61 (392, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 1.11102));

  MassPoint mpisr62 (400, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 1.12242));

  MassPoint mpisr63 (400, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 1.12138));

  MassPoint mpisr64 (400, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 1.12198));

  MassPoint mpisr65 (400, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 1.12093));

  MassPoint mpisr66 (400, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 1.11733));

  MassPoint mpisr67 (400, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 1.11907));

  MassPoint mpisr68 (400, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 1.11708));

  MassPoint mpisr69 (400, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 1.11843));

  MassPoint mpisr70 (400, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 1.11194));

  MassPoint mpisr71 (400, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 1.11043));

  MassPoint mpisr72 (400, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 1.10846));

  MassPoint mpisr73 (400, 313);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 1.10905));

  MassPoint mpisr74 (408, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 1.11916));

  MassPoint mpisr75 (417, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 1.1141));

  MassPoint mpisr76 (425, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 1.12245));

  MassPoint mpisr77 (425, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 1.12223));

  MassPoint mpisr78 (425, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 1.12179));

  MassPoint mpisr79 (425, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 1.1212));

  MassPoint mpisr80 (425, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 1.12004));

  MassPoint mpisr81 (425, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 1.11788));

  MassPoint mpisr82 (425, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 1.11047));

  MassPoint mpisr83 (425, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 1.11136));

  MassPoint mpisr84 (425, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 1.11147));

  MassPoint mpisr85 (425, 338);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 1.10975));

  MassPoint mpisr86 (433, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 1.11947));

  MassPoint mpisr87 (442, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 1.1165));

  MassPoint mpisr88 (450, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 1.12939));

  MassPoint mpisr89 (450, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 1.1298));

  MassPoint mpisr90 (450, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 1.1246));

  MassPoint mpisr91 (450, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 1.12658));

  MassPoint mpisr92 (450, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 1.12292));

  MassPoint mpisr93 (450, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 1.12345));

  MassPoint mpisr94 (450, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 1.12255));

  MassPoint mpisr95 (450, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 1.12208));

  MassPoint mpisr96 (450, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 1.11882));

  MassPoint mpisr97 (450, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 1.11625));

  MassPoint mpisr98 (450, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 1.11436));

  MassPoint mpisr99 (450, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 1.11378));

  MassPoint mpisr100 (450, 363);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 1.113));

  MassPoint mpisr101 (458, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 1.12095));

  MassPoint mpisr102 (467, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 1.11719));

  MassPoint mpisr103 (475, 175);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 1.12955));

  MassPoint mpisr104 (475, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 1.12545));

  MassPoint mpisr105 (475, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 1.13005));

  MassPoint mpisr106 (475, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 1.12246));

  MassPoint mpisr107 (475, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 1.12453));

  MassPoint mpisr108 (475, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 1.12444));

  MassPoint mpisr109 (475, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 1.11947));

  MassPoint mpisr110 (475, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 1.11748));

  MassPoint mpisr111 (475, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 1.11656));

  MassPoint mpisr112 (475, 388);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 1.11581));

  MassPoint mpisr113 (483, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 1.12339));

  MassPoint mpisr114 (492, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 1.12258));

  MassPoint mpisr115 (500, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 1.13682));

  MassPoint mpisr116 (500, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 1.13149));

  MassPoint mpisr117 (500, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 1.13353));

  MassPoint mpisr118 (500, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 1.13513));

  MassPoint mpisr119 (500, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 1.13079));

  MassPoint mpisr120 (500, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 1.12813));

  MassPoint mpisr121 (500, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 1.12317));

  MassPoint mpisr122 (500, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 1.12434));

  MassPoint mpisr123 (500, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 1.12236));

  MassPoint mpisr124 (500, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 1.12149));

  MassPoint mpisr125 (500, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 1.11871));

  MassPoint mpisr126 (500, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 1.12142));

  MassPoint mpisr127 (500, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 1.11723));

  MassPoint mpisr128 (500, 413);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 1.11538));

  MassPoint mpisr129 (508, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 1.12803));

  MassPoint mpisr130 (517, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 1.11866));

  MassPoint mpisr131 (525, 225);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 1.13105));

  MassPoint mpisr132 (525, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 1.12959));

  MassPoint mpisr133 (525, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 1.1259));

  MassPoint mpisr134 (525, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 1.1258));

  MassPoint mpisr135 (525, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 1.12817));

  MassPoint mpisr136 (525, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 1.12323));

  MassPoint mpisr137 (525, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 1.12303));

  MassPoint mpisr138 (525, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 1.12073));

  MassPoint mpisr139 (525, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 1.11615));

  MassPoint mpisr140 (525, 438);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 1.12139));

  MassPoint mpisr141 (533, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 1.12904));

  MassPoint mpisr142 (542, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 1.12538));

  MassPoint mpisr143 (550, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 1.14376));

  MassPoint mpisr144 (550, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 1.14437));

  MassPoint mpisr145 (550, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 1.14353));

  MassPoint mpisr146 (550, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 1.13843));

  MassPoint mpisr147 (550, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 1.13628));

  MassPoint mpisr148 (550, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 1.14014));

  MassPoint mpisr149 (550, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 1.13282));

  MassPoint mpisr150 (550, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 1.13635));

  MassPoint mpisr151 (550, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 1.12987));

  MassPoint mpisr152 (550, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 1.13007));

  MassPoint mpisr153 (550, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 1.12971));

  MassPoint mpisr154 (550, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 1.12552));

  MassPoint mpisr155 (550, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 1.12439));

  MassPoint mpisr156 (550, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 1.12283));

  MassPoint mpisr157 (550, 463);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 1.12003));

  MassPoint mpisr158 (558, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 1.12618));

  MassPoint mpisr159 (567, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 1.1325));

  MassPoint mpisr160 (575, 275);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 1.13582));

  MassPoint mpisr161 (575, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 1.1326));

  MassPoint mpisr162 (575, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 1.13672));

  MassPoint mpisr163 (575, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 1.13294));

  MassPoint mpisr164 (575, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 1.13459));

  MassPoint mpisr165 (575, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 1.13107));

  MassPoint mpisr166 (575, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 1.12545));

  MassPoint mpisr167 (575, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 1.12595));

  MassPoint mpisr168 (575, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 1.12067));

  MassPoint mpisr169 (575, 488);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 1.12097));

  MassPoint mpisr170 (583, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 1.13288));

  MassPoint mpisr171 (592, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 1.1246));

  MassPoint mpisr172 (600, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 1.14754));

  MassPoint mpisr173 (600, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 1.15362));

  MassPoint mpisr174 (600, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 1.15038));

  MassPoint mpisr175 (600, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 1.14482));

  MassPoint mpisr176 (600, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 1.14897));

  MassPoint mpisr177 (600, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 1.13337));

  MassPoint mpisr178 (600, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 1.14057));

  MassPoint mpisr179 (600, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 1.13166));

  MassPoint mpisr180 (600, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 1.1415));

  MassPoint mpisr181 (600, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 1.1313));

  MassPoint mpisr182 (600, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 1.13188));

  MassPoint mpisr183 (600, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 1.13149));

  MassPoint mpisr184 (600, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr184, 1.12344));

  MassPoint mpisr185 (600, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr185, 1.12318));

  MassPoint mpisr186 (600, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr186, 1.114));

  MassPoint mpisr187 (600, 513);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr187, 1.12465));

  MassPoint mpisr188 (608, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr188, 1.13318));

  MassPoint mpisr189 (617, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr189, 1.13585));

  MassPoint mpisr190 (625, 325);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr190, 1.13639));

  MassPoint mpisr191 (625, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr191, 1.12981));

  MassPoint mpisr192 (625, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr192, 1.13574));

  MassPoint mpisr193 (625, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr193, 1.13988));

  MassPoint mpisr194 (625, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr194, 1.12688));

  MassPoint mpisr195 (625, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr195, 1.13626));

  MassPoint mpisr196 (625, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr196, 1.12071));

  MassPoint mpisr197 (625, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr197, 1.12731));

  MassPoint mpisr198 (625, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr198, 1.12857));

  MassPoint mpisr199 (625, 538);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr199, 1.12557));

  MassPoint mpisr200 (633, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr200, 1.13177));

  MassPoint mpisr201 (642, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr201, 1.13209));

  MassPoint mpisr202 (650, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr202, 1.1482));

  MassPoint mpisr203 (650, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr203, 1.15763));

  MassPoint mpisr204 (650, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr204, 1.15855));

  MassPoint mpisr205 (650, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr205, 1.14548));

  MassPoint mpisr206 (650, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr206, 1.15617));

  MassPoint mpisr207 (650, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr207, 1.14977));

  MassPoint mpisr208 (650, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr208, 1.13456));

  MassPoint mpisr209 (650, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr209, 1.1338));

  MassPoint mpisr210 (650, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr210, 1.12472));

  MassPoint mpisr211 (650, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr211, 1.13228));

  MassPoint mpisr212 (650, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr212, 1.14667));

  MassPoint mpisr213 (650, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr213, 1.13262));

  MassPoint mpisr214 (650, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr214, 1.14043));

  MassPoint mpisr215 (650, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr215, 1.12905));

  MassPoint mpisr216 (650, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr216, 1.13324));

  MassPoint mpisr217 (650, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr217, 1.11932));

  MassPoint mpisr218 (650, 563);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr218, 1.12924));

  MassPoint mpisr219 (658, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr219, 1.12514));

  MassPoint mpisr220 (667, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr220, 1.12435));

  MassPoint mpisr221 (675, 375);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr221, 1.14882));

  MassPoint mpisr222 (675, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr222, 1.14118));

  MassPoint mpisr223 (675, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr223, 1.13847));

  MassPoint mpisr224 (675, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr224, 1.13794));

  MassPoint mpisr225 (675, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr225, 1.14758));

  MassPoint mpisr226 (675, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr226, 1.14628));

  MassPoint mpisr227 (675, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr227, 1.13982));

  MassPoint mpisr228 (675, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr228, 1.13435));

  MassPoint mpisr229 (675, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr229, 1.12391));

  MassPoint mpisr230 (675, 588);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr230, 1.12392));

  MassPoint mpisr231 (683, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr231, 1.13949));

  MassPoint mpisr232 (692, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr232, 1.13009));

  MassPoint mpisr233 (700, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr233, 1.15994));

  MassPoint mpisr234 (700, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr234, 1.16807));

  MassPoint mpisr235 (700, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr235, 1.15958));

  MassPoint mpisr236 (700, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr236, 1.16041));

  MassPoint mpisr237 (700, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr237, 1.15311));

  MassPoint mpisr238 (700, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr238, 1.155));

  MassPoint mpisr239 (700, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr239, 1.1535));

  MassPoint mpisr240 (700, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr240, 1.15872));

  MassPoint mpisr241 (700, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr241, 1.14083));

  MassPoint mpisr242 (700, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr242, 1.14703));

  MassPoint mpisr243 (700, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr243, 1.14662));

  MassPoint mpisr244 (700, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr244, 1.14751));

  MassPoint mpisr245 (700, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr245, 1.14056));

  MassPoint mpisr246 (700, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr246, 1.13689));

  MassPoint mpisr247 (700, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr247, 1.14499));

  MassPoint mpisr248 (700, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr248, 1.12886));

  MassPoint mpisr249 (700, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr249, 1.12144));

  MassPoint mpisr250 (700, 613);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr250, 1.12191));

  MassPoint mpisr251 (708, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr251, 1.13039));

  MassPoint mpisr252 (717, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr252, 1.14802));

  MassPoint mpisr253 (725, 425);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr253, 1.14836));

  MassPoint mpisr254 (725, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr254, 1.14807));

  MassPoint mpisr255 (725, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr255, 1.13741));

  MassPoint mpisr256 (725, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr256, 1.13104));

  MassPoint mpisr257 (725, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr257, 1.14591));

  MassPoint mpisr258 (725, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr258, 1.13439));

  MassPoint mpisr259 (725, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr259, 1.11445));

  MassPoint mpisr260 (725, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr260, 1.12962));

  MassPoint mpisr261 (725, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr261, 1.13948));

  MassPoint mpisr262 (725, 638);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr262, 1.12652));

  MassPoint mpisr263 (733, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr263, 1.14281));

  MassPoint mpisr264 (742, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr264, 1.12871));

  MassPoint mpisr265 (750, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr265, 1.16636));

  MassPoint mpisr266 (750, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr266, 1.15857));

  MassPoint mpisr267 (750, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr267, 1.16633));

  MassPoint mpisr268 (750, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr268, 1.16172));

  MassPoint mpisr269 (750, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr269, 1.1725));

  MassPoint mpisr270 (750, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr270, 1.1429));

  MassPoint mpisr271 (750, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr271, 1.14705));

  MassPoint mpisr272 (750, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr272, 1.15204));

  MassPoint mpisr273 (750, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr273, 1.12263));

  MassPoint mpisr274 (750, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr274, 1.15549));

  MassPoint mpisr275 (750, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr275, 1.1495));

  MassPoint mpisr276 (750, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr276, 1.15792));

  MassPoint mpisr277 (750, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr277, 1.14426));

  MassPoint mpisr278 (750, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr278, 1.135));

  MassPoint mpisr279 (750, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr279, 1.11677));

  MassPoint mpisr280 (750, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr280, 1.14233));

  MassPoint mpisr281 (750, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr281, 1.1338));

  MassPoint mpisr282 (750, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr282, 1.12112));

  MassPoint mpisr283 (758, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr283, 1.13855));

  MassPoint mpisr284 (767, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr284, 1.14584));

  MassPoint mpisr285 (775, 475);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr285, 1.15004));

  MassPoint mpisr286 (775, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr286, 1.13012));

  MassPoint mpisr287 (775, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr287, 1.13234));

  MassPoint mpisr288 (775, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr288, 1.14058));

  MassPoint mpisr289 (775, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr289, 1.1484));

  MassPoint mpisr290 (775, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr290, 1.14797));

  MassPoint mpisr291 (775, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr291, 1.12148));

  MassPoint mpisr292 (775, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr292, 1.13118));

  MassPoint mpisr293 (783, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr293, 1.13886));

  MassPoint mpisr294 (792, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr294, 1.12784));

  MassPoint mpisr295 (800, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr295, 1.16574));

  MassPoint mpisr296 (800, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr296, 1.16605));

  MassPoint mpisr297 (800, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr297, 1.16356));

  MassPoint mpisr298 (800, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr298, 1.15644));

  MassPoint mpisr299 (800, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr299, 1.16301));

  MassPoint mpisr300 (800, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr300, 1.16884));

  MassPoint mpisr301 (800, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr301, 1.16093));

  MassPoint mpisr302 (800, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr302, 1.16127));

  MassPoint mpisr303 (800, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr303, 1.15514));

  MassPoint mpisr304 (800, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr304, 1.15149));

  MassPoint mpisr305 (800, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr305, 1.14955));

  MassPoint mpisr306 (800, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr306, 1.1405));

  MassPoint mpisr307 (800, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr307, 1.14245));

  MassPoint mpisr308 (800, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr308, 1.14363));

  MassPoint mpisr309 (800, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr309, 1.13627));

  MassPoint mpisr310 (800, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr310, 1.13621));

  MassPoint mpisr311 (800, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr311, 1.1292));

  MassPoint mpisr312 (808, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr312, 1.14083));

  MassPoint mpisr313 (817, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr313, 1.13463));

  MassPoint mpisr314 (825, 525);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr314, 1.15137));

  MassPoint mpisr315 (825, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr315, 1.14363));

  MassPoint mpisr316 (825, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr316, 1.13486));

  MassPoint mpisr317 (825, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr317, 1.13348));

  MassPoint mpisr318 (825, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr318, 1.14001));

  MassPoint mpisr319 (825, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr319, 1.1343));

  MassPoint mpisr320 (833, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr320, 1.1546));

  MassPoint mpisr321 (850, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr321, 1.18319));

  MassPoint mpisr322 (850, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr322, 1.19036));

  MassPoint mpisr323 (850, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr323, 1.16609));

  MassPoint mpisr324 (850, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr324, 1.17082));

  MassPoint mpisr325 (850, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr325, 1.1754));

  MassPoint mpisr326 (850, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr326, 1.15958));

  MassPoint mpisr327 (850, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr327, 1.15632));

  MassPoint mpisr328 (850, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr328, 1.15144));

  MassPoint mpisr329 (850, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr329, 1.15964));

  MassPoint mpisr330 (850, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr330, 1.15231));

  MassPoint mpisr331 (850, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr331, 1.1398));

  MassPoint mpisr332 (850, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr332, 1.12633));

  MassPoint mpisr333 (850, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr333, 1.13893));

  MassPoint mpisr334 (850, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr334, 1.13216));

  MassPoint mpisr335 (850, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr335, 1.14319));

  MassPoint mpisr336 (850, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr336, 1.13978));

  MassPoint mpisr337 (875, 575);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr337, 1.1396));

  MassPoint mpisr338 (875, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr338, 1.13458));

  MassPoint mpisr339 (875, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr339, 1.1422));

  MassPoint mpisr340 (875, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr340, 1.14109));

  MassPoint mpisr341 (900, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr341, 1.17442));

  MassPoint mpisr342 (900, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr342, 1.17724));

  MassPoint mpisr343 (900, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr343, 1.17204));

  MassPoint mpisr344 (900, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr344, 1.17348));

  MassPoint mpisr345 (900, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr345, 1.1739));

  MassPoint mpisr346 (900, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr346, 1.15122));

  MassPoint mpisr347 (900, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr347, 1.16623));

  MassPoint mpisr348 (900, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr348, 1.16449));

  MassPoint mpisr349 (900, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr349, 1.1783));

  MassPoint mpisr350 (900, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr350, 1.1727));

  MassPoint mpisr351 (900, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr351, 1.15646));

  MassPoint mpisr352 (900, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr352, 1.16534));

  MassPoint mpisr353 (900, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr353, 1.14733));

  MassPoint mpisr354 (900, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr354, 1.15658));

  MassPoint mpisr355 (900, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr355, 1.15216));

  MassPoint mpisr356 (925, 625);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr356, 1.14785));

  MassPoint mpisr357 (925, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr357, 1.13313));

  MassPoint mpisr358 (950, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr358, 1.16598));

  MassPoint mpisr359 (950, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr359, 1.17097));

  MassPoint mpisr360 (950, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr360, 1.17542));

  MassPoint mpisr361 (950, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr361, 1.17886));

  MassPoint mpisr362 (950, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr362, 1.16387));

  MassPoint mpisr363 (950, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr363, 1.17795));

  MassPoint mpisr364 (950, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr364, 1.16374));

  MassPoint mpisr365 (950, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr365, 1.17061));

  MassPoint mpisr366 (950, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr366, 1.17313));

  MassPoint mpisr367 (950, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr367, 1.15202));

  MassPoint mpisr368 (950, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr368, 1.16556));

  MassPoint mpisr369 (950, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr369, 1.15359));

  MassPoint mpisr370 (950, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr370, 1.15367));

  MassPoint mpisr371 (950, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr371, 1.15128));

  MassPoint mpisr372 (1000, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr372, 1.16677));

  MassPoint mpisr373 (1000, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr373, 1.17969));

  MassPoint mpisr374 (1000, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr374, 1.19146));

  MassPoint mpisr375 (1000, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr375, 1.18331));

  MassPoint mpisr376 (1000, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr376, 1.19863));

  MassPoint mpisr377 (1000, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr377, 1.15597));

  MassPoint mpisr378 (1000, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr378, 1.15939));

  MassPoint mpisr379 (1000, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr379, 1.18102));

  MassPoint mpisr380 (1000, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr380, 1.17657));

  MassPoint mpisr381 (1000, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr381, 1.17589));

  MassPoint mpisr382 (1000, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr382, 1.15259));

  MassPoint mpisr383 (1000, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr383, 1.17282));

  MassPoint mpisr384 (1000, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr384, 1.15516));

  MassPoint mpisr385 (1000, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr385, 1.15701));

  MassPoint mpisr386 (1050, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr386, 1.19495));

  MassPoint mpisr387 (1050, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr387, 1.19304));

  MassPoint mpisr388 (1050, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr388, 1.18343));

  MassPoint mpisr389 (1050, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr389, 1.18851));

  MassPoint mpisr390 (1050, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr390, 1.19881));

  MassPoint mpisr391 (1050, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr391, 1.18418));

  MassPoint mpisr392 (1050, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr392, 1.19096));

  MassPoint mpisr393 (1050, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr393, 1.20406));

  MassPoint mpisr394 (1050, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr394, 1.17712));

  MassPoint mpisr395 (1050, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr395, 1.1965));

  MassPoint mpisr396 (1050, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr396, 1.16857));

  MassPoint mpisr397 (1050, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr397, 1.16531));

  MassPoint mpisr398 (1050, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr398, 1.16034));

  MassPoint mpisr399 (1050, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr399, 1.15379));

  MassPoint mpisr400 (1100, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr400, 1.1881));

  MassPoint mpisr401 (1100, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr401, 1.2007));

  MassPoint mpisr402 (1100, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr402, 1.17227));

  MassPoint mpisr403 (1100, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr403, 1.18213));

  MassPoint mpisr404 (1100, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr404, 1.18482));

  MassPoint mpisr405 (1100, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr405, 1.20425));

  MassPoint mpisr406 (1100, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr406, 1.20723));

  MassPoint mpisr407 (1100, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr407, 1.18643));

  MassPoint mpisr408 (1100, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr408, 1.20526));

  MassPoint mpisr409 (1100, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr409, 1.19195));

  MassPoint mpisr410 (1100, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr410, 1.16528));

  MassPoint mpisr411 (1100, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr411, 1.17079));

  MassPoint mpisr412 (1100, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr412, 1.16766));

  MassPoint mpisr413 (1100, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr413, 1.15882));

  MassPoint mpisr414 (1150, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr414, 1.20455));

  MassPoint mpisr415 (1150, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr415, 1.17532));

  MassPoint mpisr416 (1150, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr416, 1.20228));

  MassPoint mpisr417 (1150, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr417, 1.19624));

  MassPoint mpisr418 (1150, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr418, 1.1929));

  MassPoint mpisr419 (1150, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr419, 1.20319));

  MassPoint mpisr420 (1150, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr420, 1.1762));

  MassPoint mpisr421 (1150, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr421, 1.18633));

  MassPoint mpisr422 (1150, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr422, 1.18039));

  MassPoint mpisr423 (1150, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr423, 1.19445));

  MassPoint mpisr424 (1150, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr424, 1.17712));

  MassPoint mpisr425 (1150, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr425, 1.17515));

  MassPoint mpisr426 (1150, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr426, 1.17276));

  MassPoint mpisr427 (1150, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr427, 1.15744));

  MassPoint mpisr428 (1200, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr428, 1.19212));

  MassPoint mpisr429 (1200, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr429, 1.1862));

  MassPoint mpisr430 (1200, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr430, 1.20723));

  MassPoint mpisr431 (1200, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr431, 1.18632));

  MassPoint mpisr432 (1200, 200);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr432, 1.19046));

  MassPoint mpisr433 (1200, 250);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr433, 1.20552));

  MassPoint mpisr434 (1200, 300);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr434, 1.19503));

  MassPoint mpisr435 (1200, 350);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr435, 1.20144));

  MassPoint mpisr436 (1200, 400);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr436, 1.16597));

  MassPoint mpisr437 (1200, 450);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr437, 1.17837));

  MassPoint mpisr438 (1200, 500);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr438, 1.17703));

  MassPoint mpisr439 (1200, 550);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr439, 1.19067));

  MassPoint mpisr440 (1200, 600);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr440, 1.18409));

  MassPoint mpisr441 (1200, 650);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr441, 1.16849));

  MassPoint mpisr442 (150, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 1.06748));

  MassPoint mpisr443 (150, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 1.06533));

  MassPoint mpisr444 (150, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 1.06588));

  MassPoint mpisr445 (150, 63);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 1.06486));

  MassPoint mpisr446 (167, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 1.07194));

  MassPoint mpisr447 (175, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 1.07908));

  MassPoint mpisr448 (175, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 1.07409));

  MassPoint mpisr449 (175, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 1.07029));

  MassPoint mpisr450 (175, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 1.07287));

  MassPoint mpisr451 (175, 88);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 1.07242));

  MassPoint mpisr452 (183, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 1.08198));

  MassPoint mpisr453 (192, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 1.07791));

  MassPoint mpisr454 (200, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 1.08562));

  MassPoint mpisr455 (200, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 1.08602));

  MassPoint mpisr456 (200, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 1.0791));

  MassPoint mpisr457 (200, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 1.07631));

  MassPoint mpisr458 (200, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 1.07787));

  MassPoint mpisr459 (200, 113);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 1.07844));

  MassPoint mpisr460 (208, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 1.08555));

  MassPoint mpisr461 (217, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 1.08409));

  MassPoint mpisr462 (225, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 1.09003));

  MassPoint mpisr463 (225, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 1.09003));

  MassPoint mpisr464 (225, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 1.08924));

  MassPoint mpisr465 (225, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 1.08491));

  MassPoint mpisr466 (225, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 1.08162));

  MassPoint mpisr467 (225, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 1.08207));

  MassPoint mpisr468 (225, 138);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 1.08303));

  MassPoint mpisr469 (233, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 1.09247));

  MassPoint mpisr470 (242, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 1.08935));

  MassPoint mpisr471 (250, 1);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 1.09542));

  MassPoint mpisr472 (250, 25);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 1.09653));

  MassPoint mpisr473 (250, 50);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 1.09457));

  MassPoint mpisr474 (250, 75);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 1.09351));

  MassPoint mpisr475 (250, 100);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 1.08955));

  MassPoint mpisr476 (250, 125);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 1.08645));

  MassPoint mpisr477 (250, 150);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 1.08635));

  MassPoint mpisr478 (250, 163);
  StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 1.08735));

  } else if (SUSYProductionProcess=="T2bW") {

    MassPoint mp0 (200, 1);
  StopCrossSection cs0 (64.5085, 9.29555);
  MassPointParameters mpp0 (cs0, 1303484);
  StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

  MassPoint mp1 (200, 25);
  StopCrossSection cs1 (64.5085, 9.29555);
  MassPointParameters mpp1 (cs1, 1288848);
  StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

  MassPoint mp2 (250, 1);
  StopCrossSection cs2 (21.5949, 3.03613);
  MassPointParameters mpp2 (cs2, 436238);
  StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

  MassPoint mp3 (250, 50);
  StopCrossSection cs3 (21.5949, 3.03613);
  MassPointParameters mpp3 (cs3, 449820);
  StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

  MassPoint mp4 (250, 75);
  StopCrossSection cs4 (21.5949, 3.03613);
  MassPointParameters mpp4 (cs4, 427512);
  StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

  MassPoint mp5 (300, 1);
  StopCrossSection cs5 (8.51615, 1.18564);
  MassPointParameters mpp5 (cs5, 484233);
  StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

  MassPoint mp6 (300, 50);
  StopCrossSection cs6 (8.51615, 1.18564);
  MassPointParameters mpp6 (cs6, 494887);
  StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

  MassPoint mp7 (300, 100);
  StopCrossSection cs7 (8.51615, 1.18564);
  MassPointParameters mpp7 (cs7, 499325);
  StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

  MassPoint mp8 (300, 125);
  StopCrossSection cs8 (8.51615, 1.18564);
  MassPointParameters mpp8 (cs8, 506831);
  StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

  MassPoint mp9 (350, 1);
  StopCrossSection cs9 (3.78661, 0.5183);
  MassPointParameters mpp9 (cs9, 451790);
  StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

  MassPoint mp10 (350, 50);
  StopCrossSection cs10 (3.78661, 0.5183);
  MassPointParameters mpp10 (cs10, 466618);
  StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

  MassPoint mp11 (350, 100);
  StopCrossSection cs11 (3.78661, 0.5183);
  MassPointParameters mpp11 (cs11, 466000);
  StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

  MassPoint mp12 (350, 150);
  StopCrossSection cs12 (3.78661, 0.5183);
  MassPointParameters mpp12 (cs12, 454791);
  StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

  MassPoint mp13 (350, 175);
  StopCrossSection cs13 (3.78661, 0.5183);
  MassPointParameters mpp13 (cs13, 472065);
  StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

  MassPoint mp14 (400, 1);
  StopCrossSection cs14 (1.83537, 0.251418);
  MassPointParameters mpp14 (cs14, 473145);
  StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

  MassPoint mp15 (400, 25);
  StopCrossSection cs15 (1.83537, 0.251418);
  MassPointParameters mpp15 (cs15, 490867);
  StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

  MassPoint mp16 (400, 50);
  StopCrossSection cs16 (1.83537, 0.251418);
  MassPointParameters mpp16 (cs16, 471815);
  StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

  MassPoint mp17 (400, 100);
  StopCrossSection cs17 (1.83537, 0.251418);
  MassPointParameters mpp17 (cs17, 490933);
  StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

  MassPoint mp18 (400, 150);
  StopCrossSection cs18 (1.83537, 0.251418);
  MassPointParameters mpp18 (cs18, 488429);
  StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

  MassPoint mp19 (400, 175);
  StopCrossSection cs19 (1.83537, 0.251418);
  MassPointParameters mpp19 (cs19, 479593);
  StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

  MassPoint mp20 (400, 200);
  StopCrossSection cs20 (1.83537, 0.251418);
  MassPointParameters mpp20 (cs20, 488018);
  StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

  MassPoint mp21 (400, 225);
  StopCrossSection cs21 (1.83537, 0.251418);
  MassPointParameters mpp21 (cs21, 485918);
  StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

  MassPoint mp22 (425, 150);
  StopCrossSection cs22 (1.31169, 0.177095);
  MassPointParameters mpp22 (cs22, 483238);
  StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

  MassPoint mp23 (425, 175);
  StopCrossSection cs23 (1.31169, 0.177095);
  MassPointParameters mpp23 (cs23, 481465);
  StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

  MassPoint mp24 (425, 200);
  StopCrossSection cs24 (1.31169, 0.177095);
  MassPointParameters mpp24 (cs24, 488517);
  StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

  MassPoint mp25 (425, 225);
  StopCrossSection cs25 (1.31169, 0.177095);
  MassPointParameters mpp25 (cs25, 489062);
  StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

  MassPoint mp26 (425, 250);
  StopCrossSection cs26 (1.31169, 0.177095);
  MassPointParameters mpp26 (cs26, 477411);
  StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

  MassPoint mp27 (450, 1);
  StopCrossSection cs27 (0.948333, 0.127607);
  MassPointParameters mpp27 (cs27, 356218);
  StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

  MassPoint mp28 (450, 25);
  StopCrossSection cs28 (0.948333, 0.127607);
  MassPointParameters mpp28 (cs28, 362006);
  StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

  MassPoint mp29 (450, 50);
  StopCrossSection cs29 (0.948333, 0.127607);
  MassPointParameters mpp29 (cs29, 355832);
  StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

  MassPoint mp30 (450, 100);
  StopCrossSection cs30 (0.948333, 0.127607);
  MassPointParameters mpp30 (cs30, 358993);
  StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

  MassPoint mp31 (450, 150);
  StopCrossSection cs31 (0.948333, 0.127607);
  MassPointParameters mpp31 (cs31, 340703);
  StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

  MassPoint mp32 (450, 175);
  StopCrossSection cs32 (0.948333, 0.127607);
  MassPointParameters mpp32 (cs32, 350016);
  StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

  MassPoint mp33 (450, 200);
  StopCrossSection cs33 (0.948333, 0.127607);
  MassPointParameters mpp33 (cs33, 368308);
  StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

  MassPoint mp34 (450, 225);
  StopCrossSection cs34 (0.948333, 0.127607);
  MassPointParameters mpp34 (cs34, 355880);
  StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

  MassPoint mp35 (450, 250);
  StopCrossSection cs35 (0.948333, 0.127607);
  MassPointParameters mpp35 (cs35, 365649);
  StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

  MassPoint mp36 (450, 275);
  StopCrossSection cs36 (0.948333, 0.127607);
  MassPointParameters mpp36 (cs36, 358347);
  StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

  MassPoint mp37 (475, 175);
  StopCrossSection cs37 (0.697075, 0.0933565);
  MassPointParameters mpp37 (cs37, 247325);
  StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

  MassPoint mp38 (475, 200);
  StopCrossSection cs38 (0.697075, 0.0933565);
  MassPointParameters mpp38 (cs38, 262692);
  StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

  MassPoint mp39 (475, 225);
  StopCrossSection cs39 (0.697075, 0.0933565);
  MassPointParameters mpp39 (cs39, 262318);
  StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

  MassPoint mp40 (475, 250);
  StopCrossSection cs40 (0.697075, 0.0933565);
  MassPointParameters mpp40 (cs40, 255054);
  StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

  MassPoint mp41 (475, 275);
  StopCrossSection cs41 (0.697075, 0.0933565);
  MassPointParameters mpp41 (cs41, 244972);
  StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

  MassPoint mp42 (475, 300);
  StopCrossSection cs42 (0.697075, 0.0933565);
  MassPointParameters mpp42 (cs42, 243417);
  StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

  MassPoint mp43 (500, 1);
  StopCrossSection cs43 (0.51848, 0.0693711);
  MassPointParameters mpp43 (cs43, 201307);
  StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

  MassPoint mp44 (500, 25);
  StopCrossSection cs44 (0.51848, 0.0693711);
  MassPointParameters mpp44 (cs44, 212317);
  StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

  MassPoint mp45 (500, 50);
  StopCrossSection cs45 (0.51848, 0.0693711);
  MassPointParameters mpp45 (cs45, 212826);
  StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

  MassPoint mp46 (500, 100);
  StopCrossSection cs46 (0.51848, 0.0693711);
  MassPointParameters mpp46 (cs46, 200785);
  StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

  MassPoint mp47 (500, 150);
  StopCrossSection cs47 (0.51848, 0.0693711);
  MassPointParameters mpp47 (cs47, 200723);
  StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

  MassPoint mp48 (500, 200);
  StopCrossSection cs48 (0.51848, 0.0693711);
  MassPointParameters mpp48 (cs48, 199262);
  StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

  MassPoint mp49 (500, 225);
  StopCrossSection cs49 (0.51848, 0.0693711);
  MassPointParameters mpp49 (cs49, 204058);
  StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

  MassPoint mp50 (500, 250);
  StopCrossSection cs50 (0.51848, 0.0693711);
  MassPointParameters mpp50 (cs50, 204871);
  StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

  MassPoint mp51 (500, 275);
  StopCrossSection cs51 (0.51848, 0.0693711);
  MassPointParameters mpp51 (cs51, 208231);
  StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

  MassPoint mp52 (500, 300);
  StopCrossSection cs52 (0.51848, 0.0693711);
  MassPointParameters mpp52 (cs52, 203628);
  StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

  MassPoint mp53 (500, 325);
  StopCrossSection cs53 (0.51848, 0.0693711);
  MassPointParameters mpp53 (cs53, 200916);
  StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

  MassPoint mp54 (525, 225);
  StopCrossSection cs54 (0.390303, 0.0520832);
  MassPointParameters mpp54 (cs54, 154451);
  StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

  MassPoint mp55 (525, 250);
  StopCrossSection cs55 (0.390303, 0.0520832);
  MassPointParameters mpp55 (cs55, 153176);
  StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

  MassPoint mp56 (525, 275);
  StopCrossSection cs56 (0.390303, 0.0520832);
  MassPointParameters mpp56 (cs56, 150789);
  StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

  MassPoint mp57 (525, 300);
  StopCrossSection cs57 (0.390303, 0.0520832);
  MassPointParameters mpp57 (cs57, 155321);
  StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

  MassPoint mp58 (525, 325);
  StopCrossSection cs58 (0.390303, 0.0520832);
  MassPointParameters mpp58 (cs58, 149900);
  StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

  MassPoint mp59 (525, 350);
  StopCrossSection cs59 (0.390303, 0.0520832);
  MassPointParameters mpp59 (cs59, 149457);
  StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

  MassPoint mp60 (550, 1);
  StopCrossSection cs60 (0.296128, 0.0392923);
  MassPointParameters mpp60 (cs60, 116898);
  StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

  MassPoint mp61 (550, 25);
  StopCrossSection cs61 (0.296128, 0.0392923);
  MassPointParameters mpp61 (cs61, 109757);
  StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

  MassPoint mp62 (550, 50);
  StopCrossSection cs62 (0.296128, 0.0392923);
  MassPointParameters mpp62 (cs62, 116641);
  StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

  MassPoint mp63 (550, 100);
  StopCrossSection cs63 (0.296128, 0.0392923);
  MassPointParameters mpp63 (cs63, 116513);
  StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

  MassPoint mp64 (550, 150);
  StopCrossSection cs64 (0.296128, 0.0392923);
  MassPointParameters mpp64 (cs64, 101581);
  StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

  MassPoint mp65 (550, 200);
  StopCrossSection cs65 (0.296128, 0.0392923);
  MassPointParameters mpp65 (cs65, 119930);
  StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

  MassPoint mp66 (550, 250);
  StopCrossSection cs66 (0.296128, 0.0392923);
  MassPointParameters mpp66 (cs66, 103959);
  StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

  MassPoint mp67 (550, 275);
  StopCrossSection cs67 (0.296128, 0.0392923);
  MassPointParameters mpp67 (cs67, 106931);
  StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

  MassPoint mp68 (550, 300);
  StopCrossSection cs68 (0.296128, 0.0392923);
  MassPointParameters mpp68 (cs68, 121645);
  StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

  MassPoint mp69 (550, 325);
  StopCrossSection cs69 (0.296128, 0.0392923);
  MassPointParameters mpp69 (cs69, 107690);
  StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

  MassPoint mp70 (550, 350);
  StopCrossSection cs70 (0.296128, 0.0392923);
  MassPointParameters mpp70 (cs70, 119989);
  StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

  MassPoint mp71 (550, 375);
  StopCrossSection cs71 (0.296128, 0.0392923);
  MassPointParameters mpp71 (cs71, 114977);
  StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

  MassPoint mp72 (575, 275);
  StopCrossSection cs72 (0.226118, 0.0300151);
  MassPointParameters mpp72 (cs72, 84656);
  StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

  MassPoint mp73 (575, 300);
  StopCrossSection cs73 (0.226118, 0.0300151);
  MassPointParameters mpp73 (cs73, 91742);
  StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

  MassPoint mp74 (575, 325);
  StopCrossSection cs74 (0.226118, 0.0300151);
  MassPointParameters mpp74 (cs74, 86053);
  StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

  MassPoint mp75 (575, 350);
  StopCrossSection cs75 (0.226118, 0.0300151);
  MassPointParameters mpp75 (cs75, 90118);
  StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

  MassPoint mp76 (575, 375);
  StopCrossSection cs76 (0.226118, 0.0300151);
  MassPointParameters mpp76 (cs76, 94295);
  StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

  MassPoint mp77 (575, 400);
  StopCrossSection cs77 (0.226118, 0.0300151);
  MassPointParameters mpp77 (cs77, 87651);
  StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

  MassPoint mp78 (600, 1);
  StopCrossSection cs78 (0.174599, 0.02306);
  MassPointParameters mpp78 (cs78, 75309);
  StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

  MassPoint mp79 (600, 25);
  StopCrossSection cs79 (0.174599, 0.02306);
  MassPointParameters mpp79 (cs79, 61972);
  StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

  MassPoint mp80 (600, 50);
  StopCrossSection cs80 (0.174599, 0.02306);
  MassPointParameters mpp80 (cs80, 72420);
  StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

  MassPoint mp81 (600, 100);
  StopCrossSection cs81 (0.174599, 0.02306);
  MassPointParameters mpp81 (cs81, 70813);
  StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

  MassPoint mp82 (600, 150);
  StopCrossSection cs82 (0.174599, 0.02306);
  MassPointParameters mpp82 (cs82, 68613);
  StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

  MassPoint mp83 (600, 200);
  StopCrossSection cs83 (0.174599, 0.02306);
  MassPointParameters mpp83 (cs83, 72541);
  StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

  MassPoint mp84 (600, 250);
  StopCrossSection cs84 (0.174599, 0.02306);
  MassPointParameters mpp84 (cs84, 79567);
  StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

  MassPoint mp85 (600, 300);
  StopCrossSection cs85 (0.174599, 0.02306);
  MassPointParameters mpp85 (cs85, 75914);
  StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

  MassPoint mp86 (600, 325);
  StopCrossSection cs86 (0.174599, 0.02306);
  MassPointParameters mpp86 (cs86, 77149);
  StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

  MassPoint mp87 (600, 350);
  StopCrossSection cs87 (0.174599, 0.02306);
  MassPointParameters mpp87 (cs87, 68390);
  StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

  MassPoint mp88 (600, 375);
  StopCrossSection cs88 (0.174599, 0.02306);
  MassPointParameters mpp88 (cs88, 76063);
  StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

  MassPoint mp89 (600, 400);
  StopCrossSection cs89 (0.174599, 0.02306);
  MassPointParameters mpp89 (cs89, 58228);
  StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

  MassPoint mp90 (600, 425);
  StopCrossSection cs90 (0.174599, 0.02306);
  MassPointParameters mpp90 (cs90, 68447);
  StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

  MassPoint mp91 (625, 325);
  StopCrossSection cs91 (0.136372, 0.0174504);
  MassPointParameters mpp91 (cs91, 51529);
  StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

  MassPoint mp92 (625, 350);
  StopCrossSection cs92 (0.136372, 0.0174504);
  MassPointParameters mpp92 (cs92, 57452);
  StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

  MassPoint mp93 (625, 375);
  StopCrossSection cs93 (0.136372, 0.0174504);
  MassPointParameters mpp93 (cs93, 53942);
  StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

  MassPoint mp94 (625, 400);
  StopCrossSection cs94 (0.136372, 0.0174504);
  MassPointParameters mpp94 (cs94, 59477);
  StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

  MassPoint mp95 (625, 425);
  StopCrossSection cs95 (0.136372, 0.0174504);
  MassPointParameters mpp95 (cs95, 60591);
  StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

  MassPoint mp96 (625, 450);
  StopCrossSection cs96 (0.136372, 0.0174504);
  MassPointParameters mpp96 (cs96, 55706);
  StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

  MassPoint mp97 (650, 1);
  StopCrossSection cs97 (0.107045, 0.0138336);
  MassPointParameters mpp97 (cs97, 40792);
  StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

  MassPoint mp98 (650, 25);
  StopCrossSection cs98 (0.107045, 0.0138336);
  MassPointParameters mpp98 (cs98, 40510);
  StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

  MassPoint mp99 (650, 50);
  StopCrossSection cs99 (0.107045, 0.0138336);
  MassPointParameters mpp99 (cs99, 42398);
  StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

  MassPoint mp100 (650, 100);
  StopCrossSection cs100 (0.107045, 0.0138336);
  MassPointParameters mpp100 (cs100, 41429);
  StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

  MassPoint mp101 (650, 150);
  StopCrossSection cs101 (0.107045, 0.0138336);
  MassPointParameters mpp101 (cs101, 40640);
  StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

  MassPoint mp102 (650, 200);
  StopCrossSection cs102 (0.107045, 0.0138336);
  MassPointParameters mpp102 (cs102, 45144);
  StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

  MassPoint mp103 (650, 250);
  StopCrossSection cs103 (0.107045, 0.0138336);
  MassPointParameters mpp103 (cs103, 48738);
  StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

  MassPoint mp104 (650, 300);
  StopCrossSection cs104 (0.107045, 0.0138336);
  MassPointParameters mpp104 (cs104, 44408);
  StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

  MassPoint mp105 (650, 350);
  StopCrossSection cs105 (0.107045, 0.0138336);
  MassPointParameters mpp105 (cs105, 41797);
  StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

  MassPoint mp106 (650, 375);
  StopCrossSection cs106 (0.107045, 0.0138336);
  MassPointParameters mpp106 (cs106, 43091);
  StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

  MassPoint mp107 (650, 400);
  StopCrossSection cs107 (0.107045, 0.0138336);
  MassPointParameters mpp107 (cs107, 39667);
  StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

  MassPoint mp108 (650, 425);
  StopCrossSection cs108 (0.107045, 0.0138336);
  MassPointParameters mpp108 (cs108, 46521);
  StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

  MassPoint mp109 (650, 450);
  StopCrossSection cs109 (0.107045, 0.0138336);
  MassPointParameters mpp109 (cs109, 41411);
  StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

  MassPoint mp110 (650, 475);
  StopCrossSection cs110 (0.107045, 0.0138336);
  MassPointParameters mpp110 (cs110, 41492);
  StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

  MassPoint mp111 (675, 375);
  StopCrossSection cs111 (0.0844877, 0.0110428);
  MassPointParameters mpp111 (cs111, 32940);
  StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

  MassPoint mp112 (675, 400);
  StopCrossSection cs112 (0.0844877, 0.0110428);
  MassPointParameters mpp112 (cs112, 33239);
  StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

  MassPoint mp113 (675, 425);
  StopCrossSection cs113 (0.0844877, 0.0110428);
  MassPointParameters mpp113 (cs113, 33173);
  StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

  MassPoint mp114 (675, 450);
  StopCrossSection cs114 (0.0844877, 0.0110428);
  MassPointParameters mpp114 (cs114, 34046);
  StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

  MassPoint mp115 (675, 475);
  StopCrossSection cs115 (0.0844877, 0.0110428);
  MassPointParameters mpp115 (cs115, 31403);
  StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

  MassPoint mp116 (675, 500);
  StopCrossSection cs116 (0.0844877, 0.0110428);
  MassPointParameters mpp116 (cs116, 27505);
  StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

  MassPoint mp117 (700, 1);
  StopCrossSection cs117 (0.0670476, 0.00894609);
  MassPointParameters mpp117 (cs117, 25235);
  StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

  MassPoint mp118 (700, 25);
  StopCrossSection cs118 (0.0670476, 0.00894609);
  MassPointParameters mpp118 (cs118, 26325);
  StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

  MassPoint mp119 (700, 50);
  StopCrossSection cs119 (0.0670476, 0.00894609);
  MassPointParameters mpp119 (cs119, 27536);
  StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

  MassPoint mp120 (700, 100);
  StopCrossSection cs120 (0.0670476, 0.00894609);
  MassPointParameters mpp120 (cs120, 26978);
  StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

  MassPoint mp121 (700, 150);
  StopCrossSection cs121 (0.0670476, 0.00894609);
  MassPointParameters mpp121 (cs121, 25047);
  StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

  MassPoint mp122 (700, 200);
  StopCrossSection cs122 (0.0670476, 0.00894609);
  MassPointParameters mpp122 (cs122, 28454);
  StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

  MassPoint mp123 (700, 250);
  StopCrossSection cs123 (0.0670476, 0.00894609);
  MassPointParameters mpp123 (cs123, 24587);
  StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

  MassPoint mp124 (700, 300);
  StopCrossSection cs124 (0.0670476, 0.00894609);
  MassPointParameters mpp124 (cs124, 24040);
  StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

  MassPoint mp125 (700, 350);
  StopCrossSection cs125 (0.0670476, 0.00894609);
  MassPointParameters mpp125 (cs125, 25284);
  StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

  MassPoint mp126 (700, 400);
  StopCrossSection cs126 (0.0670476, 0.00894609);
  MassPointParameters mpp126 (cs126, 25034);
  StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

  MassPoint mp127 (700, 425);
  StopCrossSection cs127 (0.0670476, 0.00894609);
  MassPointParameters mpp127 (cs127, 30318);
  StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

  MassPoint mp128 (700, 450);
  StopCrossSection cs128 (0.0670476, 0.00894609);
  MassPointParameters mpp128 (cs128, 29963);
  StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

  MassPoint mp129 (700, 475);
  StopCrossSection cs129 (0.0670476, 0.00894609);
  MassPointParameters mpp129 (cs129, 30219);
  StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

  MassPoint mp130 (700, 500);
  StopCrossSection cs130 (0.0670476, 0.00894609);
  MassPointParameters mpp130 (cs130, 25703);
  StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

  MassPoint mp131 (700, 525);
  StopCrossSection cs131 (0.0670476, 0.00894609);
  MassPointParameters mpp131 (cs131, 26718);
  StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

  MassPoint mp132 (725, 425);
  StopCrossSection cs132 (0.0536438, 0.00728504);
  MassPointParameters mpp132 (cs132, 23102);
  StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

  MassPoint mp133 (725, 450);
  StopCrossSection cs133 (0.0536438, 0.00728504);
  MassPointParameters mpp133 (cs133, 23361);
  StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

  MassPoint mp134 (725, 475);
  StopCrossSection cs134 (0.0536438, 0.00728504);
  MassPointParameters mpp134 (cs134, 19971);
  StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

  MassPoint mp135 (725, 500);
  StopCrossSection cs135 (0.0536438, 0.00728504);
  MassPointParameters mpp135 (cs135, 19884);
  StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

  MassPoint mp136 (725, 525);
  StopCrossSection cs136 (0.0536438, 0.00728504);
  MassPointParameters mpp136 (cs136, 22133);
  StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

  MassPoint mp137 (725, 550);
  StopCrossSection cs137 (0.0536438, 0.00728504);
  MassPointParameters mpp137 (cs137, 20228);
  StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

  MassPoint mp138 (750, 1);
  StopCrossSection cs138 (0.0431418, 0.00593006);
  MassPointParameters mpp138 (cs138, 16536);
  StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

  MassPoint mp139 (750, 25);
  StopCrossSection cs139 (0.0431418, 0.00593006);
  MassPointParameters mpp139 (cs139, 19935);
  StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

  MassPoint mp140 (750, 50);
  StopCrossSection cs140 (0.0431418, 0.00593006);
  MassPointParameters mpp140 (cs140, 17075);
  StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

  MassPoint mp141 (750, 100);
  StopCrossSection cs141 (0.0431418, 0.00593006);
  MassPointParameters mpp141 (cs141, 19905);
  StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

  MassPoint mp142 (750, 150);
  StopCrossSection cs142 (0.0431418, 0.00593006);
  MassPointParameters mpp142 (cs142, 18883);
  StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

  MassPoint mp143 (750, 200);
  StopCrossSection cs143 (0.0431418, 0.00593006);
  MassPointParameters mpp143 (cs143, 19627);
  StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

  MassPoint mp144 (750, 250);
  StopCrossSection cs144 (0.0431418, 0.00593006);
  MassPointParameters mpp144 (cs144, 21012);
  StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

  MassPoint mp145 (750, 300);
  StopCrossSection cs145 (0.0431418, 0.00593006);
  MassPointParameters mpp145 (cs145, 21077);
  StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

  MassPoint mp146 (750, 350);
  StopCrossSection cs146 (0.0431418, 0.00593006);
  MassPointParameters mpp146 (cs146, 16686);
  StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

  MassPoint mp147 (750, 400);
  StopCrossSection cs147 (0.0431418, 0.00593006);
  MassPointParameters mpp147 (cs147, 19529);
  StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

  MassPoint mp148 (750, 450);
  StopCrossSection cs148 (0.0431418, 0.00593006);
  MassPointParameters mpp148 (cs148, 20074);
  StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

  MassPoint mp149 (750, 475);
  StopCrossSection cs149 (0.0431418, 0.00593006);
  MassPointParameters mpp149 (cs149, 18360);
  StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

  MassPoint mp150 (750, 500);
  StopCrossSection cs150 (0.0431418, 0.00593006);
  MassPointParameters mpp150 (cs150, 21264);
  StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

  MassPoint mp151 (750, 525);
  StopCrossSection cs151 (0.0431418, 0.00593006);
  MassPointParameters mpp151 (cs151, 20903);
  StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

  MassPoint mp152 (750, 550);
  StopCrossSection cs152 (0.0431418, 0.00593006);
  MassPointParameters mpp152 (cs152, 19826);
  StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

  MassPoint mp153 (750, 575);
  StopCrossSection cs153 (0.0431418, 0.00593006);
  MassPointParameters mpp153 (cs153, 22460);
  StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

  MassPoint mp154 (775, 475);
  StopCrossSection cs154 (0.0348796, 0.00486909);
  MassPointParameters mpp154 (cs154, 21766);
  StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

  MassPoint mp155 (775, 500);
  StopCrossSection cs155 (0.0348796, 0.00486909);
  MassPointParameters mpp155 (cs155, 20114);
  StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

  MassPoint mp156 (775, 525);
  StopCrossSection cs156 (0.0348796, 0.00486909);
  MassPointParameters mpp156 (cs156, 23384);
  StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

  MassPoint mp157 (775, 550);
  StopCrossSection cs157 (0.0348796, 0.00486909);
  MassPointParameters mpp157 (cs157, 19371);
  StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

  MassPoint mp158 (775, 575);
  StopCrossSection cs158 (0.0348796, 0.00486909);
  MassPointParameters mpp158 (cs158, 18251);
  StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

  MassPoint mp159 (775, 600);
  StopCrossSection cs159 (0.0348796, 0.00486909);
  MassPointParameters mpp159 (cs159, 19486);
  StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

  MassPoint mp160 (800, 1);
  StopCrossSection cs160 (0.0283338, 0.00401518);
  MassPointParameters mpp160 (cs160, 17504);
  StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

  MassPoint mp161 (800, 25);
  StopCrossSection cs161 (0.0283338, 0.00401518);
  MassPointParameters mpp161 (cs161, 19155);
  StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

  MassPoint mp162 (800, 50);
  StopCrossSection cs162 (0.0283338, 0.00401518);
  MassPointParameters mpp162 (cs162, 22741);
  StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

  MassPoint mp163 (800, 100);
  StopCrossSection cs163 (0.0283338, 0.00401518);
  MassPointParameters mpp163 (cs163, 21713);
  StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

  MassPoint mp164 (800, 150);
  StopCrossSection cs164 (0.0283338, 0.00401518);
  MassPointParameters mpp164 (cs164, 19666);
  StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

  MassPoint mp165 (800, 200);
  StopCrossSection cs165 (0.0283338, 0.00401518);
  MassPointParameters mpp165 (cs165, 20759);
  StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

  MassPoint mp166 (800, 250);
  StopCrossSection cs166 (0.0283338, 0.00401518);
  MassPointParameters mpp166 (cs166, 18986);
  StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

  MassPoint mp167 (800, 300);
  StopCrossSection cs167 (0.0283338, 0.00401518);
  MassPointParameters mpp167 (cs167, 17993);
  StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

  MassPoint mp168 (800, 350);
  StopCrossSection cs168 (0.0283338, 0.00401518);
  MassPointParameters mpp168 (cs168, 21190);
  StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

  MassPoint mp169 (800, 400);
  StopCrossSection cs169 (0.0283338, 0.00401518);
  MassPointParameters mpp169 (cs169, 20812);
  StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

  MassPoint mp170 (800, 450);
  StopCrossSection cs170 (0.0283338, 0.00401518);
  MassPointParameters mpp170 (cs170, 16138);
  StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

  MassPoint mp171 (800, 500);
  StopCrossSection cs171 (0.0283338, 0.00401518);
  MassPointParameters mpp171 (cs171, 19483);
  StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

  MassPoint mp172 (800, 525);
  StopCrossSection cs172 (0.0283338, 0.00401518);
  MassPointParameters mpp172 (cs172, 19655);
  StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

  MassPoint mp173 (800, 550);
  StopCrossSection cs173 (0.0283338, 0.00401518);
  MassPointParameters mpp173 (cs173, 17671);
  StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

  MassPoint mp174 (800, 575);
  StopCrossSection cs174 (0.0283338, 0.00401518);
  MassPointParameters mpp174 (cs174, 20756);
  StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

  MassPoint mp175 (800, 600);
  StopCrossSection cs175 (0.0283338, 0.00401518);
  MassPointParameters mpp175 (cs175, 19047);
  StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

  MassPoint mp176 (800, 625);
  StopCrossSection cs176 (0.0283338, 0.00401518);
  MassPointParameters mpp176 (cs176, 18849);
  StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

  MassPoint mp177 (825, 525);
  StopCrossSection cs177 (0.0230866, 0.00333435);
  MassPointParameters mpp177 (cs177, 20604);
  StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

  MassPoint mp178 (825, 550);
  StopCrossSection cs178 (0.0230866, 0.00333435);
  MassPointParameters mpp178 (cs178, 20741);
  StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

  MassPoint mp179 (825, 575);
  StopCrossSection cs179 (0.0230866, 0.00333435);
  MassPointParameters mpp179 (cs179, 20101);
  StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

  MassPoint mp180 (825, 600);
  StopCrossSection cs180 (0.0230866, 0.00333435);
  MassPointParameters mpp180 (cs180, 19675);
  StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

  MassPoint mp181 (825, 625);
  StopCrossSection cs181 (0.0230866, 0.00333435);
  MassPointParameters mpp181 (cs181, 18161);
  StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

  MassPoint mp182 (825, 650);
  StopCrossSection cs182 (0.0230866, 0.00333435);
  MassPointParameters mpp182 (cs182, 21313);
  StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

  MassPoint mp183 (850, 1);
  StopCrossSection cs183 (0.0189612, 0.00278768);
  MassPointParameters mpp183 (cs183, 19645);
  StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

  MassPoint mp184 (850, 25);
  StopCrossSection cs184 (0.0189612, 0.00278768);
  MassPointParameters mpp184 (cs184, 17931);
  StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

  MassPoint mp185 (850, 50);
  StopCrossSection cs185 (0.0189612, 0.00278768);
  MassPointParameters mpp185 (cs185, 24577);
  StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

  MassPoint mp186 (850, 100);
  StopCrossSection cs186 (0.0189612, 0.00278768);
  MassPointParameters mpp186 (cs186, 19734);
  StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

  MassPoint mp187 (850, 150);
  StopCrossSection cs187 (0.0189612, 0.00278768);
  MassPointParameters mpp187 (cs187, 19040);
  StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

  MassPoint mp188 (850, 200);
  StopCrossSection cs188 (0.0189612, 0.00278768);
  MassPointParameters mpp188 (cs188, 18390);
  StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

  MassPoint mp189 (850, 250);
  StopCrossSection cs189 (0.0189612, 0.00278768);
  MassPointParameters mpp189 (cs189, 17991);
  StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

  MassPoint mp190 (850, 300);
  StopCrossSection cs190 (0.0189612, 0.00278768);
  MassPointParameters mpp190 (cs190, 24589);
  StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

  MassPoint mp191 (850, 350);
  StopCrossSection cs191 (0.0189612, 0.00278768);
  MassPointParameters mpp191 (cs191, 22171);
  StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

  MassPoint mp192 (850, 400);
  StopCrossSection cs192 (0.0189612, 0.00278768);
  MassPointParameters mpp192 (cs192, 19987);
  StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

  MassPoint mp193 (850, 450);
  StopCrossSection cs193 (0.0189612, 0.00278768);
  MassPointParameters mpp193 (cs193, 22249);
  StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

  MassPoint mp194 (850, 500);
  StopCrossSection cs194 (0.0189612, 0.00278768);
  MassPointParameters mpp194 (cs194, 16688);
  StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

  MassPoint mp195 (850, 550);
  StopCrossSection cs195 (0.0189612, 0.00278768);
  MassPointParameters mpp195 (cs195, 19826);
  StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

  MassPoint mp196 (850, 575);
  StopCrossSection cs196 (0.0189612, 0.00278768);
  MassPointParameters mpp196 (cs196, 18930);
  StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

  MassPoint mp197 (850, 600);
  StopCrossSection cs197 (0.0189612, 0.00278768);
  MassPointParameters mpp197 (cs197, 21920);
  StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

  MassPoint mp198 (850, 625);
  StopCrossSection cs198 (0.0189612, 0.00278768);
  MassPointParameters mpp198 (cs198, 18747);
  StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

  MassPoint mp199 (850, 650);
  StopCrossSection cs199 (0.0189612, 0.00278768);
  MassPointParameters mpp199 (cs199, 18113);
  StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

  MassPoint mp200 (875, 575);
  StopCrossSection cs200 (0.015625, 0.00233698);
  MassPointParameters mpp200 (cs200, 18254);
  StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

  MassPoint mp201 (875, 600);
  StopCrossSection cs201 (0.015625, 0.00233698);
  MassPointParameters mpp201 (cs201, 19643);
  StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

  MassPoint mp202 (875, 625);
  StopCrossSection cs202 (0.015625, 0.00233698);
  MassPointParameters mpp202 (cs202, 19809);
  StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

  MassPoint mp203 (875, 650);
  StopCrossSection cs203 (0.015625, 0.00233698);
  MassPointParameters mpp203 (cs203, 20077);
  StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

  MassPoint mp204 (900, 1);
  StopCrossSection cs204 (0.0128895, 0.00195954);
  MassPointParameters mpp204 (cs204, 20895);
  StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

  MassPoint mp205 (900, 25);
  StopCrossSection cs205 (0.0128895, 0.00195954);
  MassPointParameters mpp205 (cs205, 18456);
  StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

  MassPoint mp206 (900, 50);
  StopCrossSection cs206 (0.0128895, 0.00195954);
  MassPointParameters mpp206 (cs206, 19265);
  StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

  MassPoint mp207 (900, 100);
  StopCrossSection cs207 (0.0128895, 0.00195954);
  MassPointParameters mpp207 (cs207, 19411);
  StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

  MassPoint mp208 (900, 150);
  StopCrossSection cs208 (0.0128895, 0.00195954);
  MassPointParameters mpp208 (cs208, 18558);
  StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

  MassPoint mp209 (900, 200);
  StopCrossSection cs209 (0.0128895, 0.00195954);
  MassPointParameters mpp209 (cs209, 19514);
  StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

  MassPoint mp210 (900, 250);
  StopCrossSection cs210 (0.0128895, 0.00195954);
  MassPointParameters mpp210 (cs210, 19523);
  StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

  MassPoint mp211 (900, 300);
  StopCrossSection cs211 (0.0128895, 0.00195954);
  MassPointParameters mpp211 (cs211, 20078);
  StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

  MassPoint mp212 (900, 350);
  StopCrossSection cs212 (0.0128895, 0.00195954);
  MassPointParameters mpp212 (cs212, 14159);
  StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

  MassPoint mp213 (900, 400);
  StopCrossSection cs213 (0.0128895, 0.00195954);
  MassPointParameters mpp213 (cs213, 17589);
  StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

  MassPoint mp214 (900, 450);
  StopCrossSection cs214 (0.0128895, 0.00195954);
  MassPointParameters mpp214 (cs214, 19315);
  StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

  MassPoint mp215 (900, 500);
  StopCrossSection cs215 (0.0128895, 0.00195954);
  MassPointParameters mpp215 (cs215, 18717);
  StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

  MassPoint mp216 (900, 550);
  StopCrossSection cs216 (0.0128895, 0.00195954);
  MassPointParameters mpp216 (cs216, 18629);
  StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

  MassPoint mp217 (900, 600);
  StopCrossSection cs217 (0.0128895, 0.00195954);
  MassPointParameters mpp217 (cs217, 16608);
  StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

  MassPoint mp218 (900, 625);
  StopCrossSection cs218 (0.0128895, 0.00195954);
  MassPointParameters mpp218 (cs218, 21621);
  StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

  MassPoint mp219 (900, 650);
  StopCrossSection cs219 (0.0128895, 0.00195954);
  MassPointParameters mpp219 (cs219, 23983);
  StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

  MassPoint mp220 (925, 625);
  StopCrossSection cs220 (0.0106631, 0.00165071);
  MassPointParameters mpp220 (cs220, 15711);
  StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

  MassPoint mp221 (925, 650);
  StopCrossSection cs221 (0.0106631, 0.00165071);
  MassPointParameters mpp221 (cs221, 19644);
  StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

  MassPoint mp222 (950, 1);
  StopCrossSection cs222 (0.00883465, 0.0013886);
  MassPointParameters mpp222 (cs222, 22178);
  StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

  MassPoint mp223 (950, 25);
  StopCrossSection cs223 (0.00883465, 0.0013886);
  MassPointParameters mpp223 (cs223, 16618);
  StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

  MassPoint mp224 (950, 50);
  StopCrossSection cs224 (0.00883465, 0.0013886);
  MassPointParameters mpp224 (cs224, 19464);
  StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

  MassPoint mp225 (950, 100);
  StopCrossSection cs225 (0.00883465, 0.0013886);
  MassPointParameters mpp225 (cs225, 16899);
  StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

  MassPoint mp226 (950, 150);
  StopCrossSection cs226 (0.00883465, 0.0013886);
  MassPointParameters mpp226 (cs226, 22678);
  StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

  MassPoint mp227 (950, 200);
  StopCrossSection cs227 (0.00883465, 0.0013886);
  MassPointParameters mpp227 (cs227, 23038);
  StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

  MassPoint mp228 (950, 250);
  StopCrossSection cs228 (0.00883465, 0.0013886);
  MassPointParameters mpp228 (cs228, 18749);
  StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

  MassPoint mp229 (950, 300);
  StopCrossSection cs229 (0.00883465, 0.0013886);
  MassPointParameters mpp229 (cs229, 18127);
  StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

  MassPoint mp230 (950, 350);
  StopCrossSection cs230 (0.00883465, 0.0013886);
  MassPointParameters mpp230 (cs230, 19831);
  StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

  MassPoint mp231 (950, 400);
  StopCrossSection cs231 (0.00883465, 0.0013886);
  MassPointParameters mpp231 (cs231, 18863);
  StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

  MassPoint mp232 (950, 450);
  StopCrossSection cs232 (0.00883465, 0.0013886);
  MassPointParameters mpp232 (cs232, 19757);
  StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

  MassPoint mp233 (950, 500);
  StopCrossSection cs233 (0.00883465, 0.0013886);
  MassPointParameters mpp233 (cs233, 18387);
  StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

  MassPoint mp234 (950, 550);
  StopCrossSection cs234 (0.00883465, 0.0013886);
  MassPointParameters mpp234 (cs234, 17909);
  StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

  MassPoint mp235 (950, 600);
  StopCrossSection cs235 (0.00883465, 0.0013886);
  MassPointParameters mpp235 (cs235, 19142);
  StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

  MassPoint mp236 (950, 650);
  StopCrossSection cs236 (0.00883465, 0.0013886);
  MassPointParameters mpp236 (cs236, 20012);
  StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

  MassPoint mp237 (1000, 1);
  StopCrossSection cs237 (0.00615134, 0.00100238);
  MassPointParameters mpp237 (cs237, 19680);
  StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

  MassPoint mp238 (1000, 25);
  StopCrossSection cs238 (0.00615134, 0.00100238);
  MassPointParameters mpp238 (cs238, 20929);
  StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

  MassPoint mp239 (1000, 50);
  StopCrossSection cs239 (0.00615134, 0.00100238);
  MassPointParameters mpp239 (cs239, 17755);
  StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

  MassPoint mp240 (1000, 100);
  StopCrossSection cs240 (0.00615134, 0.00100238);
  MassPointParameters mpp240 (cs240, 21340);
  StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

  MassPoint mp241 (1000, 150);
  StopCrossSection cs241 (0.00615134, 0.00100238);
  MassPointParameters mpp241 (cs241, 22839);
  StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

  MassPoint mp242 (1000, 200);
  StopCrossSection cs242 (0.00615134, 0.00100238);
  MassPointParameters mpp242 (cs242, 19978);
  StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

  MassPoint mp243 (1000, 250);
  StopCrossSection cs243 (0.00615134, 0.00100238);
  MassPointParameters mpp243 (cs243, 20383);
  StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

  MassPoint mp244 (1000, 300);
  StopCrossSection cs244 (0.00615134, 0.00100238);
  MassPointParameters mpp244 (cs244, 20197);
  StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

  MassPoint mp245 (1000, 350);
  StopCrossSection cs245 (0.00615134, 0.00100238);
  MassPointParameters mpp245 (cs245, 21754);
  StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

  MassPoint mp246 (1000, 400);
  StopCrossSection cs246 (0.00615134, 0.00100238);
  MassPointParameters mpp246 (cs246, 19783);
  StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

  MassPoint mp247 (1000, 450);
  StopCrossSection cs247 (0.00615134, 0.00100238);
  MassPointParameters mpp247 (cs247, 23064);
  StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

  MassPoint mp248 (1000, 500);
  StopCrossSection cs248 (0.00615134, 0.00100238);
  MassPointParameters mpp248 (cs248, 18958);
  StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

  MassPoint mp249 (1000, 550);
  StopCrossSection cs249 (0.00615134, 0.00100238);
  MassPointParameters mpp249 (cs249, 18070);
  StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

  MassPoint mp250 (1000, 600);
  StopCrossSection cs250 (0.00615134, 0.00100238);
  MassPointParameters mpp250 (cs250, 17746);
  StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

  MassPoint mp251 (1000, 650);
  StopCrossSection cs251 (0.00615134, 0.00100238);
  MassPointParameters mpp251 (cs251, 18073);
  StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

  MassPoint mp252 (1050, 1);
  StopCrossSection cs252 (0.00432261, 0.000725589);
  MassPointParameters mpp252 (cs252, 18108);
  StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

  MassPoint mp253 (1050, 25);
  StopCrossSection cs253 (0.00432261, 0.000725589);
  MassPointParameters mpp253 (cs253, 26158);
  StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

  MassPoint mp254 (1050, 50);
  StopCrossSection cs254 (0.00432261, 0.000725589);
  MassPointParameters mpp254 (cs254, 20229);
  StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

  MassPoint mp255 (1050, 100);
  StopCrossSection cs255 (0.00432261, 0.000725589);
  MassPointParameters mpp255 (cs255, 19412);
  StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

  MassPoint mp256 (1050, 150);
  StopCrossSection cs256 (0.00432261, 0.000725589);
  MassPointParameters mpp256 (cs256, 19516);
  StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

  MassPoint mp257 (1050, 200);
  StopCrossSection cs257 (0.00432261, 0.000725589);
  MassPointParameters mpp257 (cs257, 22046);
  StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

  MassPoint mp258 (1050, 250);
  StopCrossSection cs258 (0.00432261, 0.000725589);
  MassPointParameters mpp258 (cs258, 19558);
  StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

  MassPoint mp259 (1050, 300);
  StopCrossSection cs259 (0.00432261, 0.000725589);
  MassPointParameters mpp259 (cs259, 20622);
  StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

  MassPoint mp260 (1050, 350);
  StopCrossSection cs260 (0.00432261, 0.000725589);
  MassPointParameters mpp260 (cs260, 19063);
  StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

  MassPoint mp261 (1050, 400);
  StopCrossSection cs261 (0.00432261, 0.000725589);
  MassPointParameters mpp261 (cs261, 19124);
  StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

  MassPoint mp262 (1050, 450);
  StopCrossSection cs262 (0.00432261, 0.000725589);
  MassPointParameters mpp262 (cs262, 18550);
  StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

  MassPoint mp263 (1050, 500);
  StopCrossSection cs263 (0.00432261, 0.000725589);
  MassPointParameters mpp263 (cs263, 22141);
  StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

  MassPoint mp264 (1050, 550);
  StopCrossSection cs264 (0.00432261, 0.000725589);
  MassPointParameters mpp264 (cs264, 19868);
  StopNeutralinoMap.insert(std::make_pair(mp264, mpp264));

  MassPoint mp265 (1050, 600);
  StopCrossSection cs265 (0.00432261, 0.000725589);
  MassPointParameters mpp265 (cs265, 16778);
  StopNeutralinoMap.insert(std::make_pair(mp265, mpp265));

  MassPoint mp266 (1050, 650);
  StopCrossSection cs266 (0.00432261, 0.000725589);
  MassPointParameters mpp266 (cs266, 21079);
  StopNeutralinoMap.insert(std::make_pair(mp266, mpp266));

  MassPoint mp267 (1100, 1);
  StopCrossSection cs267 (0.00307413, 0.000532983);
  MassPointParameters mpp267 (cs267, 20351);
  StopNeutralinoMap.insert(std::make_pair(mp267, mpp267));

  MassPoint mp268 (1100, 25);
  StopCrossSection cs268 (0.00307413, 0.000532983);
  MassPointParameters mpp268 (cs268, 21437);
  StopNeutralinoMap.insert(std::make_pair(mp268, mpp268));

  MassPoint mp269 (1100, 50);
  StopCrossSection cs269 (0.00307413, 0.000532983);
  MassPointParameters mpp269 (cs269, 19931);
  StopNeutralinoMap.insert(std::make_pair(mp269, mpp269));

  MassPoint mp270 (1100, 100);
  StopCrossSection cs270 (0.00307413, 0.000532983);
  MassPointParameters mpp270 (cs270, 22164);
  StopNeutralinoMap.insert(std::make_pair(mp270, mpp270));

  MassPoint mp271 (1100, 150);
  StopCrossSection cs271 (0.00307413, 0.000532983);
  MassPointParameters mpp271 (cs271, 17929);
  StopNeutralinoMap.insert(std::make_pair(mp271, mpp271));

  MassPoint mp272 (1100, 200);
  StopCrossSection cs272 (0.00307413, 0.000532983);
  MassPointParameters mpp272 (cs272, 18326);
  StopNeutralinoMap.insert(std::make_pair(mp272, mpp272));

  MassPoint mp273 (1100, 250);
  StopCrossSection cs273 (0.00307413, 0.000532983);
  MassPointParameters mpp273 (cs273, 20307);
  StopNeutralinoMap.insert(std::make_pair(mp273, mpp273));

  MassPoint mp274 (1100, 300);
  StopCrossSection cs274 (0.00307413, 0.000532983);
  MassPointParameters mpp274 (cs274, 21970);
  StopNeutralinoMap.insert(std::make_pair(mp274, mpp274));

  MassPoint mp275 (1100, 350);
  StopCrossSection cs275 (0.00307413, 0.000532983);
  MassPointParameters mpp275 (cs275, 18356);
  StopNeutralinoMap.insert(std::make_pair(mp275, mpp275));

  MassPoint mp276 (1100, 400);
  StopCrossSection cs276 (0.00307413, 0.000532983);
  MassPointParameters mpp276 (cs276, 22306);
  StopNeutralinoMap.insert(std::make_pair(mp276, mpp276));

  MassPoint mp277 (1100, 450);
  StopCrossSection cs277 (0.00307413, 0.000532983);
  MassPointParameters mpp277 (cs277, 21513);
  StopNeutralinoMap.insert(std::make_pair(mp277, mpp277));

  MassPoint mp278 (1100, 500);
  StopCrossSection cs278 (0.00307413, 0.000532983);
  MassPointParameters mpp278 (cs278, 19872);
  StopNeutralinoMap.insert(std::make_pair(mp278, mpp278));

  MassPoint mp279 (1100, 550);
  StopCrossSection cs279 (0.00307413, 0.000532983);
  MassPointParameters mpp279 (cs279, 18031);
  StopNeutralinoMap.insert(std::make_pair(mp279, mpp279));

  MassPoint mp280 (1100, 600);
  StopCrossSection cs280 (0.00307413, 0.000532983);
  MassPointParameters mpp280 (cs280, 18318);
  StopNeutralinoMap.insert(std::make_pair(mp280, mpp280));

  MassPoint mp281 (1100, 650);
  StopCrossSection cs281 (0.00307413, 0.000532983);
  MassPointParameters mpp281 (cs281, 17744);
  StopNeutralinoMap.insert(std::make_pair(mp281, mpp281));

  MassPoint mp282 (1150, 1);
  StopCrossSection cs282 (0.00221047, 0.000396247);
  MassPointParameters mpp282 (cs282, 21332);
  StopNeutralinoMap.insert(std::make_pair(mp282, mpp282));

  MassPoint mp283 (1150, 25);
  StopCrossSection cs283 (0.00221047, 0.000396247);
  MassPointParameters mpp283 (cs283, 20912);
  StopNeutralinoMap.insert(std::make_pair(mp283, mpp283));

  MassPoint mp284 (1150, 50);
  StopCrossSection cs284 (0.00221047, 0.000396247);
  MassPointParameters mpp284 (cs284, 20187);
  StopNeutralinoMap.insert(std::make_pair(mp284, mpp284));

  MassPoint mp285 (1150, 100);
  StopCrossSection cs285 (0.00221047, 0.000396247);
  MassPointParameters mpp285 (cs285, 20362);
  StopNeutralinoMap.insert(std::make_pair(mp285, mpp285));

  MassPoint mp286 (1150, 150);
  StopCrossSection cs286 (0.00221047, 0.000396247);
  MassPointParameters mpp286 (cs286, 20587);
  StopNeutralinoMap.insert(std::make_pair(mp286, mpp286));

  MassPoint mp287 (1150, 200);
  StopCrossSection cs287 (0.00221047, 0.000396247);
  MassPointParameters mpp287 (cs287, 18911);
  StopNeutralinoMap.insert(std::make_pair(mp287, mpp287));

  MassPoint mp288 (1150, 250);
  StopCrossSection cs288 (0.00221047, 0.000396247);
  MassPointParameters mpp288 (cs288, 17993);
  StopNeutralinoMap.insert(std::make_pair(mp288, mpp288));

  MassPoint mp289 (1150, 300);
  StopCrossSection cs289 (0.00221047, 0.000396247);
  MassPointParameters mpp289 (cs289, 18089);
  StopNeutralinoMap.insert(std::make_pair(mp289, mpp289));

  MassPoint mp290 (1150, 350);
  StopCrossSection cs290 (0.00221047, 0.000396247);
  MassPointParameters mpp290 (cs290, 20620);
  StopNeutralinoMap.insert(std::make_pair(mp290, mpp290));

  MassPoint mp291 (1150, 400);
  StopCrossSection cs291 (0.00221047, 0.000396247);
  MassPointParameters mpp291 (cs291, 16184);
  StopNeutralinoMap.insert(std::make_pair(mp291, mpp291));

  MassPoint mp292 (1150, 450);
  StopCrossSection cs292 (0.00221047, 0.000396247);
  MassPointParameters mpp292 (cs292, 16188);
  StopNeutralinoMap.insert(std::make_pair(mp292, mpp292));

  MassPoint mp293 (1150, 500);
  StopCrossSection cs293 (0.00221047, 0.000396247);
  MassPointParameters mpp293 (cs293, 20069);
  StopNeutralinoMap.insert(std::make_pair(mp293, mpp293));

  MassPoint mp294 (1150, 550);
  StopCrossSection cs294 (0.00221047, 0.000396247);
  MassPointParameters mpp294 (cs294, 21673);
  StopNeutralinoMap.insert(std::make_pair(mp294, mpp294));

  MassPoint mp295 (1150, 600);
  StopCrossSection cs295 (0.00221047, 0.000396247);
  MassPointParameters mpp295 (cs295, 20847);
  StopNeutralinoMap.insert(std::make_pair(mp295, mpp295));

  MassPoint mp296 (1150, 650);
  StopCrossSection cs296 (0.00221047, 0.000396247);
  MassPointParameters mpp296 (cs296, 19458);
  StopNeutralinoMap.insert(std::make_pair(mp296, mpp296));

  MassPoint mp297 (1200, 1);
  StopCrossSection cs297 (0.00159844, 0.000296045);
  MassPointParameters mpp297 (cs297, 21109);
  StopNeutralinoMap.insert(std::make_pair(mp297, mpp297));

  MassPoint mp298 (1200, 25);
  StopCrossSection cs298 (0.00159844, 0.000296045);
  MassPointParameters mpp298 (cs298, 19763);
  StopNeutralinoMap.insert(std::make_pair(mp298, mpp298));

  MassPoint mp299 (1200, 50);
  StopCrossSection cs299 (0.00159844, 0.000296045);
  MassPointParameters mpp299 (cs299, 19394);
  StopNeutralinoMap.insert(std::make_pair(mp299, mpp299));

  MassPoint mp300 (1200, 100);
  StopCrossSection cs300 (0.00159844, 0.000296045);
  MassPointParameters mpp300 (cs300, 18532);
  StopNeutralinoMap.insert(std::make_pair(mp300, mpp300));

  MassPoint mp301 (1200, 150);
  StopCrossSection cs301 (0.00159844, 0.000296045);
  MassPointParameters mpp301 (cs301, 21660);
  StopNeutralinoMap.insert(std::make_pair(mp301, mpp301));

  MassPoint mp302 (1200, 200);
  StopCrossSection cs302 (0.00159844, 0.000296045);
  MassPointParameters mpp302 (cs302, 20458);
  StopNeutralinoMap.insert(std::make_pair(mp302, mpp302));

  MassPoint mp303 (1200, 250);
  StopCrossSection cs303 (0.00159844, 0.000296045);
  MassPointParameters mpp303 (cs303, 20529);
  StopNeutralinoMap.insert(std::make_pair(mp303, mpp303));

  MassPoint mp304 (1200, 300);
  StopCrossSection cs304 (0.00159844, 0.000296045);
  MassPointParameters mpp304 (cs304, 19541);
  StopNeutralinoMap.insert(std::make_pair(mp304, mpp304));

  MassPoint mp305 (1200, 350);
  StopCrossSection cs305 (0.00159844, 0.000296045);
  MassPointParameters mpp305 (cs305, 18092);
  StopNeutralinoMap.insert(std::make_pair(mp305, mpp305));

  MassPoint mp306 (1200, 400);
  StopCrossSection cs306 (0.00159844, 0.000296045);
  MassPointParameters mpp306 (cs306, 22464);
  StopNeutralinoMap.insert(std::make_pair(mp306, mpp306));

  MassPoint mp307 (1200, 450);
  StopCrossSection cs307 (0.00159844, 0.000296045);
  MassPointParameters mpp307 (cs307, 18052);
  StopNeutralinoMap.insert(std::make_pair(mp307, mpp307));

  MassPoint mp308 (1200, 500);
  StopCrossSection cs308 (0.00159844, 0.000296045);
  MassPointParameters mpp308 (cs308, 20497);
  StopNeutralinoMap.insert(std::make_pair(mp308, mpp308));

  MassPoint mp309 (1200, 550);
  StopCrossSection cs309 (0.00159844, 0.000296045);
  MassPointParameters mpp309 (cs309, 19295);
  StopNeutralinoMap.insert(std::make_pair(mp309, mpp309));

  MassPoint mp310 (1200, 600);
  StopCrossSection cs310 (0.00159844, 0.000296045);
  MassPointParameters mpp310 (cs310, 17849);
  StopNeutralinoMap.insert(std::make_pair(mp310, mpp310));

  MassPoint mp311 (1200, 650);
  StopCrossSection cs311 (0.00159844, 0.000296045);
  MassPointParameters mpp311 (cs311, 19401);
  StopNeutralinoMap.insert(std::make_pair(mp311, mpp311));

    MassPoint mpisr0 (200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 1.06891));

    MassPoint mpisr1 (200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 1.06833));

    MassPoint mpisr2 (250, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 1.07577));

    MassPoint mpisr3 (250, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 1.07685));

    MassPoint mpisr4 (250, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 1.07761));

    MassPoint mpisr5 (300, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 1.08275));

    MassPoint mpisr6 (300, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 1.08274));

    MassPoint mpisr7 (300, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 1.08331));

    MassPoint mpisr8 (300, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 1.08483));

    MassPoint mpisr9 (350, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 1.08788));

    MassPoint mpisr10 (350, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 1.08871));

    MassPoint mpisr11 (350, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 1.0886));

    MassPoint mpisr12 (350, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 1.08749));

    MassPoint mpisr13 (350, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 1.08628));

    MassPoint mpisr14 (400, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 1.09468));

    MassPoint mpisr15 (400, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 1.09268));

    MassPoint mpisr16 (400, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 1.0914));

    MassPoint mpisr17 (400, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 1.09572));

    MassPoint mpisr18 (400, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 1.09346));

    MassPoint mpisr19 (400, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 1.09401));

    MassPoint mpisr20 (400, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 1.09285));

    MassPoint mpisr21 (400, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 1.09248));

    MassPoint mpisr22 (425, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 1.09384));

    MassPoint mpisr23 (425, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 1.09428));

    MassPoint mpisr24 (425, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 1.09588));

    MassPoint mpisr25 (425, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 1.09747));

    MassPoint mpisr26 (425, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 1.0952));

    MassPoint mpisr27 (450, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 1.09855));

    MassPoint mpisr28 (450, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 1.09658));

    MassPoint mpisr29 (450, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 1.09719));

    MassPoint mpisr30 (450, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 1.09859));

    MassPoint mpisr31 (450, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 1.09509));

    MassPoint mpisr32 (450, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 1.09524));

    MassPoint mpisr33 (450, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 1.09607));

    MassPoint mpisr34 (450, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 1.09668));

    MassPoint mpisr35 (450, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 1.09805));

    MassPoint mpisr36 (450, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 1.09655));

    MassPoint mpisr37 (475, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 1.0959));

    MassPoint mpisr38 (475, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 1.09817));

    MassPoint mpisr39 (475, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 1.10221));

    MassPoint mpisr40 (475, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 1.09938));

    MassPoint mpisr41 (475, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 1.10033));

    MassPoint mpisr42 (475, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 1.09846));

    MassPoint mpisr43 (500, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 1.10133));

    MassPoint mpisr44 (500, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 1.10302));

    MassPoint mpisr45 (500, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 1.10132));

    MassPoint mpisr46 (500, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 1.09943));

    MassPoint mpisr47 (500, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 1.09906));

    MassPoint mpisr48 (500, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 1.10149));

    MassPoint mpisr49 (500, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 1.10156));

    MassPoint mpisr50 (500, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 1.09747));

    MassPoint mpisr51 (500, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 1.1011));

    MassPoint mpisr52 (500, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 1.09801));

    MassPoint mpisr53 (500, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 1.10045));

    MassPoint mpisr54 (525, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 1.10148));

    MassPoint mpisr55 (525, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 1.1036));

    MassPoint mpisr56 (525, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 1.10016));

    MassPoint mpisr57 (525, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 1.09976));

    MassPoint mpisr58 (525, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 1.09895));

    MassPoint mpisr59 (525, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 1.1069));

    MassPoint mpisr60 (550, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 1.1011));

    MassPoint mpisr61 (550, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 1.10576));

    MassPoint mpisr62 (550, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 1.10098));

    MassPoint mpisr63 (550, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 1.10257));

    MassPoint mpisr64 (550, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 1.10358));

    MassPoint mpisr65 (550, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 1.10942));

    MassPoint mpisr66 (550, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 1.09922));

    MassPoint mpisr67 (550, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 1.1068));

    MassPoint mpisr68 (550, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 1.10184));

    MassPoint mpisr69 (550, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 1.1025));

    MassPoint mpisr70 (550, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 1.10124));

    MassPoint mpisr71 (550, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 1.105));

    MassPoint mpisr72 (575, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 1.10694));

    MassPoint mpisr73 (575, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 1.10335));

    MassPoint mpisr74 (575, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 1.10283));

    MassPoint mpisr75 (575, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 1.09966));

    MassPoint mpisr76 (575, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 1.09923));

    MassPoint mpisr77 (575, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 1.1081));

    MassPoint mpisr78 (600, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 1.0989));

    MassPoint mpisr79 (600, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 1.10556));

    MassPoint mpisr80 (600, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 1.10378));

    MassPoint mpisr81 (600, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 1.10212));

    MassPoint mpisr82 (600, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 1.11175));

    MassPoint mpisr83 (600, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 1.10552));

    MassPoint mpisr84 (600, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 1.10447));

    MassPoint mpisr85 (600, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 1.10438));

    MassPoint mpisr86 (600, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 1.10199));

    MassPoint mpisr87 (600, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 1.10425));

    MassPoint mpisr88 (600, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 1.0994));

    MassPoint mpisr89 (600, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 1.10096));

    MassPoint mpisr90 (600, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 1.10817));

    MassPoint mpisr91 (625, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 1.10715));

    MassPoint mpisr92 (625, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 1.10648));

    MassPoint mpisr93 (625, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 1.10532));

    MassPoint mpisr94 (625, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 1.10082));

    MassPoint mpisr95 (625, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 1.10645));

    MassPoint mpisr96 (625, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 1.0998));

    MassPoint mpisr97 (650, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 1.1015));

    MassPoint mpisr98 (650, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 1.11128));

    MassPoint mpisr99 (650, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 1.10242));

    MassPoint mpisr100 (650, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 1.10727));

    MassPoint mpisr101 (650, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 1.09914));

    MassPoint mpisr102 (650, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 1.11084));

    MassPoint mpisr103 (650, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 1.1057));

    MassPoint mpisr104 (650, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 1.11338));

    MassPoint mpisr105 (650, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 1.10254));

    MassPoint mpisr106 (650, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 1.11241));

    MassPoint mpisr107 (650, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 1.10113));

    MassPoint mpisr108 (650, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 1.10837));

    MassPoint mpisr109 (650, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 1.10878));

    MassPoint mpisr110 (650, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 1.09664));

    MassPoint mpisr111 (675, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 1.10065));

    MassPoint mpisr112 (675, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 1.11031));

    MassPoint mpisr113 (675, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 1.10741));

    MassPoint mpisr114 (675, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 1.11915));

    MassPoint mpisr115 (675, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 1.10184));

    MassPoint mpisr116 (675, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 1.10389));

    MassPoint mpisr117 (700, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 1.10023));

    MassPoint mpisr118 (700, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 1.10682));

    MassPoint mpisr119 (700, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 1.11211));

    MassPoint mpisr120 (700, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 1.10465));

    MassPoint mpisr121 (700, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 1.10433));

    MassPoint mpisr122 (700, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 1.10279));

    MassPoint mpisr123 (700, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 1.11108));

    MassPoint mpisr124 (700, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 1.10591));

    MassPoint mpisr125 (700, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 1.10949));

    MassPoint mpisr126 (700, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 1.10493));

    MassPoint mpisr127 (700, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 1.11462));

    MassPoint mpisr128 (700, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 1.10401));

    MassPoint mpisr129 (700, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 1.11676));

    MassPoint mpisr130 (700, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 1.11259));

    MassPoint mpisr131 (700, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 1.11505));

    MassPoint mpisr132 (725, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 1.1283));

    MassPoint mpisr133 (725, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 1.10693));

    MassPoint mpisr134 (725, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 1.10782));

    MassPoint mpisr135 (725, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 1.11533));

    MassPoint mpisr136 (725, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 1.10444));

    MassPoint mpisr137 (725, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 1.10349));

    MassPoint mpisr138 (750, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 1.11104));

    MassPoint mpisr139 (750, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 1.11135));

    MassPoint mpisr140 (750, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 1.11771));

    MassPoint mpisr141 (750, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 1.11833));

    MassPoint mpisr142 (750, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 1.10662));

    MassPoint mpisr143 (750, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 1.10525));

    MassPoint mpisr144 (750, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 1.11096));

    MassPoint mpisr145 (750, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 1.10524));

    MassPoint mpisr146 (750, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 1.11939));

    MassPoint mpisr147 (750, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 1.11093));

    MassPoint mpisr148 (750, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 1.11445));

    MassPoint mpisr149 (750, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 1.1121));

    MassPoint mpisr150 (750, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 1.10568));

    MassPoint mpisr151 (750, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 1.10609));

    MassPoint mpisr152 (750, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 1.10599));

    MassPoint mpisr153 (750, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 1.10598));

    MassPoint mpisr154 (775, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 1.10851));

    MassPoint mpisr155 (775, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 1.11621));

    MassPoint mpisr156 (775, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 1.10974));

    MassPoint mpisr157 (775, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 1.11084));

    MassPoint mpisr158 (775, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 1.10703));

    MassPoint mpisr159 (775, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 1.09569));

    MassPoint mpisr160 (800, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 1.10751));

    MassPoint mpisr161 (800, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 1.13048));

    MassPoint mpisr162 (800, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 1.11047));

    MassPoint mpisr163 (800, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 1.10482));

    MassPoint mpisr164 (800, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 1.11245));

    MassPoint mpisr165 (800, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 1.12227));

    MassPoint mpisr166 (800, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 1.11512));

    MassPoint mpisr167 (800, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 1.11202));

    MassPoint mpisr168 (800, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 1.10707));

    MassPoint mpisr169 (800, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 1.10553));

    MassPoint mpisr170 (800, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 1.11134));

    MassPoint mpisr171 (800, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 1.11618));

    MassPoint mpisr172 (800, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 1.11498));

    MassPoint mpisr173 (800, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 1.11406));

    MassPoint mpisr174 (800, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 1.11213));

    MassPoint mpisr175 (800, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 1.11513));

    MassPoint mpisr176 (800, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 1.10763));

    MassPoint mpisr177 (825, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 1.11562));

    MassPoint mpisr178 (825, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 1.11946));

    MassPoint mpisr179 (825, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 1.11145));

    MassPoint mpisr180 (825, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 1.09832));

    MassPoint mpisr181 (825, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 1.10955));

    MassPoint mpisr182 (825, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 1.10486));

    MassPoint mpisr183 (850, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 1.12687));

    MassPoint mpisr184 (850, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr184, 1.12646));

    MassPoint mpisr185 (850, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr185, 1.10606));

    MassPoint mpisr186 (850, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr186, 1.10663));

    MassPoint mpisr187 (850, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr187, 1.11226));

    MassPoint mpisr188 (850, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr188, 1.11794));

    MassPoint mpisr189 (850, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr189, 1.10477));

    MassPoint mpisr190 (850, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr190, 1.11613));

    MassPoint mpisr191 (850, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr191, 1.11178));

    MassPoint mpisr192 (850, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr192, 1.11974));

    MassPoint mpisr193 (850, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr193, 1.11545));

    MassPoint mpisr194 (850, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr194, 1.1101));

    MassPoint mpisr195 (850, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr195, 1.11822));

    MassPoint mpisr196 (850, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr196, 1.10295));

    MassPoint mpisr197 (850, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr197, 1.11063));

    MassPoint mpisr198 (850, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr198, 1.10976));

    MassPoint mpisr199 (850, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr199, 1.10164));

    MassPoint mpisr200 (875, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr200, 1.11328));

    MassPoint mpisr201 (875, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr201, 1.11539));

    MassPoint mpisr202 (875, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr202, 1.11903));

    MassPoint mpisr203 (875, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr203, 1.10545));

    MassPoint mpisr204 (900, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr204, 1.10812));

    MassPoint mpisr205 (900, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr205, 1.11282));

    MassPoint mpisr206 (900, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr206, 1.10565));

    MassPoint mpisr207 (900, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr207, 1.11757));

    MassPoint mpisr208 (900, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr208, 1.10693));

    MassPoint mpisr209 (900, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr209, 1.11296));

    MassPoint mpisr210 (900, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr210, 1.11575));

    MassPoint mpisr211 (900, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr211, 1.11291));

    MassPoint mpisr212 (900, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr212, 1.12047));

    MassPoint mpisr213 (900, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr213, 1.10508));

    MassPoint mpisr214 (900, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr214, 1.10314));

    MassPoint mpisr215 (900, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr215, 1.10586));

    MassPoint mpisr216 (900, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr216, 1.11914));

    MassPoint mpisr217 (900, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr217, 1.11769));

    MassPoint mpisr218 (900, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr218, 1.10472));

    MassPoint mpisr219 (900, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr219, 1.10832));

    MassPoint mpisr220 (925, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr220, 1.11251));

    MassPoint mpisr221 (925, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr221, 1.10787));

    MassPoint mpisr222 (950, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr222, 1.11151));

    MassPoint mpisr223 (950, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr223, 1.12045));

    MassPoint mpisr224 (950, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr224, 1.11982));

    MassPoint mpisr225 (950, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr225, 1.11178));

    MassPoint mpisr226 (950, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr226, 1.10774));

    MassPoint mpisr227 (950, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr227, 1.11249));

    MassPoint mpisr228 (950, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr228, 1.11641));

    MassPoint mpisr229 (950, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr229, 1.11859));

    MassPoint mpisr230 (950, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr230, 1.10866));

    MassPoint mpisr231 (950, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr231, 1.10442));

    MassPoint mpisr232 (950, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr232, 1.10805));

    MassPoint mpisr233 (950, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr233, 1.12345));

    MassPoint mpisr234 (950, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr234, 1.10288));

    MassPoint mpisr235 (950, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr235, 1.10638));

    MassPoint mpisr236 (950, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr236, 1.11091));

    MassPoint mpisr237 (1000, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr237, 1.12683));

    MassPoint mpisr238 (1000, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr238, 1.11317));

    MassPoint mpisr239 (1000, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr239, 1.11588));

    MassPoint mpisr240 (1000, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr240, 1.12554));

    MassPoint mpisr241 (1000, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr241, 1.11334));

    MassPoint mpisr242 (1000, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr242, 1.12326));

    MassPoint mpisr243 (1000, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr243, 1.11818));

    MassPoint mpisr244 (1000, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr244, 1.1037));

    MassPoint mpisr245 (1000, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr245, 1.11078));

    MassPoint mpisr246 (1000, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr246, 1.1188));

    MassPoint mpisr247 (1000, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr247, 1.11305));

    MassPoint mpisr248 (1000, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr248, 1.10636));

    MassPoint mpisr249 (1000, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr249, 1.10231));

    MassPoint mpisr250 (1000, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr250, 1.11214));

    MassPoint mpisr251 (1000, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr251, 1.11435));

    MassPoint mpisr252 (1050, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr252, 1.11981));

    MassPoint mpisr253 (1050, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr253, 1.11172));

    MassPoint mpisr254 (1050, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr254, 1.11414));

    MassPoint mpisr255 (1050, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr255, 1.11949));

    MassPoint mpisr256 (1050, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr256, 1.11711));

    MassPoint mpisr257 (1050, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr257, 1.11621));

    MassPoint mpisr258 (1050, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr258, 1.11607));

    MassPoint mpisr259 (1050, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr259, 1.1081));

    MassPoint mpisr260 (1050, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr260, 1.12002));

    MassPoint mpisr261 (1050, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr261, 1.11774));

    MassPoint mpisr262 (1050, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr262, 1.12748));

    MassPoint mpisr263 (1050, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr263, 1.1134));

    MassPoint mpisr264 (1050, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr264, 1.10302));

    MassPoint mpisr265 (1050, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr265, 1.11142));

    MassPoint mpisr266 (1050, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr266, 1.12287));

    MassPoint mpisr267 (1100, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr267, 1.11723));

    MassPoint mpisr268 (1100, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr268, 1.11225));

    MassPoint mpisr269 (1100, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr269, 1.10703));

    MassPoint mpisr270 (1100, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr270, 1.1141));

    MassPoint mpisr271 (1100, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr271, 1.10852));

    MassPoint mpisr272 (1100, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr272, 1.12228));

    MassPoint mpisr273 (1100, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr273, 1.11737));

    MassPoint mpisr274 (1100, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr274, 1.1148));

    MassPoint mpisr275 (1100, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr275, 1.11497));

    MassPoint mpisr276 (1100, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr276, 1.11278));

    MassPoint mpisr277 (1100, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr277, 1.10245));

    MassPoint mpisr278 (1100, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr278, 1.11723));

    MassPoint mpisr279 (1100, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr279, 1.1123));

    MassPoint mpisr280 (1100, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr280, 1.12391));

    MassPoint mpisr281 (1100, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr281, 1.10755));

    MassPoint mpisr282 (1150, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr282, 1.11707));

    MassPoint mpisr283 (1150, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr283, 1.1163));

    MassPoint mpisr284 (1150, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr284, 1.1306));

    MassPoint mpisr285 (1150, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr285, 1.11461));

    MassPoint mpisr286 (1150, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr286, 1.10805));

    MassPoint mpisr287 (1150, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr287, 1.12104));

    MassPoint mpisr288 (1150, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr288, 1.11341));

    MassPoint mpisr289 (1150, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr289, 1.10857));

    MassPoint mpisr290 (1150, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr290, 1.11329));

    MassPoint mpisr291 (1150, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr291, 1.12935));

    MassPoint mpisr292 (1150, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr292, 1.11576));

    MassPoint mpisr293 (1150, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr293, 1.13024));

    MassPoint mpisr294 (1150, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr294, 1.10307));

    MassPoint mpisr295 (1150, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr295, 1.11961));

    MassPoint mpisr296 (1150, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr296, 1.12111));

    MassPoint mpisr297 (1200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr297, 1.11539));

    MassPoint mpisr298 (1200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr298, 1.11261));

    MassPoint mpisr299 (1200, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr299, 1.11725));

    MassPoint mpisr300 (1200, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr300, 1.12438));

    MassPoint mpisr301 (1200, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr301, 1.11309));

    MassPoint mpisr302 (1200, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr302, 1.11823));

    MassPoint mpisr303 (1200, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr303, 1.11205));

    MassPoint mpisr304 (1200, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr304, 1.11474));

    MassPoint mpisr305 (1200, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr305, 1.11198));

    MassPoint mpisr306 (1200, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr306, 1.10503));

    MassPoint mpisr307 (1200, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr307, 1.12286));

    MassPoint mpisr308 (1200, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr308, 1.11219));

    MassPoint mpisr309 (1200, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr309, 1.11328));

    MassPoint mpisr310 (1200, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr310, 1.11423));

    MassPoint mpisr311 (1200, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr311, 1.11246));

  } else if (SUSYProductionProcess=="T2bt") {

     MassPoint mp0 (200, 1);
  StopCrossSection cs0 (64.5085, 9.29555);
  MassPointParameters mpp0 (cs0, 1248580);
  StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

  MassPoint mp1 (225, 25);
  StopCrossSection cs1 (36.3818, 5.17309);
  MassPointParameters mpp1 (cs1, 692125);
  StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

  MassPoint mp2 (250, 1);
  StopCrossSection cs2 (21.5949, 3.03613);
  MassPointParameters mpp2 (cs2, 434391);
  StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

  MassPoint mp3 (250, 25);
  StopCrossSection cs3 (21.5949, 3.03613);
  MassPointParameters mpp3 (cs3, 451646);
  StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

  MassPoint mp4 (250, 50);
  StopCrossSection cs4 (21.5949, 3.03613);
  MassPointParameters mpp4 (cs4, 444437);
  StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

  MassPoint mp5 (275, 75);
  StopCrossSection cs5 (13.3231, 1.89919);
  MassPointParameters mpp5 (cs5, 438466);
  StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

  MassPoint mp6 (300, 1);
  StopCrossSection cs6 (8.51615, 1.18564);
  MassPointParameters mpp6 (cs6, 479507);
  StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

  MassPoint mp7 (300, 25);
  StopCrossSection cs7 (8.51615, 1.18564);
  MassPointParameters mpp7 (cs7, 491569);
  StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

  MassPoint mp8 (300, 50);
  StopCrossSection cs8 (8.51615, 1.18564);
  MassPointParameters mpp8 (cs8, 489655);
  StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

  MassPoint mp9 (300, 75);
  StopCrossSection cs9 (8.51615, 1.18564);
  MassPointParameters mpp9 (cs9, 491821);
  StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

  MassPoint mp10 (300, 100);
  StopCrossSection cs10 (8.51615, 1.18564);
  MassPointParameters mpp10 (cs10, 504788);
  StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

  MassPoint mp11 (325, 75);
  StopCrossSection cs11 (5.60471, 0.774257);
  MassPointParameters mpp11 (cs11, 473601);
  StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

  MassPoint mp12 (325, 100);
  StopCrossSection cs12 (5.60471, 0.774257);
  MassPointParameters mpp12 (cs12, 472959);
  StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

  MassPoint mp13 (325, 125);
  StopCrossSection cs13 (5.60471, 0.774257);
  MassPointParameters mpp13 (cs13, 474782);
  StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

  MassPoint mp14 (350, 1);
  StopCrossSection cs14 (3.78661, 0.5183);
  MassPointParameters mpp14 (cs14, 463542);
  StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

  MassPoint mp15 (350, 25);
  StopCrossSection cs15 (3.78661, 0.5183);
  MassPointParameters mpp15 (cs15, 467417);
  StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

  MassPoint mp16 (350, 50);
  StopCrossSection cs16 (3.78661, 0.5183);
  MassPointParameters mpp16 (cs16, 467109);
  StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

  MassPoint mp17 (350, 75);
  StopCrossSection cs17 (3.78661, 0.5183);
  MassPointParameters mpp17 (cs17, 474349);
  StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

  MassPoint mp18 (350, 100);
  StopCrossSection cs18 (3.78661, 0.5183);
  MassPointParameters mpp18 (cs18, 456313);
  StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

  MassPoint mp19 (350, 125);
  StopCrossSection cs19 (3.78661, 0.5183);
  MassPointParameters mpp19 (cs19, 451862);
  StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

  MassPoint mp20 (350, 150);
  StopCrossSection cs20 (3.78661, 0.5183);
  MassPointParameters mpp20 (cs20, 466737);
  StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

  MassPoint mp21 (375, 75);
  StopCrossSection cs21 (2.61162, 0.361649);
  MassPointParameters mpp21 (cs21, 446309);
  StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

  MassPoint mp22 (375, 100);
  StopCrossSection cs22 (2.61162, 0.361649);
  MassPointParameters mpp22 (cs22, 445491);
  StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

  MassPoint mp23 (375, 125);
  StopCrossSection cs23 (2.61162, 0.361649);
  MassPointParameters mpp23 (cs23, 443873);
  StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

  MassPoint mp24 (375, 150);
  StopCrossSection cs24 (2.61162, 0.361649);
  MassPointParameters mpp24 (cs24, 454022);
  StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

  MassPoint mp25 (375, 175);
  StopCrossSection cs25 (2.61162, 0.361649);
  MassPointParameters mpp25 (cs25, 472272);
  StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

  MassPoint mp26 (400, 1);
  StopCrossSection cs26 (1.83537, 0.251418);
  MassPointParameters mpp26 (cs26, 467276);
  StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

  MassPoint mp27 (400, 25);
  StopCrossSection cs27 (1.83537, 0.251418);
  MassPointParameters mpp27 (cs27, 470207);
  StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

  MassPoint mp28 (400, 50);
  StopCrossSection cs28 (1.83537, 0.251418);
  MassPointParameters mpp28 (cs28, 481250);
  StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

  MassPoint mp29 (400, 100);
  StopCrossSection cs29 (1.83537, 0.251418);
  MassPointParameters mpp29 (cs29, 487553);
  StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

  MassPoint mp30 (400, 125);
  StopCrossSection cs30 (1.83537, 0.251418);
  MassPointParameters mpp30 (cs30, 489058);
  StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

  MassPoint mp31 (400, 150);
  StopCrossSection cs31 (1.83537, 0.251418);
  MassPointParameters mpp31 (cs31, 488799);
  StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

  MassPoint mp32 (400, 175);
  StopCrossSection cs32 (1.83537, 0.251418);
  MassPointParameters mpp32 (cs32, 482168);
  StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

  MassPoint mp33 (400, 200);
  StopCrossSection cs33 (1.83537, 0.251418);
  MassPointParameters mpp33 (cs33, 487766);
  StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

  MassPoint mp34 (425, 125);
  StopCrossSection cs34 (1.31169, 0.177095);
  MassPointParameters mpp34 (cs34, 478131);
  StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

  MassPoint mp35 (425, 150);
  StopCrossSection cs35 (1.31169, 0.177095);
  MassPointParameters mpp35 (cs35, 480046);
  StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

  MassPoint mp36 (425, 175);
  StopCrossSection cs36 (1.31169, 0.177095);
  MassPointParameters mpp36 (cs36, 455404);
  StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

  MassPoint mp37 (425, 200);
  StopCrossSection cs37 (1.31169, 0.177095);
  MassPointParameters mpp37 (cs37, 466226);
  StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

  MassPoint mp38 (425, 225);
  StopCrossSection cs38 (1.31169, 0.177095);
  MassPointParameters mpp38 (cs38, 476639);
  StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

  MassPoint mp39 (450, 1);
  StopCrossSection cs39 (0.948333, 0.127607);
  MassPointParameters mpp39 (cs39, 348296);
  StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

  MassPoint mp40 (450, 25);
  StopCrossSection cs40 (0.948333, 0.127607);
  MassPointParameters mpp40 (cs40, 347277);
  StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

  MassPoint mp41 (450, 50);
  StopCrossSection cs41 (0.948333, 0.127607);
  MassPointParameters mpp41 (cs41, 353544);
  StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

  MassPoint mp42 (450, 100);
  StopCrossSection cs42 (0.948333, 0.127607);
  MassPointParameters mpp42 (cs42, 360239);
  StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

  MassPoint mp43 (450, 150);
  StopCrossSection cs43 (0.948333, 0.127607);
  MassPointParameters mpp43 (cs43, 352870);
  StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

  MassPoint mp44 (450, 175);
  StopCrossSection cs44 (0.948333, 0.127607);
  MassPointParameters mpp44 (cs44, 348622);
  StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

  MassPoint mp45 (450, 200);
  StopCrossSection cs45 (0.948333, 0.127607);
  MassPointParameters mpp45 (cs45, 356432);
  StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

  MassPoint mp46 (450, 225);
  StopCrossSection cs46 (0.948333, 0.127607);
  MassPointParameters mpp46 (cs46, 363885);
  StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

  MassPoint mp47 (450, 250);
  StopCrossSection cs47 (0.948333, 0.127607);
  MassPointParameters mpp47 (cs47, 352189);
  StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

  MassPoint mp48 (475, 175);
  StopCrossSection cs48 (0.697075, 0.0933565);
  MassPointParameters mpp48 (cs48, 264823);
  StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

  MassPoint mp49 (475, 200);
  StopCrossSection cs49 (0.697075, 0.0933565);
  MassPointParameters mpp49 (cs49, 253858);
  StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

  MassPoint mp50 (475, 225);
  StopCrossSection cs50 (0.697075, 0.0933565);
  MassPointParameters mpp50 (cs50, 253433);
  StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

  MassPoint mp51 (475, 250);
  StopCrossSection cs51 (0.697075, 0.0933565);
  MassPointParameters mpp51 (cs51, 264092);
  StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

  MassPoint mp52 (475, 275);
  StopCrossSection cs52 (0.697075, 0.0933565);
  MassPointParameters mpp52 (cs52, 249571);
  StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

  MassPoint mp53 (500, 1);
  StopCrossSection cs53 (0.51848, 0.0693711);
  MassPointParameters mpp53 (cs53, 213165);
  StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

  MassPoint mp54 (500, 25);
  StopCrossSection cs54 (0.51848, 0.0693711);
  MassPointParameters mpp54 (cs54, 206314);
  StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

  MassPoint mp55 (500, 50);
  StopCrossSection cs55 (0.51848, 0.0693711);
  MassPointParameters mpp55 (cs55, 206920);
  StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

  MassPoint mp56 (500, 100);
  StopCrossSection cs56 (0.51848, 0.0693711);
  MassPointParameters mpp56 (cs56, 196813);
  StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

  MassPoint mp57 (500, 150);
  StopCrossSection cs57 (0.51848, 0.0693711);
  MassPointParameters mpp57 (cs57, 202245);
  StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

  MassPoint mp58 (500, 200);
  StopCrossSection cs58 (0.51848, 0.0693711);
  MassPointParameters mpp58 (cs58, 199928);
  StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

  MassPoint mp59 (500, 225);
  StopCrossSection cs59 (0.51848, 0.0693711);
  MassPointParameters mpp59 (cs59, 195953);
  StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

  MassPoint mp60 (500, 250);
  StopCrossSection cs60 (0.51848, 0.0693711);
  MassPointParameters mpp60 (cs60, 196844);
  StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

  MassPoint mp61 (500, 275);
  StopCrossSection cs61 (0.51848, 0.0693711);
  MassPointParameters mpp61 (cs61, 196339);
  StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

  MassPoint mp62 (500, 300);
  StopCrossSection cs62 (0.51848, 0.0693711);
  MassPointParameters mpp62 (cs62, 208625);
  StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

  MassPoint mp63 (525, 225);
  StopCrossSection cs63 (0.390303, 0.0520832);
  MassPointParameters mpp63 (cs63, 150541);
  StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

  MassPoint mp64 (525, 250);
  StopCrossSection cs64 (0.390303, 0.0520832);
  MassPointParameters mpp64 (cs64, 143785);
  StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

  MassPoint mp65 (525, 275);
  StopCrossSection cs65 (0.390303, 0.0520832);
  MassPointParameters mpp65 (cs65, 147121);
  StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

  MassPoint mp66 (525, 300);
  StopCrossSection cs66 (0.390303, 0.0520832);
  MassPointParameters mpp66 (cs66, 141636);
  StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

  MassPoint mp67 (525, 325);
  StopCrossSection cs67 (0.390303, 0.0520832);
  MassPointParameters mpp67 (cs67, 147778);
  StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

  MassPoint mp68 (550, 1);
  StopCrossSection cs68 (0.296128, 0.0392923);
  MassPointParameters mpp68 (cs68, 115656);
  StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

  MassPoint mp69 (550, 25);
  StopCrossSection cs69 (0.296128, 0.0392923);
  MassPointParameters mpp69 (cs69, 110772);
  StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

  MassPoint mp70 (550, 50);
  StopCrossSection cs70 (0.296128, 0.0392923);
  MassPointParameters mpp70 (cs70, 112813);
  StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

  MassPoint mp71 (550, 100);
  StopCrossSection cs71 (0.296128, 0.0392923);
  MassPointParameters mpp71 (cs71, 106423);
  StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

  MassPoint mp72 (550, 150);
  StopCrossSection cs72 (0.296128, 0.0392923);
  MassPointParameters mpp72 (cs72, 120964);
  StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

  MassPoint mp73 (550, 200);
  StopCrossSection cs73 (0.296128, 0.0392923);
  MassPointParameters mpp73 (cs73, 103880);
  StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

  MassPoint mp74 (550, 250);
  StopCrossSection cs74 (0.296128, 0.0392923);
  MassPointParameters mpp74 (cs74, 115835);
  StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

  MassPoint mp75 (550, 275);
  StopCrossSection cs75 (0.296128, 0.0392923);
  MassPointParameters mpp75 (cs75, 109044);
  StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

  MassPoint mp76 (550, 300);
  StopCrossSection cs76 (0.296128, 0.0392923);
  MassPointParameters mpp76 (cs76, 108508);
  StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

  MassPoint mp77 (550, 325);
  StopCrossSection cs77 (0.296128, 0.0392923);
  MassPointParameters mpp77 (cs77, 115623);
  StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

  MassPoint mp78 (550, 350);
  StopCrossSection cs78 (0.296128, 0.0392923);
  MassPointParameters mpp78 (cs78, 110447);
  StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

  MassPoint mp79 (575, 275);
  StopCrossSection cs79 (0.226118, 0.0300151);
  MassPointParameters mpp79 (cs79, 86721);
  StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

  MassPoint mp80 (575, 300);
  StopCrossSection cs80 (0.226118, 0.0300151);
  MassPointParameters mpp80 (cs80, 86492);
  StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

  MassPoint mp81 (575, 325);
  StopCrossSection cs81 (0.226118, 0.0300151);
  MassPointParameters mpp81 (cs81, 89492);
  StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

  MassPoint mp82 (575, 350);
  StopCrossSection cs82 (0.226118, 0.0300151);
  MassPointParameters mpp82 (cs82, 90411);
  StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

  MassPoint mp83 (575, 375);
  StopCrossSection cs83 (0.226118, 0.0300151);
  MassPointParameters mpp83 (cs83, 83086);
  StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

  MassPoint mp84 (600, 1);
  StopCrossSection cs84 (0.174599, 0.02306);
  MassPointParameters mpp84 (cs84, 70857);
  StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

  MassPoint mp85 (600, 25);
  StopCrossSection cs85 (0.174599, 0.02306);
  MassPointParameters mpp85 (cs85, 74502);
  StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

  MassPoint mp86 (600, 50);
  StopCrossSection cs86 (0.174599, 0.02306);
  MassPointParameters mpp86 (cs86, 67479);
  StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

  MassPoint mp87 (600, 100);
  StopCrossSection cs87 (0.174599, 0.02306);
  MassPointParameters mpp87 (cs87, 71085);
  StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

  MassPoint mp88 (600, 150);
  StopCrossSection cs88 (0.174599, 0.02306);
  MassPointParameters mpp88 (cs88, 73574);
  StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

  MassPoint mp89 (600, 200);
  StopCrossSection cs89 (0.174599, 0.02306);
  MassPointParameters mpp89 (cs89, 70645);
  StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

  MassPoint mp90 (600, 250);
  StopCrossSection cs90 (0.174599, 0.02306);
  MassPointParameters mpp90 (cs90, 74241);
  StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

  MassPoint mp91 (600, 300);
  StopCrossSection cs91 (0.174599, 0.02306);
  MassPointParameters mpp91 (cs91, 65374);
  StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

  MassPoint mp92 (600, 325);
  StopCrossSection cs92 (0.174599, 0.02306);
  MassPointParameters mpp92 (cs92, 77056);
  StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

  MassPoint mp93 (600, 350);
  StopCrossSection cs93 (0.174599, 0.02306);
  MassPointParameters mpp93 (cs93, 67981);
  StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

  MassPoint mp94 (600, 375);
  StopCrossSection cs94 (0.174599, 0.02306);
  MassPointParameters mpp94 (cs94, 69909);
  StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

  MassPoint mp95 (600, 400);
  StopCrossSection cs95 (0.174599, 0.02306);
  MassPointParameters mpp95 (cs95, 69659);
  StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

  MassPoint mp96 (625, 325);
  StopCrossSection cs96 (0.136372, 0.0174504);
  MassPointParameters mpp96 (cs96, 50316);
  StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

  MassPoint mp97 (625, 350);
  StopCrossSection cs97 (0.136372, 0.0174504);
  MassPointParameters mpp97 (cs97, 56045);
  StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

  MassPoint mp98 (625, 375);
  StopCrossSection cs98 (0.136372, 0.0174504);
  MassPointParameters mpp98 (cs98, 57829);
  StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

  MassPoint mp99 (625, 400);
  StopCrossSection cs99 (0.136372, 0.0174504);
  MassPointParameters mpp99 (cs99, 55234);
  StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

  MassPoint mp100 (625, 425);
  StopCrossSection cs100 (0.136372, 0.0174504);
  MassPointParameters mpp100 (cs100, 56227);
  StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

  MassPoint mp101 (650, 1);
  StopCrossSection cs101 (0.107045, 0.0138336);
  MassPointParameters mpp101 (cs101, 43111);
  StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

  MassPoint mp102 (650, 25);
  StopCrossSection cs102 (0.107045, 0.0138336);
  MassPointParameters mpp102 (cs102, 36459);
  StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

  MassPoint mp103 (650, 50);
  StopCrossSection cs103 (0.107045, 0.0138336);
  MassPointParameters mpp103 (cs103, 45457);
  StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

  MassPoint mp104 (650, 100);
  StopCrossSection cs104 (0.107045, 0.0138336);
  MassPointParameters mpp104 (cs104, 39817);
  StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

  MassPoint mp105 (650, 150);
  StopCrossSection cs105 (0.107045, 0.0138336);
  MassPointParameters mpp105 (cs105, 44410);
  StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

  MassPoint mp106 (650, 200);
  StopCrossSection cs106 (0.107045, 0.0138336);
  MassPointParameters mpp106 (cs106, 47166);
  StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

  MassPoint mp107 (650, 250);
  StopCrossSection cs107 (0.107045, 0.0138336);
  MassPointParameters mpp107 (cs107, 44391);
  StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

  MassPoint mp108 (650, 300);
  StopCrossSection cs108 (0.107045, 0.0138336);
  MassPointParameters mpp108 (cs108, 48558);
  StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

  MassPoint mp109 (650, 350);
  StopCrossSection cs109 (0.107045, 0.0138336);
  MassPointParameters mpp109 (cs109, 41259);
  StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

  MassPoint mp110 (650, 375);
  StopCrossSection cs110 (0.107045, 0.0138336);
  MassPointParameters mpp110 (cs110, 39497);
  StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

  MassPoint mp111 (650, 400);
  StopCrossSection cs111 (0.107045, 0.0138336);
  MassPointParameters mpp111 (cs111, 45281);
  StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

  MassPoint mp112 (650, 425);
  StopCrossSection cs112 (0.107045, 0.0138336);
  MassPointParameters mpp112 (cs112, 41907);
  StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

  MassPoint mp113 (650, 450);
  StopCrossSection cs113 (0.107045, 0.0138336);
  MassPointParameters mpp113 (cs113, 47530);
  StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

  MassPoint mp114 (675, 375);
  StopCrossSection cs114 (0.0844877, 0.0110428);
  MassPointParameters mpp114 (cs114, 36856);
  StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

  MassPoint mp115 (675, 400);
  StopCrossSection cs115 (0.0844877, 0.0110428);
  MassPointParameters mpp115 (cs115, 30903);
  StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

  MassPoint mp116 (675, 425);
  StopCrossSection cs116 (0.0844877, 0.0110428);
  MassPointParameters mpp116 (cs116, 31436);
  StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

  MassPoint mp117 (675, 450);
  StopCrossSection cs117 (0.0844877, 0.0110428);
  MassPointParameters mpp117 (cs117, 28420);
  StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

  MassPoint mp118 (675, 475);
  StopCrossSection cs118 (0.0844877, 0.0110428);
  MassPointParameters mpp118 (cs118, 34934);
  StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

  MassPoint mp119 (700, 1);
  StopCrossSection cs119 (0.0670476, 0.00894609);
  MassPointParameters mpp119 (cs119, 28579);
  StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

  MassPoint mp120 (700, 25);
  StopCrossSection cs120 (0.0670476, 0.00894609);
  MassPointParameters mpp120 (cs120, 29813);
  StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

  MassPoint mp121 (700, 50);
  StopCrossSection cs121 (0.0670476, 0.00894609);
  MassPointParameters mpp121 (cs121, 26450);
  StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

  MassPoint mp122 (700, 100);
  StopCrossSection cs122 (0.0670476, 0.00894609);
  MassPointParameters mpp122 (cs122, 24337);
  StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

  MassPoint mp123 (700, 150);
  StopCrossSection cs123 (0.0670476, 0.00894609);
  MassPointParameters mpp123 (cs123, 27818);
  StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

  MassPoint mp124 (700, 200);
  StopCrossSection cs124 (0.0670476, 0.00894609);
  MassPointParameters mpp124 (cs124, 26811);
  StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

  MassPoint mp125 (700, 250);
  StopCrossSection cs125 (0.0670476, 0.00894609);
  MassPointParameters mpp125 (cs125, 29133);
  StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

  MassPoint mp126 (700, 300);
  StopCrossSection cs126 (0.0670476, 0.00894609);
  MassPointParameters mpp126 (cs126, 27832);
  StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

  MassPoint mp127 (700, 350);
  StopCrossSection cs127 (0.0670476, 0.00894609);
  MassPointParameters mpp127 (cs127, 27677);
  StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

  MassPoint mp128 (700, 400);
  StopCrossSection cs128 (0.0670476, 0.00894609);
  MassPointParameters mpp128 (cs128, 23446);
  StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

  MassPoint mp129 (700, 425);
  StopCrossSection cs129 (0.0670476, 0.00894609);
  MassPointParameters mpp129 (cs129, 28963);
  StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

  MassPoint mp130 (700, 450);
  StopCrossSection cs130 (0.0670476, 0.00894609);
  MassPointParameters mpp130 (cs130, 24059);
  StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

  MassPoint mp131 (700, 475);
  StopCrossSection cs131 (0.0670476, 0.00894609);
  MassPointParameters mpp131 (cs131, 21220);
  StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

  MassPoint mp132 (700, 500);
  StopCrossSection cs132 (0.0670476, 0.00894609);
  MassPointParameters mpp132 (cs132, 24525);
  StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

  MassPoint mp133 (725, 425);
  StopCrossSection cs133 (0.0536438, 0.00728504);
  MassPointParameters mpp133 (cs133, 18241);
  StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

  MassPoint mp134 (725, 450);
  StopCrossSection cs134 (0.0536438, 0.00728504);
  MassPointParameters mpp134 (cs134, 20691);
  StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

  MassPoint mp135 (725, 475);
  StopCrossSection cs135 (0.0536438, 0.00728504);
  MassPointParameters mpp135 (cs135, 23664);
  StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

  MassPoint mp136 (725, 500);
  StopCrossSection cs136 (0.0536438, 0.00728504);
  MassPointParameters mpp136 (cs136, 17125);
  StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

  MassPoint mp137 (725, 525);
  StopCrossSection cs137 (0.0536438, 0.00728504);
  MassPointParameters mpp137 (cs137, 22286);
  StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

  MassPoint mp138 (750, 1);
  StopCrossSection cs138 (0.0431418, 0.00593006);
  MassPointParameters mpp138 (cs138, 17531);
  StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

  MassPoint mp139 (750, 25);
  StopCrossSection cs139 (0.0431418, 0.00593006);
  MassPointParameters mpp139 (cs139, 15989);
  StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

  MassPoint mp140 (750, 50);
  StopCrossSection cs140 (0.0431418, 0.00593006);
  MassPointParameters mpp140 (cs140, 18866);
  StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

  MassPoint mp141 (750, 100);
  StopCrossSection cs141 (0.0431418, 0.00593006);
  MassPointParameters mpp141 (cs141, 19025);
  StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

  MassPoint mp142 (750, 150);
  StopCrossSection cs142 (0.0431418, 0.00593006);
  MassPointParameters mpp142 (cs142, 17424);
  StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

  MassPoint mp143 (750, 200);
  StopCrossSection cs143 (0.0431418, 0.00593006);
  MassPointParameters mpp143 (cs143, 16025);
  StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

  MassPoint mp144 (750, 250);
  StopCrossSection cs144 (0.0431418, 0.00593006);
  MassPointParameters mpp144 (cs144, 20643);
  StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

  MassPoint mp145 (750, 300);
  StopCrossSection cs145 (0.0431418, 0.00593006);
  MassPointParameters mpp145 (cs145, 17964);
  StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

  MassPoint mp146 (750, 350);
  StopCrossSection cs146 (0.0431418, 0.00593006);
  MassPointParameters mpp146 (cs146, 19301);
  StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

  MassPoint mp147 (750, 400);
  StopCrossSection cs147 (0.0431418, 0.00593006);
  MassPointParameters mpp147 (cs147, 21827);
  StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

  MassPoint mp148 (750, 450);
  StopCrossSection cs148 (0.0431418, 0.00593006);
  MassPointParameters mpp148 (cs148, 17734);
  StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

  MassPoint mp149 (750, 475);
  StopCrossSection cs149 (0.0431418, 0.00593006);
  MassPointParameters mpp149 (cs149, 20241);
  StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

  MassPoint mp150 (750, 500);
  StopCrossSection cs150 (0.0431418, 0.00593006);
  MassPointParameters mpp150 (cs150, 19566);
  StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

  MassPoint mp151 (750, 525);
  StopCrossSection cs151 (0.0431418, 0.00593006);
  MassPointParameters mpp151 (cs151, 22390);
  StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

  MassPoint mp152 (750, 550);
  StopCrossSection cs152 (0.0431418, 0.00593006);
  MassPointParameters mpp152 (cs152, 20620);
  StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

  MassPoint mp153 (775, 475);
  StopCrossSection cs153 (0.0348796, 0.00486909);
  MassPointParameters mpp153 (cs153, 21293);
  StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

  MassPoint mp154 (775, 500);
  StopCrossSection cs154 (0.0348796, 0.00486909);
  MassPointParameters mpp154 (cs154, 19164);
  StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

  MassPoint mp155 (775, 525);
  StopCrossSection cs155 (0.0348796, 0.00486909);
  MassPointParameters mpp155 (cs155, 21173);
  StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

  MassPoint mp156 (775, 550);
  StopCrossSection cs156 (0.0348796, 0.00486909);
  MassPointParameters mpp156 (cs156, 17186);
  StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

  MassPoint mp157 (775, 575);
  StopCrossSection cs157 (0.0348796, 0.00486909);
  MassPointParameters mpp157 (cs157, 17992);
  StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

  MassPoint mp158 (800, 1);
  StopCrossSection cs158 (0.0283338, 0.00401518);
  MassPointParameters mpp158 (cs158, 19843);
  StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

  MassPoint mp159 (800, 25);
  StopCrossSection cs159 (0.0283338, 0.00401518);
  MassPointParameters mpp159 (cs159, 17185);
  StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

  MassPoint mp160 (800, 50);
  StopCrossSection cs160 (0.0283338, 0.00401518);
  MassPointParameters mpp160 (cs160, 20162);
  StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

  MassPoint mp161 (800, 100);
  StopCrossSection cs161 (0.0283338, 0.00401518);
  MassPointParameters mpp161 (cs161, 20534);
  StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

  MassPoint mp162 (800, 150);
  StopCrossSection cs162 (0.0283338, 0.00401518);
  MassPointParameters mpp162 (cs162, 18402);
  StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

  MassPoint mp163 (800, 200);
  StopCrossSection cs163 (0.0283338, 0.00401518);
  MassPointParameters mpp163 (cs163, 18425);
  StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

  MassPoint mp164 (800, 250);
  StopCrossSection cs164 (0.0283338, 0.00401518);
  MassPointParameters mpp164 (cs164, 23249);
  StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

  MassPoint mp165 (800, 300);
  StopCrossSection cs165 (0.0283338, 0.00401518);
  MassPointParameters mpp165 (cs165, 18036);
  StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

  MassPoint mp166 (800, 350);
  StopCrossSection cs166 (0.0283338, 0.00401518);
  MassPointParameters mpp166 (cs166, 21135);
  StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

  MassPoint mp167 (800, 400);
  StopCrossSection cs167 (0.0283338, 0.00401518);
  MassPointParameters mpp167 (cs167, 16452);
  StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

  MassPoint mp168 (800, 450);
  StopCrossSection cs168 (0.0283338, 0.00401518);
  MassPointParameters mpp168 (cs168, 20477);
  StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

  MassPoint mp169 (800, 500);
  StopCrossSection cs169 (0.0283338, 0.00401518);
  MassPointParameters mpp169 (cs169, 19287);
  StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

  MassPoint mp170 (800, 525);
  StopCrossSection cs170 (0.0283338, 0.00401518);
  MassPointParameters mpp170 (cs170, 18729);
  StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

  MassPoint mp171 (800, 550);
  StopCrossSection cs171 (0.0283338, 0.00401518);
  MassPointParameters mpp171 (cs171, 18576);
  StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

  MassPoint mp172 (800, 575);
  StopCrossSection cs172 (0.0283338, 0.00401518);
  MassPointParameters mpp172 (cs172, 16614);
  StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

  MassPoint mp173 (800, 600);
  StopCrossSection cs173 (0.0283338, 0.00401518);
  MassPointParameters mpp173 (cs173, 20812);
  StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

  MassPoint mp174 (825, 525);
  StopCrossSection cs174 (0.0230866, 0.00333435);
  MassPointParameters mpp174 (cs174, 18806);
  StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

  MassPoint mp175 (825, 550);
  StopCrossSection cs175 (0.0230866, 0.00333435);
  MassPointParameters mpp175 (cs175, 17147);
  StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

  MassPoint mp176 (825, 575);
  StopCrossSection cs176 (0.0230866, 0.00333435);
  MassPointParameters mpp176 (cs176, 16702);
  StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

  MassPoint mp177 (825, 600);
  StopCrossSection cs177 (0.0230866, 0.00333435);
  MassPointParameters mpp177 (cs177, 20559);
  StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

  MassPoint mp178 (825, 625);
  StopCrossSection cs178 (0.0230866, 0.00333435);
  MassPointParameters mpp178 (cs178, 18118);
  StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

  MassPoint mp179 (850, 1);
  StopCrossSection cs179 (0.0189612, 0.00278768);
  MassPointParameters mpp179 (cs179, 18604);
  StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

  MassPoint mp180 (850, 25);
  StopCrossSection cs180 (0.0189612, 0.00278768);
  MassPointParameters mpp180 (cs180, 19759);
  StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

  MassPoint mp181 (850, 50);
  StopCrossSection cs181 (0.0189612, 0.00278768);
  MassPointParameters mpp181 (cs181, 17194);
  StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

  MassPoint mp182 (850, 100);
  StopCrossSection cs182 (0.0189612, 0.00278768);
  MassPointParameters mpp182 (cs182, 21961);
  StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

  MassPoint mp183 (850, 150);
  StopCrossSection cs183 (0.0189612, 0.00278768);
  MassPointParameters mpp183 (cs183, 19253);
  StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

  MassPoint mp184 (850, 200);
  StopCrossSection cs184 (0.0189612, 0.00278768);
  MassPointParameters mpp184 (cs184, 19385);
  StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

  MassPoint mp185 (850, 250);
  StopCrossSection cs185 (0.0189612, 0.00278768);
  MassPointParameters mpp185 (cs185, 15605);
  StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

  MassPoint mp186 (850, 300);
  StopCrossSection cs186 (0.0189612, 0.00278768);
  MassPointParameters mpp186 (cs186, 17491);
  StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

  MassPoint mp187 (850, 350);
  StopCrossSection cs187 (0.0189612, 0.00278768);
  MassPointParameters mpp187 (cs187, 19659);
  StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

  MassPoint mp188 (850, 400);
  StopCrossSection cs188 (0.0189612, 0.00278768);
  MassPointParameters mpp188 (cs188, 22296);
  StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

  MassPoint mp189 (850, 450);
  StopCrossSection cs189 (0.0189612, 0.00278768);
  MassPointParameters mpp189 (cs189, 19968);
  StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

  MassPoint mp190 (850, 500);
  StopCrossSection cs190 (0.0189612, 0.00278768);
  MassPointParameters mpp190 (cs190, 16118);
  StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

  MassPoint mp191 (850, 550);
  StopCrossSection cs191 (0.0189612, 0.00278768);
  MassPointParameters mpp191 (cs191, 17261);
  StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

  MassPoint mp192 (850, 575);
  StopCrossSection cs192 (0.0189612, 0.00278768);
  MassPointParameters mpp192 (cs192, 19260);
  StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

  MassPoint mp193 (850, 600);
  StopCrossSection cs193 (0.0189612, 0.00278768);
  MassPointParameters mpp193 (cs193, 20181);
  StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

  MassPoint mp194 (850, 625);
  StopCrossSection cs194 (0.0189612, 0.00278768);
  MassPointParameters mpp194 (cs194, 15327);
  StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

  MassPoint mp195 (850, 650);
  StopCrossSection cs195 (0.0189612, 0.00278768);
  MassPointParameters mpp195 (cs195, 22532);
  StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

  MassPoint mp196 (875, 575);
  StopCrossSection cs196 (0.015625, 0.00233698);
  MassPointParameters mpp196 (cs196, 19785);
  StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

  MassPoint mp197 (875, 600);
  StopCrossSection cs197 (0.015625, 0.00233698);
  MassPointParameters mpp197 (cs197, 20076);
  StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

  MassPoint mp198 (875, 625);
  StopCrossSection cs198 (0.015625, 0.00233698);
  MassPointParameters mpp198 (cs198, 18636);
  StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

  MassPoint mp199 (875, 650);
  StopCrossSection cs199 (0.015625, 0.00233698);
  MassPointParameters mpp199 (cs199, 20249);
  StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

  MassPoint mp200 (900, 1);
  StopCrossSection cs200 (0.0128895, 0.00195954);
  MassPointParameters mpp200 (cs200, 20534);
  StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

  MassPoint mp201 (900, 25);
  StopCrossSection cs201 (0.0128895, 0.00195954);
  MassPointParameters mpp201 (cs201, 18365);
  StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

  MassPoint mp202 (900, 50);
  StopCrossSection cs202 (0.0128895, 0.00195954);
  MassPointParameters mpp202 (cs202, 15322);
  StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

  MassPoint mp203 (900, 100);
  StopCrossSection cs203 (0.0128895, 0.00195954);
  MassPointParameters mpp203 (cs203, 19210);
  StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

  MassPoint mp204 (900, 150);
  StopCrossSection cs204 (0.0128895, 0.00195954);
  MassPointParameters mpp204 (cs204, 20162);
  StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

  MassPoint mp205 (900, 200);
  StopCrossSection cs205 (0.0128895, 0.00195954);
  MassPointParameters mpp205 (cs205, 18150);
  StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

  MassPoint mp206 (900, 250);
  StopCrossSection cs206 (0.0128895, 0.00195954);
  MassPointParameters mpp206 (cs206, 16871);
  StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

  MassPoint mp207 (900, 300);
  StopCrossSection cs207 (0.0128895, 0.00195954);
  MassPointParameters mpp207 (cs207, 20477);
  StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

  MassPoint mp208 (900, 350);
  StopCrossSection cs208 (0.0128895, 0.00195954);
  MassPointParameters mpp208 (cs208, 22252);
  StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

  MassPoint mp209 (900, 400);
  StopCrossSection cs209 (0.0128895, 0.00195954);
  MassPointParameters mpp209 (cs209, 19879);
  StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

  MassPoint mp210 (900, 450);
  StopCrossSection cs210 (0.0128895, 0.00195954);
  MassPointParameters mpp210 (cs210, 23165);
  StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

  MassPoint mp211 (900, 500);
  StopCrossSection cs211 (0.0128895, 0.00195954);
  MassPointParameters mpp211 (cs211, 18700);
  StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

  MassPoint mp212 (900, 550);
  StopCrossSection cs212 (0.0128895, 0.00195954);
  MassPointParameters mpp212 (cs212, 17235);
  StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

  MassPoint mp213 (900, 600);
  StopCrossSection cs213 (0.0128895, 0.00195954);
  MassPointParameters mpp213 (cs213, 20553);
  StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

  MassPoint mp214 (900, 625);
  StopCrossSection cs214 (0.0128895, 0.00195954);
  MassPointParameters mpp214 (cs214, 18591);
  StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

  MassPoint mp215 (900, 650);
  StopCrossSection cs215 (0.0128895, 0.00195954);
  MassPointParameters mpp215 (cs215, 18949);
  StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

  MassPoint mp216 (925, 625);
  StopCrossSection cs216 (0.0106631, 0.00165071);
  MassPointParameters mpp216 (cs216, 17086);
  StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

  MassPoint mp217 (925, 650);
  StopCrossSection cs217 (0.0106631, 0.00165071);
  MassPointParameters mpp217 (cs217, 20343);
  StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

  MassPoint mp218 (950, 1);
  StopCrossSection cs218 (0.00883465, 0.0013886);
  MassPointParameters mpp218 (cs218, 16187);
  StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

  MassPoint mp219 (950, 25);
  StopCrossSection cs219 (0.00883465, 0.0013886);
  MassPointParameters mpp219 (cs219, 17633);
  StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

  MassPoint mp220 (950, 50);
  StopCrossSection cs220 (0.00883465, 0.0013886);
  MassPointParameters mpp220 (cs220, 20421);
  StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

  MassPoint mp221 (950, 100);
  StopCrossSection cs221 (0.00883465, 0.0013886);
  MassPointParameters mpp221 (cs221, 19636);
  StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

  MassPoint mp222 (950, 150);
  StopCrossSection cs222 (0.00883465, 0.0013886);
  MassPointParameters mpp222 (cs222, 18950);
  StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

  MassPoint mp223 (950, 200);
  StopCrossSection cs223 (0.00883465, 0.0013886);
  MassPointParameters mpp223 (cs223, 20194);
  StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

  MassPoint mp224 (950, 250);
  StopCrossSection cs224 (0.00883465, 0.0013886);
  MassPointParameters mpp224 (cs224, 19800);
  StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

  MassPoint mp225 (950, 300);
  StopCrossSection cs225 (0.00883465, 0.0013886);
  MassPointParameters mpp225 (cs225, 20428);
  StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

  MassPoint mp226 (950, 350);
  StopCrossSection cs226 (0.00883465, 0.0013886);
  MassPointParameters mpp226 (cs226, 19136);
  StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

  MassPoint mp227 (950, 400);
  StopCrossSection cs227 (0.00883465, 0.0013886);
  MassPointParameters mpp227 (cs227, 20865);
  StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

  MassPoint mp228 (950, 450);
  StopCrossSection cs228 (0.00883465, 0.0013886);
  MassPointParameters mpp228 (cs228, 19950);
  StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

  MassPoint mp229 (950, 500);
  StopCrossSection cs229 (0.00883465, 0.0013886);
  MassPointParameters mpp229 (cs229, 21176);
  StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

  MassPoint mp230 (950, 550);
  StopCrossSection cs230 (0.00883465, 0.0013886);
  MassPointParameters mpp230 (cs230, 22074);
  StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

  MassPoint mp231 (950, 600);
  StopCrossSection cs231 (0.00883465, 0.0013886);
  MassPointParameters mpp231 (cs231, 19550);
  StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

  MassPoint mp232 (950, 650);
  StopCrossSection cs232 (0.00883465, 0.0013886);
  MassPointParameters mpp232 (cs232, 22255);
  StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

  MassPoint mp233 (1000, 1);
  StopCrossSection cs233 (0.00615134, 0.00100238);
  MassPointParameters mpp233 (cs233, 21092);
  StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

  MassPoint mp234 (1000, 25);
  StopCrossSection cs234 (0.00615134, 0.00100238);
  MassPointParameters mpp234 (cs234, 17999);
  StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

  MassPoint mp235 (1000, 50);
  StopCrossSection cs235 (0.00615134, 0.00100238);
  MassPointParameters mpp235 (cs235, 15539);
  StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

  MassPoint mp236 (1000, 100);
  StopCrossSection cs236 (0.00615134, 0.00100238);
  MassPointParameters mpp236 (cs236, 20869);
  StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

  MassPoint mp237 (1000, 150);
  StopCrossSection cs237 (0.00615134, 0.00100238);
  MassPointParameters mpp237 (cs237, 18984);
  StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

  MassPoint mp238 (1000, 200);
  StopCrossSection cs238 (0.00615134, 0.00100238);
  MassPointParameters mpp238 (cs238, 17474);
  StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

  MassPoint mp239 (1000, 250);
  StopCrossSection cs239 (0.00615134, 0.00100238);
  MassPointParameters mpp239 (cs239, 19055);
  StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

  MassPoint mp240 (1000, 300);
  StopCrossSection cs240 (0.00615134, 0.00100238);
  MassPointParameters mpp240 (cs240, 18716);
  StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

  MassPoint mp241 (1000, 350);
  StopCrossSection cs241 (0.00615134, 0.00100238);
  MassPointParameters mpp241 (cs241, 15660);
  StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

  MassPoint mp242 (1000, 400);
  StopCrossSection cs242 (0.00615134, 0.00100238);
  MassPointParameters mpp242 (cs242, 20424);
  StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

  MassPoint mp243 (1000, 450);
  StopCrossSection cs243 (0.00615134, 0.00100238);
  MassPointParameters mpp243 (cs243, 22462);
  StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

  MassPoint mp244 (1000, 500);
  StopCrossSection cs244 (0.00615134, 0.00100238);
  MassPointParameters mpp244 (cs244, 21309);
  StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

  MassPoint mp245 (1000, 550);
  StopCrossSection cs245 (0.00615134, 0.00100238);
  MassPointParameters mpp245 (cs245, 17671);
  StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

  MassPoint mp246 (1000, 600);
  StopCrossSection cs246 (0.00615134, 0.00100238);
  MassPointParameters mpp246 (cs246, 20393);
  StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

  MassPoint mp247 (1000, 650);
  StopCrossSection cs247 (0.00615134, 0.00100238);
  MassPointParameters mpp247 (cs247, 17726);
  StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

  MassPoint mp248 (1050, 1);
  StopCrossSection cs248 (0.00432261, 0.000725589);
  MassPointParameters mpp248 (cs248, 19562);
  StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

  MassPoint mp249 (1050, 25);
  StopCrossSection cs249 (0.00432261, 0.000725589);
  MassPointParameters mpp249 (cs249, 17507);
  StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

  MassPoint mp250 (1050, 50);
  StopCrossSection cs250 (0.00432261, 0.000725589);
  MassPointParameters mpp250 (cs250, 20717);
  StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

  MassPoint mp251 (1050, 100);
  StopCrossSection cs251 (0.00432261, 0.000725589);
  MassPointParameters mpp251 (cs251, 15745);
  StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

  MassPoint mp252 (1050, 150);
  StopCrossSection cs252 (0.00432261, 0.000725589);
  MassPointParameters mpp252 (cs252, 20857);
  StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

  MassPoint mp253 (1050, 200);
  StopCrossSection cs253 (0.00432261, 0.000725589);
  MassPointParameters mpp253 (cs253, 18603);
  StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

  MassPoint mp254 (1050, 250);
  StopCrossSection cs254 (0.00432261, 0.000725589);
  MassPointParameters mpp254 (cs254, 18054);
  StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

  MassPoint mp255 (1050, 300);
  StopCrossSection cs255 (0.00432261, 0.000725589);
  MassPointParameters mpp255 (cs255, 19381);
  StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

  MassPoint mp256 (1050, 350);
  StopCrossSection cs256 (0.00432261, 0.000725589);
  MassPointParameters mpp256 (cs256, 19238);
  StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

  MassPoint mp257 (1050, 400);
  StopCrossSection cs257 (0.00432261, 0.000725589);
  MassPointParameters mpp257 (cs257, 20859);
  StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

  MassPoint mp258 (1050, 450);
  StopCrossSection cs258 (0.00432261, 0.000725589);
  MassPointParameters mpp258 (cs258, 21656);
  StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

  MassPoint mp259 (1050, 500);
  StopCrossSection cs259 (0.00432261, 0.000725589);
  MassPointParameters mpp259 (cs259, 20780);
  StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

  MassPoint mp260 (1050, 550);
  StopCrossSection cs260 (0.00432261, 0.000725589);
  MassPointParameters mpp260 (cs260, 18268);
  StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

  MassPoint mp261 (1050, 600);
  StopCrossSection cs261 (0.00432261, 0.000725589);
  MassPointParameters mpp261 (cs261, 15805);
  StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

  MassPoint mp262 (1050, 650);
  StopCrossSection cs262 (0.00432261, 0.000725589);
  MassPointParameters mpp262 (cs262, 18991);
  StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

  MassPoint mp263 (1100, 1);
  StopCrossSection cs263 (0.00307413, 0.000532983);
  MassPointParameters mpp263 (cs263, 18885);
  StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

  MassPoint mp264 (1100, 25);
  StopCrossSection cs264 (0.00307413, 0.000532983);
  MassPointParameters mpp264 (cs264, 19779);
  StopNeutralinoMap.insert(std::make_pair(mp264, mpp264));

  MassPoint mp265 (1100, 50);
  StopCrossSection cs265 (0.00307413, 0.000532983);
  MassPointParameters mpp265 (cs265, 22552);
  StopNeutralinoMap.insert(std::make_pair(mp265, mpp265));

  MassPoint mp266 (1100, 100);
  StopCrossSection cs266 (0.00307413, 0.000532983);
  MassPointParameters mpp266 (cs266, 17366);
  StopNeutralinoMap.insert(std::make_pair(mp266, mpp266));

  MassPoint mp267 (1100, 150);
  StopCrossSection cs267 (0.00307413, 0.000532983);
  MassPointParameters mpp267 (cs267, 16108);
  StopNeutralinoMap.insert(std::make_pair(mp267, mpp267));

  MassPoint mp268 (1100, 200);
  StopCrossSection cs268 (0.00307413, 0.000532983);
  MassPointParameters mpp268 (cs268, 15763);
  StopNeutralinoMap.insert(std::make_pair(mp268, mpp268));

  MassPoint mp269 (1100, 250);
  StopCrossSection cs269 (0.00307413, 0.000532983);
  MassPointParameters mpp269 (cs269, 16791);
  StopNeutralinoMap.insert(std::make_pair(mp269, mpp269));

  MassPoint mp270 (1100, 300);
  StopCrossSection cs270 (0.00307413, 0.000532983);
  MassPointParameters mpp270 (cs270, 22229);
  StopNeutralinoMap.insert(std::make_pair(mp270, mpp270));

  MassPoint mp271 (1100, 350);
  StopCrossSection cs271 (0.00307413, 0.000532983);
  MassPointParameters mpp271 (cs271, 20562);
  StopNeutralinoMap.insert(std::make_pair(mp271, mpp271));

  MassPoint mp272 (1100, 400);
  StopCrossSection cs272 (0.00307413, 0.000532983);
  MassPointParameters mpp272 (cs272, 19404);
  StopNeutralinoMap.insert(std::make_pair(mp272, mpp272));

  MassPoint mp273 (1100, 450);
  StopCrossSection cs273 (0.00307413, 0.000532983);
  MassPointParameters mpp273 (cs273, 23138);
  StopNeutralinoMap.insert(std::make_pair(mp273, mpp273));

  MassPoint mp274 (1100, 500);
  StopCrossSection cs274 (0.00307413, 0.000532983);
  MassPointParameters mpp274 (cs274, 21000);
  StopNeutralinoMap.insert(std::make_pair(mp274, mpp274));

  MassPoint mp275 (1100, 550);
  StopCrossSection cs275 (0.00307413, 0.000532983);
  MassPointParameters mpp275 (cs275, 16980);
  StopNeutralinoMap.insert(std::make_pair(mp275, mpp275));

  MassPoint mp276 (1100, 600);
  StopCrossSection cs276 (0.00307413, 0.000532983);
  MassPointParameters mpp276 (cs276, 19039);
  StopNeutralinoMap.insert(std::make_pair(mp276, mpp276));

  MassPoint mp277 (1100, 650);
  StopCrossSection cs277 (0.00307413, 0.000532983);
  MassPointParameters mpp277 (cs277, 18708);
  StopNeutralinoMap.insert(std::make_pair(mp277, mpp277));

  MassPoint mp278 (1150, 1);
  StopCrossSection cs278 (0.00221047, 0.000396247);
  MassPointParameters mpp278 (cs278, 17675);
  StopNeutralinoMap.insert(std::make_pair(mp278, mpp278));

  MassPoint mp279 (1150, 25);
  StopCrossSection cs279 (0.00221047, 0.000396247);
  MassPointParameters mpp279 (cs279, 21435);
  StopNeutralinoMap.insert(std::make_pair(mp279, mpp279));

  MassPoint mp280 (1150, 50);
  StopCrossSection cs280 (0.00221047, 0.000396247);
  MassPointParameters mpp280 (cs280, 19505);
  StopNeutralinoMap.insert(std::make_pair(mp280, mpp280));

  MassPoint mp281 (1150, 100);
  StopCrossSection cs281 (0.00221047, 0.000396247);
  MassPointParameters mpp281 (cs281, 18396);
  StopNeutralinoMap.insert(std::make_pair(mp281, mpp281));

  MassPoint mp282 (1150, 150);
  StopCrossSection cs282 (0.00221047, 0.000396247);
  MassPointParameters mpp282 (cs282, 16143);
  StopNeutralinoMap.insert(std::make_pair(mp282, mpp282));

  MassPoint mp283 (1150, 200);
  StopCrossSection cs283 (0.00221047, 0.000396247);
  MassPointParameters mpp283 (cs283, 19285);
  StopNeutralinoMap.insert(std::make_pair(mp283, mpp283));

  MassPoint mp284 (1150, 250);
  StopCrossSection cs284 (0.00221047, 0.000396247);
  MassPointParameters mpp284 (cs284, 21000);
  StopNeutralinoMap.insert(std::make_pair(mp284, mpp284));

  MassPoint mp285 (1150, 300);
  StopCrossSection cs285 (0.00221047, 0.000396247);
  MassPointParameters mpp285 (cs285, 19163);
  StopNeutralinoMap.insert(std::make_pair(mp285, mpp285));

  MassPoint mp286 (1150, 350);
  StopCrossSection cs286 (0.00221047, 0.000396247);
  MassPointParameters mpp286 (cs286, 18627);
  StopNeutralinoMap.insert(std::make_pair(mp286, mpp286));

  MassPoint mp287 (1150, 400);
  StopCrossSection cs287 (0.00221047, 0.000396247);
  MassPointParameters mpp287 (cs287, 18949);
  StopNeutralinoMap.insert(std::make_pair(mp287, mpp287));

  MassPoint mp288 (1150, 450);
  StopCrossSection cs288 (0.00221047, 0.000396247);
  MassPointParameters mpp288 (cs288, 18411);
  StopNeutralinoMap.insert(std::make_pair(mp288, mpp288));

  MassPoint mp289 (1150, 500);
  StopCrossSection cs289 (0.00221047, 0.000396247);
  MassPointParameters mpp289 (cs289, 21965);
  StopNeutralinoMap.insert(std::make_pair(mp289, mpp289));

  MassPoint mp290 (1150, 550);
  StopCrossSection cs290 (0.00221047, 0.000396247);
  MassPointParameters mpp290 (cs290, 19979);
  StopNeutralinoMap.insert(std::make_pair(mp290, mpp290));

  MassPoint mp291 (1150, 600);
  StopCrossSection cs291 (0.00221047, 0.000396247);
  MassPointParameters mpp291 (cs291, 17004);
  StopNeutralinoMap.insert(std::make_pair(mp291, mpp291));

  MassPoint mp292 (1150, 650);
  StopCrossSection cs292 (0.00221047, 0.000396247);
  MassPointParameters mpp292 (cs292, 20394);
  StopNeutralinoMap.insert(std::make_pair(mp292, mpp292));

  MassPoint mp293 (1200, 1);
  StopCrossSection cs293 (0.00159844, 0.000296045);
  MassPointParameters mpp293 (cs293, 19941);
  StopNeutralinoMap.insert(std::make_pair(mp293, mpp293));

  MassPoint mp294 (1200, 25);
  StopCrossSection cs294 (0.00159844, 0.000296045);
  MassPointParameters mpp294 (cs294, 20714);
  StopNeutralinoMap.insert(std::make_pair(mp294, mpp294));

  MassPoint mp295 (1200, 50);
  StopCrossSection cs295 (0.00159844, 0.000296045);
  MassPointParameters mpp295 (cs295, 18719);
  StopNeutralinoMap.insert(std::make_pair(mp295, mpp295));

  MassPoint mp296 (1200, 100);
  StopCrossSection cs296 (0.00159844, 0.000296045);
  MassPointParameters mpp296 (cs296, 19162);
  StopNeutralinoMap.insert(std::make_pair(mp296, mpp296));

  MassPoint mp297 (1200, 150);
  StopCrossSection cs297 (0.00159844, 0.000296045);
  MassPointParameters mpp297 (cs297, 22916);
  StopNeutralinoMap.insert(std::make_pair(mp297, mpp297));

  MassPoint mp298 (1200, 200);
  StopCrossSection cs298 (0.00159844, 0.000296045);
  MassPointParameters mpp298 (cs298, 22456);
  StopNeutralinoMap.insert(std::make_pair(mp298, mpp298));

  MassPoint mp299 (1200, 250);
  StopCrossSection cs299 (0.00159844, 0.000296045);
  MassPointParameters mpp299 (cs299, 20413);
  StopNeutralinoMap.insert(std::make_pair(mp299, mpp299));

  MassPoint mp300 (1200, 300);
  StopCrossSection cs300 (0.00159844, 0.000296045);
  MassPointParameters mpp300 (cs300, 16368);
  StopNeutralinoMap.insert(std::make_pair(mp300, mpp300));

  MassPoint mp301 (1200, 350);
  StopCrossSection cs301 (0.00159844, 0.000296045);
  MassPointParameters mpp301 (cs301, 23483);
  StopNeutralinoMap.insert(std::make_pair(mp301, mpp301));

  MassPoint mp302 (1200, 400);
  StopCrossSection cs302 (0.00159844, 0.000296045);
  MassPointParameters mpp302 (cs302, 17650);
  StopNeutralinoMap.insert(std::make_pair(mp302, mpp302));

  MassPoint mp303 (1200, 450);
  StopCrossSection cs303 (0.00159844, 0.000296045);
  MassPointParameters mpp303 (cs303, 17878);
  StopNeutralinoMap.insert(std::make_pair(mp303, mpp303));

  MassPoint mp304 (1200, 500);
  StopCrossSection cs304 (0.00159844, 0.000296045);
  MassPointParameters mpp304 (cs304, 22338);
  StopNeutralinoMap.insert(std::make_pair(mp304, mpp304));

  MassPoint mp305 (1200, 550);
  StopCrossSection cs305 (0.00159844, 0.000296045);
  MassPointParameters mpp305 (cs305, 18919);
  StopNeutralinoMap.insert(std::make_pair(mp305, mpp305));

  MassPoint mp306 (1200, 600);
  StopCrossSection cs306 (0.00159844, 0.000296045);
  MassPointParameters mpp306 (cs306, 18137);
  StopNeutralinoMap.insert(std::make_pair(mp306, mpp306));

  MassPoint mp307 (1200, 650);
  StopCrossSection cs307 (0.00159844, 0.000296045);
  MassPointParameters mpp307 (cs307, 19046);
  StopNeutralinoMap.insert(std::make_pair(mp307, mpp307)); 

  } else if (SUSYProductionProcess=="T2tt_dM10to80") {

      MassPoint mp0 (250, 170);
  StopCrossSection cs0 (2.26682, 0.318703);
  MassPointParameters mpp0 (cs0, 95436);
  StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

  MassPoint mp1 (250, 180);
  StopCrossSection cs1 (2.26682, 0.318703);
  MassPointParameters mpp1 (cs1, 91978);
  StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

  MassPoint mp2 (250, 190);
  StopCrossSection cs2 (2.26682, 0.318703);
  MassPointParameters mpp2 (cs2, 91004);
  StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

  MassPoint mp3 (250, 200);
  StopCrossSection cs3 (2.26682, 0.318703);
  MassPointParameters mpp3 (cs3, 91824);
  StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

  MassPoint mp4 (250, 210);
  StopCrossSection cs4 (2.26682, 0.318703);
  MassPointParameters mpp4 (cs4, 100564);
  StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

  MassPoint mp5 (250, 220);
  StopCrossSection cs5 (2.26682, 0.318703);
  MassPointParameters mpp5 (cs5, 101832);
  StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

  MassPoint mp6 (250, 230);
  StopCrossSection cs6 (2.26682, 0.318703);
  MassPointParameters mpp6 (cs6, 88158);
  StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

  MassPoint mp7 (250, 240);
  StopCrossSection cs7 (2.26682, 0.318703);
  MassPointParameters mpp7 (cs7, 91853);
  StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

  MassPoint mp8 (275, 195);
  StopCrossSection cs8 (1.39853, 0.199358);
  MassPointParameters mpp8 (cs8, 89132);
  StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

  MassPoint mp9 (275, 205);
  StopCrossSection cs9 (1.39853, 0.199358);
  MassPointParameters mpp9 (cs9, 84020);
  StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

  MassPoint mp10 (275, 215);
  StopCrossSection cs10 (1.39853, 0.199358);
  MassPointParameters mpp10 (cs10, 85779);
  StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

  MassPoint mp11 (275, 225);
  StopCrossSection cs11 (1.39853, 0.199358);
  MassPointParameters mpp11 (cs11, 89781);
  StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

  MassPoint mp12 (275, 235);
  StopCrossSection cs12 (1.39853, 0.199358);
  MassPointParameters mpp12 (cs12, 81058);
  StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

  MassPoint mp13 (275, 245);
  StopCrossSection cs13 (1.39853, 0.199358);
  MassPointParameters mpp13 (cs13, 85145);
  StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

  MassPoint mp14 (275, 255);
  StopCrossSection cs14 (1.39853, 0.199358);
  MassPointParameters mpp14 (cs14, 84281);
  StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

  MassPoint mp15 (275, 265);
  StopCrossSection cs15 (1.39853, 0.199358);
  MassPointParameters mpp15 (cs15, 81770);
  StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

  MassPoint mp16 (300, 220);
  StopCrossSection cs16 (0.89394, 0.124457);
  MassPointParameters mpp16 (cs16, 96542);
  StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

  MassPoint mp17 (300, 230);
  StopCrossSection cs17 (0.89394, 0.124457);
  MassPointParameters mpp17 (cs17, 93347);
  StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

  MassPoint mp18 (300, 240);
  StopCrossSection cs18 (0.89394, 0.124457);
  MassPointParameters mpp18 (cs18, 98572);
  StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

  MassPoint mp19 (300, 250);
  StopCrossSection cs19 (0.89394, 0.124457);
  MassPointParameters mpp19 (cs19, 93643);
  StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

  MassPoint mp20 (300, 260);
  StopCrossSection cs20 (0.89394, 0.124457);
  MassPointParameters mpp20 (cs20, 95497);
  StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

  MassPoint mp21 (300, 270);
  StopCrossSection cs21 (0.89394, 0.124457);
  MassPointParameters mpp21 (cs21, 102308);
  StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

  MassPoint mp22 (300, 280);
  StopCrossSection cs22 (0.89394, 0.124457);
  MassPointParameters mpp22 (cs22, 94657);
  StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

  MassPoint mp23 (300, 290);
  StopCrossSection cs23 (0.89394, 0.124457);
  MassPointParameters mpp23 (cs23, 102608);
  StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

  MassPoint mp24 (325, 245);
  StopCrossSection cs24 (0.588326, 0.0812738);
  MassPointParameters mpp24 (cs24, 94846);
  StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

  MassPoint mp25 (325, 255);
  StopCrossSection cs25 (0.588326, 0.0812738);
  MassPointParameters mpp25 (cs25, 95462);
  StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

  MassPoint mp26 (325, 265);
  StopCrossSection cs26 (0.588326, 0.0812738);
  MassPointParameters mpp26 (cs26, 95331);
  StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

  MassPoint mp27 (325, 275);
  StopCrossSection cs27 (0.588326, 0.0812738);
  MassPointParameters mpp27 (cs27, 97742);
  StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

  MassPoint mp28 (325, 285);
  StopCrossSection cs28 (0.588326, 0.0812738);
  MassPointParameters mpp28 (cs28, 96751);
  StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

  MassPoint mp29 (325, 295);
  StopCrossSection cs29 (0.588326, 0.0812738);
  MassPointParameters mpp29 (cs29, 95524);
  StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

  MassPoint mp30 (325, 305);
  StopCrossSection cs30 (0.588326, 0.0812738);
  MassPointParameters mpp30 (cs30, 95966);
  StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

  MassPoint mp31 (325, 315);
  StopCrossSection cs31 (0.588326, 0.0812738);
  MassPointParameters mpp31 (cs31, 89832);
  StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

  MassPoint mp32 (350, 270);
  StopCrossSection cs32 (0.39748, 0.0544059);
  MassPointParameters mpp32 (cs32, 87849);
  StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

  MassPoint mp33 (350, 280);
  StopCrossSection cs33 (0.39748, 0.0544059);
  MassPointParameters mpp33 (cs33, 95563);
  StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

  MassPoint mp34 (350, 290);
  StopCrossSection cs34 (0.39748, 0.0544059);
  MassPointParameters mpp34 (cs34, 91776);
  StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

  MassPoint mp35 (350, 300);
  StopCrossSection cs35 (0.39748, 0.0544059);
  MassPointParameters mpp35 (cs35, 94245);
  StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

  MassPoint mp36 (350, 310);
  StopCrossSection cs36 (0.39748, 0.0544059);
  MassPointParameters mpp36 (cs36, 89530);
  StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

  MassPoint mp37 (350, 320);
  StopCrossSection cs37 (0.39748, 0.0544059);
  MassPointParameters mpp37 (cs37, 91496);
  StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

  MassPoint mp38 (350, 330);
  StopCrossSection cs38 (0.39748, 0.0544059);
  MassPointParameters mpp38 (cs38, 99128);
  StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

  MassPoint mp39 (350, 340);
  StopCrossSection cs39 (0.39748, 0.0544059);
  MassPointParameters mpp39 (cs39, 88049);
  StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

  MassPoint mp40 (375, 295);
  StopCrossSection cs40 (0.274142, 0.0379623);
  MassPointParameters mpp40 (cs40, 90692);
  StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

  MassPoint mp41 (375, 305);
  StopCrossSection cs41 (0.274142, 0.0379623);
  MassPointParameters mpp41 (cs41, 89846);
  StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

  MassPoint mp42 (375, 315);
  StopCrossSection cs42 (0.274142, 0.0379623);
  MassPointParameters mpp42 (cs42, 96882);
  StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

  MassPoint mp43 (375, 325);
  StopCrossSection cs43 (0.274142, 0.0379623);
  MassPointParameters mpp43 (cs43, 87138);
  StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

  MassPoint mp44 (375, 335);
  StopCrossSection cs44 (0.274142, 0.0379623);
  MassPointParameters mpp44 (cs44, 89899);
  StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

  MassPoint mp45 (375, 345);
  StopCrossSection cs45 (0.274142, 0.0379623);
  MassPointParameters mpp45 (cs45, 90039);
  StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

  MassPoint mp46 (375, 355);
  StopCrossSection cs46 (0.274142, 0.0379623);
  MassPointParameters mpp46 (cs46, 90712);
  StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

  MassPoint mp47 (375, 365);
  StopCrossSection cs47 (0.274142, 0.0379623);
  MassPointParameters mpp47 (cs47, 90014);
  StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

  MassPoint mp48 (400, 320);
  StopCrossSection cs48 (0.192659, 0.0263914);
  MassPointParameters mpp48 (cs48, 66572);
  StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

  MassPoint mp49 (400, 330);
  StopCrossSection cs49 (0.192659, 0.0263914);
  MassPointParameters mpp49 (cs49, 71479);
  StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

  MassPoint mp50 (400, 340);
  StopCrossSection cs50 (0.192659, 0.0263914);
  MassPointParameters mpp50 (cs50, 75022);
  StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

  MassPoint mp51 (400, 350);
  StopCrossSection cs51 (0.192659, 0.0263914);
  MassPointParameters mpp51 (cs51, 73463);
  StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

  MassPoint mp52 (400, 360);
  StopCrossSection cs52 (0.192659, 0.0263914);
  MassPointParameters mpp52 (cs52, 71290);
  StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

  MassPoint mp53 (400, 370);
  StopCrossSection cs53 (0.192659, 0.0263914);
  MassPointParameters mpp53 (cs53, 74126);
  StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

  MassPoint mp54 (400, 380);
  StopCrossSection cs54 (0.192659, 0.0263914);
  MassPointParameters mpp54 (cs54, 69293);
  StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

  MassPoint mp55 (400, 390);
  StopCrossSection cs55 (0.192659, 0.0263914);
  MassPointParameters mpp55 (cs55, 71957);
  StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

  MassPoint mp56 (425, 345);
  StopCrossSection cs56 (0.137688, 0.0185897);
  MassPointParameters mpp56 (cs56, 51989);
  StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

  MassPoint mp57 (425, 355);
  StopCrossSection cs57 (0.137688, 0.0185897);
  MassPointParameters mpp57 (cs57, 48788);
  StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

  MassPoint mp58 (425, 365);
  StopCrossSection cs58 (0.137688, 0.0185897);
  MassPointParameters mpp58 (cs58, 50397);
  StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

  MassPoint mp59 (425, 375);
  StopCrossSection cs59 (0.137688, 0.0185897);
  MassPointParameters mpp59 (cs59, 51463);
  StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

  MassPoint mp60 (425, 385);
  StopCrossSection cs60 (0.137688, 0.0185897);
  MassPointParameters mpp60 (cs60, 52567);
  StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

  MassPoint mp61 (425, 395);
  StopCrossSection cs61 (0.137688, 0.0185897);
  MassPointParameters mpp61 (cs61, 51532);
  StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

  MassPoint mp62 (425, 405);
  StopCrossSection cs62 (0.137688, 0.0185897);
  MassPointParameters mpp62 (cs62, 54543);
  StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

  MassPoint mp63 (425, 415);
  StopCrossSection cs63 (0.137688, 0.0185897);
  MassPointParameters mpp63 (cs63, 44811);
  StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

  MassPoint mp64 (450, 370);
  StopCrossSection cs64 (0.0995465, 0.0133949);
  MassPointParameters mpp64 (cs64, 34040);
  StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

  MassPoint mp65 (450, 380);
  StopCrossSection cs65 (0.0995465, 0.0133949);
  MassPointParameters mpp65 (cs65, 32869);
  StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

  MassPoint mp66 (450, 390);
  StopCrossSection cs66 (0.0995465, 0.0133949);
  MassPointParameters mpp66 (cs66, 35807);
  StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

  MassPoint mp67 (450, 400);
  StopCrossSection cs67 (0.0995465, 0.0133949);
  MassPointParameters mpp67 (cs67, 37792);
  StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

  MassPoint mp68 (450, 410);
  StopCrossSection cs68 (0.0995465, 0.0133949);
  MassPointParameters mpp68 (cs68, 34377);
  StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

  MassPoint mp69 (450, 420);
  StopCrossSection cs69 (0.0995465, 0.0133949);
  MassPointParameters mpp69 (cs69, 34935);
  StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

  MassPoint mp70 (450, 430);
  StopCrossSection cs70 (0.0995465, 0.0133949);
  MassPointParameters mpp70 (cs70, 33601);
  StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

  MassPoint mp71 (450, 440);
  StopCrossSection cs71 (0.0995465, 0.0133949);
  MassPointParameters mpp71 (cs71, 35024);
  StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

  MassPoint mp72 (475, 395);
  StopCrossSection cs72 (0.073172, 0.00979963);
  MassPointParameters mpp72 (cs72, 26368);
  StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

  MassPoint mp73 (475, 405);
  StopCrossSection cs73 (0.073172, 0.00979963);
  MassPointParameters mpp73 (cs73, 26976);
  StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

  MassPoint mp74 (475, 415);
  StopCrossSection cs74 (0.073172, 0.00979963);
  MassPointParameters mpp74 (cs74, 26169);
  StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

  MassPoint mp75 (475, 425);
  StopCrossSection cs75 (0.073172, 0.00979963);
  MassPointParameters mpp75 (cs75, 23895);
  StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

  MassPoint mp76 (475, 435);
  StopCrossSection cs76 (0.073172, 0.00979963);
  MassPointParameters mpp76 (cs76, 27329);
  StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

  MassPoint mp77 (475, 445);
  StopCrossSection cs77 (0.073172, 0.00979963);
  MassPointParameters mpp77 (cs77, 24556);
  StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

  MassPoint mp78 (475, 455);
  StopCrossSection cs78 (0.073172, 0.00979963);
  MassPointParameters mpp78 (cs78, 25150);
  StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

  MassPoint mp79 (475, 465);
  StopCrossSection cs79 (0.073172, 0.00979963);
  MassPointParameters mpp79 (cs79, 25689);
  StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

  MassPoint mp80 (500, 420);
  StopCrossSection cs80 (0.0544248, 0.00728188);
  MassPointParameters mpp80 (cs80, 23764);
  StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

  MassPoint mp81 (500, 430);
  StopCrossSection cs81 (0.0544248, 0.00728188);
  MassPointParameters mpp81 (cs81, 19921);
  StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

  MassPoint mp82 (500, 440);
  StopCrossSection cs82 (0.0544248, 0.00728188);
  MassPointParameters mpp82 (cs82, 18723);
  StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

  MassPoint mp83 (500, 450);
  StopCrossSection cs83 (0.0544248, 0.00728188);
  MassPointParameters mpp83 (cs83, 23976);
  StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

  MassPoint mp84 (500, 460);
  StopCrossSection cs84 (0.0544248, 0.00728188);
  MassPointParameters mpp84 (cs84, 18363);
  StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

  MassPoint mp85 (500, 470);
  StopCrossSection cs85 (0.0544248, 0.00728188);
  MassPointParameters mpp85 (cs85, 25042);
  StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

  MassPoint mp86 (500, 480);
  StopCrossSection cs86 (0.0544248, 0.00728188);
  MassPointParameters mpp86 (cs86, 20462);
  StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

  MassPoint mp87 (500, 490);
  StopCrossSection cs87 (0.0544248, 0.00728188);
  MassPointParameters mpp87 (cs87, 20792);
  StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

  MassPoint mp88 (525, 445);
  StopCrossSection cs88 (0.0409701, 0.00546717);
  MassPointParameters mpp88 (cs88, 14432);
  StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

  MassPoint mp89 (525, 455);
  StopCrossSection cs89 (0.0409701, 0.00546717);
  MassPointParameters mpp89 (cs89, 14956);
  StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

  MassPoint mp90 (525, 465);
  StopCrossSection cs90 (0.0409701, 0.00546717);
  MassPointParameters mpp90 (cs90, 16337);
  StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

  MassPoint mp91 (525, 475);
  StopCrossSection cs91 (0.0409701, 0.00546717);
  MassPointParameters mpp91 (cs91, 14904);
  StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

  MassPoint mp92 (525, 485);
  StopCrossSection cs92 (0.0409701, 0.00546717);
  MassPointParameters mpp92 (cs92, 18261);
  StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

  MassPoint mp93 (525, 495);
  StopCrossSection cs93 (0.0409701, 0.00546717);
  MassPointParameters mpp93 (cs93, 14707);
  StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

  MassPoint mp94 (525, 505);
  StopCrossSection cs94 (0.0409701, 0.00546717);
  MassPointParameters mpp94 (cs94, 15106);
  StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

  MassPoint mp95 (525, 515);
  StopCrossSection cs95 (0.0409701, 0.00546717);
  MassPointParameters mpp95 (cs95, 14954);
  StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

  MassPoint mp96 (550, 470);
  StopCrossSection cs96 (0.0310846, 0.00412452);
  MassPointParameters mpp96 (cs96, 11150);
  StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

  MassPoint mp97 (550, 480);
  StopCrossSection cs97 (0.0310846, 0.00412452);
  MassPointParameters mpp97 (cs97, 11139);
  StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

  MassPoint mp98 (550, 490);
  StopCrossSection cs98 (0.0310846, 0.00412452);
  MassPointParameters mpp98 (cs98, 11039);
  StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

  MassPoint mp99 (550, 500);
  StopCrossSection cs99 (0.0310846, 0.00412452);
  MassPointParameters mpp99 (cs99, 10528);
  StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

  MassPoint mp100 (550, 510);
  StopCrossSection cs100 (0.0310846, 0.00412452);
  MassPointParameters mpp100 (cs100, 10116);
  StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

  MassPoint mp101 (550, 520);
  StopCrossSection cs101 (0.0310846, 0.00412452);
  MassPointParameters mpp101 (cs101, 9730);
  StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

  MassPoint mp102 (550, 530);
  StopCrossSection cs102 (0.0310846, 0.00412452);
  MassPointParameters mpp102 (cs102, 11692);
  StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

  MassPoint mp103 (550, 540);
  StopCrossSection cs103 (0.0310846, 0.00412452);
  MassPointParameters mpp103 (cs103, 12097);
  StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

  MassPoint mp104 (575, 495);
  StopCrossSection cs104 (0.0237356, 0.00315069);
  MassPointParameters mpp104 (cs104, 8357);
  StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

  MassPoint mp105 (575, 505);
  StopCrossSection cs105 (0.0237356, 0.00315069);
  MassPointParameters mpp105 (cs105, 8136);
  StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

  MassPoint mp106 (575, 515);
  StopCrossSection cs106 (0.0237356, 0.00315069);
  MassPointParameters mpp106 (cs106, 9021);
  StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

  MassPoint mp107 (575, 525);
  StopCrossSection cs107 (0.0237356, 0.00315069);
  MassPointParameters mpp107 (cs107, 9636);
  StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

  MassPoint mp108 (575, 535);
  StopCrossSection cs108 (0.0237356, 0.00315069);
  MassPointParameters mpp108 (cs108, 7998);
  StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

  MassPoint mp109 (575, 545);
  StopCrossSection cs109 (0.0237356, 0.00315069);
  MassPointParameters mpp109 (cs109, 7718);
  StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

  MassPoint mp110 (575, 555);
  StopCrossSection cs110 (0.0237356, 0.00315069);
  MassPointParameters mpp110 (cs110, 7390);
  StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

  MassPoint mp111 (575, 565);
  StopCrossSection cs111 (0.0237356, 0.00315069);
  MassPointParameters mpp111 (cs111, 9702);
  StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

  MassPoint mp112 (600, 520);
  StopCrossSection cs112 (0.0183277, 0.00242061);
  MassPointParameters mpp112 (cs112, 6349);
  StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

  MassPoint mp113 (600, 530);
  StopCrossSection cs113 (0.0183277, 0.00242061);
  MassPointParameters mpp113 (cs113, 9048);
  StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

  MassPoint mp114 (600, 540);
  StopCrossSection cs114 (0.0183277, 0.00242061);
  MassPointParameters mpp114 (cs114, 7185);
  StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

  MassPoint mp115 (600, 550);
  StopCrossSection cs115 (0.0183277, 0.00242061);
  MassPointParameters mpp115 (cs115, 5744);
  StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

  MassPoint mp116 (600, 560);
  StopCrossSection cs116 (0.0183277, 0.00242061);
  MassPointParameters mpp116 (cs116, 8702);
  StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

  MassPoint mp117 (600, 570);
  StopCrossSection cs117 (0.0183277, 0.00242061);
  MassPointParameters mpp117 (cs117, 7448);
  StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

  MassPoint mp118 (600, 580);
  StopCrossSection cs118 (0.0183277, 0.00242061);
  MassPointParameters mpp118 (cs118, 8031);
  StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

  MassPoint mp119 (600, 590);
  StopCrossSection cs119 (0.0183277, 0.00242061);
  MassPointParameters mpp119 (cs119, 6549);
  StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

  MassPoint mp120 (625, 545);
  StopCrossSection cs120 (0.014315, 0.00183177);
  MassPointParameters mpp120 (cs120, 5638);
  StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

  MassPoint mp121 (625, 555);
  StopCrossSection cs121 (0.014315, 0.00183177);
  MassPointParameters mpp121 (cs121, 5012);
  StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

  MassPoint mp122 (625, 565);
  StopCrossSection cs122 (0.014315, 0.00183177);
  MassPointParameters mpp122 (cs122, 5341);
  StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

  MassPoint mp123 (625, 575);
  StopCrossSection cs123 (0.014315, 0.00183177);
  MassPointParameters mpp123 (cs123, 5232);
  StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

  MassPoint mp124 (625, 585);
  StopCrossSection cs124 (0.014315, 0.00183177);
  MassPointParameters mpp124 (cs124, 5828);
  StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

  MassPoint mp125 (625, 595);
  StopCrossSection cs125 (0.014315, 0.00183177);
  MassPointParameters mpp125 (cs125, 6161);
  StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

  MassPoint mp126 (625, 605);
  StopCrossSection cs126 (0.014315, 0.00183177);
  MassPointParameters mpp126 (cs126, 7245);
  StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

  MassPoint mp127 (625, 615);
  StopCrossSection cs127 (0.014315, 0.00183177);
  MassPointParameters mpp127 (cs127, 6701);
  StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

  MassPoint mp128 (650, 570);
  StopCrossSection cs128 (0.0112365, 0.00145212);
  MassPointParameters mpp128 (cs128, 4476);
  StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

  MassPoint mp129 (650, 580);
  StopCrossSection cs129 (0.0112365, 0.00145212);
  MassPointParameters mpp129 (cs129, 5149);
  StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

  MassPoint mp130 (650, 590);
  StopCrossSection cs130 (0.0112365, 0.00145212);
  MassPointParameters mpp130 (cs130, 4424);
  StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

  MassPoint mp131 (650, 600);
  StopCrossSection cs131 (0.0112365, 0.00145212);
  MassPointParameters mpp131 (cs131, 5439);
  StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

  MassPoint mp132 (650, 610);
  StopCrossSection cs132 (0.0112365, 0.00145212);
  MassPointParameters mpp132 (cs132, 5604);
  StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

  MassPoint mp133 (650, 620);
  StopCrossSection cs133 (0.0112365, 0.00145212);
  MassPointParameters mpp133 (cs133, 4554);
  StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

  MassPoint mp134 (650, 630);
  StopCrossSection cs134 (0.0112365, 0.00145212);
  MassPointParameters mpp134 (cs134, 3536);
  StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

  MassPoint mp135 (650, 640);
  StopCrossSection cs135 (0.0112365, 0.00145212);
  MassPointParameters mpp135 (cs135, 5063);
  StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

  MassPoint mp136 (675, 595);
  StopCrossSection cs136 (0.00886867, 0.00115916);
  MassPointParameters mpp136 (cs136, 4522);
  StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

  MassPoint mp137 (675, 605);
  StopCrossSection cs137 (0.00886867, 0.00115916);
  MassPointParameters mpp137 (cs137, 2964);
  StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

  MassPoint mp138 (675, 615);
  StopCrossSection cs138 (0.00886867, 0.00115916);
  MassPointParameters mpp138 (cs138, 4586);
  StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

  MassPoint mp139 (675, 625);
  StopCrossSection cs139 (0.00886867, 0.00115916);
  MassPointParameters mpp139 (cs139, 5479);
  StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

  MassPoint mp140 (675, 635);
  StopCrossSection cs140 (0.00886867, 0.00115916);
  MassPointParameters mpp140 (cs140, 4258);
  StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

  MassPoint mp141 (675, 645);
  StopCrossSection cs141 (0.00886867, 0.00115916);
  MassPointParameters mpp141 (cs141, 4391);
  StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

  MassPoint mp142 (675, 655);
  StopCrossSection cs142 (0.00886867, 0.00115916);
  MassPointParameters mpp142 (cs142, 5943);
  StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

  MassPoint mp143 (675, 665);
  StopCrossSection cs143 (0.00886867, 0.00115916);
  MassPointParameters mpp143 (cs143, 5106);
  StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

  MassPoint mp144 (700, 620);
  StopCrossSection cs144 (0.00703799, 0.000939072);
  MassPointParameters mpp144 (cs144, 6403);
  StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

  MassPoint mp145 (700, 630);
  StopCrossSection cs145 (0.00703799, 0.000939072);
  MassPointParameters mpp145 (cs145, 5467);
  StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

  MassPoint mp146 (700, 640);
  StopCrossSection cs146 (0.00703799, 0.000939072);
  MassPointParameters mpp146 (cs146, 4510);
  StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

  MassPoint mp147 (700, 650);
  StopCrossSection cs147 (0.00703799, 0.000939072);
  MassPointParameters mpp147 (cs147, 4710);
  StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

  MassPoint mp148 (700, 660);
  StopCrossSection cs148 (0.00703799, 0.000939072);
  MassPointParameters mpp148 (cs148, 5143);
  StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

  MassPoint mp149 (700, 670);
  StopCrossSection cs149 (0.00703799, 0.000939072);
  MassPointParameters mpp149 (cs149, 6092);
  StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

  MassPoint mp150 (700, 680);
  StopCrossSection cs150 (0.00703799, 0.000939072);
  MassPointParameters mpp150 (cs150, 4422);
  StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

  MassPoint mp151 (700, 690);
  StopCrossSection cs151 (0.00703799, 0.000939072);
  MassPointParameters mpp151 (cs151, 4798);
  StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

  MassPoint mp152 (725, 645);
  StopCrossSection cs152 (0.00563099, 0.000764711);
  MassPointParameters mpp152 (cs152, 4559);
  StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

  MassPoint mp153 (725, 655);
  StopCrossSection cs153 (0.00563099, 0.000764711);
  MassPointParameters mpp153 (cs153, 6432);
  StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

  MassPoint mp154 (725, 665);
  StopCrossSection cs154 (0.00563099, 0.000764711);
  MassPointParameters mpp154 (cs154, 3755);
  StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

  MassPoint mp155 (725, 675);
  StopCrossSection cs155 (0.00563099, 0.000764711);
  MassPointParameters mpp155 (cs155, 4455);
  StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

  MassPoint mp156 (725, 685);
  StopCrossSection cs156 (0.00563099, 0.000764711);
  MassPointParameters mpp156 (cs156, 4395);
  StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

  MassPoint mp157 (725, 695);
  StopCrossSection cs157 (0.00563099, 0.000764711);
  MassPointParameters mpp157 (cs157, 3838);
  StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

  MassPoint mp158 (725, 705);
  StopCrossSection cs158 (0.00563099, 0.000764711);
  MassPointParameters mpp158 (cs158, 5015);
  StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

  MassPoint mp159 (725, 715);
  StopCrossSection cs159 (0.00563099, 0.000764711);
  MassPointParameters mpp159 (cs159, 4430);
  StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

  MassPoint mp160 (750, 670);
  StopCrossSection cs160 (0.00452859, 0.000622478);
  MassPointParameters mpp160 (cs160, 4655);
  StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

  MassPoint mp161 (750, 680);
  StopCrossSection cs161 (0.00452859, 0.000622478);
  MassPointParameters mpp161 (cs161, 5476);
  StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

  MassPoint mp162 (750, 690);
  StopCrossSection cs162 (0.00452859, 0.000622478);
  MassPointParameters mpp162 (cs162, 5010);
  StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

  MassPoint mp163 (750, 700);
  StopCrossSection cs163 (0.00452859, 0.000622478);
  MassPointParameters mpp163 (cs163, 4744);
  StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

  MassPoint mp164 (750, 710);
  StopCrossSection cs164 (0.00452859, 0.000622478);
  MassPointParameters mpp164 (cs164, 6815);
  StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

  MassPoint mp165 (750, 720);
  StopCrossSection cs165 (0.00452859, 0.000622478);
  MassPointParameters mpp165 (cs165, 5300);
  StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

  MassPoint mp166 (750, 730);
  StopCrossSection cs166 (0.00452859, 0.000622478);
  MassPointParameters mpp166 (cs166, 4679);
  StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

  MassPoint mp167 (750, 740);
  StopCrossSection cs167 (0.00452859, 0.000622478);
  MassPointParameters mpp167 (cs167, 6185);
  StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

  MassPoint mp168 (775, 695);
  StopCrossSection cs168 (0.00366131, 0.000511108);
  MassPointParameters mpp168 (cs168, 5177);
  StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

  MassPoint mp169 (775, 705);
  StopCrossSection cs169 (0.00366131, 0.000511108);
  MassPointParameters mpp169 (cs169, 4634);
  StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

  MassPoint mp170 (775, 715);
  StopCrossSection cs170 (0.00366131, 0.000511108);
  MassPointParameters mpp170 (cs170, 4192);
  StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

  MassPoint mp171 (775, 725);
  StopCrossSection cs171 (0.00366131, 0.000511108);
  MassPointParameters mpp171 (cs171, 5002);
  StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

  MassPoint mp172 (775, 735);
  StopCrossSection cs172 (0.00366131, 0.000511108);
  MassPointParameters mpp172 (cs172, 7060);
  StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

  MassPoint mp173 (775, 745);
  StopCrossSection cs173 (0.00366131, 0.000511108);
  MassPointParameters mpp173 (cs173, 3574);
  StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

  MassPoint mp174 (775, 755);
  StopCrossSection cs174 (0.00366131, 0.000511108);
  MassPointParameters mpp174 (cs174, 4180);
  StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

  MassPoint mp175 (775, 765);
  StopCrossSection cs175 (0.00366131, 0.000511108);
  MassPointParameters mpp175 (cs175, 3970);
  StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

  MassPoint mp176 (800, 720);
  StopCrossSection cs176 (0.0029742, 0.000421474);
  MassPointParameters mpp176 (cs176, 5515);
  StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

  MassPoint mp177 (800, 730);
  StopCrossSection cs177 (0.0029742, 0.000421474);
  MassPointParameters mpp177 (cs177, 4158);
  StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

  MassPoint mp178 (800, 740);
  StopCrossSection cs178 (0.0029742, 0.000421474);
  MassPointParameters mpp178 (cs178, 4079);
  StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

  MassPoint mp179 (800, 750);
  StopCrossSection cs179 (0.0029742, 0.000421474);
  MassPointParameters mpp179 (cs179, 5820);
  StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

  MassPoint mp180 (800, 760);
  StopCrossSection cs180 (0.0029742, 0.000421474);
  MassPointParameters mpp180 (cs180, 4377);
  StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

  MassPoint mp181 (800, 770);
  StopCrossSection cs181 (0.0029742, 0.000421474);
  MassPointParameters mpp181 (cs181, 4446);
  StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

  MassPoint mp182 (800, 780);
  StopCrossSection cs182 (0.0029742, 0.000421474);
  MassPointParameters mpp182 (cs182, 3409);
  StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

  MassPoint mp183 (800, 790);
  StopCrossSection cs183 (0.0029742, 0.000421474);
  MassPointParameters mpp183 (cs183, 4296);
  StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

    MassPoint mpisr0 (250, 170);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 1.08804));

    MassPoint mpisr1 (250, 180);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 1.08978));

    MassPoint mpisr2 (250, 190);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 1.08777));

    MassPoint mpisr3 (250, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 1.08765));

    MassPoint mpisr4 (250, 210);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 1.09132));

    MassPoint mpisr5 (250, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 1.09642));

    MassPoint mpisr6 (250, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 1.09773));

    MassPoint mpisr7 (250, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 1.09266));

    MassPoint mpisr8 (275, 195);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 1.09126));

    MassPoint mpisr9 (275, 205);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 1.0915));

    MassPoint mpisr10 (275, 215);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 1.09312));

    MassPoint mpisr11 (275, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 1.09188));

    MassPoint mpisr12 (275, 235);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 1.09601));

    MassPoint mpisr13 (275, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 1.09509));

    MassPoint mpisr14 (275, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 1.11015));

    MassPoint mpisr15 (275, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 1.09729));

    MassPoint mpisr16 (300, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 1.09549));

    MassPoint mpisr17 (300, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 1.09394));

    MassPoint mpisr18 (300, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 1.0951));

    MassPoint mpisr19 (300, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 1.09355));

    MassPoint mpisr20 (300, 260);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 1.09842));

    MassPoint mpisr21 (300, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 1.10406));

    MassPoint mpisr22 (300, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 1.10514));

    MassPoint mpisr23 (300, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 1.09951));

    MassPoint mpisr24 (325, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 1.0995));

    MassPoint mpisr25 (325, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 1.09865));

    MassPoint mpisr26 (325, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 1.09628));

    MassPoint mpisr27 (325, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 1.09828));

    MassPoint mpisr28 (325, 285);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 1.09882));

    MassPoint mpisr29 (325, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 1.10676));

    MassPoint mpisr30 (325, 305);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 1.09828));

    MassPoint mpisr31 (325, 315);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 1.10053));

    MassPoint mpisr32 (350, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 1.10384));

    MassPoint mpisr33 (350, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 1.10194));

    MassPoint mpisr34 (350, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 1.10312));

    MassPoint mpisr35 (350, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 1.10124));

    MassPoint mpisr36 (350, 310);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 1.10483));

    MassPoint mpisr37 (350, 320);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 1.10759));

    MassPoint mpisr38 (350, 330);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 1.10926));

    MassPoint mpisr39 (350, 340);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 1.09704));

    MassPoint mpisr40 (375, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 1.10555));

    MassPoint mpisr41 (375, 305);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 1.10516));

    MassPoint mpisr42 (375, 315);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 1.10435));

    MassPoint mpisr43 (375, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 1.1027));

    MassPoint mpisr44 (375, 335);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 1.10666));

    MassPoint mpisr45 (375, 345);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 1.10822));

    MassPoint mpisr46 (375, 355);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 1.10522));

    MassPoint mpisr47 (375, 365);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 1.101));

    MassPoint mpisr48 (400, 320);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 1.10678));

    MassPoint mpisr49 (400, 330);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 1.11024));

    MassPoint mpisr50 (400, 340);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 1.10686));

    MassPoint mpisr51 (400, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 1.10909));

    MassPoint mpisr52 (400, 360);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 1.10846));

    MassPoint mpisr53 (400, 370);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 1.10616));

    MassPoint mpisr54 (400, 380);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 1.1078));

    MassPoint mpisr55 (400, 390);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 1.09846));

    MassPoint mpisr56 (425, 345);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 1.11375));

    MassPoint mpisr57 (425, 355);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 1.10726));

    MassPoint mpisr58 (425, 365);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 1.11124));

    MassPoint mpisr59 (425, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 1.10937));

    MassPoint mpisr60 (425, 385);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 1.10536));

    MassPoint mpisr61 (425, 395);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 1.10899));

    MassPoint mpisr62 (425, 405);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 1.11036));

    MassPoint mpisr63 (425, 415);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 1.09619));

    MassPoint mpisr64 (450, 370);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 1.11356));

    MassPoint mpisr65 (450, 380);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 1.10915));

    MassPoint mpisr66 (450, 390);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 1.1077));

    MassPoint mpisr67 (450, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 1.11319));

    MassPoint mpisr68 (450, 410);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 1.11408));

    MassPoint mpisr69 (450, 420);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 1.11554));

    MassPoint mpisr70 (450, 430);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 1.10361));

    MassPoint mpisr71 (450, 440);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 1.09186));

    MassPoint mpisr72 (475, 395);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 1.11843));

    MassPoint mpisr73 (475, 405);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 1.11806));

    MassPoint mpisr74 (475, 415);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 1.11355));

    MassPoint mpisr75 (475, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 1.11386));

    MassPoint mpisr76 (475, 435);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 1.1104));

    MassPoint mpisr77 (475, 445);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 1.10262));

    MassPoint mpisr78 (475, 455);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 1.10949));

    MassPoint mpisr79 (475, 465);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 1.08611));

    MassPoint mpisr80 (500, 420);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 1.11964));

    MassPoint mpisr81 (500, 430);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 1.11323));

    MassPoint mpisr82 (500, 440);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 1.11675));

    MassPoint mpisr83 (500, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 1.11589));

    MassPoint mpisr84 (500, 460);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 1.12506));

    MassPoint mpisr85 (500, 470);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 1.11073));

    MassPoint mpisr86 (500, 480);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 1.12785));

    MassPoint mpisr87 (500, 490);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 1.09421));

    MassPoint mpisr88 (525, 445);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 1.12114));

    MassPoint mpisr89 (525, 455);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 1.11815));

    MassPoint mpisr90 (525, 465);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 1.11832));

    MassPoint mpisr91 (525, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 1.12254));

    MassPoint mpisr92 (525, 485);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 1.11582));

    MassPoint mpisr93 (525, 495);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 1.11543));

    MassPoint mpisr94 (525, 505);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 1.10735));

    MassPoint mpisr95 (525, 515);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 1.07599));

    MassPoint mpisr96 (550, 470);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 1.125));

    MassPoint mpisr97 (550, 480);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 1.11435));

    MassPoint mpisr98 (550, 490);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 1.11737));

    MassPoint mpisr99 (550, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 1.12473));

    MassPoint mpisr100 (550, 510);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 1.11266));

    MassPoint mpisr101 (550, 520);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 1.11193));

    MassPoint mpisr102 (550, 530);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 1.09595));

    MassPoint mpisr103 (550, 540);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 1.10684));

    MassPoint mpisr104 (575, 495);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 1.12715));

    MassPoint mpisr105 (575, 505);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 1.11361));

    MassPoint mpisr106 (575, 515);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 1.1216));

    MassPoint mpisr107 (575, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 1.11988));

    MassPoint mpisr108 (575, 535);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 1.12597));

    MassPoint mpisr109 (575, 545);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 1.11273));

    MassPoint mpisr110 (575, 555);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 1.08039));

    MassPoint mpisr111 (575, 565);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 1.0852));

    MassPoint mpisr112 (600, 520);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 1.1235));

    MassPoint mpisr113 (600, 530);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 1.11675));

    MassPoint mpisr114 (600, 540);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 1.118));

    MassPoint mpisr115 (600, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 1.13499));

    MassPoint mpisr116 (600, 560);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 1.1236));

    MassPoint mpisr117 (600, 570);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 1.13368));

    MassPoint mpisr118 (600, 580);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 1.12476));

    MassPoint mpisr119 (600, 590);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 1.10153));

    MassPoint mpisr120 (625, 545);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 1.1344));

    MassPoint mpisr121 (625, 555);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 1.12242));

    MassPoint mpisr122 (625, 565);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 1.12846));

    MassPoint mpisr123 (625, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 1.123));

    MassPoint mpisr124 (625, 585);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 1.1144));

    MassPoint mpisr125 (625, 595);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 1.11987));

    MassPoint mpisr126 (625, 605);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 1.07413));

    MassPoint mpisr127 (625, 615);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 1.11942));

    MassPoint mpisr128 (650, 570);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 1.13852));

    MassPoint mpisr129 (650, 580);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 1.13365));

    MassPoint mpisr130 (650, 590);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 1.12655));

    MassPoint mpisr131 (650, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 1.11674));

    MassPoint mpisr132 (650, 610);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 1.1215));

    MassPoint mpisr133 (650, 620);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 1.1188));

    MassPoint mpisr134 (650, 630);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 1.07199));

    MassPoint mpisr135 (650, 640);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 1.09518));

    MassPoint mpisr136 (675, 595);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 1.12964));

    MassPoint mpisr137 (675, 605);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 1.1174));

    MassPoint mpisr138 (675, 615);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 1.12533));

    MassPoint mpisr139 (675, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 1.11918));

    MassPoint mpisr140 (675, 635);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 1.12641));

    MassPoint mpisr141 (675, 645);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 1.13555));

    MassPoint mpisr142 (675, 655);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 1.14813));

    MassPoint mpisr143 (675, 665);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 1.10022));

    MassPoint mpisr144 (700, 620);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 1.13289));

    MassPoint mpisr145 (700, 630);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 1.13528));

    MassPoint mpisr146 (700, 640);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 1.11819));

    MassPoint mpisr147 (700, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 1.12001));

    MassPoint mpisr148 (700, 660);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 1.12116));

    MassPoint mpisr149 (700, 670);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 1.13309));

    MassPoint mpisr150 (700, 680);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 1.10314));

    MassPoint mpisr151 (700, 690);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 1.1197));

    MassPoint mpisr152 (725, 645);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 1.13072));

    MassPoint mpisr153 (725, 655);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 1.12496));

    MassPoint mpisr154 (725, 665);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 1.1363));

    MassPoint mpisr155 (725, 675);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 1.13207));

    MassPoint mpisr156 (725, 685);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 1.1498));

    MassPoint mpisr157 (725, 695);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 1.12125));

    MassPoint mpisr158 (725, 705);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 1.12742));

    MassPoint mpisr159 (725, 715);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 1.08159));

    MassPoint mpisr160 (750, 670);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 1.12547));

    MassPoint mpisr161 (750, 680);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 1.1256));

    MassPoint mpisr162 (750, 690);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 1.12983));

    MassPoint mpisr163 (750, 700);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 1.12136));

    MassPoint mpisr164 (750, 710);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 1.1403));

    MassPoint mpisr165 (750, 720);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 1.1143));

    MassPoint mpisr166 (750, 730);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 1.08711));

    MassPoint mpisr167 (750, 740);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 1.08817));

    MassPoint mpisr168 (775, 695);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 1.13168));

    MassPoint mpisr169 (775, 705);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 1.12821));

    MassPoint mpisr170 (775, 715);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 1.12534));

    MassPoint mpisr171 (775, 725);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 1.11686));

    MassPoint mpisr172 (775, 735);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 1.12553));

    MassPoint mpisr173 (775, 745);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 1.12157));

    MassPoint mpisr174 (775, 755);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 1.07711));

    MassPoint mpisr175 (775, 765);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 1.11223));

    MassPoint mpisr176 (800, 720);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 1.13877));

    MassPoint mpisr177 (800, 730);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 1.13434));

    MassPoint mpisr178 (800, 740);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 1.13937));

    MassPoint mpisr179 (800, 750);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 1.11316));

    MassPoint mpisr180 (800, 760);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 1.1286));

    MassPoint mpisr181 (800, 770);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 1.12292));

    MassPoint mpisr182 (800, 780);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 1.12931));

    MassPoint mpisr183 (800, 790);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 1.07682));

  } else if (SUSYProductionProcess=="T2bW_dM10to80") {

      MassPoint mp0 (250, 170);
  StopCrossSection cs0 (2.26682, 0.318703);
  MassPointParameters mpp0 (cs0, 100742);
  StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

  MassPoint mp1 (250, 180);
  StopCrossSection cs1 (2.26682, 0.318703);
  MassPointParameters mpp1 (cs1, 101467);
  StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

  MassPoint mp2 (250, 190);
  StopCrossSection cs2 (2.26682, 0.318703);
  MassPointParameters mpp2 (cs2, 92237);
  StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

  MassPoint mp3 (250, 200);
  StopCrossSection cs3 (2.26682, 0.318703);
  MassPointParameters mpp3 (cs3, 95963);
  StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

  MassPoint mp4 (250, 210);
  StopCrossSection cs4 (2.26682, 0.318703);
  MassPointParameters mpp4 (cs4, 96919);
  StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

  MassPoint mp5 (250, 220);
  StopCrossSection cs5 (2.26682, 0.318703);
  MassPointParameters mpp5 (cs5, 95434);
  StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

  MassPoint mp6 (250, 230);
  StopCrossSection cs6 (2.26682, 0.318703);
  MassPointParameters mpp6 (cs6, 94169);
  StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

  MassPoint mp7 (250, 240);
  StopCrossSection cs7 (2.26682, 0.318703);
  MassPointParameters mpp7 (cs7, 90539);
  StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

  MassPoint mp8 (275, 195);
  StopCrossSection cs8 (1.39853, 0.199358);
  MassPointParameters mpp8 (cs8, 87584);
  StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

  MassPoint mp9 (275, 205);
  StopCrossSection cs9 (1.39853, 0.199358);
  MassPointParameters mpp9 (cs9, 80180);
  StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

  MassPoint mp10 (275, 215);
  StopCrossSection cs10 (1.39853, 0.199358);
  MassPointParameters mpp10 (cs10, 92129);
  StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

  MassPoint mp11 (275, 225);
  StopCrossSection cs11 (1.39853, 0.199358);
  MassPointParameters mpp11 (cs11, 91877);
  StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

  MassPoint mp12 (275, 235);
  StopCrossSection cs12 (1.39853, 0.199358);
  MassPointParameters mpp12 (cs12, 82581);
  StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

  MassPoint mp13 (275, 245);
  StopCrossSection cs13 (1.39853, 0.199358);
  MassPointParameters mpp13 (cs13, 88021);
  StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

  MassPoint mp14 (275, 255);
  StopCrossSection cs14 (1.39853, 0.199358);
  MassPointParameters mpp14 (cs14, 86185);
  StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

  MassPoint mp15 (275, 265);
  StopCrossSection cs15 (1.39853, 0.199358);
  MassPointParameters mpp15 (cs15, 85306);
  StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

  MassPoint mp16 (300, 220);
  StopCrossSection cs16 (0.89394, 0.124457);
  MassPointParameters mpp16 (cs16, 94655);
  StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

  MassPoint mp17 (300, 230);
  StopCrossSection cs17 (0.89394, 0.124457);
  MassPointParameters mpp17 (cs17, 106028);
  StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

  MassPoint mp18 (300, 240);
  StopCrossSection cs18 (0.89394, 0.124457);
  MassPointParameters mpp18 (cs18, 102263);
  StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

  MassPoint mp19 (300, 250);
  StopCrossSection cs19 (0.89394, 0.124457);
  MassPointParameters mpp19 (cs19, 95677);
  StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

  MassPoint mp20 (300, 260);
  StopCrossSection cs20 (0.89394, 0.124457);
  MassPointParameters mpp20 (cs20, 98657);
  StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

  MassPoint mp21 (300, 270);
  StopCrossSection cs21 (0.89394, 0.124457);
  MassPointParameters mpp21 (cs21, 92277);
  StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

  MassPoint mp22 (300, 280);
  StopCrossSection cs22 (0.89394, 0.124457);
  MassPointParameters mpp22 (cs22, 95333);
  StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

  MassPoint mp23 (300, 290);
  StopCrossSection cs23 (0.89394, 0.124457);
  MassPointParameters mpp23 (cs23, 101376);
  StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

  MassPoint mp24 (325, 245);
  StopCrossSection cs24 (0.588326, 0.0812738);
  MassPointParameters mpp24 (cs24, 91491);
  StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

  MassPoint mp25 (325, 255);
  StopCrossSection cs25 (0.588326, 0.0812738);
  MassPointParameters mpp25 (cs25, 104371);
  StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

  MassPoint mp26 (325, 265);
  StopCrossSection cs26 (0.588326, 0.0812738);
  MassPointParameters mpp26 (cs26, 99250);
  StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

  MassPoint mp27 (325, 275);
  StopCrossSection cs27 (0.588326, 0.0812738);
  MassPointParameters mpp27 (cs27, 94450);
  StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

  MassPoint mp28 (325, 285);
  StopCrossSection cs28 (0.588326, 0.0812738);
  MassPointParameters mpp28 (cs28, 97836);
  StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

  MassPoint mp29 (325, 295);
  StopCrossSection cs29 (0.588326, 0.0812738);
  MassPointParameters mpp29 (cs29, 89118);
  StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

  MassPoint mp30 (325, 305);
  StopCrossSection cs30 (0.588326, 0.0812738);
  MassPointParameters mpp30 (cs30, 94327);
  StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

  MassPoint mp31 (325, 315);
  StopCrossSection cs31 (0.588326, 0.0812738);
  MassPointParameters mpp31 (cs31, 94698);
  StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

  MassPoint mp32 (350, 270);
  StopCrossSection cs32 (0.39748, 0.0544059);
  MassPointParameters mpp32 (cs32, 98516);
  StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

  MassPoint mp33 (350, 280);
  StopCrossSection cs33 (0.39748, 0.0544059);
  MassPointParameters mpp33 (cs33, 94678);
  StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

  MassPoint mp34 (350, 290);
  StopCrossSection cs34 (0.39748, 0.0544059);
  MassPointParameters mpp34 (cs34, 92260);
  StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

  MassPoint mp35 (350, 300);
  StopCrossSection cs35 (0.39748, 0.0544059);
  MassPointParameters mpp35 (cs35, 90474);
  StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

  MassPoint mp36 (350, 310);
  StopCrossSection cs36 (0.39748, 0.0544059);
  MassPointParameters mpp36 (cs36, 92340);
  StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

  MassPoint mp37 (350, 320);
  StopCrossSection cs37 (0.39748, 0.0544059);
  MassPointParameters mpp37 (cs37, 97778);
  StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

  MassPoint mp38 (350, 330);
  StopCrossSection cs38 (0.39748, 0.0544059);
  MassPointParameters mpp38 (cs38, 93123);
  StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

  MassPoint mp39 (350, 340);
  StopCrossSection cs39 (0.39748, 0.0544059);
  MassPointParameters mpp39 (cs39, 94925);
  StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

  MassPoint mp40 (375, 295);
  StopCrossSection cs40 (0.274142, 0.0379623);
  MassPointParameters mpp40 (cs40, 89326);
  StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

  MassPoint mp41 (375, 305);
  StopCrossSection cs41 (0.274142, 0.0379623);
  MassPointParameters mpp41 (cs41, 93868);
  StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

  MassPoint mp42 (375, 315);
  StopCrossSection cs42 (0.274142, 0.0379623);
  MassPointParameters mpp42 (cs42, 90673);
  StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

  MassPoint mp43 (375, 325);
  StopCrossSection cs43 (0.274142, 0.0379623);
  MassPointParameters mpp43 (cs43, 91507);
  StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

  MassPoint mp44 (375, 335);
  StopCrossSection cs44 (0.274142, 0.0379623);
  MassPointParameters mpp44 (cs44, 93793);
  StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

  MassPoint mp45 (375, 345);
  StopCrossSection cs45 (0.274142, 0.0379623);
  MassPointParameters mpp45 (cs45, 91700);
  StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

  MassPoint mp46 (375, 355);
  StopCrossSection cs46 (0.274142, 0.0379623);
  MassPointParameters mpp46 (cs46, 93032);
  StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

  MassPoint mp47 (375, 365);
  StopCrossSection cs47 (0.274142, 0.0379623);
  MassPointParameters mpp47 (cs47, 96022);
  StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

  MassPoint mp48 (400, 320);
  StopCrossSection cs48 (0.192659, 0.0263914);
  MassPointParameters mpp48 (cs48, 67175);
  StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

  MassPoint mp49 (400, 330);
  StopCrossSection cs49 (0.192659, 0.0263914);
  MassPointParameters mpp49 (cs49, 76298);
  StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

  MassPoint mp50 (400, 340);
  StopCrossSection cs50 (0.192659, 0.0263914);
  MassPointParameters mpp50 (cs50, 72466);
  StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

  MassPoint mp51 (400, 350);
  StopCrossSection cs51 (0.192659, 0.0263914);
  MassPointParameters mpp51 (cs51, 71548);
  StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

  MassPoint mp52 (400, 360);
  StopCrossSection cs52 (0.192659, 0.0263914);
  MassPointParameters mpp52 (cs52, 69838);
  StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

  MassPoint mp53 (400, 370);
  StopCrossSection cs53 (0.192659, 0.0263914);
  MassPointParameters mpp53 (cs53, 76312);
  StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

  MassPoint mp54 (400, 380);
  StopCrossSection cs54 (0.192659, 0.0263914);
  MassPointParameters mpp54 (cs54, 73463);
  StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

  MassPoint mp55 (400, 390);
  StopCrossSection cs55 (0.192659, 0.0263914);
  MassPointParameters mpp55 (cs55, 75725);
  StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

  MassPoint mp56 (425, 345);
  StopCrossSection cs56 (0.137688, 0.0185897);
  MassPointParameters mpp56 (cs56, 52271);
  StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

  MassPoint mp57 (425, 355);
  StopCrossSection cs57 (0.137688, 0.0185897);
  MassPointParameters mpp57 (cs57, 51883);
  StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

  MassPoint mp58 (425, 365);
  StopCrossSection cs58 (0.137688, 0.0185897);
  MassPointParameters mpp58 (cs58, 43879);
  StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

  MassPoint mp59 (425, 375);
  StopCrossSection cs59 (0.137688, 0.0185897);
  MassPointParameters mpp59 (cs59, 57487);
  StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

  MassPoint mp60 (425, 385);
  StopCrossSection cs60 (0.137688, 0.0185897);
  MassPointParameters mpp60 (cs60, 50116);
  StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

  MassPoint mp61 (425, 395);
  StopCrossSection cs61 (0.137688, 0.0185897);
  MassPointParameters mpp61 (cs61, 53885);
  StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

  MassPoint mp62 (425, 405);
  StopCrossSection cs62 (0.137688, 0.0185897);
  MassPointParameters mpp62 (cs62, 56107);
  StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

  MassPoint mp63 (425, 415);
  StopCrossSection cs63 (0.137688, 0.0185897);
  MassPointParameters mpp63 (cs63, 51452);
  StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

  MassPoint mp64 (450, 370);
  StopCrossSection cs64 (0.0995465, 0.0133949);
  MassPointParameters mpp64 (cs64, 37084);
  StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

  MassPoint mp65 (450, 380);
  StopCrossSection cs65 (0.0995465, 0.0133949);
  MassPointParameters mpp65 (cs65, 36237);
  StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

  MassPoint mp66 (450, 390);
  StopCrossSection cs66 (0.0995465, 0.0133949);
  MassPointParameters mpp66 (cs66, 37626);
  StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

  MassPoint mp67 (450, 400);
  StopCrossSection cs67 (0.0995465, 0.0133949);
  MassPointParameters mpp67 (cs67, 30615);
  StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

  MassPoint mp68 (450, 410);
  StopCrossSection cs68 (0.0995465, 0.0133949);
  MassPointParameters mpp68 (cs68, 36564);
  StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

  MassPoint mp69 (450, 420);
  StopCrossSection cs69 (0.0995465, 0.0133949);
  MassPointParameters mpp69 (cs69, 38685);
  StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

  MassPoint mp70 (450, 430);
  StopCrossSection cs70 (0.0995465, 0.0133949);
  MassPointParameters mpp70 (cs70, 35981);
  StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

  MassPoint mp71 (450, 440);
  StopCrossSection cs71 (0.0995465, 0.0133949);
  MassPointParameters mpp71 (cs71, 30664);
  StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

  MassPoint mp72 (475, 395);
  StopCrossSection cs72 (0.073172, 0.00979963);
  MassPointParameters mpp72 (cs72, 25947);
  StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

  MassPoint mp73 (475, 405);
  StopCrossSection cs73 (0.073172, 0.00979963);
  MassPointParameters mpp73 (cs73, 19095);
  StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

  MassPoint mp74 (475, 415);
  StopCrossSection cs74 (0.073172, 0.00979963);
  MassPointParameters mpp74 (cs74, 26414);
  StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

  MassPoint mp75 (475, 425);
  StopCrossSection cs75 (0.073172, 0.00979963);
  MassPointParameters mpp75 (cs75, 29578);
  StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

  MassPoint mp76 (475, 435);
  StopCrossSection cs76 (0.073172, 0.00979963);
  MassPointParameters mpp76 (cs76, 26269);
  StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

  MassPoint mp77 (475, 445);
  StopCrossSection cs77 (0.073172, 0.00979963);
  MassPointParameters mpp77 (cs77, 25370);
  StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

  MassPoint mp78 (475, 455);
  StopCrossSection cs78 (0.073172, 0.00979963);
  MassPointParameters mpp78 (cs78, 26254);
  StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

  MassPoint mp79 (475, 465);
  StopCrossSection cs79 (0.073172, 0.00979963);
  MassPointParameters mpp79 (cs79, 24310);
  StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

  MassPoint mp80 (500, 420);
  StopCrossSection cs80 (0.0544248, 0.00728188);
  MassPointParameters mpp80 (cs80, 20681);
  StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

  MassPoint mp81 (500, 430);
  StopCrossSection cs81 (0.0544248, 0.00728188);
  MassPointParameters mpp81 (cs81, 18965);
  StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

  MassPoint mp82 (500, 440);
  StopCrossSection cs82 (0.0544248, 0.00728188);
  MassPointParameters mpp82 (cs82, 17933);
  StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

  MassPoint mp83 (500, 450);
  StopCrossSection cs83 (0.0544248, 0.00728188);
  MassPointParameters mpp83 (cs83, 20755);
  StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

  MassPoint mp84 (500, 460);
  StopCrossSection cs84 (0.0544248, 0.00728188);
  MassPointParameters mpp84 (cs84, 19930);
  StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

  MassPoint mp85 (500, 470);
  StopCrossSection cs85 (0.0544248, 0.00728188);
  MassPointParameters mpp85 (cs85, 21938);
  StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

  MassPoint mp86 (500, 480);
  StopCrossSection cs86 (0.0544248, 0.00728188);
  MassPointParameters mpp86 (cs86, 20529);
  StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

  MassPoint mp87 (500, 490);
  StopCrossSection cs87 (0.0544248, 0.00728188);
  MassPointParameters mpp87 (cs87, 20851);
  StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

  MassPoint mp88 (525, 445);
  StopCrossSection cs88 (0.0409701, 0.00546717);
  MassPointParameters mpp88 (cs88, 17222);
  StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

  MassPoint mp89 (525, 455);
  StopCrossSection cs89 (0.0409701, 0.00546717);
  MassPointParameters mpp89 (cs89, 17456);
  StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

  MassPoint mp90 (525, 465);
  StopCrossSection cs90 (0.0409701, 0.00546717);
  MassPointParameters mpp90 (cs90, 16765);
  StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

  MassPoint mp91 (525, 475);
  StopCrossSection cs91 (0.0409701, 0.00546717);
  MassPointParameters mpp91 (cs91, 15376);
  StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

  MassPoint mp92 (525, 485);
  StopCrossSection cs92 (0.0409701, 0.00546717);
  MassPointParameters mpp92 (cs92, 15990);
  StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

  MassPoint mp93 (525, 495);
  StopCrossSection cs93 (0.0409701, 0.00546717);
  MassPointParameters mpp93 (cs93, 15970);
  StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

  MassPoint mp94 (525, 505);
  StopCrossSection cs94 (0.0409701, 0.00546717);
  MassPointParameters mpp94 (cs94, 13708);
  StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

  MassPoint mp95 (525, 515);
  StopCrossSection cs95 (0.0409701, 0.00546717);
  MassPointParameters mpp95 (cs95, 14555);
  StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

  MassPoint mp96 (550, 470);
  StopCrossSection cs96 (0.0310846, 0.00412452);
  MassPointParameters mpp96 (cs96, 14093);
  StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

  MassPoint mp97 (550, 480);
  StopCrossSection cs97 (0.0310846, 0.00412452);
  MassPointParameters mpp97 (cs97, 12448);
  StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

  MassPoint mp98 (550, 490);
  StopCrossSection cs98 (0.0310846, 0.00412452);
  MassPointParameters mpp98 (cs98, 11235);
  StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

  MassPoint mp99 (550, 500);
  StopCrossSection cs99 (0.0310846, 0.00412452);
  MassPointParameters mpp99 (cs99, 11669);
  StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

  MassPoint mp100 (550, 510);
  StopCrossSection cs100 (0.0310846, 0.00412452);
  MassPointParameters mpp100 (cs100, 10863);
  StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

  MassPoint mp101 (550, 520);
  StopCrossSection cs101 (0.0310846, 0.00412452);
  MassPointParameters mpp101 (cs101, 12237);
  StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

  MassPoint mp102 (550, 530);
  StopCrossSection cs102 (0.0310846, 0.00412452);
  MassPointParameters mpp102 (cs102, 12750);
  StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

  MassPoint mp103 (550, 540);
  StopCrossSection cs103 (0.0310846, 0.00412452);
  MassPointParameters mpp103 (cs103, 11141);
  StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

  MassPoint mp104 (575, 495);
  StopCrossSection cs104 (0.0237356, 0.00315069);
  MassPointParameters mpp104 (cs104, 6806);
  StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

  MassPoint mp105 (575, 505);
  StopCrossSection cs105 (0.0237356, 0.00315069);
  MassPointParameters mpp105 (cs105, 9317);
  StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

  MassPoint mp106 (575, 515);
  StopCrossSection cs106 (0.0237356, 0.00315069);
  MassPointParameters mpp106 (cs106, 9978);
  StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

  MassPoint mp107 (575, 525);
  StopCrossSection cs107 (0.0237356, 0.00315069);
  MassPointParameters mpp107 (cs107, 7403);
  StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

  MassPoint mp108 (575, 535);
  StopCrossSection cs108 (0.0237356, 0.00315069);
  MassPointParameters mpp108 (cs108, 9209);
  StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

  MassPoint mp109 (575, 545);
  StopCrossSection cs109 (0.0237356, 0.00315069);
  MassPointParameters mpp109 (cs109, 7472);
  StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

  MassPoint mp110 (575, 555);
  StopCrossSection cs110 (0.0237356, 0.00315069);
  MassPointParameters mpp110 (cs110, 10148);
  StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

  MassPoint mp111 (575, 565);
  StopCrossSection cs111 (0.0237356, 0.00315069);
  MassPointParameters mpp111 (cs111, 8846);
  StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

  MassPoint mp112 (600, 520);
  StopCrossSection cs112 (0.0183277, 0.00242061);
  MassPointParameters mpp112 (cs112, 6731);
  StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

  MassPoint mp113 (600, 530);
  StopCrossSection cs113 (0.0183277, 0.00242061);
  MassPointParameters mpp113 (cs113, 7373);
  StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

  MassPoint mp114 (600, 540);
  StopCrossSection cs114 (0.0183277, 0.00242061);
  MassPointParameters mpp114 (cs114, 7015);
  StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

  MassPoint mp115 (600, 550);
  StopCrossSection cs115 (0.0183277, 0.00242061);
  MassPointParameters mpp115 (cs115, 9219);
  StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

  MassPoint mp116 (600, 560);
  StopCrossSection cs116 (0.0183277, 0.00242061);
  MassPointParameters mpp116 (cs116, 6575);
  StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

  MassPoint mp117 (600, 570);
  StopCrossSection cs117 (0.0183277, 0.00242061);
  MassPointParameters mpp117 (cs117, 8442);
  StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

  MassPoint mp118 (600, 580);
  StopCrossSection cs118 (0.0183277, 0.00242061);
  MassPointParameters mpp118 (cs118, 7737);
  StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

  MassPoint mp119 (600, 590);
  StopCrossSection cs119 (0.0183277, 0.00242061);
  MassPointParameters mpp119 (cs119, 7306);
  StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

  MassPoint mp120 (625, 545);
  StopCrossSection cs120 (0.014315, 0.00183177);
  MassPointParameters mpp120 (cs120, 6206);
  StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

  MassPoint mp121 (625, 555);
  StopCrossSection cs121 (0.014315, 0.00183177);
  MassPointParameters mpp121 (cs121, 6296);
  StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

  MassPoint mp122 (625, 565);
  StopCrossSection cs122 (0.014315, 0.00183177);
  MassPointParameters mpp122 (cs122, 5107);
  StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

  MassPoint mp123 (625, 575);
  StopCrossSection cs123 (0.014315, 0.00183177);
  MassPointParameters mpp123 (cs123, 5529);
  StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

  MassPoint mp124 (625, 585);
  StopCrossSection cs124 (0.014315, 0.00183177);
  MassPointParameters mpp124 (cs124, 6950);
  StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

  MassPoint mp125 (625, 595);
  StopCrossSection cs125 (0.014315, 0.00183177);
  MassPointParameters mpp125 (cs125, 4196);
  StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

  MassPoint mp126 (625, 605);
  StopCrossSection cs126 (0.014315, 0.00183177);
  MassPointParameters mpp126 (cs126, 5199);
  StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

  MassPoint mp127 (625, 615);
  StopCrossSection cs127 (0.014315, 0.00183177);
  MassPointParameters mpp127 (cs127, 4751);
  StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

  MassPoint mp128 (650, 570);
  StopCrossSection cs128 (0.0112365, 0.00145212);
  MassPointParameters mpp128 (cs128, 3964);
  StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

  MassPoint mp129 (650, 580);
  StopCrossSection cs129 (0.0112365, 0.00145212);
  MassPointParameters mpp129 (cs129, 4435);
  StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

  MassPoint mp130 (650, 590);
  StopCrossSection cs130 (0.0112365, 0.00145212);
  MassPointParameters mpp130 (cs130, 4254);
  StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

  MassPoint mp131 (650, 600);
  StopCrossSection cs131 (0.0112365, 0.00145212);
  MassPointParameters mpp131 (cs131, 3923);
  StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

  MassPoint mp132 (650, 610);
  StopCrossSection cs132 (0.0112365, 0.00145212);
  MassPointParameters mpp132 (cs132, 3629);
  StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

  MassPoint mp133 (650, 620);
  StopCrossSection cs133 (0.0112365, 0.00145212);
  MassPointParameters mpp133 (cs133, 3548);
  StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

  MassPoint mp134 (650, 630);
  StopCrossSection cs134 (0.0112365, 0.00145212);
  MassPointParameters mpp134 (cs134, 5131);
  StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

  MassPoint mp135 (650, 640);
  StopCrossSection cs135 (0.0112365, 0.00145212);
  MassPointParameters mpp135 (cs135, 4457);
  StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

  MassPoint mp136 (675, 595);
  StopCrossSection cs136 (0.00886867, 0.00115916);
  MassPointParameters mpp136 (cs136, 4555);
  StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

  MassPoint mp137 (675, 605);
  StopCrossSection cs137 (0.00886867, 0.00115916);
  MassPointParameters mpp137 (cs137, 3788);
  StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

  MassPoint mp138 (675, 615);
  StopCrossSection cs138 (0.00886867, 0.00115916);
  MassPointParameters mpp138 (cs138, 4948);
  StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

  MassPoint mp139 (675, 625);
  StopCrossSection cs139 (0.00886867, 0.00115916);
  MassPointParameters mpp139 (cs139, 3439);
  StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

  MassPoint mp140 (675, 635);
  StopCrossSection cs140 (0.00886867, 0.00115916);
  MassPointParameters mpp140 (cs140, 3527);
  StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

  MassPoint mp141 (675, 645);
  StopCrossSection cs141 (0.00886867, 0.00115916);
  MassPointParameters mpp141 (cs141, 3568);
  StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

  MassPoint mp142 (675, 655);
  StopCrossSection cs142 (0.00886867, 0.00115916);
  MassPointParameters mpp142 (cs142, 2435);
  StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

  MassPoint mp143 (675, 665);
  StopCrossSection cs143 (0.00886867, 0.00115916);
  MassPointParameters mpp143 (cs143, 4493);
  StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

  MassPoint mp144 (700, 620);
  StopCrossSection cs144 (0.00703799, 0.000939072);
  MassPointParameters mpp144 (cs144, 3909);
  StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

  MassPoint mp145 (700, 630);
  StopCrossSection cs145 (0.00703799, 0.000939072);
  MassPointParameters mpp145 (cs145, 3415);
  StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

  MassPoint mp146 (700, 640);
  StopCrossSection cs146 (0.00703799, 0.000939072);
  MassPointParameters mpp146 (cs146, 4365);
  StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

  MassPoint mp147 (700, 650);
  StopCrossSection cs147 (0.00703799, 0.000939072);
  MassPointParameters mpp147 (cs147, 2652);
  StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

  MassPoint mp148 (700, 660);
  StopCrossSection cs148 (0.00703799, 0.000939072);
  MassPointParameters mpp148 (cs148, 5575);
  StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

  MassPoint mp149 (700, 670);
  StopCrossSection cs149 (0.00703799, 0.000939072);
  MassPointParameters mpp149 (cs149, 4830);
  StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

  MassPoint mp150 (700, 680);
  StopCrossSection cs150 (0.00703799, 0.000939072);
  MassPointParameters mpp150 (cs150, 3316);
  StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

  MassPoint mp151 (700, 690);
  StopCrossSection cs151 (0.00703799, 0.000939072);
  MassPointParameters mpp151 (cs151, 4901);
  StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

  MassPoint mp152 (725, 645);
  StopCrossSection cs152 (0.00563099, 0.000764711);
  MassPointParameters mpp152 (cs152, 3097);
  StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

  MassPoint mp153 (725, 655);
  StopCrossSection cs153 (0.00563099, 0.000764711);
  MassPointParameters mpp153 (cs153, 4584);
  StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

  MassPoint mp154 (725, 665);
  StopCrossSection cs154 (0.00563099, 0.000764711);
  MassPointParameters mpp154 (cs154, 2932);
  StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

  MassPoint mp155 (725, 675);
  StopCrossSection cs155 (0.00563099, 0.000764711);
  MassPointParameters mpp155 (cs155, 3338);
  StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

  MassPoint mp156 (725, 685);
  StopCrossSection cs156 (0.00563099, 0.000764711);
  MassPointParameters mpp156 (cs156, 3832);
  StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

  MassPoint mp157 (725, 695);
  StopCrossSection cs157 (0.00563099, 0.000764711);
  MassPointParameters mpp157 (cs157, 2472);
  StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

  MassPoint mp158 (725, 705);
  StopCrossSection cs158 (0.00563099, 0.000764711);
  MassPointParameters mpp158 (cs158, 4703);
  StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

  MassPoint mp159 (725, 715);
  StopCrossSection cs159 (0.00563099, 0.000764711);
  MassPointParameters mpp159 (cs159, 4955);
  StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

  MassPoint mp160 (750, 670);
  StopCrossSection cs160 (0.00452859, 0.000622478);
  MassPointParameters mpp160 (cs160, 3618);
  StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

  MassPoint mp161 (750, 680);
  StopCrossSection cs161 (0.00452859, 0.000622478);
  MassPointParameters mpp161 (cs161, 3754);
  StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

  MassPoint mp162 (750, 690);
  StopCrossSection cs162 (0.00452859, 0.000622478);
  MassPointParameters mpp162 (cs162, 5450);
  StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

  MassPoint mp163 (750, 700);
  StopCrossSection cs163 (0.00452859, 0.000622478);
  MassPointParameters mpp163 (cs163, 4344);
  StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

  MassPoint mp164 (750, 710);
  StopCrossSection cs164 (0.00452859, 0.000622478);
  MassPointParameters mpp164 (cs164, 2952);
  StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

  MassPoint mp165 (750, 720);
  StopCrossSection cs165 (0.00452859, 0.000622478);
  MassPointParameters mpp165 (cs165, 4611);
  StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

  MassPoint mp166 (750, 730);
  StopCrossSection cs166 (0.00452859, 0.000622478);
  MassPointParameters mpp166 (cs166, 5604);
  StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

  MassPoint mp167 (750, 740);
  StopCrossSection cs167 (0.00452859, 0.000622478);
  MassPointParameters mpp167 (cs167, 4525);
  StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

  MassPoint mp168 (775, 695);
  StopCrossSection cs168 (0.00366131, 0.000511108);
  MassPointParameters mpp168 (cs168, 3559);
  StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

  MassPoint mp169 (775, 705);
  StopCrossSection cs169 (0.00366131, 0.000511108);
  MassPointParameters mpp169 (cs169, 3996);
  StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

  MassPoint mp170 (775, 715);
  StopCrossSection cs170 (0.00366131, 0.000511108);
  MassPointParameters mpp170 (cs170, 2943);
  StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

  MassPoint mp171 (775, 725);
  StopCrossSection cs171 (0.00366131, 0.000511108);
  MassPointParameters mpp171 (cs171, 5219);
  StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

  MassPoint mp172 (775, 735);
  StopCrossSection cs172 (0.00366131, 0.000511108);
  MassPointParameters mpp172 (cs172, 3012);
  StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

  MassPoint mp173 (775, 745);
  StopCrossSection cs173 (0.00366131, 0.000511108);
  MassPointParameters mpp173 (cs173, 3844);
  StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

  MassPoint mp174 (775, 755);
  StopCrossSection cs174 (0.00366131, 0.000511108);
  MassPointParameters mpp174 (cs174, 3763);
  StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

  MassPoint mp175 (775, 765);
  StopCrossSection cs175 (0.00366131, 0.000511108);
  MassPointParameters mpp175 (cs175, 4145);
  StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

  MassPoint mp176 (800, 720);
  StopCrossSection cs176 (0.0029742, 0.000421474);
  MassPointParameters mpp176 (cs176, 3834);
  StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

  MassPoint mp177 (800, 730);
  StopCrossSection cs177 (0.0029742, 0.000421474);
  MassPointParameters mpp177 (cs177, 3488);
  StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

  MassPoint mp178 (800, 740);
  StopCrossSection cs178 (0.0029742, 0.000421474);
  MassPointParameters mpp178 (cs178, 5081);
  StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

  MassPoint mp179 (800, 750);
  StopCrossSection cs179 (0.0029742, 0.000421474);
  MassPointParameters mpp179 (cs179, 3928);
  StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

  MassPoint mp180 (800, 760);
  StopCrossSection cs180 (0.0029742, 0.000421474);
  MassPointParameters mpp180 (cs180, 4641);
  StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

  MassPoint mp181 (800, 770);
  StopCrossSection cs181 (0.0029742, 0.000421474);
  MassPointParameters mpp181 (cs181, 3199);
  StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

  MassPoint mp182 (800, 780);
  StopCrossSection cs182 (0.0029742, 0.000421474);
  MassPointParameters mpp182 (cs182, 3580);
  StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

  MassPoint mp183 (800, 790);
  StopCrossSection cs183 (0.0029742, 0.000421474);
  MassPointParameters mpp183 (cs183, 3312);
  StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

    MassPoint mpisr0 (250, 170);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 1.08158));

    MassPoint mpisr1 (250, 180);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 1.08166));

    MassPoint mpisr2 (250, 190);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 1.08566));

    MassPoint mpisr3 (250, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 1.08455));

    MassPoint mpisr4 (250, 210);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 1.08754));

    MassPoint mpisr5 (250, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 1.09412));

    MassPoint mpisr6 (250, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 1.10574));

    MassPoint mpisr7 (250, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 1.09891));

    MassPoint mpisr8 (275, 195);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 1.08508));

    MassPoint mpisr9 (275, 205);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 1.08633));

    MassPoint mpisr10 (275, 215);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 1.08656));

    MassPoint mpisr11 (275, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 1.09132));

    MassPoint mpisr12 (275, 235);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 1.09632));

    MassPoint mpisr13 (275, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 1.0955));

    MassPoint mpisr14 (275, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 1.09997));

    MassPoint mpisr15 (275, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 1.09354));

    MassPoint mpisr16 (300, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 1.08542));

    MassPoint mpisr17 (300, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 1.08676));

    MassPoint mpisr18 (300, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 1.09143));

    MassPoint mpisr19 (300, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 1.09106));

    MassPoint mpisr20 (300, 260);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 1.0953));

    MassPoint mpisr21 (300, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 1.1019));

    MassPoint mpisr22 (300, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 1.10497));

    MassPoint mpisr23 (300, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 1.09336));

    MassPoint mpisr24 (325, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 1.08724));

    MassPoint mpisr25 (325, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 1.09022));

    MassPoint mpisr26 (325, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 1.09166));

    MassPoint mpisr27 (325, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 1.09687));

    MassPoint mpisr28 (325, 285);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 1.10077));

    MassPoint mpisr29 (325, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 1.10821));

    MassPoint mpisr30 (325, 305);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 1.10459));

    MassPoint mpisr31 (325, 315);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 1.09395));

    MassPoint mpisr32 (350, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 1.09442));

    MassPoint mpisr33 (350, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 1.0922));

    MassPoint mpisr34 (350, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 1.09587));

    MassPoint mpisr35 (350, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 1.09387));

    MassPoint mpisr36 (350, 310);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 1.10222));

    MassPoint mpisr37 (350, 320);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 1.10323));

    MassPoint mpisr38 (350, 330);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 1.11104));

    MassPoint mpisr39 (350, 340);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 1.09708));

    MassPoint mpisr40 (375, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 1.09635));

    MassPoint mpisr41 (375, 305);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 1.09766));

    MassPoint mpisr42 (375, 315);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 1.09889));

    MassPoint mpisr43 (375, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 1.10073));

    MassPoint mpisr44 (375, 335);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 1.09887));

    MassPoint mpisr45 (375, 345);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 1.10888));

    MassPoint mpisr46 (375, 355);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 1.10516));

    MassPoint mpisr47 (375, 365);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 1.10858));

    MassPoint mpisr48 (400, 320);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 1.09969));

    MassPoint mpisr49 (400, 330);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 1.09751));

    MassPoint mpisr50 (400, 340);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 1.09988));

    MassPoint mpisr51 (400, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 1.10064));

    MassPoint mpisr52 (400, 360);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 1.1111));

    MassPoint mpisr53 (400, 370);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 1.10501));

    MassPoint mpisr54 (400, 380);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 1.1088));

    MassPoint mpisr55 (400, 390);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 1.0916));

    MassPoint mpisr56 (425, 345);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 1.10181));

    MassPoint mpisr57 (425, 355);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 1.10152));

    MassPoint mpisr58 (425, 365);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 1.1015));

    MassPoint mpisr59 (425, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 1.10349));

    MassPoint mpisr60 (425, 385);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 1.1028));

    MassPoint mpisr61 (425, 395);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 1.10537));

    MassPoint mpisr62 (425, 405);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 1.10482));

    MassPoint mpisr63 (425, 415);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 1.10031));

    MassPoint mpisr64 (450, 370);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 1.1013));

    MassPoint mpisr65 (450, 380);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 1.10133));

    MassPoint mpisr66 (450, 390);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 1.09837));

    MassPoint mpisr67 (450, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 1.10554));

    MassPoint mpisr68 (450, 410);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 1.11908));

    MassPoint mpisr69 (450, 420);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 1.10628));

    MassPoint mpisr70 (450, 430);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 1.10637));

    MassPoint mpisr71 (450, 440);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 1.09586));

    MassPoint mpisr72 (475, 395);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 1.10044));

    MassPoint mpisr73 (475, 405);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 1.1011));

    MassPoint mpisr74 (475, 415);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 1.10867));

    MassPoint mpisr75 (475, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 1.11382));

    MassPoint mpisr76 (475, 435);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 1.11314));

    MassPoint mpisr77 (475, 445);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 1.10462));

    MassPoint mpisr78 (475, 455);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 1.10112));

    MassPoint mpisr79 (475, 465);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 1.10307));

    MassPoint mpisr80 (500, 420);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 1.10391));

    MassPoint mpisr81 (500, 430);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 1.10723));

    MassPoint mpisr82 (500, 440);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 1.10087));

    MassPoint mpisr83 (500, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 1.10421));

    MassPoint mpisr84 (500, 460);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 1.11266));

    MassPoint mpisr85 (500, 470);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 1.11778));

    MassPoint mpisr86 (500, 480);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 1.11041));

    MassPoint mpisr87 (500, 490);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 1.0895));

    MassPoint mpisr88 (525, 445);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 1.10493));

    MassPoint mpisr89 (525, 455);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 1.10485));

    MassPoint mpisr90 (525, 465);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 1.10811));

    MassPoint mpisr91 (525, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 1.12256));

    MassPoint mpisr92 (525, 485);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 1.11394));

    MassPoint mpisr93 (525, 495);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 1.11617));

    MassPoint mpisr94 (525, 505);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 1.10873));

    MassPoint mpisr95 (525, 515);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 1.11315));

    MassPoint mpisr96 (550, 470);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 1.11107));

    MassPoint mpisr97 (550, 480);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 1.11389));

    MassPoint mpisr98 (550, 490);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 1.11121));

    MassPoint mpisr99 (550, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 1.11533));

    MassPoint mpisr100 (550, 510);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 1.09908));

    MassPoint mpisr101 (550, 520);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 1.1086));

    MassPoint mpisr102 (550, 530);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 1.11236));

    MassPoint mpisr103 (550, 540);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 1.11776));

    MassPoint mpisr104 (575, 495);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 1.103));

    MassPoint mpisr105 (575, 505);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 1.10957));

    MassPoint mpisr106 (575, 515);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 1.1142));

    MassPoint mpisr107 (575, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 1.10931));

    MassPoint mpisr108 (575, 535);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 1.09947));

    MassPoint mpisr109 (575, 545);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 1.11581));

    MassPoint mpisr110 (575, 555);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 1.13398));

    MassPoint mpisr111 (575, 565);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 1.10856));

    MassPoint mpisr112 (600, 520);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 1.11167));

    MassPoint mpisr113 (600, 530);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 1.11947));

    MassPoint mpisr114 (600, 540);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 1.12451));

    MassPoint mpisr115 (600, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 1.10787));

    MassPoint mpisr116 (600, 560);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 1.12015));

    MassPoint mpisr117 (600, 570);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 1.09304));

    MassPoint mpisr118 (600, 580);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 1.12323));

    MassPoint mpisr119 (600, 590);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 1.10676));

    MassPoint mpisr120 (625, 545);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 1.10411));

    MassPoint mpisr121 (625, 555);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 1.11528));

    MassPoint mpisr122 (625, 565);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 1.11709));

    MassPoint mpisr123 (625, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 1.10349));

    MassPoint mpisr124 (625, 585);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 1.13929));

    MassPoint mpisr125 (625, 595);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 1.10873));

    MassPoint mpisr126 (625, 605);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 1.09062));

    MassPoint mpisr127 (625, 615);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 1.10637));

    MassPoint mpisr128 (650, 570);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 1.1359));

    MassPoint mpisr129 (650, 580);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 1.11268));

    MassPoint mpisr130 (650, 590);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 1.10882));

    MassPoint mpisr131 (650, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 1.13406));

    MassPoint mpisr132 (650, 610);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 1.11648));

    MassPoint mpisr133 (650, 620);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 1.10802));

    MassPoint mpisr134 (650, 630);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 1.09903));

    MassPoint mpisr135 (650, 640);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 1.10213));

    MassPoint mpisr136 (675, 595);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 1.11537));

    MassPoint mpisr137 (675, 605);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 1.11534));

    MassPoint mpisr138 (675, 615);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 1.11263));

    MassPoint mpisr139 (675, 625);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 1.1053));

    MassPoint mpisr140 (675, 635);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 1.11221));

    MassPoint mpisr141 (675, 645);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 1.11569));

    MassPoint mpisr142 (675, 655);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 1.10957));

    MassPoint mpisr143 (675, 665);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 1.09085));

    MassPoint mpisr144 (700, 620);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 1.11947));

    MassPoint mpisr145 (700, 630);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 1.09326));

    MassPoint mpisr146 (700, 640);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 1.12441));

    MassPoint mpisr147 (700, 650);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 1.13781));

    MassPoint mpisr148 (700, 660);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 1.12405));

    MassPoint mpisr149 (700, 670);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 1.08469));

    MassPoint mpisr150 (700, 680);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 1.12055));

    MassPoint mpisr151 (700, 690);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 1.08751));

    MassPoint mpisr152 (725, 645);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 1.11073));

    MassPoint mpisr153 (725, 655);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 1.1149));

    MassPoint mpisr154 (725, 665);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 1.10102));

    MassPoint mpisr155 (725, 675);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 1.14248));

    MassPoint mpisr156 (725, 685);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 1.11597));

    MassPoint mpisr157 (725, 695);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 1.10771));

    MassPoint mpisr158 (725, 705);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 1.13643));

    MassPoint mpisr159 (725, 715);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 1.15252));

    MassPoint mpisr160 (750, 670);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 1.10376));

    MassPoint mpisr161 (750, 680);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 1.10571));

    MassPoint mpisr162 (750, 690);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 1.11908));

    MassPoint mpisr163 (750, 700);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 1.12139));

    MassPoint mpisr164 (750, 710);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 1.12358));

    MassPoint mpisr165 (750, 720);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 1.0925));

    MassPoint mpisr166 (750, 730);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 1.10219));

    MassPoint mpisr167 (750, 740);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 1.08507));

    MassPoint mpisr168 (775, 695);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 1.09527));

    MassPoint mpisr169 (775, 705);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 1.11478));

    MassPoint mpisr170 (775, 715);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 1.1017));

    MassPoint mpisr171 (775, 725);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 1.12332));

    MassPoint mpisr172 (775, 735);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 1.10247));

    MassPoint mpisr173 (775, 745);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 1.16125));

    MassPoint mpisr174 (775, 755);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 1.09552));

    MassPoint mpisr175 (775, 765);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 1.09196));

    MassPoint mpisr176 (800, 720);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 1.10733));

    MassPoint mpisr177 (800, 730);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 1.10749));

    MassPoint mpisr178 (800, 740);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 1.10998));

    MassPoint mpisr179 (800, 750);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 1.10816));

    MassPoint mpisr180 (800, 760);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 1.11668));

    MassPoint mpisr181 (800, 770);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 1.11462));

    MassPoint mpisr182 (800, 780);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 1.10283));

    MassPoint mpisr183 (800, 790);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 1.08246));

  } else  if (SUSYProductionProcess=="TChiWW") {

    MassPoint mp0 (100, 1);
    StopCrossSection cs0 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp0 (cs0, 9420);
    StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

    MassPoint mp1 (100, 10);
    StopCrossSection cs1 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp1 (cs1, 12284);
    StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

    MassPoint mp2 (100, 20);
    StopCrossSection cs2 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp2 (cs2, 10380);
    StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

    MassPoint mp3 (100, 30);
    StopCrossSection cs3 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp3 (cs3, 9929);
    StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

    MassPoint mp4 (100, 40);
    StopCrossSection cs4 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp4 (cs4, 10266);
    StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

    MassPoint mp5 (100, 50);
    StopCrossSection cs5 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp5 (cs5, 8878);
    StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

    MassPoint mp6 (100, 60);
    StopCrossSection cs6 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp6 (cs6, 14979);
    StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

    MassPoint mp7 (100, 70);
    StopCrossSection cs7 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp7 (cs7, 8922);
    StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

    MassPoint mp8 (100, 80);
    StopCrossSection cs8 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp8 (cs8, 9563);
    StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

    MassPoint mp9 (100, 90);
    StopCrossSection cs9 (11.6119*0.10497000068426132, 0.518613*0.10497000068426132);
    MassPointParameters mpp9 (cs9, 9979);
    StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

    MassPoint mp10 (125, 1);
    StopCrossSection cs10 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp10 (cs10, 9697);
    StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

    MassPoint mp11 (125, 25);
    StopCrossSection cs11 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp11 (cs11, 9402);
    StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

    MassPoint mp12 (125, 35);
    StopCrossSection cs12 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp12 (cs12, 8741);
    StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

    MassPoint mp13 (125, 45);
    StopCrossSection cs13 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp13 (cs13, 9925);
    StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

    MassPoint mp14 (125, 55);
    StopCrossSection cs14 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp14 (cs14, 10919);
    StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

    MassPoint mp15 (125, 65);
    StopCrossSection cs15 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp15 (cs15, 11634);
    StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

    MassPoint mp16 (125, 75);
    StopCrossSection cs16 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp16 (cs16, 9531);
    StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

    MassPoint mp17 (125, 85);
    StopCrossSection cs17 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp17 (cs17, 10189);
    StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

    MassPoint mp18 (125, 95);
    StopCrossSection cs18 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp18 (cs18, 10010);
    StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

    MassPoint mp19 (125, 105);
    StopCrossSection cs19 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp19 (cs19, 6833);
    StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

    MassPoint mp20 (125, 115);
    StopCrossSection cs20 (5.09052*0.10497000068426132, 0.249469*0.10497000068426132);
    MassPointParameters mpp20 (cs20, 9177);
    StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

    MassPoint mp21 (150, 1);
    StopCrossSection cs21 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp21 (cs21, 8712);
    StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

    MassPoint mp22 (150, 25);
    StopCrossSection cs22 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp22 (cs22, 9998);
    StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

    MassPoint mp23 (150, 50);
    StopCrossSection cs23 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp23 (cs23, 9000);
    StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

    MassPoint mp24 (150, 60);
    StopCrossSection cs24 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp24 (cs24, 10011);
    StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

    MassPoint mp25 (150, 70);
    StopCrossSection cs25 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp25 (cs25, 10082);
    StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

    MassPoint mp26 (150, 80);
    StopCrossSection cs26 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp26 (cs26, 10965);
    StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

    MassPoint mp27 (150, 90);
    StopCrossSection cs27 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp27 (cs27, 9410);
    StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

    MassPoint mp28 (150, 100);
    StopCrossSection cs28 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp28 (cs28, 10688);
    StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

    MassPoint mp29 (150, 110);
    StopCrossSection cs29 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp29 (cs29, 9391);
    StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

    MassPoint mp30 (150, 120);
    StopCrossSection cs30 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp30 (cs30, 12049);
    StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

    MassPoint mp31 (150, 130);
    StopCrossSection cs31 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp31 (cs31, 12286);
    StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

    MassPoint mp32 (150, 140);
    StopCrossSection cs32 (2.61231*0.10497000068426132, 0.138156*0.10497000068426132);
    MassPointParameters mpp32 (cs32, 9822);
    StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

    MassPoint mp33 (175, 1);
    StopCrossSection cs33 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp33 (cs33, 11653);
    StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

    MassPoint mp34 (175, 25);
    StopCrossSection cs34 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp34 (cs34, 9207);
    StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

    MassPoint mp35 (175, 50);
    StopCrossSection cs35 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp35 (cs35, 13018);
    StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

    MassPoint mp36 (175, 75);
    StopCrossSection cs36 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp36 (cs36, 9655);
    StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

    MassPoint mp37 (175, 85);
    StopCrossSection cs37 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp37 (cs37, 9239);
    StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

    MassPoint mp38 (175, 95);
    StopCrossSection cs38 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp38 (cs38, 7660);
    StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

    MassPoint mp39 (175, 105);
    StopCrossSection cs39 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp39 (cs39, 10277);
    StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

    MassPoint mp40 (175, 115);
    StopCrossSection cs40 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp40 (cs40, 11310);
    StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

    MassPoint mp41 (175, 125);
    StopCrossSection cs41 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp41 (cs41, 9089);
    StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

    MassPoint mp42 (175, 135);
    StopCrossSection cs42 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp42 (cs42, 7766);
    StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

    MassPoint mp43 (175, 145);
    StopCrossSection cs43 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp43 (cs43, 9057);
    StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

    MassPoint mp44 (175, 155);
    StopCrossSection cs44 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp44 (cs44, 10360);
    StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

    MassPoint mp45 (175, 165);
    StopCrossSection cs45 (1.48242*0.10497000068426132, 0.0832672*0.10497000068426132);
    MassPointParameters mpp45 (cs45, 9230);
    StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

    MassPoint mp46 (200, 1);
    StopCrossSection cs46 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp46 (cs46, 9331);
    StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

    MassPoint mp47 (200, 25);
    StopCrossSection cs47 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp47 (cs47, 8336);
    StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

    MassPoint mp48 (200, 50);
    StopCrossSection cs48 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp48 (cs48, 6692);
    StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

    MassPoint mp49 (200, 75);
    StopCrossSection cs49 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp49 (cs49, 9266);
    StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

    MassPoint mp50 (200, 100);
    StopCrossSection cs50 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp50 (cs50, 9616);
    StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

    MassPoint mp51 (200, 110);
    StopCrossSection cs51 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp51 (cs51, 9147);
    StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

    MassPoint mp52 (200, 120);
    StopCrossSection cs52 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp52 (cs52, 7513);
    StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

    MassPoint mp53 (200, 130);
    StopCrossSection cs53 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp53 (cs53, 9170);
    StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

    MassPoint mp54 (200, 140);
    StopCrossSection cs54 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp54 (cs54, 8397);
    StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

    MassPoint mp55 (200, 150);
    StopCrossSection cs55 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp55 (cs55, 8211);
    StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

    MassPoint mp56 (200, 160);
    StopCrossSection cs56 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp56 (cs56, 11887);
    StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

    MassPoint mp57 (200, 170);
    StopCrossSection cs57 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp57 (cs57, 8654);
    StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

    MassPoint mp58 (200, 180);
    StopCrossSection cs58 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp58 (cs58, 10365);
    StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

    MassPoint mp59 (200, 190);
    StopCrossSection cs59 (0.902569*0.10497000068426132, 0.0537411*0.10497000068426132);
    MassPointParameters mpp59 (cs59, 7246);
    StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

    MassPoint mp60 (225, 1);
    StopCrossSection cs60 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp60 (cs60, 10607);
    StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

    MassPoint mp61 (225, 25);
    StopCrossSection cs61 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp61 (cs61, 10514);
    StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

    MassPoint mp62 (225, 50);
    StopCrossSection cs62 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp62 (cs62, 10030);
    StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

    MassPoint mp63 (225, 75);
    StopCrossSection cs63 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp63 (cs63, 7324);
    StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

    MassPoint mp64 (225, 100);
    StopCrossSection cs64 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp64 (cs64, 11695);
    StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

    MassPoint mp65 (225, 125);
    StopCrossSection cs65 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp65 (cs65, 8181);
    StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

    MassPoint mp66 (225, 135);
    StopCrossSection cs66 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp66 (cs66, 11536);
    StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

    MassPoint mp67 (225, 145);
    StopCrossSection cs67 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp67 (cs67, 11047);
    StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

    MassPoint mp68 (225, 155);
    StopCrossSection cs68 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp68 (cs68, 11502);
    StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

    MassPoint mp69 (225, 165);
    StopCrossSection cs69 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp69 (cs69, 8948);
    StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

    MassPoint mp70 (225, 175);
    StopCrossSection cs70 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp70 (cs70, 9969);
    StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

    MassPoint mp71 (225, 185);
    StopCrossSection cs71 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp71 (cs71, 8481);
    StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

    MassPoint mp72 (225, 195);
    StopCrossSection cs72 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp72 (cs72, 10373);
    StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

    MassPoint mp73 (225, 205);
    StopCrossSection cs73 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp73 (cs73, 10074);
    StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

    MassPoint mp74 (225, 215);
    StopCrossSection cs74 (0.579564*0.10497000068426132, 0.0360699*0.10497000068426132);
    MassPointParameters mpp74 (cs74, 12059);
    StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

    MassPoint mp75 (250, 1);
    StopCrossSection cs75 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp75 (cs75, 10086);
    StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

    MassPoint mp76 (250, 25);
    StopCrossSection cs76 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp76 (cs76, 9025);
    StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

    MassPoint mp77 (250, 50);
    StopCrossSection cs77 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp77 (cs77, 10152);
    StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

    MassPoint mp78 (250, 75);
    StopCrossSection cs78 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp78 (cs78, 10108);
    StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

    MassPoint mp79 (250, 100);
    StopCrossSection cs79 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp79 (cs79, 9166);
    StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

    MassPoint mp80 (250, 125);
    StopCrossSection cs80 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp80 (cs80, 8417);
    StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

    MassPoint mp81 (250, 150);
    StopCrossSection cs81 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp81 (cs81, 9471);
    StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

    MassPoint mp82 (250, 160);
    StopCrossSection cs82 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp82 (cs82, 10351);
    StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

    MassPoint mp83 (250, 170);
    StopCrossSection cs83 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp83 (cs83, 6554);
    StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

    MassPoint mp84 (250, 180);
    StopCrossSection cs84 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp84 (cs84, 10184);
    StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

    MassPoint mp85 (250, 190);
    StopCrossSection cs85 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp85 (cs85, 10063);
    StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

    MassPoint mp86 (250, 200);
    StopCrossSection cs86 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp86 (cs86, 10064);
    StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

    MassPoint mp87 (250, 210);
    StopCrossSection cs87 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp87 (cs87, 8928);
    StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

    MassPoint mp88 (250, 220);
    StopCrossSection cs88 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp88 (cs88, 10957);
    StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

    MassPoint mp89 (250, 230);
    StopCrossSection cs89 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp89 (cs89, 10765);
    StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

    MassPoint mp90 (250, 240);
    StopCrossSection cs90 (0.387534*0.10497000068426132, 0.0253131*0.10497000068426132);
    MassPointParameters mpp90 (cs90, 10244);
    StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

    MassPoint mp91 (275, 1);
    StopCrossSection cs91 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp91 (cs91, 9317);
    StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

    MassPoint mp92 (275, 25);
    StopCrossSection cs92 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp92 (cs92, 10812);
    StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

    MassPoint mp93 (275, 50);
    StopCrossSection cs93 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp93 (cs93, 8649);
    StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

    MassPoint mp94 (275, 75);
    StopCrossSection cs94 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp94 (cs94, 9468);
    StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

    MassPoint mp95 (275, 100);
    StopCrossSection cs95 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp95 (cs95, 13109);
    StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

    MassPoint mp96 (275, 125);
    StopCrossSection cs96 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp96 (cs96, 10211);
    StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

    MassPoint mp97 (275, 150);
    StopCrossSection cs97 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp97 (cs97, 10414);
    StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

    MassPoint mp98 (275, 175);
    StopCrossSection cs98 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp98 (cs98, 9190);
    StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

    MassPoint mp99 (275, 185);
    StopCrossSection cs99 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp99 (cs99, 10527);
    StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

    MassPoint mp100 (275, 195);
    StopCrossSection cs100 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp100 (cs100, 9795);
    StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

    MassPoint mp101 (275, 205);
    StopCrossSection cs101 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp101 (cs101, 10410);
    StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

    MassPoint mp102 (275, 215);
    StopCrossSection cs102 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp102 (cs102, 7931);
    StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

    MassPoint mp103 (275, 225);
    StopCrossSection cs103 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp103 (cs103, 9466);
    StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

    MassPoint mp104 (275, 235);
    StopCrossSection cs104 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp104 (cs104, 6385);
    StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

    MassPoint mp105 (275, 245);
    StopCrossSection cs105 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp105 (cs105, 9462);
    StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

    MassPoint mp106 (275, 255);
    StopCrossSection cs106 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp106 (cs106, 7884);
    StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

    MassPoint mp107 (275, 265);
    StopCrossSection cs107 (0.267786*0.10497000068426132, 0.0182886*0.10497000068426132);
    MassPointParameters mpp107 (cs107, 7044);
    StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

    MassPoint mp108 (300, 1);
    StopCrossSection cs108 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp108 (cs108, 12852);
    StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

    MassPoint mp109 (300, 25);
    StopCrossSection cs109 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp109 (cs109, 12218);
    StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

    MassPoint mp110 (300, 50);
    StopCrossSection cs110 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp110 (cs110, 8669);
    StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

    MassPoint mp111 (300, 75);
    StopCrossSection cs111 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp111 (cs111, 10154);
    StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

    MassPoint mp112 (300, 100);
    StopCrossSection cs112 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp112 (cs112, 11072);
    StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

    MassPoint mp113 (300, 125);
    StopCrossSection cs113 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp113 (cs113, 10010);
    StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

    MassPoint mp114 (300, 150);
    StopCrossSection cs114 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp114 (cs114, 13387);
    StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

    MassPoint mp115 (300, 175);
    StopCrossSection cs115 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp115 (cs115, 8317);
    StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

    MassPoint mp116 (300, 200);
    StopCrossSection cs116 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp116 (cs116, 12459);
    StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

    MassPoint mp117 (300, 210);
    StopCrossSection cs117 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp117 (cs117, 9943);
    StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

    MassPoint mp118 (300, 220);
    StopCrossSection cs118 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp118 (cs118, 12212);
    StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

    MassPoint mp119 (300, 230);
    StopCrossSection cs119 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp119 (cs119, 10049);
    StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

    MassPoint mp120 (300, 240);
    StopCrossSection cs120 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp120 (cs120, 9910);
    StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

    MassPoint mp121 (300, 250);
    StopCrossSection cs121 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp121 (cs121, 7641);
    StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

    MassPoint mp122 (300, 260);
    StopCrossSection cs122 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp122 (cs122, 10489);
    StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

    MassPoint mp123 (300, 270);
    StopCrossSection cs123 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp123 (cs123, 10869);
    StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

    MassPoint mp124 (300, 280);
    StopCrossSection cs124 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp124 (cs124, 8367);
    StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

    MassPoint mp125 (300, 290);
    StopCrossSection cs125 (0.190159*0.10497000068426132, 0.0134438*0.10497000068426132);
    MassPointParameters mpp125 (cs125, 9204);
    StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

    MassPoint mp126 (325, 1);
    StopCrossSection cs126 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp126 (cs126, 10764);
    StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

    MassPoint mp127 (325, 25);
    StopCrossSection cs127 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp127 (cs127, 8367);
    StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

    MassPoint mp128 (325, 50);
    StopCrossSection cs128 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp128 (cs128, 7790);
    StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

    MassPoint mp129 (325, 75);
    StopCrossSection cs129 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp129 (cs129, 9075);
    StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

    MassPoint mp130 (325, 100);
    StopCrossSection cs130 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp130 (cs130, 10646);
    StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

    MassPoint mp131 (325, 125);
    StopCrossSection cs131 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp131 (cs131, 8765);
    StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

    MassPoint mp132 (325, 150);
    StopCrossSection cs132 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp132 (cs132, 7582);
    StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

    MassPoint mp133 (325, 175);
    StopCrossSection cs133 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp133 (cs133, 9223);
    StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

    MassPoint mp134 (325, 200);
    StopCrossSection cs134 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp134 (cs134, 11785);
    StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

    MassPoint mp135 (325, 225);
    StopCrossSection cs135 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp135 (cs135, 9481);
    StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

    MassPoint mp136 (325, 235);
    StopCrossSection cs136 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp136 (cs136, 8498);
    StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

    MassPoint mp137 (325, 245);
    StopCrossSection cs137 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp137 (cs137, 10412);
    StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

    MassPoint mp138 (325, 255);
    StopCrossSection cs138 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp138 (cs138, 10513);
    StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

    MassPoint mp139 (325, 265);
    StopCrossSection cs139 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp139 (cs139, 12402);
    StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

    MassPoint mp140 (325, 275);
    StopCrossSection cs140 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp140 (cs140, 9491);
    StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

    MassPoint mp141 (325, 285);
    StopCrossSection cs141 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp141 (cs141, 13370);
    StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

    MassPoint mp142 (325, 295);
    StopCrossSection cs142 (0.138086*0.10497000068426132, 0.0101835*0.10497000068426132);
    MassPointParameters mpp142 (cs142, 10203);
    StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

    MassPoint mp143 (350, 1);
    StopCrossSection cs143 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp143 (cs143, 10587);
    StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

    MassPoint mp144 (350, 25);
    StopCrossSection cs144 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp144 (cs144, 8897);
    StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

    MassPoint mp145 (350, 50);
    StopCrossSection cs145 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp145 (cs145, 12844);
    StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

    MassPoint mp146 (350, 75);
    StopCrossSection cs146 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp146 (cs146, 10488);
    StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

    MassPoint mp147 (350, 100);
    StopCrossSection cs147 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp147 (cs147, 11414);
    StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

    MassPoint mp148 (350, 125);
    StopCrossSection cs148 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp148 (cs148, 10060);
    StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

    MassPoint mp149 (350, 150);
    StopCrossSection cs149 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp149 (cs149, 10025);
    StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

    MassPoint mp150 (350, 175);
    StopCrossSection cs150 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp150 (cs150, 8639);
    StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

    MassPoint mp151 (350, 200);
    StopCrossSection cs151 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp151 (cs151, 10403);
    StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

    MassPoint mp152 (350, 225);
    StopCrossSection cs152 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp152 (cs152, 9279);
    StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

    MassPoint mp153 (350, 250);
    StopCrossSection cs153 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp153 (cs153, 9705);
    StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

    MassPoint mp154 (350, 260);
    StopCrossSection cs154 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp154 (cs154, 8503);
    StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

    MassPoint mp155 (350, 270);
    StopCrossSection cs155 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp155 (cs155, 9666);
    StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

    MassPoint mp156 (350, 280);
    StopCrossSection cs156 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp156 (cs156, 11032);
    StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

    MassPoint mp157 (350, 290);
    StopCrossSection cs157 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp157 (cs157, 8563);
    StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

    MassPoint mp158 (350, 300);
    StopCrossSection cs158 (0.102199*0.10497000068426132, 0.00775261*0.10497000068426132);
    MassPointParameters mpp158 (cs158, 9625);
    StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

    MassPoint mp159 (375, 1);
    StopCrossSection cs159 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp159 (cs159, 10547);
    StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

    MassPoint mp160 (375, 25);
    StopCrossSection cs160 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp160 (cs160, 7438);
    StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

    MassPoint mp161 (375, 50);
    StopCrossSection cs161 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp161 (cs161, 9898);
    StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

    MassPoint mp162 (375, 75);
    StopCrossSection cs162 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp162 (cs162, 11821);
    StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

    MassPoint mp163 (375, 100);
    StopCrossSection cs163 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp163 (cs163, 10532);
    StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

    MassPoint mp164 (375, 125);
    StopCrossSection cs164 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp164 (cs164, 9018);
    StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

    MassPoint mp165 (375, 150);
    StopCrossSection cs165 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp165 (cs165, 11477);
    StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

    MassPoint mp166 (375, 175);
    StopCrossSection cs166 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp166 (cs166, 10054);
    StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

    MassPoint mp167 (375, 200);
    StopCrossSection cs167 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp167 (cs167, 9650);
    StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

    MassPoint mp168 (375, 225);
    StopCrossSection cs168 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp168 (cs168, 8401);
    StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

    MassPoint mp169 (375, 250);
    StopCrossSection cs169 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp169 (cs169, 8332);
    StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

    MassPoint mp170 (375, 275);
    StopCrossSection cs170 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp170 (cs170, 9073);
    StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

    MassPoint mp171 (375, 285);
    StopCrossSection cs171 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp171 (cs171, 10285);
    StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

    MassPoint mp172 (375, 295);
    StopCrossSection cs172 (0.0768342*0.10497000068426132, 0.00602606*0.10497000068426132);
    MassPointParameters mpp172 (cs172, 10839);
    StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

    MassPoint mp173 (400, 1);
    StopCrossSection cs173 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp173 (cs173, 10555);
    StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

    MassPoint mp174 (400, 25);
    StopCrossSection cs174 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp174 (cs174, 9777);
    StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

    MassPoint mp175 (400, 50);
    StopCrossSection cs175 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp175 (cs175, 12478);
    StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

    MassPoint mp176 (400, 75);
    StopCrossSection cs176 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp176 (cs176, 9871);
    StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

    MassPoint mp177 (400, 100);
    StopCrossSection cs177 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp177 (cs177, 11894);
    StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

    MassPoint mp178 (400, 125);
    StopCrossSection cs178 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp178 (cs178, 14373);
    StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

    MassPoint mp179 (400, 150);
    StopCrossSection cs179 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp179 (cs179, 9099);
    StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

    MassPoint mp180 (400, 175);
    StopCrossSection cs180 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp180 (cs180, 8715);
    StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

    MassPoint mp181 (400, 200);
    StopCrossSection cs181 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp181 (cs181, 9846);
    StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

    MassPoint mp182 (400, 225);
    StopCrossSection cs182 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp182 (cs182, 9546);
    StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

    MassPoint mp183 (400, 250);
    StopCrossSection cs183 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp183 (cs183, 10563);
    StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

    MassPoint mp184 (400, 275);
    StopCrossSection cs184 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp184 (cs184, 12184);
    StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

    MassPoint mp185 (400, 300);
    StopCrossSection cs185 (0.0586311*0.10497000068426132, 0.0047276*0.10497000068426132);
    MassPointParameters mpp185 (cs185, 9075);
    StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

    MassPoint mp186 (425, 1);
    StopCrossSection cs186 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp186 (cs186, 11567);
    StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

    MassPoint mp187 (425, 25);
    StopCrossSection cs187 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp187 (cs187, 12659);
    StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

    MassPoint mp188 (425, 50);
    StopCrossSection cs188 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp188 (cs188, 10392);
    StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

    MassPoint mp189 (425, 75);
    StopCrossSection cs189 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp189 (cs189, 10451);
    StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

    MassPoint mp190 (425, 100);
    StopCrossSection cs190 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp190 (cs190, 9876);
    StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

    MassPoint mp191 (425, 125);
    StopCrossSection cs191 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp191 (cs191, 11265);
    StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

    MassPoint mp192 (425, 150);
    StopCrossSection cs192 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp192 (cs192, 11434);
    StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

    MassPoint mp193 (425, 175);
    StopCrossSection cs193 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp193 (cs193, 10008);
    StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

    MassPoint mp194 (425, 200);
    StopCrossSection cs194 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp194 (cs194, 8629);
    StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

    MassPoint mp195 (425, 225);
    StopCrossSection cs195 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp195 (cs195, 8239);
    StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

    MassPoint mp196 (425, 250);
    StopCrossSection cs196 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp196 (cs196, 11700);
    StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

    MassPoint mp197 (425, 275);
    StopCrossSection cs197 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp197 (cs197, 11621);
    StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

    MassPoint mp198 (425, 300);
    StopCrossSection cs198 (0.0452189*0.10497000068426132, 0.00371547*0.10497000068426132);
    MassPointParameters mpp198 (cs198, 8663);
    StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

    MassPoint mp199 (450, 1);
    StopCrossSection cs199 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp199 (cs199, 10265);
    StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

    MassPoint mp200 (450, 25);
    StopCrossSection cs200 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp200 (cs200, 9405);
    StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

    MassPoint mp201 (450, 50);
    StopCrossSection cs201 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp201 (cs201, 12097);
    StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

    MassPoint mp202 (450, 75);
    StopCrossSection cs202 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp202 (cs202, 7898);
    StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

    MassPoint mp203 (450, 100);
    StopCrossSection cs203 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp203 (cs203, 10924);
    StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

    MassPoint mp204 (450, 125);
    StopCrossSection cs204 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp204 (cs204, 9094);
    StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

    MassPoint mp205 (450, 150);
    StopCrossSection cs205 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp205 (cs205, 10905);
    StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

    MassPoint mp206 (450, 175);
    StopCrossSection cs206 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp206 (cs206, 10757);
    StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

    MassPoint mp207 (450, 200);
    StopCrossSection cs207 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp207 (cs207, 7059);
    StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

    MassPoint mp208 (450, 225);
    StopCrossSection cs208 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp208 (cs208, 9364);
    StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

    MassPoint mp209 (450, 250);
    StopCrossSection cs209 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp209 (cs209, 11258);
    StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

    MassPoint mp210 (450, 275);
    StopCrossSection cs210 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp210 (cs210, 9498);
    StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

    MassPoint mp211 (450, 300);
    StopCrossSection cs211 (0.0353143*0.10497000068426132, 0.00297283*0.10497000068426132);
    MassPointParameters mpp211 (cs211, 10269);
    StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

    MassPoint mp212 (475, 1);
    StopCrossSection cs212 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp212 (cs212, 9549);
    StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

    MassPoint mp213 (475, 25);
    StopCrossSection cs213 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp213 (cs213, 10266);
    StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

    MassPoint mp214 (475, 50);
    StopCrossSection cs214 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp214 (cs214, 6495);
    StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

    MassPoint mp215 (475, 75);
    StopCrossSection cs215 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp215 (cs215, 8915);
    StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

    MassPoint mp216 (475, 100);
    StopCrossSection cs216 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp216 (cs216, 10572);
    StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

    MassPoint mp217 (475, 125);
    StopCrossSection cs217 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp217 (cs217, 7900);
    StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

    MassPoint mp218 (475, 150);
    StopCrossSection cs218 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp218 (cs218, 12024);
    StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

    MassPoint mp219 (475, 175);
    StopCrossSection cs219 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp219 (cs219, 10933);
    StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

    MassPoint mp220 (475, 200);
    StopCrossSection cs220 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp220 (cs220, 9902);
    StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

    MassPoint mp221 (475, 225);
    StopCrossSection cs221 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp221 (cs221, 9349);
    StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

    MassPoint mp222 (475, 250);
    StopCrossSection cs222 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp222 (cs222, 9808);
    StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

    MassPoint mp223 (475, 275);
    StopCrossSection cs223 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp223 (cs223, 9735);
    StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

    MassPoint mp224 (475, 300);
    StopCrossSection cs224 (0.0278342*0.10497000068426132, 0.00241224*0.10497000068426132);
    MassPointParameters mpp224 (cs224, 9611);
    StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

    MassPoint mp225 (500, 1);
    StopCrossSection cs225 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp225 (cs225, 9875);
    StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

    MassPoint mp226 (500, 25);
    StopCrossSection cs226 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp226 (cs226, 8958);
    StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

    MassPoint mp227 (500, 50);
    StopCrossSection cs227 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp227 (cs227, 11410);
    StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

    MassPoint mp228 (500, 75);
    StopCrossSection cs228 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp228 (cs228, 10445);
    StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

    MassPoint mp229 (500, 100);
    StopCrossSection cs229 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp229 (cs229, 10248);
    StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

    MassPoint mp230 (500, 125);
    StopCrossSection cs230 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp230 (cs230, 8968);
    StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

    MassPoint mp231 (500, 150);
    StopCrossSection cs231 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp231 (cs231, 8903);
    StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

    MassPoint mp232 (500, 175);
    StopCrossSection cs232 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp232 (cs232, 10196);
    StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

    MassPoint mp233 (500, 200);
    StopCrossSection cs233 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp233 (cs233, 8635);
    StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

    MassPoint mp234 (500, 225);
    StopCrossSection cs234 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp234 (cs234, 8907);
    StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

    MassPoint mp235 (500, 250);
    StopCrossSection cs235 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp235 (cs235, 9941);
    StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

    MassPoint mp236 (500, 275);
    StopCrossSection cs236 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp236 (cs236, 9858);
    StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

    MassPoint mp237 (500, 300);
    StopCrossSection cs237 (0.0221265*0.10497000068426132, 0.00194904*0.10497000068426132);
    MassPointParameters mpp237 (cs237, 8842);
    StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

    MassPoint mp238 (525, 1);
    StopCrossSection cs238 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp238 (cs238, 11236);
    StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

    MassPoint mp239 (525, 25);
    StopCrossSection cs239 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp239 (cs239, 10433);
    StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

    MassPoint mp240 (525, 50);
    StopCrossSection cs240 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp240 (cs240, 11863);
    StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

    MassPoint mp241 (525, 75);
    StopCrossSection cs241 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp241 (cs241, 9039);
    StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

    MassPoint mp242 (525, 100);
    StopCrossSection cs242 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp242 (cs242, 10290);
    StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

    MassPoint mp243 (525, 125);
    StopCrossSection cs243 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp243 (cs243, 11212);
    StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

    MassPoint mp244 (525, 150);
    StopCrossSection cs244 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp244 (cs244, 10205);
    StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

    MassPoint mp245 (525, 175);
    StopCrossSection cs245 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp245 (cs245, 12148);
    StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

    MassPoint mp246 (525, 200);
    StopCrossSection cs246 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp246 (cs246, 9229);
    StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

    MassPoint mp247 (525, 225);
    StopCrossSection cs247 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp247 (cs247, 8604);
    StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

    MassPoint mp248 (525, 250);
    StopCrossSection cs248 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp248 (cs248, 9957);
    StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

    MassPoint mp249 (525, 275);
    StopCrossSection cs249 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp249 (cs249, 10899);
    StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

    MassPoint mp250 (525, 300);
    StopCrossSection cs250 (0.0177394*0.10497000068426132, 0.0015992*0.10497000068426132);
    MassPointParameters mpp250 (cs250, 9637);
    StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

    MassPoint mp251 (550, 1);
    StopCrossSection cs251 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp251 (cs251, 11607);
    StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

    MassPoint mp252 (550, 25);
    StopCrossSection cs252 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp252 (cs252, 9464);
    StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

    MassPoint mp253 (550, 50);
    StopCrossSection cs253 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp253 (cs253, 12204);
    StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

    MassPoint mp254 (550, 75);
    StopCrossSection cs254 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp254 (cs254, 8792);
    StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

    MassPoint mp255 (550, 100);
    StopCrossSection cs255 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp255 (cs255, 10208);
    StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

    MassPoint mp256 (550, 125);
    StopCrossSection cs256 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp256 (cs256, 11061);
    StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

    MassPoint mp257 (550, 150);
    StopCrossSection cs257 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp257 (cs257, 9578);
    StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

    MassPoint mp258 (550, 175);
    StopCrossSection cs258 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp258 (cs258, 6966);
    StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

    MassPoint mp259 (550, 200);
    StopCrossSection cs259 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp259 (cs259, 10114);
    StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

    MassPoint mp260 (550, 225);
    StopCrossSection cs260 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp260 (cs260, 9440);
    StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

    MassPoint mp261 (550, 250);
    StopCrossSection cs261 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp261 (cs261, 9892);
    StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

    MassPoint mp262 (550, 275);
    StopCrossSection cs262 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp262 (cs262, 10043);
    StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

    MassPoint mp263 (550, 300);
    StopCrossSection cs263 (0.0143134*0.10497000068426132, 0.00132368*0.10497000068426132);
    MassPointParameters mpp263 (cs263, 10804);
    StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

    MassPoint mpisr0 (100, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 0.979482));

    MassPoint mpisr1 (100, 10);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 0.979291));

    MassPoint mpisr2 (100, 20);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 0.97965));

    MassPoint mpisr3 (100, 30);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 0.979744));

    MassPoint mpisr4 (100, 40);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 0.980261));

    MassPoint mpisr5 (100, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 0.979351));

    MassPoint mpisr6 (100, 60);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 0.979077));

    MassPoint mpisr7 (100, 70);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 0.979729));

    MassPoint mpisr8 (100, 80);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 0.978412));

    MassPoint mpisr9 (100, 90);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 0.978541));

    MassPoint mpisr10 (125, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 0.976435));

    MassPoint mpisr11 (125, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 0.976273));

    MassPoint mpisr12 (125, 35);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 0.976526));

    MassPoint mpisr13 (125, 45);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 0.976108));

    MassPoint mpisr14 (125, 55);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 0.975815));

    MassPoint mpisr15 (125, 65);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 0.97593));

    MassPoint mpisr16 (125, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 0.976072));

    MassPoint mpisr17 (125, 85);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 0.975889));

    MassPoint mpisr18 (125, 95);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 0.976967));

    MassPoint mpisr19 (125, 105);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 0.975936));

    MassPoint mpisr20 (125, 115);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 0.975948));

    MassPoint mpisr21 (150, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 0.973206));

    MassPoint mpisr22 (150, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 0.974135));

    MassPoint mpisr23 (150, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 0.974379));

    MassPoint mpisr24 (150, 60);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 0.974667));

    MassPoint mpisr25 (150, 70);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 0.974343));

    MassPoint mpisr26 (150, 80);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 0.973899));

    MassPoint mpisr27 (150, 90);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 0.973572));

    MassPoint mpisr28 (150, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 0.974067));

    MassPoint mpisr29 (150, 110);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 0.974052));

    MassPoint mpisr30 (150, 120);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 0.974155));

    MassPoint mpisr31 (150, 130);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 0.973967));

    MassPoint mpisr32 (150, 140);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 0.974501));

    MassPoint mpisr33 (175, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 0.971943));

    MassPoint mpisr34 (175, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 0.972609));

    MassPoint mpisr35 (175, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 0.971971));

    MassPoint mpisr36 (175, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 0.971853));

    MassPoint mpisr37 (175, 85);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 0.971327));

    MassPoint mpisr38 (175, 95);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 0.972113));

    MassPoint mpisr39 (175, 105);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 0.972523));

    MassPoint mpisr40 (175, 115);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 0.972248));

    MassPoint mpisr41 (175, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 0.972579));

    MassPoint mpisr42 (175, 135);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 0.971835));

    MassPoint mpisr43 (175, 145);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 0.97206));

    MassPoint mpisr44 (175, 155);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 0.972175));

    MassPoint mpisr45 (175, 165);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 0.971584));

    MassPoint mpisr46 (200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 0.971014));

    MassPoint mpisr47 (200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 0.97096));

    MassPoint mpisr48 (200, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 0.971121));

    MassPoint mpisr49 (200, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 0.969838));

    MassPoint mpisr50 (200, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 0.971022));

    MassPoint mpisr51 (200, 110);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 0.970827));

    MassPoint mpisr52 (200, 120);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 0.971459));

    MassPoint mpisr53 (200, 130);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 0.971269));

    MassPoint mpisr54 (200, 140);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 0.970483));

    MassPoint mpisr55 (200, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 0.97047));

    MassPoint mpisr56 (200, 160);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 0.970795));

    MassPoint mpisr57 (200, 170);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 0.971016));

    MassPoint mpisr58 (200, 180);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 0.971545));

    MassPoint mpisr59 (200, 190);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 0.970704));

    MassPoint mpisr60 (225, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 0.968939));

    MassPoint mpisr61 (225, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 0.969527));

    MassPoint mpisr62 (225, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 0.969568));

    MassPoint mpisr63 (225, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 0.969396));

    MassPoint mpisr64 (225, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 0.970176));

    MassPoint mpisr65 (225, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 0.969956));

    MassPoint mpisr66 (225, 135);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 0.969819));

    MassPoint mpisr67 (225, 145);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 0.969421));

    MassPoint mpisr68 (225, 155);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 0.970204));

    MassPoint mpisr69 (225, 165);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 0.97036));

    MassPoint mpisr70 (225, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 0.96964));

    MassPoint mpisr71 (225, 185);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 0.969453));

    MassPoint mpisr72 (225, 195);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 0.969386));

    MassPoint mpisr73 (225, 205);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 0.969487));

    MassPoint mpisr74 (225, 215);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 0.970003));

    MassPoint mpisr75 (250, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 0.969872));

    MassPoint mpisr76 (250, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 0.969372));

    MassPoint mpisr77 (250, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 0.968423));

    MassPoint mpisr78 (250, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 0.968808));

    MassPoint mpisr79 (250, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 0.969167));

    MassPoint mpisr80 (250, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 0.968867));

    MassPoint mpisr81 (250, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 0.969023));

    MassPoint mpisr82 (250, 160);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 0.969224));

    MassPoint mpisr83 (250, 170);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 0.968948));

    MassPoint mpisr84 (250, 180);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 0.968238));

    MassPoint mpisr85 (250, 190);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 0.969943));

    MassPoint mpisr86 (250, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 0.968186));

    MassPoint mpisr87 (250, 210);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 0.96977));

    MassPoint mpisr88 (250, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 0.969791));

    MassPoint mpisr89 (250, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 0.969142));

    MassPoint mpisr90 (250, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 0.969201));

    MassPoint mpisr91 (275, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 0.968658));

    MassPoint mpisr92 (275, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 0.96799));

    MassPoint mpisr93 (275, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 0.968896));

    MassPoint mpisr94 (275, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 0.969008));

    MassPoint mpisr95 (275, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 0.967596));

    MassPoint mpisr96 (275, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 0.967917));

    MassPoint mpisr97 (275, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 0.968283));

    MassPoint mpisr98 (275, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 0.968738));

    MassPoint mpisr99 (275, 185);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 0.96805));

    MassPoint mpisr100 (275, 195);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 0.967991));

    MassPoint mpisr101 (275, 205);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 0.968154));

    MassPoint mpisr102 (275, 215);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 0.968284));

    MassPoint mpisr103 (275, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 0.968724));

    MassPoint mpisr104 (275, 235);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 0.968599));

    MassPoint mpisr105 (275, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 0.967843));

    MassPoint mpisr106 (275, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 0.96842));

    MassPoint mpisr107 (275, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 0.96851));

    MassPoint mpisr108 (300, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 0.969111));

    MassPoint mpisr109 (300, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 0.967404));

    MassPoint mpisr110 (300, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 0.967042));

    MassPoint mpisr111 (300, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 0.967121));

    MassPoint mpisr112 (300, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 0.967979));

    MassPoint mpisr113 (300, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 0.967942));

    MassPoint mpisr114 (300, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 0.967327));

    MassPoint mpisr115 (300, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 0.968522));

    MassPoint mpisr116 (300, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 0.967));

    MassPoint mpisr117 (300, 210);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 0.967836));

    MassPoint mpisr118 (300, 220);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 0.967979));

    MassPoint mpisr119 (300, 230);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 0.968546));

    MassPoint mpisr120 (300, 240);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 0.967271));

    MassPoint mpisr121 (300, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 0.967847));

    MassPoint mpisr122 (300, 260);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 0.967843));

    MassPoint mpisr123 (300, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 0.967436));

    MassPoint mpisr124 (300, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 0.967264));

    MassPoint mpisr125 (300, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 0.967542));

    MassPoint mpisr126 (325, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 0.96702));

    MassPoint mpisr127 (325, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 0.966497));

    MassPoint mpisr128 (325, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 0.966322));

    MassPoint mpisr129 (325, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 0.966762));

    MassPoint mpisr130 (325, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 0.967293));

    MassPoint mpisr131 (325, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 0.967777));

    MassPoint mpisr132 (325, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 0.966301));

    MassPoint mpisr133 (325, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 0.966903));

    MassPoint mpisr134 (325, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 0.967456));

    MassPoint mpisr135 (325, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 0.967486));

    MassPoint mpisr136 (325, 235);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 0.96619));

    MassPoint mpisr137 (325, 245);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 0.96826));

    MassPoint mpisr138 (325, 255);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 0.968469));

    MassPoint mpisr139 (325, 265);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 0.967264));

    MassPoint mpisr140 (325, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 0.967516));

    MassPoint mpisr141 (325, 285);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 0.967752));

    MassPoint mpisr142 (325, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 0.968066));

    MassPoint mpisr143 (350, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 0.967278));

    MassPoint mpisr144 (350, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 0.968112));

    MassPoint mpisr145 (350, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 0.967276));

    MassPoint mpisr146 (350, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 0.966347));

    MassPoint mpisr147 (350, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 0.966795));

    MassPoint mpisr148 (350, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 0.966992));

    MassPoint mpisr149 (350, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 0.967448));

    MassPoint mpisr150 (350, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 0.967539));

    MassPoint mpisr151 (350, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 0.966998));

    MassPoint mpisr152 (350, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 0.966831));

    MassPoint mpisr153 (350, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 0.96626));

    MassPoint mpisr154 (350, 260);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 0.968198));

    MassPoint mpisr155 (350, 270);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 0.967095));

    MassPoint mpisr156 (350, 280);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 0.967125));

    MassPoint mpisr157 (350, 290);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 0.967827));

    MassPoint mpisr158 (350, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 0.967704));

    MassPoint mpisr159 (375, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 0.967122));

    MassPoint mpisr160 (375, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 0.966776));

    MassPoint mpisr161 (375, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 0.966362));

    MassPoint mpisr162 (375, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 0.966277));

    MassPoint mpisr163 (375, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 0.966935));

    MassPoint mpisr164 (375, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 0.966855));

    MassPoint mpisr165 (375, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 0.966982));

    MassPoint mpisr166 (375, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 0.9664));

    MassPoint mpisr167 (375, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 0.966695));

    MassPoint mpisr168 (375, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 0.96743));

    MassPoint mpisr169 (375, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 0.966327));

    MassPoint mpisr170 (375, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 0.965897));

    MassPoint mpisr171 (375, 285);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 0.967533));

    MassPoint mpisr172 (375, 295);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 0.966041));

    MassPoint mpisr173 (400, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 0.967399));

    MassPoint mpisr174 (400, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 0.966809));

    MassPoint mpisr175 (400, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 0.966637));

    MassPoint mpisr176 (400, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 0.966052));

    MassPoint mpisr177 (400, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 0.96504));

    MassPoint mpisr178 (400, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 0.966968));

    MassPoint mpisr179 (400, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 0.965626));

    MassPoint mpisr180 (400, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 0.966012));

    MassPoint mpisr181 (400, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 0.966552));

    MassPoint mpisr182 (400, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 0.966385));

    MassPoint mpisr183 (400, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 0.966515));

    MassPoint mpisr184 (400, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr184, 0.966889));

    MassPoint mpisr185 (400, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr185, 0.965676));

    MassPoint mpisr186 (425, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr186, 0.966367));

    MassPoint mpisr187 (425, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr187, 0.966412));

    MassPoint mpisr188 (425, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr188, 0.966341));

    MassPoint mpisr189 (425, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr189, 0.966677));

    MassPoint mpisr190 (425, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr190, 0.96607));

    MassPoint mpisr191 (425, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr191, 0.965802));

    MassPoint mpisr192 (425, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr192, 0.966712));

    MassPoint mpisr193 (425, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr193, 0.964715));

    MassPoint mpisr194 (425, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr194, 0.967058));

    MassPoint mpisr195 (425, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr195, 0.966075));

    MassPoint mpisr196 (425, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr196, 0.966253));

    MassPoint mpisr197 (425, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr197, 0.966277));

    MassPoint mpisr198 (425, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr198, 0.9648));

    MassPoint mpisr199 (450, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr199, 0.965784));

    MassPoint mpisr200 (450, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr200, 0.966933));

    MassPoint mpisr201 (450, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr201, 0.966534));

    MassPoint mpisr202 (450, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr202, 0.966135));

    MassPoint mpisr203 (450, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr203, 0.966227));

    MassPoint mpisr204 (450, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr204, 0.965503));

    MassPoint mpisr205 (450, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr205, 0.966469));

    MassPoint mpisr206 (450, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr206, 0.966023));

    MassPoint mpisr207 (450, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr207, 0.965813));

    MassPoint mpisr208 (450, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr208, 0.964627));

    MassPoint mpisr209 (450, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr209, 0.966138));

    MassPoint mpisr210 (450, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr210, 0.964627));

    MassPoint mpisr211 (450, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr211, 0.964879));

    MassPoint mpisr212 (475, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr212, 0.965821));

    MassPoint mpisr213 (475, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr213, 0.964233));

    MassPoint mpisr214 (475, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr214, 0.966288));

    MassPoint mpisr215 (475, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr215, 0.965493));

    MassPoint mpisr216 (475, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr216, 0.96523));

    MassPoint mpisr217 (475, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr217, 0.964697));

    MassPoint mpisr218 (475, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr218, 0.964927));

    MassPoint mpisr219 (475, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr219, 0.965595));

    MassPoint mpisr220 (475, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr220, 0.965318));

    MassPoint mpisr221 (475, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr221, 0.965805));

    MassPoint mpisr222 (475, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr222, 0.96491));

    MassPoint mpisr223 (475, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr223, 0.966017));

    MassPoint mpisr224 (475, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr224, 0.965332));

    MassPoint mpisr225 (500, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr225, 0.965965));

    MassPoint mpisr226 (500, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr226, 0.964581));

    MassPoint mpisr227 (500, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr227, 0.964036));

    MassPoint mpisr228 (500, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr228, 0.965495));

    MassPoint mpisr229 (500, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr229, 0.965576));

    MassPoint mpisr230 (500, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr230, 0.964829));

    MassPoint mpisr231 (500, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr231, 0.964861));

    MassPoint mpisr232 (500, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr232, 0.965929));

    MassPoint mpisr233 (500, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr233, 0.965567));

    MassPoint mpisr234 (500, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr234, 0.963732));

    MassPoint mpisr235 (500, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr235, 0.965099));

    MassPoint mpisr236 (500, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr236, 0.965363));

    MassPoint mpisr237 (500, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr237, 0.965687));

    MassPoint mpisr238 (525, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr238, 0.965232));

    MassPoint mpisr239 (525, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr239, 0.965107));

    MassPoint mpisr240 (525, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr240, 0.965822));

    MassPoint mpisr241 (525, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr241, 0.966384));

    MassPoint mpisr242 (525, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr242, 0.964546));

    MassPoint mpisr243 (525, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr243, 0.965466));

    MassPoint mpisr244 (525, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr244, 0.964069));

    MassPoint mpisr245 (525, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr245, 0.963896));

    MassPoint mpisr246 (525, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr246, 0.965183));

    MassPoint mpisr247 (525, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr247, 0.965134));

    MassPoint mpisr248 (525, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr248, 0.965367));

    MassPoint mpisr249 (525, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr249, 0.965549));

    MassPoint mpisr250 (525, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr250, 0.96399));

    MassPoint mpisr251 (550, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr251, 0.964865));

    MassPoint mpisr252 (550, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr252, 0.965371));

    MassPoint mpisr253 (550, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr253, 0.964985));

    MassPoint mpisr254 (550, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr254, 0.964544));

    MassPoint mpisr255 (550, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr255, 0.964609));

    MassPoint mpisr256 (550, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr256, 0.964677));

    MassPoint mpisr257 (550, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr257, 0.965843));

    MassPoint mpisr258 (550, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr258, 0.964733));

    MassPoint mpisr259 (550, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr259, 0.965905));

    MassPoint mpisr260 (550, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr260, 0.965081));

    MassPoint mpisr261 (550, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr261, 0.96486));

    MassPoint mpisr262 (550, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr262, 0.965365));

    MassPoint mpisr263 (550, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr263, 0.96469));

  } else  if (SUSYProductionProcess=="TChiStau") {

    MassPoint mp0 (100, 1);
    StopCrossSection cs0 (11.6119, 0.518613);
    MassPointParameters mpp0 (cs0, 44168);
    StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

    MassPoint mp1 (100, 25);
    StopCrossSection cs1 (11.6119, 0.518613);
    MassPointParameters mpp1 (cs1, 44704);
    StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

    MassPoint mp2 (100, 50);
    StopCrossSection cs2 (11.6119, 0.518613);
    MassPointParameters mpp2 (cs2, 49651);
    StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

    MassPoint mp3 (125, 1);
    StopCrossSection cs3 (5.09052, 0.249469);
    MassPointParameters mpp3 (cs3, 43146);
    StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

    MassPoint mp4 (125, 25);
    StopCrossSection cs4 (5.09052, 0.249469);
    MassPointParameters mpp4 (cs4, 41854);
    StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

    MassPoint mp5 (125, 50);
    StopCrossSection cs5 (5.09052, 0.249469);
    MassPointParameters mpp5 (cs5, 52042);
    StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

    MassPoint mp6 (125, 75);
    StopCrossSection cs6 (5.09052, 0.249469);
    MassPointParameters mpp6 (cs6, 48671);
    StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

    MassPoint mp7 (150, 1);
    StopCrossSection cs7 (2.61231, 0.138156);
    MassPointParameters mpp7 (cs7, 54802);
    StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

    MassPoint mp8 (150, 25);
    StopCrossSection cs8 (2.61231, 0.138156);
    MassPointParameters mpp8 (cs8, 52085);
    StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

    MassPoint mp9 (150, 50);
    StopCrossSection cs9 (2.61231, 0.138156);
    MassPointParameters mpp9 (cs9, 54363);
    StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

    MassPoint mp10 (150, 75);
    StopCrossSection cs10 (2.61231, 0.138156);
    MassPointParameters mpp10 (cs10, 45246);
    StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

    MassPoint mp11 (150, 100);
    StopCrossSection cs11 (2.61231, 0.138156);
    MassPointParameters mpp11 (cs11, 55831);
    StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

    MassPoint mp12 (175, 1);
    StopCrossSection cs12 (1.48242, 0.0832672);
    MassPointParameters mpp12 (cs12, 50049);
    StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

    MassPoint mp13 (175, 25);
    StopCrossSection cs13 (1.48242, 0.0832672);
    MassPointParameters mpp13 (cs13, 40788);
    StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

    MassPoint mp14 (175, 50);
    StopCrossSection cs14 (1.48242, 0.0832672);
    MassPointParameters mpp14 (cs14, 46623);
    StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

    MassPoint mp15 (175, 75);
    StopCrossSection cs15 (1.48242, 0.0832672);
    MassPointParameters mpp15 (cs15, 47727);
    StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

    MassPoint mp16 (175, 100);
    StopCrossSection cs16 (1.48242, 0.0832672);
    MassPointParameters mpp16 (cs16, 50555);
    StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

    MassPoint mp17 (175, 125);
    StopCrossSection cs17 (1.48242, 0.0832672);
    MassPointParameters mpp17 (cs17, 49645);
    StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

    MassPoint mp18 (200, 1);
    StopCrossSection cs18 (0.902569, 0.0537411);
    MassPointParameters mpp18 (cs18, 46045);
    StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

    MassPoint mp19 (200, 25);
    StopCrossSection cs19 (0.902569, 0.0537411);
    MassPointParameters mpp19 (cs19, 43907);
    StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

    MassPoint mp20 (200, 50);
    StopCrossSection cs20 (0.902569, 0.0537411);
    MassPointParameters mpp20 (cs20, 54253);
    StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

    MassPoint mp21 (200, 75);
    StopCrossSection cs21 (0.902569, 0.0537411);
    MassPointParameters mpp21 (cs21, 53546);
    StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

    MassPoint mp22 (200, 100);
    StopCrossSection cs22 (0.902569, 0.0537411);
    MassPointParameters mpp22 (cs22, 48295);
    StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

    MassPoint mp23 (200, 125);
    StopCrossSection cs23 (0.902569, 0.0537411);
    MassPointParameters mpp23 (cs23, 49652);
    StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

    MassPoint mp24 (200, 150);
    StopCrossSection cs24 (0.902569, 0.0537411);
    MassPointParameters mpp24 (cs24, 45246);
    StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

    MassPoint mp25 (225, 1);
    StopCrossSection cs25 (0.579564, 0.0360699);
    MassPointParameters mpp25 (cs25, 48295);
    StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

    MassPoint mp26 (225, 25);
    StopCrossSection cs26 (0.579564, 0.0360699);
    MassPointParameters mpp26 (cs26, 47199);
    StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

    MassPoint mp27 (225, 50);
    StopCrossSection cs27 (0.579564, 0.0360699);
    MassPointParameters mpp27 (cs27, 47222);
    StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

    MassPoint mp28 (225, 75);
    StopCrossSection cs28 (0.579564, 0.0360699);
    MassPointParameters mpp28 (cs28, 49390);
    StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

    MassPoint mp29 (225, 100);
    StopCrossSection cs29 (0.579564, 0.0360699);
    MassPointParameters mpp29 (cs29, 48253);
    StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

    MassPoint mp30 (225, 125);
    StopCrossSection cs30 (0.579564, 0.0360699);
    MassPointParameters mpp30 (cs30, 50548);
    StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

    MassPoint mp31 (225, 150);
    StopCrossSection cs31 (0.579564, 0.0360699);
    MassPointParameters mpp31 (cs31, 52865);
    StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

    MassPoint mp32 (225, 175);
    StopCrossSection cs32 (0.579564, 0.0360699);
    MassPointParameters mpp32 (cs32, 54711);
    StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

    MassPoint mp33 (250, 1);
    StopCrossSection cs33 (0.387534, 0.0253131);
    MassPointParameters mpp33 (cs33, 48833);
    StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

    MassPoint mp34 (250, 25);
    StopCrossSection cs34 (0.387534, 0.0253131);
    MassPointParameters mpp34 (cs34, 50294);
    StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

    MassPoint mp35 (250, 50);
    StopCrossSection cs35 (0.387534, 0.0253131);
    MassPointParameters mpp35 (cs35, 48679);
    StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

    MassPoint mp36 (250, 75);
    StopCrossSection cs36 (0.387534, 0.0253131);
    MassPointParameters mpp36 (cs36, 52960);
    StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

    MassPoint mp37 (250, 100);
    StopCrossSection cs37 (0.387534, 0.0253131);
    MassPointParameters mpp37 (cs37, 47810);
    StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

    MassPoint mp38 (250, 125);
    StopCrossSection cs38 (0.387534, 0.0253131);
    MassPointParameters mpp38 (cs38, 41500);
    StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

    MassPoint mp39 (250, 150);
    StopCrossSection cs39 (0.387534, 0.0253131);
    MassPointParameters mpp39 (cs39, 51971);
    StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

    MassPoint mp40 (250, 175);
    StopCrossSection cs40 (0.387534, 0.0253131);
    MassPointParameters mpp40 (cs40, 44111);
    StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

    MassPoint mp41 (250, 200);
    StopCrossSection cs41 (0.387534, 0.0253131);
    MassPointParameters mpp41 (cs41, 53180);
    StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

    MassPoint mp42 (275, 1);
    StopCrossSection cs42 (0.267786, 0.0182886);
    MassPointParameters mpp42 (cs42, 50229);
    StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

    MassPoint mp43 (275, 25);
    StopCrossSection cs43 (0.267786, 0.0182886);
    MassPointParameters mpp43 (cs43, 47630);
    StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

    MassPoint mp44 (275, 50);
    StopCrossSection cs44 (0.267786, 0.0182886);
    MassPointParameters mpp44 (cs44, 46932);
    StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

    MassPoint mp45 (275, 75);
    StopCrossSection cs45 (0.267786, 0.0182886);
    MassPointParameters mpp45 (cs45, 47800);
    StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

    MassPoint mp46 (275, 100);
    StopCrossSection cs46 (0.267786, 0.0182886);
    MassPointParameters mpp46 (cs46, 50566);
    StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

    MassPoint mp47 (275, 125);
    StopCrossSection cs47 (0.267786, 0.0182886);
    MassPointParameters mpp47 (cs47, 50265);
    StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

    MassPoint mp48 (275, 150);
    StopCrossSection cs48 (0.267786, 0.0182886);
    MassPointParameters mpp48 (cs48, 47796);
    StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

    MassPoint mp49 (275, 175);
    StopCrossSection cs49 (0.267786, 0.0182886);
    MassPointParameters mpp49 (cs49, 44594);
    StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

    MassPoint mp50 (275, 200);
    StopCrossSection cs50 (0.267786, 0.0182886);
    MassPointParameters mpp50 (cs50, 44726);
    StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

    MassPoint mp51 (275, 225);
    StopCrossSection cs51 (0.267786, 0.0182886);
    MassPointParameters mpp51 (cs51, 49753);
    StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

    MassPoint mp52 (300, 1);
    StopCrossSection cs52 (0.190159, 0.0134438);
    MassPointParameters mpp52 (cs52, 53138);
    StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

    MassPoint mp53 (300, 25);
    StopCrossSection cs53 (0.190159, 0.0134438);
    MassPointParameters mpp53 (cs53, 55234);
    StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

    MassPoint mp54 (300, 50);
    StopCrossSection cs54 (0.190159, 0.0134438);
    MassPointParameters mpp54 (cs54, 46042);
    StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

    MassPoint mp55 (300, 75);
    StopCrossSection cs55 (0.190159, 0.0134438);
    MassPointParameters mpp55 (cs55, 50565);
    StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

    MassPoint mp56 (300, 100);
    StopCrossSection cs56 (0.190159, 0.0134438);
    MassPointParameters mpp56 (cs56, 48047);
    StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

    MassPoint mp57 (300, 125);
    StopCrossSection cs57 (0.190159, 0.0134438);
    MassPointParameters mpp57 (cs57, 48436);
    StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

    MassPoint mp58 (300, 150);
    StopCrossSection cs58 (0.190159, 0.0134438);
    MassPointParameters mpp58 (cs58, 52201);
    StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

    MassPoint mp59 (300, 175);
    StopCrossSection cs59 (0.190159, 0.0134438);
    MassPointParameters mpp59 (cs59, 45733);
    StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

    MassPoint mp60 (300, 200);
    StopCrossSection cs60 (0.190159, 0.0134438);
    MassPointParameters mpp60 (cs60, 51335);
    StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

    MassPoint mp61 (300, 225);
    StopCrossSection cs61 (0.190159, 0.0134438);
    MassPointParameters mpp61 (cs61, 53053);
    StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

    MassPoint mp62 (300, 250);
    StopCrossSection cs62 (0.190159, 0.0134438);
    MassPointParameters mpp62 (cs62, 54844);
    StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

    MassPoint mp63 (325, 1);
    StopCrossSection cs63 (0.138086, 0.0101835);
    MassPointParameters mpp63 (cs63, 50587);
    StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

    MassPoint mp64 (325, 25);
    StopCrossSection cs64 (0.138086, 0.0101835);
    MassPointParameters mpp64 (cs64, 43017);
    StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

    MassPoint mp65 (325, 50);
    StopCrossSection cs65 (0.138086, 0.0101835);
    MassPointParameters mpp65 (cs65, 46860);
    StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

    MassPoint mp66 (325, 75);
    StopCrossSection cs66 (0.138086, 0.0101835);
    MassPointParameters mpp66 (cs66, 50615);
    StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

    MassPoint mp67 (325, 100);
    StopCrossSection cs67 (0.138086, 0.0101835);
    MassPointParameters mpp67 (cs67, 47278);
    StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

    MassPoint mp68 (325, 125);
    StopCrossSection cs68 (0.138086, 0.0101835);
    MassPointParameters mpp68 (cs68, 51544);
    StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

    MassPoint mp69 (325, 150);
    StopCrossSection cs69 (0.138086, 0.0101835);
    MassPointParameters mpp69 (cs69, 51708);
    StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

    MassPoint mp70 (325, 175);
    StopCrossSection cs70 (0.138086, 0.0101835);
    MassPointParameters mpp70 (cs70, 51115);
    StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

    MassPoint mp71 (325, 200);
    StopCrossSection cs71 (0.138086, 0.0101835);
    MassPointParameters mpp71 (cs71, 48942);
    StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

    MassPoint mp72 (325, 225);
    StopCrossSection cs72 (0.138086, 0.0101835);
    MassPointParameters mpp72 (cs72, 50650);
    StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

    MassPoint mp73 (325, 250);
    StopCrossSection cs73 (0.138086, 0.0101835);
    MassPointParameters mpp73 (cs73, 48847);
    StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

    MassPoint mp74 (325, 275);
    StopCrossSection cs74 (0.138086, 0.0101835);
    MassPointParameters mpp74 (cs74, 55897);
    StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

    MassPoint mp75 (350, 1);
    StopCrossSection cs75 (0.102199, 0.00775261);
    MassPointParameters mpp75 (cs75, 42248);
    StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

    MassPoint mp76 (350, 25);
    StopCrossSection cs76 (0.102199, 0.00775261);
    MassPointParameters mpp76 (cs76, 46438);
    StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

    MassPoint mp77 (350, 50);
    StopCrossSection cs77 (0.102199, 0.00775261);
    MassPointParameters mpp77 (cs77, 51544);
    StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

    MassPoint mp78 (350, 75);
    StopCrossSection cs78 (0.102199, 0.00775261);
    MassPointParameters mpp78 (cs78, 48376);
    StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

    MassPoint mp79 (350, 100);
    StopCrossSection cs79 (0.102199, 0.00775261);
    MassPointParameters mpp79 (cs79, 42591);
    StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

    MassPoint mp80 (350, 125);
    StopCrossSection cs80 (0.102199, 0.00775261);
    MassPointParameters mpp80 (cs80, 46458);
    StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

    MassPoint mp81 (350, 150);
    StopCrossSection cs81 (0.102199, 0.00775261);
    MassPointParameters mpp81 (cs81, 43171);
    StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

    MassPoint mp82 (350, 175);
    StopCrossSection cs82 (0.102199, 0.00775261);
    MassPointParameters mpp82 (cs82, 49760);
    StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

    MassPoint mp83 (350, 200);
    StopCrossSection cs83 (0.102199, 0.00775261);
    MassPointParameters mpp83 (cs83, 51197);
    StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

    MassPoint mp84 (350, 225);
    StopCrossSection cs84 (0.102199, 0.00775261);
    MassPointParameters mpp84 (cs84, 41253);
    StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

    MassPoint mp85 (350, 250);
    StopCrossSection cs85 (0.102199, 0.00775261);
    MassPointParameters mpp85 (cs85, 50708);
    StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

    MassPoint mp86 (350, 275);
    StopCrossSection cs86 (0.102199, 0.00775261);
    MassPointParameters mpp86 (cs86, 49405);
    StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

    MassPoint mp87 (350, 300);
    StopCrossSection cs87 (0.102199, 0.00775261);
    MassPointParameters mpp87 (cs87, 48136);
    StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

    MassPoint mp88 (375, 1);
    StopCrossSection cs88 (0.0768342, 0.00602606);
    MassPointParameters mpp88 (cs88, 51375);
    StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

    MassPoint mp89 (375, 25);
    StopCrossSection cs89 (0.0768342, 0.00602606);
    MassPointParameters mpp89 (cs89, 46669);
    StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

    MassPoint mp90 (375, 50);
    StopCrossSection cs90 (0.0768342, 0.00602606);
    MassPointParameters mpp90 (cs90, 47021);
    StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

    MassPoint mp91 (375, 75);
    StopCrossSection cs91 (0.0768342, 0.00602606);
    MassPointParameters mpp91 (cs91, 46646);
    StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

    MassPoint mp92 (375, 100);
    StopCrossSection cs92 (0.0768342, 0.00602606);
    MassPointParameters mpp92 (cs92, 47367);
    StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

    MassPoint mp93 (375, 125);
    StopCrossSection cs93 (0.0768342, 0.00602606);
    MassPointParameters mpp93 (cs93, 47243);
    StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

    MassPoint mp94 (375, 150);
    StopCrossSection cs94 (0.0768342, 0.00602606);
    MassPointParameters mpp94 (cs94, 51769);
    StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

    MassPoint mp95 (375, 175);
    StopCrossSection cs95 (0.0768342, 0.00602606);
    MassPointParameters mpp95 (cs95, 46694);
    StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

    MassPoint mp96 (375, 200);
    StopCrossSection cs96 (0.0768342, 0.00602606);
    MassPointParameters mpp96 (cs96, 46383);
    StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

    MassPoint mp97 (375, 225);
    StopCrossSection cs97 (0.0768342, 0.00602606);
    MassPointParameters mpp97 (cs97, 46487);
    StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

    MassPoint mp98 (375, 250);
    StopCrossSection cs98 (0.0768342, 0.00602606);
    MassPointParameters mpp98 (cs98, 44148);
    StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

    MassPoint mp99 (375, 275);
    StopCrossSection cs99 (0.0768342, 0.00602606);
    MassPointParameters mpp99 (cs99, 48977);
    StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

    MassPoint mp100 (375, 300);
    StopCrossSection cs100 (0.0768342, 0.00602606);
    MassPointParameters mpp100 (cs100, 55103);
    StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

    MassPoint mp101 (375, 325);
    StopCrossSection cs101 (0.0768342, 0.00602606);
    MassPointParameters mpp101 (cs101, 44515);
    StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

    MassPoint mp102 (400, 1);
    StopCrossSection cs102 (0.0586311, 0.0047276);
    MassPointParameters mpp102 (cs102, 49305);
    StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

    MassPoint mp103 (400, 25);
    StopCrossSection cs103 (0.0586311, 0.0047276);
    MassPointParameters mpp103 (cs103, 44793);
    StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

    MassPoint mp104 (400, 50);
    StopCrossSection cs104 (0.0586311, 0.0047276);
    MassPointParameters mpp104 (cs104, 48746);
    StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

    MassPoint mp105 (400, 75);
    StopCrossSection cs105 (0.0586311, 0.0047276);
    MassPointParameters mpp105 (cs105, 50661);
    StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

    MassPoint mp106 (400, 100);
    StopCrossSection cs106 (0.0586311, 0.0047276);
    MassPointParameters mpp106 (cs106, 49967);
    StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

    MassPoint mp107 (400, 125);
    StopCrossSection cs107 (0.0586311, 0.0047276);
    MassPointParameters mpp107 (cs107, 42918);
    StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

    MassPoint mp108 (400, 150);
    StopCrossSection cs108 (0.0586311, 0.0047276);
    MassPointParameters mpp108 (cs108, 45949);
    StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

    MassPoint mp109 (400, 175);
    StopCrossSection cs109 (0.0586311, 0.0047276);
    MassPointParameters mpp109 (cs109, 55332);
    StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

    MassPoint mp110 (400, 200);
    StopCrossSection cs110 (0.0586311, 0.0047276);
    MassPointParameters mpp110 (cs110, 50534);
    StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

    MassPoint mp111 (400, 225);
    StopCrossSection cs111 (0.0586311, 0.0047276);
    MassPointParameters mpp111 (cs111, 50573);
    StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

    MassPoint mp112 (400, 250);
    StopCrossSection cs112 (0.0586311, 0.0047276);
    MassPointParameters mpp112 (cs112, 50402);
    StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

    MassPoint mp113 (400, 275);
    StopCrossSection cs113 (0.0586311, 0.0047276);
    MassPointParameters mpp113 (cs113, 49443);
    StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

    MassPoint mp114 (400, 300);
    StopCrossSection cs114 (0.0586311, 0.0047276);
    MassPointParameters mpp114 (cs114, 47760);
    StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

    MassPoint mp115 (400, 325);
    StopCrossSection cs115 (0.0586311, 0.0047276);
    MassPointParameters mpp115 (cs115, 50429);
    StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

    MassPoint mp116 (400, 350);
    StopCrossSection cs116 (0.0586311, 0.0047276);
    MassPointParameters mpp116 (cs116, 51051);
    StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

    MassPoint mp117 (425, 1);
    StopCrossSection cs117 (0.0452189, 0.00371547);
    MassPointParameters mpp117 (cs117, 47877);
    StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

    MassPoint mp118 (425, 25);
    StopCrossSection cs118 (0.0452189, 0.00371547);
    MassPointParameters mpp118 (cs118, 46323);
    StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

    MassPoint mp119 (425, 50);
    StopCrossSection cs119 (0.0452189, 0.00371547);
    MassPointParameters mpp119 (cs119, 45129);
    StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

    MassPoint mp120 (425, 75);
    StopCrossSection cs120 (0.0452189, 0.00371547);
    MassPointParameters mpp120 (cs120, 48708);
    StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

    MassPoint mp121 (425, 100);
    StopCrossSection cs121 (0.0452189, 0.00371547);
    MassPointParameters mpp121 (cs121, 51087);
    StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

    MassPoint mp122 (425, 125);
    StopCrossSection cs122 (0.0452189, 0.00371547);
    MassPointParameters mpp122 (cs122, 48445);
    StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

    MassPoint mp123 (425, 150);
    StopCrossSection cs123 (0.0452189, 0.00371547);
    MassPointParameters mpp123 (cs123, 54166);
    StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

    MassPoint mp124 (425, 175);
    StopCrossSection cs124 (0.0452189, 0.00371547);
    MassPointParameters mpp124 (cs124, 50746);
    StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

    MassPoint mp125 (425, 200);
    StopCrossSection cs125 (0.0452189, 0.00371547);
    MassPointParameters mpp125 (cs125, 44941);
    StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

    MassPoint mp126 (425, 225);
    StopCrossSection cs126 (0.0452189, 0.00371547);
    MassPointParameters mpp126 (cs126, 46926);
    StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

    MassPoint mp127 (425, 250);
    StopCrossSection cs127 (0.0452189, 0.00371547);
    MassPointParameters mpp127 (cs127, 52769);
    StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

    MassPoint mp128 (425, 275);
    StopCrossSection cs128 (0.0452189, 0.00371547);
    MassPointParameters mpp128 (cs128, 45435);
    StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

    MassPoint mp129 (425, 300);
    StopCrossSection cs129 (0.0452189, 0.00371547);
    MassPointParameters mpp129 (cs129, 50006);
    StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

    MassPoint mp130 (425, 325);
    StopCrossSection cs130 (0.0452189, 0.00371547);
    MassPointParameters mpp130 (cs130, 46501);
    StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

    MassPoint mp131 (425, 350);
    StopCrossSection cs131 (0.0452189, 0.00371547);
    MassPointParameters mpp131 (cs131, 52746);
    StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

    MassPoint mp132 (425, 375);
    StopCrossSection cs132 (0.0452189, 0.00371547);
    MassPointParameters mpp132 (cs132, 47817);
    StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

    MassPoint mp133 (450, 1);
    StopCrossSection cs133 (0.0353143, 0.00297283);
    MassPointParameters mpp133 (cs133, 49478);
    StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

    MassPoint mp134 (450, 25);
    StopCrossSection cs134 (0.0353143, 0.00297283);
    MassPointParameters mpp134 (cs134, 49225);
    StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

    MassPoint mp135 (450, 50);
    StopCrossSection cs135 (0.0353143, 0.00297283);
    MassPointParameters mpp135 (cs135, 51626);
    StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

    MassPoint mp136 (450, 75);
    StopCrossSection cs136 (0.0353143, 0.00297283);
    MassPointParameters mpp136 (cs136, 49536);
    StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

    MassPoint mp137 (450, 100);
    StopCrossSection cs137 (0.0353143, 0.00297283);
    MassPointParameters mpp137 (cs137, 43933);
    StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

    MassPoint mp138 (450, 125);
    StopCrossSection cs138 (0.0353143, 0.00297283);
    MassPointParameters mpp138 (cs138, 45297);
    StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

    MassPoint mp139 (450, 150);
    StopCrossSection cs139 (0.0353143, 0.00297283);
    MassPointParameters mpp139 (cs139, 54727);
    StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

    MassPoint mp140 (450, 175);
    StopCrossSection cs140 (0.0353143, 0.00297283);
    MassPointParameters mpp140 (cs140, 47776);
    StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

    MassPoint mp141 (450, 200);
    StopCrossSection cs141 (0.0353143, 0.00297283);
    MassPointParameters mpp141 (cs141, 47759);
    StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

    MassPoint mp142 (450, 225);
    StopCrossSection cs142 (0.0353143, 0.00297283);
    MassPointParameters mpp142 (cs142, 52321);
    StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

    MassPoint mp143 (450, 250);
    StopCrossSection cs143 (0.0353143, 0.00297283);
    MassPointParameters mpp143 (cs143, 48819);
    StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

    MassPoint mp144 (450, 275);
    StopCrossSection cs144 (0.0353143, 0.00297283);
    MassPointParameters mpp144 (cs144, 49917);
    StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

    MassPoint mp145 (450, 300);
    StopCrossSection cs145 (0.0353143, 0.00297283);
    MassPointParameters mpp145 (cs145, 52125);
    StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

    MassPoint mp146 (450, 325);
    StopCrossSection cs146 (0.0353143, 0.00297283);
    MassPointParameters mpp146 (cs146, 47004);
    StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

    MassPoint mp147 (450, 350);
    StopCrossSection cs147 (0.0353143, 0.00297283);
    MassPointParameters mpp147 (cs147, 50835);
    StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

    MassPoint mp148 (450, 375);
    StopCrossSection cs148 (0.0353143, 0.00297283);
    MassPointParameters mpp148 (cs148, 48461);
    StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

    MassPoint mp149 (450, 400);
    StopCrossSection cs149 (0.0353143, 0.00297283);
    MassPointParameters mpp149 (cs149, 48588);
    StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

    MassPoint mp150 (475, 1);
    StopCrossSection cs150 (0.0278342, 0.00241224);
    MassPointParameters mpp150 (cs150, 43453);
    StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

    MassPoint mp151 (475, 25);
    StopCrossSection cs151 (0.0278342, 0.00241224);
    MassPointParameters mpp151 (cs151, 47063);
    StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

    MassPoint mp152 (475, 50);
    StopCrossSection cs152 (0.0278342, 0.00241224);
    MassPointParameters mpp152 (cs152, 50748);
    StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

    MassPoint mp153 (475, 75);
    StopCrossSection cs153 (0.0278342, 0.00241224);
    MassPointParameters mpp153 (cs153, 51314);
    StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

    MassPoint mp154 (475, 100);
    StopCrossSection cs154 (0.0278342, 0.00241224);
    MassPointParameters mpp154 (cs154, 52818);
    StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

    MassPoint mp155 (475, 125);
    StopCrossSection cs155 (0.0278342, 0.00241224);
    MassPointParameters mpp155 (cs155, 47298);
    StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

    MassPoint mp156 (475, 150);
    StopCrossSection cs156 (0.0278342, 0.00241224);
    MassPointParameters mpp156 (cs156, 47951);
    StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

    MassPoint mp157 (475, 175);
    StopCrossSection cs157 (0.0278342, 0.00241224);
    MassPointParameters mpp157 (cs157, 47872);
    StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

    MassPoint mp158 (475, 200);
    StopCrossSection cs158 (0.0278342, 0.00241224);
    MassPointParameters mpp158 (cs158, 52102);
    StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

    MassPoint mp159 (475, 225);
    StopCrossSection cs159 (0.0278342, 0.00241224);
    MassPointParameters mpp159 (cs159, 45317);
    StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

    MassPoint mp160 (475, 250);
    StopCrossSection cs160 (0.0278342, 0.00241224);
    MassPointParameters mpp160 (cs160, 47713);
    StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

    MassPoint mp161 (475, 275);
    StopCrossSection cs161 (0.0278342, 0.00241224);
    MassPointParameters mpp161 (cs161, 48251);
    StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

    MassPoint mp162 (475, 300);
    StopCrossSection cs162 (0.0278342, 0.00241224);
    MassPointParameters mpp162 (cs162, 50812);
    StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

    MassPoint mp163 (475, 325);
    StopCrossSection cs163 (0.0278342, 0.00241224);
    MassPointParameters mpp163 (cs163, 49146);
    StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

    MassPoint mp164 (475, 350);
    StopCrossSection cs164 (0.0278342, 0.00241224);
    MassPointParameters mpp164 (cs164, 47464);
    StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

    MassPoint mp165 (475, 375);
    StopCrossSection cs165 (0.0278342, 0.00241224);
    MassPointParameters mpp165 (cs165, 45184);
    StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

    MassPoint mp166 (475, 400);
    StopCrossSection cs166 (0.0278342, 0.00241224);
    MassPointParameters mpp166 (cs166, 40669);
    StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

    MassPoint mp167 (475, 425);
    StopCrossSection cs167 (0.0278342, 0.00241224);
    MassPointParameters mpp167 (cs167, 49499);
    StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

    MassPoint mp168 (500, 1);
    StopCrossSection cs168 (0.0221265, 0.00194904);
    MassPointParameters mpp168 (cs168, 54644);
    StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

    MassPoint mp169 (500, 25);
    StopCrossSection cs169 (0.0221265, 0.00194904);
    MassPointParameters mpp169 (cs169, 50474);
    StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

    MassPoint mp170 (500, 50);
    StopCrossSection cs170 (0.0221265, 0.00194904);
    MassPointParameters mpp170 (cs170, 46442);
    StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

    MassPoint mp171 (500, 75);
    StopCrossSection cs171 (0.0221265, 0.00194904);
    MassPointParameters mpp171 (cs171, 41924);
    StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

    MassPoint mp172 (500, 100);
    StopCrossSection cs172 (0.0221265, 0.00194904);
    MassPointParameters mpp172 (cs172, 46435);
    StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

    MassPoint mp173 (500, 125);
    StopCrossSection cs173 (0.0221265, 0.00194904);
    MassPointParameters mpp173 (cs173, 47343);
    StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

    MassPoint mp174 (500, 150);
    StopCrossSection cs174 (0.0221265, 0.00194904);
    MassPointParameters mpp174 (cs174, 48051);
    StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

    MassPoint mp175 (500, 175);
    StopCrossSection cs175 (0.0221265, 0.00194904);
    MassPointParameters mpp175 (cs175, 51888);
    StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

    MassPoint mp176 (500, 200);
    StopCrossSection cs176 (0.0221265, 0.00194904);
    MassPointParameters mpp176 (cs176, 42944);
    StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

    MassPoint mp177 (500, 225);
    StopCrossSection cs177 (0.0221265, 0.00194904);
    MassPointParameters mpp177 (cs177, 47612);
    StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

    MassPoint mp178 (500, 250);
    StopCrossSection cs178 (0.0221265, 0.00194904);
    MassPointParameters mpp178 (cs178, 48769);
    StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

    MassPoint mp179 (500, 275);
    StopCrossSection cs179 (0.0221265, 0.00194904);
    MassPointParameters mpp179 (cs179, 47555);
    StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

    MassPoint mp180 (500, 300);
    StopCrossSection cs180 (0.0221265, 0.00194904);
    MassPointParameters mpp180 (cs180, 50845);
    StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

    MassPoint mp181 (500, 325);
    StopCrossSection cs181 (0.0221265, 0.00194904);
    MassPointParameters mpp181 (cs181, 45177);
    StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

    MassPoint mp182 (500, 350);
    StopCrossSection cs182 (0.0221265, 0.00194904);
    MassPointParameters mpp182 (cs182, 51275);
    StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

    MassPoint mp183 (500, 375);
    StopCrossSection cs183 (0.0221265, 0.00194904);
    MassPointParameters mpp183 (cs183, 48370);
    StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

    MassPoint mp184 (500, 400);
    StopCrossSection cs184 (0.0221265, 0.00194904);
    MassPointParameters mpp184 (cs184, 45208);
    StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

    MassPoint mp185 (500, 425);
    StopCrossSection cs185 (0.0221265, 0.00194904);
    MassPointParameters mpp185 (cs185, 50458);
    StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

    MassPoint mp186 (500, 450);
    StopCrossSection cs186 (0.0221265, 0.00194904);
    MassPointParameters mpp186 (cs186, 45081);
    StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

    MassPoint mp187 (525, 1);
    StopCrossSection cs187 (0.0177394, 0.0015992);
    MassPointParameters mpp187 (cs187, 44598);
    StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

    MassPoint mp188 (525, 25);
    StopCrossSection cs188 (0.0177394, 0.0015992);
    MassPointParameters mpp188 (cs188, 47186);
    StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

    MassPoint mp189 (525, 50);
    StopCrossSection cs189 (0.0177394, 0.0015992);
    MassPointParameters mpp189 (cs189, 52792);
    StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

    MassPoint mp190 (525, 75);
    StopCrossSection cs190 (0.0177394, 0.0015992);
    MassPointParameters mpp190 (cs190, 53837);
    StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

    MassPoint mp191 (525, 100);
    StopCrossSection cs191 (0.0177394, 0.0015992);
    MassPointParameters mpp191 (cs191, 52956);
    StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

    MassPoint mp192 (525, 125);
    StopCrossSection cs192 (0.0177394, 0.0015992);
    MassPointParameters mpp192 (cs192, 53612);
    StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

    MassPoint mp193 (525, 150);
    StopCrossSection cs193 (0.0177394, 0.0015992);
    MassPointParameters mpp193 (cs193, 47889);
    StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

    MassPoint mp194 (525, 175);
    StopCrossSection cs194 (0.0177394, 0.0015992);
    MassPointParameters mpp194 (cs194, 50426);
    StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

    MassPoint mp195 (525, 200);
    StopCrossSection cs195 (0.0177394, 0.0015992);
    MassPointParameters mpp195 (cs195, 48427);
    StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

    MassPoint mp196 (525, 225);
    StopCrossSection cs196 (0.0177394, 0.0015992);
    MassPointParameters mpp196 (cs196, 48810);
    StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

    MassPoint mp197 (525, 250);
    StopCrossSection cs197 (0.0177394, 0.0015992);
    MassPointParameters mpp197 (cs197, 47635);
    StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

    MassPoint mp198 (525, 275);
    StopCrossSection cs198 (0.0177394, 0.0015992);
    MassPointParameters mpp198 (cs198, 49730);
    StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

    MassPoint mp199 (525, 300);
    StopCrossSection cs199 (0.0177394, 0.0015992);
    MassPointParameters mpp199 (cs199, 50812);
    StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

    MassPoint mp200 (525, 325);
    StopCrossSection cs200 (0.0177394, 0.0015992);
    MassPointParameters mpp200 (cs200, 47594);
    StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

    MassPoint mp201 (525, 350);
    StopCrossSection cs201 (0.0177394, 0.0015992);
    MassPointParameters mpp201 (cs201, 48992);
    StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

    MassPoint mp202 (525, 375);
    StopCrossSection cs202 (0.0177394, 0.0015992);
    MassPointParameters mpp202 (cs202, 52826);
    StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

    MassPoint mp203 (525, 400);
    StopCrossSection cs203 (0.0177394, 0.0015992);
    MassPointParameters mpp203 (cs203, 46197);
    StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

    MassPoint mp204 (525, 425);
    StopCrossSection cs204 (0.0177394, 0.0015992);
    MassPointParameters mpp204 (cs204, 50043);
    StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

    MassPoint mp205 (525, 450);
    StopCrossSection cs205 (0.0177394, 0.0015992);
    MassPointParameters mpp205 (cs205, 52260);
    StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

    MassPoint mp206 (525, 475);
    StopCrossSection cs206 (0.0177394, 0.0015992);
    MassPointParameters mpp206 (cs206, 48969);
    StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

    MassPoint mp207 (550, 1);
    StopCrossSection cs207 (0.0143134, 0.00132368);
    MassPointParameters mpp207 (cs207, 49488);
    StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

    MassPoint mp208 (550, 25);
    StopCrossSection cs208 (0.0143134, 0.00132368);
    MassPointParameters mpp208 (cs208, 49327);
    StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

    MassPoint mp209 (550, 50);
    StopCrossSection cs209 (0.0143134, 0.00132368);
    MassPointParameters mpp209 (cs209, 51714);
    StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

    MassPoint mp210 (550, 75);
    StopCrossSection cs210 (0.0143134, 0.00132368);
    MassPointParameters mpp210 (cs210, 53344);
    StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

    MassPoint mp211 (550, 100);
    StopCrossSection cs211 (0.0143134, 0.00132368);
    MassPointParameters mpp211 (cs211, 48333);
    StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

    MassPoint mp212 (550, 125);
    StopCrossSection cs212 (0.0143134, 0.00132368);
    MassPointParameters mpp212 (cs212, 53841);
    StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

    MassPoint mp213 (550, 150);
    StopCrossSection cs213 (0.0143134, 0.00132368);
    MassPointParameters mpp213 (cs213, 50576);
    StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

    MassPoint mp214 (550, 175);
    StopCrossSection cs214 (0.0143134, 0.00132368);
    MassPointParameters mpp214 (cs214, 47784);
    StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

    MassPoint mp215 (550, 200);
    StopCrossSection cs215 (0.0143134, 0.00132368);
    MassPointParameters mpp215 (cs215, 54623);
    StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

    MassPoint mp216 (550, 225);
    StopCrossSection cs216 (0.0143134, 0.00132368);
    MassPointParameters mpp216 (cs216, 49907);
    StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

    MassPoint mp217 (550, 250);
    StopCrossSection cs217 (0.0143134, 0.00132368);
    MassPointParameters mpp217 (cs217, 53260);
    StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

    MassPoint mp218 (550, 275);
    StopCrossSection cs218 (0.0143134, 0.00132368);
    MassPointParameters mpp218 (cs218, 52684);
    StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

    MassPoint mp219 (550, 300);
    StopCrossSection cs219 (0.0143134, 0.00132368);
    MassPointParameters mpp219 (cs219, 50095);
    StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

    MassPoint mp220 (550, 325);
    StopCrossSection cs220 (0.0143134, 0.00132368);
    MassPointParameters mpp220 (cs220, 46097);
    StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

    MassPoint mp221 (550, 350);
    StopCrossSection cs221 (0.0143134, 0.00132368);
    MassPointParameters mpp221 (cs221, 53728);
    StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

    MassPoint mp222 (550, 375);
    StopCrossSection cs222 (0.0143134, 0.00132368);
    MassPointParameters mpp222 (cs222, 51250);
    StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

    MassPoint mp223 (550, 400);
    StopCrossSection cs223 (0.0143134, 0.00132368);
    MassPointParameters mpp223 (cs223, 53630);
    StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

    MassPoint mp224 (550, 425);
    StopCrossSection cs224 (0.0143134, 0.00132368);
    MassPointParameters mpp224 (cs224, 49600);
    StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

    MassPoint mp225 (550, 450);
    StopCrossSection cs225 (0.0143134, 0.00132368);
    MassPointParameters mpp225 (cs225, 52892);
    StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

    MassPoint mp226 (550, 475);
    StopCrossSection cs226 (0.0143134, 0.00132368);
    MassPointParameters mpp226 (cs226, 52504);
    StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

    MassPoint mp227 (550, 500);
    StopCrossSection cs227 (0.0143134, 0.00132368);
    MassPointParameters mpp227 (cs227, 50430);
    StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

    MassPoint mp228 (575, 1);
    StopCrossSection cs228 (0.0116266, 0.00109669);
    MassPointParameters mpp228 (cs228, 47849);
    StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

    MassPoint mp229 (575, 25);
    StopCrossSection cs229 (0.0116266, 0.00109669);
    MassPointParameters mpp229 (cs229, 52509);
    StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

    MassPoint mp230 (575, 50);
    StopCrossSection cs230 (0.0116266, 0.00109669);
    MassPointParameters mpp230 (cs230, 51128);
    StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

    MassPoint mp231 (575, 75);
    StopCrossSection cs231 (0.0116266, 0.00109669);
    MassPointParameters mpp231 (cs231, 49426);
    StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

    MassPoint mp232 (575, 100);
    StopCrossSection cs232 (0.0116266, 0.00109669);
    MassPointParameters mpp232 (cs232, 52263);
    StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

    MassPoint mp233 (575, 125);
    StopCrossSection cs233 (0.0116266, 0.00109669);
    MassPointParameters mpp233 (cs233, 46747);
    StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

    MassPoint mp234 (575, 150);
    StopCrossSection cs234 (0.0116266, 0.00109669);
    MassPointParameters mpp234 (cs234, 48846);
    StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

    MassPoint mp235 (575, 175);
    StopCrossSection cs235 (0.0116266, 0.00109669);
    MassPointParameters mpp235 (cs235, 50822);
    StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

    MassPoint mp236 (575, 200);
    StopCrossSection cs236 (0.0116266, 0.00109669);
    MassPointParameters mpp236 (cs236, 53514);
    StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

    MassPoint mp237 (575, 225);
    StopCrossSection cs237 (0.0116266, 0.00109669);
    MassPointParameters mpp237 (cs237, 50054);
    StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

    MassPoint mp238 (575, 250);
    StopCrossSection cs238 (0.0116266, 0.00109669);
    MassPointParameters mpp238 (cs238, 52018);
    StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

    MassPoint mp239 (575, 275);
    StopCrossSection cs239 (0.0116266, 0.00109669);
    MassPointParameters mpp239 (cs239, 51749);
    StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

    MassPoint mp240 (575, 300);
    StopCrossSection cs240 (0.0116266, 0.00109669);
    MassPointParameters mpp240 (cs240, 49912);
    StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

    MassPoint mp241 (575, 325);
    StopCrossSection cs241 (0.0116266, 0.00109669);
    MassPointParameters mpp241 (cs241, 46171);
    StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

    MassPoint mp242 (575, 350);
    StopCrossSection cs242 (0.0116266, 0.00109669);
    MassPointParameters mpp242 (cs242, 51370);
    StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

    MassPoint mp243 (575, 375);
    StopCrossSection cs243 (0.0116266, 0.00109669);
    MassPointParameters mpp243 (cs243, 55962);
    StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

    MassPoint mp244 (575, 400);
    StopCrossSection cs244 (0.0116266, 0.00109669);
    MassPointParameters mpp244 (cs244, 51237);
    StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

    MassPoint mp245 (575, 425);
    StopCrossSection cs245 (0.0116266, 0.00109669);
    MassPointParameters mpp245 (cs245, 52559);
    StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

    MassPoint mp246 (575, 450);
    StopCrossSection cs246 (0.0116266, 0.00109669);
    MassPointParameters mpp246 (cs246, 49699);
    StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

    MassPoint mp247 (575, 475);
    StopCrossSection cs247 (0.0116266, 0.00109669);
    MassPointParameters mpp247 (cs247, 46908);
    StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

    MassPoint mp248 (575, 500);
    StopCrossSection cs248 (0.0116266, 0.00109669);
    MassPointParameters mpp248 (cs248, 44929);
    StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

    MassPoint mp249 (575, 525);
    StopCrossSection cs249 (0.0116266, 0.00109669);
    MassPointParameters mpp249 (cs249, 50646);
    StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

    MassPoint mp250 (600, 1);
    StopCrossSection cs250 (0.00949913, 0.000912324);
    MassPointParameters mpp250 (cs250, 53613);
    StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

    MassPoint mp251 (600, 25);
    StopCrossSection cs251 (0.00949913, 0.000912324);
    MassPointParameters mpp251 (cs251, 51880);
    StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

    MassPoint mp252 (600, 50);
    StopCrossSection cs252 (0.00949913, 0.000912324);
    MassPointParameters mpp252 (cs252, 48985);
    StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

    MassPoint mp253 (600, 75);
    StopCrossSection cs253 (0.00949913, 0.000912324);
    MassPointParameters mpp253 (cs253, 52556);
    StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

    MassPoint mp254 (600, 100);
    StopCrossSection cs254 (0.00949913, 0.000912324);
    MassPointParameters mpp254 (cs254, 47154);
    StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

    MassPoint mp255 (600, 125);
    StopCrossSection cs255 (0.00949913, 0.000912324);
    MassPointParameters mpp255 (cs255, 51322);
    StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

    MassPoint mp256 (600, 150);
    StopCrossSection cs256 (0.00949913, 0.000912324);
    MassPointParameters mpp256 (cs256, 49298);
    StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

    MassPoint mp257 (600, 175);
    StopCrossSection cs257 (0.00949913, 0.000912324);
    MassPointParameters mpp257 (cs257, 52058);
    StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

    MassPoint mp258 (600, 200);
    StopCrossSection cs258 (0.00949913, 0.000912324);
    MassPointParameters mpp258 (cs258, 51161);
    StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

    MassPoint mp259 (600, 225);
    StopCrossSection cs259 (0.00949913, 0.000912324);
    MassPointParameters mpp259 (cs259, 53003);
    StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

    MassPoint mp260 (600, 250);
    StopCrossSection cs260 (0.00949913, 0.000912324);
    MassPointParameters mpp260 (cs260, 50757);
    StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

    MassPoint mp261 (600, 275);
    StopCrossSection cs261 (0.00949913, 0.000912324);
    MassPointParameters mpp261 (cs261, 41548);
    StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

    MassPoint mp262 (600, 300);
    StopCrossSection cs262 (0.00949913, 0.000912324);
    MassPointParameters mpp262 (cs262, 46591);
    StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

    MassPoint mp263 (600, 325);
    StopCrossSection cs263 (0.00949913, 0.000912324);
    MassPointParameters mpp263 (cs263, 47768);
    StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

    MassPoint mp264 (600, 350);
    StopCrossSection cs264 (0.00949913, 0.000912324);
    MassPointParameters mpp264 (cs264, 50587);
    StopNeutralinoMap.insert(std::make_pair(mp264, mpp264));

    MassPoint mp265 (600, 375);
    StopCrossSection cs265 (0.00949913, 0.000912324);
    MassPointParameters mpp265 (cs265, 51760);
    StopNeutralinoMap.insert(std::make_pair(mp265, mpp265));

    MassPoint mp266 (600, 400);
    StopCrossSection cs266 (0.00949913, 0.000912324);
    MassPointParameters mpp266 (cs266, 48210);
    StopNeutralinoMap.insert(std::make_pair(mp266, mpp266));

    MassPoint mp267 (600, 425);
    StopCrossSection cs267 (0.00949913, 0.000912324);
    MassPointParameters mpp267 (cs267, 49047);
    StopNeutralinoMap.insert(std::make_pair(mp267, mpp267));

    MassPoint mp268 (600, 450);
    StopCrossSection cs268 (0.00949913, 0.000912324);
    MassPointParameters mpp268 (cs268, 50952);
    StopNeutralinoMap.insert(std::make_pair(mp268, mpp268));

    MassPoint mp269 (600, 475);
    StopCrossSection cs269 (0.00949913, 0.000912324);
    MassPointParameters mpp269 (cs269, 50018);
    StopNeutralinoMap.insert(std::make_pair(mp269, mpp269));

    MassPoint mp270 (600, 500);
    StopCrossSection cs270 (0.00949913, 0.000912324);
    MassPointParameters mpp270 (cs270, 45452);
    StopNeutralinoMap.insert(std::make_pair(mp270, mpp270));

    MassPoint mp271 (600, 525);
    StopCrossSection cs271 (0.00949913, 0.000912324);
    MassPointParameters mpp271 (cs271, 50522);
    StopNeutralinoMap.insert(std::make_pair(mp271, mpp271));

    MassPoint mp272 (600, 550);
    StopCrossSection cs272 (0.00949913, 0.000912324);
    MassPointParameters mpp272 (cs272, 47966);
    StopNeutralinoMap.insert(std::make_pair(mp272, mpp272));

    MassPoint mp273 (625, 1);
    StopCrossSection cs273 (0.00780081, 0.000768988);
    MassPointParameters mpp273 (cs273, 50628);
    StopNeutralinoMap.insert(std::make_pair(mp273, mpp273));

    MassPoint mp274 (625, 25);
    StopCrossSection cs274 (0.00780081, 0.000768988);
    MassPointParameters mpp274 (cs274, 51150);
    StopNeutralinoMap.insert(std::make_pair(mp274, mpp274));

    MassPoint mp275 (625, 50);
    StopCrossSection cs275 (0.00780081, 0.000768988);
    MassPointParameters mpp275 (cs275, 50141);
    StopNeutralinoMap.insert(std::make_pair(mp275, mpp275));

    MassPoint mp276 (625, 75);
    StopCrossSection cs276 (0.00780081, 0.000768988);
    MassPointParameters mpp276 (cs276, 40673);
    StopNeutralinoMap.insert(std::make_pair(mp276, mpp276));

    MassPoint mp277 (625, 100);
    StopCrossSection cs277 (0.00780081, 0.000768988);
    MassPointParameters mpp277 (cs277, 44829);
    StopNeutralinoMap.insert(std::make_pair(mp277, mpp277));

    MassPoint mp278 (625, 125);
    StopCrossSection cs278 (0.00780081, 0.000768988);
    MassPointParameters mpp278 (cs278, 48033);
    StopNeutralinoMap.insert(std::make_pair(mp278, mpp278));

    MassPoint mp279 (625, 150);
    StopCrossSection cs279 (0.00780081, 0.000768988);
    MassPointParameters mpp279 (cs279, 50014);
    StopNeutralinoMap.insert(std::make_pair(mp279, mpp279));

    MassPoint mp280 (625, 175);
    StopCrossSection cs280 (0.00780081, 0.000768988);
    MassPointParameters mpp280 (cs280, 44179);
    StopNeutralinoMap.insert(std::make_pair(mp280, mpp280));

    MassPoint mp281 (625, 200);
    StopCrossSection cs281 (0.00780081, 0.000768988);
    MassPointParameters mpp281 (cs281, 50260);
    StopNeutralinoMap.insert(std::make_pair(mp281, mpp281));

    MassPoint mp282 (625, 225);
    StopCrossSection cs282 (0.00780081, 0.000768988);
    MassPointParameters mpp282 (cs282, 50056);
    StopNeutralinoMap.insert(std::make_pair(mp282, mpp282));

    MassPoint mp283 (625, 250);
    StopCrossSection cs283 (0.00780081, 0.000768988);
    MassPointParameters mpp283 (cs283, 48346);
    StopNeutralinoMap.insert(std::make_pair(mp283, mpp283));

    MassPoint mp284 (625, 275);
    StopCrossSection cs284 (0.00780081, 0.000768988);
    MassPointParameters mpp284 (cs284, 51145);
    StopNeutralinoMap.insert(std::make_pair(mp284, mpp284));

    MassPoint mp285 (625, 300);
    StopCrossSection cs285 (0.00780081, 0.000768988);
    MassPointParameters mpp285 (cs285, 51222);
    StopNeutralinoMap.insert(std::make_pair(mp285, mpp285));

    MassPoint mp286 (625, 325);
    StopCrossSection cs286 (0.00780081, 0.000768988);
    MassPointParameters mpp286 (cs286, 53082);
    StopNeutralinoMap.insert(std::make_pair(mp286, mpp286));

    MassPoint mp287 (625, 350);
    StopCrossSection cs287 (0.00780081, 0.000768988);
    MassPointParameters mpp287 (cs287, 48483);
    StopNeutralinoMap.insert(std::make_pair(mp287, mpp287));

    MassPoint mp288 (625, 375);
    StopCrossSection cs288 (0.00780081, 0.000768988);
    MassPointParameters mpp288 (cs288, 49404);
    StopNeutralinoMap.insert(std::make_pair(mp288, mpp288));

    MassPoint mp289 (625, 400);
    StopCrossSection cs289 (0.00780081, 0.000768988);
    MassPointParameters mpp289 (cs289, 49466);
    StopNeutralinoMap.insert(std::make_pair(mp289, mpp289));

    MassPoint mp290 (625, 425);
    StopCrossSection cs290 (0.00780081, 0.000768988);
    MassPointParameters mpp290 (cs290, 49460);
    StopNeutralinoMap.insert(std::make_pair(mp290, mpp290));

    MassPoint mp291 (625, 450);
    StopCrossSection cs291 (0.00780081, 0.000768988);
    MassPointParameters mpp291 (cs291, 52267);
    StopNeutralinoMap.insert(std::make_pair(mp291, mpp291));

    MassPoint mp292 (625, 475);
    StopCrossSection cs292 (0.00780081, 0.000768988);
    MassPointParameters mpp292 (cs292, 48652);
    StopNeutralinoMap.insert(std::make_pair(mp292, mpp292));

    MassPoint mp293 (625, 500);
    StopCrossSection cs293 (0.00780081, 0.000768988);
    MassPointParameters mpp293 (cs293, 44152);
    StopNeutralinoMap.insert(std::make_pair(mp293, mpp293));

    MassPoint mp294 (625, 525);
    StopCrossSection cs294 (0.00780081, 0.000768988);
    MassPointParameters mpp294 (cs294, 50810);
    StopNeutralinoMap.insert(std::make_pair(mp294, mpp294));

    MassPoint mp295 (625, 550);
    StopCrossSection cs295 (0.00780081, 0.000768988);
    MassPointParameters mpp295 (cs295, 46070);
    StopNeutralinoMap.insert(std::make_pair(mp295, mpp295));

    MassPoint mp296 (625, 575);
    StopCrossSection cs296 (0.00780081, 0.000768988);
    MassPointParameters mpp296 (cs296, 50942);
    StopNeutralinoMap.insert(std::make_pair(mp296, mpp296));

    MassPoint mp297 (650, 1);
    StopCrossSection cs297 (0.00643244, 0.000638889);
    MassPointParameters mpp297 (cs297, 45847);
    StopNeutralinoMap.insert(std::make_pair(mp297, mpp297));

    MassPoint mp298 (650, 25);
    StopCrossSection cs298 (0.00643244, 0.000638889);
    MassPointParameters mpp298 (cs298, 53208);
    StopNeutralinoMap.insert(std::make_pair(mp298, mpp298));

    MassPoint mp299 (650, 50);
    StopCrossSection cs299 (0.00643244, 0.000638889);
    MassPointParameters mpp299 (cs299, 41186);
    StopNeutralinoMap.insert(std::make_pair(mp299, mpp299));

    MassPoint mp300 (650, 75);
    StopCrossSection cs300 (0.00643244, 0.000638889);
    MassPointParameters mpp300 (cs300, 43007);
    StopNeutralinoMap.insert(std::make_pair(mp300, mpp300));

    MassPoint mp301 (650, 100);
    StopCrossSection cs301 (0.00643244, 0.000638889);
    MassPointParameters mpp301 (cs301, 52006);
    StopNeutralinoMap.insert(std::make_pair(mp301, mpp301));

    MassPoint mp302 (650, 125);
    StopCrossSection cs302 (0.00643244, 0.000638889);
    MassPointParameters mpp302 (cs302, 48256);
    StopNeutralinoMap.insert(std::make_pair(mp302, mpp302));

    MassPoint mp303 (650, 150);
    StopCrossSection cs303 (0.00643244, 0.000638889);
    MassPointParameters mpp303 (cs303, 47606);
    StopNeutralinoMap.insert(std::make_pair(mp303, mpp303));

    MassPoint mp304 (650, 175);
    StopCrossSection cs304 (0.00643244, 0.000638889);
    MassPointParameters mpp304 (cs304, 47625);
    StopNeutralinoMap.insert(std::make_pair(mp304, mpp304));

    MassPoint mp305 (650, 200);
    StopCrossSection cs305 (0.00643244, 0.000638889);
    MassPointParameters mpp305 (cs305, 45379);
    StopNeutralinoMap.insert(std::make_pair(mp305, mpp305));

    MassPoint mp306 (650, 225);
    StopCrossSection cs306 (0.00643244, 0.000638889);
    MassPointParameters mpp306 (cs306, 54097);
    StopNeutralinoMap.insert(std::make_pair(mp306, mpp306));

    MassPoint mp307 (650, 250);
    StopCrossSection cs307 (0.00643244, 0.000638889);
    MassPointParameters mpp307 (cs307, 52543);
    StopNeutralinoMap.insert(std::make_pair(mp307, mpp307));

    MassPoint mp308 (650, 275);
    StopCrossSection cs308 (0.00643244, 0.000638889);
    MassPointParameters mpp308 (cs308, 50666);
    StopNeutralinoMap.insert(std::make_pair(mp308, mpp308));

    MassPoint mp309 (650, 300);
    StopCrossSection cs309 (0.00643244, 0.000638889);
    MassPointParameters mpp309 (cs309, 46964);
    StopNeutralinoMap.insert(std::make_pair(mp309, mpp309));

    MassPoint mp310 (650, 325);
    StopCrossSection cs310 (0.00643244, 0.000638889);
    MassPointParameters mpp310 (cs310, 43878);
    StopNeutralinoMap.insert(std::make_pair(mp310, mpp310));

    MassPoint mp311 (650, 350);
    StopCrossSection cs311 (0.00643244, 0.000638889);
    MassPointParameters mpp311 (cs311, 48494);
    StopNeutralinoMap.insert(std::make_pair(mp311, mpp311));

    MassPoint mp312 (650, 375);
    StopCrossSection cs312 (0.00643244, 0.000638889);
    MassPointParameters mpp312 (cs312, 52693);
    StopNeutralinoMap.insert(std::make_pair(mp312, mpp312));

    MassPoint mp313 (650, 400);
    StopCrossSection cs313 (0.00643244, 0.000638889);
    MassPointParameters mpp313 (cs313, 46197);
    StopNeutralinoMap.insert(std::make_pair(mp313, mpp313));

    MassPoint mp314 (650, 425);
    StopCrossSection cs314 (0.00643244, 0.000638889);
    MassPointParameters mpp314 (cs314, 49678);
    StopNeutralinoMap.insert(std::make_pair(mp314, mpp314));

    MassPoint mp315 (650, 450);
    StopCrossSection cs315 (0.00643244, 0.000638889);
    MassPointParameters mpp315 (cs315, 51026);
    StopNeutralinoMap.insert(std::make_pair(mp315, mpp315));

    MassPoint mp316 (650, 475);
    StopCrossSection cs316 (0.00643244, 0.000638889);
    MassPointParameters mpp316 (cs316, 42480);
    StopNeutralinoMap.insert(std::make_pair(mp316, mpp316));

    MassPoint mp317 (650, 500);
    StopCrossSection cs317 (0.00643244, 0.000638889);
    MassPointParameters mpp317 (cs317, 51976);
    StopNeutralinoMap.insert(std::make_pair(mp317, mpp317));

    MassPoint mp318 (650, 525);
    StopCrossSection cs318 (0.00643244, 0.000638889);
    MassPointParameters mpp318 (cs318, 43650);
    StopNeutralinoMap.insert(std::make_pair(mp318, mpp318));

    MassPoint mp319 (650, 550);
    StopCrossSection cs319 (0.00643244, 0.000638889);
    MassPointParameters mpp319 (cs319, 48155);
    StopNeutralinoMap.insert(std::make_pair(mp319, mpp319));

    MassPoint mp320 (650, 575);
    StopCrossSection cs320 (0.00643244, 0.000638889);
    MassPointParameters mpp320 (cs320, 50770);
    StopNeutralinoMap.insert(std::make_pair(mp320, mpp320));

    MassPoint mp321 (650, 600);
    StopCrossSection cs321 (0.00643244, 0.000638889);
    MassPointParameters mpp321 (cs321, 50832);
    StopNeutralinoMap.insert(std::make_pair(mp321, mpp321));

    MassPoint mp322 (675, 1);
    StopCrossSection cs322 (0.00533642, 0.000541519);
    MassPointParameters mpp322 (cs322, 47268);
    StopNeutralinoMap.insert(std::make_pair(mp322, mpp322));

    MassPoint mp323 (675, 25);
    StopCrossSection cs323 (0.00533642, 0.000541519);
    MassPointParameters mpp323 (cs323, 47286);
    StopNeutralinoMap.insert(std::make_pair(mp323, mpp323));

    MassPoint mp324 (675, 50);
    StopCrossSection cs324 (0.00533642, 0.000541519);
    MassPointParameters mpp324 (cs324, 50388);
    StopNeutralinoMap.insert(std::make_pair(mp324, mpp324));

    MassPoint mp325 (675, 75);
    StopCrossSection cs325 (0.00533642, 0.000541519);
    MassPointParameters mpp325 (cs325, 47859);
    StopNeutralinoMap.insert(std::make_pair(mp325, mpp325));

    MassPoint mp326 (675, 100);
    StopCrossSection cs326 (0.00533642, 0.000541519);
    MassPointParameters mpp326 (cs326, 49018);
    StopNeutralinoMap.insert(std::make_pair(mp326, mpp326));

    MassPoint mp327 (675, 125);
    StopCrossSection cs327 (0.00533642, 0.000541519);
    MassPointParameters mpp327 (cs327, 42176);
    StopNeutralinoMap.insert(std::make_pair(mp327, mpp327));

    MassPoint mp328 (675, 150);
    StopCrossSection cs328 (0.00533642, 0.000541519);
    MassPointParameters mpp328 (cs328, 49883);
    StopNeutralinoMap.insert(std::make_pair(mp328, mpp328));

    MassPoint mp329 (675, 175);
    StopCrossSection cs329 (0.00533642, 0.000541519);
    MassPointParameters mpp329 (cs329, 46186);
    StopNeutralinoMap.insert(std::make_pair(mp329, mpp329));

    MassPoint mp330 (675, 200);
    StopCrossSection cs330 (0.00533642, 0.000541519);
    MassPointParameters mpp330 (cs330, 51010);
    StopNeutralinoMap.insert(std::make_pair(mp330, mpp330));

    MassPoint mp331 (675, 225);
    StopCrossSection cs331 (0.00533642, 0.000541519);
    MassPointParameters mpp331 (cs331, 47734);
    StopNeutralinoMap.insert(std::make_pair(mp331, mpp331));

    MassPoint mp332 (675, 250);
    StopCrossSection cs332 (0.00533642, 0.000541519);
    MassPointParameters mpp332 (cs332, 45208);
    StopNeutralinoMap.insert(std::make_pair(mp332, mpp332));

    MassPoint mp333 (675, 275);
    StopCrossSection cs333 (0.00533642, 0.000541519);
    MassPointParameters mpp333 (cs333, 48475);
    StopNeutralinoMap.insert(std::make_pair(mp333, mpp333));

    MassPoint mp334 (675, 300);
    StopCrossSection cs334 (0.00533642, 0.000541519);
    MassPointParameters mpp334 (cs334, 54957);
    StopNeutralinoMap.insert(std::make_pair(mp334, mpp334));

    MassPoint mp335 (675, 325);
    StopCrossSection cs335 (0.00533642, 0.000541519);
    MassPointParameters mpp335 (cs335, 47373);
    StopNeutralinoMap.insert(std::make_pair(mp335, mpp335));

    MassPoint mp336 (675, 350);
    StopCrossSection cs336 (0.00533642, 0.000541519);
    MassPointParameters mpp336 (cs336, 50675);
    StopNeutralinoMap.insert(std::make_pair(mp336, mpp336));

    MassPoint mp337 (675, 375);
    StopCrossSection cs337 (0.00533642, 0.000541519);
    MassPointParameters mpp337 (cs337, 41396);
    StopNeutralinoMap.insert(std::make_pair(mp337, mpp337));

    MassPoint mp338 (675, 400);
    StopCrossSection cs338 (0.00533642, 0.000541519);
    MassPointParameters mpp338 (cs338, 51192);
    StopNeutralinoMap.insert(std::make_pair(mp338, mpp338));

    MassPoint mp339 (675, 425);
    StopCrossSection cs339 (0.00533642, 0.000541519);
    MassPointParameters mpp339 (cs339, 47809);
    StopNeutralinoMap.insert(std::make_pair(mp339, mpp339));

    MassPoint mp340 (675, 450);
    StopCrossSection cs340 (0.00533642, 0.000541519);
    MassPointParameters mpp340 (cs340, 47596);
    StopNeutralinoMap.insert(std::make_pair(mp340, mpp340));

    MassPoint mp341 (675, 475);
    StopCrossSection cs341 (0.00533642, 0.000541519);
    MassPointParameters mpp341 (cs341, 53884);
    StopNeutralinoMap.insert(std::make_pair(mp341, mpp341));

    MassPoint mp342 (675, 500);
    StopCrossSection cs342 (0.00533642, 0.000541519);
    MassPointParameters mpp342 (cs342, 50410);
    StopNeutralinoMap.insert(std::make_pair(mp342, mpp342));

    MassPoint mp343 (675, 525);
    StopCrossSection cs343 (0.00533642, 0.000541519);
    MassPointParameters mpp343 (cs343, 44762);
    StopNeutralinoMap.insert(std::make_pair(mp343, mpp343));

    MassPoint mp344 (675, 550);
    StopCrossSection cs344 (0.00533642, 0.000541519);
    MassPointParameters mpp344 (cs344, 46476);
    StopNeutralinoMap.insert(std::make_pair(mp344, mpp344));

    MassPoint mp345 (675, 575);
    StopCrossSection cs345 (0.00533642, 0.000541519);
    MassPointParameters mpp345 (cs345, 47077);
    StopNeutralinoMap.insert(std::make_pair(mp345, mpp345));

    MassPoint mp346 (675, 600);
    StopCrossSection cs346 (0.00533642, 0.000541519);
    MassPointParameters mpp346 (cs346, 51380);
    StopNeutralinoMap.insert(std::make_pair(mp346, mpp346));

    MassPoint mp347 (700, 1);
    StopCrossSection cs347 (0.0044387, 0.000457123);
    MassPointParameters mpp347 (cs347, 48092);
    StopNeutralinoMap.insert(std::make_pair(mp347, mpp347));

    MassPoint mp348 (700, 25);
    StopCrossSection cs348 (0.0044387, 0.000457123);
    MassPointParameters mpp348 (cs348, 50206);
    StopNeutralinoMap.insert(std::make_pair(mp348, mpp348));

    MassPoint mp349 (700, 50);
    StopCrossSection cs349 (0.0044387, 0.000457123);
    MassPointParameters mpp349 (cs349, 44996);
    StopNeutralinoMap.insert(std::make_pair(mp349, mpp349));

    MassPoint mp350 (700, 75);
    StopCrossSection cs350 (0.0044387, 0.000457123);
    MassPointParameters mpp350 (cs350, 43303);
    StopNeutralinoMap.insert(std::make_pair(mp350, mpp350));

    MassPoint mp351 (700, 100);
    StopCrossSection cs351 (0.0044387, 0.000457123);
    MassPointParameters mpp351 (cs351, 52318);
    StopNeutralinoMap.insert(std::make_pair(mp351, mpp351));

    MassPoint mp352 (700, 125);
    StopCrossSection cs352 (0.0044387, 0.000457123);
    MassPointParameters mpp352 (cs352, 52874);
    StopNeutralinoMap.insert(std::make_pair(mp352, mpp352));

    MassPoint mp353 (700, 150);
    StopCrossSection cs353 (0.0044387, 0.000457123);
    MassPointParameters mpp353 (cs353, 50128);
    StopNeutralinoMap.insert(std::make_pair(mp353, mpp353));

    MassPoint mp354 (700, 175);
    StopCrossSection cs354 (0.0044387, 0.000457123);
    MassPointParameters mpp354 (cs354, 49619);
    StopNeutralinoMap.insert(std::make_pair(mp354, mpp354));

    MassPoint mp355 (700, 200);
    StopCrossSection cs355 (0.0044387, 0.000457123);
    MassPointParameters mpp355 (cs355, 49873);
    StopNeutralinoMap.insert(std::make_pair(mp355, mpp355));

    MassPoint mp356 (700, 225);
    StopCrossSection cs356 (0.0044387, 0.000457123);
    MassPointParameters mpp356 (cs356, 44017);
    StopNeutralinoMap.insert(std::make_pair(mp356, mpp356));

    MassPoint mp357 (700, 250);
    StopCrossSection cs357 (0.0044387, 0.000457123);
    MassPointParameters mpp357 (cs357, 48215);
    StopNeutralinoMap.insert(std::make_pair(mp357, mpp357));

    MassPoint mp358 (700, 275);
    StopCrossSection cs358 (0.0044387, 0.000457123);
    MassPointParameters mpp358 (cs358, 51344);
    StopNeutralinoMap.insert(std::make_pair(mp358, mpp358));

    MassPoint mp359 (700, 300);
    StopCrossSection cs359 (0.0044387, 0.000457123);
    MassPointParameters mpp359 (cs359, 54461);
    StopNeutralinoMap.insert(std::make_pair(mp359, mpp359));

    MassPoint mp360 (700, 325);
    StopCrossSection cs360 (0.0044387, 0.000457123);
    MassPointParameters mpp360 (cs360, 47025);
    StopNeutralinoMap.insert(std::make_pair(mp360, mpp360));

    MassPoint mp361 (700, 350);
    StopCrossSection cs361 (0.0044387, 0.000457123);
    MassPointParameters mpp361 (cs361, 53220);
    StopNeutralinoMap.insert(std::make_pair(mp361, mpp361));

    MassPoint mp362 (700, 375);
    StopCrossSection cs362 (0.0044387, 0.000457123);
    MassPointParameters mpp362 (cs362, 44983);
    StopNeutralinoMap.insert(std::make_pair(mp362, mpp362));

    MassPoint mp363 (700, 400);
    StopCrossSection cs363 (0.0044387, 0.000457123);
    MassPointParameters mpp363 (cs363, 40253);
    StopNeutralinoMap.insert(std::make_pair(mp363, mpp363));

    MassPoint mp364 (700, 425);
    StopCrossSection cs364 (0.0044387, 0.000457123);
    MassPointParameters mpp364 (cs364, 50819);
    StopNeutralinoMap.insert(std::make_pair(mp364, mpp364));

    MassPoint mp365 (700, 450);
    StopCrossSection cs365 (0.0044387, 0.000457123);
    MassPointParameters mpp365 (cs365, 48140);
    StopNeutralinoMap.insert(std::make_pair(mp365, mpp365));

    MassPoint mp366 (700, 475);
    StopCrossSection cs366 (0.0044387, 0.000457123);
    MassPointParameters mpp366 (cs366, 50851);
    StopNeutralinoMap.insert(std::make_pair(mp366, mpp366));

    MassPoint mp367 (700, 500);
    StopCrossSection cs367 (0.0044387, 0.000457123);
    MassPointParameters mpp367 (cs367, 44310);
    StopNeutralinoMap.insert(std::make_pair(mp367, mpp367));

    MassPoint mp368 (700, 525);
    StopCrossSection cs368 (0.0044387, 0.000457123);
    MassPointParameters mpp368 (cs368, 50292);
    StopNeutralinoMap.insert(std::make_pair(mp368, mpp368));

    MassPoint mp369 (700, 550);
    StopCrossSection cs369 (0.0044387, 0.000457123);
    MassPointParameters mpp369 (cs369, 50226);
    StopNeutralinoMap.insert(std::make_pair(mp369, mpp369));

    MassPoint mp370 (700, 575);
    StopCrossSection cs370 (0.0044387, 0.000457123);
    MassPointParameters mpp370 (cs370, 50771);
    StopNeutralinoMap.insert(std::make_pair(mp370, mpp370));

    MassPoint mp371 (700, 600);
    StopCrossSection cs371 (0.0044387, 0.000457123);
    MassPointParameters mpp371 (cs371, 49120);
    StopNeutralinoMap.insert(std::make_pair(mp371, mpp371));

    MassPoint mp372 (725, 1);
    StopCrossSection cs372 (0.00370675, 0.000385799);
    MassPointParameters mpp372 (cs372, 49654);
    StopNeutralinoMap.insert(std::make_pair(mp372, mpp372));

    MassPoint mp373 (725, 25);
    StopCrossSection cs373 (0.00370675, 0.000385799);
    MassPointParameters mpp373 (cs373, 54279);
    StopNeutralinoMap.insert(std::make_pair(mp373, mpp373));

    MassPoint mp374 (725, 50);
    StopCrossSection cs374 (0.00370675, 0.000385799);
    MassPointParameters mpp374 (cs374, 43911);
    StopNeutralinoMap.insert(std::make_pair(mp374, mpp374));

    MassPoint mp375 (725, 75);
    StopCrossSection cs375 (0.00370675, 0.000385799);
    MassPointParameters mpp375 (cs375, 48816);
    StopNeutralinoMap.insert(std::make_pair(mp375, mpp375));

    MassPoint mp376 (725, 100);
    StopCrossSection cs376 (0.00370675, 0.000385799);
    MassPointParameters mpp376 (cs376, 52643);
    StopNeutralinoMap.insert(std::make_pair(mp376, mpp376));

    MassPoint mp377 (725, 125);
    StopCrossSection cs377 (0.00370675, 0.000385799);
    MassPointParameters mpp377 (cs377, 50389);
    StopNeutralinoMap.insert(std::make_pair(mp377, mpp377));

    MassPoint mp378 (725, 150);
    StopCrossSection cs378 (0.00370675, 0.000385799);
    MassPointParameters mpp378 (cs378, 53972);
    StopNeutralinoMap.insert(std::make_pair(mp378, mpp378));

    MassPoint mp379 (725, 175);
    StopCrossSection cs379 (0.00370675, 0.000385799);
    MassPointParameters mpp379 (cs379, 52198);
    StopNeutralinoMap.insert(std::make_pair(mp379, mpp379));

    MassPoint mp380 (725, 200);
    StopCrossSection cs380 (0.00370675, 0.000385799);
    MassPointParameters mpp380 (cs380, 47309);
    StopNeutralinoMap.insert(std::make_pair(mp380, mpp380));

    MassPoint mp381 (725, 225);
    StopCrossSection cs381 (0.00370675, 0.000385799);
    MassPointParameters mpp381 (cs381, 51338);
    StopNeutralinoMap.insert(std::make_pair(mp381, mpp381));

    MassPoint mp382 (725, 250);
    StopCrossSection cs382 (0.00370675, 0.000385799);
    MassPointParameters mpp382 (cs382, 46512);
    StopNeutralinoMap.insert(std::make_pair(mp382, mpp382));

    MassPoint mp383 (725, 275);
    StopCrossSection cs383 (0.00370675, 0.000385799);
    MassPointParameters mpp383 (cs383, 44440);
    StopNeutralinoMap.insert(std::make_pair(mp383, mpp383));

    MassPoint mp384 (725, 300);
    StopCrossSection cs384 (0.00370675, 0.000385799);
    MassPointParameters mpp384 (cs384, 47354);
    StopNeutralinoMap.insert(std::make_pair(mp384, mpp384));

    MassPoint mp385 (725, 325);
    StopCrossSection cs385 (0.00370675, 0.000385799);
    MassPointParameters mpp385 (cs385, 53844);
    StopNeutralinoMap.insert(std::make_pair(mp385, mpp385));

    MassPoint mp386 (725, 350);
    StopCrossSection cs386 (0.00370675, 0.000385799);
    MassPointParameters mpp386 (cs386, 50538);
    StopNeutralinoMap.insert(std::make_pair(mp386, mpp386));

    MassPoint mp387 (725, 375);
    StopCrossSection cs387 (0.00370675, 0.000385799);
    MassPointParameters mpp387 (cs387, 51728);
    StopNeutralinoMap.insert(std::make_pair(mp387, mpp387));

    MassPoint mp388 (725, 400);
    StopCrossSection cs388 (0.00370675, 0.000385799);
    MassPointParameters mpp388 (cs388, 46897);
    StopNeutralinoMap.insert(std::make_pair(mp388, mpp388));

    MassPoint mp389 (725, 425);
    StopCrossSection cs389 (0.00370675, 0.000385799);
    MassPointParameters mpp389 (cs389, 49736);
    StopNeutralinoMap.insert(std::make_pair(mp389, mpp389));

    MassPoint mp390 (725, 450);
    StopCrossSection cs390 (0.00370675, 0.000385799);
    MassPointParameters mpp390 (cs390, 56034);
    StopNeutralinoMap.insert(std::make_pair(mp390, mpp390));

    MassPoint mp391 (725, 475);
    StopCrossSection cs391 (0.00370675, 0.000385799);
    MassPointParameters mpp391 (cs391, 54963);
    StopNeutralinoMap.insert(std::make_pair(mp391, mpp391));

    MassPoint mp392 (725, 500);
    StopCrossSection cs392 (0.00370675, 0.000385799);
    MassPointParameters mpp392 (cs392, 47607);
    StopNeutralinoMap.insert(std::make_pair(mp392, mpp392));

    MassPoint mp393 (725, 525);
    StopCrossSection cs393 (0.00370675, 0.000385799);
    MassPointParameters mpp393 (cs393, 51313);
    StopNeutralinoMap.insert(std::make_pair(mp393, mpp393));

    MassPoint mp394 (725, 550);
    StopCrossSection cs394 (0.00370675, 0.000385799);
    MassPointParameters mpp394 (cs394, 51166);
    StopNeutralinoMap.insert(std::make_pair(mp394, mpp394));

    MassPoint mp395 (725, 575);
    StopCrossSection cs395 (0.00370675, 0.000385799);
    MassPointParameters mpp395 (cs395, 54147);
    StopNeutralinoMap.insert(std::make_pair(mp395, mpp395));

    MassPoint mp396 (725, 600);
    StopCrossSection cs396 (0.00370675, 0.000385799);
    MassPointParameters mpp396 (cs396, 47308);
    StopNeutralinoMap.insert(std::make_pair(mp396, mpp396));

    MassPoint mp397 (750, 1);
    StopCrossSection cs397 (0.00310861, 0.000330353);
    MassPointParameters mpp397 (cs397, 44838);
    StopNeutralinoMap.insert(std::make_pair(mp397, mpp397));

    MassPoint mp398 (750, 25);
    StopCrossSection cs398 (0.00310861, 0.000330353);
    MassPointParameters mpp398 (cs398, 46064);
    StopNeutralinoMap.insert(std::make_pair(mp398, mpp398));

    MassPoint mp399 (750, 50);
    StopCrossSection cs399 (0.00310861, 0.000330353);
    MassPointParameters mpp399 (cs399, 50667);
    StopNeutralinoMap.insert(std::make_pair(mp399, mpp399));

    MassPoint mp400 (750, 75);
    StopCrossSection cs400 (0.00310861, 0.000330353);
    MassPointParameters mpp400 (cs400, 49931);
    StopNeutralinoMap.insert(std::make_pair(mp400, mpp400));

    MassPoint mp401 (750, 100);
    StopCrossSection cs401 (0.00310861, 0.000330353);
    MassPointParameters mpp401 (cs401, 47500);
    StopNeutralinoMap.insert(std::make_pair(mp401, mpp401));

    MassPoint mp402 (750, 125);
    StopCrossSection cs402 (0.00310861, 0.000330353);
    MassPointParameters mpp402 (cs402, 47999);
    StopNeutralinoMap.insert(std::make_pair(mp402, mpp402));

    MassPoint mp403 (750, 150);
    StopCrossSection cs403 (0.00310861, 0.000330353);
    MassPointParameters mpp403 (cs403, 54054);
    StopNeutralinoMap.insert(std::make_pair(mp403, mpp403));

    MassPoint mp404 (750, 175);
    StopCrossSection cs404 (0.00310861, 0.000330353);
    MassPointParameters mpp404 (cs404, 51175);
    StopNeutralinoMap.insert(std::make_pair(mp404, mpp404));

    MassPoint mp405 (750, 200);
    StopCrossSection cs405 (0.00310861, 0.000330353);
    MassPointParameters mpp405 (cs405, 49615);
    StopNeutralinoMap.insert(std::make_pair(mp405, mpp405));

    MassPoint mp406 (750, 225);
    StopCrossSection cs406 (0.00310861, 0.000330353);
    MassPointParameters mpp406 (cs406, 44755);
    StopNeutralinoMap.insert(std::make_pair(mp406, mpp406));

    MassPoint mp407 (750, 250);
    StopCrossSection cs407 (0.00310861, 0.000330353);
    MassPointParameters mpp407 (cs407, 47874);
    StopNeutralinoMap.insert(std::make_pair(mp407, mpp407));

    MassPoint mp408 (750, 275);
    StopCrossSection cs408 (0.00310861, 0.000330353);
    MassPointParameters mpp408 (cs408, 49319);
    StopNeutralinoMap.insert(std::make_pair(mp408, mpp408));

    MassPoint mp409 (750, 300);
    StopCrossSection cs409 (0.00310861, 0.000330353);
    MassPointParameters mpp409 (cs409, 51223);
    StopNeutralinoMap.insert(std::make_pair(mp409, mpp409));

    MassPoint mp410 (750, 325);
    StopCrossSection cs410 (0.00310861, 0.000330353);
    MassPointParameters mpp410 (cs410, 46208);
    StopNeutralinoMap.insert(std::make_pair(mp410, mpp410));

    MassPoint mp411 (750, 350);
    StopCrossSection cs411 (0.00310861, 0.000330353);
    MassPointParameters mpp411 (cs411, 51153);
    StopNeutralinoMap.insert(std::make_pair(mp411, mpp411));

    MassPoint mp412 (750, 375);
    StopCrossSection cs412 (0.00310861, 0.000330353);
    MassPointParameters mpp412 (cs412, 54547);
    StopNeutralinoMap.insert(std::make_pair(mp412, mpp412));

    MassPoint mp413 (750, 400);
    StopCrossSection cs413 (0.00310861, 0.000330353);
    MassPointParameters mpp413 (cs413, 45297);
    StopNeutralinoMap.insert(std::make_pair(mp413, mpp413));

    MassPoint mp414 (750, 425);
    StopCrossSection cs414 (0.00310861, 0.000330353);
    MassPointParameters mpp414 (cs414, 48519);
    StopNeutralinoMap.insert(std::make_pair(mp414, mpp414));

    MassPoint mp415 (750, 450);
    StopCrossSection cs415 (0.00310861, 0.000330353);
    MassPointParameters mpp415 (cs415, 42164);
    StopNeutralinoMap.insert(std::make_pair(mp415, mpp415));

    MassPoint mp416 (750, 475);
    StopCrossSection cs416 (0.00310861, 0.000330353);
    MassPointParameters mpp416 (cs416, 46163);
    StopNeutralinoMap.insert(std::make_pair(mp416, mpp416));

    MassPoint mp417 (750, 500);
    StopCrossSection cs417 (0.00310861, 0.000330353);
    MassPointParameters mpp417 (cs417, 48119);
    StopNeutralinoMap.insert(std::make_pair(mp417, mpp417));

    MassPoint mp418 (750, 525);
    StopCrossSection cs418 (0.00310861, 0.000330353);
    MassPointParameters mpp418 (cs418, 52644);
    StopNeutralinoMap.insert(std::make_pair(mp418, mpp418));

    MassPoint mp419 (750, 550);
    StopCrossSection cs419 (0.00310861, 0.000330353);
    MassPointParameters mpp419 (cs419, 51787);
    StopNeutralinoMap.insert(std::make_pair(mp419, mpp419));

    MassPoint mp420 (750, 575);
    StopCrossSection cs420 (0.00310861, 0.000330353);
    MassPointParameters mpp420 (cs420, 52876);
    StopNeutralinoMap.insert(std::make_pair(mp420, mpp420));

    MassPoint mp421 (750, 600);
    StopCrossSection cs421 (0.00310861, 0.000330353);
    MassPointParameters mpp421 (cs421, 50870);
    StopNeutralinoMap.insert(std::make_pair(mp421, mpp421));

    MassPoint mp422 (775, 1);
    StopCrossSection cs422 (0.00261656, 0.000283139);
    MassPointParameters mpp422 (cs422, 53007);
    StopNeutralinoMap.insert(std::make_pair(mp422, mpp422));

    MassPoint mp423 (775, 25);
    StopCrossSection cs423 (0.00261656, 0.000283139);
    MassPointParameters mpp423 (cs423, 46246);
    StopNeutralinoMap.insert(std::make_pair(mp423, mpp423));

    MassPoint mp424 (775, 50);
    StopCrossSection cs424 (0.00261656, 0.000283139);
    MassPointParameters mpp424 (cs424, 47856);
    StopNeutralinoMap.insert(std::make_pair(mp424, mpp424));

    MassPoint mp425 (775, 75);
    StopCrossSection cs425 (0.00261656, 0.000283139);
    MassPointParameters mpp425 (cs425, 48268);
    StopNeutralinoMap.insert(std::make_pair(mp425, mpp425));

    MassPoint mp426 (775, 100);
    StopCrossSection cs426 (0.00261656, 0.000283139);
    MassPointParameters mpp426 (cs426, 51294);
    StopNeutralinoMap.insert(std::make_pair(mp426, mpp426));

    MassPoint mp427 (775, 125);
    StopCrossSection cs427 (0.00261656, 0.000283139);
    MassPointParameters mpp427 (cs427, 48224);
    StopNeutralinoMap.insert(std::make_pair(mp427, mpp427));

    MassPoint mp428 (775, 150);
    StopCrossSection cs428 (0.00261656, 0.000283139);
    MassPointParameters mpp428 (cs428, 48545);
    StopNeutralinoMap.insert(std::make_pair(mp428, mpp428));

    MassPoint mp429 (775, 175);
    StopCrossSection cs429 (0.00261656, 0.000283139);
    MassPointParameters mpp429 (cs429, 50181);
    StopNeutralinoMap.insert(std::make_pair(mp429, mpp429));

    MassPoint mp430 (775, 200);
    StopCrossSection cs430 (0.00261656, 0.000283139);
    MassPointParameters mpp430 (cs430, 49232);
    StopNeutralinoMap.insert(std::make_pair(mp430, mpp430));

    MassPoint mp431 (775, 225);
    StopCrossSection cs431 (0.00261656, 0.000283139);
    MassPointParameters mpp431 (cs431, 48471);
    StopNeutralinoMap.insert(std::make_pair(mp431, mpp431));

    MassPoint mp432 (775, 250);
    StopCrossSection cs432 (0.00261656, 0.000283139);
    MassPointParameters mpp432 (cs432, 47661);
    StopNeutralinoMap.insert(std::make_pair(mp432, mpp432));

    MassPoint mp433 (775, 275);
    StopCrossSection cs433 (0.00261656, 0.000283139);
    MassPointParameters mpp433 (cs433, 53054);
    StopNeutralinoMap.insert(std::make_pair(mp433, mpp433));

    MassPoint mp434 (775, 300);
    StopCrossSection cs434 (0.00261656, 0.000283139);
    MassPointParameters mpp434 (cs434, 57102);
    StopNeutralinoMap.insert(std::make_pair(mp434, mpp434));

    MassPoint mp435 (775, 325);
    StopCrossSection cs435 (0.00261656, 0.000283139);
    MassPointParameters mpp435 (cs435, 54810);
    StopNeutralinoMap.insert(std::make_pair(mp435, mpp435));

    MassPoint mp436 (775, 350);
    StopCrossSection cs436 (0.00261656, 0.000283139);
    MassPointParameters mpp436 (cs436, 49821);
    StopNeutralinoMap.insert(std::make_pair(mp436, mpp436));

    MassPoint mp437 (775, 375);
    StopCrossSection cs437 (0.00261656, 0.000283139);
    MassPointParameters mpp437 (cs437, 46639);
    StopNeutralinoMap.insert(std::make_pair(mp437, mpp437));

    MassPoint mp438 (775, 400);
    StopCrossSection cs438 (0.00261656, 0.000283139);
    MassPointParameters mpp438 (cs438, 45677);
    StopNeutralinoMap.insert(std::make_pair(mp438, mpp438));

    MassPoint mp439 (775, 425);
    StopCrossSection cs439 (0.00261656, 0.000283139);
    MassPointParameters mpp439 (cs439, 48957);
    StopNeutralinoMap.insert(std::make_pair(mp439, mpp439));

    MassPoint mp440 (775, 450);
    StopCrossSection cs440 (0.00261656, 0.000283139);
    MassPointParameters mpp440 (cs440, 50332);
    StopNeutralinoMap.insert(std::make_pair(mp440, mpp440));

    MassPoint mp441 (775, 475);
    StopCrossSection cs441 (0.00261656, 0.000283139);
    MassPointParameters mpp441 (cs441, 42853);
    StopNeutralinoMap.insert(std::make_pair(mp441, mpp441));

    MassPoint mp442 (775, 500);
    StopCrossSection cs442 (0.00261656, 0.000283139);
    MassPointParameters mpp442 (cs442, 48317);
    StopNeutralinoMap.insert(std::make_pair(mp442, mpp442));

    MassPoint mp443 (775, 525);
    StopCrossSection cs443 (0.00261656, 0.000283139);
    MassPointParameters mpp443 (cs443, 51120);
    StopNeutralinoMap.insert(std::make_pair(mp443, mpp443));

    MassPoint mp444 (775, 550);
    StopCrossSection cs444 (0.00261656, 0.000283139);
    MassPointParameters mpp444 (cs444, 44471);
    StopNeutralinoMap.insert(std::make_pair(mp444, mpp444));

    MassPoint mp445 (775, 575);
    StopCrossSection cs445 (0.00261656, 0.000283139);
    MassPointParameters mpp445 (cs445, 48917);
    StopNeutralinoMap.insert(std::make_pair(mp445, mpp445));

    MassPoint mp446 (775, 600);
    StopCrossSection cs446 (0.00261656, 0.000283139);
    MassPointParameters mpp446 (cs446, 54290);
    StopNeutralinoMap.insert(std::make_pair(mp446, mpp446));

    MassPoint mp447 (800, 1);
    StopCrossSection cs447 (0.00221197, 0.000245196);
    MassPointParameters mpp447 (cs447, 49747);
    StopNeutralinoMap.insert(std::make_pair(mp447, mpp447));

    MassPoint mp448 (800, 25);
    StopCrossSection cs448 (0.00221197, 0.000245196);
    MassPointParameters mpp448 (cs448, 47776);
    StopNeutralinoMap.insert(std::make_pair(mp448, mpp448));

    MassPoint mp449 (800, 50);
    StopCrossSection cs449 (0.00221197, 0.000245196);
    MassPointParameters mpp449 (cs449, 51222);
    StopNeutralinoMap.insert(std::make_pair(mp449, mpp449));

    MassPoint mp450 (800, 75);
    StopCrossSection cs450 (0.00221197, 0.000245196);
    MassPointParameters mpp450 (cs450, 48427);
    StopNeutralinoMap.insert(std::make_pair(mp450, mpp450));

    MassPoint mp451 (800, 100);
    StopCrossSection cs451 (0.00221197, 0.000245196);
    MassPointParameters mpp451 (cs451, 47709);
    StopNeutralinoMap.insert(std::make_pair(mp451, mpp451));

    MassPoint mp452 (800, 125);
    StopCrossSection cs452 (0.00221197, 0.000245196);
    MassPointParameters mpp452 (cs452, 47551);
    StopNeutralinoMap.insert(std::make_pair(mp452, mpp452));

    MassPoint mp453 (800, 150);
    StopCrossSection cs453 (0.00221197, 0.000245196);
    MassPointParameters mpp453 (cs453, 50106);
    StopNeutralinoMap.insert(std::make_pair(mp453, mpp453));

    MassPoint mp454 (800, 175);
    StopCrossSection cs454 (0.00221197, 0.000245196);
    MassPointParameters mpp454 (cs454, 55232);
    StopNeutralinoMap.insert(std::make_pair(mp454, mpp454));

    MassPoint mp455 (800, 200);
    StopCrossSection cs455 (0.00221197, 0.000245196);
    MassPointParameters mpp455 (cs455, 47616);
    StopNeutralinoMap.insert(std::make_pair(mp455, mpp455));

    MassPoint mp456 (800, 225);
    StopCrossSection cs456 (0.00221197, 0.000245196);
    MassPointParameters mpp456 (cs456, 54843);
    StopNeutralinoMap.insert(std::make_pair(mp456, mpp456));

    MassPoint mp457 (800, 250);
    StopCrossSection cs457 (0.00221197, 0.000245196);
    MassPointParameters mpp457 (cs457, 53711);
    StopNeutralinoMap.insert(std::make_pair(mp457, mpp457));

    MassPoint mp458 (800, 275);
    StopCrossSection cs458 (0.00221197, 0.000245196);
    MassPointParameters mpp458 (cs458, 48662);
    StopNeutralinoMap.insert(std::make_pair(mp458, mpp458));

    MassPoint mp459 (800, 300);
    StopCrossSection cs459 (0.00221197, 0.000245196);
    MassPointParameters mpp459 (cs459, 47712);
    StopNeutralinoMap.insert(std::make_pair(mp459, mpp459));

    MassPoint mp460 (800, 325);
    StopCrossSection cs460 (0.00221197, 0.000245196);
    MassPointParameters mpp460 (cs460, 45909);
    StopNeutralinoMap.insert(std::make_pair(mp460, mpp460));

    MassPoint mp461 (800, 350);
    StopCrossSection cs461 (0.00221197, 0.000245196);
    MassPointParameters mpp461 (cs461, 50099);
    StopNeutralinoMap.insert(std::make_pair(mp461, mpp461));

    MassPoint mp462 (800, 375);
    StopCrossSection cs462 (0.00221197, 0.000245196);
    MassPointParameters mpp462 (cs462, 52501);
    StopNeutralinoMap.insert(std::make_pair(mp462, mpp462));

    MassPoint mp463 (800, 400);
    StopCrossSection cs463 (0.00221197, 0.000245196);
    MassPointParameters mpp463 (cs463, 41825);
    StopNeutralinoMap.insert(std::make_pair(mp463, mpp463));

    MassPoint mp464 (800, 425);
    StopCrossSection cs464 (0.00221197, 0.000245196);
    MassPointParameters mpp464 (cs464, 49529);
    StopNeutralinoMap.insert(std::make_pair(mp464, mpp464));

    MassPoint mp465 (800, 450);
    StopCrossSection cs465 (0.00221197, 0.000245196);
    MassPointParameters mpp465 (cs465, 52564);
    StopNeutralinoMap.insert(std::make_pair(mp465, mpp465));

    MassPoint mp466 (800, 475);
    StopCrossSection cs466 (0.00221197, 0.000245196);
    MassPointParameters mpp466 (cs466, 48811);
    StopNeutralinoMap.insert(std::make_pair(mp466, mpp466));

    MassPoint mp467 (800, 500);
    StopCrossSection cs467 (0.00221197, 0.000245196);
    MassPointParameters mpp467 (cs467, 49303);
    StopNeutralinoMap.insert(std::make_pair(mp467, mpp467));

    MassPoint mp468 (800, 525);
    StopCrossSection cs468 (0.00221197, 0.000245196);
    MassPointParameters mpp468 (cs468, 44363);
    StopNeutralinoMap.insert(std::make_pair(mp468, mpp468));

    MassPoint mp469 (800, 550);
    StopCrossSection cs469 (0.00221197, 0.000245196);
    MassPointParameters mpp469 (cs469, 51626);
    StopNeutralinoMap.insert(std::make_pair(mp469, mpp469));

    MassPoint mp470 (800, 575);
    StopCrossSection cs470 (0.00221197, 0.000245196);
    MassPointParameters mpp470 (cs470, 51956);
    StopNeutralinoMap.insert(std::make_pair(mp470, mpp470));

    MassPoint mp471 (800, 600);
    StopCrossSection cs471 (0.00221197, 0.000245196);
    MassPointParameters mpp471 (cs471, 48690);
    StopNeutralinoMap.insert(std::make_pair(mp471, mpp471));

    MassPoint mpisr0 (100, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 0.979484));

    MassPoint mpisr1 (100, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 0.979383));

    MassPoint mpisr2 (100, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 0.979387));

    MassPoint mpisr3 (125, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 0.976002));

    MassPoint mpisr4 (125, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 0.976321));

    MassPoint mpisr5 (125, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 0.976085));

    MassPoint mpisr6 (125, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 0.976627));

    MassPoint mpisr7 (150, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 0.974469));

    MassPoint mpisr8 (150, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 0.974205));

    MassPoint mpisr9 (150, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 0.973861));

    MassPoint mpisr10 (150, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 0.974392));

    MassPoint mpisr11 (150, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 0.974044));

    MassPoint mpisr12 (175, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 0.972008));

    MassPoint mpisr13 (175, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 0.971918));

    MassPoint mpisr14 (175, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 0.97162));

    MassPoint mpisr15 (175, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 0.971737));

    MassPoint mpisr16 (175, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 0.971931));

    MassPoint mpisr17 (175, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 0.972129));

    MassPoint mpisr18 (200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 0.970492));

    MassPoint mpisr19 (200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 0.970458));

    MassPoint mpisr20 (200, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 0.970609));

    MassPoint mpisr21 (200, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 0.970983));

    MassPoint mpisr22 (200, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 0.970525));

    MassPoint mpisr23 (200, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 0.970592));

    MassPoint mpisr24 (200, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 0.970573));

    MassPoint mpisr25 (225, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 0.969327));

    MassPoint mpisr26 (225, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 0.969657));

    MassPoint mpisr27 (225, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 0.970022));

    MassPoint mpisr28 (225, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 0.969366));

    MassPoint mpisr29 (225, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 0.969095));

    MassPoint mpisr30 (225, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 0.96997));

    MassPoint mpisr31 (225, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 0.969576));

    MassPoint mpisr32 (225, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 0.969564));

    MassPoint mpisr33 (250, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 0.968625));

    MassPoint mpisr34 (250, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 0.968521));

    MassPoint mpisr35 (250, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 0.968577));

    MassPoint mpisr36 (250, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 0.968843));

    MassPoint mpisr37 (250, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 0.969211));

    MassPoint mpisr38 (250, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 0.968396));

    MassPoint mpisr39 (250, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 0.968383));

    MassPoint mpisr40 (250, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 0.968677));

    MassPoint mpisr41 (250, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 0.968671));

    MassPoint mpisr42 (275, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 0.967603));

    MassPoint mpisr43 (275, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 0.968007));

    MassPoint mpisr44 (275, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 0.967501));

    MassPoint mpisr45 (275, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 0.968058));

    MassPoint mpisr46 (275, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 0.967587));

    MassPoint mpisr47 (275, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 0.967957));

    MassPoint mpisr48 (275, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 0.968358));

    MassPoint mpisr49 (275, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 0.96754));

    MassPoint mpisr50 (275, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 0.967782));

    MassPoint mpisr51 (275, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 0.967845));

    MassPoint mpisr52 (300, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 0.967236));

    MassPoint mpisr53 (300, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 0.966929));

    MassPoint mpisr54 (300, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 0.967609));

    MassPoint mpisr55 (300, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 0.967228));

    MassPoint mpisr56 (300, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 0.967423));

    MassPoint mpisr57 (300, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 0.967233));

    MassPoint mpisr58 (300, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 0.967102));

    MassPoint mpisr59 (300, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 0.967076));

    MassPoint mpisr60 (300, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 0.967169));

    MassPoint mpisr61 (300, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 0.96742));

    MassPoint mpisr62 (300, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 0.967101));

    MassPoint mpisr63 (325, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 0.96674));

    MassPoint mpisr64 (325, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 0.967027));

    MassPoint mpisr65 (325, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 0.966995));

    MassPoint mpisr66 (325, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 0.966774));

    MassPoint mpisr67 (325, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 0.966837));

    MassPoint mpisr68 (325, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 0.966847));

    MassPoint mpisr69 (325, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 0.966536));

    MassPoint mpisr70 (325, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 0.966709));

    MassPoint mpisr71 (325, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 0.967417));

    MassPoint mpisr72 (325, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 0.966682));

    MassPoint mpisr73 (325, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 0.966689));

    MassPoint mpisr74 (325, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 0.966755));

    MassPoint mpisr75 (350, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 0.966317));

    MassPoint mpisr76 (350, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 0.966768));

    MassPoint mpisr77 (350, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 0.966021));

    MassPoint mpisr78 (350, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 0.966345));

    MassPoint mpisr79 (350, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 0.966563));

    MassPoint mpisr80 (350, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 0.96658));

    MassPoint mpisr81 (350, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 0.966744));

    MassPoint mpisr82 (350, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 0.966812));

    MassPoint mpisr83 (350, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 0.966515));

    MassPoint mpisr84 (350, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 0.966357));

    MassPoint mpisr85 (350, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 0.966626));

    MassPoint mpisr86 (350, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 0.966774));

    MassPoint mpisr87 (350, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 0.966369));

    MassPoint mpisr88 (375, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 0.965879));

    MassPoint mpisr89 (375, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 0.965847));

    MassPoint mpisr90 (375, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 0.96669));

    MassPoint mpisr91 (375, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 0.966481));

    MassPoint mpisr92 (375, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 0.966439));

    MassPoint mpisr93 (375, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 0.966181));

    MassPoint mpisr94 (375, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 0.966429));

    MassPoint mpisr95 (375, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 0.966081));

    MassPoint mpisr96 (375, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 0.96611));

    MassPoint mpisr97 (375, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 0.966361));

    MassPoint mpisr98 (375, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 0.966157));

    MassPoint mpisr99 (375, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 0.966436));

    MassPoint mpisr100 (375, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 0.966332));

    MassPoint mpisr101 (375, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 0.966608));

    MassPoint mpisr102 (400, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 0.966496));

    MassPoint mpisr103 (400, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 0.96609));

    MassPoint mpisr104 (400, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 0.965977));

    MassPoint mpisr105 (400, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 0.96626));

    MassPoint mpisr106 (400, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 0.96623));

    MassPoint mpisr107 (400, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 0.966493));

    MassPoint mpisr108 (400, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 0.966043));

    MassPoint mpisr109 (400, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 0.966067));

    MassPoint mpisr110 (400, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 0.966145));

    MassPoint mpisr111 (400, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 0.966132));

    MassPoint mpisr112 (400, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 0.966293));

    MassPoint mpisr113 (400, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 0.965803));

    MassPoint mpisr114 (400, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 0.966021));

    MassPoint mpisr115 (400, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 0.966072));

    MassPoint mpisr116 (400, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 0.966014));

    MassPoint mpisr117 (425, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 0.966005));

    MassPoint mpisr118 (425, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 0.965973));

    MassPoint mpisr119 (425, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 0.96547));

    MassPoint mpisr120 (425, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 0.966185));

    MassPoint mpisr121 (425, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 0.965972));

    MassPoint mpisr122 (425, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 0.966));

    MassPoint mpisr123 (425, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 0.966132));

    MassPoint mpisr124 (425, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 0.965792));

    MassPoint mpisr125 (425, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 0.965788));

    MassPoint mpisr126 (425, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 0.966256));

    MassPoint mpisr127 (425, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 0.966428));

    MassPoint mpisr128 (425, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 0.96607));

    MassPoint mpisr129 (425, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 0.965932));

    MassPoint mpisr130 (425, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 0.965888));

    MassPoint mpisr131 (425, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 0.965665));

    MassPoint mpisr132 (425, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 0.965852));

    MassPoint mpisr133 (450, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 0.965124));

    MassPoint mpisr134 (450, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 0.965034));

    MassPoint mpisr135 (450, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 0.965785));

    MassPoint mpisr136 (450, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 0.965857));

    MassPoint mpisr137 (450, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 0.965423));

    MassPoint mpisr138 (450, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 0.96578));

    MassPoint mpisr139 (450, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 0.965651));

    MassPoint mpisr140 (450, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 0.965669));

    MassPoint mpisr141 (450, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 0.965455));

    MassPoint mpisr142 (450, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 0.964995));

    MassPoint mpisr143 (450, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 0.965442));

    MassPoint mpisr144 (450, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 0.965574));

    MassPoint mpisr145 (450, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 0.965733));

    MassPoint mpisr146 (450, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 0.96561));

    MassPoint mpisr147 (450, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 0.965278));

    MassPoint mpisr148 (450, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 0.965607));

    MassPoint mpisr149 (450, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 0.965313));

    MassPoint mpisr150 (475, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 0.965659));

    MassPoint mpisr151 (475, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 0.965653));

    MassPoint mpisr152 (475, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 0.965217));

    MassPoint mpisr153 (475, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 0.965224));

    MassPoint mpisr154 (475, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 0.965494));

    MassPoint mpisr155 (475, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 0.965574));

    MassPoint mpisr156 (475, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 0.964985));

    MassPoint mpisr157 (475, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 0.965599));

    MassPoint mpisr158 (475, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 0.965721));

    MassPoint mpisr159 (475, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 0.965445));

    MassPoint mpisr160 (475, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 0.96546));

    MassPoint mpisr161 (475, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 0.96575));

    MassPoint mpisr162 (475, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 0.96567));

    MassPoint mpisr163 (475, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 0.965686));

    MassPoint mpisr164 (475, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 0.965457));

    MassPoint mpisr165 (475, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 0.965506));

    MassPoint mpisr166 (475, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 0.965763));

    MassPoint mpisr167 (475, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 0.964972));

    MassPoint mpisr168 (500, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 0.964948));

    MassPoint mpisr169 (500, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 0.964783));

    MassPoint mpisr170 (500, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 0.965831));

    MassPoint mpisr171 (500, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 0.965024));

    MassPoint mpisr172 (500, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 0.965301));

    MassPoint mpisr173 (500, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 0.965168));

    MassPoint mpisr174 (500, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 0.966132));

    MassPoint mpisr175 (500, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 0.965145));

    MassPoint mpisr176 (500, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 0.965042));

    MassPoint mpisr177 (500, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 0.964716));

    MassPoint mpisr178 (500, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 0.965531));

    MassPoint mpisr179 (500, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 0.965182));

    MassPoint mpisr180 (500, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 0.965293));

    MassPoint mpisr181 (500, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 0.964981));

    MassPoint mpisr182 (500, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 0.96531));

    MassPoint mpisr183 (500, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 0.965257));

    MassPoint mpisr184 (500, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr184, 0.965123));

    MassPoint mpisr185 (500, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr185, 0.964918));

    MassPoint mpisr186 (500, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr186, 0.965527));

    MassPoint mpisr187 (525, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr187, 0.964922));

    MassPoint mpisr188 (525, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr188, 0.964928));

    MassPoint mpisr189 (525, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr189, 0.965274));

    MassPoint mpisr190 (525, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr190, 0.965074));

    MassPoint mpisr191 (525, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr191, 0.965122));

    MassPoint mpisr192 (525, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr192, 0.96509));

    MassPoint mpisr193 (525, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr193, 0.965308));

    MassPoint mpisr194 (525, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr194, 0.965149));

    MassPoint mpisr195 (525, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr195, 0.965235));

    MassPoint mpisr196 (525, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr196, 0.964805));

    MassPoint mpisr197 (525, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr197, 0.965248));

    MassPoint mpisr198 (525, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr198, 0.965477));

    MassPoint mpisr199 (525, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr199, 0.964826));

    MassPoint mpisr200 (525, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr200, 0.965222));

    MassPoint mpisr201 (525, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr201, 0.965147));

    MassPoint mpisr202 (525, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr202, 0.964869));

    MassPoint mpisr203 (525, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr203, 0.964919));

    MassPoint mpisr204 (525, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr204, 0.965457));

    MassPoint mpisr205 (525, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr205, 0.965312));

    MassPoint mpisr206 (525, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr206, 0.964839));

    MassPoint mpisr207 (550, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr207, 0.965208));

    MassPoint mpisr208 (550, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr208, 0.965097));

    MassPoint mpisr209 (550, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr209, 0.965227));

    MassPoint mpisr210 (550, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr210, 0.964997));

    MassPoint mpisr211 (550, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr211, 0.964968));

    MassPoint mpisr212 (550, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr212, 0.965001));

    MassPoint mpisr213 (550, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr213, 0.96504));

    MassPoint mpisr214 (550, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr214, 0.964894));

    MassPoint mpisr215 (550, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr215, 0.965021));

    MassPoint mpisr216 (550, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr216, 0.965095));

    MassPoint mpisr217 (550, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr217, 0.964938));

    MassPoint mpisr218 (550, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr218, 0.965106));

    MassPoint mpisr219 (550, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr219, 0.965344));

    MassPoint mpisr220 (550, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr220, 0.965022));

    MassPoint mpisr221 (550, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr221, 0.965071));

    MassPoint mpisr222 (550, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr222, 0.965061));

    MassPoint mpisr223 (550, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr223, 0.964456));

    MassPoint mpisr224 (550, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr224, 0.965089));

    MassPoint mpisr225 (550, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr225, 0.964748));

    MassPoint mpisr226 (550, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr226, 0.965034));

    MassPoint mpisr227 (550, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr227, 0.964716));

    MassPoint mpisr228 (575, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr228, 0.964988));

    MassPoint mpisr229 (575, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr229, 0.964761));

    MassPoint mpisr230 (575, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr230, 0.964973));

    MassPoint mpisr231 (575, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr231, 0.965011));

    MassPoint mpisr232 (575, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr232, 0.964705));

    MassPoint mpisr233 (575, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr233, 0.964968));

    MassPoint mpisr234 (575, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr234, 0.964993));

    MassPoint mpisr235 (575, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr235, 0.964824));

    MassPoint mpisr236 (575, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr236, 0.964591));

    MassPoint mpisr237 (575, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr237, 0.964654));

    MassPoint mpisr238 (575, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr238, 0.965003));

    MassPoint mpisr239 (575, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr239, 0.965242));

    MassPoint mpisr240 (575, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr240, 0.964495));

    MassPoint mpisr241 (575, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr241, 0.965039));

    MassPoint mpisr242 (575, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr242, 0.965245));

    MassPoint mpisr243 (575, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr243, 0.965096));

    MassPoint mpisr244 (575, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr244, 0.965515));

    MassPoint mpisr245 (575, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr245, 0.964991));

    MassPoint mpisr246 (575, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr246, 0.964293));

    MassPoint mpisr247 (575, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr247, 0.96422));

    MassPoint mpisr248 (575, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr248, 0.964796));

    MassPoint mpisr249 (575, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr249, 0.965348));

    MassPoint mpisr250 (600, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr250, 0.964441));

    MassPoint mpisr251 (600, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr251, 0.965345));

    MassPoint mpisr252 (600, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr252, 0.964989));

    MassPoint mpisr253 (600, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr253, 0.965339));

    MassPoint mpisr254 (600, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr254, 0.964377));

    MassPoint mpisr255 (600, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr255, 0.964319));

    MassPoint mpisr256 (600, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr256, 0.965322));

    MassPoint mpisr257 (600, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr257, 0.96528));

    MassPoint mpisr258 (600, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr258, 0.964956));

    MassPoint mpisr259 (600, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr259, 0.965391));

    MassPoint mpisr260 (600, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr260, 0.965519));

    MassPoint mpisr261 (600, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr261, 0.965275));

    MassPoint mpisr262 (600, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr262, 0.964568));

    MassPoint mpisr263 (600, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr263, 0.96512));

    MassPoint mpisr264 (600, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr264, 0.964936));

    MassPoint mpisr265 (600, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr265, 0.964805));

    MassPoint mpisr266 (600, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr266, 0.964709));

    MassPoint mpisr267 (600, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr267, 0.9649));

    MassPoint mpisr268 (600, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr268, 0.964549));

    MassPoint mpisr269 (600, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr269, 0.965087));

    MassPoint mpisr270 (600, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr270, 0.964494));

    MassPoint mpisr271 (600, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr271, 0.964574));

    MassPoint mpisr272 (600, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr272, 0.964701));

    MassPoint mpisr273 (625, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr273, 0.964636));

    MassPoint mpisr274 (625, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr274, 0.964714));

    MassPoint mpisr275 (625, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr275, 0.96499));

    MassPoint mpisr276 (625, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr276, 0.964554));

    MassPoint mpisr277 (625, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr277, 0.964262));

    MassPoint mpisr278 (625, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr278, 0.96444));

    MassPoint mpisr279 (625, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr279, 0.964461));

    MassPoint mpisr280 (625, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr280, 0.963986));

    MassPoint mpisr281 (625, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr281, 0.96485));

    MassPoint mpisr282 (625, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr282, 0.964591));

    MassPoint mpisr283 (625, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr283, 0.964522));

    MassPoint mpisr284 (625, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr284, 0.964121));

    MassPoint mpisr285 (625, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr285, 0.964559));

    MassPoint mpisr286 (625, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr286, 0.965675));

    MassPoint mpisr287 (625, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr287, 0.964294));

    MassPoint mpisr288 (625, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr288, 0.965385));

    MassPoint mpisr289 (625, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr289, 0.965163));

    MassPoint mpisr290 (625, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr290, 0.964557));

    MassPoint mpisr291 (625, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr291, 0.964939));

    MassPoint mpisr292 (625, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr292, 0.964737));

    MassPoint mpisr293 (625, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr293, 0.964734));

    MassPoint mpisr294 (625, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr294, 0.964672));

    MassPoint mpisr295 (625, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr295, 0.964878));

    MassPoint mpisr296 (625, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr296, 0.964864));

    MassPoint mpisr297 (650, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr297, 0.964757));

    MassPoint mpisr298 (650, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr298, 0.964961));

    MassPoint mpisr299 (650, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr299, 0.964444));

    MassPoint mpisr300 (650, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr300, 0.964555));

    MassPoint mpisr301 (650, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr301, 0.96473));

    MassPoint mpisr302 (650, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr302, 0.964208));

    MassPoint mpisr303 (650, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr303, 0.965532));

    MassPoint mpisr304 (650, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr304, 0.964927));

    MassPoint mpisr305 (650, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr305, 0.964705));

    MassPoint mpisr306 (650, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr306, 0.964535));

    MassPoint mpisr307 (650, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr307, 0.964714));

    MassPoint mpisr308 (650, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr308, 0.964371));

    MassPoint mpisr309 (650, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr309, 0.964959));

    MassPoint mpisr310 (650, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr310, 0.964313));

    MassPoint mpisr311 (650, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr311, 0.964554));

    MassPoint mpisr312 (650, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr312, 0.964661));

    MassPoint mpisr313 (650, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr313, 0.964696));

    MassPoint mpisr314 (650, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr314, 0.964713));

    MassPoint mpisr315 (650, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr315, 0.96442));

    MassPoint mpisr316 (650, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr316, 0.964668));

    MassPoint mpisr317 (650, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr317, 0.964927));

    MassPoint mpisr318 (650, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr318, 0.964607));

    MassPoint mpisr319 (650, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr319, 0.964556));

    MassPoint mpisr320 (650, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr320, 0.964665));

    MassPoint mpisr321 (650, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr321, 0.964466));

    MassPoint mpisr322 (675, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr322, 0.964728));

    MassPoint mpisr323 (675, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr323, 0.964321));

    MassPoint mpisr324 (675, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr324, 0.964057));

    MassPoint mpisr325 (675, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr325, 0.964098));

    MassPoint mpisr326 (675, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr326, 0.96457));

    MassPoint mpisr327 (675, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr327, 0.964511));

    MassPoint mpisr328 (675, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr328, 0.964766));

    MassPoint mpisr329 (675, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr329, 0.964394));

    MassPoint mpisr330 (675, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr330, 0.964795));

    MassPoint mpisr331 (675, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr331, 0.964441));

    MassPoint mpisr332 (675, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr332, 0.964387));

    MassPoint mpisr333 (675, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr333, 0.964767));

    MassPoint mpisr334 (675, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr334, 0.96427));

    MassPoint mpisr335 (675, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr335, 0.965319));

    MassPoint mpisr336 (675, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr336, 0.964827));

    MassPoint mpisr337 (675, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr337, 0.964891));

    MassPoint mpisr338 (675, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr338, 0.96523));

    MassPoint mpisr339 (675, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr339, 0.964899));

    MassPoint mpisr340 (675, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr340, 0.9648));

    MassPoint mpisr341 (675, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr341, 0.96408));

    MassPoint mpisr342 (675, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr342, 0.964322));

    MassPoint mpisr343 (675, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr343, 0.964621));

    MassPoint mpisr344 (675, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr344, 0.964884));

    MassPoint mpisr345 (675, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr345, 0.964801));

    MassPoint mpisr346 (675, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr346, 0.964839));

    MassPoint mpisr347 (700, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr347, 0.964625));

    MassPoint mpisr348 (700, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr348, 0.964269));

    MassPoint mpisr349 (700, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr349, 0.964682));

    MassPoint mpisr350 (700, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr350, 0.964421));

    MassPoint mpisr351 (700, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr351, 0.964646));

    MassPoint mpisr352 (700, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr352, 0.965158));

    MassPoint mpisr353 (700, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr353, 0.964861));

    MassPoint mpisr354 (700, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr354, 0.964884));

    MassPoint mpisr355 (700, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr355, 0.964323));

    MassPoint mpisr356 (700, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr356, 0.964602));

    MassPoint mpisr357 (700, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr357, 0.964591));

    MassPoint mpisr358 (700, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr358, 0.963935));

    MassPoint mpisr359 (700, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr359, 0.964367));

    MassPoint mpisr360 (700, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr360, 0.964426));

    MassPoint mpisr361 (700, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr361, 0.964791));

    MassPoint mpisr362 (700, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr362, 0.964878));

    MassPoint mpisr363 (700, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr363, 0.96459));

    MassPoint mpisr364 (700, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr364, 0.96459));

    MassPoint mpisr365 (700, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr365, 0.96435));

    MassPoint mpisr366 (700, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr366, 0.964525));

    MassPoint mpisr367 (700, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr367, 0.963882));

    MassPoint mpisr368 (700, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr368, 0.964665));

    MassPoint mpisr369 (700, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr369, 0.964329));

    MassPoint mpisr370 (700, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr370, 0.964872));

    MassPoint mpisr371 (700, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr371, 0.964974));

    MassPoint mpisr372 (725, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr372, 0.964265));

    MassPoint mpisr373 (725, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr373, 0.964692));

    MassPoint mpisr374 (725, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr374, 0.964861));

    MassPoint mpisr375 (725, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr375, 0.964442));

    MassPoint mpisr376 (725, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr376, 0.964371));

    MassPoint mpisr377 (725, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr377, 0.964634));

    MassPoint mpisr378 (725, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr378, 0.964275));

    MassPoint mpisr379 (725, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr379, 0.964636));

    MassPoint mpisr380 (725, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr380, 0.964142));

    MassPoint mpisr381 (725, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr381, 0.964366));

    MassPoint mpisr382 (725, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr382, 0.964723));

    MassPoint mpisr383 (725, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr383, 0.964329));

    MassPoint mpisr384 (725, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr384, 0.964974));

    MassPoint mpisr385 (725, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr385, 0.964726));

    MassPoint mpisr386 (725, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr386, 0.964529));

    MassPoint mpisr387 (725, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr387, 0.964612));

    MassPoint mpisr388 (725, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr388, 0.964815));

    MassPoint mpisr389 (725, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr389, 0.964411));

    MassPoint mpisr390 (725, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr390, 0.964926));

    MassPoint mpisr391 (725, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr391, 0.964624));

    MassPoint mpisr392 (725, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr392, 0.964818));

    MassPoint mpisr393 (725, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr393, 0.964734));

    MassPoint mpisr394 (725, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr394, 0.96438));

    MassPoint mpisr395 (725, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr395, 0.964613));

    MassPoint mpisr396 (725, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr396, 0.964597));

    MassPoint mpisr397 (750, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr397, 0.964778));

    MassPoint mpisr398 (750, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr398, 0.964296));

    MassPoint mpisr399 (750, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr399, 0.964373));

    MassPoint mpisr400 (750, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr400, 0.964574));

    MassPoint mpisr401 (750, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr401, 0.965004));

    MassPoint mpisr402 (750, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr402, 0.964968));

    MassPoint mpisr403 (750, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr403, 0.964209));

    MassPoint mpisr404 (750, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr404, 0.9646));

    MassPoint mpisr405 (750, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr405, 0.964802));

    MassPoint mpisr406 (750, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr406, 0.964423));

    MassPoint mpisr407 (750, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr407, 0.964794));

    MassPoint mpisr408 (750, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr408, 0.964728));

    MassPoint mpisr409 (750, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr409, 0.964093));

    MassPoint mpisr410 (750, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr410, 0.96493));

    MassPoint mpisr411 (750, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr411, 0.964512));

    MassPoint mpisr412 (750, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr412, 0.964627));

    MassPoint mpisr413 (750, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr413, 0.96479));

    MassPoint mpisr414 (750, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr414, 0.964886));

    MassPoint mpisr415 (750, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr415, 0.964173));

    MassPoint mpisr416 (750, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr416, 0.964544));

    MassPoint mpisr417 (750, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr417, 0.96435));

    MassPoint mpisr418 (750, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr418, 0.964751));

    MassPoint mpisr419 (750, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr419, 0.964493));

    MassPoint mpisr420 (750, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr420, 0.964721));

    MassPoint mpisr421 (750, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr421, 0.964076));

    MassPoint mpisr422 (775, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr422, 0.964999));

    MassPoint mpisr423 (775, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr423, 0.965114));

    MassPoint mpisr424 (775, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr424, 0.964482));

    MassPoint mpisr425 (775, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr425, 0.964651));

    MassPoint mpisr426 (775, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr426, 0.964172));

    MassPoint mpisr427 (775, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr427, 0.964685));

    MassPoint mpisr428 (775, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr428, 0.964766));

    MassPoint mpisr429 (775, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr429, 0.964474));

    MassPoint mpisr430 (775, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr430, 0.964852));

    MassPoint mpisr431 (775, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr431, 0.96443));

    MassPoint mpisr432 (775, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr432, 0.964672));

    MassPoint mpisr433 (775, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr433, 0.964417));

    MassPoint mpisr434 (775, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr434, 0.964258));

    MassPoint mpisr435 (775, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr435, 0.964752));

    MassPoint mpisr436 (775, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr436, 0.964186));

    MassPoint mpisr437 (775, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr437, 0.964416));

    MassPoint mpisr438 (775, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr438, 0.964365));

    MassPoint mpisr439 (775, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr439, 0.964396));

    MassPoint mpisr440 (775, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr440, 0.964473));

    MassPoint mpisr441 (775, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr441, 0.964295));

    MassPoint mpisr442 (775, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr442, 0.96428));

    MassPoint mpisr443 (775, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr443, 0.964453));

    MassPoint mpisr444 (775, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr444, 0.964037));

    MassPoint mpisr445 (775, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr445, 0.964883));

    MassPoint mpisr446 (775, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr446, 0.96431));

    MassPoint mpisr447 (800, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr447, 0.964929));

    MassPoint mpisr448 (800, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr448, 0.964337));

    MassPoint mpisr449 (800, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr449, 0.964475));

    MassPoint mpisr450 (800, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr450, 0.964118));

    MassPoint mpisr451 (800, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr451, 0.964236));

    MassPoint mpisr452 (800, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr452, 0.964344));

    MassPoint mpisr453 (800, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr453, 0.964407));

    MassPoint mpisr454 (800, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr454, 0.965062));

    MassPoint mpisr455 (800, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr455, 0.964501));

    MassPoint mpisr456 (800, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr456, 0.964749));

    MassPoint mpisr457 (800, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr457, 0.964539));

    MassPoint mpisr458 (800, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr458, 0.964277));

    MassPoint mpisr459 (800, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr459, 0.964511));

    MassPoint mpisr460 (800, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr460, 0.964709));

    MassPoint mpisr461 (800, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr461, 0.964341));

    MassPoint mpisr462 (800, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr462, 0.964656));

    MassPoint mpisr463 (800, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr463, 0.964026));

    MassPoint mpisr464 (800, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr464, 0.964551));

    MassPoint mpisr465 (800, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr465, 0.964847));

    MassPoint mpisr466 (800, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr466, 0.964131));

    MassPoint mpisr467 (800, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr467, 0.9647));

    MassPoint mpisr468 (800, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr468, 0.964169));

    MassPoint mpisr469 (800, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr469, 0.964785));

    MassPoint mpisr470 (800, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr470, 0.964601));

    MassPoint mpisr471 (800, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr471, 0.964478)); 

  } else  if (SUSYProductionProcess=="TChiSlep") {

    MassPoint mp0 (100, 1);
    StopCrossSection cs0 (11.6119, 0.518613);
    MassPointParameters mpp0 (cs0, 8625);
    StopNeutralinoMap.insert(std::make_pair(mp0, mpp0));

    MassPoint mp1 (100, 25);
    StopCrossSection cs1 (11.6119, 0.518613);
    MassPointParameters mpp1 (cs1, 11890);
    StopNeutralinoMap.insert(std::make_pair(mp1, mpp1));

    MassPoint mp2 (100, 50);
    StopCrossSection cs2 (11.6119, 0.518613);
    MassPointParameters mpp2 (cs2, 12845);
    StopNeutralinoMap.insert(std::make_pair(mp2, mpp2));

    MassPoint mp3 (125, 1);
    StopCrossSection cs3 (5.09052, 0.249469);
    MassPointParameters mpp3 (cs3, 9201);
    StopNeutralinoMap.insert(std::make_pair(mp3, mpp3));

    MassPoint mp4 (125, 25);
    StopCrossSection cs4 (5.09052, 0.249469);
    MassPointParameters mpp4 (cs4, 8450);
    StopNeutralinoMap.insert(std::make_pair(mp4, mpp4));

    MassPoint mp5 (125, 50);
    StopCrossSection cs5 (5.09052, 0.249469);
    MassPointParameters mpp5 (cs5, 12360);
    StopNeutralinoMap.insert(std::make_pair(mp5, mpp5));

    MassPoint mp6 (125, 75);
    StopCrossSection cs6 (5.09052, 0.249469);
    MassPointParameters mpp6 (cs6, 9903);
    StopNeutralinoMap.insert(std::make_pair(mp6, mpp6));

    MassPoint mp7 (150, 1);
    StopCrossSection cs7 (2.61231, 0.138156);
    MassPointParameters mpp7 (cs7, 11132);
    StopNeutralinoMap.insert(std::make_pair(mp7, mpp7));

    MassPoint mp8 (150, 25);
    StopCrossSection cs8 (2.61231, 0.138156);
    MassPointParameters mpp8 (cs8, 9169);
    StopNeutralinoMap.insert(std::make_pair(mp8, mpp8));

    MassPoint mp9 (150, 50);
    StopCrossSection cs9 (2.61231, 0.138156);
    MassPointParameters mpp9 (cs9, 10942);
    StopNeutralinoMap.insert(std::make_pair(mp9, mpp9));

    MassPoint mp10 (150, 75);
    StopCrossSection cs10 (2.61231, 0.138156);
    MassPointParameters mpp10 (cs10, 12012);
    StopNeutralinoMap.insert(std::make_pair(mp10, mpp10));

    MassPoint mp11 (150, 100);
    StopCrossSection cs11 (2.61231, 0.138156);
    MassPointParameters mpp11 (cs11, 10779);
    StopNeutralinoMap.insert(std::make_pair(mp11, mpp11));

    MassPoint mp12 (175, 1);
    StopCrossSection cs12 (1.48242, 0.0832672);
    MassPointParameters mpp12 (cs12, 9601);
    StopNeutralinoMap.insert(std::make_pair(mp12, mpp12));

    MassPoint mp13 (175, 25);
    StopCrossSection cs13 (1.48242, 0.0832672);
    MassPointParameters mpp13 (cs13, 11529);
    StopNeutralinoMap.insert(std::make_pair(mp13, mpp13));

    MassPoint mp14 (175, 50);
    StopCrossSection cs14 (1.48242, 0.0832672);
    MassPointParameters mpp14 (cs14, 7725);
    StopNeutralinoMap.insert(std::make_pair(mp14, mpp14));

    MassPoint mp15 (175, 75);
    StopCrossSection cs15 (1.48242, 0.0832672);
    MassPointParameters mpp15 (cs15, 8560);
    StopNeutralinoMap.insert(std::make_pair(mp15, mpp15));

    MassPoint mp16 (175, 100);
    StopCrossSection cs16 (1.48242, 0.0832672);
    MassPointParameters mpp16 (cs16, 6674);
    StopNeutralinoMap.insert(std::make_pair(mp16, mpp16));

    MassPoint mp17 (175, 125);
    StopCrossSection cs17 (1.48242, 0.0832672);
    MassPointParameters mpp17 (cs17, 11973);
    StopNeutralinoMap.insert(std::make_pair(mp17, mpp17));

    MassPoint mp18 (200, 1);
    StopCrossSection cs18 (0.902569, 0.0537411);
    MassPointParameters mpp18 (cs18, 10065);
    StopNeutralinoMap.insert(std::make_pair(mp18, mpp18));

    MassPoint mp19 (200, 25);
    StopCrossSection cs19 (0.902569, 0.0537411);
    MassPointParameters mpp19 (cs19, 8416);
    StopNeutralinoMap.insert(std::make_pair(mp19, mpp19));

    MassPoint mp20 (200, 50);
    StopCrossSection cs20 (0.902569, 0.0537411);
    MassPointParameters mpp20 (cs20, 10885);
    StopNeutralinoMap.insert(std::make_pair(mp20, mpp20));

    MassPoint mp21 (200, 75);
    StopCrossSection cs21 (0.902569, 0.0537411);
    MassPointParameters mpp21 (cs21, 7196);
    StopNeutralinoMap.insert(std::make_pair(mp21, mpp21));

    MassPoint mp22 (200, 100);
    StopCrossSection cs22 (0.902569, 0.0537411);
    MassPointParameters mpp22 (cs22, 11505);
    StopNeutralinoMap.insert(std::make_pair(mp22, mpp22));

    MassPoint mp23 (200, 125);
    StopCrossSection cs23 (0.902569, 0.0537411);
    MassPointParameters mpp23 (cs23, 11256);
    StopNeutralinoMap.insert(std::make_pair(mp23, mpp23));

    MassPoint mp24 (200, 150);
    StopCrossSection cs24 (0.902569, 0.0537411);
    MassPointParameters mpp24 (cs24, 9382);
    StopNeutralinoMap.insert(std::make_pair(mp24, mpp24));

    MassPoint mp25 (225, 1);
    StopCrossSection cs25 (0.579564, 0.0360699);
    MassPointParameters mpp25 (cs25, 11009);
    StopNeutralinoMap.insert(std::make_pair(mp25, mpp25));

    MassPoint mp26 (225, 25);
    StopCrossSection cs26 (0.579564, 0.0360699);
    MassPointParameters mpp26 (cs26, 7996);
    StopNeutralinoMap.insert(std::make_pair(mp26, mpp26));

    MassPoint mp27 (225, 50);
    StopCrossSection cs27 (0.579564, 0.0360699);
    MassPointParameters mpp27 (cs27, 10446);
    StopNeutralinoMap.insert(std::make_pair(mp27, mpp27));

    MassPoint mp28 (225, 75);
    StopCrossSection cs28 (0.579564, 0.0360699);
    MassPointParameters mpp28 (cs28, 9081);
    StopNeutralinoMap.insert(std::make_pair(mp28, mpp28));

    MassPoint mp29 (225, 100);
    StopCrossSection cs29 (0.579564, 0.0360699);
    MassPointParameters mpp29 (cs29, 13155);
    StopNeutralinoMap.insert(std::make_pair(mp29, mpp29));

    MassPoint mp30 (225, 125);
    StopCrossSection cs30 (0.579564, 0.0360699);
    MassPointParameters mpp30 (cs30, 11469);
    StopNeutralinoMap.insert(std::make_pair(mp30, mpp30));

    MassPoint mp31 (225, 150);
    StopCrossSection cs31 (0.579564, 0.0360699);
    MassPointParameters mpp31 (cs31, 9419);
    StopNeutralinoMap.insert(std::make_pair(mp31, mpp31));

    MassPoint mp32 (225, 175);
    StopCrossSection cs32 (0.579564, 0.0360699);
    MassPointParameters mpp32 (cs32, 10153);
    StopNeutralinoMap.insert(std::make_pair(mp32, mpp32));

    MassPoint mp33 (250, 1);
    StopCrossSection cs33 (0.387534, 0.0253131);
    MassPointParameters mpp33 (cs33, 11869);
    StopNeutralinoMap.insert(std::make_pair(mp33, mpp33));

    MassPoint mp34 (250, 25);
    StopCrossSection cs34 (0.387534, 0.0253131);
    MassPointParameters mpp34 (cs34, 8842);
    StopNeutralinoMap.insert(std::make_pair(mp34, mpp34));

    MassPoint mp35 (250, 50);
    StopCrossSection cs35 (0.387534, 0.0253131);
    MassPointParameters mpp35 (cs35, 9855);
    StopNeutralinoMap.insert(std::make_pair(mp35, mpp35));

    MassPoint mp36 (250, 75);
    StopCrossSection cs36 (0.387534, 0.0253131);
    MassPointParameters mpp36 (cs36, 10466);
    StopNeutralinoMap.insert(std::make_pair(mp36, mpp36));

    MassPoint mp37 (250, 100);
    StopCrossSection cs37 (0.387534, 0.0253131);
    MassPointParameters mpp37 (cs37, 9187);
    StopNeutralinoMap.insert(std::make_pair(mp37, mpp37));

    MassPoint mp38 (250, 125);
    StopCrossSection cs38 (0.387534, 0.0253131);
    MassPointParameters mpp38 (cs38, 11896);
    StopNeutralinoMap.insert(std::make_pair(mp38, mpp38));

    MassPoint mp39 (250, 150);
    StopCrossSection cs39 (0.387534, 0.0253131);
    MassPointParameters mpp39 (cs39, 7147);
    StopNeutralinoMap.insert(std::make_pair(mp39, mpp39));

    MassPoint mp40 (250, 175);
    StopCrossSection cs40 (0.387534, 0.0253131);
    MassPointParameters mpp40 (cs40, 7674);
    StopNeutralinoMap.insert(std::make_pair(mp40, mpp40));

    MassPoint mp41 (250, 200);
    StopCrossSection cs41 (0.387534, 0.0253131);
    MassPointParameters mpp41 (cs41, 9130);
    StopNeutralinoMap.insert(std::make_pair(mp41, mpp41));

    MassPoint mp42 (275, 1);
    StopCrossSection cs42 (0.267786, 0.0182886);
    MassPointParameters mpp42 (cs42, 7574);
    StopNeutralinoMap.insert(std::make_pair(mp42, mpp42));

    MassPoint mp43 (275, 25);
    StopCrossSection cs43 (0.267786, 0.0182886);
    MassPointParameters mpp43 (cs43, 9284);
    StopNeutralinoMap.insert(std::make_pair(mp43, mpp43));

    MassPoint mp44 (275, 50);
    StopCrossSection cs44 (0.267786, 0.0182886);
    MassPointParameters mpp44 (cs44, 11034);
    StopNeutralinoMap.insert(std::make_pair(mp44, mpp44));

    MassPoint mp45 (275, 75);
    StopCrossSection cs45 (0.267786, 0.0182886);
    MassPointParameters mpp45 (cs45, 8287);
    StopNeutralinoMap.insert(std::make_pair(mp45, mpp45));

    MassPoint mp46 (275, 100);
    StopCrossSection cs46 (0.267786, 0.0182886);
    MassPointParameters mpp46 (cs46, 10626);
    StopNeutralinoMap.insert(std::make_pair(mp46, mpp46));

    MassPoint mp47 (275, 125);
    StopCrossSection cs47 (0.267786, 0.0182886);
    MassPointParameters mpp47 (cs47, 9369);
    StopNeutralinoMap.insert(std::make_pair(mp47, mpp47));

    MassPoint mp48 (275, 150);
    StopCrossSection cs48 (0.267786, 0.0182886);
    MassPointParameters mpp48 (cs48, 10156);
    StopNeutralinoMap.insert(std::make_pair(mp48, mpp48));

    MassPoint mp49 (275, 175);
    StopCrossSection cs49 (0.267786, 0.0182886);
    MassPointParameters mpp49 (cs49, 9102);
    StopNeutralinoMap.insert(std::make_pair(mp49, mpp49));

    MassPoint mp50 (275, 200);
    StopCrossSection cs50 (0.267786, 0.0182886);
    MassPointParameters mpp50 (cs50, 11781);
    StopNeutralinoMap.insert(std::make_pair(mp50, mpp50));

    MassPoint mp51 (275, 225);
    StopCrossSection cs51 (0.267786, 0.0182886);
    MassPointParameters mpp51 (cs51, 9505);
    StopNeutralinoMap.insert(std::make_pair(mp51, mpp51));

    MassPoint mp52 (300, 1);
    StopCrossSection cs52 (0.190159, 0.0134438);
    MassPointParameters mpp52 (cs52, 10052);
    StopNeutralinoMap.insert(std::make_pair(mp52, mpp52));

    MassPoint mp53 (300, 25);
    StopCrossSection cs53 (0.190159, 0.0134438);
    MassPointParameters mpp53 (cs53, 10917);
    StopNeutralinoMap.insert(std::make_pair(mp53, mpp53));

    MassPoint mp54 (300, 50);
    StopCrossSection cs54 (0.190159, 0.0134438);
    MassPointParameters mpp54 (cs54, 8040);
    StopNeutralinoMap.insert(std::make_pair(mp54, mpp54));

    MassPoint mp55 (300, 75);
    StopCrossSection cs55 (0.190159, 0.0134438);
    MassPointParameters mpp55 (cs55, 7470);
    StopNeutralinoMap.insert(std::make_pair(mp55, mpp55));

    MassPoint mp56 (300, 100);
    StopCrossSection cs56 (0.190159, 0.0134438);
    MassPointParameters mpp56 (cs56, 11334);
    StopNeutralinoMap.insert(std::make_pair(mp56, mpp56));

    MassPoint mp57 (300, 125);
    StopCrossSection cs57 (0.190159, 0.0134438);
    MassPointParameters mpp57 (cs57, 9540);
    StopNeutralinoMap.insert(std::make_pair(mp57, mpp57));

    MassPoint mp58 (300, 150);
    StopCrossSection cs58 (0.190159, 0.0134438);
    MassPointParameters mpp58 (cs58, 11137);
    StopNeutralinoMap.insert(std::make_pair(mp58, mpp58));

    MassPoint mp59 (300, 175);
    StopCrossSection cs59 (0.190159, 0.0134438);
    MassPointParameters mpp59 (cs59, 11029);
    StopNeutralinoMap.insert(std::make_pair(mp59, mpp59));

    MassPoint mp60 (300, 200);
    StopCrossSection cs60 (0.190159, 0.0134438);
    MassPointParameters mpp60 (cs60, 11164);
    StopNeutralinoMap.insert(std::make_pair(mp60, mpp60));

    MassPoint mp61 (300, 225);
    StopCrossSection cs61 (0.190159, 0.0134438);
    MassPointParameters mpp61 (cs61, 9929);
    StopNeutralinoMap.insert(std::make_pair(mp61, mpp61));

    MassPoint mp62 (300, 250);
    StopCrossSection cs62 (0.190159, 0.0134438);
    MassPointParameters mpp62 (cs62, 10078);
    StopNeutralinoMap.insert(std::make_pair(mp62, mpp62));

    MassPoint mp63 (325, 1);
    StopCrossSection cs63 (0.138086, 0.0101835);
    MassPointParameters mpp63 (cs63, 10713);
    StopNeutralinoMap.insert(std::make_pair(mp63, mpp63));

    MassPoint mp64 (325, 25);
    StopCrossSection cs64 (0.138086, 0.0101835);
    MassPointParameters mpp64 (cs64, 8988);
    StopNeutralinoMap.insert(std::make_pair(mp64, mpp64));

    MassPoint mp65 (325, 50);
    StopCrossSection cs65 (0.138086, 0.0101835);
    MassPointParameters mpp65 (cs65, 7569);
    StopNeutralinoMap.insert(std::make_pair(mp65, mpp65));

    MassPoint mp66 (325, 75);
    StopCrossSection cs66 (0.138086, 0.0101835);
    MassPointParameters mpp66 (cs66, 6793);
    StopNeutralinoMap.insert(std::make_pair(mp66, mpp66));

    MassPoint mp67 (325, 100);
    StopCrossSection cs67 (0.138086, 0.0101835);
    MassPointParameters mpp67 (cs67, 7608);
    StopNeutralinoMap.insert(std::make_pair(mp67, mpp67));

    MassPoint mp68 (325, 125);
    StopCrossSection cs68 (0.138086, 0.0101835);
    MassPointParameters mpp68 (cs68, 11778);
    StopNeutralinoMap.insert(std::make_pair(mp68, mpp68));

    MassPoint mp69 (325, 150);
    StopCrossSection cs69 (0.138086, 0.0101835);
    MassPointParameters mpp69 (cs69, 9547);
    StopNeutralinoMap.insert(std::make_pair(mp69, mpp69));

    MassPoint mp70 (325, 175);
    StopCrossSection cs70 (0.138086, 0.0101835);
    MassPointParameters mpp70 (cs70, 11230);
    StopNeutralinoMap.insert(std::make_pair(mp70, mpp70));

    MassPoint mp71 (325, 200);
    StopCrossSection cs71 (0.138086, 0.0101835);
    MassPointParameters mpp71 (cs71, 10468);
    StopNeutralinoMap.insert(std::make_pair(mp71, mpp71));

    MassPoint mp72 (325, 225);
    StopCrossSection cs72 (0.138086, 0.0101835);
    MassPointParameters mpp72 (cs72, 8122);
    StopNeutralinoMap.insert(std::make_pair(mp72, mpp72));

    MassPoint mp73 (325, 250);
    StopCrossSection cs73 (0.138086, 0.0101835);
    MassPointParameters mpp73 (cs73, 8664);
    StopNeutralinoMap.insert(std::make_pair(mp73, mpp73));

    MassPoint mp74 (325, 275);
    StopCrossSection cs74 (0.138086, 0.0101835);
    MassPointParameters mpp74 (cs74, 9403);
    StopNeutralinoMap.insert(std::make_pair(mp74, mpp74));

    MassPoint mp75 (350, 1);
    StopCrossSection cs75 (0.102199, 0.00775261);
    MassPointParameters mpp75 (cs75, 10723);
    StopNeutralinoMap.insert(std::make_pair(mp75, mpp75));

    MassPoint mp76 (350, 25);
    StopCrossSection cs76 (0.102199, 0.00775261);
    MassPointParameters mpp76 (cs76, 6873);
    StopNeutralinoMap.insert(std::make_pair(mp76, mpp76));

    MassPoint mp77 (350, 50);
    StopCrossSection cs77 (0.102199, 0.00775261);
    MassPointParameters mpp77 (cs77, 7917);
    StopNeutralinoMap.insert(std::make_pair(mp77, mpp77));

    MassPoint mp78 (350, 75);
    StopCrossSection cs78 (0.102199, 0.00775261);
    MassPointParameters mpp78 (cs78, 9548);
    StopNeutralinoMap.insert(std::make_pair(mp78, mpp78));

    MassPoint mp79 (350, 100);
    StopCrossSection cs79 (0.102199, 0.00775261);
    MassPointParameters mpp79 (cs79, 9267);
    StopNeutralinoMap.insert(std::make_pair(mp79, mpp79));

    MassPoint mp80 (350, 125);
    StopCrossSection cs80 (0.102199, 0.00775261);
    MassPointParameters mpp80 (cs80, 10030);
    StopNeutralinoMap.insert(std::make_pair(mp80, mpp80));

    MassPoint mp81 (350, 150);
    StopCrossSection cs81 (0.102199, 0.00775261);
    MassPointParameters mpp81 (cs81, 12319);
    StopNeutralinoMap.insert(std::make_pair(mp81, mpp81));

    MassPoint mp82 (350, 175);
    StopCrossSection cs82 (0.102199, 0.00775261);
    MassPointParameters mpp82 (cs82, 12558);
    StopNeutralinoMap.insert(std::make_pair(mp82, mpp82));

    MassPoint mp83 (350, 200);
    StopCrossSection cs83 (0.102199, 0.00775261);
    MassPointParameters mpp83 (cs83, 9694);
    StopNeutralinoMap.insert(std::make_pair(mp83, mpp83));

    MassPoint mp84 (350, 225);
    StopCrossSection cs84 (0.102199, 0.00775261);
    MassPointParameters mpp84 (cs84, 6913);
    StopNeutralinoMap.insert(std::make_pair(mp84, mpp84));

    MassPoint mp85 (350, 250);
    StopCrossSection cs85 (0.102199, 0.00775261);
    MassPointParameters mpp85 (cs85, 10531);
    StopNeutralinoMap.insert(std::make_pair(mp85, mpp85));

    MassPoint mp86 (350, 275);
    StopCrossSection cs86 (0.102199, 0.00775261);
    MassPointParameters mpp86 (cs86, 9286);
    StopNeutralinoMap.insert(std::make_pair(mp86, mpp86));

    MassPoint mp87 (350, 300);
    StopCrossSection cs87 (0.102199, 0.00775261);
    MassPointParameters mpp87 (cs87, 9862);
    StopNeutralinoMap.insert(std::make_pair(mp87, mpp87));

    MassPoint mp88 (375, 1);
    StopCrossSection cs88 (0.0768342, 0.00602606);
    MassPointParameters mpp88 (cs88, 10185);
    StopNeutralinoMap.insert(std::make_pair(mp88, mpp88));

    MassPoint mp89 (375, 25);
    StopCrossSection cs89 (0.0768342, 0.00602606);
    MassPointParameters mpp89 (cs89, 8200);
    StopNeutralinoMap.insert(std::make_pair(mp89, mpp89));

    MassPoint mp90 (375, 50);
    StopCrossSection cs90 (0.0768342, 0.00602606);
    MassPointParameters mpp90 (cs90, 10921);
    StopNeutralinoMap.insert(std::make_pair(mp90, mpp90));

    MassPoint mp91 (375, 75);
    StopCrossSection cs91 (0.0768342, 0.00602606);
    MassPointParameters mpp91 (cs91, 9279);
    StopNeutralinoMap.insert(std::make_pair(mp91, mpp91));

    MassPoint mp92 (375, 100);
    StopCrossSection cs92 (0.0768342, 0.00602606);
    MassPointParameters mpp92 (cs92, 6808);
    StopNeutralinoMap.insert(std::make_pair(mp92, mpp92));

    MassPoint mp93 (375, 125);
    StopCrossSection cs93 (0.0768342, 0.00602606);
    MassPointParameters mpp93 (cs93, 8529);
    StopNeutralinoMap.insert(std::make_pair(mp93, mpp93));

    MassPoint mp94 (375, 150);
    StopCrossSection cs94 (0.0768342, 0.00602606);
    MassPointParameters mpp94 (cs94, 11022);
    StopNeutralinoMap.insert(std::make_pair(mp94, mpp94));

    MassPoint mp95 (375, 175);
    StopCrossSection cs95 (0.0768342, 0.00602606);
    MassPointParameters mpp95 (cs95, 8076);
    StopNeutralinoMap.insert(std::make_pair(mp95, mpp95));

    MassPoint mp96 (375, 200);
    StopCrossSection cs96 (0.0768342, 0.00602606);
    MassPointParameters mpp96 (cs96, 9394);
    StopNeutralinoMap.insert(std::make_pair(mp96, mpp96));

    MassPoint mp97 (375, 225);
    StopCrossSection cs97 (0.0768342, 0.00602606);
    MassPointParameters mpp97 (cs97, 9647);
    StopNeutralinoMap.insert(std::make_pair(mp97, mpp97));

    MassPoint mp98 (375, 250);
    StopCrossSection cs98 (0.0768342, 0.00602606);
    MassPointParameters mpp98 (cs98, 8084);
    StopNeutralinoMap.insert(std::make_pair(mp98, mpp98));

    MassPoint mp99 (375, 275);
    StopCrossSection cs99 (0.0768342, 0.00602606);
    MassPointParameters mpp99 (cs99, 10570);
    StopNeutralinoMap.insert(std::make_pair(mp99, mpp99));

    MassPoint mp100 (375, 300);
    StopCrossSection cs100 (0.0768342, 0.00602606);
    MassPointParameters mpp100 (cs100, 8454);
    StopNeutralinoMap.insert(std::make_pair(mp100, mpp100));

    MassPoint mp101 (375, 325);
    StopCrossSection cs101 (0.0768342, 0.00602606);
    MassPointParameters mpp101 (cs101, 11182);
    StopNeutralinoMap.insert(std::make_pair(mp101, mpp101));

    MassPoint mp102 (400, 1);
    StopCrossSection cs102 (0.0586311, 0.0047276);
    MassPointParameters mpp102 (cs102, 8180);
    StopNeutralinoMap.insert(std::make_pair(mp102, mpp102));

    MassPoint mp103 (400, 25);
    StopCrossSection cs103 (0.0586311, 0.0047276);
    MassPointParameters mpp103 (cs103, 9144);
    StopNeutralinoMap.insert(std::make_pair(mp103, mpp103));

    MassPoint mp104 (400, 50);
    StopCrossSection cs104 (0.0586311, 0.0047276);
    MassPointParameters mpp104 (cs104, 14234);
    StopNeutralinoMap.insert(std::make_pair(mp104, mpp104));

    MassPoint mp105 (400, 75);
    StopCrossSection cs105 (0.0586311, 0.0047276);
    MassPointParameters mpp105 (cs105, 10180);
    StopNeutralinoMap.insert(std::make_pair(mp105, mpp105));

    MassPoint mp106 (400, 100);
    StopCrossSection cs106 (0.0586311, 0.0047276);
    MassPointParameters mpp106 (cs106, 8656);
    StopNeutralinoMap.insert(std::make_pair(mp106, mpp106));

    MassPoint mp107 (400, 125);
    StopCrossSection cs107 (0.0586311, 0.0047276);
    MassPointParameters mpp107 (cs107, 11250);
    StopNeutralinoMap.insert(std::make_pair(mp107, mpp107));

    MassPoint mp108 (400, 150);
    StopCrossSection cs108 (0.0586311, 0.0047276);
    MassPointParameters mpp108 (cs108, 8933);
    StopNeutralinoMap.insert(std::make_pair(mp108, mpp108));

    MassPoint mp109 (400, 175);
    StopCrossSection cs109 (0.0586311, 0.0047276);
    MassPointParameters mpp109 (cs109, 12403);
    StopNeutralinoMap.insert(std::make_pair(mp109, mpp109));

    MassPoint mp110 (400, 200);
    StopCrossSection cs110 (0.0586311, 0.0047276);
    MassPointParameters mpp110 (cs110, 7706);
    StopNeutralinoMap.insert(std::make_pair(mp110, mpp110));

    MassPoint mp111 (400, 225);
    StopCrossSection cs111 (0.0586311, 0.0047276);
    MassPointParameters mpp111 (cs111, 8795);
    StopNeutralinoMap.insert(std::make_pair(mp111, mpp111));

    MassPoint mp112 (400, 250);
    StopCrossSection cs112 (0.0586311, 0.0047276);
    MassPointParameters mpp112 (cs112, 11721);
    StopNeutralinoMap.insert(std::make_pair(mp112, mpp112));

    MassPoint mp113 (400, 275);
    StopCrossSection cs113 (0.0586311, 0.0047276);
    MassPointParameters mpp113 (cs113, 11593);
    StopNeutralinoMap.insert(std::make_pair(mp113, mpp113));

    MassPoint mp114 (400, 300);
    StopCrossSection cs114 (0.0586311, 0.0047276);
    MassPointParameters mpp114 (cs114, 10036);
    StopNeutralinoMap.insert(std::make_pair(mp114, mpp114));

    MassPoint mp115 (400, 325);
    StopCrossSection cs115 (0.0586311, 0.0047276);
    MassPointParameters mpp115 (cs115, 10870);
    StopNeutralinoMap.insert(std::make_pair(mp115, mpp115));

    MassPoint mp116 (400, 350);
    StopCrossSection cs116 (0.0586311, 0.0047276);
    MassPointParameters mpp116 (cs116, 10291);
    StopNeutralinoMap.insert(std::make_pair(mp116, mpp116));

    MassPoint mp117 (425, 1);
    StopCrossSection cs117 (0.0452189, 0.00371547);
    MassPointParameters mpp117 (cs117, 10002);
    StopNeutralinoMap.insert(std::make_pair(mp117, mpp117));

    MassPoint mp118 (425, 25);
    StopCrossSection cs118 (0.0452189, 0.00371547);
    MassPointParameters mpp118 (cs118, 12401);
    StopNeutralinoMap.insert(std::make_pair(mp118, mpp118));

    MassPoint mp119 (425, 50);
    StopCrossSection cs119 (0.0452189, 0.00371547);
    MassPointParameters mpp119 (cs119, 8267);
    StopNeutralinoMap.insert(std::make_pair(mp119, mpp119));

    MassPoint mp120 (425, 75);
    StopCrossSection cs120 (0.0452189, 0.00371547);
    MassPointParameters mpp120 (cs120, 11919);
    StopNeutralinoMap.insert(std::make_pair(mp120, mpp120));

    MassPoint mp121 (425, 100);
    StopCrossSection cs121 (0.0452189, 0.00371547);
    MassPointParameters mpp121 (cs121, 9785);
    StopNeutralinoMap.insert(std::make_pair(mp121, mpp121));

    MassPoint mp122 (425, 125);
    StopCrossSection cs122 (0.0452189, 0.00371547);
    MassPointParameters mpp122 (cs122, 10192);
    StopNeutralinoMap.insert(std::make_pair(mp122, mpp122));

    MassPoint mp123 (425, 150);
    StopCrossSection cs123 (0.0452189, 0.00371547);
    MassPointParameters mpp123 (cs123, 11793);
    StopNeutralinoMap.insert(std::make_pair(mp123, mpp123));

    MassPoint mp124 (425, 175);
    StopCrossSection cs124 (0.0452189, 0.00371547);
    MassPointParameters mpp124 (cs124, 8942);
    StopNeutralinoMap.insert(std::make_pair(mp124, mpp124));

    MassPoint mp125 (425, 200);
    StopCrossSection cs125 (0.0452189, 0.00371547);
    MassPointParameters mpp125 (cs125, 9816);
    StopNeutralinoMap.insert(std::make_pair(mp125, mpp125));

    MassPoint mp126 (425, 225);
    StopCrossSection cs126 (0.0452189, 0.00371547);
    MassPointParameters mpp126 (cs126, 8097);
    StopNeutralinoMap.insert(std::make_pair(mp126, mpp126));

    MassPoint mp127 (425, 250);
    StopCrossSection cs127 (0.0452189, 0.00371547);
    MassPointParameters mpp127 (cs127, 10061);
    StopNeutralinoMap.insert(std::make_pair(mp127, mpp127));

    MassPoint mp128 (425, 275);
    StopCrossSection cs128 (0.0452189, 0.00371547);
    MassPointParameters mpp128 (cs128, 7061);
    StopNeutralinoMap.insert(std::make_pair(mp128, mpp128));

    MassPoint mp129 (425, 300);
    StopCrossSection cs129 (0.0452189, 0.00371547);
    MassPointParameters mpp129 (cs129, 12116);
    StopNeutralinoMap.insert(std::make_pair(mp129, mpp129));

    MassPoint mp130 (425, 325);
    StopCrossSection cs130 (0.0452189, 0.00371547);
    MassPointParameters mpp130 (cs130, 7767);
    StopNeutralinoMap.insert(std::make_pair(mp130, mpp130));

    MassPoint mp131 (425, 350);
    StopCrossSection cs131 (0.0452189, 0.00371547);
    MassPointParameters mpp131 (cs131, 8984);
    StopNeutralinoMap.insert(std::make_pair(mp131, mpp131));

    MassPoint mp132 (425, 375);
    StopCrossSection cs132 (0.0452189, 0.00371547);
    MassPointParameters mpp132 (cs132, 10291);
    StopNeutralinoMap.insert(std::make_pair(mp132, mpp132));

    MassPoint mp133 (450, 1);
    StopCrossSection cs133 (0.0353143, 0.00297283);
    MassPointParameters mpp133 (cs133, 9496);
    StopNeutralinoMap.insert(std::make_pair(mp133, mpp133));

    MassPoint mp134 (450, 25);
    StopCrossSection cs134 (0.0353143, 0.00297283);
    MassPointParameters mpp134 (cs134, 9840);
    StopNeutralinoMap.insert(std::make_pair(mp134, mpp134));

    MassPoint mp135 (450, 50);
    StopCrossSection cs135 (0.0353143, 0.00297283);
    MassPointParameters mpp135 (cs135, 8628);
    StopNeutralinoMap.insert(std::make_pair(mp135, mpp135));

    MassPoint mp136 (450, 75);
    StopCrossSection cs136 (0.0353143, 0.00297283);
    MassPointParameters mpp136 (cs136, 8139);
    StopNeutralinoMap.insert(std::make_pair(mp136, mpp136));

    MassPoint mp137 (450, 100);
    StopCrossSection cs137 (0.0353143, 0.00297283);
    MassPointParameters mpp137 (cs137, 10629);
    StopNeutralinoMap.insert(std::make_pair(mp137, mpp137));

    MassPoint mp138 (450, 125);
    StopCrossSection cs138 (0.0353143, 0.00297283);
    MassPointParameters mpp138 (cs138, 8914);
    StopNeutralinoMap.insert(std::make_pair(mp138, mpp138));

    MassPoint mp139 (450, 150);
    StopCrossSection cs139 (0.0353143, 0.00297283);
    MassPointParameters mpp139 (cs139, 6636);
    StopNeutralinoMap.insert(std::make_pair(mp139, mpp139));

    MassPoint mp140 (450, 175);
    StopCrossSection cs140 (0.0353143, 0.00297283);
    MassPointParameters mpp140 (cs140, 9281);
    StopNeutralinoMap.insert(std::make_pair(mp140, mpp140));

    MassPoint mp141 (450, 200);
    StopCrossSection cs141 (0.0353143, 0.00297283);
    MassPointParameters mpp141 (cs141, 10347);
    StopNeutralinoMap.insert(std::make_pair(mp141, mpp141));

    MassPoint mp142 (450, 225);
    StopCrossSection cs142 (0.0353143, 0.00297283);
    MassPointParameters mpp142 (cs142, 9520);
    StopNeutralinoMap.insert(std::make_pair(mp142, mpp142));

    MassPoint mp143 (450, 250);
    StopCrossSection cs143 (0.0353143, 0.00297283);
    MassPointParameters mpp143 (cs143, 9963);
    StopNeutralinoMap.insert(std::make_pair(mp143, mpp143));

    MassPoint mp144 (450, 275);
    StopCrossSection cs144 (0.0353143, 0.00297283);
    MassPointParameters mpp144 (cs144, 11475);
    StopNeutralinoMap.insert(std::make_pair(mp144, mpp144));

    MassPoint mp145 (450, 300);
    StopCrossSection cs145 (0.0353143, 0.00297283);
    MassPointParameters mpp145 (cs145, 8759);
    StopNeutralinoMap.insert(std::make_pair(mp145, mpp145));

    MassPoint mp146 (450, 325);
    StopCrossSection cs146 (0.0353143, 0.00297283);
    MassPointParameters mpp146 (cs146, 10892);
    StopNeutralinoMap.insert(std::make_pair(mp146, mpp146));

    MassPoint mp147 (450, 350);
    StopCrossSection cs147 (0.0353143, 0.00297283);
    MassPointParameters mpp147 (cs147, 11120);
    StopNeutralinoMap.insert(std::make_pair(mp147, mpp147));

    MassPoint mp148 (450, 375);
    StopCrossSection cs148 (0.0353143, 0.00297283);
    MassPointParameters mpp148 (cs148, 8672);
    StopNeutralinoMap.insert(std::make_pair(mp148, mpp148));

    MassPoint mp149 (450, 400);
    StopCrossSection cs149 (0.0353143, 0.00297283);
    MassPointParameters mpp149 (cs149, 10438);
    StopNeutralinoMap.insert(std::make_pair(mp149, mpp149));

    MassPoint mp150 (475, 1);
    StopCrossSection cs150 (0.0278342, 0.00241224);
    MassPointParameters mpp150 (cs150, 7249);
    StopNeutralinoMap.insert(std::make_pair(mp150, mpp150));

    MassPoint mp151 (475, 25);
    StopCrossSection cs151 (0.0278342, 0.00241224);
    MassPointParameters mpp151 (cs151, 10422);
    StopNeutralinoMap.insert(std::make_pair(mp151, mpp151));

    MassPoint mp152 (475, 50);
    StopCrossSection cs152 (0.0278342, 0.00241224);
    MassPointParameters mpp152 (cs152, 11239);
    StopNeutralinoMap.insert(std::make_pair(mp152, mpp152));

    MassPoint mp153 (475, 75);
    StopCrossSection cs153 (0.0278342, 0.00241224);
    MassPointParameters mpp153 (cs153, 11610);
    StopNeutralinoMap.insert(std::make_pair(mp153, mpp153));

    MassPoint mp154 (475, 100);
    StopCrossSection cs154 (0.0278342, 0.00241224);
    MassPointParameters mpp154 (cs154, 9638);
    StopNeutralinoMap.insert(std::make_pair(mp154, mpp154));

    MassPoint mp155 (475, 125);
    StopCrossSection cs155 (0.0278342, 0.00241224);
    MassPointParameters mpp155 (cs155, 9129);
    StopNeutralinoMap.insert(std::make_pair(mp155, mpp155));

    MassPoint mp156 (475, 150);
    StopCrossSection cs156 (0.0278342, 0.00241224);
    MassPointParameters mpp156 (cs156, 10317);
    StopNeutralinoMap.insert(std::make_pair(mp156, mpp156));

    MassPoint mp157 (475, 175);
    StopCrossSection cs157 (0.0278342, 0.00241224);
    MassPointParameters mpp157 (cs157, 8420);
    StopNeutralinoMap.insert(std::make_pair(mp157, mpp157));

    MassPoint mp158 (475, 200);
    StopCrossSection cs158 (0.0278342, 0.00241224);
    MassPointParameters mpp158 (cs158, 11173);
    StopNeutralinoMap.insert(std::make_pair(mp158, mpp158));

    MassPoint mp159 (475, 225);
    StopCrossSection cs159 (0.0278342, 0.00241224);
    MassPointParameters mpp159 (cs159, 9265);
    StopNeutralinoMap.insert(std::make_pair(mp159, mpp159));

    MassPoint mp160 (475, 250);
    StopCrossSection cs160 (0.0278342, 0.00241224);
    MassPointParameters mpp160 (cs160, 8553);
    StopNeutralinoMap.insert(std::make_pair(mp160, mpp160));

    MassPoint mp161 (475, 275);
    StopCrossSection cs161 (0.0278342, 0.00241224);
    MassPointParameters mpp161 (cs161, 10697);
    StopNeutralinoMap.insert(std::make_pair(mp161, mpp161));

    MassPoint mp162 (475, 300);
    StopCrossSection cs162 (0.0278342, 0.00241224);
    MassPointParameters mpp162 (cs162, 9275);
    StopNeutralinoMap.insert(std::make_pair(mp162, mpp162));

    MassPoint mp163 (475, 325);
    StopCrossSection cs163 (0.0278342, 0.00241224);
    MassPointParameters mpp163 (cs163, 11217);
    StopNeutralinoMap.insert(std::make_pair(mp163, mpp163));

    MassPoint mp164 (475, 350);
    StopCrossSection cs164 (0.0278342, 0.00241224);
    MassPointParameters mpp164 (cs164, 9039);
    StopNeutralinoMap.insert(std::make_pair(mp164, mpp164));

    MassPoint mp165 (475, 375);
    StopCrossSection cs165 (0.0278342, 0.00241224);
    MassPointParameters mpp165 (cs165, 10089);
    StopNeutralinoMap.insert(std::make_pair(mp165, mpp165));

    MassPoint mp166 (475, 400);
    StopCrossSection cs166 (0.0278342, 0.00241224);
    MassPointParameters mpp166 (cs166, 8235);
    StopNeutralinoMap.insert(std::make_pair(mp166, mpp166));

    MassPoint mp167 (475, 425);
    StopCrossSection cs167 (0.0278342, 0.00241224);
    MassPointParameters mpp167 (cs167, 11470);
    StopNeutralinoMap.insert(std::make_pair(mp167, mpp167));

    MassPoint mp168 (500, 1);
    StopCrossSection cs168 (0.0221265, 0.00194904);
    MassPointParameters mpp168 (cs168, 10959);
    StopNeutralinoMap.insert(std::make_pair(mp168, mpp168));

    MassPoint mp169 (500, 25);
    StopCrossSection cs169 (0.0221265, 0.00194904);
    MassPointParameters mpp169 (cs169, 8939);
    StopNeutralinoMap.insert(std::make_pair(mp169, mpp169));

    MassPoint mp170 (500, 50);
    StopCrossSection cs170 (0.0221265, 0.00194904);
    MassPointParameters mpp170 (cs170, 12301);
    StopNeutralinoMap.insert(std::make_pair(mp170, mpp170));

    MassPoint mp171 (500, 75);
    StopCrossSection cs171 (0.0221265, 0.00194904);
    MassPointParameters mpp171 (cs171, 9597);
    StopNeutralinoMap.insert(std::make_pair(mp171, mpp171));

    MassPoint mp172 (500, 100);
    StopCrossSection cs172 (0.0221265, 0.00194904);
    MassPointParameters mpp172 (cs172, 8051);
    StopNeutralinoMap.insert(std::make_pair(mp172, mpp172));

    MassPoint mp173 (500, 125);
    StopCrossSection cs173 (0.0221265, 0.00194904);
    MassPointParameters mpp173 (cs173, 8768);
    StopNeutralinoMap.insert(std::make_pair(mp173, mpp173));

    MassPoint mp174 (500, 150);
    StopCrossSection cs174 (0.0221265, 0.00194904);
    MassPointParameters mpp174 (cs174, 9610);
    StopNeutralinoMap.insert(std::make_pair(mp174, mpp174));

    MassPoint mp175 (500, 175);
    StopCrossSection cs175 (0.0221265, 0.00194904);
    MassPointParameters mpp175 (cs175, 9347);
    StopNeutralinoMap.insert(std::make_pair(mp175, mpp175));

    MassPoint mp176 (500, 200);
    StopCrossSection cs176 (0.0221265, 0.00194904);
    MassPointParameters mpp176 (cs176, 10709);
    StopNeutralinoMap.insert(std::make_pair(mp176, mpp176));

    MassPoint mp177 (500, 225);
    StopCrossSection cs177 (0.0221265, 0.00194904);
    MassPointParameters mpp177 (cs177, 10744);
    StopNeutralinoMap.insert(std::make_pair(mp177, mpp177));

    MassPoint mp178 (500, 250);
    StopCrossSection cs178 (0.0221265, 0.00194904);
    MassPointParameters mpp178 (cs178, 7921);
    StopNeutralinoMap.insert(std::make_pair(mp178, mpp178));

    MassPoint mp179 (500, 275);
    StopCrossSection cs179 (0.0221265, 0.00194904);
    MassPointParameters mpp179 (cs179, 10542);
    StopNeutralinoMap.insert(std::make_pair(mp179, mpp179));

    MassPoint mp180 (500, 300);
    StopCrossSection cs180 (0.0221265, 0.00194904);
    MassPointParameters mpp180 (cs180, 11685);
    StopNeutralinoMap.insert(std::make_pair(mp180, mpp180));

    MassPoint mp181 (500, 325);
    StopCrossSection cs181 (0.0221265, 0.00194904);
    MassPointParameters mpp181 (cs181, 9177);
    StopNeutralinoMap.insert(std::make_pair(mp181, mpp181));

    MassPoint mp182 (500, 350);
    StopCrossSection cs182 (0.0221265, 0.00194904);
    MassPointParameters mpp182 (cs182, 11439);
    StopNeutralinoMap.insert(std::make_pair(mp182, mpp182));

    MassPoint mp183 (500, 375);
    StopCrossSection cs183 (0.0221265, 0.00194904);
    MassPointParameters mpp183 (cs183, 8450);
    StopNeutralinoMap.insert(std::make_pair(mp183, mpp183));

    MassPoint mp184 (500, 400);
    StopCrossSection cs184 (0.0221265, 0.00194904);
    MassPointParameters mpp184 (cs184, 8119);
    StopNeutralinoMap.insert(std::make_pair(mp184, mpp184));

    MassPoint mp185 (500, 425);
    StopCrossSection cs185 (0.0221265, 0.00194904);
    MassPointParameters mpp185 (cs185, 9269);
    StopNeutralinoMap.insert(std::make_pair(mp185, mpp185));

    MassPoint mp186 (500, 450);
    StopCrossSection cs186 (0.0221265, 0.00194904);
    MassPointParameters mpp186 (cs186, 9182);
    StopNeutralinoMap.insert(std::make_pair(mp186, mpp186));

    MassPoint mp187 (525, 1);
    StopCrossSection cs187 (0.0177394, 0.0015992);
    MassPointParameters mpp187 (cs187, 8608);
    StopNeutralinoMap.insert(std::make_pair(mp187, mpp187));

    MassPoint mp188 (525, 25);
    StopCrossSection cs188 (0.0177394, 0.0015992);
    MassPointParameters mpp188 (cs188, 7459);
    StopNeutralinoMap.insert(std::make_pair(mp188, mpp188));

    MassPoint mp189 (525, 50);
    StopCrossSection cs189 (0.0177394, 0.0015992);
    MassPointParameters mpp189 (cs189, 9103);
    StopNeutralinoMap.insert(std::make_pair(mp189, mpp189));

    MassPoint mp190 (525, 75);
    StopCrossSection cs190 (0.0177394, 0.0015992);
    MassPointParameters mpp190 (cs190, 9180);
    StopNeutralinoMap.insert(std::make_pair(mp190, mpp190));

    MassPoint mp191 (525, 100);
    StopCrossSection cs191 (0.0177394, 0.0015992);
    MassPointParameters mpp191 (cs191, 10620);
    StopNeutralinoMap.insert(std::make_pair(mp191, mpp191));

    MassPoint mp192 (525, 125);
    StopCrossSection cs192 (0.0177394, 0.0015992);
    MassPointParameters mpp192 (cs192, 9554);
    StopNeutralinoMap.insert(std::make_pair(mp192, mpp192));

    MassPoint mp193 (525, 150);
    StopCrossSection cs193 (0.0177394, 0.0015992);
    MassPointParameters mpp193 (cs193, 9098);
    StopNeutralinoMap.insert(std::make_pair(mp193, mpp193));

    MassPoint mp194 (525, 175);
    StopCrossSection cs194 (0.0177394, 0.0015992);
    MassPointParameters mpp194 (cs194, 13788);
    StopNeutralinoMap.insert(std::make_pair(mp194, mpp194));

    MassPoint mp195 (525, 200);
    StopCrossSection cs195 (0.0177394, 0.0015992);
    MassPointParameters mpp195 (cs195, 8656);
    StopNeutralinoMap.insert(std::make_pair(mp195, mpp195));

    MassPoint mp196 (525, 225);
    StopCrossSection cs196 (0.0177394, 0.0015992);
    MassPointParameters mpp196 (cs196, 9977);
    StopNeutralinoMap.insert(std::make_pair(mp196, mpp196));

    MassPoint mp197 (525, 250);
    StopCrossSection cs197 (0.0177394, 0.0015992);
    MassPointParameters mpp197 (cs197, 7321);
    StopNeutralinoMap.insert(std::make_pair(mp197, mpp197));

    MassPoint mp198 (525, 275);
    StopCrossSection cs198 (0.0177394, 0.0015992);
    MassPointParameters mpp198 (cs198, 10389);
    StopNeutralinoMap.insert(std::make_pair(mp198, mpp198));

    MassPoint mp199 (525, 300);
    StopCrossSection cs199 (0.0177394, 0.0015992);
    MassPointParameters mpp199 (cs199, 11743);
    StopNeutralinoMap.insert(std::make_pair(mp199, mpp199));

    MassPoint mp200 (525, 325);
    StopCrossSection cs200 (0.0177394, 0.0015992);
    MassPointParameters mpp200 (cs200, 10820);
    StopNeutralinoMap.insert(std::make_pair(mp200, mpp200));

    MassPoint mp201 (525, 350);
    StopCrossSection cs201 (0.0177394, 0.0015992);
    MassPointParameters mpp201 (cs201, 7900);
    StopNeutralinoMap.insert(std::make_pair(mp201, mpp201));

    MassPoint mp202 (525, 375);
    StopCrossSection cs202 (0.0177394, 0.0015992);
    MassPointParameters mpp202 (cs202, 11208);
    StopNeutralinoMap.insert(std::make_pair(mp202, mpp202));

    MassPoint mp203 (525, 400);
    StopCrossSection cs203 (0.0177394, 0.0015992);
    MassPointParameters mpp203 (cs203, 12322);
    StopNeutralinoMap.insert(std::make_pair(mp203, mpp203));

    MassPoint mp204 (525, 425);
    StopCrossSection cs204 (0.0177394, 0.0015992);
    MassPointParameters mpp204 (cs204, 12071);
    StopNeutralinoMap.insert(std::make_pair(mp204, mpp204));

    MassPoint mp205 (525, 450);
    StopCrossSection cs205 (0.0177394, 0.0015992);
    MassPointParameters mpp205 (cs205, 10739);
    StopNeutralinoMap.insert(std::make_pair(mp205, mpp205));

    MassPoint mp206 (525, 475);
    StopCrossSection cs206 (0.0177394, 0.0015992);
    MassPointParameters mpp206 (cs206, 10219);
    StopNeutralinoMap.insert(std::make_pair(mp206, mpp206));

    MassPoint mp207 (550, 1);
    StopCrossSection cs207 (0.0143134, 0.00132368);
    MassPointParameters mpp207 (cs207, 10881);
    StopNeutralinoMap.insert(std::make_pair(mp207, mpp207));

    MassPoint mp208 (550, 25);
    StopCrossSection cs208 (0.0143134, 0.00132368);
    MassPointParameters mpp208 (cs208, 11462);
    StopNeutralinoMap.insert(std::make_pair(mp208, mpp208));

    MassPoint mp209 (550, 50);
    StopCrossSection cs209 (0.0143134, 0.00132368);
    MassPointParameters mpp209 (cs209, 10430);
    StopNeutralinoMap.insert(std::make_pair(mp209, mpp209));

    MassPoint mp210 (550, 75);
    StopCrossSection cs210 (0.0143134, 0.00132368);
    MassPointParameters mpp210 (cs210, 9232);
    StopNeutralinoMap.insert(std::make_pair(mp210, mpp210));

    MassPoint mp211 (550, 100);
    StopCrossSection cs211 (0.0143134, 0.00132368);
    MassPointParameters mpp211 (cs211, 11523);
    StopNeutralinoMap.insert(std::make_pair(mp211, mpp211));

    MassPoint mp212 (550, 125);
    StopCrossSection cs212 (0.0143134, 0.00132368);
    MassPointParameters mpp212 (cs212, 8869);
    StopNeutralinoMap.insert(std::make_pair(mp212, mpp212));

    MassPoint mp213 (550, 150);
    StopCrossSection cs213 (0.0143134, 0.00132368);
    MassPointParameters mpp213 (cs213, 9229);
    StopNeutralinoMap.insert(std::make_pair(mp213, mpp213));

    MassPoint mp214 (550, 175);
    StopCrossSection cs214 (0.0143134, 0.00132368);
    MassPointParameters mpp214 (cs214, 11196);
    StopNeutralinoMap.insert(std::make_pair(mp214, mpp214));

    MassPoint mp215 (550, 200);
    StopCrossSection cs215 (0.0143134, 0.00132368);
    MassPointParameters mpp215 (cs215, 12472);
    StopNeutralinoMap.insert(std::make_pair(mp215, mpp215));

    MassPoint mp216 (550, 225);
    StopCrossSection cs216 (0.0143134, 0.00132368);
    MassPointParameters mpp216 (cs216, 8643);
    StopNeutralinoMap.insert(std::make_pair(mp216, mpp216));

    MassPoint mp217 (550, 250);
    StopCrossSection cs217 (0.0143134, 0.00132368);
    MassPointParameters mpp217 (cs217, 9345);
    StopNeutralinoMap.insert(std::make_pair(mp217, mpp217));

    MassPoint mp218 (550, 275);
    StopCrossSection cs218 (0.0143134, 0.00132368);
    MassPointParameters mpp218 (cs218, 8314);
    StopNeutralinoMap.insert(std::make_pair(mp218, mpp218));

    MassPoint mp219 (550, 300);
    StopCrossSection cs219 (0.0143134, 0.00132368);
    MassPointParameters mpp219 (cs219, 8842);
    StopNeutralinoMap.insert(std::make_pair(mp219, mpp219));

    MassPoint mp220 (550, 325);
    StopCrossSection cs220 (0.0143134, 0.00132368);
    MassPointParameters mpp220 (cs220, 10318);
    StopNeutralinoMap.insert(std::make_pair(mp220, mpp220));

    MassPoint mp221 (550, 350);
    StopCrossSection cs221 (0.0143134, 0.00132368);
    MassPointParameters mpp221 (cs221, 10006);
    StopNeutralinoMap.insert(std::make_pair(mp221, mpp221));

    MassPoint mp222 (550, 375);
    StopCrossSection cs222 (0.0143134, 0.00132368);
    MassPointParameters mpp222 (cs222, 10379);
    StopNeutralinoMap.insert(std::make_pair(mp222, mpp222));

    MassPoint mp223 (550, 400);
    StopCrossSection cs223 (0.0143134, 0.00132368);
    MassPointParameters mpp223 (cs223, 9028);
    StopNeutralinoMap.insert(std::make_pair(mp223, mpp223));

    MassPoint mp224 (550, 425);
    StopCrossSection cs224 (0.0143134, 0.00132368);
    MassPointParameters mpp224 (cs224, 9750);
    StopNeutralinoMap.insert(std::make_pair(mp224, mpp224));

    MassPoint mp225 (550, 450);
    StopCrossSection cs225 (0.0143134, 0.00132368);
    MassPointParameters mpp225 (cs225, 9125);
    StopNeutralinoMap.insert(std::make_pair(mp225, mpp225));

    MassPoint mp226 (550, 475);
    StopCrossSection cs226 (0.0143134, 0.00132368);
    MassPointParameters mpp226 (cs226, 10044);
    StopNeutralinoMap.insert(std::make_pair(mp226, mpp226));

    MassPoint mp227 (550, 500);
    StopCrossSection cs227 (0.0143134, 0.00132368);
    MassPointParameters mpp227 (cs227, 10171);
    StopNeutralinoMap.insert(std::make_pair(mp227, mpp227));

    MassPoint mp228 (575, 1);
    StopCrossSection cs228 (0.0116266, 0.00109669);
    MassPointParameters mpp228 (cs228, 10829);
    StopNeutralinoMap.insert(std::make_pair(mp228, mpp228));

    MassPoint mp229 (575, 25);
    StopCrossSection cs229 (0.0116266, 0.00109669);
    MassPointParameters mpp229 (cs229, 11052);
    StopNeutralinoMap.insert(std::make_pair(mp229, mpp229));

    MassPoint mp230 (575, 50);
    StopCrossSection cs230 (0.0116266, 0.00109669);
    MassPointParameters mpp230 (cs230, 10915);
    StopNeutralinoMap.insert(std::make_pair(mp230, mpp230));

    MassPoint mp231 (575, 75);
    StopCrossSection cs231 (0.0116266, 0.00109669);
    MassPointParameters mpp231 (cs231, 10234);
    StopNeutralinoMap.insert(std::make_pair(mp231, mpp231));

    MassPoint mp232 (575, 100);
    StopCrossSection cs232 (0.0116266, 0.00109669);
    MassPointParameters mpp232 (cs232, 7214);
    StopNeutralinoMap.insert(std::make_pair(mp232, mpp232));

    MassPoint mp233 (575, 125);
    StopCrossSection cs233 (0.0116266, 0.00109669);
    MassPointParameters mpp233 (cs233, 11342);
    StopNeutralinoMap.insert(std::make_pair(mp233, mpp233));

    MassPoint mp234 (575, 150);
    StopCrossSection cs234 (0.0116266, 0.00109669);
    MassPointParameters mpp234 (cs234, 10086);
    StopNeutralinoMap.insert(std::make_pair(mp234, mpp234));

    MassPoint mp235 (575, 175);
    StopCrossSection cs235 (0.0116266, 0.00109669);
    MassPointParameters mpp235 (cs235, 10994);
    StopNeutralinoMap.insert(std::make_pair(mp235, mpp235));

    MassPoint mp236 (575, 200);
    StopCrossSection cs236 (0.0116266, 0.00109669);
    MassPointParameters mpp236 (cs236, 9767);
    StopNeutralinoMap.insert(std::make_pair(mp236, mpp236));

    MassPoint mp237 (575, 225);
    StopCrossSection cs237 (0.0116266, 0.00109669);
    MassPointParameters mpp237 (cs237, 8807);
    StopNeutralinoMap.insert(std::make_pair(mp237, mpp237));

    MassPoint mp238 (575, 250);
    StopCrossSection cs238 (0.0116266, 0.00109669);
    MassPointParameters mpp238 (cs238, 10007);
    StopNeutralinoMap.insert(std::make_pair(mp238, mpp238));

    MassPoint mp239 (575, 275);
    StopCrossSection cs239 (0.0116266, 0.00109669);
    MassPointParameters mpp239 (cs239, 9540);
    StopNeutralinoMap.insert(std::make_pair(mp239, mpp239));

    MassPoint mp240 (575, 300);
    StopCrossSection cs240 (0.0116266, 0.00109669);
    MassPointParameters mpp240 (cs240, 9325);
    StopNeutralinoMap.insert(std::make_pair(mp240, mpp240));

    MassPoint mp241 (575, 325);
    StopCrossSection cs241 (0.0116266, 0.00109669);
    MassPointParameters mpp241 (cs241, 10815);
    StopNeutralinoMap.insert(std::make_pair(mp241, mpp241));

    MassPoint mp242 (575, 350);
    StopCrossSection cs242 (0.0116266, 0.00109669);
    MassPointParameters mpp242 (cs242, 7039);
    StopNeutralinoMap.insert(std::make_pair(mp242, mpp242));

    MassPoint mp243 (575, 375);
    StopCrossSection cs243 (0.0116266, 0.00109669);
    MassPointParameters mpp243 (cs243, 10835);
    StopNeutralinoMap.insert(std::make_pair(mp243, mpp243));

    MassPoint mp244 (575, 400);
    StopCrossSection cs244 (0.0116266, 0.00109669);
    MassPointParameters mpp244 (cs244, 8270);
    StopNeutralinoMap.insert(std::make_pair(mp244, mpp244));

    MassPoint mp245 (575, 425);
    StopCrossSection cs245 (0.0116266, 0.00109669);
    MassPointParameters mpp245 (cs245, 9965);
    StopNeutralinoMap.insert(std::make_pair(mp245, mpp245));

    MassPoint mp246 (575, 450);
    StopCrossSection cs246 (0.0116266, 0.00109669);
    MassPointParameters mpp246 (cs246, 8679);
    StopNeutralinoMap.insert(std::make_pair(mp246, mpp246));

    MassPoint mp247 (575, 475);
    StopCrossSection cs247 (0.0116266, 0.00109669);
    MassPointParameters mpp247 (cs247, 10320);
    StopNeutralinoMap.insert(std::make_pair(mp247, mpp247));

    MassPoint mp248 (575, 500);
    StopCrossSection cs248 (0.0116266, 0.00109669);
    MassPointParameters mpp248 (cs248, 11639);
    StopNeutralinoMap.insert(std::make_pair(mp248, mpp248));

    MassPoint mp249 (575, 525);
    StopCrossSection cs249 (0.0116266, 0.00109669);
    MassPointParameters mpp249 (cs249, 7738);
    StopNeutralinoMap.insert(std::make_pair(mp249, mpp249));

    MassPoint mp250 (600, 1);
    StopCrossSection cs250 (0.00949913, 0.000912324);
    MassPointParameters mpp250 (cs250, 10527);
    StopNeutralinoMap.insert(std::make_pair(mp250, mpp250));

    MassPoint mp251 (600, 25);
    StopCrossSection cs251 (0.00949913, 0.000912324);
    MassPointParameters mpp251 (cs251, 13538);
    StopNeutralinoMap.insert(std::make_pair(mp251, mpp251));

    MassPoint mp252 (600, 50);
    StopCrossSection cs252 (0.00949913, 0.000912324);
    MassPointParameters mpp252 (cs252, 9501);
    StopNeutralinoMap.insert(std::make_pair(mp252, mpp252));

    MassPoint mp253 (600, 75);
    StopCrossSection cs253 (0.00949913, 0.000912324);
    MassPointParameters mpp253 (cs253, 10430);
    StopNeutralinoMap.insert(std::make_pair(mp253, mpp253));

    MassPoint mp254 (600, 100);
    StopCrossSection cs254 (0.00949913, 0.000912324);
    MassPointParameters mpp254 (cs254, 10030);
    StopNeutralinoMap.insert(std::make_pair(mp254, mpp254));

    MassPoint mp255 (600, 125);
    StopCrossSection cs255 (0.00949913, 0.000912324);
    MassPointParameters mpp255 (cs255, 10586);
    StopNeutralinoMap.insert(std::make_pair(mp255, mpp255));

    MassPoint mp256 (600, 150);
    StopCrossSection cs256 (0.00949913, 0.000912324);
    MassPointParameters mpp256 (cs256, 11202);
    StopNeutralinoMap.insert(std::make_pair(mp256, mpp256));

    MassPoint mp257 (600, 175);
    StopCrossSection cs257 (0.00949913, 0.000912324);
    MassPointParameters mpp257 (cs257, 10170);
    StopNeutralinoMap.insert(std::make_pair(mp257, mpp257));

    MassPoint mp258 (600, 200);
    StopCrossSection cs258 (0.00949913, 0.000912324);
    MassPointParameters mpp258 (cs258, 8913);
    StopNeutralinoMap.insert(std::make_pair(mp258, mpp258));

    MassPoint mp259 (600, 225);
    StopCrossSection cs259 (0.00949913, 0.000912324);
    MassPointParameters mpp259 (cs259, 8247);
    StopNeutralinoMap.insert(std::make_pair(mp259, mpp259));

    MassPoint mp260 (600, 250);
    StopCrossSection cs260 (0.00949913, 0.000912324);
    MassPointParameters mpp260 (cs260, 12189);
    StopNeutralinoMap.insert(std::make_pair(mp260, mpp260));

    MassPoint mp261 (600, 275);
    StopCrossSection cs261 (0.00949913, 0.000912324);
    MassPointParameters mpp261 (cs261, 10313);
    StopNeutralinoMap.insert(std::make_pair(mp261, mpp261));

    MassPoint mp262 (600, 300);
    StopCrossSection cs262 (0.00949913, 0.000912324);
    MassPointParameters mpp262 (cs262, 8271);
    StopNeutralinoMap.insert(std::make_pair(mp262, mpp262));

    MassPoint mp263 (600, 325);
    StopCrossSection cs263 (0.00949913, 0.000912324);
    MassPointParameters mpp263 (cs263, 11083);
    StopNeutralinoMap.insert(std::make_pair(mp263, mpp263));

    MassPoint mp264 (600, 350);
    StopCrossSection cs264 (0.00949913, 0.000912324);
    MassPointParameters mpp264 (cs264, 9276);
    StopNeutralinoMap.insert(std::make_pair(mp264, mpp264));

    MassPoint mp265 (600, 375);
    StopCrossSection cs265 (0.00949913, 0.000912324);
    MassPointParameters mpp265 (cs265, 8672);
    StopNeutralinoMap.insert(std::make_pair(mp265, mpp265));

    MassPoint mp266 (600, 400);
    StopCrossSection cs266 (0.00949913, 0.000912324);
    MassPointParameters mpp266 (cs266, 9720);
    StopNeutralinoMap.insert(std::make_pair(mp266, mpp266));

    MassPoint mp267 (600, 425);
    StopCrossSection cs267 (0.00949913, 0.000912324);
    MassPointParameters mpp267 (cs267, 6899);
    StopNeutralinoMap.insert(std::make_pair(mp267, mpp267));

    MassPoint mp268 (600, 450);
    StopCrossSection cs268 (0.00949913, 0.000912324);
    MassPointParameters mpp268 (cs268, 9472);
    StopNeutralinoMap.insert(std::make_pair(mp268, mpp268));

    MassPoint mp269 (600, 475);
    StopCrossSection cs269 (0.00949913, 0.000912324);
    MassPointParameters mpp269 (cs269, 14832);
    StopNeutralinoMap.insert(std::make_pair(mp269, mpp269));

    MassPoint mp270 (600, 500);
    StopCrossSection cs270 (0.00949913, 0.000912324);
    MassPointParameters mpp270 (cs270, 10062);
    StopNeutralinoMap.insert(std::make_pair(mp270, mpp270));

    MassPoint mp271 (600, 525);
    StopCrossSection cs271 (0.00949913, 0.000912324);
    MassPointParameters mpp271 (cs271, 10175);
    StopNeutralinoMap.insert(std::make_pair(mp271, mpp271));

    MassPoint mp272 (600, 550);
    StopCrossSection cs272 (0.00949913, 0.000912324);
    MassPointParameters mpp272 (cs272, 9599);
    StopNeutralinoMap.insert(std::make_pair(mp272, mpp272));

    MassPoint mp273 (625, 1);
    StopCrossSection cs273 (0.00780081, 0.000768988);
    MassPointParameters mpp273 (cs273, 10223);
    StopNeutralinoMap.insert(std::make_pair(mp273, mpp273));

    MassPoint mp274 (625, 25);
    StopCrossSection cs274 (0.00780081, 0.000768988);
    MassPointParameters mpp274 (cs274, 9798);
    StopNeutralinoMap.insert(std::make_pair(mp274, mpp274));

    MassPoint mp275 (625, 50);
    StopCrossSection cs275 (0.00780081, 0.000768988);
    MassPointParameters mpp275 (cs275, 8478);
    StopNeutralinoMap.insert(std::make_pair(mp275, mpp275));

    MassPoint mp276 (625, 75);
    StopCrossSection cs276 (0.00780081, 0.000768988);
    MassPointParameters mpp276 (cs276, 8564);
    StopNeutralinoMap.insert(std::make_pair(mp276, mpp276));

    MassPoint mp277 (625, 100);
    StopCrossSection cs277 (0.00780081, 0.000768988);
    MassPointParameters mpp277 (cs277, 7793);
    StopNeutralinoMap.insert(std::make_pair(mp277, mpp277));

    MassPoint mp278 (625, 125);
    StopCrossSection cs278 (0.00780081, 0.000768988);
    MassPointParameters mpp278 (cs278, 9113);
    StopNeutralinoMap.insert(std::make_pair(mp278, mpp278));

    MassPoint mp279 (625, 150);
    StopCrossSection cs279 (0.00780081, 0.000768988);
    MassPointParameters mpp279 (cs279, 9434);
    StopNeutralinoMap.insert(std::make_pair(mp279, mpp279));

    MassPoint mp280 (625, 175);
    StopCrossSection cs280 (0.00780081, 0.000768988);
    MassPointParameters mpp280 (cs280, 10009);
    StopNeutralinoMap.insert(std::make_pair(mp280, mpp280));

    MassPoint mp281 (625, 200);
    StopCrossSection cs281 (0.00780081, 0.000768988);
    MassPointParameters mpp281 (cs281, 7091);
    StopNeutralinoMap.insert(std::make_pair(mp281, mpp281));

    MassPoint mp282 (625, 225);
    StopCrossSection cs282 (0.00780081, 0.000768988);
    MassPointParameters mpp282 (cs282, 9006);
    StopNeutralinoMap.insert(std::make_pair(mp282, mpp282));

    MassPoint mp283 (625, 250);
    StopCrossSection cs283 (0.00780081, 0.000768988);
    MassPointParameters mpp283 (cs283, 8233);
    StopNeutralinoMap.insert(std::make_pair(mp283, mpp283));

    MassPoint mp284 (625, 275);
    StopCrossSection cs284 (0.00780081, 0.000768988);
    MassPointParameters mpp284 (cs284, 10373);
    StopNeutralinoMap.insert(std::make_pair(mp284, mpp284));

    MassPoint mp285 (625, 300);
    StopCrossSection cs285 (0.00780081, 0.000768988);
    MassPointParameters mpp285 (cs285, 9558);
    StopNeutralinoMap.insert(std::make_pair(mp285, mpp285));

    MassPoint mp286 (625, 325);
    StopCrossSection cs286 (0.00780081, 0.000768988);
    MassPointParameters mpp286 (cs286, 10825);
    StopNeutralinoMap.insert(std::make_pair(mp286, mpp286));

    MassPoint mp287 (625, 350);
    StopCrossSection cs287 (0.00780081, 0.000768988);
    MassPointParameters mpp287 (cs287, 10064);
    StopNeutralinoMap.insert(std::make_pair(mp287, mpp287));

    MassPoint mp288 (625, 375);
    StopCrossSection cs288 (0.00780081, 0.000768988);
    MassPointParameters mpp288 (cs288, 10663);
    StopNeutralinoMap.insert(std::make_pair(mp288, mpp288));

    MassPoint mp289 (625, 400);
    StopCrossSection cs289 (0.00780081, 0.000768988);
    MassPointParameters mpp289 (cs289, 8142);
    StopNeutralinoMap.insert(std::make_pair(mp289, mpp289));

    MassPoint mp290 (625, 425);
    StopCrossSection cs290 (0.00780081, 0.000768988);
    MassPointParameters mpp290 (cs290, 8315);
    StopNeutralinoMap.insert(std::make_pair(mp290, mpp290));

    MassPoint mp291 (625, 450);
    StopCrossSection cs291 (0.00780081, 0.000768988);
    MassPointParameters mpp291 (cs291, 8839);
    StopNeutralinoMap.insert(std::make_pair(mp291, mpp291));

    MassPoint mp292 (625, 475);
    StopCrossSection cs292 (0.00780081, 0.000768988);
    MassPointParameters mpp292 (cs292, 8679);
    StopNeutralinoMap.insert(std::make_pair(mp292, mpp292));

    MassPoint mp293 (625, 500);
    StopCrossSection cs293 (0.00780081, 0.000768988);
    MassPointParameters mpp293 (cs293, 9274);
    StopNeutralinoMap.insert(std::make_pair(mp293, mpp293));

    MassPoint mp294 (625, 525);
    StopCrossSection cs294 (0.00780081, 0.000768988);
    MassPointParameters mpp294 (cs294, 11146);
    StopNeutralinoMap.insert(std::make_pair(mp294, mpp294));

    MassPoint mp295 (625, 550);
    StopCrossSection cs295 (0.00780081, 0.000768988);
    MassPointParameters mpp295 (cs295, 6998);
    StopNeutralinoMap.insert(std::make_pair(mp295, mpp295));

    MassPoint mp296 (625, 575);
    StopCrossSection cs296 (0.00780081, 0.000768988);
    MassPointParameters mpp296 (cs296, 9727);
    StopNeutralinoMap.insert(std::make_pair(mp296, mpp296));

    MassPoint mp297 (650, 1);
    StopCrossSection cs297 (0.00643244, 0.000638889);
    MassPointParameters mpp297 (cs297, 8402);
    StopNeutralinoMap.insert(std::make_pair(mp297, mpp297));

    MassPoint mp298 (650, 25);
    StopCrossSection cs298 (0.00643244, 0.000638889);
    MassPointParameters mpp298 (cs298, 8728);
    StopNeutralinoMap.insert(std::make_pair(mp298, mpp298));

    MassPoint mp299 (650, 50);
    StopCrossSection cs299 (0.00643244, 0.000638889);
    MassPointParameters mpp299 (cs299, 9169);
    StopNeutralinoMap.insert(std::make_pair(mp299, mpp299));

    MassPoint mp300 (650, 75);
    StopCrossSection cs300 (0.00643244, 0.000638889);
    MassPointParameters mpp300 (cs300, 8729);
    StopNeutralinoMap.insert(std::make_pair(mp300, mpp300));

    MassPoint mp301 (650, 100);
    StopCrossSection cs301 (0.00643244, 0.000638889);
    MassPointParameters mpp301 (cs301, 10072);
    StopNeutralinoMap.insert(std::make_pair(mp301, mpp301));

    MassPoint mp302 (650, 125);
    StopCrossSection cs302 (0.00643244, 0.000638889);
    MassPointParameters mpp302 (cs302, 8031);
    StopNeutralinoMap.insert(std::make_pair(mp302, mpp302));

    MassPoint mp303 (650, 150);
    StopCrossSection cs303 (0.00643244, 0.000638889);
    MassPointParameters mpp303 (cs303, 8850);
    StopNeutralinoMap.insert(std::make_pair(mp303, mpp303));

    MassPoint mp304 (650, 175);
    StopCrossSection cs304 (0.00643244, 0.000638889);
    MassPointParameters mpp304 (cs304, 9667);
    StopNeutralinoMap.insert(std::make_pair(mp304, mpp304));

    MassPoint mp305 (650, 200);
    StopCrossSection cs305 (0.00643244, 0.000638889);
    MassPointParameters mpp305 (cs305, 11324);
    StopNeutralinoMap.insert(std::make_pair(mp305, mpp305));

    MassPoint mp306 (650, 225);
    StopCrossSection cs306 (0.00643244, 0.000638889);
    MassPointParameters mpp306 (cs306, 8605);
    StopNeutralinoMap.insert(std::make_pair(mp306, mpp306));

    MassPoint mp307 (650, 250);
    StopCrossSection cs307 (0.00643244, 0.000638889);
    MassPointParameters mpp307 (cs307, 10657);
    StopNeutralinoMap.insert(std::make_pair(mp307, mpp307));

    MassPoint mp308 (650, 275);
    StopCrossSection cs308 (0.00643244, 0.000638889);
    MassPointParameters mpp308 (cs308, 12414);
    StopNeutralinoMap.insert(std::make_pair(mp308, mpp308));

    MassPoint mp309 (650, 300);
    StopCrossSection cs309 (0.00643244, 0.000638889);
    MassPointParameters mpp309 (cs309, 7558);
    StopNeutralinoMap.insert(std::make_pair(mp309, mpp309));

    MassPoint mp310 (650, 325);
    StopCrossSection cs310 (0.00643244, 0.000638889);
    MassPointParameters mpp310 (cs310, 10587);
    StopNeutralinoMap.insert(std::make_pair(mp310, mpp310));

    MassPoint mp311 (650, 350);
    StopCrossSection cs311 (0.00643244, 0.000638889);
    MassPointParameters mpp311 (cs311, 8277);
    StopNeutralinoMap.insert(std::make_pair(mp311, mpp311));

    MassPoint mp312 (650, 375);
    StopCrossSection cs312 (0.00643244, 0.000638889);
    MassPointParameters mpp312 (cs312, 10160);
    StopNeutralinoMap.insert(std::make_pair(mp312, mpp312));

    MassPoint mp313 (650, 400);
    StopCrossSection cs313 (0.00643244, 0.000638889);
    MassPointParameters mpp313 (cs313, 11045);
    StopNeutralinoMap.insert(std::make_pair(mp313, mpp313));

    MassPoint mp314 (650, 425);
    StopCrossSection cs314 (0.00643244, 0.000638889);
    MassPointParameters mpp314 (cs314, 10574);
    StopNeutralinoMap.insert(std::make_pair(mp314, mpp314));

    MassPoint mp315 (650, 450);
    StopCrossSection cs315 (0.00643244, 0.000638889);
    MassPointParameters mpp315 (cs315, 11600);
    StopNeutralinoMap.insert(std::make_pair(mp315, mpp315));

    MassPoint mp316 (650, 475);
    StopCrossSection cs316 (0.00643244, 0.000638889);
    MassPointParameters mpp316 (cs316, 10152);
    StopNeutralinoMap.insert(std::make_pair(mp316, mpp316));

    MassPoint mp317 (650, 500);
    StopCrossSection cs317 (0.00643244, 0.000638889);
    MassPointParameters mpp317 (cs317, 9881);
    StopNeutralinoMap.insert(std::make_pair(mp317, mpp317));

    MassPoint mp318 (650, 525);
    StopCrossSection cs318 (0.00643244, 0.000638889);
    MassPointParameters mpp318 (cs318, 11143);
    StopNeutralinoMap.insert(std::make_pair(mp318, mpp318));

    MassPoint mp319 (650, 550);
    StopCrossSection cs319 (0.00643244, 0.000638889);
    MassPointParameters mpp319 (cs319, 8953);
    StopNeutralinoMap.insert(std::make_pair(mp319, mpp319));

    MassPoint mp320 (650, 575);
    StopCrossSection cs320 (0.00643244, 0.000638889);
    MassPointParameters mpp320 (cs320, 10989);
    StopNeutralinoMap.insert(std::make_pair(mp320, mpp320));

    MassPoint mp321 (650, 600);
    StopCrossSection cs321 (0.00643244, 0.000638889);
    MassPointParameters mpp321 (cs321, 11589);
    StopNeutralinoMap.insert(std::make_pair(mp321, mpp321));

    MassPoint mp322 (675, 1);
    StopCrossSection cs322 (0.00533642, 0.000541519);
    MassPointParameters mpp322 (cs322, 8289);
    StopNeutralinoMap.insert(std::make_pair(mp322, mpp322));

    MassPoint mp323 (675, 25);
    StopCrossSection cs323 (0.00533642, 0.000541519);
    MassPointParameters mpp323 (cs323, 9664);
    StopNeutralinoMap.insert(std::make_pair(mp323, mpp323));

    MassPoint mp324 (675, 50);
    StopCrossSection cs324 (0.00533642, 0.000541519);
    MassPointParameters mpp324 (cs324, 9956);
    StopNeutralinoMap.insert(std::make_pair(mp324, mpp324));

    MassPoint mp325 (675, 75);
    StopCrossSection cs325 (0.00533642, 0.000541519);
    MassPointParameters mpp325 (cs325, 10132);
    StopNeutralinoMap.insert(std::make_pair(mp325, mpp325));

    MassPoint mp326 (675, 100);
    StopCrossSection cs326 (0.00533642, 0.000541519);
    MassPointParameters mpp326 (cs326, 10810);
    StopNeutralinoMap.insert(std::make_pair(mp326, mpp326));

    MassPoint mp327 (675, 125);
    StopCrossSection cs327 (0.00533642, 0.000541519);
    MassPointParameters mpp327 (cs327, 10863);
    StopNeutralinoMap.insert(std::make_pair(mp327, mpp327));

    MassPoint mp328 (675, 150);
    StopCrossSection cs328 (0.00533642, 0.000541519);
    MassPointParameters mpp328 (cs328, 10569);
    StopNeutralinoMap.insert(std::make_pair(mp328, mpp328));

    MassPoint mp329 (675, 175);
    StopCrossSection cs329 (0.00533642, 0.000541519);
    MassPointParameters mpp329 (cs329, 7268);
    StopNeutralinoMap.insert(std::make_pair(mp329, mpp329));

    MassPoint mp330 (675, 200);
    StopCrossSection cs330 (0.00533642, 0.000541519);
    MassPointParameters mpp330 (cs330, 10793);
    StopNeutralinoMap.insert(std::make_pair(mp330, mpp330));

    MassPoint mp331 (675, 225);
    StopCrossSection cs331 (0.00533642, 0.000541519);
    MassPointParameters mpp331 (cs331, 12232);
    StopNeutralinoMap.insert(std::make_pair(mp331, mpp331));

    MassPoint mp332 (675, 250);
    StopCrossSection cs332 (0.00533642, 0.000541519);
    MassPointParameters mpp332 (cs332, 11222);
    StopNeutralinoMap.insert(std::make_pair(mp332, mpp332));

    MassPoint mp333 (675, 275);
    StopCrossSection cs333 (0.00533642, 0.000541519);
    MassPointParameters mpp333 (cs333, 8038);
    StopNeutralinoMap.insert(std::make_pair(mp333, mpp333));

    MassPoint mp334 (675, 300);
    StopCrossSection cs334 (0.00533642, 0.000541519);
    MassPointParameters mpp334 (cs334, 9343);
    StopNeutralinoMap.insert(std::make_pair(mp334, mpp334));

    MassPoint mp335 (675, 325);
    StopCrossSection cs335 (0.00533642, 0.000541519);
    MassPointParameters mpp335 (cs335, 7430);
    StopNeutralinoMap.insert(std::make_pair(mp335, mpp335));

    MassPoint mp336 (675, 350);
    StopCrossSection cs336 (0.00533642, 0.000541519);
    MassPointParameters mpp336 (cs336, 8925);
    StopNeutralinoMap.insert(std::make_pair(mp336, mpp336));

    MassPoint mp337 (675, 375);
    StopCrossSection cs337 (0.00533642, 0.000541519);
    MassPointParameters mpp337 (cs337, 10322);
    StopNeutralinoMap.insert(std::make_pair(mp337, mpp337));

    MassPoint mp338 (675, 400);
    StopCrossSection cs338 (0.00533642, 0.000541519);
    MassPointParameters mpp338 (cs338, 9039);
    StopNeutralinoMap.insert(std::make_pair(mp338, mpp338));

    MassPoint mp339 (675, 425);
    StopCrossSection cs339 (0.00533642, 0.000541519);
    MassPointParameters mpp339 (cs339, 9165);
    StopNeutralinoMap.insert(std::make_pair(mp339, mpp339));

    MassPoint mp340 (675, 450);
    StopCrossSection cs340 (0.00533642, 0.000541519);
    MassPointParameters mpp340 (cs340, 8510);
    StopNeutralinoMap.insert(std::make_pair(mp340, mpp340));

    MassPoint mp341 (675, 475);
    StopCrossSection cs341 (0.00533642, 0.000541519);
    MassPointParameters mpp341 (cs341, 8536);
    StopNeutralinoMap.insert(std::make_pair(mp341, mpp341));

    MassPoint mp342 (675, 500);
    StopCrossSection cs342 (0.00533642, 0.000541519);
    MassPointParameters mpp342 (cs342, 10106);
    StopNeutralinoMap.insert(std::make_pair(mp342, mpp342));

    MassPoint mp343 (675, 525);
    StopCrossSection cs343 (0.00533642, 0.000541519);
    MassPointParameters mpp343 (cs343, 10390);
    StopNeutralinoMap.insert(std::make_pair(mp343, mpp343));

    MassPoint mp344 (675, 550);
    StopCrossSection cs344 (0.00533642, 0.000541519);
    MassPointParameters mpp344 (cs344, 8160);
    StopNeutralinoMap.insert(std::make_pair(mp344, mpp344));

    MassPoint mp345 (675, 575);
    StopCrossSection cs345 (0.00533642, 0.000541519);
    MassPointParameters mpp345 (cs345, 9999);
    StopNeutralinoMap.insert(std::make_pair(mp345, mpp345));

    MassPoint mp346 (675, 600);
    StopCrossSection cs346 (0.00533642, 0.000541519);
    MassPointParameters mpp346 (cs346, 11630);
    StopNeutralinoMap.insert(std::make_pair(mp346, mpp346));

    MassPoint mp347 (700, 1);
    StopCrossSection cs347 (0.0044387, 0.000457123);
    MassPointParameters mpp347 (cs347, 9518);
    StopNeutralinoMap.insert(std::make_pair(mp347, mpp347));

    MassPoint mp348 (700, 25);
    StopCrossSection cs348 (0.0044387, 0.000457123);
    MassPointParameters mpp348 (cs348, 8403);
    StopNeutralinoMap.insert(std::make_pair(mp348, mpp348));

    MassPoint mp349 (700, 50);
    StopCrossSection cs349 (0.0044387, 0.000457123);
    MassPointParameters mpp349 (cs349, 11377);
    StopNeutralinoMap.insert(std::make_pair(mp349, mpp349));

    MassPoint mp350 (700, 75);
    StopCrossSection cs350 (0.0044387, 0.000457123);
    MassPointParameters mpp350 (cs350, 11289);
    StopNeutralinoMap.insert(std::make_pair(mp350, mpp350));

    MassPoint mp351 (700, 100);
    StopCrossSection cs351 (0.0044387, 0.000457123);
    MassPointParameters mpp351 (cs351, 7889);
    StopNeutralinoMap.insert(std::make_pair(mp351, mpp351));

    MassPoint mp352 (700, 125);
    StopCrossSection cs352 (0.0044387, 0.000457123);
    MassPointParameters mpp352 (cs352, 10619);
    StopNeutralinoMap.insert(std::make_pair(mp352, mpp352));

    MassPoint mp353 (700, 150);
    StopCrossSection cs353 (0.0044387, 0.000457123);
    MassPointParameters mpp353 (cs353, 9114);
    StopNeutralinoMap.insert(std::make_pair(mp353, mpp353));

    MassPoint mp354 (700, 175);
    StopCrossSection cs354 (0.0044387, 0.000457123);
    MassPointParameters mpp354 (cs354, 9259);
    StopNeutralinoMap.insert(std::make_pair(mp354, mpp354));

    MassPoint mp355 (700, 200);
    StopCrossSection cs355 (0.0044387, 0.000457123);
    MassPointParameters mpp355 (cs355, 10139);
    StopNeutralinoMap.insert(std::make_pair(mp355, mpp355));

    MassPoint mp356 (700, 225);
    StopCrossSection cs356 (0.0044387, 0.000457123);
    MassPointParameters mpp356 (cs356, 10465);
    StopNeutralinoMap.insert(std::make_pair(mp356, mpp356));

    MassPoint mp357 (700, 250);
    StopCrossSection cs357 (0.0044387, 0.000457123);
    MassPointParameters mpp357 (cs357, 9865);
    StopNeutralinoMap.insert(std::make_pair(mp357, mpp357));

    MassPoint mp358 (700, 275);
    StopCrossSection cs358 (0.0044387, 0.000457123);
    MassPointParameters mpp358 (cs358, 10170);
    StopNeutralinoMap.insert(std::make_pair(mp358, mpp358));

    MassPoint mp359 (700, 300);
    StopCrossSection cs359 (0.0044387, 0.000457123);
    MassPointParameters mpp359 (cs359, 11538);
    StopNeutralinoMap.insert(std::make_pair(mp359, mpp359));

    MassPoint mp360 (700, 325);
    StopCrossSection cs360 (0.0044387, 0.000457123);
    MassPointParameters mpp360 (cs360, 10110);
    StopNeutralinoMap.insert(std::make_pair(mp360, mpp360));

    MassPoint mp361 (700, 350);
    StopCrossSection cs361 (0.0044387, 0.000457123);
    MassPointParameters mpp361 (cs361, 9014);
    StopNeutralinoMap.insert(std::make_pair(mp361, mpp361));

    MassPoint mp362 (700, 375);
    StopCrossSection cs362 (0.0044387, 0.000457123);
    MassPointParameters mpp362 (cs362, 9739);
    StopNeutralinoMap.insert(std::make_pair(mp362, mpp362));

    MassPoint mp363 (700, 400);
    StopCrossSection cs363 (0.0044387, 0.000457123);
    MassPointParameters mpp363 (cs363, 7566);
    StopNeutralinoMap.insert(std::make_pair(mp363, mpp363));

    MassPoint mp364 (700, 425);
    StopCrossSection cs364 (0.0044387, 0.000457123);
    MassPointParameters mpp364 (cs364, 8382);
    StopNeutralinoMap.insert(std::make_pair(mp364, mpp364));

    MassPoint mp365 (700, 450);
    StopCrossSection cs365 (0.0044387, 0.000457123);
    MassPointParameters mpp365 (cs365, 9861);
    StopNeutralinoMap.insert(std::make_pair(mp365, mpp365));

    MassPoint mp366 (700, 475);
    StopCrossSection cs366 (0.0044387, 0.000457123);
    MassPointParameters mpp366 (cs366, 10184);
    StopNeutralinoMap.insert(std::make_pair(mp366, mpp366));

    MassPoint mp367 (700, 500);
    StopCrossSection cs367 (0.0044387, 0.000457123);
    MassPointParameters mpp367 (cs367, 8202);
    StopNeutralinoMap.insert(std::make_pair(mp367, mpp367));

    MassPoint mp368 (700, 525);
    StopCrossSection cs368 (0.0044387, 0.000457123);
    MassPointParameters mpp368 (cs368, 12038);
    StopNeutralinoMap.insert(std::make_pair(mp368, mpp368));

    MassPoint mp369 (700, 550);
    StopCrossSection cs369 (0.0044387, 0.000457123);
    MassPointParameters mpp369 (cs369, 8475);
    StopNeutralinoMap.insert(std::make_pair(mp369, mpp369));

    MassPoint mp370 (700, 575);
    StopCrossSection cs370 (0.0044387, 0.000457123);
    MassPointParameters mpp370 (cs370, 8612);
    StopNeutralinoMap.insert(std::make_pair(mp370, mpp370));

    MassPoint mp371 (700, 600);
    StopCrossSection cs371 (0.0044387, 0.000457123);
    MassPointParameters mpp371 (cs371, 11479);
    StopNeutralinoMap.insert(std::make_pair(mp371, mpp371));

    MassPoint mp372 (725, 1);
    StopCrossSection cs372 (0.00370675, 0.000385799);
    MassPointParameters mpp372 (cs372, 9683);
    StopNeutralinoMap.insert(std::make_pair(mp372, mpp372));

    MassPoint mp373 (725, 25);
    StopCrossSection cs373 (0.00370675, 0.000385799);
    MassPointParameters mpp373 (cs373, 11784);
    StopNeutralinoMap.insert(std::make_pair(mp373, mpp373));

    MassPoint mp374 (725, 50);
    StopCrossSection cs374 (0.00370675, 0.000385799);
    MassPointParameters mpp374 (cs374, 6940);
    StopNeutralinoMap.insert(std::make_pair(mp374, mpp374));

    MassPoint mp375 (725, 75);
    StopCrossSection cs375 (0.00370675, 0.000385799);
    MassPointParameters mpp375 (cs375, 13035);
    StopNeutralinoMap.insert(std::make_pair(mp375, mpp375));

    MassPoint mp376 (725, 100);
    StopCrossSection cs376 (0.00370675, 0.000385799);
    MassPointParameters mpp376 (cs376, 11261);
    StopNeutralinoMap.insert(std::make_pair(mp376, mpp376));

    MassPoint mp377 (725, 125);
    StopCrossSection cs377 (0.00370675, 0.000385799);
    MassPointParameters mpp377 (cs377, 10616);
    StopNeutralinoMap.insert(std::make_pair(mp377, mpp377));

    MassPoint mp378 (725, 150);
    StopCrossSection cs378 (0.00370675, 0.000385799);
    MassPointParameters mpp378 (cs378, 8528);
    StopNeutralinoMap.insert(std::make_pair(mp378, mpp378));

    MassPoint mp379 (725, 175);
    StopCrossSection cs379 (0.00370675, 0.000385799);
    MassPointParameters mpp379 (cs379, 11044);
    StopNeutralinoMap.insert(std::make_pair(mp379, mpp379));

    MassPoint mp380 (725, 200);
    StopCrossSection cs380 (0.00370675, 0.000385799);
    MassPointParameters mpp380 (cs380, 10797);
    StopNeutralinoMap.insert(std::make_pair(mp380, mpp380));

    MassPoint mp381 (725, 225);
    StopCrossSection cs381 (0.00370675, 0.000385799);
    MassPointParameters mpp381 (cs381, 8092);
    StopNeutralinoMap.insert(std::make_pair(mp381, mpp381));

    MassPoint mp382 (725, 250);
    StopCrossSection cs382 (0.00370675, 0.000385799);
    MassPointParameters mpp382 (cs382, 10422);
    StopNeutralinoMap.insert(std::make_pair(mp382, mpp382));

    MassPoint mp383 (725, 275);
    StopCrossSection cs383 (0.00370675, 0.000385799);
    MassPointParameters mpp383 (cs383, 9075);
    StopNeutralinoMap.insert(std::make_pair(mp383, mpp383));

    MassPoint mp384 (725, 300);
    StopCrossSection cs384 (0.00370675, 0.000385799);
    MassPointParameters mpp384 (cs384, 8673);
    StopNeutralinoMap.insert(std::make_pair(mp384, mpp384));

    MassPoint mp385 (725, 325);
    StopCrossSection cs385 (0.00370675, 0.000385799);
    MassPointParameters mpp385 (cs385, 8416);
    StopNeutralinoMap.insert(std::make_pair(mp385, mpp385));

    MassPoint mp386 (725, 350);
    StopCrossSection cs386 (0.00370675, 0.000385799);
    MassPointParameters mpp386 (cs386, 11066);
    StopNeutralinoMap.insert(std::make_pair(mp386, mpp386));

    MassPoint mp387 (725, 375);
    StopCrossSection cs387 (0.00370675, 0.000385799);
    MassPointParameters mpp387 (cs387, 11390);
    StopNeutralinoMap.insert(std::make_pair(mp387, mpp387));

    MassPoint mp388 (725, 400);
    StopCrossSection cs388 (0.00370675, 0.000385799);
    MassPointParameters mpp388 (cs388, 11453);
    StopNeutralinoMap.insert(std::make_pair(mp388, mpp388));

    MassPoint mp389 (725, 425);
    StopCrossSection cs389 (0.00370675, 0.000385799);
    MassPointParameters mpp389 (cs389, 10094);
    StopNeutralinoMap.insert(std::make_pair(mp389, mpp389));

    MassPoint mp390 (725, 450);
    StopCrossSection cs390 (0.00370675, 0.000385799);
    MassPointParameters mpp390 (cs390, 7249);
    StopNeutralinoMap.insert(std::make_pair(mp390, mpp390));

    MassPoint mp391 (725, 475);
    StopCrossSection cs391 (0.00370675, 0.000385799);
    MassPointParameters mpp391 (cs391, 8634);
    StopNeutralinoMap.insert(std::make_pair(mp391, mpp391));

    MassPoint mp392 (725, 500);
    StopCrossSection cs392 (0.00370675, 0.000385799);
    MassPointParameters mpp392 (cs392, 12076);
    StopNeutralinoMap.insert(std::make_pair(mp392, mpp392));

    MassPoint mp393 (725, 525);
    StopCrossSection cs393 (0.00370675, 0.000385799);
    MassPointParameters mpp393 (cs393, 10605);
    StopNeutralinoMap.insert(std::make_pair(mp393, mpp393));

    MassPoint mp394 (725, 550);
    StopCrossSection cs394 (0.00370675, 0.000385799);
    MassPointParameters mpp394 (cs394, 9683);
    StopNeutralinoMap.insert(std::make_pair(mp394, mpp394));

    MassPoint mp395 (725, 575);
    StopCrossSection cs395 (0.00370675, 0.000385799);
    MassPointParameters mpp395 (cs395, 10696);
    StopNeutralinoMap.insert(std::make_pair(mp395, mpp395));

    MassPoint mp396 (725, 600);
    StopCrossSection cs396 (0.00370675, 0.000385799);
    MassPointParameters mpp396 (cs396, 10683);
    StopNeutralinoMap.insert(std::make_pair(mp396, mpp396));

    MassPoint mp397 (750, 1);
    StopCrossSection cs397 (0.00310861, 0.000330353);
    MassPointParameters mpp397 (cs397, 11724);
    StopNeutralinoMap.insert(std::make_pair(mp397, mpp397));

    MassPoint mp398 (750, 25);
    StopCrossSection cs398 (0.00310861, 0.000330353);
    MassPointParameters mpp398 (cs398, 11565);
    StopNeutralinoMap.insert(std::make_pair(mp398, mpp398));

    MassPoint mp399 (750, 50);
    StopCrossSection cs399 (0.00310861, 0.000330353);
    MassPointParameters mpp399 (cs399, 11761);
    StopNeutralinoMap.insert(std::make_pair(mp399, mpp399));

    MassPoint mp400 (750, 75);
    StopCrossSection cs400 (0.00310861, 0.000330353);
    MassPointParameters mpp400 (cs400, 11726);
    StopNeutralinoMap.insert(std::make_pair(mp400, mpp400));

    MassPoint mp401 (750, 100);
    StopCrossSection cs401 (0.00310861, 0.000330353);
    MassPointParameters mpp401 (cs401, 8057);
    StopNeutralinoMap.insert(std::make_pair(mp401, mpp401));

    MassPoint mp402 (750, 125);
    StopCrossSection cs402 (0.00310861, 0.000330353);
    MassPointParameters mpp402 (cs402, 11231);
    StopNeutralinoMap.insert(std::make_pair(mp402, mpp402));

    MassPoint mp403 (750, 150);
    StopCrossSection cs403 (0.00310861, 0.000330353);
    MassPointParameters mpp403 (cs403, 9646);
    StopNeutralinoMap.insert(std::make_pair(mp403, mpp403));

    MassPoint mp404 (750, 175);
    StopCrossSection cs404 (0.00310861, 0.000330353);
    MassPointParameters mpp404 (cs404, 11141);
    StopNeutralinoMap.insert(std::make_pair(mp404, mpp404));

    MassPoint mp405 (750, 200);
    StopCrossSection cs405 (0.00310861, 0.000330353);
    MassPointParameters mpp405 (cs405, 10005);
    StopNeutralinoMap.insert(std::make_pair(mp405, mpp405));

    MassPoint mp406 (750, 225);
    StopCrossSection cs406 (0.00310861, 0.000330353);
    MassPointParameters mpp406 (cs406, 9547);
    StopNeutralinoMap.insert(std::make_pair(mp406, mpp406));

    MassPoint mp407 (750, 250);
    StopCrossSection cs407 (0.00310861, 0.000330353);
    MassPointParameters mpp407 (cs407, 10304);
    StopNeutralinoMap.insert(std::make_pair(mp407, mpp407));

    MassPoint mp408 (750, 275);
    StopCrossSection cs408 (0.00310861, 0.000330353);
    MassPointParameters mpp408 (cs408, 10703);
    StopNeutralinoMap.insert(std::make_pair(mp408, mpp408));

    MassPoint mp409 (750, 300);
    StopCrossSection cs409 (0.00310861, 0.000330353);
    MassPointParameters mpp409 (cs409, 11342);
    StopNeutralinoMap.insert(std::make_pair(mp409, mpp409));

    MassPoint mp410 (750, 325);
    StopCrossSection cs410 (0.00310861, 0.000330353);
    MassPointParameters mpp410 (cs410, 8466);
    StopNeutralinoMap.insert(std::make_pair(mp410, mpp410));

    MassPoint mp411 (750, 350);
    StopCrossSection cs411 (0.00310861, 0.000330353);
    MassPointParameters mpp411 (cs411, 12166);
    StopNeutralinoMap.insert(std::make_pair(mp411, mpp411));

    MassPoint mp412 (750, 375);
    StopCrossSection cs412 (0.00310861, 0.000330353);
    MassPointParameters mpp412 (cs412, 9765);
    StopNeutralinoMap.insert(std::make_pair(mp412, mpp412));

    MassPoint mp413 (750, 400);
    StopCrossSection cs413 (0.00310861, 0.000330353);
    MassPointParameters mpp413 (cs413, 10085);
    StopNeutralinoMap.insert(std::make_pair(mp413, mpp413));

    MassPoint mp414 (750, 425);
    StopCrossSection cs414 (0.00310861, 0.000330353);
    MassPointParameters mpp414 (cs414, 10825);
    StopNeutralinoMap.insert(std::make_pair(mp414, mpp414));

    MassPoint mp415 (750, 450);
    StopCrossSection cs415 (0.00310861, 0.000330353);
    MassPointParameters mpp415 (cs415, 10469);
    StopNeutralinoMap.insert(std::make_pair(mp415, mpp415));

    MassPoint mp416 (750, 475);
    StopCrossSection cs416 (0.00310861, 0.000330353);
    MassPointParameters mpp416 (cs416, 10614);
    StopNeutralinoMap.insert(std::make_pair(mp416, mpp416));

    MassPoint mp417 (750, 500);
    StopCrossSection cs417 (0.00310861, 0.000330353);
    MassPointParameters mpp417 (cs417, 10561);
    StopNeutralinoMap.insert(std::make_pair(mp417, mpp417));

    MassPoint mp418 (750, 525);
    StopCrossSection cs418 (0.00310861, 0.000330353);
    MassPointParameters mpp418 (cs418, 10217);
    StopNeutralinoMap.insert(std::make_pair(mp418, mpp418));

    MassPoint mp419 (750, 550);
    StopCrossSection cs419 (0.00310861, 0.000330353);
    MassPointParameters mpp419 (cs419, 8478);
    StopNeutralinoMap.insert(std::make_pair(mp419, mpp419));

    MassPoint mp420 (750, 575);
    StopCrossSection cs420 (0.00310861, 0.000330353);
    MassPointParameters mpp420 (cs420, 8197);
    StopNeutralinoMap.insert(std::make_pair(mp420, mpp420));

    MassPoint mp421 (750, 600);
    StopCrossSection cs421 (0.00310861, 0.000330353);
    MassPointParameters mpp421 (cs421, 10734);
    StopNeutralinoMap.insert(std::make_pair(mp421, mpp421));

    MassPoint mp422 (775, 1);
    StopCrossSection cs422 (0.00261656, 0.000283139);
    MassPointParameters mpp422 (cs422, 10233);
    StopNeutralinoMap.insert(std::make_pair(mp422, mpp422));

    MassPoint mp423 (775, 25);
    StopCrossSection cs423 (0.00261656, 0.000283139);
    MassPointParameters mpp423 (cs423, 12032);
    StopNeutralinoMap.insert(std::make_pair(mp423, mpp423));

    MassPoint mp424 (775, 50);
    StopCrossSection cs424 (0.00261656, 0.000283139);
    MassPointParameters mpp424 (cs424, 10309);
    StopNeutralinoMap.insert(std::make_pair(mp424, mpp424));

    MassPoint mp425 (775, 75);
    StopCrossSection cs425 (0.00261656, 0.000283139);
    MassPointParameters mpp425 (cs425, 10780);
    StopNeutralinoMap.insert(std::make_pair(mp425, mpp425));

    MassPoint mp426 (775, 100);
    StopCrossSection cs426 (0.00261656, 0.000283139);
    MassPointParameters mpp426 (cs426, 8577);
    StopNeutralinoMap.insert(std::make_pair(mp426, mpp426));

    MassPoint mp427 (775, 125);
    StopCrossSection cs427 (0.00261656, 0.000283139);
    MassPointParameters mpp427 (cs427, 11512);
    StopNeutralinoMap.insert(std::make_pair(mp427, mpp427));

    MassPoint mp428 (775, 150);
    StopCrossSection cs428 (0.00261656, 0.000283139);
    MassPointParameters mpp428 (cs428, 9715);
    StopNeutralinoMap.insert(std::make_pair(mp428, mpp428));

    MassPoint mp429 (775, 175);
    StopCrossSection cs429 (0.00261656, 0.000283139);
    MassPointParameters mpp429 (cs429, 10039);
    StopNeutralinoMap.insert(std::make_pair(mp429, mpp429));

    MassPoint mp430 (775, 200);
    StopCrossSection cs430 (0.00261656, 0.000283139);
    MassPointParameters mpp430 (cs430, 9076);
    StopNeutralinoMap.insert(std::make_pair(mp430, mpp430));

    MassPoint mp431 (775, 225);
    StopCrossSection cs431 (0.00261656, 0.000283139);
    MassPointParameters mpp431 (cs431, 11571);
    StopNeutralinoMap.insert(std::make_pair(mp431, mpp431));

    MassPoint mp432 (775, 250);
    StopCrossSection cs432 (0.00261656, 0.000283139);
    MassPointParameters mpp432 (cs432, 8055);
    StopNeutralinoMap.insert(std::make_pair(mp432, mpp432));

    MassPoint mp433 (775, 275);
    StopCrossSection cs433 (0.00261656, 0.000283139);
    MassPointParameters mpp433 (cs433, 9001);
    StopNeutralinoMap.insert(std::make_pair(mp433, mpp433));

    MassPoint mp434 (775, 300);
    StopCrossSection cs434 (0.00261656, 0.000283139);
    MassPointParameters mpp434 (cs434, 11357);
    StopNeutralinoMap.insert(std::make_pair(mp434, mpp434));

    MassPoint mp435 (775, 325);
    StopCrossSection cs435 (0.00261656, 0.000283139);
    MassPointParameters mpp435 (cs435, 11478);
    StopNeutralinoMap.insert(std::make_pair(mp435, mpp435));

    MassPoint mp436 (775, 350);
    StopCrossSection cs436 (0.00261656, 0.000283139);
    MassPointParameters mpp436 (cs436, 10798);
    StopNeutralinoMap.insert(std::make_pair(mp436, mpp436));

    MassPoint mp437 (775, 375);
    StopCrossSection cs437 (0.00261656, 0.000283139);
    MassPointParameters mpp437 (cs437, 10849);
    StopNeutralinoMap.insert(std::make_pair(mp437, mpp437));

    MassPoint mp438 (775, 400);
    StopCrossSection cs438 (0.00261656, 0.000283139);
    MassPointParameters mpp438 (cs438, 11336);
    StopNeutralinoMap.insert(std::make_pair(mp438, mpp438));

    MassPoint mp439 (775, 425);
    StopCrossSection cs439 (0.00261656, 0.000283139);
    MassPointParameters mpp439 (cs439, 11183);
    StopNeutralinoMap.insert(std::make_pair(mp439, mpp439));

    MassPoint mp440 (775, 450);
    StopCrossSection cs440 (0.00261656, 0.000283139);
    MassPointParameters mpp440 (cs440, 11130);
    StopNeutralinoMap.insert(std::make_pair(mp440, mpp440));

    MassPoint mp441 (775, 475);
    StopCrossSection cs441 (0.00261656, 0.000283139);
    MassPointParameters mpp441 (cs441, 8683);
    StopNeutralinoMap.insert(std::make_pair(mp441, mpp441));

    MassPoint mp442 (775, 500);
    StopCrossSection cs442 (0.00261656, 0.000283139);
    MassPointParameters mpp442 (cs442, 10186);
    StopNeutralinoMap.insert(std::make_pair(mp442, mpp442));

    MassPoint mp443 (775, 525);
    StopCrossSection cs443 (0.00261656, 0.000283139);
    MassPointParameters mpp443 (cs443, 8313);
    StopNeutralinoMap.insert(std::make_pair(mp443, mpp443));

    MassPoint mp444 (775, 550);
    StopCrossSection cs444 (0.00261656, 0.000283139);
    MassPointParameters mpp444 (cs444, 8554);
    StopNeutralinoMap.insert(std::make_pair(mp444, mpp444));

    MassPoint mp445 (775, 575);
    StopCrossSection cs445 (0.00261656, 0.000283139);
    MassPointParameters mpp445 (cs445, 12557);
    StopNeutralinoMap.insert(std::make_pair(mp445, mpp445));

    MassPoint mp446 (775, 600);
    StopCrossSection cs446 (0.00261656, 0.000283139);
    MassPointParameters mpp446 (cs446, 8718);
    StopNeutralinoMap.insert(std::make_pair(mp446, mpp446));

    MassPoint mp447 (800, 1);
    StopCrossSection cs447 (0.00221197, 0.000245196);
    MassPointParameters mpp447 (cs447, 7401);
    StopNeutralinoMap.insert(std::make_pair(mp447, mpp447));

    MassPoint mp448 (800, 25);
    StopCrossSection cs448 (0.00221197, 0.000245196);
    MassPointParameters mpp448 (cs448, 10303);
    StopNeutralinoMap.insert(std::make_pair(mp448, mpp448));

    MassPoint mp449 (800, 50);
    StopCrossSection cs449 (0.00221197, 0.000245196);
    MassPointParameters mpp449 (cs449, 10522);
    StopNeutralinoMap.insert(std::make_pair(mp449, mpp449));

    MassPoint mp450 (800, 75);
    StopCrossSection cs450 (0.00221197, 0.000245196);
    MassPointParameters mpp450 (cs450, 11842);
    StopNeutralinoMap.insert(std::make_pair(mp450, mpp450));

    MassPoint mp451 (800, 100);
    StopCrossSection cs451 (0.00221197, 0.000245196);
    MassPointParameters mpp451 (cs451, 11364);
    StopNeutralinoMap.insert(std::make_pair(mp451, mpp451));

    MassPoint mp452 (800, 125);
    StopCrossSection cs452 (0.00221197, 0.000245196);
    MassPointParameters mpp452 (cs452, 7996);
    StopNeutralinoMap.insert(std::make_pair(mp452, mpp452));

    MassPoint mp453 (800, 150);
    StopCrossSection cs453 (0.00221197, 0.000245196);
    MassPointParameters mpp453 (cs453, 8674);
    StopNeutralinoMap.insert(std::make_pair(mp453, mpp453));

    MassPoint mp454 (800, 175);
    StopCrossSection cs454 (0.00221197, 0.000245196);
    MassPointParameters mpp454 (cs454, 9150);
    StopNeutralinoMap.insert(std::make_pair(mp454, mpp454));

    MassPoint mp455 (800, 200);
    StopCrossSection cs455 (0.00221197, 0.000245196);
    MassPointParameters mpp455 (cs455, 10402);
    StopNeutralinoMap.insert(std::make_pair(mp455, mpp455));

    MassPoint mp456 (800, 225);
    StopCrossSection cs456 (0.00221197, 0.000245196);
    MassPointParameters mpp456 (cs456, 9645);
    StopNeutralinoMap.insert(std::make_pair(mp456, mpp456));

    MassPoint mp457 (800, 250);
    StopCrossSection cs457 (0.00221197, 0.000245196);
    MassPointParameters mpp457 (cs457, 10344);
    StopNeutralinoMap.insert(std::make_pair(mp457, mpp457));

    MassPoint mp458 (800, 275);
    StopCrossSection cs458 (0.00221197, 0.000245196);
    MassPointParameters mpp458 (cs458, 7010);
    StopNeutralinoMap.insert(std::make_pair(mp458, mpp458));

    MassPoint mp459 (800, 300);
    StopCrossSection cs459 (0.00221197, 0.000245196);
    MassPointParameters mpp459 (cs459, 10259);
    StopNeutralinoMap.insert(std::make_pair(mp459, mpp459));

    MassPoint mp460 (800, 325);
    StopCrossSection cs460 (0.00221197, 0.000245196);
    MassPointParameters mpp460 (cs460, 11380);
    StopNeutralinoMap.insert(std::make_pair(mp460, mpp460));

    MassPoint mp461 (800, 350);
    StopCrossSection cs461 (0.00221197, 0.000245196);
    MassPointParameters mpp461 (cs461, 10624);
    StopNeutralinoMap.insert(std::make_pair(mp461, mpp461));

    MassPoint mp462 (800, 375);
    StopCrossSection cs462 (0.00221197, 0.000245196);
    MassPointParameters mpp462 (cs462, 12881);
    StopNeutralinoMap.insert(std::make_pair(mp462, mpp462));

    MassPoint mp463 (800, 400);
    StopCrossSection cs463 (0.00221197, 0.000245196);
    MassPointParameters mpp463 (cs463, 9584);
    StopNeutralinoMap.insert(std::make_pair(mp463, mpp463));

    MassPoint mp464 (800, 425);
    StopCrossSection cs464 (0.00221197, 0.000245196);
    MassPointParameters mpp464 (cs464, 9235);
    StopNeutralinoMap.insert(std::make_pair(mp464, mpp464));

    MassPoint mp465 (800, 450);
    StopCrossSection cs465 (0.00221197, 0.000245196);
    MassPointParameters mpp465 (cs465, 12168);
    StopNeutralinoMap.insert(std::make_pair(mp465, mpp465));

    MassPoint mp466 (800, 475);
    StopCrossSection cs466 (0.00221197, 0.000245196);
    MassPointParameters mpp466 (cs466, 8508);
    StopNeutralinoMap.insert(std::make_pair(mp466, mpp466));

    MassPoint mp467 (800, 500);
    StopCrossSection cs467 (0.00221197, 0.000245196);
    MassPointParameters mpp467 (cs467, 10758);
    StopNeutralinoMap.insert(std::make_pair(mp467, mpp467));

    MassPoint mp468 (800, 525);
    StopCrossSection cs468 (0.00221197, 0.000245196);
    MassPointParameters mpp468 (cs468, 11118);
    StopNeutralinoMap.insert(std::make_pair(mp468, mpp468));

    MassPoint mp469 (800, 550);
    StopCrossSection cs469 (0.00221197, 0.000245196);
    MassPointParameters mpp469 (cs469, 9167);
    StopNeutralinoMap.insert(std::make_pair(mp469, mpp469));

    MassPoint mp470 (800, 575);
    StopCrossSection cs470 (0.00221197, 0.000245196);
    MassPointParameters mpp470 (cs470, 9986);
    StopNeutralinoMap.insert(std::make_pair(mp470, mpp470));

    MassPoint mp471 (800, 600);
    StopCrossSection cs471 (0.00221197, 0.000245196);
    MassPointParameters mpp471 (cs471, 10662);
    StopNeutralinoMap.insert(std::make_pair(mp471, mpp471));

    MassPoint mpext0 (825, 1);
    StopCrossSection csext0 (0.00186142, 0.000201762);
    MassPointParameters mppext0 (csext0, 8916);
    StopNeutralinoMap.insert(std::make_pair(mpext0, mppext0));

    MassPoint mpext1 (825, 25);
    StopCrossSection csext1 (0.00186142, 0.000201762);
    MassPointParameters mppext1 (csext1, 9070);
    StopNeutralinoMap.insert(std::make_pair(mpext1, mppext1));

    MassPoint mpext2 (825, 50);
    StopCrossSection csext2 (0.00186142, 0.000201762);
    MassPointParameters mppext2 (csext2, 13241);
    StopNeutralinoMap.insert(std::make_pair(mpext2, mppext2));

    MassPoint mpext3 (825, 75);
    StopCrossSection csext3 (0.00186142, 0.000201762);
    MassPointParameters mppext3 (csext3, 8117);
    StopNeutralinoMap.insert(std::make_pair(mpext3, mppext3));

    MassPoint mpext4 (825, 100);
    StopCrossSection csext4 (0.00186142, 0.000201762);
    MassPointParameters mppext4 (csext4, 11753);
    StopNeutralinoMap.insert(std::make_pair(mpext4, mppext4));

    MassPoint mpext5 (825, 125);
    StopCrossSection csext5 (0.00186142, 0.000201762);
    MassPointParameters mppext5 (csext5, 9669);
    StopNeutralinoMap.insert(std::make_pair(mpext5, mppext5));

    MassPoint mpext6 (825, 150);
    StopCrossSection csext6 (0.00186142, 0.000201762);
    MassPointParameters mppext6 (csext6, 11899);
    StopNeutralinoMap.insert(std::make_pair(mpext6, mppext6));

    MassPoint mpext7 (825, 175);
    StopCrossSection csext7 (0.00186142, 0.000201762);
    MassPointParameters mppext7 (csext7, 9897);
    StopNeutralinoMap.insert(std::make_pair(mpext7, mppext7));

    MassPoint mpext8 (825, 200);
    StopCrossSection csext8 (0.00186142, 0.000201762);
    MassPointParameters mppext8 (csext8, 9368);
    StopNeutralinoMap.insert(std::make_pair(mpext8, mppext8));

    MassPoint mpext9 (825, 225);
    StopCrossSection csext9 (0.00186142, 0.000201762);
    MassPointParameters mppext9 (csext9, 11492);
    StopNeutralinoMap.insert(std::make_pair(mpext9, mppext9));

    MassPoint mpext10 (825, 250);
    StopCrossSection csext10 (0.00186142, 0.000201762);
    MassPointParameters mppext10 (csext10, 11350);
    StopNeutralinoMap.insert(std::make_pair(mpext10, mppext10));

    MassPoint mpext11 (825, 275);
    StopCrossSection csext11 (0.00186142, 0.000201762);
    MassPointParameters mppext11 (csext11, 8913);
    StopNeutralinoMap.insert(std::make_pair(mpext11, mppext11));

    MassPoint mpext12 (825, 300);
    StopCrossSection csext12 (0.00186142, 0.000201762);
    MassPointParameters mppext12 (csext12, 11326);
    StopNeutralinoMap.insert(std::make_pair(mpext12, mppext12));

    MassPoint mpext13 (825, 325);
    StopCrossSection csext13 (0.00186142, 0.000201762);
    MassPointParameters mppext13 (csext13, 12785);
    StopNeutralinoMap.insert(std::make_pair(mpext13, mppext13));

    MassPoint mpext14 (825, 350);
    StopCrossSection csext14 (0.00186142, 0.000201762);
    MassPointParameters mppext14 (csext14, 10873);
    StopNeutralinoMap.insert(std::make_pair(mpext14, mppext14));

    MassPoint mpext15 (825, 375);
    StopCrossSection csext15 (0.00186142, 0.000201762);
    MassPointParameters mppext15 (csext15, 11987);
    StopNeutralinoMap.insert(std::make_pair(mpext15, mppext15));

    MassPoint mpext16 (825, 400);
    StopCrossSection csext16 (0.00186142, 0.000201762);
    MassPointParameters mppext16 (csext16, 11298);
    StopNeutralinoMap.insert(std::make_pair(mpext16, mppext16));

    MassPoint mpext17 (825, 425);
    StopCrossSection csext17 (0.00186142, 0.000201762);
    MassPointParameters mppext17 (csext17, 8799);
    StopNeutralinoMap.insert(std::make_pair(mpext17, mppext17));

    MassPoint mpext18 (825, 450);
    StopCrossSection csext18 (0.00186142, 0.000201762);
    MassPointParameters mppext18 (csext18, 8809);
    StopNeutralinoMap.insert(std::make_pair(mpext18, mppext18));

    MassPoint mpext19 (825, 475);
    StopCrossSection csext19 (0.00186142, 0.000201762);
    MassPointParameters mppext19 (csext19, 9519);
    StopNeutralinoMap.insert(std::make_pair(mpext19, mppext19));

    MassPoint mpext20 (825, 500);
    StopCrossSection csext20 (0.00186142, 0.000201762);
    MassPointParameters mppext20 (csext20, 8177);
    StopNeutralinoMap.insert(std::make_pair(mpext20, mppext20));

    MassPoint mpext21 (825, 525);
    StopCrossSection csext21 (0.00186142, 0.000201762);
    MassPointParameters mppext21 (csext21, 10502);
    StopNeutralinoMap.insert(std::make_pair(mpext21, mppext21));

    MassPoint mpext22 (825, 550);
    StopCrossSection csext22 (0.00186142, 0.000201762);
    MassPointParameters mppext22 (csext22, 8388);
    StopNeutralinoMap.insert(std::make_pair(mpext22, mppext22));

    MassPoint mpext23 (825, 575);
    StopCrossSection csext23 (0.00186142, 0.000201762);
    MassPointParameters mppext23 (csext23, 9804);
    StopNeutralinoMap.insert(std::make_pair(mpext23, mppext23));

    MassPoint mpext24 (825, 600);
    StopCrossSection csext24 (0.00186142, 0.000201762);
    MassPointParameters mppext24 (csext24, 8710);
    StopNeutralinoMap.insert(std::make_pair(mpext24, mppext24));

    MassPoint mpext25 (850, 1);
    StopCrossSection csext25 (0.00158356, 0.000177806);
    MassPointParameters mppext25 (csext25, 10942);
    StopNeutralinoMap.insert(std::make_pair(mpext25, mppext25));

    MassPoint mpext26 (850, 25);
    StopCrossSection csext26 (0.00158356, 0.000177806);
    MassPointParameters mppext26 (csext26, 8745);
    StopNeutralinoMap.insert(std::make_pair(mpext26, mppext26));

    MassPoint mpext27 (850, 50);
    StopCrossSection csext27 (0.00158356, 0.000177806);
    MassPointParameters mppext27 (csext27, 8888);
    StopNeutralinoMap.insert(std::make_pair(mpext27, mppext27));

    MassPoint mpext28 (850, 75);
    StopCrossSection csext28 (0.00158356, 0.000177806);
    MassPointParameters mppext28 (csext28, 8299);
    StopNeutralinoMap.insert(std::make_pair(mpext28, mppext28));

    MassPoint mpext29 (850, 100);
    StopCrossSection csext29 (0.00158356, 0.000177806);
    MassPointParameters mppext29 (csext29, 12612);
    StopNeutralinoMap.insert(std::make_pair(mpext29, mppext29));

    MassPoint mpext30 (850, 125);
    StopCrossSection csext30 (0.00158356, 0.000177806);
    MassPointParameters mppext30 (csext30, 10722);
    StopNeutralinoMap.insert(std::make_pair(mpext30, mppext30));

    MassPoint mpext31 (850, 150);
    StopCrossSection csext31 (0.00158356, 0.000177806);
    MassPointParameters mppext31 (csext31, 10393);
    StopNeutralinoMap.insert(std::make_pair(mpext31, mppext31));

    MassPoint mpext32 (850, 175);
    StopCrossSection csext32 (0.00158356, 0.000177806);
    MassPointParameters mppext32 (csext32, 9052);
    StopNeutralinoMap.insert(std::make_pair(mpext32, mppext32));

    MassPoint mpext33 (850, 200);
    StopCrossSection csext33 (0.00158356, 0.000177806);
    MassPointParameters mppext33 (csext33, 9551);
    StopNeutralinoMap.insert(std::make_pair(mpext33, mppext33));

    MassPoint mpext34 (850, 225);
    StopCrossSection csext34 (0.00158356, 0.000177806);
    MassPointParameters mppext34 (csext34, 10612);
    StopNeutralinoMap.insert(std::make_pair(mpext34, mppext34));

    MassPoint mpext35 (850, 250);
    StopCrossSection csext35 (0.00158356, 0.000177806);
    MassPointParameters mppext35 (csext35, 9933);
    StopNeutralinoMap.insert(std::make_pair(mpext35, mppext35));

    MassPoint mpext36 (850, 275);
    StopCrossSection csext36 (0.00158356, 0.000177806);
    MassPointParameters mppext36 (csext36, 11430);
    StopNeutralinoMap.insert(std::make_pair(mpext36, mppext36));

    MassPoint mpext37 (850, 300);
    StopCrossSection csext37 (0.00158356, 0.000177806);
    MassPointParameters mppext37 (csext37, 9245);
    StopNeutralinoMap.insert(std::make_pair(mpext37, mppext37));

    MassPoint mpext38 (850, 325);
    StopCrossSection csext38 (0.00158356, 0.000177806);
    MassPointParameters mppext38 (csext38, 11135);
    StopNeutralinoMap.insert(std::make_pair(mpext38, mppext38));

    MassPoint mpext39 (850, 350);
    StopCrossSection csext39 (0.00158356, 0.000177806);
    MassPointParameters mppext39 (csext39, 9137);
    StopNeutralinoMap.insert(std::make_pair(mpext39, mppext39));

    MassPoint mpext40 (850, 375);
    StopCrossSection csext40 (0.00158356, 0.000177806);
    MassPointParameters mppext40 (csext40, 10509);
    StopNeutralinoMap.insert(std::make_pair(mpext40, mppext40));

    MassPoint mpext41 (850, 400);
    StopCrossSection csext41 (0.00158356, 0.000177806);
    MassPointParameters mppext41 (csext41, 9920);
    StopNeutralinoMap.insert(std::make_pair(mpext41, mppext41));

    MassPoint mpext42 (850, 425);
    StopCrossSection csext42 (0.00158356, 0.000177806);
    MassPointParameters mppext42 (csext42, 10869);
    StopNeutralinoMap.insert(std::make_pair(mpext42, mppext42));

    MassPoint mpext43 (850, 450);
    StopCrossSection csext43 (0.00158356, 0.000177806);
    MassPointParameters mppext43 (csext43, 13727);
    StopNeutralinoMap.insert(std::make_pair(mpext43, mppext43));

    MassPoint mpext44 (850, 475);
    StopCrossSection csext44 (0.00158356, 0.000177806);
    MassPointParameters mppext44 (csext44, 11411);
    StopNeutralinoMap.insert(std::make_pair(mpext44, mppext44));

    MassPoint mpext45 (850, 500);
    StopCrossSection csext45 (0.00158356, 0.000177806);
    MassPointParameters mppext45 (csext45, 10155);
    StopNeutralinoMap.insert(std::make_pair(mpext45, mppext45));

    MassPoint mpext46 (850, 525);
    StopCrossSection csext46 (0.00158356, 0.000177806);
    MassPointParameters mppext46 (csext46, 9237);
    StopNeutralinoMap.insert(std::make_pair(mpext46, mppext46));

    MassPoint mpext47 (850, 550);
    StopCrossSection csext47 (0.00158356, 0.000177806);
    MassPointParameters mppext47 (csext47, 9936);
    StopNeutralinoMap.insert(std::make_pair(mpext47, mppext47));

    MassPoint mpext48 (850, 575);
    StopCrossSection csext48 (0.00158356, 0.000177806);
    MassPointParameters mppext48 (csext48, 9449);
    StopNeutralinoMap.insert(std::make_pair(mpext48, mppext48));

    MassPoint mpext49 (850, 600);
    StopCrossSection csext49 (0.00158356, 0.000177806);
    MassPointParameters mppext49 (csext49, 9047);
    StopNeutralinoMap.insert(std::make_pair(mpext49, mppext49));

    MassPoint mpext50 (875, 1);
    StopCrossSection csext50 (0.00134699, 0.000150075);
    MassPointParameters mppext50 (csext50, 9791);
    StopNeutralinoMap.insert(std::make_pair(mpext50, mppext50));

    MassPoint mpext51 (875, 25);
    StopCrossSection csext51 (0.00134699, 0.000150075);
    MassPointParameters mppext51 (csext51, 12365);
    StopNeutralinoMap.insert(std::make_pair(mpext51, mppext51));

    MassPoint mpext52 (875, 50);
    StopCrossSection csext52 (0.00134699, 0.000150075);
    MassPointParameters mppext52 (csext52, 11819);
    StopNeutralinoMap.insert(std::make_pair(mpext52, mppext52));

    MassPoint mpext53 (875, 75);
    StopCrossSection csext53 (0.00134699, 0.000150075);
    MassPointParameters mppext53 (csext53, 8680);
    StopNeutralinoMap.insert(std::make_pair(mpext53, mppext53));

    MassPoint mpext54 (875, 100);
    StopCrossSection csext54 (0.00134699, 0.000150075);
    MassPointParameters mppext54 (csext54, 9542);
    StopNeutralinoMap.insert(std::make_pair(mpext54, mppext54));

    MassPoint mpext55 (875, 125);
    StopCrossSection csext55 (0.00134699, 0.000150075);
    MassPointParameters mppext55 (csext55, 10652);
    StopNeutralinoMap.insert(std::make_pair(mpext55, mppext55));

    MassPoint mpext56 (875, 150);
    StopCrossSection csext56 (0.00134699, 0.000150075);
    MassPointParameters mppext56 (csext56, 6869);
    StopNeutralinoMap.insert(std::make_pair(mpext56, mppext56));

    MassPoint mpext57 (875, 175);
    StopCrossSection csext57 (0.00134699, 0.000150075);
    MassPointParameters mppext57 (csext57, 10580);
    StopNeutralinoMap.insert(std::make_pair(mpext57, mppext57));

    MassPoint mpext58 (875, 200);
    StopCrossSection csext58 (0.00134699, 0.000150075);
    MassPointParameters mppext58 (csext58, 11055);
    StopNeutralinoMap.insert(std::make_pair(mpext58, mppext58));

    MassPoint mpext59 (875, 225);
    StopCrossSection csext59 (0.00134699, 0.000150075);
    MassPointParameters mppext59 (csext59, 10806);
    StopNeutralinoMap.insert(std::make_pair(mpext59, mppext59));

    MassPoint mpext60 (875, 250);
    StopCrossSection csext60 (0.00134699, 0.000150075);
    MassPointParameters mppext60 (csext60, 10163);
    StopNeutralinoMap.insert(std::make_pair(mpext60, mppext60));

    MassPoint mpext61 (875, 275);
    StopCrossSection csext61 (0.00134699, 0.000150075);
    MassPointParameters mppext61 (csext61, 10809);
    StopNeutralinoMap.insert(std::make_pair(mpext61, mppext61));

    MassPoint mpext62 (875, 300);
    StopCrossSection csext62 (0.00134699, 0.000150075);
    MassPointParameters mppext62 (csext62, 11076);
    StopNeutralinoMap.insert(std::make_pair(mpext62, mppext62));

    MassPoint mpext63 (875, 325);
    StopCrossSection csext63 (0.00134699, 0.000150075);
    MassPointParameters mppext63 (csext63, 12265);
    StopNeutralinoMap.insert(std::make_pair(mpext63, mppext63));

    MassPoint mpext64 (875, 350);
    StopCrossSection csext64 (0.00134699, 0.000150075);
    MassPointParameters mppext64 (csext64, 10143);
    StopNeutralinoMap.insert(std::make_pair(mpext64, mppext64));

    MassPoint mpext65 (875, 375);
    StopCrossSection csext65 (0.00134699, 0.000150075);
    MassPointParameters mppext65 (csext65, 9047);
    StopNeutralinoMap.insert(std::make_pair(mpext65, mppext65));

    MassPoint mpext66 (875, 400);
    StopCrossSection csext66 (0.00134699, 0.000150075);
    MassPointParameters mppext66 (csext66, 11825);
    StopNeutralinoMap.insert(std::make_pair(mpext66, mppext66));

    MassPoint mpext67 (875, 425);
    StopCrossSection csext67 (0.00134699, 0.000150075);
    MassPointParameters mppext67 (csext67, 11518);
    StopNeutralinoMap.insert(std::make_pair(mpext67, mppext67));

    MassPoint mpext68 (875, 450);
    StopCrossSection csext68 (0.00134699, 0.000150075);
    MassPointParameters mppext68 (csext68, 8572);
    StopNeutralinoMap.insert(std::make_pair(mpext68, mppext68));

    MassPoint mpext69 (875, 475);
    StopCrossSection csext69 (0.00134699, 0.000150075);
    MassPointParameters mppext69 (csext69, 9515);
    StopNeutralinoMap.insert(std::make_pair(mpext69, mppext69));

    MassPoint mpext70 (875, 500);
    StopCrossSection csext70 (0.00134699, 0.000150075);
    MassPointParameters mppext70 (csext70, 9687);
    StopNeutralinoMap.insert(std::make_pair(mpext70, mppext70));

    MassPoint mpext71 (875, 525);
    StopCrossSection csext71 (0.00134699, 0.000150075);
    MassPointParameters mppext71 (csext71, 10418);
    StopNeutralinoMap.insert(std::make_pair(mpext71, mppext71));

    MassPoint mpext72 (875, 550);
    StopCrossSection csext72 (0.00134699, 0.000150075);
    MassPointParameters mppext72 (csext72, 7877);
    StopNeutralinoMap.insert(std::make_pair(mpext72, mppext72));

    MassPoint mpext73 (875, 575);
    StopCrossSection csext73 (0.00134699, 0.000150075);
    MassPointParameters mppext73 (csext73, 10652);
    StopNeutralinoMap.insert(std::make_pair(mpext73, mppext73));

    MassPoint mpext74 (875, 600);
    StopCrossSection csext74 (0.00134699, 0.000150075);
    MassPointParameters mppext74 (csext74, 9831);
    StopNeutralinoMap.insert(std::make_pair(mpext74, mppext74));

    MassPoint mpext75 (900, 1);
    StopCrossSection csext75 (0.00115301, 0.000135822);
    MassPointParameters mppext75 (csext75, 6790);
    StopNeutralinoMap.insert(std::make_pair(mpext75, mppext75));

    MassPoint mpext76 (900, 25);
    StopCrossSection csext76 (0.00115301, 0.000135822);
    MassPointParameters mppext76 (csext76, 10766);
    StopNeutralinoMap.insert(std::make_pair(mpext76, mppext76));

    MassPoint mpext77 (900, 50);
    StopCrossSection csext77 (0.00115301, 0.000135822);
    MassPointParameters mppext77 (csext77, 11396);
    StopNeutralinoMap.insert(std::make_pair(mpext77, mppext77));

    MassPoint mpext78 (900, 75);
    StopCrossSection csext78 (0.00115301, 0.000135822);
    MassPointParameters mppext78 (csext78, 10038);
    StopNeutralinoMap.insert(std::make_pair(mpext78, mppext78));

    MassPoint mpext79 (900, 100);
    StopCrossSection csext79 (0.00115301, 0.000135822);
    MassPointParameters mppext79 (csext79, 9856);
    StopNeutralinoMap.insert(std::make_pair(mpext79, mppext79));

    MassPoint mpext80 (900, 125);
    StopCrossSection csext80 (0.00115301, 0.000135822);
    MassPointParameters mppext80 (csext80, 7403);
    StopNeutralinoMap.insert(std::make_pair(mpext80, mppext80));

    MassPoint mpext81 (900, 150);
    StopCrossSection csext81 (0.00115301, 0.000135822);
    MassPointParameters mppext81 (csext81, 9749);
    StopNeutralinoMap.insert(std::make_pair(mpext81, mppext81));

    MassPoint mpext82 (900, 175);
    StopCrossSection csext82 (0.00115301, 0.000135822);
    MassPointParameters mppext82 (csext82, 8267);
    StopNeutralinoMap.insert(std::make_pair(mpext82, mppext82));

    MassPoint mpext83 (900, 200);
    StopCrossSection csext83 (0.00115301, 0.000135822);
    MassPointParameters mppext83 (csext83, 11446);
    StopNeutralinoMap.insert(std::make_pair(mpext83, mppext83));

    MassPoint mpext84 (900, 225);
    StopCrossSection csext84 (0.00115301, 0.000135822);
    MassPointParameters mppext84 (csext84, 10870);
    StopNeutralinoMap.insert(std::make_pair(mpext84, mppext84));

    MassPoint mpext85 (900, 250);
    StopCrossSection csext85 (0.00115301, 0.000135822);
    MassPointParameters mppext85 (csext85, 7312);
    StopNeutralinoMap.insert(std::make_pair(mpext85, mppext85));

    MassPoint mpext86 (900, 275);
    StopCrossSection csext86 (0.00115301, 0.000135822);
    MassPointParameters mppext86 (csext86, 8646);
    StopNeutralinoMap.insert(std::make_pair(mpext86, mppext86));

    MassPoint mpext87 (900, 300);
    StopCrossSection csext87 (0.00115301, 0.000135822);
    MassPointParameters mppext87 (csext87, 10252);
    StopNeutralinoMap.insert(std::make_pair(mpext87, mppext87));

    MassPoint mpext88 (900, 325);
    StopCrossSection csext88 (0.00115301, 0.000135822);
    MassPointParameters mppext88 (csext88, 8515);
    StopNeutralinoMap.insert(std::make_pair(mpext88, mppext88));

    MassPoint mpext89 (900, 350);
    StopCrossSection csext89 (0.00115301, 0.000135822);
    MassPointParameters mppext89 (csext89, 9660);
    StopNeutralinoMap.insert(std::make_pair(mpext89, mppext89));

    MassPoint mpext90 (900, 375);
    StopCrossSection csext90 (0.00115301, 0.000135822);
    MassPointParameters mppext90 (csext90, 8347);
    StopNeutralinoMap.insert(std::make_pair(mpext90, mppext90));

    MassPoint mpext91 (900, 400);
    StopCrossSection csext91 (0.00115301, 0.000135822);
    MassPointParameters mppext91 (csext91, 7756);
    StopNeutralinoMap.insert(std::make_pair(mpext91, mppext91));

    MassPoint mpext92 (900, 425);
    StopCrossSection csext92 (0.00115301, 0.000135822);
    MassPointParameters mppext92 (csext92, 10342);
    StopNeutralinoMap.insert(std::make_pair(mpext92, mppext92));

    MassPoint mpext93 (900, 450);
    StopCrossSection csext93 (0.00115301, 0.000135822);
    MassPointParameters mppext93 (csext93, 10366);
    StopNeutralinoMap.insert(std::make_pair(mpext93, mppext93));

    MassPoint mpext94 (900, 475);
    StopCrossSection csext94 (0.00115301, 0.000135822);
    MassPointParameters mppext94 (csext94, 9846);
    StopNeutralinoMap.insert(std::make_pair(mpext94, mppext94));

    MassPoint mpext95 (900, 500);
    StopCrossSection csext95 (0.00115301, 0.000135822);
    MassPointParameters mppext95 (csext95, 8284);
    StopNeutralinoMap.insert(std::make_pair(mpext95, mppext95));

    MassPoint mpext96 (900, 525);
    StopCrossSection csext96 (0.00115301, 0.000135822);
    MassPointParameters mppext96 (csext96, 8463);
    StopNeutralinoMap.insert(std::make_pair(mpext96, mppext96));

    MassPoint mpext97 (900, 550);
    StopCrossSection csext97 (0.00115301, 0.000135822);
    MassPointParameters mppext97 (csext97, 9575);
    StopNeutralinoMap.insert(std::make_pair(mpext97, mppext97));

    MassPoint mpext98 (900, 575);
    StopCrossSection csext98 (0.00115301, 0.000135822);
    MassPointParameters mppext98 (csext98, 10610);
    StopNeutralinoMap.insert(std::make_pair(mpext98, mppext98));

    MassPoint mpext99 (900, 600);
    StopCrossSection csext99 (0.00115301, 0.000135822);
    MassPointParameters mppext99 (csext99, 11412);
    StopNeutralinoMap.insert(std::make_pair(mpext99, mppext99));

    MassPoint mpext100 (925, 1);
    StopCrossSection csext100 (0.000981406, 0.000114539);
    MassPointParameters mppext100 (csext100, 10398);
    StopNeutralinoMap.insert(std::make_pair(mpext100, mppext100));

    MassPoint mpext101 (925, 25);
    StopCrossSection csext101 (0.000981406, 0.000114539);
    MassPointParameters mppext101 (csext101, 7579);
    StopNeutralinoMap.insert(std::make_pair(mpext101, mppext101));

    MassPoint mpext102 (925, 50);
    StopCrossSection csext102 (0.000981406, 0.000114539);
    MassPointParameters mppext102 (csext102, 10408);
    StopNeutralinoMap.insert(std::make_pair(mpext102, mppext102));

    MassPoint mpext103 (925, 75);
    StopCrossSection csext103 (0.000981406, 0.000114539);
    MassPointParameters mppext103 (csext103, 9682);
    StopNeutralinoMap.insert(std::make_pair(mpext103, mppext103));

    MassPoint mpext104 (925, 100);
    StopCrossSection csext104 (0.000981406, 0.000114539);
    MassPointParameters mppext104 (csext104, 11661);
    StopNeutralinoMap.insert(std::make_pair(mpext104, mppext104));

    MassPoint mpext105 (925, 125);
    StopCrossSection csext105 (0.000981406, 0.000114539);
    MassPointParameters mppext105 (csext105, 12416);
    StopNeutralinoMap.insert(std::make_pair(mpext105, mppext105));

    MassPoint mpext106 (925, 150);
    StopCrossSection csext106 (0.000981406, 0.000114539);
    MassPointParameters mppext106 (csext106, 9033);
    StopNeutralinoMap.insert(std::make_pair(mpext106, mppext106));

    MassPoint mpext107 (925, 175);
    StopCrossSection csext107 (0.000981406, 0.000114539);
    MassPointParameters mppext107 (csext107, 10037);
    StopNeutralinoMap.insert(std::make_pair(mpext107, mppext107));

    MassPoint mpext108 (925, 200);
    StopCrossSection csext108 (0.000981406, 0.000114539);
    MassPointParameters mppext108 (csext108, 9191);
    StopNeutralinoMap.insert(std::make_pair(mpext108, mppext108));

    MassPoint mpext109 (925, 225);
    StopCrossSection csext109 (0.000981406, 0.000114539);
    MassPointParameters mppext109 (csext109, 9674);
    StopNeutralinoMap.insert(std::make_pair(mpext109, mppext109));

    MassPoint mpext110 (925, 250);
    StopCrossSection csext110 (0.000981406, 0.000114539);
    MassPointParameters mppext110 (csext110, 6284);
    StopNeutralinoMap.insert(std::make_pair(mpext110, mppext110));

    MassPoint mpext111 (925, 275);
    StopCrossSection csext111 (0.000981406, 0.000114539);
    MassPointParameters mppext111 (csext111, 10927);
    StopNeutralinoMap.insert(std::make_pair(mpext111, mppext111));

    MassPoint mpext112 (925, 300);
    StopCrossSection csext112 (0.000981406, 0.000114539);
    MassPointParameters mppext112 (csext112, 10023);
    StopNeutralinoMap.insert(std::make_pair(mpext112, mppext112));

    MassPoint mpext113 (925, 325);
    StopCrossSection csext113 (0.000981406, 0.000114539);
    MassPointParameters mppext113 (csext113, 7857);
    StopNeutralinoMap.insert(std::make_pair(mpext113, mppext113));

    MassPoint mpext114 (925, 350);
    StopCrossSection csext114 (0.000981406, 0.000114539);
    MassPointParameters mppext114 (csext114, 8275);
    StopNeutralinoMap.insert(std::make_pair(mpext114, mppext114));

    MassPoint mpext115 (925, 375);
    StopCrossSection csext115 (0.000981406, 0.000114539);
    MassPointParameters mppext115 (csext115, 9336);
    StopNeutralinoMap.insert(std::make_pair(mpext115, mppext115));

    MassPoint mpext116 (925, 400);
    StopCrossSection csext116 (0.000981406, 0.000114539);
    MassPointParameters mppext116 (csext116, 9101);
    StopNeutralinoMap.insert(std::make_pair(mpext116, mppext116));

    MassPoint mpext117 (925, 425);
    StopCrossSection csext117 (0.000981406, 0.000114539);
    MassPointParameters mppext117 (csext117, 9155);
    StopNeutralinoMap.insert(std::make_pair(mpext117, mppext117));

    MassPoint mpext118 (925, 450);
    StopCrossSection csext118 (0.000981406, 0.000114539);
    MassPointParameters mppext118 (csext118, 11149);
    StopNeutralinoMap.insert(std::make_pair(mpext118, mppext118));

    MassPoint mpext119 (925, 475);
    StopCrossSection csext119 (0.000981406, 0.000114539);
    MassPointParameters mppext119 (csext119, 9630);
    StopNeutralinoMap.insert(std::make_pair(mpext119, mppext119));

    MassPoint mpext120 (925, 500);
    StopCrossSection csext120 (0.000981406, 0.000114539);
    MassPointParameters mppext120 (csext120, 8055);
    StopNeutralinoMap.insert(std::make_pair(mpext120, mppext120));

    MassPoint mpext121 (925, 525);
    StopCrossSection csext121 (0.000981406, 0.000114539);
    MassPointParameters mppext121 (csext121, 9426);
    StopNeutralinoMap.insert(std::make_pair(mpext121, mppext121));

    MassPoint mpext122 (925, 550);
    StopCrossSection csext122 (0.000981406, 0.000114539);
    MassPointParameters mppext122 (csext122, 11088);
    StopNeutralinoMap.insert(std::make_pair(mpext122, mppext122));

    MassPoint mpext123 (925, 575);
    StopCrossSection csext123 (0.000981406, 0.000114539);
    MassPointParameters mppext123 (csext123, 11436);
    StopNeutralinoMap.insert(std::make_pair(mpext123, mppext123));

    MassPoint mpext124 (925, 600);
    StopCrossSection csext124 (0.000981406, 0.000114539);
    MassPointParameters mppext124 (csext124, 11362);
    StopNeutralinoMap.insert(std::make_pair(mpext124, mppext124));

    MassPoint mpext125 (950, 1);
    StopCrossSection csext125 (0.000842779, 0.000102086);
    MassPointParameters mppext125 (csext125, 9436);
    StopNeutralinoMap.insert(std::make_pair(mpext125, mppext125));

    MassPoint mpext126 (950, 25);
    StopCrossSection csext126 (0.000842779, 0.000102086);
    MassPointParameters mppext126 (csext126, 6634);
    StopNeutralinoMap.insert(std::make_pair(mpext126, mppext126));

    MassPoint mpext127 (950, 50);
    StopCrossSection csext127 (0.000842779, 0.000102086);
    MassPointParameters mppext127 (csext127, 10563);
    StopNeutralinoMap.insert(std::make_pair(mpext127, mppext127));

    MassPoint mpext128 (950, 75);
    StopCrossSection csext128 (0.000842779, 0.000102086);
    MassPointParameters mppext128 (csext128, 11722);
    StopNeutralinoMap.insert(std::make_pair(mpext128, mppext128));

    MassPoint mpext129 (950, 100);
    StopCrossSection csext129 (0.000842779, 0.000102086);
    MassPointParameters mppext129 (csext129, 12044);
    StopNeutralinoMap.insert(std::make_pair(mpext129, mppext129));

    MassPoint mpext130 (950, 125);
    StopCrossSection csext130 (0.000842779, 0.000102086);
    MassPointParameters mppext130 (csext130, 8644);
    StopNeutralinoMap.insert(std::make_pair(mpext130, mppext130));

    MassPoint mpext131 (950, 150);
    StopCrossSection csext131 (0.000842779, 0.000102086);
    MassPointParameters mppext131 (csext131, 11906);
    StopNeutralinoMap.insert(std::make_pair(mpext131, mppext131));

    MassPoint mpext132 (950, 175);
    StopCrossSection csext132 (0.000842779, 0.000102086);
    MassPointParameters mppext132 (csext132, 9332);
    StopNeutralinoMap.insert(std::make_pair(mpext132, mppext132));

    MassPoint mpext133 (950, 200);
    StopCrossSection csext133 (0.000842779, 0.000102086);
    MassPointParameters mppext133 (csext133, 8493);
    StopNeutralinoMap.insert(std::make_pair(mpext133, mppext133));

    MassPoint mpext134 (950, 225);
    StopCrossSection csext134 (0.000842779, 0.000102086);
    MassPointParameters mppext134 (csext134, 8154);
    StopNeutralinoMap.insert(std::make_pair(mpext134, mppext134));

    MassPoint mpext135 (950, 250);
    StopCrossSection csext135 (0.000842779, 0.000102086);
    MassPointParameters mppext135 (csext135, 8354);
    StopNeutralinoMap.insert(std::make_pair(mpext135, mppext135));

    MassPoint mpext136 (950, 275);
    StopCrossSection csext136 (0.000842779, 0.000102086);
    MassPointParameters mppext136 (csext136, 9339);
    StopNeutralinoMap.insert(std::make_pair(mpext136, mppext136));

    MassPoint mpext137 (950, 300);
    StopCrossSection csext137 (0.000842779, 0.000102086);
    MassPointParameters mppext137 (csext137, 9293);
    StopNeutralinoMap.insert(std::make_pair(mpext137, mppext137));

    MassPoint mpext138 (950, 325);
    StopCrossSection csext138 (0.000842779, 0.000102086);
    MassPointParameters mppext138 (csext138, 9403);
    StopNeutralinoMap.insert(std::make_pair(mpext138, mppext138));

    MassPoint mpext139 (950, 350);
    StopCrossSection csext139 (0.000842779, 0.000102086);
    MassPointParameters mppext139 (csext139, 12125);
    StopNeutralinoMap.insert(std::make_pair(mpext139, mppext139));

    MassPoint mpext140 (950, 375);
    StopCrossSection csext140 (0.000842779, 0.000102086);
    MassPointParameters mppext140 (csext140, 11423);
    StopNeutralinoMap.insert(std::make_pair(mpext140, mppext140));

    MassPoint mpext141 (950, 400);
    StopCrossSection csext141 (0.000842779, 0.000102086);
    MassPointParameters mppext141 (csext141, 11618);
    StopNeutralinoMap.insert(std::make_pair(mpext141, mppext141));

    MassPoint mpext142 (950, 425);
    StopCrossSection csext142 (0.000842779, 0.000102086);
    MassPointParameters mppext142 (csext142, 8794);
    StopNeutralinoMap.insert(std::make_pair(mpext142, mppext142));

    MassPoint mpext143 (950, 450);
    StopCrossSection csext143 (0.000842779, 0.000102086);
    MassPointParameters mppext143 (csext143, 8153);
    StopNeutralinoMap.insert(std::make_pair(mpext143, mppext143));

    MassPoint mpext144 (950, 475);
    StopCrossSection csext144 (0.000842779, 0.000102086);
    MassPointParameters mppext144 (csext144, 10589);
    StopNeutralinoMap.insert(std::make_pair(mpext144, mppext144));

    MassPoint mpext145 (950, 500);
    StopCrossSection csext145 (0.000842779, 0.000102086);
    MassPointParameters mppext145 (csext145, 9449);
    StopNeutralinoMap.insert(std::make_pair(mpext145, mppext145));

    MassPoint mpext146 (950, 525);
    StopCrossSection csext146 (0.000842779, 0.000102086);
    MassPointParameters mppext146 (csext146, 7878);
    StopNeutralinoMap.insert(std::make_pair(mpext146, mppext146));

    MassPoint mpext147 (950, 550);
    StopCrossSection csext147 (0.000842779, 0.000102086);
    MassPointParameters mppext147 (csext147, 9004);
    StopNeutralinoMap.insert(std::make_pair(mpext147, mppext147));

    MassPoint mpext148 (950, 575);
    StopCrossSection csext148 (0.000842779, 0.000102086);
    MassPointParameters mppext148 (csext148, 6614);
    StopNeutralinoMap.insert(std::make_pair(mpext148, mppext148));

    MassPoint mpext149 (950, 600);
    StopCrossSection csext149 (0.000842779, 0.000102086);
    MassPointParameters mppext149 (csext149, 8445);
    StopNeutralinoMap.insert(std::make_pair(mpext149, mppext149));

    MassPoint mpext150 (975, 1);
    StopCrossSection csext150 (0.000713432, 7.79702e-05);
    MassPointParameters mppext150 (csext150, 10244);
    StopNeutralinoMap.insert(std::make_pair(mpext150, mppext150));

    MassPoint mpext151 (975, 25);
    StopCrossSection csext151 (0.000713432, 7.79702e-05);
    MassPointParameters mppext151 (csext151, 12411);
    StopNeutralinoMap.insert(std::make_pair(mpext151, mppext151));

    MassPoint mpext152 (975, 50);
    StopCrossSection csext152 (0.000713432, 7.79702e-05);
    MassPointParameters mppext152 (csext152, 10250);
    StopNeutralinoMap.insert(std::make_pair(mpext152, mppext152));

    MassPoint mpext153 (975, 75);
    StopCrossSection csext153 (0.000713432, 7.79702e-05);
    MassPointParameters mppext153 (csext153, 8546);
    StopNeutralinoMap.insert(std::make_pair(mpext153, mppext153));

    MassPoint mpext154 (975, 100);
    StopCrossSection csext154 (0.000713432, 7.79702e-05);
    MassPointParameters mppext154 (csext154, 7674);
    StopNeutralinoMap.insert(std::make_pair(mpext154, mppext154));

    MassPoint mpext155 (975, 125);
    StopCrossSection csext155 (0.000713432, 7.79702e-05);
    MassPointParameters mppext155 (csext155, 10585);
    StopNeutralinoMap.insert(std::make_pair(mpext155, mppext155));

    MassPoint mpext156 (975, 150);
    StopCrossSection csext156 (0.000713432, 7.79702e-05);
    MassPointParameters mppext156 (csext156, 9386);
    StopNeutralinoMap.insert(std::make_pair(mpext156, mppext156));

    MassPoint mpext157 (975, 175);
    StopCrossSection csext157 (0.000713432, 7.79702e-05);
    MassPointParameters mppext157 (csext157, 9535);
    StopNeutralinoMap.insert(std::make_pair(mpext157, mppext157));

    MassPoint mpext158 (975, 200);
    StopCrossSection csext158 (0.000713432, 7.79702e-05);
    MassPointParameters mppext158 (csext158, 8901);
    StopNeutralinoMap.insert(std::make_pair(mpext158, mppext158));

    MassPoint mpext159 (975, 225);
    StopCrossSection csext159 (0.000713432, 7.79702e-05);
    MassPointParameters mppext159 (csext159, 10233);
    StopNeutralinoMap.insert(std::make_pair(mpext159, mppext159));

    MassPoint mpext160 (975, 250);
    StopCrossSection csext160 (0.000713432, 7.79702e-05);
    MassPointParameters mppext160 (csext160, 9334);
    StopNeutralinoMap.insert(std::make_pair(mpext160, mppext160));

    MassPoint mpext161 (975, 275);
    StopCrossSection csext161 (0.000713432, 7.79702e-05);
    MassPointParameters mppext161 (csext161, 10234);
    StopNeutralinoMap.insert(std::make_pair(mpext161, mppext161));

    MassPoint mpext162 (975, 300);
    StopCrossSection csext162 (0.000713432, 7.79702e-05);
    MassPointParameters mppext162 (csext162, 10524);
    StopNeutralinoMap.insert(std::make_pair(mpext162, mppext162));

    MassPoint mpext163 (975, 325);
    StopCrossSection csext163 (0.000713432, 7.79702e-05);
    MassPointParameters mppext163 (csext163, 8362);
    StopNeutralinoMap.insert(std::make_pair(mpext163, mppext163));

    MassPoint mpext164 (975, 350);
    StopCrossSection csext164 (0.000713432, 7.79702e-05);
    MassPointParameters mppext164 (csext164, 10017);
    StopNeutralinoMap.insert(std::make_pair(mpext164, mppext164));

    MassPoint mpext165 (975, 375);
    StopCrossSection csext165 (0.000713432, 7.79702e-05);
    MassPointParameters mppext165 (csext165, 10744);
    StopNeutralinoMap.insert(std::make_pair(mpext165, mppext165));

    MassPoint mpext166 (975, 400);
    StopCrossSection csext166 (0.000713432, 7.79702e-05);
    MassPointParameters mppext166 (csext166, 9911);
    StopNeutralinoMap.insert(std::make_pair(mpext166, mppext166));

    MassPoint mpext167 (975, 425);
    StopCrossSection csext167 (0.000713432, 7.79702e-05);
    MassPointParameters mppext167 (csext167, 9444);
    StopNeutralinoMap.insert(std::make_pair(mpext167, mppext167));

    MassPoint mpext168 (975, 450);
    StopCrossSection csext168 (0.000713432, 7.79702e-05);
    MassPointParameters mppext168 (csext168, 10392);
    StopNeutralinoMap.insert(std::make_pair(mpext168, mppext168));

    MassPoint mpext169 (975, 475);
    StopCrossSection csext169 (0.000713432, 7.79702e-05);
    MassPointParameters mppext169 (csext169, 10054);
    StopNeutralinoMap.insert(std::make_pair(mpext169, mppext169));

    MassPoint mpext170 (975, 500);
    StopCrossSection csext170 (0.000713432, 7.79702e-05);
    MassPointParameters mppext170 (csext170, 8084);
    StopNeutralinoMap.insert(std::make_pair(mpext170, mppext170));

    MassPoint mpext171 (975, 525);
    StopCrossSection csext171 (0.000713432, 7.79702e-05);
    MassPointParameters mppext171 (csext171, 7420);
    StopNeutralinoMap.insert(std::make_pair(mpext171, mppext171));

    MassPoint mpext172 (975, 550);
    StopCrossSection csext172 (0.000713432, 7.79702e-05);
    MassPointParameters mppext172 (csext172, 9563);
    StopNeutralinoMap.insert(std::make_pair(mpext172, mppext172));

    MassPoint mpext173 (975, 575);
    StopCrossSection csext173 (0.000713432, 7.79702e-05);
    MassPointParameters mppext173 (csext173, 8492);
    StopNeutralinoMap.insert(std::make_pair(mpext173, mppext173));

    MassPoint mpext174 (975, 600);
    StopCrossSection csext174 (0.000713432, 7.79702e-05);
    MassPointParameters mppext174 (csext174, 7691);
    StopNeutralinoMap.insert(std::make_pair(mpext174, mppext174));

    MassPoint mpext175 (1000, 1);
    StopCrossSection csext175 (0.000621866, 7.71005e-05);
    MassPointParameters mppext175 (csext175, 10809);
    StopNeutralinoMap.insert(std::make_pair(mpext175, mppext175));

    MassPoint mpext176 (1000, 25);
    StopCrossSection csext176 (0.000621866, 7.71005e-05);
    MassPointParameters mppext176 (csext176, 10115);
    StopNeutralinoMap.insert(std::make_pair(mpext176, mppext176));

    MassPoint mpext177 (1000, 50);
    StopCrossSection csext177 (0.000621866, 7.71005e-05);
    MassPointParameters mppext177 (csext177, 10214);
    StopNeutralinoMap.insert(std::make_pair(mpext177, mppext177));

    MassPoint mpext178 (1000, 75);
    StopCrossSection csext178 (0.000621866, 7.71005e-05);
    MassPointParameters mppext178 (csext178, 10123);
    StopNeutralinoMap.insert(std::make_pair(mpext178, mppext178));

    MassPoint mpext179 (1000, 100);
    StopCrossSection csext179 (0.000621866, 7.71005e-05);
    MassPointParameters mppext179 (csext179, 10889);
    StopNeutralinoMap.insert(std::make_pair(mpext179, mppext179));

    MassPoint mpext180 (1000, 125);
    StopCrossSection csext180 (0.000621866, 7.71005e-05);
    MassPointParameters mppext180 (csext180, 12454);
    StopNeutralinoMap.insert(std::make_pair(mpext180, mppext180));

    MassPoint mpext181 (1000, 150);
    StopCrossSection csext181 (0.000621866, 7.71005e-05);
    MassPointParameters mppext181 (csext181, 8607);
    StopNeutralinoMap.insert(std::make_pair(mpext181, mppext181));

    MassPoint mpext182 (1000, 175);
    StopCrossSection csext182 (0.000621866, 7.71005e-05);
    MassPointParameters mppext182 (csext182, 7451);
    StopNeutralinoMap.insert(std::make_pair(mpext182, mppext182));

    MassPoint mpext183 (1000, 200);
    StopCrossSection csext183 (0.000621866, 7.71005e-05);
    MassPointParameters mppext183 (csext183, 10304);
    StopNeutralinoMap.insert(std::make_pair(mpext183, mppext183));

    MassPoint mpext184 (1000, 225);
    StopCrossSection csext184 (0.000621866, 7.71005e-05);
    MassPointParameters mppext184 (csext184, 9579);
    StopNeutralinoMap.insert(std::make_pair(mpext184, mppext184));

    MassPoint mpext185 (1000, 250);
    StopCrossSection csext185 (0.000621866, 7.71005e-05);
    MassPointParameters mppext185 (csext185, 10702);
    StopNeutralinoMap.insert(std::make_pair(mpext185, mppext185));

    MassPoint mpext186 (1000, 275);
    StopCrossSection csext186 (0.000621866, 7.71005e-05);
    MassPointParameters mppext186 (csext186, 9615);
    StopNeutralinoMap.insert(std::make_pair(mpext186, mppext186));

    MassPoint mpext187 (1000, 300);
    StopCrossSection csext187 (0.000621866, 7.71005e-05);
    MassPointParameters mppext187 (csext187, 10469);
    StopNeutralinoMap.insert(std::make_pair(mpext187, mppext187));

    MassPoint mpext188 (1000, 325);
    StopCrossSection csext188 (0.000621866, 7.71005e-05);
    MassPointParameters mppext188 (csext188, 11374);
    StopNeutralinoMap.insert(std::make_pair(mpext188, mppext188));

    MassPoint mpext189 (1000, 350);
    StopCrossSection csext189 (0.000621866, 7.71005e-05);
    MassPointParameters mppext189 (csext189, 9928);
    StopNeutralinoMap.insert(std::make_pair(mpext189, mppext189));

    MassPoint mpext190 (1000, 375);
    StopCrossSection csext190 (0.000621866, 7.71005e-05);
    MassPointParameters mppext190 (csext190, 10025);
    StopNeutralinoMap.insert(std::make_pair(mpext190, mppext190));

    MassPoint mpext191 (1000, 400);
    StopCrossSection csext191 (0.000621866, 7.71005e-05);
    MassPointParameters mppext191 (csext191, 11254);
    StopNeutralinoMap.insert(std::make_pair(mpext191, mppext191));

    MassPoint mpext192 (1000, 425);
    StopCrossSection csext192 (0.000621866, 7.71005e-05);
    MassPointParameters mppext192 (csext192, 8227);
    StopNeutralinoMap.insert(std::make_pair(mpext192, mppext192));

    MassPoint mpext193 (1000, 450);
    StopCrossSection csext193 (0.000621866, 7.71005e-05);
    MassPointParameters mppext193 (csext193, 8365);
    StopNeutralinoMap.insert(std::make_pair(mpext193, mppext193));

    MassPoint mpext194 (1000, 475);
    StopCrossSection csext194 (0.000621866, 7.71005e-05);
    MassPointParameters mppext194 (csext194, 9310);
    StopNeutralinoMap.insert(std::make_pair(mpext194, mppext194));

    MassPoint mpext195 (1000, 500);
    StopCrossSection csext195 (0.000621866, 7.71005e-05);
    MassPointParameters mppext195 (csext195, 7832);
    StopNeutralinoMap.insert(std::make_pair(mpext195, mppext195));

    MassPoint mpext196 (1000, 525);
    StopCrossSection csext196 (0.000621866, 7.71005e-05);
    MassPointParameters mppext196 (csext196, 11109);
    StopNeutralinoMap.insert(std::make_pair(mpext196, mppext196));

    MassPoint mpext197 (1000, 550);
    StopCrossSection csext197 (0.000621866, 7.71005e-05);
    MassPointParameters mppext197 (csext197, 9124);
    StopNeutralinoMap.insert(std::make_pair(mpext197, mppext197));

    MassPoint mpext198 (1000, 575);
    StopCrossSection csext198 (0.000621866, 7.71005e-05);
    MassPointParameters mppext198 (csext198, 11023);
    StopNeutralinoMap.insert(std::make_pair(mpext198, mppext198));

    MassPoint mpext199 (1000, 600);
    StopCrossSection csext199 (0.000621866, 7.71005e-05);
    MassPointParameters mppext199 (csext199, 12278);
    StopNeutralinoMap.insert(std::make_pair(mpext199, mppext199));

    MassPoint mpext200 (1025, 1);
    StopCrossSection csext200 (0.000535563, 6.67594e-05);
    MassPointParameters mppext200 (csext200, 8067);
    StopNeutralinoMap.insert(std::make_pair(mpext200, mppext200));

    MassPoint mpext201 (1025, 25);
    StopCrossSection csext201 (0.000535563, 6.67594e-05);
    MassPointParameters mppext201 (csext201, 10949);
    StopNeutralinoMap.insert(std::make_pair(mpext201, mppext201));

    MassPoint mpext202 (1025, 50);
    StopCrossSection csext202 (0.000535563, 6.67594e-05);
    MassPointParameters mppext202 (csext202, 8913);
    StopNeutralinoMap.insert(std::make_pair(mpext202, mppext202));

    MassPoint mpext203 (1025, 75);
    StopCrossSection csext203 (0.000535563, 6.67594e-05);
    MassPointParameters mppext203 (csext203, 8209);
    StopNeutralinoMap.insert(std::make_pair(mpext203, mppext203));

    MassPoint mpext204 (1025, 100);
    StopCrossSection csext204 (0.000535563, 6.67594e-05);
    MassPointParameters mppext204 (csext204, 9162);
    StopNeutralinoMap.insert(std::make_pair(mpext204, mppext204));

    MassPoint mpext205 (1025, 125);
    StopCrossSection csext205 (0.000535563, 6.67594e-05);
    MassPointParameters mppext205 (csext205, 10853);
    StopNeutralinoMap.insert(std::make_pair(mpext205, mppext205));

    MassPoint mpext206 (1025, 150);
    StopCrossSection csext206 (0.000535563, 6.67594e-05);
    MassPointParameters mppext206 (csext206, 8173);
    StopNeutralinoMap.insert(std::make_pair(mpext206, mppext206));

    MassPoint mpext207 (1025, 175);
    StopCrossSection csext207 (0.000535563, 6.67594e-05);
    MassPointParameters mppext207 (csext207, 10316);
    StopNeutralinoMap.insert(std::make_pair(mpext207, mppext207));

    MassPoint mpext208 (1025, 200);
    StopCrossSection csext208 (0.000535563, 6.67594e-05);
    MassPointParameters mppext208 (csext208, 8712);
    StopNeutralinoMap.insert(std::make_pair(mpext208, mppext208));

    MassPoint mpext209 (1025, 225);
    StopCrossSection csext209 (0.000535563, 6.67594e-05);
    MassPointParameters mppext209 (csext209, 8602);
    StopNeutralinoMap.insert(std::make_pair(mpext209, mppext209));

    MassPoint mpext210 (1025, 250);
    StopCrossSection csext210 (0.000535563, 6.67594e-05);
    MassPointParameters mppext210 (csext210, 11096);
    StopNeutralinoMap.insert(std::make_pair(mpext210, mppext210));

    MassPoint mpext211 (1025, 275);
    StopCrossSection csext211 (0.000535563, 6.67594e-05);
    MassPointParameters mppext211 (csext211, 9534);
    StopNeutralinoMap.insert(std::make_pair(mpext211, mppext211));

    MassPoint mpext212 (1025, 300);
    StopCrossSection csext212 (0.000535563, 6.67594e-05);
    MassPointParameters mppext212 (csext212, 11208);
    StopNeutralinoMap.insert(std::make_pair(mpext212, mppext212));

    MassPoint mpext213 (1025, 325);
    StopCrossSection csext213 (0.000535563, 6.67594e-05);
    MassPointParameters mppext213 (csext213, 9077);
    StopNeutralinoMap.insert(std::make_pair(mpext213, mppext213));

    MassPoint mpext214 (1025, 350);
    StopCrossSection csext214 (0.000535563, 6.67594e-05);
    MassPointParameters mppext214 (csext214, 9430);
    StopNeutralinoMap.insert(std::make_pair(mpext214, mppext214));

    MassPoint mpext215 (1025, 375);
    StopCrossSection csext215 (0.000535563, 6.67594e-05);
    MassPointParameters mppext215 (csext215, 9427);
    StopNeutralinoMap.insert(std::make_pair(mpext215, mppext215));

    MassPoint mpext216 (1025, 400);
    StopCrossSection csext216 (0.000535563, 6.67594e-05);
    MassPointParameters mppext216 (csext216, 12597);
    StopNeutralinoMap.insert(std::make_pair(mpext216, mppext216));

    MassPoint mpext217 (1025, 425);
    StopCrossSection csext217 (0.000535563, 6.67594e-05);
    MassPointParameters mppext217 (csext217, 8672);
    StopNeutralinoMap.insert(std::make_pair(mpext217, mppext217));

    MassPoint mpext218 (1025, 450);
    StopCrossSection csext218 (0.000535563, 6.67594e-05);
    MassPointParameters mppext218 (csext218, 7805);
    StopNeutralinoMap.insert(std::make_pair(mpext218, mppext218));

    MassPoint mpext219 (1025, 475);
    StopCrossSection csext219 (0.000535563, 6.67594e-05);
    MassPointParameters mppext219 (csext219, 7463);
    StopNeutralinoMap.insert(std::make_pair(mpext219, mppext219));

    MassPoint mpext220 (1025, 500);
    StopCrossSection csext220 (0.000535563, 6.67594e-05);
    MassPointParameters mppext220 (csext220, 10623);
    StopNeutralinoMap.insert(std::make_pair(mpext220, mppext220));

    MassPoint mpext221 (1025, 525);
    StopCrossSection csext221 (0.000535563, 6.67594e-05);
    MassPointParameters mppext221 (csext221, 9074);
    StopNeutralinoMap.insert(std::make_pair(mpext221, mppext221));

    MassPoint mpext222 (1025, 550);
    StopCrossSection csext222 (0.000535563, 6.67594e-05);
    MassPointParameters mppext222 (csext222, 7078);
    StopNeutralinoMap.insert(std::make_pair(mpext222, mppext222));

    MassPoint mpext223 (1025, 575);
    StopCrossSection csext223 (0.000535563, 6.67594e-05);
    MassPointParameters mppext223 (csext223, 12158);
    StopNeutralinoMap.insert(std::make_pair(mpext223, mppext223));

    MassPoint mpext224 (1025, 600);
    StopCrossSection csext224 (0.000535563, 6.67594e-05);
    MassPointParameters mppext224 (csext224, 9815);
    StopNeutralinoMap.insert(std::make_pair(mpext224, mppext224));

    MassPoint mpext225 (1050, 1);
    StopCrossSection csext225 (0.000458716, 5.69349e-05);
    MassPointParameters mppext225 (csext225, 10414);
    StopNeutralinoMap.insert(std::make_pair(mpext225, mppext225));

    MassPoint mpext226 (1050, 25);
    StopCrossSection csext226 (0.000458716, 5.69349e-05);
    MassPointParameters mppext226 (csext226, 8768);
    StopNeutralinoMap.insert(std::make_pair(mpext226, mppext226));

    MassPoint mpext227 (1050, 50);
    StopCrossSection csext227 (0.000458716, 5.69349e-05);
    MassPointParameters mppext227 (csext227, 10676);
    StopNeutralinoMap.insert(std::make_pair(mpext227, mppext227));

    MassPoint mpext228 (1050, 75);
    StopCrossSection csext228 (0.000458716, 5.69349e-05);
    MassPointParameters mppext228 (csext228, 7387);
    StopNeutralinoMap.insert(std::make_pair(mpext228, mppext228));

    MassPoint mpext229 (1050, 100);
    StopCrossSection csext229 (0.000458716, 5.69349e-05);
    MassPointParameters mppext229 (csext229, 10135);
    StopNeutralinoMap.insert(std::make_pair(mpext229, mppext229));

    MassPoint mpext230 (1050, 125);
    StopCrossSection csext230 (0.000458716, 5.69349e-05);
    MassPointParameters mppext230 (csext230, 8982);
    StopNeutralinoMap.insert(std::make_pair(mpext230, mppext230));

    MassPoint mpext231 (1050, 150);
    StopCrossSection csext231 (0.000458716, 5.69349e-05);
    MassPointParameters mppext231 (csext231, 8050);
    StopNeutralinoMap.insert(std::make_pair(mpext231, mppext231));

    MassPoint mpext232 (1050, 175);
    StopCrossSection csext232 (0.000458716, 5.69349e-05);
    MassPointParameters mppext232 (csext232, 11223);
    StopNeutralinoMap.insert(std::make_pair(mpext232, mppext232));

    MassPoint mpext233 (1050, 200);
    StopCrossSection csext233 (0.000458716, 5.69349e-05);
    MassPointParameters mppext233 (csext233, 9934);
    StopNeutralinoMap.insert(std::make_pair(mpext233, mppext233));

    MassPoint mpext234 (1050, 225);
    StopCrossSection csext234 (0.000458716, 5.69349e-05);
    MassPointParameters mppext234 (csext234, 9892);
    StopNeutralinoMap.insert(std::make_pair(mpext234, mppext234));

    MassPoint mpext235 (1050, 250);
    StopCrossSection csext235 (0.000458716, 5.69349e-05);
    MassPointParameters mppext235 (csext235, 8524);
    StopNeutralinoMap.insert(std::make_pair(mpext235, mppext235));

    MassPoint mpext236 (1050, 275);
    StopCrossSection csext236 (0.000458716, 5.69349e-05);
    MassPointParameters mppext236 (csext236, 9754);
    StopNeutralinoMap.insert(std::make_pair(mpext236, mppext236));

    MassPoint mpext237 (1050, 300);
    StopCrossSection csext237 (0.000458716, 5.69349e-05);
    MassPointParameters mppext237 (csext237, 10192);
    StopNeutralinoMap.insert(std::make_pair(mpext237, mppext237));

    MassPoint mpext238 (1050, 325);
    StopCrossSection csext238 (0.000458716, 5.69349e-05);
    MassPointParameters mppext238 (csext238, 11288);
    StopNeutralinoMap.insert(std::make_pair(mpext238, mppext238));

    MassPoint mpext239 (1050, 350);
    StopCrossSection csext239 (0.000458716, 5.69349e-05);
    MassPointParameters mppext239 (csext239, 9901);
    StopNeutralinoMap.insert(std::make_pair(mpext239, mppext239));

    MassPoint mpext240 (1050, 375);
    StopCrossSection csext240 (0.000458716, 5.69349e-05);
    MassPointParameters mppext240 (csext240, 9937);
    StopNeutralinoMap.insert(std::make_pair(mpext240, mppext240));

    MassPoint mpext241 (1050, 400);
    StopCrossSection csext241 (0.000458716, 5.69349e-05);
    MassPointParameters mppext241 (csext241, 9876);
    StopNeutralinoMap.insert(std::make_pair(mpext241, mppext241));

    MassPoint mpext242 (1050, 425);
    StopCrossSection csext242 (0.000458716, 5.69349e-05);
    MassPointParameters mppext242 (csext242, 10499);
    StopNeutralinoMap.insert(std::make_pair(mpext242, mppext242));

    MassPoint mpext243 (1050, 450);
    StopCrossSection csext243 (0.000458716, 5.69349e-05);
    MassPointParameters mppext243 (csext243, 10409);
    StopNeutralinoMap.insert(std::make_pair(mpext243, mppext243));

    MassPoint mpext244 (1050, 475);
    StopCrossSection csext244 (0.000458716, 5.69349e-05);
    MassPointParameters mppext244 (csext244, 9524);
    StopNeutralinoMap.insert(std::make_pair(mpext244, mppext244));

    MassPoint mpext245 (1050, 500);
    StopCrossSection csext245 (0.000458716, 5.69349e-05);
    MassPointParameters mppext245 (csext245, 8551);
    StopNeutralinoMap.insert(std::make_pair(mpext245, mppext245));

    MassPoint mpext246 (1050, 525);
    StopCrossSection csext246 (0.000458716, 5.69349e-05);
    MassPointParameters mppext246 (csext246, 10065);
    StopNeutralinoMap.insert(std::make_pair(mpext246, mppext246));

    MassPoint mpext247 (1050, 550);
    StopCrossSection csext247 (0.000458716, 5.69349e-05);
    MassPointParameters mppext247 (csext247, 9879);
    StopNeutralinoMap.insert(std::make_pair(mpext247, mppext247));

    MassPoint mpext248 (1050, 575);
    StopCrossSection csext248 (0.000458716, 5.69349e-05);
    MassPointParameters mppext248 (csext248, 9497);
    StopNeutralinoMap.insert(std::make_pair(mpext248, mppext248));

    MassPoint mpext249 (1050, 600);
    StopCrossSection csext249 (0.000458716, 5.69349e-05);
    MassPointParameters mppext249 (csext249, 11842);
    StopNeutralinoMap.insert(std::make_pair(mpext249, mppext249));

    MassPoint mpext250 (1075, 1);
    StopCrossSection csext250 (0.000398794, 5.06191e-05);
    MassPointParameters mppext250 (csext250, 7072);
    StopNeutralinoMap.insert(std::make_pair(mpext250, mppext250));

    MassPoint mpext251 (1075, 25);
    StopCrossSection csext251 (0.000398794, 5.06191e-05);
    MassPointParameters mppext251 (csext251, 8566);
    StopNeutralinoMap.insert(std::make_pair(mpext251, mppext251));

    MassPoint mpext252 (1075, 50);
    StopCrossSection csext252 (0.000398794, 5.06191e-05);
    MassPointParameters mppext252 (csext252, 9073);
    StopNeutralinoMap.insert(std::make_pair(mpext252, mppext252));

    MassPoint mpext253 (1075, 75);
    StopCrossSection csext253 (0.000398794, 5.06191e-05);
    MassPointParameters mppext253 (csext253, 10104);
    StopNeutralinoMap.insert(std::make_pair(mpext253, mppext253));

    MassPoint mpext254 (1075, 100);
    StopCrossSection csext254 (0.000398794, 5.06191e-05);
    MassPointParameters mppext254 (csext254, 10261);
    StopNeutralinoMap.insert(std::make_pair(mpext254, mppext254));

    MassPoint mpext255 (1075, 125);
    StopCrossSection csext255 (0.000398794, 5.06191e-05);
    MassPointParameters mppext255 (csext255, 8456);
    StopNeutralinoMap.insert(std::make_pair(mpext255, mppext255));

    MassPoint mpext256 (1075, 150);
    StopCrossSection csext256 (0.000398794, 5.06191e-05);
    MassPointParameters mppext256 (csext256, 11801);
    StopNeutralinoMap.insert(std::make_pair(mpext256, mppext256));

    MassPoint mpext257 (1075, 175);
    StopCrossSection csext257 (0.000398794, 5.06191e-05);
    MassPointParameters mppext257 (csext257, 7624);
    StopNeutralinoMap.insert(std::make_pair(mpext257, mppext257));

    MassPoint mpext258 (1075, 200);
    StopCrossSection csext258 (0.000398794, 5.06191e-05);
    MassPointParameters mppext258 (csext258, 8976);
    StopNeutralinoMap.insert(std::make_pair(mpext258, mppext258));

    MassPoint mpext259 (1075, 225);
    StopCrossSection csext259 (0.000398794, 5.06191e-05);
    MassPointParameters mppext259 (csext259, 9402);
    StopNeutralinoMap.insert(std::make_pair(mpext259, mppext259));

    MassPoint mpext260 (1075, 250);
    StopCrossSection csext260 (0.000398794, 5.06191e-05);
    MassPointParameters mppext260 (csext260, 9299);
    StopNeutralinoMap.insert(std::make_pair(mpext260, mppext260));

    MassPoint mpext261 (1075, 275);
    StopCrossSection csext261 (0.000398794, 5.06191e-05);
    MassPointParameters mppext261 (csext261, 7919);
    StopNeutralinoMap.insert(std::make_pair(mpext261, mppext261));

    MassPoint mpext262 (1075, 300);
    StopCrossSection csext262 (0.000398794, 5.06191e-05);
    MassPointParameters mppext262 (csext262, 8605);
    StopNeutralinoMap.insert(std::make_pair(mpext262, mppext262));

    MassPoint mpext263 (1075, 325);
    StopCrossSection csext263 (0.000398794, 5.06191e-05);
    MassPointParameters mppext263 (csext263, 7110);
    StopNeutralinoMap.insert(std::make_pair(mpext263, mppext263));

    MassPoint mpext264 (1075, 350);
    StopCrossSection csext264 (0.000398794, 5.06191e-05);
    MassPointParameters mppext264 (csext264, 8429);
    StopNeutralinoMap.insert(std::make_pair(mpext264, mppext264));

    MassPoint mpext265 (1075, 375);
    StopCrossSection csext265 (0.000398794, 5.06191e-05);
    MassPointParameters mppext265 (csext265, 7753);
    StopNeutralinoMap.insert(std::make_pair(mpext265, mppext265));

    MassPoint mpext266 (1075, 400);
    StopCrossSection csext266 (0.000398794, 5.06191e-05);
    MassPointParameters mppext266 (csext266, 11266);
    StopNeutralinoMap.insert(std::make_pair(mpext266, mppext266));

    MassPoint mpext267 (1075, 425);
    StopCrossSection csext267 (0.000398794, 5.06191e-05);
    MassPointParameters mppext267 (csext267, 10700);
    StopNeutralinoMap.insert(std::make_pair(mpext267, mppext267));

    MassPoint mpext268 (1075, 450);
    StopCrossSection csext268 (0.000398794, 5.06191e-05);
    MassPointParameters mppext268 (csext268, 8702);
    StopNeutralinoMap.insert(std::make_pair(mpext268, mppext268));

    MassPoint mpext269 (1075, 475);
    StopCrossSection csext269 (0.000398794, 5.06191e-05);
    MassPointParameters mppext269 (csext269, 11726);
    StopNeutralinoMap.insert(std::make_pair(mpext269, mppext269));

    MassPoint mpext270 (1075, 500);
    StopCrossSection csext270 (0.000398794, 5.06191e-05);
    MassPointParameters mppext270 (csext270, 8726);
    StopNeutralinoMap.insert(std::make_pair(mpext270, mppext270));

    MassPoint mpext271 (1075, 525);
    StopCrossSection csext271 (0.000398794, 5.06191e-05);
    MassPointParameters mppext271 (csext271, 10126);
    StopNeutralinoMap.insert(std::make_pair(mpext271, mppext271));

    MassPoint mpext272 (1075, 550);
    StopCrossSection csext272 (0.000398794, 5.06191e-05);
    MassPointParameters mppext272 (csext272, 9370);
    StopNeutralinoMap.insert(std::make_pair(mpext272, mppext272));

    MassPoint mpext273 (1075, 575);
    StopCrossSection csext273 (0.000398794, 5.06191e-05);
    MassPointParameters mppext273 (csext273, 8625);
    StopNeutralinoMap.insert(std::make_pair(mpext273, mppext273));

    MassPoint mpext274 (1075, 600);
    StopCrossSection csext274 (0.000398794, 5.06191e-05);
    MassPointParameters mppext274 (csext274, 9305);
    StopNeutralinoMap.insert(std::make_pair(mpext274, mppext274));

    MassPoint mpext275 (1100, 1);
    StopCrossSection csext275 (0.000342626, 4.27672e-05);
    MassPointParameters mppext275 (csext275, 7984);
    StopNeutralinoMap.insert(std::make_pair(mpext275, mppext275));

    MassPoint mpext276 (1100, 25);
    StopCrossSection csext276 (0.000342626, 4.27672e-05);
    MassPointParameters mppext276 (csext276, 7954);
    StopNeutralinoMap.insert(std::make_pair(mpext276, mppext276));

    MassPoint mpext277 (1100, 50);
    StopCrossSection csext277 (0.000342626, 4.27672e-05);
    MassPointParameters mppext277 (csext277, 7909);
    StopNeutralinoMap.insert(std::make_pair(mpext277, mppext277));

    MassPoint mpext278 (1100, 75);
    StopCrossSection csext278 (0.000342626, 4.27672e-05);
    MassPointParameters mppext278 (csext278, 10380);
    StopNeutralinoMap.insert(std::make_pair(mpext278, mppext278));

    MassPoint mpext279 (1100, 100);
    StopCrossSection csext279 (0.000342626, 4.27672e-05);
    MassPointParameters mppext279 (csext279, 11421);
    StopNeutralinoMap.insert(std::make_pair(mpext279, mppext279));

    MassPoint mpext280 (1100, 125);
    StopCrossSection csext280 (0.000342626, 4.27672e-05);
    MassPointParameters mppext280 (csext280, 8404);
    StopNeutralinoMap.insert(std::make_pair(mpext280, mppext280));

    MassPoint mpext281 (1100, 150);
    StopCrossSection csext281 (0.000342626, 4.27672e-05);
    MassPointParameters mppext281 (csext281, 9276);
    StopNeutralinoMap.insert(std::make_pair(mpext281, mppext281));

    MassPoint mpext282 (1100, 175);
    StopCrossSection csext282 (0.000342626, 4.27672e-05);
    MassPointParameters mppext282 (csext282, 9855);
    StopNeutralinoMap.insert(std::make_pair(mpext282, mppext282));

    MassPoint mpext283 (1100, 200);
    StopCrossSection csext283 (0.000342626, 4.27672e-05);
    MassPointParameters mppext283 (csext283, 8638);
    StopNeutralinoMap.insert(std::make_pair(mpext283, mppext283));

    MassPoint mpext284 (1100, 225);
    StopCrossSection csext284 (0.000342626, 4.27672e-05);
    MassPointParameters mppext284 (csext284, 9938);
    StopNeutralinoMap.insert(std::make_pair(mpext284, mppext284));

    MassPoint mpext285 (1100, 250);
    StopCrossSection csext285 (0.000342626, 4.27672e-05);
    MassPointParameters mppext285 (csext285, 8982);
    StopNeutralinoMap.insert(std::make_pair(mpext285, mppext285));

    MassPoint mpext286 (1100, 275);
    StopCrossSection csext286 (0.000342626, 4.27672e-05);
    MassPointParameters mppext286 (csext286, 10128);
    StopNeutralinoMap.insert(std::make_pair(mpext286, mppext286));

    MassPoint mpext287 (1100, 300);
    StopCrossSection csext287 (0.000342626, 4.27672e-05);
    MassPointParameters mppext287 (csext287, 8983);
    StopNeutralinoMap.insert(std::make_pair(mpext287, mppext287));

    MassPoint mpext288 (1100, 325);
    StopCrossSection csext288 (0.000342626, 4.27672e-05);
    MassPointParameters mppext288 (csext288, 9525);
    StopNeutralinoMap.insert(std::make_pair(mpext288, mppext288));

    MassPoint mpext289 (1100, 350);
    StopCrossSection csext289 (0.000342626, 4.27672e-05);
    MassPointParameters mppext289 (csext289, 10102);
    StopNeutralinoMap.insert(std::make_pair(mpext289, mppext289));

    MassPoint mpext290 (1100, 375);
    StopCrossSection csext290 (0.000342626, 4.27672e-05);
    MassPointParameters mppext290 (csext290, 10754);
    StopNeutralinoMap.insert(std::make_pair(mpext290, mppext290));

    MassPoint mpext291 (1100, 400);
    StopCrossSection csext291 (0.000342626, 4.27672e-05);
    MassPointParameters mppext291 (csext291, 11141);
    StopNeutralinoMap.insert(std::make_pair(mpext291, mppext291));

    MassPoint mpext292 (1100, 425);
    StopCrossSection csext292 (0.000342626, 4.27672e-05);
    MassPointParameters mppext292 (csext292, 8779);
    StopNeutralinoMap.insert(std::make_pair(mpext292, mppext292));

    MassPoint mpext293 (1100, 450);
    StopCrossSection csext293 (0.000342626, 4.27672e-05);
    MassPointParameters mppext293 (csext293, 8701);
    StopNeutralinoMap.insert(std::make_pair(mpext293, mppext293));

    MassPoint mpext294 (1100, 475);
    StopCrossSection csext294 (0.000342626, 4.27672e-05);
    MassPointParameters mppext294 (csext294, 6195);
    StopNeutralinoMap.insert(std::make_pair(mpext294, mppext294));

    MassPoint mpext295 (1100, 500);
    StopCrossSection csext295 (0.000342626, 4.27672e-05);
    MassPointParameters mppext295 (csext295, 7495);
    StopNeutralinoMap.insert(std::make_pair(mpext295, mppext295));

    MassPoint mpext296 (1100, 525);
    StopCrossSection csext296 (0.000342626, 4.27672e-05);
    MassPointParameters mppext296 (csext296, 9848);
    StopNeutralinoMap.insert(std::make_pair(mpext296, mppext296));

    MassPoint mpext297 (1100, 550);
    StopCrossSection csext297 (0.000342626, 4.27672e-05);
    MassPointParameters mppext297 (csext297, 11707);
    StopNeutralinoMap.insert(std::make_pair(mpext297, mppext297));

    MassPoint mpext298 (1100, 575);
    StopCrossSection csext298 (0.000342626, 4.27672e-05);
    MassPointParameters mppext298 (csext298, 10740);
    StopNeutralinoMap.insert(std::make_pair(mpext298, mppext298));

    MassPoint mpext299 (1100, 600);
    StopCrossSection csext299 (0.000342626, 4.27672e-05);
    MassPointParameters mppext299 (csext299, 9088);
    StopNeutralinoMap.insert(std::make_pair(mpext299, mppext299));

    MassPoint mpext300 (1125, 1);
    StopCrossSection csext300 (0.000301119, 4.14674e-05);
    MassPointParameters mppext300 (csext300, 10980);
    StopNeutralinoMap.insert(std::make_pair(mpext300, mppext300));

    MassPoint mpext301 (1125, 25);
    StopCrossSection csext301 (0.000301119, 4.14674e-05);
    MassPointParameters mppext301 (csext301, 8888);
    StopNeutralinoMap.insert(std::make_pair(mpext301, mppext301));

    MassPoint mpext302 (1125, 50);
    StopCrossSection csext302 (0.000301119, 4.14674e-05);
    MassPointParameters mppext302 (csext302, 10348);
    StopNeutralinoMap.insert(std::make_pair(mpext302, mppext302));

    MassPoint mpext303 (1125, 75);
    StopCrossSection csext303 (0.000301119, 4.14674e-05);
    MassPointParameters mppext303 (csext303, 10740);
    StopNeutralinoMap.insert(std::make_pair(mpext303, mppext303));

    MassPoint mpext304 (1125, 100);
    StopCrossSection csext304 (0.000301119, 4.14674e-05);
    MassPointParameters mppext304 (csext304, 9565);
    StopNeutralinoMap.insert(std::make_pair(mpext304, mppext304));

    MassPoint mpext305 (1125, 125);
    StopCrossSection csext305 (0.000301119, 4.14674e-05);
    MassPointParameters mppext305 (csext305, 7600);
    StopNeutralinoMap.insert(std::make_pair(mpext305, mppext305));

    MassPoint mpext306 (1125, 150);
    StopCrossSection csext306 (0.000301119, 4.14674e-05);
    MassPointParameters mppext306 (csext306, 9610);
    StopNeutralinoMap.insert(std::make_pair(mpext306, mppext306));

    MassPoint mpext307 (1125, 175);
    StopCrossSection csext307 (0.000301119, 4.14674e-05);
    MassPointParameters mppext307 (csext307, 9636);
    StopNeutralinoMap.insert(std::make_pair(mpext307, mppext307));

    MassPoint mpext308 (1125, 200);
    StopCrossSection csext308 (0.000301119, 4.14674e-05);
    MassPointParameters mppext308 (csext308, 9063);
    StopNeutralinoMap.insert(std::make_pair(mpext308, mppext308));

    MassPoint mpext309 (1125, 225);
    StopCrossSection csext309 (0.000301119, 4.14674e-05);
    MassPointParameters mppext309 (csext309, 9929);
    StopNeutralinoMap.insert(std::make_pair(mpext309, mppext309));

    MassPoint mpext310 (1125, 250);
    StopCrossSection csext310 (0.000301119, 4.14674e-05);
    MassPointParameters mppext310 (csext310, 11759);
    StopNeutralinoMap.insert(std::make_pair(mpext310, mppext310));

    MassPoint mpext311 (1125, 275);
    StopCrossSection csext311 (0.000301119, 4.14674e-05);
    MassPointParameters mppext311 (csext311, 10582);
    StopNeutralinoMap.insert(std::make_pair(mpext311, mppext311));

    MassPoint mpext312 (1125, 300);
    StopCrossSection csext312 (0.000301119, 4.14674e-05);
    MassPointParameters mppext312 (csext312, 9749);
    StopNeutralinoMap.insert(std::make_pair(mpext312, mppext312));

    MassPoint mpext313 (1125, 325);
    StopCrossSection csext313 (0.000301119, 4.14674e-05);
    MassPointParameters mppext313 (csext313, 9386);
    StopNeutralinoMap.insert(std::make_pair(mpext313, mppext313));

    MassPoint mpext314 (1125, 350);
    StopCrossSection csext314 (0.000301119, 4.14674e-05);
    MassPointParameters mppext314 (csext314, 9001);
    StopNeutralinoMap.insert(std::make_pair(mpext314, mppext314));

    MassPoint mpext315 (1125, 375);
    StopCrossSection csext315 (0.000301119, 4.14674e-05);
    MassPointParameters mppext315 (csext315, 11250);
    StopNeutralinoMap.insert(std::make_pair(mpext315, mppext315));

    MassPoint mpext316 (1125, 400);
    StopCrossSection csext316 (0.000301119, 4.14674e-05);
    MassPointParameters mppext316 (csext316, 9693);
    StopNeutralinoMap.insert(std::make_pair(mpext316, mppext316));

    MassPoint mpext317 (1125, 425);
    StopCrossSection csext317 (0.000301119, 4.14674e-05);
    MassPointParameters mppext317 (csext317, 8573);
    StopNeutralinoMap.insert(std::make_pair(mpext317, mppext317));

    MassPoint mpext318 (1125, 450);
    StopCrossSection csext318 (0.000301119, 4.14674e-05);
    MassPointParameters mppext318 (csext318, 10164);
    StopNeutralinoMap.insert(std::make_pair(mpext318, mppext318));

    MassPoint mpext319 (1125, 475);
    StopCrossSection csext319 (0.000301119, 4.14674e-05);
    MassPointParameters mppext319 (csext319, 8891);
    StopNeutralinoMap.insert(std::make_pair(mpext319, mppext319));

    MassPoint mpext320 (1125, 500);
    StopCrossSection csext320 (0.000301119, 4.14674e-05);
    MassPointParameters mppext320 (csext320, 9660);
    StopNeutralinoMap.insert(std::make_pair(mpext320, mppext320));

    MassPoint mpext321 (1125, 525);
    StopCrossSection csext321 (0.000301119, 4.14674e-05);
    MassPointParameters mppext321 (csext321, 9195);
    StopNeutralinoMap.insert(std::make_pair(mpext321, mppext321));

    MassPoint mpext322 (1125, 550);
    StopCrossSection csext322 (0.000301119, 4.14674e-05);
    MassPointParameters mppext322 (csext322, 10301);
    StopNeutralinoMap.insert(std::make_pair(mpext322, mppext322));

    MassPoint mpext323 (1125, 575);
    StopCrossSection csext323 (0.000301119, 4.14674e-05);
    MassPointParameters mppext323 (csext323, 8300);
    StopNeutralinoMap.insert(std::make_pair(mpext323, mppext323));

    MassPoint mpext324 (1125, 600);
    StopCrossSection csext324 (0.000301119, 4.14674e-05);
    MassPointParameters mppext324 (csext324, 8860);
    StopNeutralinoMap.insert(std::make_pair(mpext324, mppext324));

    MassPoint mpext325 (1150, 1);
    StopCrossSection csext325 (0.000262408, 3.73521e-05);
    MassPointParameters mppext325 (csext325, 9244);
    StopNeutralinoMap.insert(std::make_pair(mpext325, mppext325));

    MassPoint mpext326 (1150, 25);
    StopCrossSection csext326 (0.000262408, 3.73521e-05);
    MassPointParameters mppext326 (csext326, 10277);
    StopNeutralinoMap.insert(std::make_pair(mpext326, mppext326));

    MassPoint mpext327 (1150, 50);
    StopCrossSection csext327 (0.000262408, 3.73521e-05);
    MassPointParameters mppext327 (csext327, 5603);
    StopNeutralinoMap.insert(std::make_pair(mpext327, mppext327));

    MassPoint mpext328 (1150, 75);
    StopCrossSection csext328 (0.000262408, 3.73521e-05);
    MassPointParameters mppext328 (csext328, 8460);
    StopNeutralinoMap.insert(std::make_pair(mpext328, mppext328));

    MassPoint mpext329 (1150, 100);
    StopCrossSection csext329 (0.000262408, 3.73521e-05);
    MassPointParameters mppext329 (csext329, 9455);
    StopNeutralinoMap.insert(std::make_pair(mpext329, mppext329));

    MassPoint mpext330 (1150, 125);
    StopCrossSection csext330 (0.000262408, 3.73521e-05);
    MassPointParameters mppext330 (csext330, 9331);
    StopNeutralinoMap.insert(std::make_pair(mpext330, mppext330));

    MassPoint mpext331 (1150, 150);
    StopCrossSection csext331 (0.000262408, 3.73521e-05);
    MassPointParameters mppext331 (csext331, 9056);
    StopNeutralinoMap.insert(std::make_pair(mpext331, mppext331));

    MassPoint mpext332 (1150, 175);
    StopCrossSection csext332 (0.000262408, 3.73521e-05);
    MassPointParameters mppext332 (csext332, 7805);
    StopNeutralinoMap.insert(std::make_pair(mpext332, mppext332));

    MassPoint mpext333 (1150, 200);
    StopCrossSection csext333 (0.000262408, 3.73521e-05);
    MassPointParameters mppext333 (csext333, 8412);
    StopNeutralinoMap.insert(std::make_pair(mpext333, mppext333));

    MassPoint mpext334 (1150, 225);
    StopCrossSection csext334 (0.000262408, 3.73521e-05);
    MassPointParameters mppext334 (csext334, 8902);
    StopNeutralinoMap.insert(std::make_pair(mpext334, mppext334));

    MassPoint mpext335 (1150, 250);
    StopCrossSection csext335 (0.000262408, 3.73521e-05);
    MassPointParameters mppext335 (csext335, 9656);
    StopNeutralinoMap.insert(std::make_pair(mpext335, mppext335));

    MassPoint mpext336 (1150, 275);
    StopCrossSection csext336 (0.000262408, 3.73521e-05);
    MassPointParameters mppext336 (csext336, 9682);
    StopNeutralinoMap.insert(std::make_pair(mpext336, mppext336));

    MassPoint mpext337 (1150, 300);
    StopCrossSection csext337 (0.000262408, 3.73521e-05);
    MassPointParameters mppext337 (csext337, 10662);
    StopNeutralinoMap.insert(std::make_pair(mpext337, mppext337));

    MassPoint mpext338 (1150, 325);
    StopCrossSection csext338 (0.000262408, 3.73521e-05);
    MassPointParameters mppext338 (csext338, 8308);
    StopNeutralinoMap.insert(std::make_pair(mpext338, mppext338));

    MassPoint mpext339 (1150, 350);
    StopCrossSection csext339 (0.000262408, 3.73521e-05);
    MassPointParameters mppext339 (csext339, 9554);
    StopNeutralinoMap.insert(std::make_pair(mpext339, mppext339));

    MassPoint mpext340 (1150, 375);
    StopCrossSection csext340 (0.000262408, 3.73521e-05);
    MassPointParameters mppext340 (csext340, 8600);
    StopNeutralinoMap.insert(std::make_pair(mpext340, mppext340));

    MassPoint mpext341 (1150, 400);
    StopCrossSection csext341 (0.000262408, 3.73521e-05);
    MassPointParameters mppext341 (csext341, 10462);
    StopNeutralinoMap.insert(std::make_pair(mpext341, mppext341));

    MassPoint mpext342 (1150, 425);
    StopCrossSection csext342 (0.000262408, 3.73521e-05);
    MassPointParameters mppext342 (csext342, 9871);
    StopNeutralinoMap.insert(std::make_pair(mpext342, mppext342));

    MassPoint mpext343 (1150, 450);
    StopCrossSection csext343 (0.000262408, 3.73521e-05);
    MassPointParameters mppext343 (csext343, 9131);
    StopNeutralinoMap.insert(std::make_pair(mpext343, mppext343));

    MassPoint mpext344 (1150, 475);
    StopCrossSection csext344 (0.000262408, 3.73521e-05);
    MassPointParameters mppext344 (csext344, 11746);
    StopNeutralinoMap.insert(std::make_pair(mpext344, mppext344));

    MassPoint mpext345 (1150, 500);
    StopCrossSection csext345 (0.000262408, 3.73521e-05);
    MassPointParameters mppext345 (csext345, 7617);
    StopNeutralinoMap.insert(std::make_pair(mpext345, mppext345));

    MassPoint mpext346 (1150, 525);
    StopCrossSection csext346 (0.000262408, 3.73521e-05);
    MassPointParameters mppext346 (csext346, 11239);
    StopNeutralinoMap.insert(std::make_pair(mpext346, mppext346));

    MassPoint mpext347 (1150, 550);
    StopCrossSection csext347 (0.000262408, 3.73521e-05);
    MassPointParameters mppext347 (csext347, 6853);
    StopNeutralinoMap.insert(std::make_pair(mpext347, mppext347));

    MassPoint mpext348 (1150, 575);
    StopCrossSection csext348 (0.000262408, 3.73521e-05);
    MassPointParameters mppext348 (csext348, 10230);
    StopNeutralinoMap.insert(std::make_pair(mpext348, mppext348));

    MassPoint mpext349 (1150, 600);
    StopCrossSection csext349 (0.000262408, 3.73521e-05);
    MassPointParameters mppext349 (csext349, 7760);
    StopNeutralinoMap.insert(std::make_pair(mpext349, mppext349));

    MassPoint mpext350 (1175, 1);
    StopCrossSection csext350 (0.000224723, 3.01438e-05);
    MassPointParameters mppext350 (csext350, 6954);
    StopNeutralinoMap.insert(std::make_pair(mpext350, mppext350));

    MassPoint mpext351 (1175, 25);
    StopCrossSection csext351 (0.000224723, 3.01438e-05);
    MassPointParameters mppext351 (csext351, 10079);
    StopNeutralinoMap.insert(std::make_pair(mpext351, mppext351));

    MassPoint mpext352 (1175, 50);
    StopCrossSection csext352 (0.000224723, 3.01438e-05);
    MassPointParameters mppext352 (csext352, 9567);
    StopNeutralinoMap.insert(std::make_pair(mpext352, mppext352));

    MassPoint mpext353 (1175, 75);
    StopCrossSection csext353 (0.000224723, 3.01438e-05);
    MassPointParameters mppext353 (csext353, 6447);
    StopNeutralinoMap.insert(std::make_pair(mpext353, mppext353));

    MassPoint mpext354 (1175, 100);
    StopCrossSection csext354 (0.000224723, 3.01438e-05);
    MassPointParameters mppext354 (csext354, 9620);
    StopNeutralinoMap.insert(std::make_pair(mpext354, mppext354));

    MassPoint mpext355 (1175, 125);
    StopCrossSection csext355 (0.000224723, 3.01438e-05);
    MassPointParameters mppext355 (csext355, 11746);
    StopNeutralinoMap.insert(std::make_pair(mpext355, mppext355));

    MassPoint mpext356 (1175, 150);
    StopCrossSection csext356 (0.000224723, 3.01438e-05);
    MassPointParameters mppext356 (csext356, 8553);
    StopNeutralinoMap.insert(std::make_pair(mpext356, mppext356));

    MassPoint mpext357 (1175, 175);
    StopCrossSection csext357 (0.000224723, 3.01438e-05);
    MassPointParameters mppext357 (csext357, 8201);
    StopNeutralinoMap.insert(std::make_pair(mpext357, mppext357));

    MassPoint mpext358 (1175, 200);
    StopCrossSection csext358 (0.000224723, 3.01438e-05);
    MassPointParameters mppext358 (csext358, 11981);
    StopNeutralinoMap.insert(std::make_pair(mpext358, mppext358));

    MassPoint mpext359 (1175, 225);
    StopCrossSection csext359 (0.000224723, 3.01438e-05);
    MassPointParameters mppext359 (csext359, 9636);
    StopNeutralinoMap.insert(std::make_pair(mpext359, mppext359));

    MassPoint mpext360 (1175, 250);
    StopCrossSection csext360 (0.000224723, 3.01438e-05);
    MassPointParameters mppext360 (csext360, 9831);
    StopNeutralinoMap.insert(std::make_pair(mpext360, mppext360));

    MassPoint mpext361 (1175, 275);
    StopCrossSection csext361 (0.000224723, 3.01438e-05);
    MassPointParameters mppext361 (csext361, 9052);
    StopNeutralinoMap.insert(std::make_pair(mpext361, mppext361));

    MassPoint mpext362 (1175, 300);
    StopCrossSection csext362 (0.000224723, 3.01438e-05);
    MassPointParameters mppext362 (csext362, 8645);
    StopNeutralinoMap.insert(std::make_pair(mpext362, mppext362));

    MassPoint mpext363 (1175, 325);
    StopCrossSection csext363 (0.000224723, 3.01438e-05);
    MassPointParameters mppext363 (csext363, 8300);
    StopNeutralinoMap.insert(std::make_pair(mpext363, mppext363));

    MassPoint mpext364 (1175, 350);
    StopCrossSection csext364 (0.000224723, 3.01438e-05);
    MassPointParameters mppext364 (csext364, 7917);
    StopNeutralinoMap.insert(std::make_pair(mpext364, mppext364));

    MassPoint mpext365 (1175, 375);
    StopCrossSection csext365 (0.000224723, 3.01438e-05);
    MassPointParameters mppext365 (csext365, 14331);
    StopNeutralinoMap.insert(std::make_pair(mpext365, mppext365));

    MassPoint mpext366 (1175, 400);
    StopCrossSection csext366 (0.000224723, 3.01438e-05);
    MassPointParameters mppext366 (csext366, 9212);
    StopNeutralinoMap.insert(std::make_pair(mpext366, mppext366));

    MassPoint mpext367 (1175, 425);
    StopCrossSection csext367 (0.000224723, 3.01438e-05);
    MassPointParameters mppext367 (csext367, 9760);
    StopNeutralinoMap.insert(std::make_pair(mpext367, mppext367));

    MassPoint mpext368 (1175, 450);
    StopCrossSection csext368 (0.000224723, 3.01438e-05);
    MassPointParameters mppext368 (csext368, 9937);
    StopNeutralinoMap.insert(std::make_pair(mpext368, mppext368));

    MassPoint mpext369 (1175, 475);
    StopCrossSection csext369 (0.000224723, 3.01438e-05);
    MassPointParameters mppext369 (csext369, 12153);
    StopNeutralinoMap.insert(std::make_pair(mpext369, mppext369));

    MassPoint mpext370 (1175, 500);
    StopCrossSection csext370 (0.000224723, 3.01438e-05);
    MassPointParameters mppext370 (csext370, 8611);
    StopNeutralinoMap.insert(std::make_pair(mpext370, mppext370));

    MassPoint mpext371 (1175, 525);
    StopCrossSection csext371 (0.000224723, 3.01438e-05);
    MassPointParameters mppext371 (csext371, 11344);
    StopNeutralinoMap.insert(std::make_pair(mpext371, mppext371));

    MassPoint mpext372 (1175, 550);
    StopCrossSection csext372 (0.000224723, 3.01438e-05);
    MassPointParameters mppext372 (csext372, 8382);
    StopNeutralinoMap.insert(std::make_pair(mpext372, mppext372));

    MassPoint mpext373 (1175, 575);
    StopCrossSection csext373 (0.000224723, 3.01438e-05);
    MassPointParameters mppext373 (csext373, 8075);
    StopNeutralinoMap.insert(std::make_pair(mpext373, mppext373));

    MassPoint mpext374 (1175, 600);
    StopCrossSection csext374 (0.000224723, 3.01438e-05);
    MassPointParameters mppext374 (csext374, 9580);
    StopNeutralinoMap.insert(std::make_pair(mpext374, mppext374));

    MassPoint mpext375 (1200, 1);
    StopCrossSection csext375 (0.000196044, 2.64135e-05);
    MassPointParameters mppext375 (csext375, 7451);
    StopNeutralinoMap.insert(std::make_pair(mpext375, mppext375));

    MassPoint mpext376 (1200, 25);
    StopCrossSection csext376 (0.000196044, 2.64135e-05);
    MassPointParameters mppext376 (csext376, 7846);
    StopNeutralinoMap.insert(std::make_pair(mpext376, mppext376));

    MassPoint mpext377 (1200, 50);
    StopCrossSection csext377 (0.000196044, 2.64135e-05);
    MassPointParameters mppext377 (csext377, 10689);
    StopNeutralinoMap.insert(std::make_pair(mpext377, mppext377));

    MassPoint mpext378 (1200, 75);
    StopCrossSection csext378 (0.000196044, 2.64135e-05);
    MassPointParameters mppext378 (csext378, 8714);
    StopNeutralinoMap.insert(std::make_pair(mpext378, mppext378));

    MassPoint mpext379 (1200, 100);
    StopCrossSection csext379 (0.000196044, 2.64135e-05);
    MassPointParameters mppext379 (csext379, 7788);
    StopNeutralinoMap.insert(std::make_pair(mpext379, mppext379));

    MassPoint mpext380 (1200, 125);
    StopCrossSection csext380 (0.000196044, 2.64135e-05);
    MassPointParameters mppext380 (csext380, 9477);
    StopNeutralinoMap.insert(std::make_pair(mpext380, mppext380));

    MassPoint mpext381 (1200, 150);
    StopCrossSection csext381 (0.000196044, 2.64135e-05);
    MassPointParameters mppext381 (csext381, 8370);
    StopNeutralinoMap.insert(std::make_pair(mpext381, mppext381));

    MassPoint mpext382 (1200, 175);
    StopCrossSection csext382 (0.000196044, 2.64135e-05);
    MassPointParameters mppext382 (csext382, 8930);
    StopNeutralinoMap.insert(std::make_pair(mpext382, mppext382));

    MassPoint mpext383 (1200, 200);
    StopCrossSection csext383 (0.000196044, 2.64135e-05);
    MassPointParameters mppext383 (csext383, 9465);
    StopNeutralinoMap.insert(std::make_pair(mpext383, mppext383));

    MassPoint mpext384 (1200, 225);
    StopCrossSection csext384 (0.000196044, 2.64135e-05);
    MassPointParameters mppext384 (csext384, 8879);
    StopNeutralinoMap.insert(std::make_pair(mpext384, mppext384));

    MassPoint mpext385 (1200, 250);
    StopCrossSection csext385 (0.000196044, 2.64135e-05);
    MassPointParameters mppext385 (csext385, 10140);
    StopNeutralinoMap.insert(std::make_pair(mpext385, mppext385));

    MassPoint mpext386 (1200, 275);
    StopCrossSection csext386 (0.000196044, 2.64135e-05);
    MassPointParameters mppext386 (csext386, 8521);
    StopNeutralinoMap.insert(std::make_pair(mpext386, mppext386));

    MassPoint mpext387 (1200, 300);
    StopCrossSection csext387 (0.000196044, 2.64135e-05);
    MassPointParameters mppext387 (csext387, 10449);
    StopNeutralinoMap.insert(std::make_pair(mpext387, mppext387));

    MassPoint mpext388 (1200, 325);
    StopCrossSection csext388 (0.000196044, 2.64135e-05);
    MassPointParameters mppext388 (csext388, 13930);
    StopNeutralinoMap.insert(std::make_pair(mpext388, mppext388));

    MassPoint mpext389 (1200, 350);
    StopCrossSection csext389 (0.000196044, 2.64135e-05);
    MassPointParameters mppext389 (csext389, 9972);
    StopNeutralinoMap.insert(std::make_pair(mpext389, mppext389));

    MassPoint mpext390 (1200, 375);
    StopCrossSection csext390 (0.000196044, 2.64135e-05);
    MassPointParameters mppext390 (csext390, 9681);
    StopNeutralinoMap.insert(std::make_pair(mpext390, mppext390));

    MassPoint mpext391 (1200, 400);
    StopCrossSection csext391 (0.000196044, 2.64135e-05);
    MassPointParameters mppext391 (csext391, 10861);
    StopNeutralinoMap.insert(std::make_pair(mpext391, mppext391));

    MassPoint mpext392 (1200, 425);
    StopCrossSection csext392 (0.000196044, 2.64135e-05);
    MassPointParameters mppext392 (csext392, 9773);
    StopNeutralinoMap.insert(std::make_pair(mpext392, mppext392));

    MassPoint mpext393 (1200, 450);
    StopCrossSection csext393 (0.000196044, 2.64135e-05);
    MassPointParameters mppext393 (csext393, 13042);
    StopNeutralinoMap.insert(std::make_pair(mpext393, mppext393));

    MassPoint mpext394 (1200, 475);
    StopCrossSection csext394 (0.000196044, 2.64135e-05);
    MassPointParameters mppext394 (csext394, 8334);
    StopNeutralinoMap.insert(std::make_pair(mpext394, mppext394));

    MassPoint mpext395 (1200, 500);
    StopCrossSection csext395 (0.000196044, 2.64135e-05);
    MassPointParameters mppext395 (csext395, 11570);
    StopNeutralinoMap.insert(std::make_pair(mpext395, mppext395));

    MassPoint mpext396 (1200, 525);
    StopCrossSection csext396 (0.000196044, 2.64135e-05);
    MassPointParameters mppext396 (csext396, 9849);
    StopNeutralinoMap.insert(std::make_pair(mpext396, mppext396));

    MassPoint mpext397 (1200, 550);
    StopCrossSection csext397 (0.000196044, 2.64135e-05);
    MassPointParameters mppext397 (csext397, 10377);
    StopNeutralinoMap.insert(std::make_pair(mpext397, mppext397));

    MassPoint mpext398 (1200, 575);
    StopCrossSection csext398 (0.000196044, 2.64135e-05);
    MassPointParameters mppext398 (csext398, 8202);
    StopNeutralinoMap.insert(std::make_pair(mpext398, mppext398));

    MassPoint mpext399 (1200, 600);
    StopCrossSection csext399 (0.000196044, 2.64135e-05);
    MassPointParameters mppext399 (csext399, 9880);
    StopNeutralinoMap.insert(std::make_pair(mpext399, mppext399));

    MassPoint mpisr0 (100, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr0, 0.979732));

    MassPoint mpisr1 (100, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr1, 0.97921));

    MassPoint mpisr2 (100, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr2, 0.979101));

    MassPoint mpisr3 (125, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr3, 0.976387));

    MassPoint mpisr4 (125, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr4, 0.976411));

    MassPoint mpisr5 (125, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr5, 0.976059));

    MassPoint mpisr6 (125, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr6, 0.976566));

    MassPoint mpisr7 (150, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr7, 0.97391));

    MassPoint mpisr8 (150, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr8, 0.974946));

    MassPoint mpisr9 (150, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr9, 0.973733));

    MassPoint mpisr10 (150, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr10, 0.973955));

    MassPoint mpisr11 (150, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr11, 0.973552));

    MassPoint mpisr12 (175, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr12, 0.97202));

    MassPoint mpisr13 (175, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr13, 0.971457));

    MassPoint mpisr14 (175, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr14, 0.972362));

    MassPoint mpisr15 (175, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr15, 0.972942));

    MassPoint mpisr16 (175, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr16, 0.97182));

    MassPoint mpisr17 (175, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr17, 0.97172));

    MassPoint mpisr18 (200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr18, 0.969949));

    MassPoint mpisr19 (200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr19, 0.970696));

    MassPoint mpisr20 (200, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr20, 0.970238));

    MassPoint mpisr21 (200, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr21, 0.970672));

    MassPoint mpisr22 (200, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr22, 0.971308));

    MassPoint mpisr23 (200, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr23, 0.971257));

    MassPoint mpisr24 (200, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr24, 0.970811));

    MassPoint mpisr25 (225, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr25, 0.969733));

    MassPoint mpisr26 (225, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr26, 0.969674));

    MassPoint mpisr27 (225, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr27, 0.970271));

    MassPoint mpisr28 (225, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr28, 0.96999));

    MassPoint mpisr29 (225, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr29, 0.969394));

    MassPoint mpisr30 (225, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr30, 0.969561));

    MassPoint mpisr31 (225, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr31, 0.970424));

    MassPoint mpisr32 (225, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr32, 0.968637));

    MassPoint mpisr33 (250, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr33, 0.968804));

    MassPoint mpisr34 (250, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr34, 0.969455));

    MassPoint mpisr35 (250, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr35, 0.968801));

    MassPoint mpisr36 (250, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr36, 0.969198));

    MassPoint mpisr37 (250, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr37, 0.969145));

    MassPoint mpisr38 (250, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr38, 0.968993));

    MassPoint mpisr39 (250, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr39, 0.969479));

    MassPoint mpisr40 (250, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr40, 0.969686));

    MassPoint mpisr41 (250, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr41, 0.967571));

    MassPoint mpisr42 (275, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr42, 0.967557));

    MassPoint mpisr43 (275, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr43, 0.968142));

    MassPoint mpisr44 (275, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr44, 0.967866));

    MassPoint mpisr45 (275, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr45, 0.96758));

    MassPoint mpisr46 (275, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr46, 0.969261));

    MassPoint mpisr47 (275, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr47, 0.967689));

    MassPoint mpisr48 (275, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr48, 0.967814));

    MassPoint mpisr49 (275, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr49, 0.969103));

    MassPoint mpisr50 (275, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr50, 0.967858));

    MassPoint mpisr51 (275, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr51, 0.968218));

    MassPoint mpisr52 (300, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr52, 0.966528));

    MassPoint mpisr53 (300, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr53, 0.968347));

    MassPoint mpisr54 (300, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr54, 0.96711));

    MassPoint mpisr55 (300, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr55, 0.968174));

    MassPoint mpisr56 (300, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr56, 0.967687));

    MassPoint mpisr57 (300, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr57, 0.966821));

    MassPoint mpisr58 (300, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr58, 0.96721));

    MassPoint mpisr59 (300, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr59, 0.967399));

    MassPoint mpisr60 (300, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr60, 0.967778));

    MassPoint mpisr61 (300, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr61, 0.967579));

    MassPoint mpisr62 (300, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr62, 0.968116));

    MassPoint mpisr63 (325, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr63, 0.967408));

    MassPoint mpisr64 (325, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr64, 0.966851));

    MassPoint mpisr65 (325, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr65, 0.96686));

    MassPoint mpisr66 (325, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr66, 0.965345));

    MassPoint mpisr67 (325, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr67, 0.96592));

    MassPoint mpisr68 (325, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr68, 0.967107));

    MassPoint mpisr69 (325, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr69, 0.965912));

    MassPoint mpisr70 (325, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr70, 0.967432));

    MassPoint mpisr71 (325, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr71, 0.967501));

    MassPoint mpisr72 (325, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr72, 0.966884));

    MassPoint mpisr73 (325, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr73, 0.966691));

    MassPoint mpisr74 (325, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr74, 0.96622));

    MassPoint mpisr75 (350, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr75, 0.966781));

    MassPoint mpisr76 (350, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr76, 0.965287));

    MassPoint mpisr77 (350, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr77, 0.965859));

    MassPoint mpisr78 (350, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr78, 0.966167));

    MassPoint mpisr79 (350, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr79, 0.966384));

    MassPoint mpisr80 (350, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr80, 0.966453));

    MassPoint mpisr81 (350, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr81, 0.966409));

    MassPoint mpisr82 (350, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr82, 0.966291));

    MassPoint mpisr83 (350, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr83, 0.966055));

    MassPoint mpisr84 (350, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr84, 0.966712));

    MassPoint mpisr85 (350, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr85, 0.965993));

    MassPoint mpisr86 (350, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr86, 0.967239));

    MassPoint mpisr87 (350, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr87, 0.967529));

    MassPoint mpisr88 (375, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr88, 0.967166));

    MassPoint mpisr89 (375, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr89, 0.966793));

    MassPoint mpisr90 (375, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr90, 0.966273));

    MassPoint mpisr91 (375, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr91, 0.96529));

    MassPoint mpisr92 (375, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr92, 0.967105));

    MassPoint mpisr93 (375, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr93, 0.965839));

    MassPoint mpisr94 (375, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr94, 0.966126));

    MassPoint mpisr95 (375, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr95, 0.966589));

    MassPoint mpisr96 (375, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr96, 0.966867));

    MassPoint mpisr97 (375, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr97, 0.966283));

    MassPoint mpisr98 (375, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr98, 0.966474));

    MassPoint mpisr99 (375, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr99, 0.967416));

    MassPoint mpisr100 (375, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr100, 0.965906));

    MassPoint mpisr101 (375, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr101, 0.966608));

    MassPoint mpisr102 (400, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr102, 0.966193));

    MassPoint mpisr103 (400, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr103, 0.965954));

    MassPoint mpisr104 (400, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr104, 0.965839));

    MassPoint mpisr105 (400, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr105, 0.966196));

    MassPoint mpisr106 (400, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr106, 0.967759));

    MassPoint mpisr107 (400, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr107, 0.96628));

    MassPoint mpisr108 (400, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr108, 0.965925));

    MassPoint mpisr109 (400, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr109, 0.966051));

    MassPoint mpisr110 (400, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr110, 0.966713));

    MassPoint mpisr111 (400, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr111, 0.966134));

    MassPoint mpisr112 (400, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr112, 0.965272));

    MassPoint mpisr113 (400, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr113, 0.966239));

    MassPoint mpisr114 (400, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr114, 0.966955));

    MassPoint mpisr115 (400, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr115, 0.965555));

    MassPoint mpisr116 (400, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr116, 0.966019));

    MassPoint mpisr117 (425, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr117, 0.966237));

    MassPoint mpisr118 (425, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr118, 0.966054));

    MassPoint mpisr119 (425, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr119, 0.966259));

    MassPoint mpisr120 (425, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr120, 0.966657));

    MassPoint mpisr121 (425, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr121, 0.966159));

    MassPoint mpisr122 (425, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr122, 0.965823));

    MassPoint mpisr123 (425, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr123, 0.967183));

    MassPoint mpisr124 (425, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr124, 0.965881));

    MassPoint mpisr125 (425, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr125, 0.966936));

    MassPoint mpisr126 (425, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr126, 0.966533));

    MassPoint mpisr127 (425, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr127, 0.966635));

    MassPoint mpisr128 (425, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr128, 0.966952));

    MassPoint mpisr129 (425, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr129, 0.966437));

    MassPoint mpisr130 (425, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr130, 0.966806));

    MassPoint mpisr131 (425, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr131, 0.965518));

    MassPoint mpisr132 (425, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr132, 0.965477));

    MassPoint mpisr133 (450, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr133, 0.964805));

    MassPoint mpisr134 (450, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr134, 0.965336));

    MassPoint mpisr135 (450, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr135, 0.96563));

    MassPoint mpisr136 (450, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr136, 0.965009));

    MassPoint mpisr137 (450, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr137, 0.965429));

    MassPoint mpisr138 (450, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr138, 0.965552));

    MassPoint mpisr139 (450, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr139, 0.96698));

    MassPoint mpisr140 (450, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr140, 0.964996));

    MassPoint mpisr141 (450, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr141, 0.965742));

    MassPoint mpisr142 (450, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr142, 0.965999));

    MassPoint mpisr143 (450, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr143, 0.966345));

    MassPoint mpisr144 (450, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr144, 0.965505));

    MassPoint mpisr145 (450, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr145, 0.964646));

    MassPoint mpisr146 (450, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr146, 0.966635));

    MassPoint mpisr147 (450, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr147, 0.96495));

    MassPoint mpisr148 (450, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr148, 0.964396));

    MassPoint mpisr149 (450, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr149, 0.966074));

    MassPoint mpisr150 (475, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr150, 0.966773));

    MassPoint mpisr151 (475, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr151, 0.965204));

    MassPoint mpisr152 (475, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr152, 0.965927));

    MassPoint mpisr153 (475, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr153, 0.965118));

    MassPoint mpisr154 (475, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr154, 0.966781));

    MassPoint mpisr155 (475, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr155, 0.96603));

    MassPoint mpisr156 (475, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr156, 0.965002));

    MassPoint mpisr157 (475, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr157, 0.964981));

    MassPoint mpisr158 (475, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr158, 0.966478));

    MassPoint mpisr159 (475, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr159, 0.96568));

    MassPoint mpisr160 (475, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr160, 0.965834));

    MassPoint mpisr161 (475, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr161, 0.965609));

    MassPoint mpisr162 (475, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr162, 0.965802));

    MassPoint mpisr163 (475, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr163, 0.964246));

    MassPoint mpisr164 (475, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr164, 0.964456));

    MassPoint mpisr165 (475, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr165, 0.965196));

    MassPoint mpisr166 (475, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr166, 0.965315));

    MassPoint mpisr167 (475, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr167, 0.965428));

    MassPoint mpisr168 (500, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr168, 0.965528));

    MassPoint mpisr169 (500, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr169, 0.964423));

    MassPoint mpisr170 (500, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr170, 0.965282));

    MassPoint mpisr171 (500, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr171, 0.96558));

    MassPoint mpisr172 (500, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr172, 0.964916));

    MassPoint mpisr173 (500, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr173, 0.965632));

    MassPoint mpisr174 (500, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr174, 0.965097));

    MassPoint mpisr175 (500, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr175, 0.96597));

    MassPoint mpisr176 (500, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr176, 0.96541));

    MassPoint mpisr177 (500, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr177, 0.964333));

    MassPoint mpisr178 (500, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr178, 0.964799));

    MassPoint mpisr179 (500, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr179, 0.965111));

    MassPoint mpisr180 (500, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr180, 0.965471));

    MassPoint mpisr181 (500, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr181, 0.965379));

    MassPoint mpisr182 (500, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr182, 0.965677));

    MassPoint mpisr183 (500, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr183, 0.964702));

    MassPoint mpisr184 (500, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr184, 0.963814));

    MassPoint mpisr185 (500, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr185, 0.965242));

    MassPoint mpisr186 (500, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr186, 0.965224));

    MassPoint mpisr187 (525, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr187, 0.965658));

    MassPoint mpisr188 (525, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr188, 0.964847));

    MassPoint mpisr189 (525, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr189, 0.964594));

    MassPoint mpisr190 (525, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr190, 0.964079));

    MassPoint mpisr191 (525, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr191, 0.964582));

    MassPoint mpisr192 (525, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr192, 0.963982));

    MassPoint mpisr193 (525, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr193, 0.965234));

    MassPoint mpisr194 (525, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr194, 0.965107));

    MassPoint mpisr195 (525, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr195, 0.965036));

    MassPoint mpisr196 (525, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr196, 0.964602));

    MassPoint mpisr197 (525, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr197, 0.964998));

    MassPoint mpisr198 (525, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr198, 0.964931));

    MassPoint mpisr199 (525, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr199, 0.965836));

    MassPoint mpisr200 (525, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr200, 0.966103));

    MassPoint mpisr201 (525, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr201, 0.966047));

    MassPoint mpisr202 (525, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr202, 0.964918));

    MassPoint mpisr203 (525, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr203, 0.964905));

    MassPoint mpisr204 (525, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr204, 0.964714));

    MassPoint mpisr205 (525, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr205, 0.964311));

    MassPoint mpisr206 (525, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr206, 0.964525));

    MassPoint mpisr207 (550, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr207, 0.965304));

    MassPoint mpisr208 (550, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr208, 0.96547));

    MassPoint mpisr209 (550, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr209, 0.965482));

    MassPoint mpisr210 (550, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr210, 0.964175));

    MassPoint mpisr211 (550, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr211, 0.96447));

    MassPoint mpisr212 (550, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr212, 0.964959));

    MassPoint mpisr213 (550, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr213, 0.965553));

    MassPoint mpisr214 (550, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr214, 0.964602));

    MassPoint mpisr215 (550, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr215, 0.966012));

    MassPoint mpisr216 (550, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr216, 0.96442));

    MassPoint mpisr217 (550, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr217, 0.966314));

    MassPoint mpisr218 (550, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr218, 0.964981));

    MassPoint mpisr219 (550, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr219, 0.965322));

    MassPoint mpisr220 (550, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr220, 0.9648));

    MassPoint mpisr221 (550, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr221, 0.964723));

    MassPoint mpisr222 (550, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr222, 0.965866));

    MassPoint mpisr223 (550, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr223, 0.966124));

    MassPoint mpisr224 (550, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr224, 0.965412));

    MassPoint mpisr225 (550, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr225, 0.965462));

    MassPoint mpisr226 (550, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr226, 0.965387));

    MassPoint mpisr227 (550, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr227, 0.964753));

    MassPoint mpisr228 (575, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr228, 0.964589));

    MassPoint mpisr229 (575, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr229, 0.964727));

    MassPoint mpisr230 (575, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr230, 0.9646));

    MassPoint mpisr231 (575, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr231, 0.965327));

    MassPoint mpisr232 (575, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr232, 0.963837));

    MassPoint mpisr233 (575, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr233, 0.96445));

    MassPoint mpisr234 (575, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr234, 0.964279));

    MassPoint mpisr235 (575, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr235, 0.964135));

    MassPoint mpisr236 (575, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr236, 0.964975));

    MassPoint mpisr237 (575, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr237, 0.965448));

    MassPoint mpisr238 (575, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr238, 0.965566));

    MassPoint mpisr239 (575, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr239, 0.964998));

    MassPoint mpisr240 (575, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr240, 0.964165));

    MassPoint mpisr241 (575, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr241, 0.966023));

    MassPoint mpisr242 (575, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr242, 0.964217));

    MassPoint mpisr243 (575, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr243, 0.965834));

    MassPoint mpisr244 (575, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr244, 0.964327));

    MassPoint mpisr245 (575, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr245, 0.965379));

    MassPoint mpisr246 (575, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr246, 0.965718));

    MassPoint mpisr247 (575, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr247, 0.963941));

    MassPoint mpisr248 (575, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr248, 0.96424));

    MassPoint mpisr249 (575, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr249, 0.965031));

    MassPoint mpisr250 (600, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr250, 0.965973));

    MassPoint mpisr251 (600, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr251, 0.964899));

    MassPoint mpisr252 (600, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr252, 0.964849));

    MassPoint mpisr253 (600, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr253, 0.965403));

    MassPoint mpisr254 (600, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr254, 0.964317));

    MassPoint mpisr255 (600, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr255, 0.96425));

    MassPoint mpisr256 (600, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr256, 0.963264));

    MassPoint mpisr257 (600, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr257, 0.964851));

    MassPoint mpisr258 (600, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr258, 0.964569));

    MassPoint mpisr259 (600, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr259, 0.964335));

    MassPoint mpisr260 (600, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr260, 0.965488));

    MassPoint mpisr261 (600, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr261, 0.964353));

    MassPoint mpisr262 (600, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr262, 0.965092));

    MassPoint mpisr263 (600, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr263, 0.963705));

    MassPoint mpisr264 (600, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr264, 0.964286));

    MassPoint mpisr265 (600, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr265, 0.965379));

    MassPoint mpisr266 (600, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr266, 0.965161));

    MassPoint mpisr267 (600, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr267, 0.965468));

    MassPoint mpisr268 (600, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr268, 0.965461));

    MassPoint mpisr269 (600, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr269, 0.964256));

    MassPoint mpisr270 (600, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr270, 0.966117));

    MassPoint mpisr271 (600, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr271, 0.966057));

    MassPoint mpisr272 (600, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr272, 0.965285));

    MassPoint mpisr273 (625, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr273, 0.964682));

    MassPoint mpisr274 (625, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr274, 0.963868));

    MassPoint mpisr275 (625, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr275, 0.963891));

    MassPoint mpisr276 (625, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr276, 0.964541));

    MassPoint mpisr277 (625, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr277, 0.964405));

    MassPoint mpisr278 (625, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr278, 0.963451));

    MassPoint mpisr279 (625, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr279, 0.964692));

    MassPoint mpisr280 (625, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr280, 0.964816));

    MassPoint mpisr281 (625, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr281, 0.965387));

    MassPoint mpisr282 (625, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr282, 0.965231));

    MassPoint mpisr283 (625, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr283, 0.964616));

    MassPoint mpisr284 (625, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr284, 0.964905));

    MassPoint mpisr285 (625, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr285, 0.963795));

    MassPoint mpisr286 (625, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr286, 0.965371));

    MassPoint mpisr287 (625, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr287, 0.965002));

    MassPoint mpisr288 (625, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr288, 0.963777));

    MassPoint mpisr289 (625, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr289, 0.964263));

    MassPoint mpisr290 (625, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr290, 0.964736));

    MassPoint mpisr291 (625, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr291, 0.963971));

    MassPoint mpisr292 (625, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr292, 0.9647));

    MassPoint mpisr293 (625, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr293, 0.964714));

    MassPoint mpisr294 (625, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr294, 0.965305));

    MassPoint mpisr295 (625, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr295, 0.964294));

    MassPoint mpisr296 (625, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr296, 0.9659));

    MassPoint mpisr297 (650, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr297, 0.964863));

    MassPoint mpisr298 (650, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr298, 0.963717));

    MassPoint mpisr299 (650, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr299, 0.965867));

    MassPoint mpisr300 (650, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr300, 0.964259));

    MassPoint mpisr301 (650, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr301, 0.96464));

    MassPoint mpisr302 (650, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr302, 0.965354));

    MassPoint mpisr303 (650, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr303, 0.963721));

    MassPoint mpisr304 (650, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr304, 0.964559));

    MassPoint mpisr305 (650, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr305, 0.9642));

    MassPoint mpisr306 (650, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr306, 0.964928));

    MassPoint mpisr307 (650, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr307, 0.964733));

    MassPoint mpisr308 (650, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr308, 0.965162));

    MassPoint mpisr309 (650, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr309, 0.965172));

    MassPoint mpisr310 (650, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr310, 0.965384));

    MassPoint mpisr311 (650, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr311, 0.964651));

    MassPoint mpisr312 (650, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr312, 0.963763));

    MassPoint mpisr313 (650, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr313, 0.964397));

    MassPoint mpisr314 (650, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr314, 0.96398));

    MassPoint mpisr315 (650, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr315, 0.965005));

    MassPoint mpisr316 (650, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr316, 0.965902));

    MassPoint mpisr317 (650, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr317, 0.964135));

    MassPoint mpisr318 (650, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr318, 0.96495));

    MassPoint mpisr319 (650, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr319, 0.964656));

    MassPoint mpisr320 (650, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr320, 0.965395));

    MassPoint mpisr321 (650, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr321, 0.964852));

    MassPoint mpisr322 (675, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr322, 0.965053));

    MassPoint mpisr323 (675, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr323, 0.965734));

    MassPoint mpisr324 (675, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr324, 0.964276));

    MassPoint mpisr325 (675, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr325, 0.96523));

    MassPoint mpisr326 (675, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr326, 0.964416));

    MassPoint mpisr327 (675, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr327, 0.96427));

    MassPoint mpisr328 (675, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr328, 0.96426));

    MassPoint mpisr329 (675, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr329, 0.965052));

    MassPoint mpisr330 (675, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr330, 0.964744));

    MassPoint mpisr331 (675, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr331, 0.964125));

    MassPoint mpisr332 (675, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr332, 0.965063));

    MassPoint mpisr333 (675, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr333, 0.964679));

    MassPoint mpisr334 (675, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr334, 0.965037));

    MassPoint mpisr335 (675, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr335, 0.965044));

    MassPoint mpisr336 (675, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr336, 0.965845));

    MassPoint mpisr337 (675, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr337, 0.965246));

    MassPoint mpisr338 (675, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr338, 0.963807));

    MassPoint mpisr339 (675, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr339, 0.963543));

    MassPoint mpisr340 (675, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr340, 0.964851));

    MassPoint mpisr341 (675, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr341, 0.964732));

    MassPoint mpisr342 (675, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr342, 0.965368));

    MassPoint mpisr343 (675, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr343, 0.964949));

    MassPoint mpisr344 (675, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr344, 0.964799));

    MassPoint mpisr345 (675, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr345, 0.964668));

    MassPoint mpisr346 (675, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr346, 0.964112));

    MassPoint mpisr347 (700, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr347, 0.964256));

    MassPoint mpisr348 (700, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr348, 0.96337));

    MassPoint mpisr349 (700, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr349, 0.964184));

    MassPoint mpisr350 (700, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr350, 0.963889));

    MassPoint mpisr351 (700, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr351, 0.963455));

    MassPoint mpisr352 (700, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr352, 0.965058));

    MassPoint mpisr353 (700, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr353, 0.96386));

    MassPoint mpisr354 (700, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr354, 0.963363));

    MassPoint mpisr355 (700, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr355, 0.963489));

    MassPoint mpisr356 (700, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr356, 0.963532));

    MassPoint mpisr357 (700, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr357, 0.96516));

    MassPoint mpisr358 (700, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr358, 0.964621));

    MassPoint mpisr359 (700, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr359, 0.963506));

    MassPoint mpisr360 (700, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr360, 0.964154));

    MassPoint mpisr361 (700, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr361, 0.964722));

    MassPoint mpisr362 (700, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr362, 0.964501));

    MassPoint mpisr363 (700, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr363, 0.964669));

    MassPoint mpisr364 (700, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr364, 0.965406));

    MassPoint mpisr365 (700, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr365, 0.965751));

    MassPoint mpisr366 (700, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr366, 0.964589));

    MassPoint mpisr367 (700, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr367, 0.96558));

    MassPoint mpisr368 (700, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr368, 0.964656));

    MassPoint mpisr369 (700, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr369, 0.965084));

    MassPoint mpisr370 (700, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr370, 0.964111));

    MassPoint mpisr371 (700, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr371, 0.964144));

    MassPoint mpisr372 (725, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr372, 0.96392));

    MassPoint mpisr373 (725, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr373, 0.965065));

    MassPoint mpisr374 (725, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr374, 0.963211));

    MassPoint mpisr375 (725, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr375, 0.96525));

    MassPoint mpisr376 (725, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr376, 0.964491));

    MassPoint mpisr377 (725, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr377, 0.965423));

    MassPoint mpisr378 (725, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr378, 0.964134));

    MassPoint mpisr379 (725, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr379, 0.964583));

    MassPoint mpisr380 (725, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr380, 0.966032));

    MassPoint mpisr381 (725, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr381, 0.964724));

    MassPoint mpisr382 (725, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr382, 0.96427));

    MassPoint mpisr383 (725, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr383, 0.96479));

    MassPoint mpisr384 (725, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr384, 0.96471));

    MassPoint mpisr385 (725, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr385, 0.964257));

    MassPoint mpisr386 (725, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr386, 0.964179));

    MassPoint mpisr387 (725, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr387, 0.964435));

    MassPoint mpisr388 (725, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr388, 0.963777));

    MassPoint mpisr389 (725, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr389, 0.964327));

    MassPoint mpisr390 (725, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr390, 0.964843));

    MassPoint mpisr391 (725, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr391, 0.963589));

    MassPoint mpisr392 (725, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr392, 0.963764));

    MassPoint mpisr393 (725, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr393, 0.963805));

    MassPoint mpisr394 (725, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr394, 0.965569));

    MassPoint mpisr395 (725, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr395, 0.964297));

    MassPoint mpisr396 (725, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr396, 0.96529));

    MassPoint mpisr397 (750, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr397, 0.964335));

    MassPoint mpisr398 (750, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr398, 0.964247));

    MassPoint mpisr399 (750, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr399, 0.965418));

    MassPoint mpisr400 (750, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr400, 0.964599));

    MassPoint mpisr401 (750, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr401, 0.964218));

    MassPoint mpisr402 (750, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr402, 0.964468));

    MassPoint mpisr403 (750, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr403, 0.965339));

    MassPoint mpisr404 (750, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr404, 0.964441));

    MassPoint mpisr405 (750, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr405, 0.965635));

    MassPoint mpisr406 (750, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr406, 0.96544));

    MassPoint mpisr407 (750, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr407, 0.965042));

    MassPoint mpisr408 (750, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr408, 0.965064));

    MassPoint mpisr409 (750, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr409, 0.964116));

    MassPoint mpisr410 (750, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr410, 0.965208));

    MassPoint mpisr411 (750, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr411, 0.966104));

    MassPoint mpisr412 (750, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr412, 0.964517));

    MassPoint mpisr413 (750, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr413, 0.964727));

    MassPoint mpisr414 (750, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr414, 0.964686));

    MassPoint mpisr415 (750, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr415, 0.964907));

    MassPoint mpisr416 (750, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr416, 0.964319));

    MassPoint mpisr417 (750, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr417, 0.963546));

    MassPoint mpisr418 (750, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr418, 0.964431));

    MassPoint mpisr419 (750, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr419, 0.964545));

    MassPoint mpisr420 (750, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr420, 0.96457));

    MassPoint mpisr421 (750, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr421, 0.964626));

    MassPoint mpisr422 (775, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr422, 0.964514));

    MassPoint mpisr423 (775, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr423, 0.964508));

    MassPoint mpisr424 (775, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr424, 0.963382));

    MassPoint mpisr425 (775, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr425, 0.965897));

    MassPoint mpisr426 (775, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr426, 0.964148));

    MassPoint mpisr427 (775, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr427, 0.964421));

    MassPoint mpisr428 (775, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr428, 0.964335));

    MassPoint mpisr429 (775, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr429, 0.964518));

    MassPoint mpisr430 (775, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr430, 0.965529));

    MassPoint mpisr431 (775, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr431, 0.963445));

    MassPoint mpisr432 (775, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr432, 0.962777));

    MassPoint mpisr433 (775, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr433, 0.964451));

    MassPoint mpisr434 (775, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr434, 0.963749));

    MassPoint mpisr435 (775, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr435, 0.965206));

    MassPoint mpisr436 (775, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr436, 0.965412));

    MassPoint mpisr437 (775, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr437, 0.963266));

    MassPoint mpisr438 (775, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr438, 0.964387));

    MassPoint mpisr439 (775, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr439, 0.965813));

    MassPoint mpisr440 (775, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr440, 0.965764));

    MassPoint mpisr441 (775, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr441, 0.963184));

    MassPoint mpisr442 (775, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr442, 0.964012));

    MassPoint mpisr443 (775, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr443, 0.964581));

    MassPoint mpisr444 (775, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr444, 0.965146));

    MassPoint mpisr445 (775, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr445, 0.965241));

    MassPoint mpisr446 (775, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr446, 0.963854));

    MassPoint mpisr447 (800, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr447, 0.966485));

    MassPoint mpisr448 (800, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr448, 0.96528));

    MassPoint mpisr449 (800, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr449, 0.963259));

    MassPoint mpisr450 (800, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr450, 0.964174));

    MassPoint mpisr451 (800, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr451, 0.964285));

    MassPoint mpisr452 (800, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr452, 0.963058));

    MassPoint mpisr453 (800, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr453, 0.962901));

    MassPoint mpisr454 (800, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr454, 0.963848));

    MassPoint mpisr455 (800, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr455, 0.965455));

    MassPoint mpisr456 (800, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr456, 0.964928));

    MassPoint mpisr457 (800, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr457, 0.964659));

    MassPoint mpisr458 (800, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr458, 0.965219));

    MassPoint mpisr459 (800, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr459, 0.963307));

    MassPoint mpisr460 (800, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr460, 0.965143));

    MassPoint mpisr461 (800, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr461, 0.965241));

    MassPoint mpisr462 (800, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr462, 0.962841));

    MassPoint mpisr463 (800, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr463, 0.963737));

    MassPoint mpisr464 (800, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr464, 0.965522));

    MassPoint mpisr465 (800, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr465, 0.964903));

    MassPoint mpisr466 (800, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr466, 0.963894));

    MassPoint mpisr467 (800, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr467, 0.963667));

    MassPoint mpisr468 (800, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr468, 0.965456));

    MassPoint mpisr469 (800, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr469, 0.964662));

    MassPoint mpisr470 (800, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr470, 0.964571));

    MassPoint mpisr471 (800, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpisr471, 0.96363)); 

    MassPoint mpextisr0 (825, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr0, 0.964168));

    MassPoint mpextisr1 (825, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr1, 0.964375));

    MassPoint mpextisr2 (825, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr2, 0.964693));

    MassPoint mpextisr3 (825, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr3, 0.963942));

    MassPoint mpextisr4 (825, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr4, 0.964529));

    MassPoint mpextisr5 (825, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr5, 0.96382));

    MassPoint mpextisr6 (825, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr6, 0.966221));

    MassPoint mpextisr7 (825, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr7, 0.965196));

    MassPoint mpextisr8 (825, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr8, 0.964467));

    MassPoint mpextisr9 (825, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr9, 0.965099));

    MassPoint mpextisr10 (825, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr10, 0.963552));

    MassPoint mpextisr11 (825, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr11, 0.965578));

    MassPoint mpextisr12 (825, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr12, 0.964495));

    MassPoint mpextisr13 (825, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr13, 0.964287));

    MassPoint mpextisr14 (825, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr14, 0.964752));

    MassPoint mpextisr15 (825, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr15, 0.964429));

    MassPoint mpextisr16 (825, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr16, 0.964602));

    MassPoint mpextisr17 (825, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr17, 0.965935));

    MassPoint mpextisr18 (825, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr18, 0.963812));

    MassPoint mpextisr19 (825, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr19, 0.963736));

    MassPoint mpextisr20 (825, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr20, 0.964048));

    MassPoint mpextisr21 (825, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr21, 0.963788));

    MassPoint mpextisr22 (825, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr22, 0.965663));

    MassPoint mpextisr23 (825, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr23, 0.964258));

    MassPoint mpextisr24 (825, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr24, 0.964277));

    MassPoint mpextisr25 (850, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr25, 0.965386));

    MassPoint mpextisr26 (850, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr26, 0.964786));

    MassPoint mpextisr27 (850, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr27, 0.96509));

    MassPoint mpextisr28 (850, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr28, 0.96487));

    MassPoint mpextisr29 (850, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr29, 0.965091));

    MassPoint mpextisr30 (850, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr30, 0.965118));

    MassPoint mpextisr31 (850, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr31, 0.96461));

    MassPoint mpextisr32 (850, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr32, 0.96533));

    MassPoint mpextisr33 (850, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr33, 0.965606));

    MassPoint mpextisr34 (850, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr34, 0.964426));

    MassPoint mpextisr35 (850, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr35, 0.964979));

    MassPoint mpextisr36 (850, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr36, 0.965498));

    MassPoint mpextisr37 (850, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr37, 0.963403));

    MassPoint mpextisr38 (850, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr38, 0.964459));

    MassPoint mpextisr39 (850, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr39, 0.964266));

    MassPoint mpextisr40 (850, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr40, 0.964628));

    MassPoint mpextisr41 (850, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr41, 0.96503));

    MassPoint mpextisr42 (850, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr42, 0.96481));

    MassPoint mpextisr43 (850, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr43, 0.964638));

    MassPoint mpextisr44 (850, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr44, 0.965222));

    MassPoint mpextisr45 (850, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr45, 0.963785));

    MassPoint mpextisr46 (850, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr46, 0.965439));

    MassPoint mpextisr47 (850, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr47, 0.965027));

    MassPoint mpextisr48 (850, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr48, 0.965048));

    MassPoint mpextisr49 (850, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr49, 0.964748));

    MassPoint mpextisr50 (875, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr50, 0.964753));

    MassPoint mpextisr51 (875, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr51, 0.963698));

    MassPoint mpextisr52 (875, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr52, 0.964047));

    MassPoint mpextisr53 (875, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr53, 0.964696));

    MassPoint mpextisr54 (875, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr54, 0.963867));

    MassPoint mpextisr55 (875, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr55, 0.964609));

    MassPoint mpextisr56 (875, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr56, 0.963602));

    MassPoint mpextisr57 (875, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr57, 0.964591));

    MassPoint mpextisr58 (875, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr58, 0.963809));

    MassPoint mpextisr59 (875, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr59, 0.963869));

    MassPoint mpextisr60 (875, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr60, 0.964988));

    MassPoint mpextisr61 (875, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr61, 0.964485));

    MassPoint mpextisr62 (875, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr62, 0.963995));

    MassPoint mpextisr63 (875, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr63, 0.964228));

    MassPoint mpextisr64 (875, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr64, 0.964354));

    MassPoint mpextisr65 (875, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr65, 0.964835));

    MassPoint mpextisr66 (875, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr66, 0.964206));

    MassPoint mpextisr67 (875, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr67, 0.9648));

    MassPoint mpextisr68 (875, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr68, 0.963796));

    MassPoint mpextisr69 (875, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr69, 0.963288));

    MassPoint mpextisr70 (875, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr70, 0.966201));

    MassPoint mpextisr71 (875, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr71, 0.964499));

    MassPoint mpextisr72 (875, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr72, 0.965403));

    MassPoint mpextisr73 (875, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr73, 0.964439));

    MassPoint mpextisr74 (875, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr74, 0.964408));

    MassPoint mpextisr75 (900, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr75, 0.964264));

    MassPoint mpextisr76 (900, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr76, 0.96345));

    MassPoint mpextisr77 (900, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr77, 0.962879));

    MassPoint mpextisr78 (900, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr78, 0.964744));

    MassPoint mpextisr79 (900, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr79, 0.965806));

    MassPoint mpextisr80 (900, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr80, 0.964038));

    MassPoint mpextisr81 (900, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr81, 0.963684));

    MassPoint mpextisr82 (900, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr82, 0.963218));

    MassPoint mpextisr83 (900, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr83, 0.964231));

    MassPoint mpextisr84 (900, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr84, 0.964518));

    MassPoint mpextisr85 (900, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr85, 0.96486));

    MassPoint mpextisr86 (900, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr86, 0.96543));

    MassPoint mpextisr87 (900, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr87, 0.963393));

    MassPoint mpextisr88 (900, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr88, 0.964422));

    MassPoint mpextisr89 (900, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr89, 0.965396));

    MassPoint mpextisr90 (900, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr90, 0.965121));

    MassPoint mpextisr91 (900, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr91, 0.964839));

    MassPoint mpextisr92 (900, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr92, 0.9648));

    MassPoint mpextisr93 (900, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr93, 0.964212));

    MassPoint mpextisr94 (900, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr94, 0.964037));

    MassPoint mpextisr95 (900, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr95, 0.963044));

    MassPoint mpextisr96 (900, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr96, 0.964376));

    MassPoint mpextisr97 (900, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr97, 0.963744));

    MassPoint mpextisr98 (900, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr98, 0.96496));

    MassPoint mpextisr99 (900, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr99, 0.964466));

    MassPoint mpextisr100 (925, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr100, 0.963484));

    MassPoint mpextisr101 (925, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr101, 0.965017));

    MassPoint mpextisr102 (925, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr102, 0.964795));

    MassPoint mpextisr103 (925, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr103, 0.964608));

    MassPoint mpextisr104 (925, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr104, 0.964735));

    MassPoint mpextisr105 (925, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr105, 0.964382));

    MassPoint mpextisr106 (925, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr106, 0.963825));

    MassPoint mpextisr107 (925, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr107, 0.965663));

    MassPoint mpextisr108 (925, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr108, 0.964633));

    MassPoint mpextisr109 (925, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr109, 0.965421));

    MassPoint mpextisr110 (925, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr110, 0.964619));

    MassPoint mpextisr111 (925, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr111, 0.965748));

    MassPoint mpextisr112 (925, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr112, 0.965657));

    MassPoint mpextisr113 (925, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr113, 0.964675));

    MassPoint mpextisr114 (925, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr114, 0.964532));

    MassPoint mpextisr115 (925, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr115, 0.965251));

    MassPoint mpextisr116 (925, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr116, 0.96413));

    MassPoint mpextisr117 (925, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr117, 0.96454));

    MassPoint mpextisr118 (925, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr118, 0.963473));

    MassPoint mpextisr119 (925, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr119, 0.964165));

    MassPoint mpextisr120 (925, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr120, 0.964769));

    MassPoint mpextisr121 (925, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr121, 0.966214));

    MassPoint mpextisr122 (925, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr122, 0.96424));

    MassPoint mpextisr123 (925, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr123, 0.964019));

    MassPoint mpextisr124 (925, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr124, 0.965463));

    MassPoint mpextisr125 (950, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr125, 0.963966));

    MassPoint mpextisr126 (950, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr126, 0.964226));

    MassPoint mpextisr127 (950, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr127, 0.964163));

    MassPoint mpextisr128 (950, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr128, 0.96472));

    MassPoint mpextisr129 (950, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr129, 0.965079));

    MassPoint mpextisr130 (950, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr130, 0.965125));

    MassPoint mpextisr131 (950, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr131, 0.96457));

    MassPoint mpextisr132 (950, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr132, 0.963848));

    MassPoint mpextisr133 (950, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr133, 0.964464));

    MassPoint mpextisr134 (950, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr134, 0.964944));

    MassPoint mpextisr135 (950, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr135, 0.9635));

    MassPoint mpextisr136 (950, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr136, 0.965847));

    MassPoint mpextisr137 (950, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr137, 0.965024));

    MassPoint mpextisr138 (950, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr138, 0.964954));

    MassPoint mpextisr139 (950, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr139, 0.963754));

    MassPoint mpextisr140 (950, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr140, 0.963899));

    MassPoint mpextisr141 (950, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr141, 0.965074));

    MassPoint mpextisr142 (950, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr142, 0.964824));

    MassPoint mpextisr143 (950, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr143, 0.964296));

    MassPoint mpextisr144 (950, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr144, 0.964113));

    MassPoint mpextisr145 (950, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr145, 0.963537));

    MassPoint mpextisr146 (950, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr146, 0.964353));

    MassPoint mpextisr147 (950, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr147, 0.963207));

    MassPoint mpextisr148 (950, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr148, 0.964251));

    MassPoint mpextisr149 (950, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr149, 0.964214));

    MassPoint mpextisr150 (975, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr150, 0.965225));

    MassPoint mpextisr151 (975, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr151, 0.964572));

    MassPoint mpextisr152 (975, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr152, 0.963685));

    MassPoint mpextisr153 (975, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr153, 0.964508));

    MassPoint mpextisr154 (975, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr154, 0.964693));

    MassPoint mpextisr155 (975, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr155, 0.964541));

    MassPoint mpextisr156 (975, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr156, 0.965355));

    MassPoint mpextisr157 (975, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr157, 0.963401));

    MassPoint mpextisr158 (975, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr158, 0.964996));

    MassPoint mpextisr159 (975, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr159, 0.964159));

    MassPoint mpextisr160 (975, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr160, 0.965023));

    MassPoint mpextisr161 (975, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr161, 0.965313));

    MassPoint mpextisr162 (975, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr162, 0.964879));

    MassPoint mpextisr163 (975, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr163, 0.963688));

    MassPoint mpextisr164 (975, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr164, 0.964481));

    MassPoint mpextisr165 (975, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr165, 0.965136));

    MassPoint mpextisr166 (975, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr166, 0.96436));

    MassPoint mpextisr167 (975, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr167, 0.964517));

    MassPoint mpextisr168 (975, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr168, 0.965676));

    MassPoint mpextisr169 (975, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr169, 0.963897));

    MassPoint mpextisr170 (975, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr170, 0.963542));

    MassPoint mpextisr171 (975, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr171, 0.963664));

    MassPoint mpextisr172 (975, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr172, 0.963173));

    MassPoint mpextisr173 (975, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr173, 0.964327));

    MassPoint mpextisr174 (975, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr174, 0.964977));

    MassPoint mpextisr175 (1000, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr175, 0.965279));

    MassPoint mpextisr176 (1000, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr176, 0.964222));

    MassPoint mpextisr177 (1000, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr177, 0.963601));

    MassPoint mpextisr178 (1000, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr178, 0.96407));

    MassPoint mpextisr179 (1000, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr179, 0.965639));

    MassPoint mpextisr180 (1000, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr180, 0.9656));

    MassPoint mpextisr181 (1000, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr181, 0.964916));

    MassPoint mpextisr182 (1000, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr182, 0.963958));

    MassPoint mpextisr183 (1000, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr183, 0.96453));

    MassPoint mpextisr184 (1000, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr184, 0.963703));

    MassPoint mpextisr185 (1000, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr185, 0.964405));

    MassPoint mpextisr186 (1000, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr186, 0.964236));

    MassPoint mpextisr187 (1000, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr187, 0.965098));

    MassPoint mpextisr188 (1000, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr188, 0.965008));

    MassPoint mpextisr189 (1000, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr189, 0.963953));

    MassPoint mpextisr190 (1000, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr190, 0.963697));

    MassPoint mpextisr191 (1000, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr191, 0.965479));

    MassPoint mpextisr192 (1000, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr192, 0.963312));

    MassPoint mpextisr193 (1000, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr193, 0.964043));

    MassPoint mpextisr194 (1000, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr194, 0.96481));

    MassPoint mpextisr195 (1000, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr195, 0.964194));

    MassPoint mpextisr196 (1000, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr196, 0.963148));

    MassPoint mpextisr197 (1000, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr197, 0.962871));

    MassPoint mpextisr198 (1000, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr198, 0.963774));

    MassPoint mpextisr199 (1000, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr199, 0.963659));

    MassPoint mpextisr200 (1025, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr200, 0.96427));

    MassPoint mpextisr201 (1025, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr201, 0.963859));

    MassPoint mpextisr202 (1025, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr202, 0.963663));

    MassPoint mpextisr203 (1025, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr203, 0.964333));

    MassPoint mpextisr204 (1025, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr204, 0.96437));

    MassPoint mpextisr205 (1025, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr205, 0.965647));

    MassPoint mpextisr206 (1025, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr206, 0.964651));

    MassPoint mpextisr207 (1025, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr207, 0.963193));

    MassPoint mpextisr208 (1025, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr208, 0.964108));

    MassPoint mpextisr209 (1025, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr209, 0.963751));

    MassPoint mpextisr210 (1025, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr210, 0.96388));

    MassPoint mpextisr211 (1025, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr211, 0.963601));

    MassPoint mpextisr212 (1025, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr212, 0.962937));

    MassPoint mpextisr213 (1025, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr213, 0.96446));

    MassPoint mpextisr214 (1025, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr214, 0.96485));

    MassPoint mpextisr215 (1025, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr215, 0.964568));

    MassPoint mpextisr216 (1025, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr216, 0.964644));

    MassPoint mpextisr217 (1025, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr217, 0.963905));

    MassPoint mpextisr218 (1025, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr218, 0.964044));

    MassPoint mpextisr219 (1025, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr219, 0.966088));

    MassPoint mpextisr220 (1025, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr220, 0.964519));

    MassPoint mpextisr221 (1025, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr221, 0.964048));

    MassPoint mpextisr222 (1025, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr222, 0.964554));

    MassPoint mpextisr223 (1025, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr223, 0.96488));

    MassPoint mpextisr224 (1025, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr224, 0.96368));

    MassPoint mpextisr225 (1050, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr225, 0.964514));

    MassPoint mpextisr226 (1050, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr226, 0.965266));

    MassPoint mpextisr227 (1050, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr227, 0.964959));

    MassPoint mpextisr228 (1050, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr228, 0.963248));

    MassPoint mpextisr229 (1050, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr229, 0.965007));

    MassPoint mpextisr230 (1050, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr230, 0.964789));

    MassPoint mpextisr231 (1050, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr231, 0.964296));

    MassPoint mpextisr232 (1050, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr232, 0.964469));

    MassPoint mpextisr233 (1050, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr233, 0.96427));

    MassPoint mpextisr234 (1050, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr234, 0.965736));

    MassPoint mpextisr235 (1050, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr235, 0.964649));

    MassPoint mpextisr236 (1050, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr236, 0.964358));

    MassPoint mpextisr237 (1050, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr237, 0.965295));

    MassPoint mpextisr238 (1050, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr238, 0.964147));

    MassPoint mpextisr239 (1050, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr239, 0.965328));

    MassPoint mpextisr240 (1050, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr240, 0.965256));

    MassPoint mpextisr241 (1050, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr241, 0.964209));

    MassPoint mpextisr242 (1050, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr242, 0.964832));

    MassPoint mpextisr243 (1050, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr243, 0.965513));

    MassPoint mpextisr244 (1050, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr244, 0.964415));

    MassPoint mpextisr245 (1050, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr245, 0.963679));

    MassPoint mpextisr246 (1050, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr246, 0.964788));

    MassPoint mpextisr247 (1050, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr247, 0.964387));

    MassPoint mpextisr248 (1050, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr248, 0.963018));

    MassPoint mpextisr249 (1050, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr249, 0.964544));

    MassPoint mpextisr250 (1075, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr250, 0.964737));

    MassPoint mpextisr251 (1075, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr251, 0.964239));

    MassPoint mpextisr252 (1075, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr252, 0.964085));

    MassPoint mpextisr253 (1075, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr253, 0.96429));

    MassPoint mpextisr254 (1075, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr254, 0.964733));

    MassPoint mpextisr255 (1075, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr255, 0.964248));

    MassPoint mpextisr256 (1075, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr256, 0.964097));

    MassPoint mpextisr257 (1075, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr257, 0.961808));

    MassPoint mpextisr258 (1075, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr258, 0.963918));

    MassPoint mpextisr259 (1075, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr259, 0.964635));

    MassPoint mpextisr260 (1075, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr260, 0.965103));

    MassPoint mpextisr261 (1075, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr261, 0.965696));

    MassPoint mpextisr262 (1075, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr262, 0.965386));

    MassPoint mpextisr263 (1075, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr263, 0.962371));

    MassPoint mpextisr264 (1075, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr264, 0.964616));

    MassPoint mpextisr265 (1075, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr265, 0.966295));

    MassPoint mpextisr266 (1075, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr266, 0.964244));

    MassPoint mpextisr267 (1075, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr267, 0.964093));

    MassPoint mpextisr268 (1075, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr268, 0.965159));

    MassPoint mpextisr269 (1075, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr269, 0.963935));

    MassPoint mpextisr270 (1075, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr270, 0.963582));

    MassPoint mpextisr271 (1075, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr271, 0.964106));

    MassPoint mpextisr272 (1075, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr272, 0.965512));

    MassPoint mpextisr273 (1075, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr273, 0.964106));

    MassPoint mpextisr274 (1075, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr274, 0.963995));

    MassPoint mpextisr275 (1100, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr275, 0.964632));

    MassPoint mpextisr276 (1100, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr276, 0.965384));

    MassPoint mpextisr277 (1100, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr277, 0.964604));

    MassPoint mpextisr278 (1100, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr278, 0.964683));

    MassPoint mpextisr279 (1100, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr279, 0.964693));

    MassPoint mpextisr280 (1100, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr280, 0.964429));

    MassPoint mpextisr281 (1100, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr281, 0.965072));

    MassPoint mpextisr282 (1100, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr282, 0.965259));

    MassPoint mpextisr283 (1100, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr283, 0.964215));

    MassPoint mpextisr284 (1100, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr284, 0.964857));

    MassPoint mpextisr285 (1100, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr285, 0.964181));

    MassPoint mpextisr286 (1100, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr286, 0.96524));

    MassPoint mpextisr287 (1100, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr287, 0.964121));

    MassPoint mpextisr288 (1100, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr288, 0.963973));

    MassPoint mpextisr289 (1100, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr289, 0.963619));

    MassPoint mpextisr290 (1100, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr290, 0.964593));

    MassPoint mpextisr291 (1100, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr291, 0.96542));

    MassPoint mpextisr292 (1100, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr292, 0.965503));

    MassPoint mpextisr293 (1100, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr293, 0.964285));

    MassPoint mpextisr294 (1100, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr294, 0.964792));

    MassPoint mpextisr295 (1100, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr295, 0.96374));

    MassPoint mpextisr296 (1100, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr296, 0.964634));

    MassPoint mpextisr297 (1100, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr297, 0.96405));

    MassPoint mpextisr298 (1100, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr298, 0.96453));

    MassPoint mpextisr299 (1100, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr299, 0.965175));

    MassPoint mpextisr300 (1125, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr300, 0.96395));

    MassPoint mpextisr301 (1125, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr301, 0.965144));

    MassPoint mpextisr302 (1125, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr302, 0.965006));

    MassPoint mpextisr303 (1125, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr303, 0.964331));

    MassPoint mpextisr304 (1125, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr304, 0.965099));

    MassPoint mpextisr305 (1125, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr305, 0.963981));

    MassPoint mpextisr306 (1125, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr306, 0.965062));

    MassPoint mpextisr307 (1125, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr307, 0.964552));

    MassPoint mpextisr308 (1125, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr308, 0.964622));

    MassPoint mpextisr309 (1125, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr309, 0.965512));

    MassPoint mpextisr310 (1125, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr310, 0.965363));

    MassPoint mpextisr311 (1125, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr311, 0.964386));

    MassPoint mpextisr312 (1125, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr312, 0.964565));

    MassPoint mpextisr313 (1125, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr313, 0.96415));

    MassPoint mpextisr314 (1125, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr314, 0.965767));

    MassPoint mpextisr315 (1125, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr315, 0.965632));

    MassPoint mpextisr316 (1125, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr316, 0.965235));

    MassPoint mpextisr317 (1125, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr317, 0.963272));

    MassPoint mpextisr318 (1125, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr318, 0.963881));

    MassPoint mpextisr319 (1125, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr319, 0.964627));

    MassPoint mpextisr320 (1125, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr320, 0.964677));

    MassPoint mpextisr321 (1125, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr321, 0.965414));

    MassPoint mpextisr322 (1125, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr322, 0.96427));

    MassPoint mpextisr323 (1125, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr323, 0.964424));

    MassPoint mpextisr324 (1125, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr324, 0.967451));

    MassPoint mpextisr325 (1150, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr325, 0.965409));

    MassPoint mpextisr326 (1150, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr326, 0.964054));

    MassPoint mpextisr327 (1150, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr327, 0.964983));

    MassPoint mpextisr328 (1150, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr328, 0.963982));

    MassPoint mpextisr329 (1150, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr329, 0.964854));

    MassPoint mpextisr330 (1150, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr330, 0.965132));

    MassPoint mpextisr331 (1150, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr331, 0.964967));

    MassPoint mpextisr332 (1150, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr332, 0.964456));

    MassPoint mpextisr333 (1150, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr333, 0.965696));

    MassPoint mpextisr334 (1150, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr334, 0.963807));

    MassPoint mpextisr335 (1150, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr335, 0.964777));

    MassPoint mpextisr336 (1150, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr336, 0.964993));

    MassPoint mpextisr337 (1150, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr337, 0.964896));

    MassPoint mpextisr338 (1150, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr338, 0.963944));

    MassPoint mpextisr339 (1150, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr339, 0.964558));

    MassPoint mpextisr340 (1150, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr340, 0.96474));

    MassPoint mpextisr341 (1150, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr341, 0.964808));

    MassPoint mpextisr342 (1150, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr342, 0.964078));

    MassPoint mpextisr343 (1150, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr343, 0.96475));

    MassPoint mpextisr344 (1150, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr344, 0.964286));

    MassPoint mpextisr345 (1150, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr345, 0.964531));

    MassPoint mpextisr346 (1150, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr346, 0.964286));

    MassPoint mpextisr347 (1150, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr347, 0.964641));

    MassPoint mpextisr348 (1150, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr348, 0.964541));

    MassPoint mpextisr349 (1150, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr349, 0.963941));

    MassPoint mpextisr350 (1175, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr350, 0.965145));

    MassPoint mpextisr351 (1175, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr351, 0.964782));

    MassPoint mpextisr352 (1175, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr352, 0.963424));

    MassPoint mpextisr353 (1175, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr353, 0.966404));

    MassPoint mpextisr354 (1175, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr354, 0.965211));

    MassPoint mpextisr355 (1175, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr355, 0.964945));

    MassPoint mpextisr356 (1175, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr356, 0.964921));

    MassPoint mpextisr357 (1175, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr357, 0.966636));

    MassPoint mpextisr358 (1175, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr358, 0.964998));

    MassPoint mpextisr359 (1175, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr359, 0.963744));

    MassPoint mpextisr360 (1175, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr360, 0.963793));

    MassPoint mpextisr361 (1175, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr361, 0.964429));

    MassPoint mpextisr362 (1175, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr362, 0.964271));

    MassPoint mpextisr363 (1175, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr363, 0.964494));

    MassPoint mpextisr364 (1175, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr364, 0.964572));

    MassPoint mpextisr365 (1175, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr365, 0.963985));

    MassPoint mpextisr366 (1175, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr366, 0.965101));

    MassPoint mpextisr367 (1175, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr367, 0.964557));

    MassPoint mpextisr368 (1175, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr368, 0.963826));

    MassPoint mpextisr369 (1175, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr369, 0.963818));

    MassPoint mpextisr370 (1175, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr370, 0.966207));

    MassPoint mpextisr371 (1175, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr371, 0.964882));

    MassPoint mpextisr372 (1175, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr372, 0.963717));

    MassPoint mpextisr373 (1175, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr373, 0.964517));

    MassPoint mpextisr374 (1175, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr374, 0.964668));

    MassPoint mpextisr375 (1200, 1);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr375, 0.964554));

    MassPoint mpextisr376 (1200, 25);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr376, 0.964235));

    MassPoint mpextisr377 (1200, 50);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr377, 0.964606));

    MassPoint mpextisr378 (1200, 75);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr378, 0.964498));

    MassPoint mpextisr379 (1200, 100);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr379, 0.963981));

    MassPoint mpextisr380 (1200, 125);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr380, 0.965309));

    MassPoint mpextisr381 (1200, 150);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr381, 0.964492));

    MassPoint mpextisr382 (1200, 175);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr382, 0.965112));

    MassPoint mpextisr383 (1200, 200);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr383, 0.966659));

    MassPoint mpextisr384 (1200, 225);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr384, 0.965315));

    MassPoint mpextisr385 (1200, 250);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr385, 0.964079));

    MassPoint mpextisr386 (1200, 275);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr386, 0.964282));

    MassPoint mpextisr387 (1200, 300);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr387, 0.964171));

    MassPoint mpextisr388 (1200, 325);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr388, 0.964664));

    MassPoint mpextisr389 (1200, 350);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr389, 0.96569));

    MassPoint mpextisr390 (1200, 375);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr390, 0.964475));

    MassPoint mpextisr391 (1200, 400);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr391, 0.965146));

    MassPoint mpextisr392 (1200, 425);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr392, 0.964558));

    MassPoint mpextisr393 (1200, 450);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr393, 0.964727));

    MassPoint mpextisr394 (1200, 475);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr394, 0.964368));

    MassPoint mpextisr395 (1200, 500);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr395, 0.965201));

    MassPoint mpextisr396 (1200, 525);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr396, 0.964522));

    MassPoint mpextisr397 (1200, 550);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr397, 0.964042));

    MassPoint mpextisr398 (1200, 575);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr398, 0.965902));

    MassPoint mpextisr399 (1200, 600);
    StopNeutralinoISRMap.insert(std::make_pair(mpextisr399, 0.964453));

  }

}

void AnalysisStop::CorrectEventWeight() {

  float EventBTagSF       = 1.;
  float EventBTagSF_Upb   = 1.;
  float EventBTagSF_Dob   = 1.;
  float EventBTagSF_Upl   = 1.;
  float EventBTagSF_Dol   = 1.;
  float EventBTagSF_UpFSb = 1.;
  float EventBTagSF_DoFSb = 1.;
  
  for (int ijet = 0; ijet<_njet; ijet++) {
    
    if (AnalysisJets[ijet].v.Pt() <= 20.) continue; 
    
    int ThisIndex = AnalysisJets[ijet].index;
    int ThisFlavour = std_vector_jet_HadronFlavour->at(ThisIndex);
    float MonteCarloEfficiency = BTagSF->JetTagEfficiency(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency       = MonteCarloEfficiency*BTagSF->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_Upb   = MonteCarloEfficiency*BTagSF_Upb->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_Dob   = MonteCarloEfficiency*BTagSF_Dob->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_Upl   = MonteCarloEfficiency*BTagSF_Upl->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_Dol   = MonteCarloEfficiency*BTagSF_Dol->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_UpFSb = MonteCarloEfficiency*BTagSF_UpFSb->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());
    float DataEfficiency_DoFSb = MonteCarloEfficiency*BTagSF_DoFSb->GetJetSF(ThisFlavour, AnalysisJets[ijet].v.Pt(), AnalysisJets[ijet].v.Eta());

    if (AnalysisJets[ijet].csvv2ivf>CSVv2M) {
      EventBTagSF       *= DataEfficiency/MonteCarloEfficiency;
      EventBTagSF_Upb   *= DataEfficiency_Upb/MonteCarloEfficiency;
      EventBTagSF_Dob   *= DataEfficiency_Dob/MonteCarloEfficiency;
      EventBTagSF_Upl   *= DataEfficiency_Upl/MonteCarloEfficiency;
      EventBTagSF_Dol   *= DataEfficiency_Dol/MonteCarloEfficiency;
      EventBTagSF_UpFSb *= DataEfficiency_UpFSb/MonteCarloEfficiency;
      EventBTagSF_DoFSb *= DataEfficiency_DoFSb/MonteCarloEfficiency;
    } else {
      EventBTagSF       *= (1. - DataEfficiency)/(1. - MonteCarloEfficiency);
      EventBTagSF_Upb   *= (1. - DataEfficiency_Upb)/(1. - MonteCarloEfficiency);
      EventBTagSF_Dob   *= (1. - DataEfficiency_Dob)/(1. - MonteCarloEfficiency);
      EventBTagSF_Upl   *= (1. - DataEfficiency_Upl)/(1. - MonteCarloEfficiency);
      EventBTagSF_Dol   *= (1. - DataEfficiency_Dol)/(1. - MonteCarloEfficiency);
      EventBTagSF_UpFSb *= (1. - DataEfficiency_UpFSb)/(1. - MonteCarloEfficiency);
      EventBTagSF_DoFSb *= (1. - DataEfficiency_DoFSb)/(1. - MonteCarloEfficiency);
    }

  }

  _event_weight            *= EventBTagSF/bPogSF_CSVM;
  _event_weight_Btagup      = _event_weight*EventBTagSF_Upb/EventBTagSF; 
  _event_weight_Btagdo      = _event_weight*EventBTagSF_Dob/EventBTagSF; 
  _event_weight_Btaglightup = _event_weight*EventBTagSF_Upl/EventBTagSF; 
  _event_weight_Btaglightdo = _event_weight*EventBTagSF_Dol/EventBTagSF; 
  _event_weight_BtagFSup    = _event_weight*EventBTagSF_UpFSb/EventBTagSF; 
  _event_weight_BtagFSdo    = _event_weight*EventBTagSF_DoFSb/EventBTagSF; 
  _event_weight_Idisoup    *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Idisodo    *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Idisoeleup *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Idisoeledo *= EventBTagSF/bPogSF_CSVM;  
  _event_weight_Idisomuup  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Idisomudo  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Triggerup  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Triggerdo  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Recoup     *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Recodo     *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Fastsimup  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Fastsimdo  *= EventBTagSF/bPogSF_CSVM; 
  _event_weight_Toppt      *= EventBTagSF/bPogSF_CSVM; 
  
  if (FastSimDataset!="") { 
    
    int XMass = (SUSYProductionProcess.Contains("T2")) ? susyMstop : susyMChargino;
    int NeutralinoMass = susyMLSP;
    MassPoint ThisMassPoint (XMass, NeutralinoMass);
    MassPointParameters TheseMassPointParameters = StopNeutralinoMap.at(ThisMassPoint);
    StopCrossSection ThisStopCrossSection = TheseMassPointParameters.first;
    int SampleSize = TheseMassPointParameters.second;
    _event_weight            *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW;
    _event_weight_Btagup     *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Btagdo     *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW;
    _event_weight_Btaglightup*= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Btaglightdo*= (1000.*ThisStopCrossSection.first/SampleSize)/baseW;
    _event_weight_BtagFSup   *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_BtagFSdo   *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisoup    *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisodo    *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisoeleup *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisoeledo *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisomuup  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Idisomudo  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Triggerup  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Triggerdo  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Recoup     *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Recodo     *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Fastsimup  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW; 
    _event_weight_Fastsimdo  *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW;
    _event_weight_Toppt      *= (1000.*ThisStopCrossSection.first/SampleSize)/baseW;
    
  }
  
}

bool AnalysisStop::PassFastsimJetsCleanup() {

  // Cleaning up of fastsim jets 
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17#Cleaning_up_of_fastsim_jets_from

  if (!_isfastsim) return true;
  
  for (int ijet = 0; ijet<_njet; ijet++) 
    if (AnalysisJets[ijet].v.Pt()>20 && fabs(AnalysisJets[ijet].v.Eta())<2.5) 
      if (std_vector_jet_chargedHadronFraction->at(AnalysisJets[ijet].index)<0.1) {

	bool isgenmatched = false;
	
	for (int gjet = 0; gjet<std_vector_jetGen_pt->size(); gjet++) {
	  
	  if (std_vector_jetGen_pt->at(gjet)<0.) continue;

	  TLorentzVector ThisGenJet; 
	  ThisGenJet.SetPtEtaPhiM(std_vector_jetGen_pt->at(gjet), 
				  std_vector_jetGen_eta->at(gjet), 
				  std_vector_jetGen_phi->at(gjet), 0.0);

	  if (fabs(AnalysisJets[ijet].v.DeltaR(ThisGenJet))<0.3) isgenmatched = true;

	}
	
	if (!isgenmatched) return false;
      
      }

  return true;
  
}

void AnalysisStop::ApplyISRReweighting() {

  int XMass; float _isr_weight = 1.;

  if (SUSYProductionProcess.Contains("T2")) {

    XMass = susyMstop;

    int _isrbin = (_nisrjet<nISRMultiplicityBins) ? _nisrjet : (nISRMultiplicityBins-1);
    _isr_weight = ISRBinWeight[_isrbin];

  } else {

    XMass = susyMChargino;

    int _isrbin = -1;
    for (int ib = 0; ib<nISRPtBins; ib++) 
      if (_isrpt>=ISRPtBinEdge[ib]) _isrbin = ib;
    _isr_weight = ISRPtBinWeight[_isrbin];

  }

  int NeutralinoMass = susyMLSP;
  MassPoint ThisMassPoint (XMass, NeutralinoMass);
  if (SUSYProductionProcess=="T2tt" && XMass <=250.) 
    _isr_weight *= 1.1;
  else
    _isr_weight *= StopNeutralinoISRMap.at(ThisMassPoint);

  float _isr_weight_up = (1. + _isr_weight)/2.;
  float _isr_weight_do = _isr_weight - (_isr_weight_up - _isr_weight);

  _event_weight            *= _isr_weight;
  _event_weight_Btagup     *= _isr_weight; 
  _event_weight_Btagdo     *= _isr_weight;
  _event_weight_Btaglightup*= _isr_weight; 
  _event_weight_Btaglightdo*= _isr_weight;
  _event_weight_BtagFSup   *= _isr_weight; 
  _event_weight_BtagFSdo   *= _isr_weight; 
  _event_weight_Idisoup    *= _isr_weight; 
  _event_weight_Idisodo    *= _isr_weight;  
  _event_weight_Idisoeleup *= _isr_weight; 
  _event_weight_Idisoeledo *= _isr_weight;  
  _event_weight_Idisomuup  *= _isr_weight; 
  _event_weight_Idisomudo  *= _isr_weight; 
  _event_weight_Triggerup  *= _isr_weight; 
  _event_weight_Triggerdo  *= _isr_weight; 
  _event_weight_Recoup     *= _isr_weight; 
  _event_weight_Recodo     *= _isr_weight; 
  _event_weight_Fastsimup  *= _isr_weight; 
  _event_weight_Fastsimdo  *= _isr_weight;
  _event_weight_Toppt      *= _isr_weight;
  _event_weight_Isrnjetup   = _event_weight*_isr_weight_up/_isr_weight;
  _event_weight_Isrnjetdo   = _event_weight*_isr_weight_do/_isr_weight;

}

void AnalysisStop::ApplyStopTriggerEfficiency() {

  float pt1 = (_lep1pt<1000) ? _lep1pt : 999.;
  float pt2 = (_lep2pt<1000) ? _lep2pt : 999.;

  float StopTriggerEfficiency = -1.;

  if (fabs(_lep1id)!=fabs(_lep2id)) {
    Int_t trgbin = TrgEff_em->FindFixBin(pt1, pt2);
    StopTriggerEfficiency = TrgEff_em->GetEfficiency(trgbin);
  } else if (fabs(_lep1id)==11) {
    Int_t trgbin = TrgEff_ee->FindFixBin(pt1, pt2);
    StopTriggerEfficiency = TrgEff_ee->GetEfficiency(trgbin);
  } else if (fabs(_lep1id)==13) {
    Int_t trgbin = TrgEff_mm->FindFixBin(pt1, pt2);
    StopTriggerEfficiency = TrgEff_mm->GetEfficiency(trgbin);
  }

  if (StopTriggerEfficiency<0.) return;

  StopTriggerEfficiency /= effTrigW;

  _event_weight *= StopTriggerEfficiency;
  _event_weight_Btagup     *= StopTriggerEfficiency; 
  _event_weight_Btagdo     *= StopTriggerEfficiency;
  _event_weight_Btaglightup*= StopTriggerEfficiency; 
  _event_weight_Btaglightdo*= StopTriggerEfficiency;
  _event_weight_BtagFSup   *= StopTriggerEfficiency; 
  _event_weight_BtagFSdo   *= StopTriggerEfficiency; 
  _event_weight_Idisoup    *= StopTriggerEfficiency; 
  _event_weight_Idisodo    *= StopTriggerEfficiency;  
  _event_weight_Idisoeleup *= StopTriggerEfficiency; 
  _event_weight_Idisoeledo *= StopTriggerEfficiency;  
  _event_weight_Idisomuup  *= StopTriggerEfficiency; 
  _event_weight_Idisomudo  *= StopTriggerEfficiency; 
  _event_weight_Triggerup  *= StopTriggerEfficiency; 
  _event_weight_Triggerdo  *= StopTriggerEfficiency; 
  _event_weight_Recoup     *= StopTriggerEfficiency; 
  _event_weight_Recodo     *= StopTriggerEfficiency; 
  _event_weight_Fastsimup  *= StopTriggerEfficiency; 
  _event_weight_Fastsimdo  *= StopTriggerEfficiency;
  _event_weight_Toppt      *= StopTriggerEfficiency;
  
}

void AnalysisStop::GetMiniTree(TFile *MiniTreeFile, TString systematic) {

  fChain = (TTree*) MiniTreeFile->Get("latino");

  fChain->SetBranchAddress("event",           &event);
  fChain->SetBranchAddress("run",             &run);
  fChain->SetBranchAddress("dyll",            &_dyll);
  fChain->SetBranchAddress("ptbll",           &_ptbll);
  fChain->SetBranchAddress("m2l",             &_m2l);
  fChain->SetBranchAddress("mllbb",           &_mllbb);
  fChain->SetBranchAddress("meff",            &_meff);
  fChain->SetBranchAddress("nvtx",            &nvtx);
  fChain->SetBranchAddress("_gen_ptll",       &_gen_ptll);

  fChain->SetBranchAddress("iRunPeriod",      &iRunPeriod);
  fChain->SetBranchAddress("trpu",            &trpu);
  fChain->SetBranchAddress("puW",             &puW);
  
  fChain->SetBranchAddress("dphill",          &dphill);
  
  fChain->SetBranchAddress("metPfType1",      &metPfType1);
  fChain->SetBranchAddress("metPfType1Phi",   &metPfType1Phi);
  fChain->SetBranchAddress("metGenpt",        &metGenpt);
  fChain->SetBranchAddress("dphillmet",       &_dphillmet);
  fChain->SetBranchAddress("dphilmet1",       &dphilmet1);
  fChain->SetBranchAddress("dphimetjet",      &_dphimetjet);
  fChain->SetBranchAddress("dphijet1met",     &_dphijet1met);
  fChain->SetBranchAddress("dphimetptbll",    &_dphimetptbll); 

  fChain->SetBranchAddress("jet1pt",          &jetpt1);
  fChain->SetBranchAddress("jet1phi",         &jetphi1);
  fChain->SetBranchAddress("jet1eta",         &jeteta1);
  fChain->SetBranchAddress("jet2pt",          &jetpt2);
  fChain->SetBranchAddress("jet2phi",         &jetphi2);
  fChain->SetBranchAddress("jet2eta",         &jeteta2);
  
  fChain->SetBranchAddress("lep1eta",         &_lep1eta); 
  fChain->SetBranchAddress("lep1phi",         &_lep1phi);  
  fChain->SetBranchAddress("lep1pt",          &_lep1pt);   
  fChain->SetBranchAddress("lep1id",          &_lep1id);   
  fChain->SetBranchAddress("lep1mid",         &_lep1mid); 
  fChain->SetBranchAddress("lep1isfake",      &_lep1isfake);
  fChain->SetBranchAddress("lep1mass",        &_lep1mass);
  fChain->SetBranchAddress("lep2eta",         &_lep2eta);
  fChain->SetBranchAddress("lep2phi",         &_lep2phi);
  fChain->SetBranchAddress("lep2pt",          &_lep2pt);
  fChain->SetBranchAddress("lep2isfake",      &_lep2isfake);  
  fChain->SetBranchAddress("lep2id",          &_lep2id);  
  fChain->SetBranchAddress("lep2mid",         &_lep2mid); 
  fChain->SetBranchAddress("lep2mass",        &_lep2mass);

  //if (_systematic.Contains("multilepton")) {
    fChain->SetBranchAddress("lep3eta",         &_lep3eta); 
    fChain->SetBranchAddress("lep3phi",         &_lep3phi);  
    fChain->SetBranchAddress("lep3pt",          &_lep3pt); 
    fChain->SetBranchAddress("lep3mass",        &_lep3mass); 
    fChain->SetBranchAddress("lep3type",        &_lep3type);
    fChain->SetBranchAddress("lep3id",          &_lep3id);
    fChain->SetBranchAddress("lep3idisoW",      &_lep3idisoW);  
    fChain->SetBranchAddress("lep3mid",         &_lep3mid); 
    fChain->SetBranchAddress("lep4eta",         &_lep4eta);
    fChain->SetBranchAddress("lep4phi",         &_lep4phi);
    fChain->SetBranchAddress("lep4pt",          &_lep4pt);
    fChain->SetBranchAddress("lep4mass",        &_lep4mass); 
    fChain->SetBranchAddress("lep4type",        &_lep4type);
    fChain->SetBranchAddress("lep4id",          &_lep4id);
    fChain->SetBranchAddress("lep4idisoW",      &_lep4idisoW);  
    fChain->SetBranchAddress("lep4mid",         &_lep4mid); 
    //}

  fChain->SetBranchAddress("ht",              &_ht);
  fChain->SetBranchAddress("htjets",          &_htjets);
  fChain->SetBranchAddress("htgen",           &_htgen);
  fChain->SetBranchAddress("htnojets",        &_htnojets);
  fChain->SetBranchAddress("htvisible",       &_htvisible);
  fChain->SetBranchAddress("nisrjet",         &_nisrjet);
  fChain->SetBranchAddress("isrpt",           &_isrpt);
  fChain->SetBranchAddress("njet",            &_njet);
  fChain->SetBranchAddress("nlatinojet",      &njet);
  fChain->SetBranchAddress("nlepton",         &_nlepton);
  fChain->SetBranchAddress("ntightlepton",    &_ntightlepton);
  
  fChain->SetBranchAddress("susyMLSP",        &susyMLSP);
  fChain->SetBranchAddress("susyMstop",       &susyMstop);
  fChain->SetBranchAddress("susyMChargino",   &susyMChargino);

  fChain->SetBranchAddress("effTrigW",        &effTrigW);

  if (systematic=="nominal") 
    fChain->SetBranchAddress("eventW",           &_event_weight);
  else if (systematic=="Btagup")
    fChain->SetBranchAddress("eventW_Btagup",    &_event_weight);
  else if (systematic=="Btagdo")
    fChain->SetBranchAddress("eventW_Btagdo",    &_event_weight);
  else if (systematic=="Btaglightup")
    fChain->SetBranchAddress("eventW_Btaglightup",    &_event_weight);
  else if (systematic=="Btaglightdo")
    fChain->SetBranchAddress("eventW_Btaglightdo",    &_event_weight);
  else if (systematic=="BtagFSup")
    fChain->SetBranchAddress("eventW_BtagFSup",  &_event_weight);
  else if (systematic=="BtagFSdo")
    fChain->SetBranchAddress("eventW_BtagFSdo",  &_event_weight);
  else if (systematic=="Idisoup")
    fChain->SetBranchAddress("eventW_Idisoup",   &_event_weight);
  else if (systematic=="Idisodo") 
    fChain->SetBranchAddress("eventW_Idisodo",   &_event_weight);
  else if (systematic=="Triggerup")
    fChain->SetBranchAddress("eventW_Triggerup", &_event_weight);
  else if (systematic=="Triggerdo")
    fChain->SetBranchAddress("eventW_Triggerdo", &_event_weight);
  else if (systematic=="Recoup")
    fChain->SetBranchAddress("eventW_Recoup",    &_event_weight);
  else if (systematic=="Recodo")
    fChain->SetBranchAddress("eventW_Recodo",    &_event_weight);
  else if (systematic=="Fastsimup")
    fChain->SetBranchAddress("eventW_Fastsimup", &_event_weight);
  else if (systematic=="Fastsimdo")
    fChain->SetBranchAddress("eventW_Fastsimdo", &_event_weight);
  else if (systematic=="Toppt" || systematic=="Topptup" || systematic=="Topptdo") {
    fChain->SetBranchAddress("eventW_Toppt",     &_event_weight);
    fChain->SetBranchAddress("eventW",           &_event_weight_Toppt);
  } else if (systematic=="Idisoeleup")
    fChain->SetBranchAddress("eventW_Idisoeleup",&_event_weight);
  else if (systematic=="Idisoeledo") 
    fChain->SetBranchAddress("eventW_Idisoeledo",&_event_weight);
  else if (systematic=="Idisomuup")
    fChain->SetBranchAddress("eventW_Idisomuup", &_event_weight);
  else if (systematic=="Idisomudo") 
    fChain->SetBranchAddress("eventW_Idisomudo", &_event_weight);
  else 
    fChain->SetBranchAddress("eventW",           &_event_weight);

  if (systematic=="nominal") {
    fChain->SetBranchAddress("eventW_Btagup",     &_event_weight_Btagup);
    fChain->SetBranchAddress("eventW_Btagdo",     &_event_weight_Btagdo);
    fChain->SetBranchAddress("eventW_Btaglightup",&_event_weight_Btaglightup);
    fChain->SetBranchAddress("eventW_Btaglightdo",&_event_weight_Btaglightdo);
    fChain->SetBranchAddress("eventW_BtagFSup",   &_event_weight_BtagFSup);
    fChain->SetBranchAddress("eventW_BtagFSdo",   &_event_weight_BtagFSdo);
    fChain->SetBranchAddress("eventW_Idisoup",    &_event_weight_Idisoup);
    fChain->SetBranchAddress("eventW_Idisodo",    &_event_weight_Idisodo);
    fChain->SetBranchAddress("eventW_Idisoeleup", &_event_weight_Idisoeleup);
    fChain->SetBranchAddress("eventW_Idisoeledo", &_event_weight_Idisoeledo);
    fChain->SetBranchAddress("eventW_Idisomuup",  &_event_weight_Idisomuup);
    fChain->SetBranchAddress("eventW_Idisomudo",  &_event_weight_Idisomudo);
    fChain->SetBranchAddress("eventW_Triggerup",  &_event_weight_Triggerup);
    fChain->SetBranchAddress("eventW_Triggerdo",  &_event_weight_Triggerdo);
    fChain->SetBranchAddress("eventW_Recoup",     &_event_weight_Recoup);
    fChain->SetBranchAddress("eventW_Recodo",     &_event_weight_Recodo);
    fChain->SetBranchAddress("eventW_Fastsimup",  &_event_weight_Fastsimup);
    fChain->SetBranchAddress("eventW_Fastsimdo",  &_event_weight_Fastsimdo);
    fChain->SetBranchAddress("eventW_Toppt",      &_event_weight_Toppt);
  }

  _event_weight_Isrnjetup = _event_weight;
  _event_weight_Isrnjetdo = _event_weight;

  fChain->SetBranchAddress("channel",         &_channel);
  fChain->SetBranchAddress("njet",            &_njet); 

  fChain->SetBranchAddress("nbjet30csvv2l",   &_nbjet30csvv2l);
  fChain->SetBranchAddress("nbjet30csvv2m",   &_nbjet30csvv2m);
  fChain->SetBranchAddress("nbjet30csvv2t",   &_nbjet30csvv2t);
  fChain->SetBranchAddress("leadingPtCSVv2M", &_leadingPtCSVv2M);
  fChain->SetBranchAddress("leadingPtCSVv2L", &_leadingPtCSVv2L);
  fChain->SetBranchAddress("leadingPtCSVv2T", &_leadingPtCSVv2T);
  fChain->SetBranchAddress("trailingPtCSVv2M", &_trailingPtCSVv2M);
  fChain->SetBranchAddress("trailingPtCSVv2L", &_trailingPtCSVv2L);
  fChain->SetBranchAddress("trailingPtCSVv2T", &_trailingPtCSVv2T);

  fChain->SetBranchAddress("mt2ll",           &_mt2ll);
  fChain->SetBranchAddress("mt2llgen",        &_mt2llgen);
  fChain->SetBranchAddress("mt2bb",           &_mt2bb);
  fChain->SetBranchAddress("mt2lblb",         &_mt2lblb);
  fChain->SetBranchAddress("mt2bbtrue",       &_mt2bbtrue);
  fChain->SetBranchAddress("mt2lblbcomb",     &_mt2lblbcomb);
  fChain->SetBranchAddress("mt2lblbtrue",     &_mt2lblbtrue);

  fChain->SetBranchAddress("mlb1",            &_mlb1);
  fChain->SetBranchAddress("mlb2",            &_mlb2);
  fChain->SetBranchAddress("mlb1comb",        &_mlb1comb);
  fChain->SetBranchAddress("mlb2comb",        &_mlb2comb);
  fChain->SetBranchAddress("mlb1true",        &_mlb1true);
  fChain->SetBranchAddress("mlb2true",        &_mlb2true);

  fChain->SetBranchAddress("bjet1pt",         &_bjet1pt);
  fChain->SetBranchAddress("bjet1eta",        &_bjet1eta);
  fChain->SetBranchAddress("bjet1phi",        &_bjet1phi);
  fChain->SetBranchAddress("bjet1mass",       &_bjet1mass);
  fChain->SetBranchAddress("bjet1csvv2ivf",   &_bjet1csvv2ivf);
  fChain->SetBranchAddress("bjet2pt",         &_bjet2pt);
  fChain->SetBranchAddress("bjet2eta",        &_bjet2eta);
  fChain->SetBranchAddress("bjet2phi",        &_bjet2phi);
  fChain->SetBranchAddress("bjet2mass",       &_bjet2mass);
  fChain->SetBranchAddress("bjet2csvv2ivf",   &_bjet2csvv2ivf);
  fChain->SetBranchAddress("tjet1pt",         &_tjet1pt);
  fChain->SetBranchAddress("tjet1eta",        &_tjet1eta);
  fChain->SetBranchAddress("tjet1phi",        &_tjet1phi);
  fChain->SetBranchAddress("tjet1mass",       &_tjet1mass);
  fChain->SetBranchAddress("tjet1csvv2ivf",   &_tjet1csvv2ivf);
  fChain->SetBranchAddress("tjet1assignment", &_tjet1assignment);
  fChain->SetBranchAddress("tjet2pt",         &_tjet2pt);
  fChain->SetBranchAddress("tjet2eta",        &_tjet2eta);
  fChain->SetBranchAddress("tjet2phi",        &_tjet2phi);
  fChain->SetBranchAddress("tjet2mass",       &_tjet2mass);
  fChain->SetBranchAddress("tjet2csvv2ivf",   &_tjet2csvv2ivf);
  fChain->SetBranchAddress("tjet2assignment", &_tjet2assignment);

  fChain->SetBranchAddress("MR",              &_MR);
  fChain->SetBranchAddress("R2",              &_R2);
  fChain->SetBranchAddress("Rpt",             &_Rpt);
  fChain->SetBranchAddress("invGamma",        &_invGamma);
  fChain->SetBranchAddress("Mdr",             &_Mdr);
  fChain->SetBranchAddress("DeltaPhiRll",     &_DeltaPhiRll);
  
  fChain->SetBranchAddress("dphijet1met", &_dphijet1met);
  fChain->SetBranchAddress("dphijet2met", &_dphijet2met);//, &b_dphijet2met);
  fChain->SetBranchAddress("dphijj", &_dphijj);//, &b_dphijj);
  fChain->SetBranchAddress("dphijjmet", &_dphijjmet);//, &b_dphijjmet);
  fChain->SetBranchAddress("dphilep1jet1", &_dphilep1jet1);//, &b_dphilep1jet1);
  fChain->SetBranchAddress("dphilep1jet2", &_dphilep1jet2);//, &b_dphilep1jet2);
  fChain->SetBranchAddress("dphilep2jet1", &_dphilep2jet1);//, &b_dphilep2jet1);
  fChain->SetBranchAddress("dphilep2jet2", &_dphilep2jet2);//, &b_dphilep2jet2);
  fChain->SetBranchAddress("dphill", &dphill);//, &b_dphill);
  fChain->SetBranchAddress("dphillstar", &_dphillstar);//, &b_dphillstar);
  fChain->SetBranchAddress("dphilmet1", &dphilmet1);//, &b_dphilmet1);
  fChain->SetBranchAddress("dphilmet2", &dphilmet2);//, &b_dphilmet2);
  fChain->SetBranchAddress("drll", &drll);//, &b_drll);
  fChain->SetBranchAddress("dphimetptbll", &_dphimetptbll);//, &b_dphimetptbll);
  fChain->SetBranchAddress("dphimetbbll", &_dphimetbbll);//, &b_dphimetbbll);
  fChain->SetBranchAddress("dphimetjet", &_dphimetjet);//, &b_dphimetjet);

  if (!_isfastsim) fChain->SetBranchAddress("LHEweight", &std_vector_LHE_weight);
  fChain->SetBranchAddress("jet_pt",  &std_vector_jet_pt);
  fChain->SetBranchAddress("jet_eta", &std_vector_jet_eta);
  fChain->SetBranchAddress("jet_phi", &std_vector_jet_phi);

  fChain->SetBranchAddress("ZZmass", &_ZZmass);
  fChain->SetBranchAddress("ZZdphi", &_ZZdphi);
  fChain->SetBranchAddress("ZZpt",   &_ZZpt);

}
    
//------------------------------------------------------------------------------
// ShapeWZtoWW
//------------------------------------------------------------------------------
bool AnalysisStop::ShapeWZtoWW()
{
  
  if (_nlepton<3) return false;
  //if (MET.Et()<100.) return false;

  int nTightLeptons = 3;
  if (_systematic.Contains("WZtoWWloose")) nTightLeptons = 2;
  if (_ntightlepton<nTightLeptons) return false;

  if (_systematic.Contains("veto") && _nlepton>3) return false;

  float DMZCut = (_systematic.Contains("noZcut")) ? 999. : 10.;
  if (_systematic.Contains("DM15")) DMZCut = 15.;

  float DMZ = 999.; int Wlep1 = -1, Wlep2 = -1, Lostlep = -1;

  for (int l1 = 0; l1<3; l1++) {
    for (int l2 = l1+1; l2<3; l2++) {
      if (AnalysisLeptons[l1].flavour*AnalysisLeptons[l2].flavour<0. && 
	  fabs(AnalysisLeptons[l1].flavour)==fabs(AnalysisLeptons[l2].flavour)) {

	TLorentzVector Zcand = AnalysisLeptons[l1].v + AnalysisLeptons[l2].v;
	float DMZcand = fabs( Zcand.M() - Z_MASS );
	
	if (DMZcand<DMZCut && _systematic.Contains("WZtoWWlooseZ") && _ntightlepton<3)
	  if (AnalysisLeptons[l1].type==1 && AnalysisLeptons[l2].type==1) return false;

	if (DMZcand<DMZ) {
	  
	  if (l1==0 && l2==1) {
	    Wlep2 = 2;
	    if (AnalysisLeptons[l1].flavour*AnalysisLeptons[Wlep2].flavour<0.) { Wlep1 = l1; Lostlep = l2; }
	    else if (AnalysisLeptons[l2].flavour*AnalysisLeptons[Wlep2].flavour<0.) { Wlep1 = l2; Lostlep = l1; }
	  } if (l1==0 && l2==2) {
	    if (AnalysisLeptons[l1].flavour*AnalysisLeptons[1].flavour<0.) { Wlep1 = l1; Wlep2 = 1; Lostlep = l2; }
	    else if (AnalysisLeptons[l2].flavour*AnalysisLeptons[1].flavour<0.) { Wlep1 = 1; Wlep2 = l2; Lostlep = l1; }
	  } if (l1==1 && l2==2) {
	    Wlep1 = 0;
	    if (AnalysisLeptons[l1].flavour*AnalysisLeptons[Wlep1].flavour<0.) { Wlep2 = l1; Lostlep = l2; }
	    else if (AnalysisLeptons[l2].flavour*AnalysisLeptons[Wlep1].flavour<0.) { Wlep2 = l2; Lostlep = l1; }
	  }

	  _massZcandidate =  Zcand.M();
	  
	  DMZ = DMZcand;
	  
	}
	
      }
    }
  }
    
  if (DMZ>=DMZCut || Wlep1<0 || Wlep2<0 || Lostlep<0) { 
    return false; 
  } else {

    if (AnalysisLeptons[Wlep1].v.Pt()<25. || AnalysisLeptons[Wlep2].v.Pt()<20.) return false;
    //if (AnalysisLeptons[Wlep1].v.Pt()<25. && AnalysisLeptons[Wlep2].v.Pt()<25.) return false;
    //if (AnalysisLeptons[Wlep1].v.Pt()<20. || AnalysisLeptons[Wlep2].v.Pt()<20.) return false;

    TVector3 NewMET; NewMET.SetXYZ(MET.Px() + AnalysisLeptons[Lostlep].v.Px(),
				   MET.Py() + AnalysisLeptons[Lostlep].v.Py(),
				   MET.Pz());
    MET.SetPtEtaPhiM(NewMET.Pt(), 0.0, NewMET.Phi(), 0.0);
    
    if (Wlep2==2) {
      
      if (Wlep1==1) AnalysisLeptons[0] = AnalysisLeptons[1]; 
      AnalysisLeptons[1] = AnalysisLeptons[2];

      Lepton1 = AnalysisLeptons[0];
      Lepton2 = AnalysisLeptons[1];

      mll = (Lepton1.v+Lepton2.v).M();
      _m2l  = mll;
      
      if (abs(Lepton1.flavour)!=abs(Lepton2.flavour)) _channel = em;
      else if (abs(Lepton1.flavour)==11) _channel = ee;
      else if (abs(Lepton1.flavour)==13) _channel = mm;
      
      _mt2ll    = ComputeMT2(Lepton1.v, Lepton2.v, MET);
      _MT2ll = (_mt2ll<140.) ? _mt2ll : 139.;
     
      _lep1pt = Lepton1.v.Pt();   
      _lep2pt = Lepton2.v.Pt();

      _lep1phi = Lepton1.v.Phi();   
      _lep2phi = Lepton2.v.Phi();  
  
      dphill = fabs((Lepton1.v).DeltaPhi(Lepton2.v));
      dphilmet1 = fabs((Lepton1.v).DeltaPhi(MET)); // this recompute the latino variable used in AnalysisCMS.C
      dphilmet2 = fabs((Lepton2.v).DeltaPhi(MET)); // this recompute the latino variable used in AnalysisCMS.C
      _dphillmet = fabs((Lepton1.v + Lepton2.v).DeltaPhi(MET));

    
    }
  
    return true;
  
  }

}

//------------------------------------------------------------------------------
// ShapeZZ
//------------------------------------------------------------------------------
bool AnalysisStop::ShapeZZ()
{
  
  if (_nlepton<4) return false;

  int nTightLeptons = 3;
  if (_systematic.Contains("ZZtight")) nTightLeptons = 4;
  if (_systematic.Contains("ZZloose")) nTightLeptons = 2;
  if (_ntightlepton<nTightLeptons) return false;
 
  float DMZ1 = 15., DMZ2 = 30., DMZT = 999.; int Z1lep1 = -1, Z1lep2 = -1, Z2lep1 = -1, Z2lep2 = -1;
  for (int l1 = 0; l1<4; l1++) {
    for (int l2 = l1+1; l2<4; l2++) {
      if (AnalysisLeptons[l1].flavour*AnalysisLeptons[l2].flavour<0. && 
	  fabs(AnalysisLeptons[l1].flavour)==fabs(AnalysisLeptons[l2].flavour)) { 
	
	TLorentzVector Zcand1 = AnalysisLeptons[l1].v + AnalysisLeptons[l2].v;
	float DMZcand1 = fabs( Zcand1.M() - Z_MASS );
	if (DMZcand1<DMZ1) {
	
	  for (int l3 = 0; l3<4; l3++) {
	    if (l3!=l1 && l3!=l2) {
	      for (int l4 = l3+1; l4<4; l4++) {
		if (l4!=l1 && l4!=l2) {
		  if (AnalysisLeptons[l3].flavour*AnalysisLeptons[l4].flavour<0. && 
		      fabs(AnalysisLeptons[l3].flavour)==fabs(AnalysisLeptons[l4].flavour)) {
		    
		    TLorentzVector Zcand2 = AnalysisLeptons[l3].v + AnalysisLeptons[l4].v;
		    float DMZcand2 = fabs( Zcand2.M() - Z_MASS );
		    if (DMZcand2<DMZ2) {
		      
		      float DMZTcand = sqrt(DMZcand1*DMZcand1 + DMZcand2*DMZcand2);
		      if (DMZTcand<DMZT) {
			
			bool EventOdd = event%2;
			if (DMZcand2>DMZ1 || EventOdd) {
			  Z1lep1 = l1;
			  Z1lep2 = l2;
			  Z2lep1 = l3;
			  Z2lep2 = l4;
			} else {
			  Z1lep1 = l3;
			  Z1lep2 = l4;
			  Z2lep1 = l1;
			  Z2lep2 = l2;
			}
		      
			DMZT = DMZTcand;
			
		      }

		    }

		  }
		}
	      }
	    }
	  }

	}

      }
    }
  }

  if (DMZT>900. || Z1lep1<0 || Z1lep2<0 || Z2lep1<0 || Z2lep2<0) { 
    return false;
  } else {

    Lepton1 = AnalysisLeptons[Z1lep1]; 
    Lepton2 = AnalysisLeptons[Z1lep2]; 
    Lepton Lepton3 = AnalysisLeptons[Z2lep1];
    Lepton Lepton4 = AnalysisLeptons[Z2lep2];

    //if (abs(Lepton1.flavour)!=abs(Lepton3.flavour)) _channel = em;
    //else if (abs(Lepton1.flavour)==11) _channel = ee;
    //else if (abs(Lepton1.flavour)==13) _channel = mm;

    int nEle = 0;
    if (abs(Lepton1.flavour)==11) nEle++;
    if (abs(Lepton2.flavour)==11) nEle++;
    if (abs(Lepton3.flavour)==11) nEle++;
    if (abs(Lepton4.flavour)==11) nEle++;
    if (nEle==0) _channel = mm;
    else if (nEle==2) _channel = em;
    else if (nEle==4) _channel = ee;
    else return false;

    /*
    if (abs(Lepton1.flavour)!=abs(Lepton3.flavour)) {
      Lepton Lepton5 = Lepton1, Lepton6, Lepton7 = Lepton2, Lepton8;
      if (Lepton2.flavour*Lepton3.flavour<0) {
	Lepton8 = Lepton3;
	Lepton6 = Lepton4;
      } else {
	Lepton8 = Lepton4;
	Lepton6 = Lepton3;
      }
      float DMZswitch1 = fabs( (Lepton5.v+Lepton6.v).M() - Z_MASS );
      float DMZswitch2 = fabs( (Lepton7.v+Lepton8.v).M() - Z_MASS );
      bool AllSwitch = true;
      if (AllSwitch || (DMZswitch1<DMZ2 && DMZswitch2<DMZ2 && (DMZswitch1<DMZ1 || DMZswitch2<DMZ1))) {
	TVector3 PreMET; PreMET.SetXYZ(MET.Px() + Lepton3.v.Px() + Lepton4.v.Px(),
				       MET.Py() + Lepton3.v.Py() + Lepton4.v.Py(),
				       MET.Pz());
	bool EventOdd = event%2;
	if ((!AllSwitch && (DMZswitch2>DMZ1 || (DMZswitch1<DMZ1 && EventOdd))) || (AllSwitch && DMZswitch1<DMZswitch2)) {
	  Lepton1 = Lepton5;
	  Lepton2 = Lepton6;
	  Lepton3 = Lepton7;
	  Lepton4 = Lepton8;
	} else {
	  Lepton1 = Lepton7;
	  Lepton2 = Lepton8;
	  Lepton3 = Lepton5;
	  Lepton4 = Lepton6;
	}
	TVector3 PostMET; PostMET.SetXYZ(MET.Px() + Lepton3.v.Px() + Lepton4.v.Px(),
					 MET.Py() + Lepton3.v.Py() + Lepton4.v.Py(),
					 MET.Pz());
	if (PreMET.Mag()>140. && PostMET.Mag()<140.)
	  cout << "  switch lost " << PreMET.Mag() << " " << PostMET.Mag() << " " << _njet << endl;
	if (PreMET.Mag()<140. && PostMET.Mag()>140.)
	  cout << "  switch got  " << PreMET.Mag() << " " << PostMET.Mag() << " " << _njet << endl;
      }
    }
    */

    if (Lepton1.v.Pt()<25. && Lepton2.v.Pt()<25.) return false;
    //if (Lepton1.v.Pt()<20. || Lepton2.v.Pt()<20. || Lepton3.v.Pt()<20. || Lepton4.v.Pt()<20.) return false;
    
    TVector3 NewMET; NewMET.SetXYZ(MET.Px() + Lepton3.v.Px() + Lepton4.v.Px(),
				   MET.Py() + Lepton3.v.Py() + Lepton4.v.Py(),
				   MET.Pz());
    MET.SetPtEtaPhiM(NewMET.Pt(), 0.0, NewMET.Phi(), 0.0);
    
    AnalysisLeptons[0] = Lepton1;
    AnalysisLeptons[1] = Lepton2;
    AnalysisLeptons[2] = Lepton3;
    AnalysisLeptons[3] = Lepton4;
    
    int _mjet = 0;
    for (int ij = 0; ij<_jet_pt.size(); ij++) {
      bool goodjet = true;
      TLorentzVector jetv; jetv.SetPtEtaPhiM(_jet_pt.at(ij), _jet_eta.at(ij), _jet_phi.at(ij), 0.);
      if (Lepton1.v.DeltaR(jetv)<0.4) goodjet = false;
      if (Lepton2.v.DeltaR(jetv)<0.4) goodjet = false;
      if (Lepton3.v.DeltaR(jetv)<0.4) goodjet = false;
      if (Lepton4.v.DeltaR(jetv)<0.4) goodjet = false;
      if (goodjet) _mjet++;
    }
    _njet = _mjet;
    
    _event_weight *= (_lep3idisoW*_lep4idisoW);
    
    mll = (Lepton1.v+Lepton2.v).M();
    _m2l  = mll;
    _m2lZ2 = (Lepton3.v+Lepton4.v).M();
    
    _mt2ll    = ComputeMT2(Lepton1.v, Lepton2.v, MET);
    _MT2ll = (_mt2ll<140.) ? _mt2ll : 139.;
  
    return true;
    
  }

}


//------------------------------------------------------------------------------
// ShapeZWtoZ
//------------------------------------------------------------------------------
bool AnalysisStop::ShapeZWtoZ()
{
  
  if (_nlepton<3) return false; 
  //if (MET.Et()<100.) return false;

  int nTightLeptons = 3;
  if (_systematic.Contains("loose")) nTightLeptons = 2;
  if (_ntightlepton<nTightLeptons) return false;

  float DMZCut = 10.;
  if (_systematic.Contains("DM15")) DMZCut = 15.;

  float DMZ = 999.; int Zlep1 = -1, Zlep2 = -1, Lostlep = -1;

  for (int l1 = 0; l1<3; l1++) {
    for (int l2 = l1+1; l2<3; l2++) {
      if (AnalysisLeptons[l1].flavour*AnalysisLeptons[l2].flavour<0. && 
	  fabs(AnalysisLeptons[l1].flavour)==fabs(AnalysisLeptons[l2].flavour)) {

	TLorentzVector Zcand = AnalysisLeptons[l1].v + AnalysisLeptons[l2].v;
	float DMZcand = fabs( Zcand.M() - Z_MASS );
	
	if (DMZcand<DMZCut && _systematic.Contains("looseZ") && _ntightlepton<3)
	  if (AnalysisLeptons[l1].type==1 && AnalysisLeptons[l2].type==1) return false;

	if (DMZcand<DMZ) {
	  
	  Zlep1 = l1; 
	  Zlep2 = l2;

	  if (l1==0 && l2==1) Lostlep = 2;
	  if (l1==0 && l2==2) Lostlep = 1;
	  if (l1==1 && l2==2) Lostlep = 0;
	  
	  DMZ = DMZcand;
	  
	}
	
      }
    }
  }
    
  if (DMZ>=DMZCut || Zlep1<0 || Zlep2<0 || Lostlep<0) { 
    return false; 
  } else {
    
    Lepton1 = AnalysisLeptons[Zlep1]; 
    Lepton2 = AnalysisLeptons[Zlep2];
    Lepton Lepton3 =  AnalysisLeptons[Lostlep];

    if (Lepton1.v.Pt()<25. && Lepton2.v.Pt()<25.) return false;

    TVector3 NewMET; NewMET.SetXYZ(MET.Px() + Lepton3.v.Px(),
				   MET.Py() + Lepton3.v.Py(),
				   MET.Pz());
    MET.SetPtEtaPhiM(NewMET.Pt(), 0.0, NewMET.Phi(), 0.0);
    
    AnalysisLeptons[0] = Lepton1;
    AnalysisLeptons[1] = Lepton2;
    AnalysisLeptons[2] = Lepton3;

    mll = (Lepton1.v+Lepton2.v).M();
    _m2l  = mll;

    if (abs(Lepton1.flavour)!=abs(Lepton3.flavour)) _channel = em;
    else if (abs(Lepton1.flavour)==11) _channel = ee;
    else if (abs(Lepton1.flavour)==13) _channel = mm;
    
    _mt2ll    = ComputeMT2(Lepton1.v, Lepton2.v, MET);
    _MT2ll = (_mt2ll<140.) ? _mt2ll : 139.;
  
    return true;
    
  }

}
