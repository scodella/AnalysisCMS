#include "HistogramReader.h"

using namespace std;


//------------------------------------------------------------------------------
// HistogramReader
//------------------------------------------------------------------------------
HistogramReader::HistogramReader(const TString& inputdir,
				 const TString& outputdir) :
  _inputdir        (inputdir),
  _outputdir       (outputdir),
  _stackoption     ("nostack,hist"),
  _title           ("inclusive"),
  _luminosity_fb   (-1),
  _lumilegend_fb   (-1),
  _datanorm        (false),
  _drawratio       (false),
  _drawsignificance(false),
  _drawyield       (false),
  _minitreebased   (false),
  _publicstyle     (false),
  _savepdf         (false),
  _savepng         (true)
{
  _mcfile.clear();
  _mccolor.clear();
  _mclabel.clear();
  _mcscale.clear();

  _datafile  = NULL;
  _datahist  = NULL;
  _allmchist = NULL;

  TH1::SetDefaultSumw2();
}


//------------------------------------------------------------------------------
// AddData
//------------------------------------------------------------------------------
void HistogramReader::AddData(const TString& filename,
			      const TString& label,
			      Color_t        color)
{
  TString fullname = _inputdir + "/" + filename + ".root";

  if (gSystem->AccessPathName(fullname))
    {
      printf(" [HistogramReader::AddData] Cannot access %s\n", fullname.Data());
      return;
    }

  TFile* file = new TFile(fullname, "update");

  _datacolor    = color;
  _datafile     = file;
  _datafilename = filename;
  _datalabel    = label;
}


//------------------------------------------------------------------------------
// AddProcess
//------------------------------------------------------------------------------
void HistogramReader::AddProcess(const TString& filename,
				 const TString& label,
				 Color_t        color,
				 Int_t          kind,
				 Float_t        scale)
{
  TString fullname = _inputdir + "/" + filename + ".root";

  if (gSystem->AccessPathName(fullname))
    {
      printf(" [HistogramReader::AddProcess] Cannot access %s\n", fullname.Data());
      return;
    }

  TFile* file = new TFile(fullname, "update");

  _mccolor.push_back(color);
  _mcfile.push_back(file);
  _mcfilename.push_back(filename); 
  _mclabel.push_back(label);
  _mcscale.push_back(scale);
  
  if (scale > 0. && scale != 1.)
    printf("\n [HistogramReader::AddProcess] Process %s will be scaled by %.2f\n\n", label.Data(), scale);

  if (kind == roc_signal)
    {
      _roc_signalfile.push_back(file);
      _roc_signalscale.push_back(scale);
    }
  else if (kind == roc_background)
    {
      _roc_backgroundfile.push_back(file);
      _roc_backgroundscale.push_back(scale);
    }
}


//------------------------------------------------------------------------------
// AddSignal
//------------------------------------------------------------------------------
void HistogramReader::AddSignal(const TString& filename,
				const TString& label,
				Color_t        color,
				Int_t          kind,
				Float_t        scale)
{
  TString fullname = _inputdir + "/" + filename + ".root";
  
  if (gSystem->AccessPathName(fullname))
    {
      printf(" [HistogramReader::AddSignal] Cannot access %s\n", fullname.Data());
      return;
    }

  TFile* file = new TFile(fullname, "update");

  _signalcolor.push_back(color);
  _signalfile.push_back(file);
  _signalfilename.push_back(filename);
  _signallabel.push_back(label);
  _signalscale.push_back(scale);

  if (scale > 0. && scale != 1.)
    printf("\n [HistogramReader::AddSignal] Process %s will be scaled by %.2f\n\n", label.Data(), scale);

  if (kind == roc_signal)
    {
      _roc_signalfile.push_back(file);
      _roc_signalscale.push_back(scale);
    }
  else if (kind == roc_background)
    {
      _roc_backgroundfile.push_back(file);
      _roc_backgroundscale.push_back(scale);
    }
}

//------------------------------------------------------------------------------
// AddSystematic
//------------------------------------------------------------------------------
void HistogramReader::AddSystematic(TString analysis, TString systematic)
{
  //_mycut    = mycut;
  _analysis = analysis;
  _systematics.push_back(systematic);
}

//------------------------------------------------------------------------------
// Draw
//------------------------------------------------------------------------------
void HistogramReader::Draw(TString hname,
			   TString xtitle,
			   Int_t   ngroup,
			   Int_t   precision,
			   TString units,
			   Bool_t  setlogy,
			   Bool_t  moveoverflow,
			   Float_t xmin,
			   Float_t xmax,
			   Float_t ymin,
			   Float_t ymax)
{
  TString cname = hname;

  if (_stackoption.Contains("nostack")) cname += "_nostack";

  if (setlogy) cname += "_log";

  _writeyields = (hname.Contains("_evolution")) ? true : false;

  if (_writeyields)
    {
      _yields_table.open(_outputdir + "/" + cname + ".txt");

      _writelabels = true;
    }


  TCanvas* canvas = NULL;

  TPad* pad1 = NULL;
  TPad* pad2 = NULL;


  // Set drawsignificance to false if drawratio is true
  if (_drawratio /*&& _datafile*/) _drawsignificance = false;


  if ((_drawratio /*&& _datafile*/) || _drawsignificance)
    {
      canvas = new TCanvas(cname, cname, 550, 720);

      pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
      pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);

      pad1->SetTopMargin   (0.08);
      pad1->SetBottomMargin(0.02);
      pad1->Draw();

      pad2->SetTopMargin   (0.08);
      pad2->SetBottomMargin(0.35);
      pad2->Draw();
    }
  else
    {
      canvas = new TCanvas(cname, cname, 550, 600);

      pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);

      pad1->Draw();
    }


  //----------------------------------------------------------------------------
  // pad1
  //----------------------------------------------------------------------------
  pad1->cd();
  
  pad1->SetLogy(setlogy);


  // Stack processes
  //----------------------------------------------------------------------------
  _mchist.clear();

  THStack* mcstack = new THStack(hname + "_mcstack", hname + "_mcstack");

  for (UInt_t i=0; i<_mcfile.size(); i++) {

    _mcfile[i]->cd();

    TH1D* dummy = GetHistogram(_mcfile[i], hname);//(TH1D*)_mcfile[i]->Get(hname);

    _mchist.push_back((TH1D*)dummy->Clone());

    if (_luminosity_fb > 0 && _mcscale[i] > -999) _mchist[i]->Scale(_luminosity_fb);

    if (_mcscale[i] > 0) _mchist[i]->Scale(_mcscale[i]);

    SetHistogram(_mchist[i], _mccolor[i], 1001, kDot, kSolid, 0, ngroup, moveoverflow, xmin, xmax);
    
    mcstack->Add(_mchist[i]);
  }


  // Stack signals
  //----------------------------------------------------------------------------
  _signalhist.clear();

  THStack* signalstack = new THStack(hname + "_signalstack", hname + "_signalstack");

  for (UInt_t i=0; i<_signalfile.size(); i++) {

    _signalfile[i]->cd();

    TH1D* dummy = GetHistogram(_signalfile[i], hname);//(TH1D*)_signalfile[i]->Get(hname);
    
    if (hname.Contains("_SR") && (hname.Contains("h_MT2ll_") || hname.Contains("h_MT2llisr_"))) {
      for (int isyst=0; isyst<_systematics.size(); isyst++) {
	if (_systematics.at(isyst)=="Metfastsim")  {
      
	  TString hnamegen = hname;
	  hnamegen.ReplaceAll("_SR1", "_SR1gen");
	  hnamegen.ReplaceAll("_SR2", "_SR2gen");
	  hnamegen.ReplaceAll("_SR3", "_SR3gen");
	  hnamegen.ReplaceAll("h_MT2ll_", "h_MT2llgen_");
	  hnamegen.ReplaceAll("h_MT2llisr_", "h_MT2llisrgen_");
	  TH1D* dummygen = GetHistogram(_signalfile[i], hnamegen);//(TH1D*)_signalfile[i]->Get(hnamegen);
	  dummy->Add(dummygen); // that's not good if errors are needed
	  dummy->Scale(0.5);
	  
	}
      }
    }

    _signalhist.push_back((TH1D*)dummy->Clone());

    if (_luminosity_fb > 0 && _signalscale[i] > -999) _signalhist[i]->Scale(_luminosity_fb);

    if (_signalscale[i] > 0) _signalhist[i]->Scale(_signalscale[i]);

    SetHistogram(_signalhist[i], _signalcolor[i], 0, kDot, kSolid, 3, ngroup, moveoverflow, xmin, xmax);
    
    signalstack->Add(_signalhist[i]);
  }


  // Get the data
  //----------------------------------------------------------------------------
  if (_datafile)
    {
      _datafile->cd();

      TH1D* dummy = GetHistogram(_datafile, hname);//(TH1D*)_datafile->Get(hname);

      _datahist = (TH1D*)dummy->Clone();
      
      SetHistogram(_datahist, kBlack, 0, kFullCircle, kSolid, 1, ngroup, moveoverflow, xmin, xmax);
    }


  // Normalize MC to data
  //----------------------------------------------------------------------------
  if (_datahist && _datanorm)
    {
      Float_t mcnorm   = Yield((TH1D*)(mcstack->GetStack()->Last()));
      Float_t datanorm = Yield(_datahist);

      for (UInt_t i=0; i<_mchist.size(); i++)
	{
	  _mchist[i]->Scale(datanorm / mcnorm);
	}

      mcstack->Modified();
    }


  // hfirst will contain the axis settings
  //----------------------------------------------------------------------------
  TH1D* hfirst = (TH1D*)_mchist[0]->Clone("hfirst");

  hfirst->Reset();

  hfirst->SetTitle("");


  // All MC
  //----------------------------------------------------------------------------
  _allmchist = (TH1D*)_mchist[0]->Clone("allmchist");

  _allmchist->SetName(_mchist[0]->GetName());

  // Possible modification (how to deal with systematic uncertainties?)
  //  _allmchist = (TH1D*)(mcstack->GetStack()->Last());


  // Include systematics
  //----------------------------------------------------------------------------
  //if (_mchist_syst.size() > 0) IncludeSystematics(hname);


  for (Int_t ibin=0; ibin<=_allmchist->GetNbinsX(); ibin++) {

    Float_t binValue = 0.;
    Float_t binError = 0.;

    for (UInt_t i=0; i<_mchist.size(); i++) {

      Float_t binContent   = _mchist[i]->GetBinContent(ibin);
      Float_t binStatError = sqrt(_mchist[i]->GetSumw2()->At(ibin));
      Float_t binSystError = (_mchist_syst.size() > 0) ? _mchist_syst[i]->GetBinContent(ibin) : 0.;

      binValue += binContent;
      binError += (binStatError * binStatError);
      binError += (binSystError * binSystError);
    }

    binError = sqrt(binError);

    _allmchist->SetBinContent(ibin, binValue);
    _allmchist->SetBinError  (ibin, binError);

  }

  _allmclabel = "stat";

  _allmchist->SetFillColor  (kGray+1);
  _allmchist->SetFillStyle  (   3345);
  _allmchist->SetLineColor  (kGray+1);
  _allmchist->SetMarkerColor(kGray+1);
  _allmchist->SetMarkerSize (      0);


  // Include systematics with TGraphAsymmErrors
  //----------------------------------------------------------------------------
  if (_systematics.size() > 0)  IncludeSystematics(hname);
  

  // Draw
  //----------------------------------------------------------------------------
  hfirst->Draw();

  mcstack->Draw(_stackoption + ",same");

  if (_systematics.size() > 0 ) { 
    _ErrorGr->Draw("e2,same");  
  } else { 
    if (!_stackoption.Contains("nostack")) _allmchist->Draw("e2,same");
  }

  if (_signalfile.size() > 0) signalstack->Draw("nostack,hist,same");

  if (_datahist) _datahist->Draw("ep,same");


  // Set xtitle and ytitle
  //----------------------------------------------------------------------------
  TString ytitle = Form("events / %s.%df", "%", precision);

  ytitle = Form(ytitle.Data(), hfirst->GetBinWidth(0));

  if (!units.Contains("NULL")) {
    xtitle = Form("%s [%s]", xtitle.Data(), units.Data());
    ytitle = Form("%s %s",   ytitle.Data(), units.Data());
  }


  // Adjust xaxis and yaxis
  //----------------------------------------------------------------------------
  hfirst->GetXaxis()->SetRangeUser(xmin, xmax);

  Float_t theMin = 0.0;

  Float_t theMax = (_datahist) ? GetMaximum(_datahist, xmin, xmax) : 0.0;
  
  Float_t theMaxMC = GetMaximum(_allmchist, xmin, xmax);

  if (_stackoption.Contains("nostack"))
    {
      for (UInt_t i=0; i<_mcfile.size(); i++)
	{
	  Float_t mchist_i_max = GetMaximum(_mchist[i], xmin, xmax, false);

	  if (mchist_i_max > theMaxMC) theMaxMC = mchist_i_max;
	}
    }
  
  if (theMaxMC > theMax) theMax = theMaxMC;

  Float_t theMaxSignal = 0.0;

  if (_signalfile.size() > 0)
    {
      for (UInt_t i=0; i<_signalfile.size(); i++)
	{
	  Float_t signalhist_i_max = GetMaximum(_signalhist[i], xmin, xmax, false);

	  if (signalhist_i_max > theMaxSignal) theMaxSignal = signalhist_i_max;
	}
    }

  if (theMaxSignal > theMax) theMax = theMaxSignal;
  
  if (pad1->GetLogy())
    {
      theMin = 1e-2; // 1e-5
      theMax = TMath::Power(10, TMath::Log10(theMax) + 4); // 6);
    }
  else if (!_stackoption.Contains("nostack"))
    {
      theMax *= 1.7;
    } else theMax = 0.8;

  hfirst->SetMinimum(theMin);
  hfirst->SetMaximum(theMax);

  if (ymin != -999) hfirst->SetMinimum(ymin);
  if (ymax != -999) hfirst->SetMaximum(ymax);


  // Legend
  //----------------------------------------------------------------------------
  Float_t x0     = 0.220;                         // x position of the data on the top left
  Float_t y0     = 0.843;                         // y position of the data on the top left
  Float_t xdelta = (_drawyield) ? 0.228 : 0.228;//0.170;  // x width between columns
  Float_t ydelta = 0.050;                         // y width between rows
  Int_t   nx     = 0;                             // column number
  Int_t   ny     = 0;                             // row    number

  TString opt = (_stackoption.Contains("nostack")) ? "l" : "f";


  // Data legend
  //----------------------------------------------------------------------------
  if (_datahist)
    {
      DrawLegend(x0, y0, _datahist, _datalabel.Data(), "lp");
      ny++;
    }


  // All MC legend
  //----------------------------------------------------------------------------
  if (!_stackoption.Contains("nostack"))
    {
      DrawLegend(x0, y0 - ny*ydelta, _allmchist, _allmclabel.Data(), opt);
      ny++;
    }


  // Standard Model processes legend
  //----------------------------------------------------------------------------
  Int_t nrow = (_mchist.size() > 10) ? 5 : 4;

  for (int i=0; i<_mchist.size(); i++)
    {
      if (ny == nrow)
	{
	  ny = 0;
	  nx++;
	}

      DrawLegend(x0 + nx*xdelta, y0 - ny*ydelta, _mchist[i], _mclabel[i].Data(), opt);
      ny++;
    }


  // Search signals legend
  //----------------------------------------------------------------------------
  nx = 0; ny = nrow;
  for (int i=0; i<_signalhist.size(); i++)
    {
      if (ny == nrow+3)
	{
	  ny =  nrow;
	  nx++;
	}
      DrawLegend(x0 + nx*(xdelta*1.6), y0 - ny*ydelta, _signalhist[i], _signallabel[i].Data(), "l");
      ny++;
    }


  // Titles
  //----------------------------------------------------------------------------
  Float_t xprelim = ((_drawratio /*&& _datafile*/) || _drawsignificance) ? 0.288 : 0.300;

  if (_title.EqualTo("inclusive"))
    {
      DrawLatex(61, 0.190,   0.945, 0.050, 11, "CMS");
      DrawLatex(52, xprelim, 0.945, 0.030, 11, "Preliminary");
    }
  else
    {
      DrawLatex(42, 0.190, 0.945, 0.050, 11, _title);
    }

  if (_luminosity_fb > 0)
    DrawLatex(42, 0.940, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV)", _lumilegend_fb));
  else
    DrawLatex(42, 0.940, 0.945, 0.050, 31, "(13TeV)");

  SetAxis(hfirst, xtitle, ytitle, 1.5, 1.8);


      
  if (hname.Contains("h_dphiLL_llaa")) {
    TFile *OUT = new TFile("DYFit.root","recreate");
    OUT->cd();
    _datahist->SetName("data");
    _datahist->Write();
    _allmchist->SetName("MC");
    _allmchist->Write(); 
    for (UInt_t i=0; i<_mcfile.size(); i++) {
      if (_mcfilename[i].Contains("07")) {
	_mcfile[i]->cd();
	TH1D* dummy = (TH1D*)_mcfile[i]->Get(hname);
	dummy->Scale(_luminosity_fb);
	dummy->SetName("DY");
	OUT->cd();
	dummy->Write();
      }
    }
    OUT->Close();
  }
      
  //----------------------------------------------------------------------------
  // pad2
  //----------------------------------------------------------------------------
  if (_drawratio /*&& _datafile*/)
    {
      pad2->cd();
     
      // This approach isn't yet working
      //      TGraphAsymmErrors* g = new TGraphAsymmErrors();
      //      g->Divide(_datahist, _allmchist, "cl=0.683 b(1,1) mode");
      //      g->SetMarkerStyle(kFullCircle);
      //      g->Draw("ap");

      TH1D* ratio = _datahist ? (TH1D*)_datahist->Clone("ratio") : (TH1D*)_allmchist->Clone("ratio");
      TH1D* uncertainty = (TH1D*)_allmchist->Clone("uncertainty");
       
      float mcrenormalization = 1.;
      //if (_datahist && (hname.Contains("h_njet20") || hname.Contains("h_njet30"))) 
      //mcrenormalization = _datahist->Integral(0, ratio->GetNbinsX()+1)/_allmchist->Integral(0, ratio->GetNbinsX()+1);

      for (Int_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {
	
	Float_t dtValue = _datahist ? _datahist->GetBinContent(ibin) : -999.;
	Float_t dtError = _datahist ? _datahist->GetBinError(ibin)   :    0.;
	
	Float_t mcValue = _allmchist->GetBinContent(ibin);
	Float_t mcError = (_systematics.size()<=0) ? sqrt(_allmchist->GetSumw2()->At(ibin)) : _ErrorGr->GetErrorY(ibin);
	
	Float_t ratioVal         = 999;
	Float_t ratioErr         = 999;
	Float_t uncertaintyError = 999;

	if (mcValue > 0)
	  {
	    ratioVal         = dtValue / mcValue / mcrenormalization;
	    ratioErr         = dtError / mcValue;
	    uncertaintyError = _datahist ? ratioVal * mcError / mcValue : mcError / mcValue;
	  }

	ratio->SetBinContent(ibin, ratioVal);
	ratio->SetBinError  (ibin, ratioErr);
	
	uncertainty->SetBinContent(ibin, 1.);
	uncertainty->SetBinError  (ibin, uncertaintyError);
      }

      ratio->SetTitle("");

      ratio->Draw("ep");

      ratio->GetXaxis()->SetRangeUser(xmin, xmax);
      ratio->GetYaxis()->SetRangeUser(0., 2.);

      uncertainty->Draw("e2,same");
      
      ratio->Draw("ep,same");
 
      SetAxis(ratio, xtitle, "data / MC", 1.4, 0.75);
    }
  else if (_drawsignificance)
    {
      // Barbara's stuff
    }


  //----------------------------------------------------------------------------
  // Save it
  //----------------------------------------------------------------------------
  canvas->cd();

  if (_savepdf) canvas->SaveAs(_outputdir + cname + ".pdf");
  if (_savepng) canvas->SaveAs(_outputdir + cname + ".png");

  if (_writeyields)
    {
      _yields_table << std::endl;
      
      _yields_table.close();
    }
}


//------------------------------------------------------------------------------
// CrossSection
//------------------------------------------------------------------------------
void HistogramReader::CrossSection(TString level,
				   TString channel,
				   TString process,
				   Float_t branchingratio,
				   TString signal1_filename,
				   Float_t signal1_xs,
				   Float_t signal1_ngen,
				   TString signal2_filename,
				   Float_t signal2_xs,
				   Float_t signal2_ngen)
{
  if (_luminosity_fb < 0)
    {
      printf("\n [HistogramReader::CrossSection] Warning: reading negative luminosity\n\n");
    }


  // Get the signal (example qqWW)
  //----------------------------------------------------------------------------
  TFile* signal1_file = new TFile(_inputdir + "/" + signal1_filename + ".root");

  float signal1_counterLum = Yield((TH1D*)signal1_file->Get(level + "/h_counterLum_" + channel));
  float signal1_counterRaw = Yield((TH1D*)signal1_file->Get(level + "/h_counterRaw_" + channel));

  float counterSignal = signal1_counterLum * _luminosity_fb;

  float efficiency = signal1_counterRaw / signal1_ngen;


  // Get the second signal (example ggWW)
  //----------------------------------------------------------------------------
  if (!signal2_filename.Contains("NULL"))
    {
      TFile* signal2_file = new TFile(_inputdir + "/" + signal2_filename + ".root");

      float signal2_counterLum = Yield((TH1D*)signal2_file->Get(level + "/h_counterLum_" + channel));
      float signal2_counterRaw = Yield((TH1D*)signal2_file->Get(level + "/h_counterRaw_" + channel));

      counterSignal += (signal2_counterLum * _luminosity_fb);

      float signal1_fraction = signal1_xs / (signal1_xs + signal2_xs);
      float signal2_fraction = 1. - signal1_fraction;

      float signal1_efficiency = signal1_counterRaw / signal1_ngen;
      float signal2_efficiency = signal2_counterRaw / signal2_ngen;

      efficiency = signal1_fraction*signal1_efficiency + signal2_fraction*signal2_efficiency;
    }


  // Get the backgrounds
  //----------------------------------------------------------------------------
  float counterBackground = 0;

  for (UInt_t i=0; i<_mcfile.size(); i++) {

    if (_mclabel[i].EqualTo(process)) continue;

    _mcfile[i]->cd();

    TH1D* dummy = (TH1D*)_mcfile[i]->Get(level + "/h_counterLum_" + channel);

    float counterDummy = Yield(dummy);

    if (_luminosity_fb > 0 && _mcscale[i] > -999) counterDummy *= _luminosity_fb;

    if (_mcscale[i] > 0) counterDummy *= _mcscale[i];

    counterBackground += counterDummy;
  }


  // Get the data
  //----------------------------------------------------------------------------
  if (_datafile)
    {
      _datafile->cd();

      TH1D* dummy = (TH1D*)_datafile->Get(level + "/h_counterLum_" + channel);

      _datahist = (TH1D*)dummy->Clone();      
    }

  float counterData = Yield(_datahist);


  // Cross-section calculation
  //----------------------------------------------------------------------------  
  float xs = (counterData - counterBackground) / (1e3 * _luminosity_fb * efficiency * branchingratio);
  float mu = (counterData - counterBackground) / (counterSignal);


  // Statistical error
  //----------------------------------------------------------------------------  
  float xsErrorStat = sqrt(counterData) / (1e3 * _luminosity_fb * efficiency * branchingratio);
  float muErrorStat = sqrt(counterData) / (counterSignal); 

 
  // Print
  //----------------------------------------------------------------------------  
  printf("      channel = %s\n", channel.Data());
  printf("        ndata = %.0f\n", counterData);
  printf("         nbkg = %.2f\n", counterBackground);
  printf(" ndata - nbkg = %.2f\n", counterData - counterBackground);
  printf("      nsignal = %.2f\n", counterSignal);
  printf("           mu = (ndata - nbkg) / nsignal = %.2f +- %.2f (stat) +- %.2f (lumi)\n", mu, muErrorStat, mu * lumi_error_percent / 1e2);
  printf("         lumi = %.0f pb\n", 1e3 * _luminosity_fb);
  printf("           br = %f\n", branchingratio);
  printf("          eff = %.4f\n", efficiency);
  printf("           xs = (ndata - nbkg) / (lumi * eff * br) = %.2f +- %.2f (stat) +- %.2f (lumi) pb\n\n", xs, xsErrorStat, xs * lumi_error_percent / 1e2);
}


//-----------------------------------------------------------------------------
// DrawLatex 
//------------------------------------------------------------------------------
void HistogramReader::DrawLatex(Font_t      tfont,
				Float_t     x,
				Float_t     y,
				Float_t     tsize,
				Short_t     align,
				const char* text,
				Bool_t      setndc)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC      (setndc);
  tl->SetTextAlign( align);
  tl->SetTextFont ( tfont);
  tl->SetTextSize ( tsize);

  tl->Draw("same");
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
TLegend* HistogramReader::DrawLegend(Float_t x1,
				     Float_t y1,
				     TH1*    hist,
				     TString label,
				     TString option,
				     Bool_t  drawyield,
				     Float_t tsize,
				     Float_t xoffset,
				     Float_t yoffset)
{
  drawyield &= (_drawyield && !_publicstyle);

  TLegend* legend = new TLegend(x1,
				y1,
				x1 + xoffset,
				y1 + yoffset);
  
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  TString final_label = Form(" %s", label.Data());

  if (drawyield)
    final_label = Form("%s (%.0f)", final_label.Data(), Yield(hist));

  if (Yield(hist) < 0)
    printf("\n [HistogramReader::DrawLegend] Warning: %s %s yield = %f\n\n",
	   label.Data(),
	   hist->GetName(),
	   Yield(hist));

  legend->AddEntry(hist, final_label.Data(), option.Data());
  legend->Draw();

  WriteYields(hist, label);

  return legend;
}


//------------------------------------------------------------------------------
// GetMaximum
//------------------------------------------------------------------------------
Float_t HistogramReader::GetMaximum(TH1*    hist,
				    Float_t xmin,
				    Float_t xmax,
				    Bool_t  binError)
{
  UInt_t nbins = hist->GetNbinsX();

  TAxis* axis = (TAxis*)hist->GetXaxis();
  
  Int_t firstBin = (xmin > -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax > -999) ? axis->FindBin(xmax) : nbins;

  if (firstBin < 1)     firstBin = 1;
  if (lastBin  > nbins) lastBin  = nbins;

  Float_t hmax = 0;

  for (Int_t i=firstBin; i<=lastBin; i++) {

    Float_t binHeight = hist->GetBinContent(i);

    if (binError) {
      if (_systematics.size()>0) 
	binHeight += _ErrorGr->GetErrorY(i);
      else
	binHeight += hist->GetBinError(i);
    }

    if (binHeight > hmax) hmax = binHeight;
  }

  return hmax;
}


//------------------------------------------------------------------------------
// MoveOverflows
//
// For all histogram types: nbins, xlow, xup
//
//   bin = 0;       underflow bin
//   bin = 1;       first bin with low-edge xlow INCLUDED
//   bin = nbins;   last bin with upper-edge xup EXCLUDED
//   bin = nbins+1; overflow bin
//
//------------------------------------------------------------------------------
void HistogramReader::MoveOverflows(TH1*    hist,
				    Float_t xmin,
				    Float_t xmax)
{
  int nentries = hist->GetEntries();
  int nbins    = hist->GetNbinsX();
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();


  // Underflow
  //----------------------------------------------------------------------------
  if (xmin != -999)
    {
      Int_t   firstBin = -1;
      Float_t firstVal = 0;
      Float_t firstErr = 0;
      
      for (Int_t i=0; i<=nbins+1; i++)
	{
	  if (xaxis->GetBinLowEdge(i) < xmin)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (firstBin == -1)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      firstBin = i;
	    }
	}

      firstErr = sqrt(firstErr);
  
      hist->SetBinContent(firstBin, firstVal);
      hist->SetBinError  (firstBin, firstErr);
    }


  // Overflow
  //----------------------------------------------------------------------------
  if (xmax != -999)
    {
      Int_t   lastBin = -1;
      Float_t lastVal = 0;
      Float_t lastErr = 0;
      
      for (Int_t i=nbins+1; i>=0; i--)
	{
	  Float_t lowEdge = xaxis->GetBinLowEdge(i);
      
	  if (lowEdge >= xmax)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (lastBin == -1)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      lastBin = i;
	    }
	}

      lastErr = sqrt(lastErr);
  
      hist->SetBinContent(lastBin, lastVal);
      hist->SetBinError  (lastBin, lastErr);
    }

  hist->SetEntries(nentries);
}


//------------------------------------------------------------------------------
// SetAxis
//------------------------------------------------------------------------------
void HistogramReader::SetAxis(TH1*    hist,
			      TString xtitle,
			      TString ytitle,
			      Float_t xoffset,
			      Float_t yoffset)
{
  gPad->cd();
  gPad->Update();

  // See https://root.cern.ch/doc/master/classTAttText.html#T4
  Float_t padw = gPad->XtoPixel(gPad->GetX2());
  Float_t padh = gPad->YtoPixel(gPad->GetY1());

  Float_t size = (padw < padh) ? padw : padh;

  size = 20. / size;  // Like this label size is always 20 pixels
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();
  
  xaxis->SetTitleOffset(xoffset);
  yaxis->SetTitleOffset(yoffset);

  xaxis->SetLabelSize(size);
  yaxis->SetLabelSize(size);
  xaxis->SetTitleSize(size);
  yaxis->SetTitleSize(size);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(ytitle);

  yaxis->CenterTitle();

  gPad->GetFrame()->DrawClone();
  gPad->RedrawAxis();
}

//------------------------------------------------------------------------------
// SetLuminosity
//------------------------------------------------------------------------------
//void HistogramReader::SetLuminosity(Float_t lumi,
//				    Bool_t  postfitplots)
//{
//  _lumilegend_fb = lumi;
//  _luminosity = postfitplots ? 1. : lumi;
//}



//------------------------------------------------------------------------------
// SetHistogram
//------------------------------------------------------------------------------
void HistogramReader::SetHistogram(TH1*     hist,
				   Color_t  color,
				   Style_t  fstyle,
				   Style_t  mstyle,
				   Style_t  lstyle,
				   Width_t  lwidth,
				   Int_t    ngroup,
				   Bool_t   moveoverflow,
				   Float_t& xmin,
				   Float_t& xmax)
{
  if (!hist)
    {
      printf("\n [HistogramReader::SetHistogram] Error: histogram does not exist\n\n");
      return;
    }

  if (xmin == -999) xmin = hist->GetXaxis()->GetXmin();
  if (xmax == -999) xmax = hist->GetXaxis()->GetXmax();

  hist->SetFillColor(color );
  hist->SetFillStyle(fstyle);

  hist->SetLineColor(color );
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);

  hist->SetMarkerColor(color );
  hist->SetMarkerStyle(mstyle);

  if (_stackoption.Contains("nostack") && Yield(hist) > 0)
    {
      hist->SetFillStyle(0);
      hist->SetLineWidth(2);

      hist->Scale(1. / Yield(hist));
    }


  // Rebin and move overflow bins
  //----------------------------------------------------------------------------
  if (ngroup > 0) hist->Rebin(ngroup);
  
  if (moveoverflow) MoveOverflows(hist, xmin, xmax);
}


//------------------------------------------------------------------------------
// Yield
//------------------------------------------------------------------------------
Float_t HistogramReader::Yield(TH1* hist)
{
  if (!hist) return 0;

  Float_t hist_yield = hist->Integral(-1, -1);

  return hist_yield;
}


//------------------------------------------------------------------------------
// Error
//------------------------------------------------------------------------------
Float_t HistogramReader::Error(TH1* hist)
{
  if (!hist) return 0;

  Float_t hist_error = sqrt(hist->GetSumw2()->GetSum());

  return hist_error;
}


//------------------------------------------------------------------------------
// EventsByCut
//------------------------------------------------------------------------------
void HistogramReader::EventsByCut(TFile*  file,
				  TString analysis,
				  TString hname)
{
  // Check if the evolution histogram already exists
  TH1D* test_hist = GetHistogram(file, analysis + "/" + hname + "_evolution");//(TH1D*)file->Get(analysis + "/" + hname + "_evolution");

  if (test_hist) return;


  // Get the number of bins
  Int_t nbins = 0;
  
  for (Int_t i=0; i<ncut; i++)
    {
      if (!scut[i].Contains(analysis + "/")) continue;

      nbins++;
    }


  // Create and fill the evolution histogram
  file->cd(analysis);

  TH1D* hist = new TH1D(hname + "_evolution", "", nbins, -0.5, nbins-0.5);

  for (Int_t i=0, bin=0; i<ncut; i++)
    {
      if (!scut[i].Contains(analysis + "/")) continue;

      TH1D* dummy = GetHistogram(file, scut[i] + "/" + hname);//(TH1D*)file->Get(scut[i] + "/" + hname);

      bin++;

      hist->SetBinContent(bin, Yield(dummy));
      hist->SetBinError  (bin, Error(dummy));


      // Change the evolution histogram x-axis labels
      TString tok, icut;

      Ssiz_t from = 0;

      while (scut[i].Tokenize(tok, from, "_")) icut = tok;

      hist->GetXaxis()->SetBinLabel(bin, icut);
    }


  // Write the evolution histogram
  hist->Write();
  file->cd();
}


//------------------------------------------------------------------------------
// LoopEventsByCut
//------------------------------------------------------------------------------
void HistogramReader::LoopEventsByCut(TString analysis, TString hname)
{
  if (_datafile) EventsByCut(_datafile, analysis, hname);

  for (UInt_t i=0; i<_mcfile.size(); i++) EventsByCut(_mcfile[i], analysis, hname);

  for (UInt_t i=0; i<_signalfile.size(); i++) EventsByCut(_signalfile[i], analysis, hname);
}


//------------------------------------------------------------------------------
// EventsByChannel
//------------------------------------------------------------------------------
void HistogramReader::EventsByChannel(TFile*  file,
				      TString level)
{
  // Check if the evolution histogram already exists
  TH1D* test_hist = GetHistogram(file, level + "/h_counterLum_evolution");//(TH1D*)file->Get(level + "/h_counterLum_evolution");

  if (test_hist) return;


  // Get the number of bins
  Int_t firstchannel = (level.Contains("WZ/")) ? eee : ee;
  Int_t lastchannel  = (level.Contains("WZ/")) ? lll : ll;
  
  Int_t nbins = 0;
  
  for (Int_t i=firstchannel; i<=lastchannel; i++) nbins++;


  // Create and fill the evolution histogram
  file->cd(level);

  TH1D* hist = new TH1D("h_counterLum_evolution", "", nbins, -0.5, nbins-0.5);

  for (Int_t i=firstchannel, bin=0; i<=lastchannel; i++)
    {
      TH1D* dummy = GetHistogram(file, level + "h_counterLum" + schannel[i]);//(TH1D*)file->Get(level + "/h_counterLum_" + schannel[i]);

      bin++;

      hist->SetBinContent(bin, Yield(dummy));
      hist->SetBinError  (bin, Error(dummy));

      hist->GetXaxis()->SetBinLabel(bin, lchannel[i]);
    }


  // Write the evolution histogram
  hist->Write();
  file->cd();
}


//------------------------------------------------------------------------------
// LoopEventsByChannel
//------------------------------------------------------------------------------
void HistogramReader::LoopEventsByChannel(TString level)
{
  if (_datafile) EventsByChannel(_datafile, level);

  for (UInt_t i=0; i<_mcfile.size(); i++) EventsByChannel(_mcfile[i], level);

  for (UInt_t i=0; i<_signalfile.size(); i++) EventsByChannel(_signalfile[i], level);
}


//------------------------------------------------------------------------------
// GetBestScoreX
//------------------------------------------------------------------------------
Float_t HistogramReader::GetBestScoreX(TH1*    sig_hist,
				       TH1*    bkg_hist,
				       TString fom)
{
  Int_t nbins = sig_hist->GetNbinsX();

  Float_t score_value = 0;
  Float_t score_x     = 0;
  Float_t sig_total   = Yield(sig_hist);


  // For The Punzi Effect
  // http://arxiv.org/pdf/physics/0308063v2.pdf
  Float_t a = 5.;
  Float_t b = 1.645;  // Corresponds to a p-value equal to 0.05


  for (UInt_t k=0; k<nbins+1; k++) {

    Float_t sig_yield = sig_hist->Integral(k, nbins+1);
    Float_t bkg_yield = bkg_hist->Integral(k, nbins+1);

    Float_t sig_eff = (sig_total > 0.) ? sig_yield / sig_total : -999;

    if (sig_yield > 0. && bkg_yield > 0.)
      {
	Float_t score = -999;

	if (fom.EqualTo("S / #sqrt{B}"))   score = sig_yield / sqrt(bkg_yield);
	if (fom.EqualTo("S / #sqrt{S+B}")) score = sig_yield / sqrt(sig_yield + bkg_yield);
	if (fom.EqualTo("S / B"))          score = sig_yield / bkg_yield;
	if (fom.EqualTo("Punzi Eq.6"))     score =   sig_eff / (b*b + 2*a*sqrt(bkg_yield) + b*sqrt(b*b + 4*a*sqrt(b) + 4*bkg_yield)); 
	if (fom.EqualTo("Punzi Eq.7"))     score =   sig_eff / (a/2 + sqrt(bkg_yield));

	if (score > score_value)
	  {
	    score_value = score;
	    score_x     = sig_hist->GetBinCenter(k);
	  }
      }
  }


  printf("\n [HistogramReader::GetBestScoreX] x = %.2f (%.2f < x < %.2f) has the best %s (%f)\n\n",
  	 score_x,
  	 sig_hist->GetXaxis()->GetXmin(),
  	 sig_hist->GetXaxis()->GetXmax(),
   	 fom.Data(),
  	 score_value);


  return score_x;
}


//------------------------------------------------------------------------------
// GetBestSignalScoreX
//------------------------------------------------------------------------------
Float_t HistogramReader::GetBestSignalScoreX(TString hname,
					     TString fom,
					     Int_t   ngroup)
{
  printf("\n [HistogramReader::GetBestSignalScoreX] Warning: reading only the first signal\n");


  // Get the signals
  //----------------------------------------------------------------------------
  _signalhist.clear();

  for (UInt_t i=0; i<_signalfile.size(); i++) {

    _signalfile[i]->cd();

    TH1D* dummy = (TH1D*)_signalfile[i]->Get(hname);

    _signalhist.push_back((TH1D*)dummy->Clone());

    if (_luminosity_fb > 0) _signalhist[i]->Scale(_luminosity_fb);

    if (ngroup > 0) _signalhist[i]->Rebin(ngroup);
  }

  
  // Get the backgrounds
  //----------------------------------------------------------------------------
  _mchist.clear();

  THStack* mcstack = new THStack(hname + "_mcstack", hname + "_mcstack");

  for (UInt_t i=0; i<_mcfile.size(); i++) {

    _mcfile[i]->cd();

    TH1D* dummy = (TH1D*)_mcfile[i]->Get(hname);

    _mchist.push_back((TH1D*)dummy->Clone());

    if (_luminosity_fb > 0 && _mcscale[i] > -999) _mchist[i]->Scale(_luminosity_fb);

    if (_mcscale[i] > 0) _mchist[i]->Scale(_mcscale[i]);

    if (ngroup > 0) _mchist[i]->Rebin(ngroup);

    mcstack->Add(_mchist[i]);
  }


  // Get the best score
  //----------------------------------------------------------------------------
  TH1D* backgroundhist = (TH1D*)(mcstack->GetStack()->Last());

  return GetBestScoreX(_signalhist[0], backgroundhist, fom);
}


//------------------------------------------------------------------------------
// WriteYields
//------------------------------------------------------------------------------
void HistogramReader::WriteYields(TH1*    hist,
				  TString label)
{
  TString hname = hist->GetName();

  if (!_writeyields) return;

  if (_writelabels)
    {
      _writelabels = false;

      _yields_table << Form("\n %14s", " ");
        
      for (int i=1; i<=hist->GetNbinsX(); i++) {

	TString binlabel = (TString)hist->GetXaxis()->GetBinLabel(i);
	    
	_yields_table << Form(" | %-32s", binlabel.Data());
      }

      _yields_table << Form("\n");
    }

  _yields_table << Form(" %14s", label.Data());

  for (int i=1; i<=hist->GetNbinsX(); i++) {

    float process_yield = hist->GetBinContent(i);
    float process_error = sqrt(hist->GetSumw2()->At(i));

    if (label.EqualTo("data"))
      {
	_yields_table << Form(" | %8.0f %14s", process_yield, " ");
      }
    else
      {
	_yields_table << Form(" | %11.2f +/- %7.2f", process_yield, process_error);
      }

    int denominator = (hname.Contains("counterLum_evolution")) ? hist->GetNbinsX() : 1;

    float process_percent = 1e2 * process_yield / hist->GetBinContent(denominator);

    _yields_table << Form(" (%5.1f%s)", process_percent, "%");
  }

  _yields_table << Form("\n");
}


//------------------------------------------------------------------------------
// Roc
//------------------------------------------------------------------------------
void HistogramReader::Roc(TString hname,
			  TString xtitle,
			  Int_t   npoints,
			  TString units,
			  Float_t xmin,
			  Float_t xmax,
			  TString fom)
{
  // Get the signal
  //----------------------------------------------------------------------------
  THStack* stack_sig = new THStack(hname + "_stack_sig", hname + "_stack_sig");

  for (int i=0; i<_roc_signalfile.size(); ++i)
    {
      _roc_signalfile[i]->cd();

      TH1D* dummy = (TH1D*)(_roc_signalfile[i]->Get(hname))->Clone();

      if (_luminosity_fb > 0 && _roc_signalscale[i] > -999) dummy->Scale(_luminosity_fb);
      
      stack_sig->Add(dummy);
    }

  TH1D* hSig = (TH1D*)(stack_sig->GetStack()->Last());


  // Get the backgrounds
  //----------------------------------------------------------------------------
  THStack* stack_bkg = new THStack(hname + "_stack_bkg", hname + "_stack_bkg");

  for (int j=0; j<_roc_backgroundfile.size(); ++j)
    {
      _roc_backgroundfile[j]->cd();

      TH1D* dummy = (TH1D*)(_roc_backgroundfile[j]->Get(hname))->Clone();

      if (_luminosity_fb > 0 && _roc_backgroundscale[j] > -999) dummy->Scale(_luminosity_fb);

      stack_bkg->Add(dummy);
    }

  TH1D* hBkg = (TH1D*)(stack_bkg->GetStack()->Last());


  // For The Punzi Effect
  // http://arxiv.org/pdf/physics/0308063v2.pdf
  Float_t a = 5.;
  Float_t b = 1.645;  // Corresponds to a p-value equal to 0.05


  // Compute ROC and significance
  //----------------------------------------------------------------------------
  float step = (xmax - xmin) / npoints;

  TGraph* rocGraph_min = new TGraph();
  TGraph* rocGraph_max = new TGraph();
  TGraph* sigGraph_min = new TGraph();
  TGraph* sigGraph_max = new TGraph();

  Float_t score_value_min = 0;
  Float_t score_value_max = 0;
  Float_t score_x_min     = 0;
  Float_t score_x_max     = 0;

  Float_t sigEff_score_x_min = -999;
  Float_t bkgEff_score_x_min = -999;
  Float_t sigEff_score_x_max = -999;
  Float_t bkgEff_score_x_max = -999;

  Float_t sigTotal = hSig->Integral(-1, -1);
  Float_t bkgTotal = hBkg->Integral(-1, -1);

  for (int s=0; s<=npoints; ++s) {

    Float_t sigYield_min = 0;
    Float_t sigYield_max = 0;
    Float_t bkgYield_min = 0;
    Float_t bkgYield_max = 0;

    sigYield_max += hSig->Integral(-1, hSig->FindBin(xmin + s*step));
    bkgYield_max += hBkg->Integral(-1, hBkg->FindBin(xmin + s*step));

    sigYield_min += hSig->Integral(hSig->FindBin(xmin + s*step), -1);
    bkgYield_min += hBkg->Integral(hBkg->FindBin(xmin + s*step), -1);

    Float_t sigEff_max = (sigTotal != 0) ? sigYield_max / sigTotal : -999;
    Float_t bkgEff_max = (bkgTotal != 0) ? bkgYield_max / bkgTotal : -999;

    Float_t sigEff_min = (sigTotal != 0) ? sigYield_min / sigTotal : -999;
    Float_t bkgEff_min = (bkgTotal != 0) ? bkgYield_min / bkgTotal : -999;

    Float_t score_min = -999;

    if (sigYield_min > 0. && bkgYield_min > 0.)
      {
        if (fom.EqualTo("S / #sqrt{B}"))   score_min = sigYield_min / sqrt(bkgYield_min);  
        if (fom.EqualTo("S / #sqrt{S+B}")) score_min = sigYield_min / sqrt(bkgYield_min + sigYield_min);
        if (fom.EqualTo("S / B"))          score_min = sigYield_min / bkgYield_min;
        if (fom.EqualTo("Punzi Eq.6"))     score_min =   sigEff_min / (b*b + 2*a*sqrt(bkgYield_min) + b*sqrt(b*b + 4*a*sqrt(b) + 4*bkgYield_min));
        if (fom.EqualTo("Punzi Eq.7"))     score_min =   sigEff_min / (a/2 + sqrt(bkgYield_min));
      }

    Float_t score_max = -999;

    if (sigYield_max > 0. && bkgYield_max > 0.)
      {
        if (fom.EqualTo("S / #sqrt{B}"))   score_max = sigYield_max / sqrt(bkgYield_max);
        if (fom.EqualTo("S / #sqrt{S+B}")) score_max = sigYield_max / sqrt(bkgYield_max + sigYield_max);
        if (fom.EqualTo("S / B"))          score_max = sigYield_max / bkgYield_max;
        if (fom.EqualTo("Punzi Eq.6"))     score_max =   sigEff_max / (b*b + 2*a*sqrt(bkgYield_max) + b*sqrt(b*b + 4*a*sqrt(b) + 4*bkgYield_max));
        if (fom.EqualTo("Punzi Eq.7"))     score_max =   sigEff_max / (a/2 + sqrt(bkgYield_max));
      }

    if (score_min > score_value_min) {
      score_value_min    = score_min;
      score_x_min        = xmin + s*step;
      sigEff_score_x_min = sigEff_min;
      bkgEff_score_x_min = bkgEff_min;
    }

    if (score_max > score_value_max) {
      score_value_max    = score_max;
      score_x_max        = xmin + s*step;
      sigEff_score_x_max = sigEff_max;
      bkgEff_score_x_max = bkgEff_max;
    }

    rocGraph_min->SetPoint(s, sigEff_min, 1 - bkgEff_min);
    rocGraph_max->SetPoint(s, sigEff_max, 1 - bkgEff_max);

    sigGraph_min->SetPoint(s, xmin + s*step, score_min);
    sigGraph_max->SetPoint(s, xmin + s*step, score_max);
  }


  printf("\n");
  printf(" [HistogramReader::Roc] Reading %s from %.2f to %.2f\n\n", hname.Data(), xmin, xmax);
  printf(" The best %s (%f) corresponds to x > %7.2f %s (S_eff = %6.2f\%, B_eff = %6.2f\%)\n",
	 fom.Data(),
	 score_value_min,
	 score_x_min,
	 units.Data(),
	 1e2 * sigEff_score_x_min,
	 1e2 * bkgEff_score_x_min);
  printf(" The best %s (%f) corresponds to x < %7.2f %s (S_eff = %6.2f\%, B_eff = %6.2f\%)\n",
	 fom.Data(),
	 score_value_max,
	 score_x_max,
	 units.Data(),
	 1e2 * sigEff_score_x_max,
	 1e2 * bkgEff_score_x_max);
  printf("\n");
  

  // Draw and save ROC
  //----------------------------------------------------------------------------
  Color_t color_min = kRed+1;
  Color_t color_max = kBlack;

  Style_t style_min = kFullCircle;
  Style_t style_max = kOpenCircle;

  TCanvas* rocCanvas = new TCanvas(hname + " ROC", hname + " ROC");

  rocGraph_min->SetMarkerColor(color_min);
  rocGraph_min->SetMarkerStyle(style_min);
  rocGraph_min->SetMarkerSize(0.5);

  rocGraph_max->SetMarkerColor(color_max);
  rocGraph_max->SetMarkerStyle(style_max);
  rocGraph_max->SetMarkerSize(0.5);

  rocGraph_min->Draw("ap");
  rocGraph_max->Draw("psame");

  rocGraph_min->GetXaxis()->SetRangeUser(0, 1);
  rocGraph_min->GetYaxis()->SetRangeUser(0, 1);

  DrawLatex(42, 0.190, 0.945, 0.050, 11, _title);

  SetAxis(rocGraph_min->GetHistogram(), xtitle + " signal efficiency", xtitle + " background rejection", 1.5, 1.8);

  if (_savepdf) rocCanvas->SaveAs(_outputdir + hname + "_ROC.pdf");
  if (_savepng) rocCanvas->SaveAs(_outputdir + hname + "_ROC.png");


  // Draw and save significance
  //----------------------------------------------------------------------------
  TCanvas *sigCanvas = new TCanvas(hname + " significance", hname + " significance");

  TString myxtitle = (!units.Contains("NULL")) ? xtitle + " [" + units + "]" : xtitle;

  sigGraph_min->SetMarkerColor(color_min);
  sigGraph_min->SetMarkerStyle(style_min);
  sigGraph_min->SetMarkerSize(0.5);

  sigGraph_max->SetMarkerColor(color_max);
  sigGraph_max->SetMarkerStyle(style_max);
  sigGraph_max->SetMarkerSize(0.5);

  sigGraph_min->Draw("ap");
  sigGraph_max->Draw("psame");

  Float_t ymax = (score_value_min > score_value_max) ? score_value_min : score_value_max;

  ymax *= 1.5;

  sigGraph_min->GetXaxis()->SetRangeUser(xmin, xmax);
  sigGraph_min->GetYaxis()->SetRangeUser(   0, ymax);

  DrawLatex(42, 0.190, 0.945, 0.050, 11, _title);

  TH1F* dummy_min = new TH1F("dummy_min", "", 1, 0, 1);
  TH1F* dummy_max = new TH1F("dummy_max", "", 1, 0, 1);

  dummy_min->SetLineColor  (color_min);
  dummy_min->SetMarkerColor(color_min);
  dummy_min->SetMarkerStyle(style_min);

  dummy_max->SetLineColor  (color_max);
  dummy_max->SetMarkerColor(color_max);
  dummy_max->SetMarkerStyle(style_max);

  DrawLegend(0.22, 0.84, dummy_min, Form("%s > x", xtitle.Data()), "lp", false);
  DrawLegend(0.22, 0.77, dummy_max, Form("%s < x", xtitle.Data()), "lp", false);

  SetAxis(sigGraph_min->GetHistogram(), myxtitle, fom, 1.5, 2.1);

  if (_savepdf) sigCanvas->SaveAs(_outputdir + hname + "_significance.pdf");
  if (_savepng) sigCanvas->SaveAs(_outputdir + hname + "_significance.png");

  dummy_min->Delete();
  dummy_max->Delete();
}


//------------------------------------------------------------------------------
// IncludeSystematics
//------------------------------------------------------------------------------
/*void HistogramReader::IncludeSystematics(TString hname)
{
  int nbins = _mchist[0]->GetNbinsX();


  // Loop over all processes
  //----------------------------------------------------------------------------
  for (int i=0; i<_mchist.size(); i++) {

    TH1D* myhisto = (TH1D*)_mchist[0]->Clone("myhisto");

    float suma[nbins+1]; 

    for (int k=0; k<=nbins; k++) suma[k] = 0;

    TFile* myfile0 = new TFile(_inputdir + "/" + _mcfilename.at(i) + ".root", "read");

    TH1D* dummy0 = (TH1D*)myfile0->Get(hname);


    // Loop over all systematics
    //--------------------------------------------------------------------------
    for (int j=0; j<_systematics.size(); j++) {

      TFile* myfile = new TFile(_inputdir + "/" + _mcfilename.at(i) + "_" + _systematics.at(j) + ".root", "read");

      TH1D* dummy = (TH1D*)myfile->Get(hname);


      // Loop over all bins
      //------------------------------------------------------------------------
      for (int k=0; k<=nbins; k++) {

	float diff = dummy->GetBinContent(k) - dummy0->GetBinContent(k);
	
	if (_mclabel[i] == "non-prompt") diff = 0; 

	suma[k] += diff*diff;
      }

      myfile->Close();
    }

    
    // Save the sum of systematic uncertainties per bin
    //--------------------------------------------------------------------------
    for (int k=0; k<=nbins; k++) { 
	
      myhisto->SetBinContent(k, sqrt(suma[k]));

      _mchist_syst.push_back(myhisto);
    }
  }
}*/ 

TH1D* SumSRHistograms(TFile*  file, TString HistogramName) {

  TH1D *SumHisto;
  if ((!HistogramName.Contains("_Veto") && !HistogramName.Contains("_SRs_"))) {
     SumHisto = (TH1D*) file->Get(HistogramName);
  } else if (HistogramName.Contains("_Veto") && !HistogramName.Contains("_SRs_")) {
    HistogramName.ReplaceAll("_Veto", "_NoTag");
    SumHisto = (TH1D*) file->Get(HistogramName);
    SumHisto->SetDirectory(0);
    HistogramName.ReplaceAll("_NoTag", "_NoJet");
    TH1D* SumHisto2 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto2);
  } else if (!HistogramName.Contains("_Veto") && HistogramName.Contains("_SRs_")) {
    HistogramName.ReplaceAll("_SRs", "_SR1");
    SumHisto = (TH1D*) file->Get(HistogramName);
    SumHisto->SetDirectory(0);
    HistogramName.ReplaceAll("_SR1", "_SR2");
    TH1D* SumHisto2 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto2);
    HistogramName.ReplaceAll("_SR2", "_SR3");
    TH1D* SumHisto3 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto3);
  } else if (HistogramName.Contains("_Veto") && HistogramName.Contains("_SRs_")) {
    HistogramName.ReplaceAll("_Veto", "_NoTag");
    HistogramName.ReplaceAll("_SRs", "_SR1");
    SumHisto = (TH1D*) file->Get(HistogramName);
    SumHisto->SetDirectory(0);
    HistogramName.ReplaceAll("_SR1", "_SR2");
    TH1D* SumHisto2 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto2);
    HistogramName.ReplaceAll("_SR2", "_SR3");
    TH1D* SumHisto3 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto3);
    HistogramName.ReplaceAll("_NoTag", "_NoJet");
    HistogramName.ReplaceAll("_SR3", "_SR1");
    TH1D* SumHisto4 = (TH1D*) file->Get(HistogramName);
    if (!SumHisto4) return SumHisto;
    SumHisto->Add(SumHisto4);
    HistogramName.ReplaceAll("_SR1", "_SR2");
    TH1D* SumHisto5 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto5);
    HistogramName.ReplaceAll("_SR2", "_SR3");
    TH1D* SumHisto6 = (TH1D*) file->Get(HistogramName);
    SumHisto->Add(SumHisto6);
  } 
  return SumHisto;

}

TH1D* HistogramReader::GetHistogram(TFile*  file, TString HistogramName) {
  
  TH1D* ThisHisto;
  if (!HistogramName.Contains("_sf") || _inputdir.Contains("Postfit")) {
    ThisHisto = SumSRHistograms(file, HistogramName);//(TH1D*) file->Get(HistogramName);
  } else {
    HistogramName.ReplaceAll("_sf", "_ee");
    ThisHisto =  SumSRHistograms(file, HistogramName);//(TH1D*) file->Get(HistogramName);
    HistogramName.ReplaceAll("_ee", "_mm");
    TH1D* ThisHisto2 =  SumSRHistograms(file, HistogramName);//(TH1D*) file->Get(HistogramName);
    ThisHisto->Add(ThisHisto2);
  }
  return ThisHisto;

}

void FormatTableYields(float *YY, float *EY) {

  if (*EY>=10.) {
    *EY = round(*EY);
    *YY = round(*YY);
  } else if (*EY>=1.) {
    *EY = round((*EY)*10)/10.;
    *YY = round((*YY)*10)/10.;
  } else {
    *EY = round((*EY)*100)/100.;
    *YY = round((*YY)*100)/100.;
  } 

  //if (*YY<0.) *YY = *EY;

}

void HistogramReader::IncludeSystematics(TString hname)
{
  bool _verbose = false, _dotable = true;

  float StatZero = 1.84102;

  int nsystematics = _systematics.size(); 
  int nbins        = _mchist[0]->GetNbinsX();
  int nprocess     = _mchist.size();
  int nsignals     = _signalfilename.size();
   
  // Table variables
  float yieldTab      [nprocess][nbins+1];
  float errBackTab_do [nprocess][nbins+1]; 
  float errBackTab_up [nprocess][nbins+1];  
  bool FirstSystematic[nprocess];
  bool _isPostfit = false;
  for (int k=0; k <nprocess; k++){
    FirstSystematic[k] = true;
    for (int j=1; j<=nbins; j++){ 
      yieldTab[k][j] = 0;
      errBackTab_do[k][j] = 0; 
      errBackTab_up[k][j] = 0;
    }
  }
  // TGraphAssymetriErrors variables
  float errSystDo     [nsystematics][nbins+1];
  float errSystUp     [nsystematics][nbins+1]; 
  for (int i=0; i < nsystematics; i++) {
    for (int j=1; j<=nbins; j++) {
      errSystDo [i][j] = 0;      
      errSystUp [i][j] = 0;      
     }
   }
  // Signal uncertainties
  float yieldSign [nsignals][nbins+1];
  float errSignUp [nsignals][nbins+1];
  float errSignDo [nsignals][nbins+1];
  for (int j=1; j<=nbins; j++) {
    for (int s=0; s<nsignals; s++) {
      errSignUp[s][j] = 0.;
      errSignDo[s][j] = 0.;
    }
  }

   // Loop over all systematics and processes
   //----------------------------------------------------------------------------
   for (int isyst=0; isyst<nsystematics; isyst++) 
    {

      if (_verbose) {
	// Print systematic name 
	printf( "                                                         \n");  
	printf( "                                                         \n");  
	printf( "systematic name %s \n", _systematics.at(isyst).Data() );
	printf( "----------------------------------\n"); 
      }

     for (int kproce=0; kproce<nprocess; kproce++) 
       {
       if (_analysis == "Stop" && _systematics.at(isyst) == "Toppt" && _mcfilename.at(kproce) != "04_TTTo2L2Nu") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "Fake" && _mcfilename.at(kproce) != "04_TTTo2L2Nu") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "PDF"   && _mcfilename.at(kproce)=="05_ST") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "Q2"    && _mcfilename.at(kproce)=="05_ST") continue;
       if (_analysis == "Stop" && (_systematics.at(isyst)=="BtagFS" || _systematics.at(isyst)=="Fastsim" || 
				   _systematics.at(isyst)=="Pileup" || _systematics.at(isyst)=="Metfastsim" ||
				   _systematics.at(isyst)=="Isrnjet")) continue;
       if (_systematics.at(isyst).Contains("MT2ll") && !hname.Contains("h_MT2ll")) continue;
       if (_systematics.at(isyst)=="MT2llTop" && _mcfilename.at(kproce)!="04_TTTo2L2Nu" && _mcfilename.at(kproce)!="05_ST") continue;
       if (_systematics.at(isyst)=="MT2llWW" && _mcfilename.at(kproce)!="06_WW") continue;
       if (_systematics.at(isyst)=="ttZSF" && _mcfilename.at(kproce)!="10_TTZ") continue;
       if (_systematics.at(isyst)=="ZZSF" && _mcfilename.at(kproce)!="03_VZ") continue;
       if (_systematics.at(isyst)=="DYSF" && !_mcfilename.at(kproce).Contains("07_ZJets")) continue;
       if (_systematics.at(isyst)=="ZZnojet" && _mcfilename.at(kproce)!="03_VZ") continue;
       if (_systematics.at(isyst)=="DYnojet" && (!_mcfilename.at(kproce).Contains("07_ZJets") || !hname.Contains("NoJet"))) continue;
       if (_systematics.at(isyst)=="DYshape" && (!_mcfilename.at(kproce).Contains("07_ZJets") || hname.Contains("NoJet"))) continue;
       if (_systematics.at(isyst)=="ZMETjet" && !_mcfilename.at(kproce).Contains("07_ZJets") && _mcfilename.at(kproce)!="03_VZ") continue;
       if (_systematics.at(isyst)=="normWZ" && _mcfilename.at(kproce)!="02_WZTo3LNu") continue;
       if (_systematics.at(isyst)=="normWW" && _mcfilename.at(kproce)!="06_WW") continue;
       if (_systematics.at(isyst)=="normTtbar" && _mcfilename.at(kproce)!="04_TTTo2L2Nu") continue;
       if (_systematics.at(isyst)=="normTW" && _mcfilename.at(kproce)!="05_ST") continue;
       if (_systematics.at(isyst)=="normTTW" && _mcfilename.at(kproce)!="09_TTW") continue;
       if (_systematics.at(isyst)=="normHWW" && _mcfilename.at(kproce)!="11_HWW") continue;
       if (_systematics.at(isyst)=="normVVV" && _mcfilename.at(kproce)!="13_VVV") continue;
       
       TFile* myfile0 = myfile0 = new TFile(_inputdir + "/" + _mcfilename.at(kproce) + ".root", "read");;
       
       TH1D* dummy0 = GetHistogram(myfile0, hname);//(TH1D*)myfile0->Get( hname );//nominal
       if (_luminosity_fb > 0 && _mcscale[kproce] > -999) dummy0->Scale(_luminosity_fb);		
       if (_mcscale[kproce] > 0) dummy0->Scale(_mcscale[kproce]);

       if (FirstSystematic[kproce]) {
	 for (int ibin = 1; ibin<=nbins; ibin++) {
	   yieldTab     [kproce][ibin] = dummy0->GetBinContent(ibin); 
	   if ( dummy0->GetBinContent(ibin)<0.) yieldTab [kproce][ibin] = dummy0->GetBinError(ibin);
	 }
	 FirstSystematic[kproce] = false;
       }

       if (_systematics.at(isyst) == "Postfit") {
	 _isPostfit = true;
	 myfile0->Close();
	 continue; 
       }

       if (_systematics.at(isyst) == "Statistics") {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float StatUncert2 = dummy0->GetSumw2()->At(ibin);
	   if (StatUncert2<0.0001) {
	     bool ApplyZeroStat = true;
	     if (_mcfilename.at(kproce).Contains("07_ZJets")) {
	       ApplyZeroStat = false;
	       if (ibin>1)
		 if (fabs(dummy0->GetBinContent(ibin-1)>0.005)) ApplyZeroStat = true;
	       if (ibin<nbins)
		 if (fabs(dummy0->GetBinContent(ibin+1)>0.005)) ApplyZeroStat = true;
	     }
	     if (ApplyZeroStat) 
	       StatUncert2 = TMath::Power(StatZero*dummy0->Integral()/dummy0->GetEntries(), 2);
	   }
	   errBackTab_up[kproce][ibin] += StatUncert2;
	   errBackTab_do[kproce][ibin] += StatUncert2;
	   errSystUp [isyst][ibin] += StatUncert2;
	   errSystDo [isyst][ibin] += StatUncert2;
	 }
	 myfile0->Close();
	 continue;
       }

       if (_systematics.at(isyst) == "Luminosity") {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float LumiBinError = yieldTab[kproce][ibin]*lumi_error_percent/1e2;
	   errBackTab_up [kproce][ibin] += LumiBinError*LumiBinError;
	   errBackTab_do [kproce][ibin] += LumiBinError*LumiBinError;
	   errSystUp [isyst][ibin] += LumiBinError;
	   errSystDo [isyst][ibin] += -LumiBinError;
	 }
	 myfile0->Close();
	 continue;
       }
       
       if (_systematics.at(isyst) == "Trigger") {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float TrigBinError = yieldTab[kproce][ibin]*2./1e2;
	   errBackTab_up [kproce][ibin] += TrigBinError*TrigBinError;
	   errBackTab_do [kproce][ibin] += TrigBinError*TrigBinError;
	   errSystUp [isyst][ibin] += TrigBinError;
	   errSystDo [isyst][ibin] += -TrigBinError;
	 }
	 myfile0->Close();
	 continue;
       }
       
       if (_systematics.at(isyst).Contains("MT2ll")) {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float relError = 0.;
	   if (ibin==4) relError = 0.05;
	   else if (ibin==5) relError = 0.10;
	   else if (ibin==6) relError = 0.20;
	   else if (ibin==7) relError = 0.30;
	   float absError = relError*yieldTab[kproce][ibin]; 
	   errBackTab_up [kproce][ibin] += absError*absError;
	   errBackTab_do [kproce][ibin] += absError*absError;
	   errSystUp [isyst][ibin] += absError;
	   errSystDo [isyst][ibin] += -absError;
	 }
	 myfile0->Close();
	 continue;
       }

       if (_systematics.at(isyst)=="ttZSF" || _systematics.at(isyst)=="ZZSF" || 
	   _systematics.at(isyst)=="DYSF" || _systematics.at(isyst)=="DYnojet") { 
	 if (_mcscale[kproce] > 0) {
	   float errSF;
	   if (_systematics.at(isyst)=="ttZSF") errSF = 0.36;
	   else if (_systematics.at(isyst)=="ZZSF") {
	     if (hname.Contains("NoJet")) errSF = 0.14;
	     else if (hname.Contains("Veto")) errSF = 0.17;
	     else errSF = 0.20;
	   } else if (_systematics.at(isyst)=="DYSF") errSF = 0.64;
	   else if (_systematics.at(isyst)=="DYnojet") errSF = 4.;
	   for (int ibin=1; ibin<=nbins; ibin++) {
	     float relErr = errSF;
	     if (_mcscale[kproce]>0.) relErr /= _mcscale[kproce];
	     float errNorm = dummy0->GetBinContent(ibin)*relErr; 
	     errBackTab_up [kproce][ibin] += errNorm*errNorm;
	     errBackTab_do [kproce][ibin] += errNorm*errNorm;
	     errSystUp [isyst][ibin] += errNorm;
	     errSystDo [isyst][ibin] += -errNorm;
	   }
	 }
	 myfile0->Close();
	 continue;
       }

       if (_systematics.at(isyst)=="normWZ" || _systematics.at(isyst)=="normWW" || _systematics.at(isyst)=="normTtbar" || 
	   _systematics.at(isyst)=="normTW" || _systematics.at(isyst)=="normTTW" || _systematics.at(isyst)=="normHWW" || 
	   _systematics.at(isyst)=="normVVV") {
	 float relErr = 0.5;
	 if (_systematics.at(isyst)=="normWZ") relErr = 0.04;
	 if (_systematics.at(isyst)=="normWW" || _systematics.at(isyst)=="normTtbar" || 
	     _systematics.at(isyst)=="normTW") relErr = 0.1;
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float errNorm = dummy0->GetBinContent(ibin)*relErr; 
	   errBackTab_up [kproce][ibin] += errNorm*errNorm;
	   errBackTab_do [kproce][ibin] += errNorm*errNorm;
	   errSystUp [isyst][ibin] += errNorm;
	   errSystDo [isyst][ibin] += -errNorm;
	 }
	 myfile0->Close();
	 continue;
       }
       
       TFile* myfile1;
       TFile* myfile2;
       
       if (_minitreebased) // This is for ttdm
	{
         if (isyst%2 != 0) continue; //Only takes the first one systematic_up(down).root as reference
	 myfile1 = new TFile(_inputdir + "/" + _mcfilename.at(kproce) + _systematics.at(isyst)   + ".root", "read");
	 myfile2 = new TFile(_inputdir + "/" + _mcfilename.at(kproce) + _systematics.at(isyst+1) + ".root", "read");
	}
       else
	 {
	   myfile1 = new TFile(_inputdir + "/../../" + _systematics.at(isyst) + "up/" + _analysis + "/" + _mcfilename.at(kproce)   + ".root", "read");//up
	   myfile2 = new TFile(_inputdir + "/../../" + _systematics.at(isyst) + "do/" + _analysis + "/" + _mcfilename.at(kproce)   + ".root", "read");//down
        }

       TH1D* dummy1 = GetHistogram(myfile1, hname);//(TH1D*)myfile1->Get( hname );//up
       TH1D* dummy2 = GetHistogram(myfile2, hname);//(TH1D*)myfile2->Get( hname );//down      

       if (_luminosity_fb > 0 && _mcscale[kproce] > -999)
        {
          dummy1->Scale(_luminosity_fb);
          dummy2->Scale(_luminosity_fb);
        }
		
       if (_mcscale[kproce] > 0)
	{
          dummy1->Scale(_mcscale[kproce]);
          dummy2->Scale(_mcscale[kproce]);
        }
       
       if (_verbose) {
	 // Print Process name       
	 printf( "                                                         \n");  
	 printf("process name %s\n", _mcfilename.at(kproce).Data()); 
	 printf( "                                                         \n");  
       }

       // Loop over all bins (Underflow is not included)
       //--------------------------------------------------------------------  
       for (int ibin=1; ibin<=nbins; ibin++) 
        {
	 if (dummy0->GetBinContent(ibin)<0.) dummy0->SetBinContent(ibin, 0.001);
	 if (dummy1->GetBinContent(ibin)<0.) dummy1->SetBinContent(ibin, 0.001);
	 if (dummy2->GetBinContent(ibin)<0.) dummy2->SetBinContent(ibin, 0.001);
         errSystUp [isyst][ibin] += (dummy1->GetBinContent(ibin) - dummy0->GetBinContent(ibin));
         errSystDo [isyst][ibin] += (dummy2->GetBinContent(ibin) - dummy0->GetBinContent(ibin)); 

	 float ErrUp = 0., ErrDo = 0.; 
	 float VarUp = (dummy1->GetBinContent(ibin) - dummy0->GetBinContent(ibin));
	 float VarDo = (dummy2->GetBinContent(ibin) - dummy0->GetBinContent(ibin)); 
	 if (VarDo<=0. && VarUp>=0.) {
	   ErrUp = VarUp;
	   ErrDo = VarDo;
	 } else if (VarDo>=0. && VarUp<=0.) {
	   ErrUp = VarDo;
	   ErrDo = VarUp;
	 } else if (VarDo>0.) {
	   ErrUp = (VarUp>VarDo) ? VarUp : VarDo;
	 } else if (VarDo<0.) {
	   ErrDo = (VarUp<VarDo) ? VarUp : VarDo;
	 }
         errBackTab_up[kproce][ibin] += ErrUp*ErrUp;
         errBackTab_do[kproce][ibin] += ErrDo*ErrDo;

	 if (_verbose) { 
	   // Print Bin Information per process       
	   printf( "Print Bin Information per process \n");
	   printf( "                                                         \n");  
	   printf( " bin number  = %i\n", ibin );  
	   printf( "                                                         \n");  
	   printf( " nominal     = %f\n", dummy0->GetBinContent(ibin));             
	   printf( " SF+ errUp   = %f\n", dummy1->GetBinContent(ibin));             
	   printf( " SF+ errDo   = %f\n", dummy2->GetBinContent(ibin));             
	   // Print Systematic Error per process
	   printf( "Print Systematic Error per process \n");
	   printf( "                                                         \n");  
	   printf( " ErrsystUp  = %.5f\n", (dummy1->GetBinContent(ibin) - dummy0->GetBinContent(ibin)) );
	   printf( " ErrsystDo  = %.5f\n", (dummy2->GetBinContent(ibin) - dummy0->GetBinContent(ibin)) );       
	   printf( " rel_ErrsystUp  = %.5f  %\n", (dummy1->GetBinContent(ibin) - dummy0->GetBinContent(ibin)) / dummy0->GetBinContent(ibin) );
	   printf( " rel_ErrsystDo  = %.5f  %\n", (dummy2->GetBinContent(ibin) - dummy0->GetBinContent(ibin)) / dummy0->GetBinContent(ibin) );       
	 }
         // 

        } 

       myfile0->Close();
       myfile1->Close();
       myfile2->Close();
     
      }
    }

   // Loop over signals
   //----------------------------------------------------------------------------
   bool _doMetFastSim = false;
   TString hnamegen = hname;
   if (hname.Contains("_SR") && (hname.Contains("h_MT2ll_") || hname.Contains("h_MT2llisr_"))) {
     for (int isyst=0; isyst<nsystematics; isyst++) {
       if (_systematics.at(isyst)=="Metfastsim")  {
	 
	 hnamegen.ReplaceAll("_SR1", "_SR1gen");
	 hnamegen.ReplaceAll("_SR2", "_SR2gen");
	 hnamegen.ReplaceAll("_SR3", "_SR3gen");
	 hnamegen.ReplaceAll("h_MT2ll_", "h_MT2llgen_");
	 hnamegen.ReplaceAll("h_MT2llisr_", "h_MT2llisrgen_");
	 _doMetFastSim = true;
	 
       }
     }
   }

   for (int kproce=0; kproce<nsignals; kproce++) {

     TFile* myfile0 = TFile::Open(_inputdir + "/" + _signalfilename.at(kproce) + ".root");

     TH1D* dummy0 = GetHistogram(myfile0, hname);   //(TH1D*)myfile0->Get( hname );
     TH1D* dummy3 = GetHistogram(myfile0, hnamegen);//(TH1D*)myfile0->Get( hnamegen );
     if (_luminosity_fb > 0) {
       dummy0->Scale(_luminosity_fb); 
       if (_doMetFastSim) dummy3->Scale(_luminosity_fb); 
     }

     for (int ibin = 1; ibin<=nbins; ibin++) {
       if (dummy0->GetBinContent(ibin)<0.) dummy0->SetBinContent(ibin, 0.001);
       if (dummy3->GetBinContent(ibin)<0.) dummy3->SetBinContent(ibin, 0.001);
       if (!_doMetFastSim) {
	 yieldSign [kproce][ibin] = dummy0->GetBinContent(ibin);
       } else {
	 yieldSign [kproce][ibin] = (dummy0->GetBinContent(ibin) + dummy3->GetBinContent(ibin))/2.;
       }
     }

     for (int isyst=0; isyst<nsystematics; isyst++) {

       if (_systematics.at(isyst) == "Postfit") continue; 
       if (_analysis == "Stop" && _systematics.at(isyst) == "Toppt") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "Fake") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "PDF") continue;
       if (_analysis == "Stop" && _systematics.at(isyst) == "Q2") continue;
       if (_systematics.at(isyst).Contains("MT2ll")) continue;
       if (_systematics.at(isyst)=="ttZSF" || _systematics.at(isyst)=="ZZSF" || _systematics.at(isyst)=="DYSF") continue;
       if (_systematics.at(isyst)=="ZZnojet" || _systematics.at(isyst)=="DYnojet" || _systematics.at(isyst)=="ZMETjet" || _systematics.at(isyst)=="DYshape") continue;
       if (_systematics.at(isyst)=="normWZ" || _systematics.at(isyst)=="normWW" || _systematics.at(isyst)=="normTtbar" || 
	   _systematics.at(isyst)=="normTW" || _systematics.at(isyst)=="normTTW" || _systematics.at(isyst)=="normHWW" || 
	   _systematics.at(isyst)=="normVVV") continue;

       if ( _systematics.at(isyst) == "Statistics") {
	 if (!_doMetFastSim) {
	   for (int ibin = 1; ibin<=nbins; ibin++) {
	     float StatUncert2 = dummy0->GetSumw2()->At(ibin);
	     if (StatUncert2<0.0001) StatUncert2 = TMath::Power(StatZero*dummy0->Integral()/dummy0->GetEntries(), 2);
	     errSignUp [kproce][ibin] += StatUncert2;
	     errSignDo [kproce][ibin] += StatUncert2;
	   }
	 } else {
	   for (int ibin = 1; ibin<=nbins; ibin++) {
	     float StatUncert2 = TMath::Power((dummy0->GetBinError(ibin)+dummy3->GetBinError(ibin))/2., 2);
	     if (StatUncert2<0.0001) StatUncert2 = TMath::Power(StatZero*(dummy0->Integral()+dummy3->Integral())/(dummy0->GetEntries()+dummy3->GetEntries()), 2);
	     errSignUp [kproce][ibin] += StatUncert2;
	     errSignDo [kproce][ibin] += StatUncert2;
	   }
	 }
	 continue;
       }

       if ( _systematics.at(isyst) == "Luminosity") {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float LumiBinError2 = TMath::Power(yieldSign[kproce][ibin]*lumi_error_percent/1e2, 2);
	   errSignUp [kproce][ibin] += LumiBinError2;
	   errSignDo [kproce][ibin] += LumiBinError2;
	 }
	 continue;
       }
      
       if ( _systematics.at(isyst) == "Trigger") {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   float TrigBinError2 = TMath::Power(yieldSign[kproce][ibin]*2./1e2, 2);
	   errSignUp [kproce][ibin] += TrigBinError2;
	   errSignDo [kproce][ibin] += TrigBinError2;
	 }
	 continue;
       }

       if (_systematics.at(isyst)=="Metfastsim")  {
	 for (int ibin=1; ibin<=nbins; ibin++) {
	   errSignUp [kproce][ibin] += TMath::Power(dummy0->GetBinContent(ibin)-yieldSign[kproce][ibin], 2);
	   errSignDo [kproce][ibin] += TMath::Power(dummy0->GetBinContent(ibin)-yieldSign[kproce][ibin], 2);
	 }
	 continue;
       }
       
       TString FileUpName, FileDoName;
       if (_minitreebased) { // This is for ttdm
	 if (isyst%2 != 0) continue; //Only takes the first one systematic_up(down).root as reference
	 FileUpName = _inputdir + "/" + _signalfilename.at(kproce) + _systematics.at(isyst)   + ".root";
	 FileDoName = _inputdir + "/" + _signalfilename.at(kproce) + _systematics.at(isyst+1) + ".root";
       } else {
	 FileUpName = _inputdir + "/../../" + _systematics.at(isyst) + "up/" + _analysis + "/" + _signalfilename.at(kproce) + ".root";
	 FileDoName = _inputdir + "/../../" + _systematics.at(isyst) + "do/" + _analysis + "/" + _signalfilename.at(kproce) + ".root";
	 }
       TFile* myfile1 = TFile::Open(FileUpName);
       TFile* myfile2 = TFile::Open(FileDoName);
       
       TH1D* dummy1 = GetHistogram(myfile1, hname);   //(TH1D*)myfile1->Get( hname );//up
       TH1D* dummy4 = GetHistogram(myfile1, hnamegen);//(TH1D*)myfile1->Get( hnamegen );//up
       TH1D* dummy2 = GetHistogram(myfile2, hname);   //(TH1D*)myfile2->Get( hname );//down 
       TH1D* dummy5 = GetHistogram(myfile2, hnamegen);//(TH1D*)myfile2->Get( hnamegen );//down      
       
       if (_luminosity_fb > 0) {
	 dummy1->Scale(_luminosity_fb);
	 if (_doMetFastSim) dummy4->Scale(_luminosity_fb);
	 dummy2->Scale(_luminosity_fb);
	 if (_doMetFastSim) dummy5->Scale(_luminosity_fb);
       }
		       
       // Loop over all bins (Underflow is not included)
       //--------------------------------------------------------------------  
       for (int ibin=1; ibin<=nbins; ibin++) {

	 if (dummy1->GetBinContent(ibin)<0.) dummy1->SetBinContent(ibin, 0.001);
	 if (dummy2->GetBinContent(ibin)<0.) dummy2->SetBinContent(ibin, 0.001);
	 if (dummy4->GetBinContent(ibin)<0.) dummy4->SetBinContent(ibin, 0.001);
	 if (dummy5->GetBinContent(ibin)<0.) dummy5->SetBinContent(ibin, 0.001);
	 float ErrUp = 0., ErrDo = 0.;
	 float VarUp, VarDo;
	 if (!_doMetFastSim) {
	   VarUp = dummy1->GetBinContent(ibin) - yieldSign[kproce][ibin];
	   VarDo = dummy2->GetBinContent(ibin) - yieldSign[kproce][ibin]; 
	 } else {
	   VarUp = (dummy1->GetBinContent(ibin)+dummy4->GetBinContent(ibin))/2. - yieldSign[kproce][ibin];
	   VarDo = (dummy2->GetBinContent(ibin)+dummy5->GetBinContent(ibin))/2. - yieldSign[kproce][ibin]; 
	 }
	 if (VarDo<=0. && VarUp>=0.) {
	   ErrUp = VarUp;
	   ErrDo = VarDo;
	 } else if (VarDo>=0. && VarUp<=0.) {
	   ErrUp = VarDo;
	   ErrDo = VarUp;
	 } else if (VarDo>0.) {
	   ErrUp = (VarUp>VarDo) ? VarUp : VarDo;
	 } else if (VarDo<0.) {
	   ErrDo = (VarUp<VarDo) ? VarUp : VarDo;
	 }
         errSignUp [kproce][ibin] += TMath::Power(ErrUp, 2);
         errSignDo [kproce][ibin] += TMath::Power(ErrDo, 2);
	 
       }
     
       myfile1->Close();
       myfile2->Close();
     
     }
     
     myfile0->Close();

   }

   // Create the TGraphAsymmErrors
   // ------------------------------
   float x       [nbins+1];
   float y       [nbins+1];
   float exl     [nbins+1];
   float eyl     [nbins+1];
   float exh     [nbins+1];
   float eyh     [nbins+1];
   float errStat [nbins+1];
   float errLumi [nbins+1];
   float errTrig [nbins+1];

   if (_verbose) printf( "--------------------- Printing Errors ----------------------------\n" ); 
   
   TH1D *dummySM;
   if (_isPostfit) {
     TFile *myfileSM = TFile::Open(_inputdir + "/99_TotalBackground.root");
     dummySM = GetHistogram(myfileSM, hname);
   }

   for (int  ibin=1; ibin<=nbins; ibin++)
     {
       if (!_isPostfit)
	 errStat [ibin] = sqrt(_allmchist ->GetSumw2()->At(ibin));
       else 
	 errStat [ibin] = dummySM->GetBinError(ibin);
       x       [ibin] = _allmchist ->GetXaxis()->GetBinCenter(ibin);
       y       [ibin] = _allmchist ->GetBinContent(ibin);
       errLumi [ibin] =  y[ibin] * lumi_error_percent/1e2;
       errTrig [ibin] =  y[ibin] * 2./1e2;
       exl     [ibin] = (_allmchist -> GetXaxis() -> GetBinWidth(ibin))/2;
       exh     [ibin] = exl[ibin];
   
       if (_verbose) {
	 //Print Stat and flat errors per bin
	 printf( "bin number = %i\n", ibin );  
	 printf( "                                                              \n");  
	 printf( "--------------------- Flat errors ----------------------------\n" ); 
	 printf( "---------------------------------------------------------\n");  
	 printf( "                                                              \n");  
	 printf( "errStat = %f\n", errStat [ibin] ); 
	 printf( "errLumi = %f\n", errLumi [ibin] );
	 printf( "errTrig = %f\n", errTrig [ibin] );
	 printf( "rel_errStat = %f %\n", errStat [ibin] / y [ibin] ); 
	 printf( "rel_errLumi = %f %\n", errLumi [ibin] / y [ibin]);
	 printf( "rel_errTrig = %f %\n", errTrig [ibin] / y [ibin]);
	 printf( "                                                              \n");  
	 printf( "--------------------------------------------------------\n" ); 
       }
 
       float systUp2  = 0;
       float systDo2  = 0;
       float systSym2 = 0;

       bool AddLuminosity = false, AddTrigger = false;

       for (int isyst =0; isyst< nsystematics; isyst++)
         {
         
	   if (_systematics.at(isyst)=="Luminosity") {
	     AddLuminosity = true;
	     continue;
	   } else if (_systematics.at(isyst)=="Trigger") {
	     AddTrigger = true;
	     continue;
	   } else if (_systematics.at(isyst)=="Postfit") {
	     continue;
	   }
	   
	   if (_verbose)
	     if ( errSystUp [isyst][ibin] * errSystDo [isyst][ibin] > 0 ) 
	       {
		 printf( "WARNING! errSystUp and errSystDo have the same sign!\n") ; 
		 printf( "The systematic is %s, the bin is %i \n\n", _systematics.at(isyst).Data(), ibin);         
		 printf( "errSystUp = %f\n errSystDo = %f\n", errSystUp [isyst][ibin], errSystDo [isyst][ibin]); 
	       }

	   if (_systematics.at(isyst)=="Statistics") {
	     errSystUp [isyst][ibin] =     sqrt( errSystUp [isyst][ibin]  );
	     errSystDo [isyst][ibin] = -1.*sqrt( errSystDo [isyst][ibin]  );
	   }
	   
	   // Assymetric errors
	   if ( errSystUp [isyst][ibin] < 0 && errSystDo [isyst][ibin] > 0 )
	     {
	       float midErrUp = errSystUp[isyst][ibin]; 
	       float midErrDo = errSystDo[isyst][ibin]; 
	       errSystUp [isyst][ibin] = midErrDo; 
	       errSystDo [isyst][ibin] = midErrUp;
	     }
	   
	   systUp2  += errSystUp [isyst][ibin] * errSystUp [isyst][ibin];    
	   systDo2  += errSystDo [isyst][ibin] * errSystDo [isyst][ibin];
	   // Symmetric errors
	   //systSym2 += sqrt( (errSystUp [isyst][ibin] + errSystDo [isyst][ibin]) * (errSystUp [isyst][ibin] + errSystDo [isyst][ibin]) )  / 2 ;   
          
	   if (_verbose) {
	     // Print Systematic Total Errors (sum over all processe---s)          
	     printf( "                                                        \n");  
	     printf( "----  %s  Systematic Total Error\n", _systematics.at(isyst).Data());  
	     //          printf( "systematic name %s \n", _systematics.at(isyst).Data() ); 
	     printf( "                                                        \n");  
	     printf( "ErrsystUp  = %.5f\n", errSystUp [isyst][ibin] );           
	     printf( "ErrsystDo  = %.5f\n", errSystDo [isyst][ibin] );
	     printf( "rel_ErrsystUp  = %.5f  %\n", errSystUp [isyst][ibin] / y [ibin] );
	     printf( "rel_ErrsystDo  = %.5f  %\n", errSystDo [isyst][ibin] / y [ibin] );       
	     printf( "ErrsystSym = %.5f\n", fabs(errSystUp [isyst][ibin])/2 + fabs(errSystDo [isyst][ibin])/2 );
	     printf( "----------------------------------------------------------\n");  
	     printf( "                                                        \n");  
	   }
	   
         }
       
       if (!AddLuminosity) errLumi[ibin] = 0.;
       if (!AddTrigger) errTrig[ibin] = 0.;
       //eyl [ibin] = sqrt( errStat[ibin]*errStat[ibin] + errLumi[ibin]*errLumi[ibin] + errTrig[ibin]*errTrig[ibin] + systDo2);
       //eyh [ibin] = sqrt( errStat[ibin]*errStat[ibin] + errLumi[ibin]*errLumi[ibin] + errTrig[ibin]*errTrig[ibin] + systUp2);
       eyl [ibin] = sqrt( errLumi[ibin]*errLumi[ibin] + errTrig[ibin]*errTrig[ibin] + systDo2);
       eyh [ibin] = sqrt( errLumi[ibin]*errLumi[ibin] + errTrig[ibin]*errTrig[ibin] + systUp2);
        
       if (_verbose) {
	 //Print Total Error per bin    
	 printf( "----------------------- Total Error -------------------\n"); 
	 printf( "err_up   = %.5f\n", eyh [ibin] );   
	 printf( "err_down = %.5f\n", eyl [ibin] );   
	 printf( "rel_err_down = %.5f %\n", eyl [ibin] / y [ibin] );   
	 printf( "rel_err_up   = %.5f %\n", eyh [ibin] / y [ibin]);   
	 printf( "----------------------------------------------------------\n");  
	  printf( "                                                        \n");  
       }
       
     }

   x[0] = -999.;
   _ErrorGr = new TGraphAsymmErrors(nbins+1,x,y,exl,exh,eyl,eyh);
   
   _ErrorGr->SetMarkerColor(kGray+1);
   _ErrorGr->SetMarkerSize (      0);
   _ErrorGr->SetLineColor  (kGray+1);
   _ErrorGr->SetFillColor  (kGray+1);
   _ErrorGr->SetFillStyle  (   3345);
   
   if (_verbose) {
     _ErrorGr->Print();
     for (int ff = 1; ff<=nbins; ff++)
       cout << "ErrorY[" << ff << "] = " << _ErrorGr->GetErrorY(ff) << endl;
   }

   if (_dotable && hname.Contains("h_MT2ll")) {
     
     TString TableFlag = hname; TableFlag.ReplaceAll("Stop/", "_"); TableFlag.ReplaceAll("/h", "");
     std::ofstream inFile("./Tables/Yields" + TableFlag + ".tex",std::ios::out);
     // Process | nbins = 7;    

     //inFile << "\\begin{table}[htb]" << endl;
     inFile << "\\tiny" << endl;
     inFile << "\\begin{center}" << endl;
     inFile << "\\begin{tabular}{|l|ccccccc|}" << endl;
     inFile << "\\hline" << endl;
     //Yield & stat_error & systematic_error & total_error" << endl;
     inFile << "\\hline" << endl;
     inFile << "$M_{T2}(ll})$ bin &";
     for (int ibin=1; ibin<nbins; ibin++) 
       inFile <<  (ibin-1)*20 <<"-"<< ibin*20 << "~\\GeV & ";
     inFile << "\\ge " << (nbins-1)*20 << "~\\GeV \\\\" << endl;
     inFile << "\\hline" << endl;
     for (int kproce=0; kproce<nprocess; kproce++) {       
       TString ThisLabel = _mclabel[kproce].Data();
       ThisLabel.ReplaceAll("#", "\\");
       inFile << "$" << ThisLabel << "$";  
       //inFile << _mclabel[kproce].Data();     
       for (int ibin=1; ibin<=nbins; ibin++) {
	 float ThisYield = yieldTab[kproce][ibin];
	 float ThisError = (sqrt(errBackTab_up[kproce][ibin])+sqrt(errBackTab_do[kproce][ibin]))/2.;
	 FormatTableYields(&ThisYield, &ThisError);
	 inFile << " & $" << ThisYield <<  " \\pm " << ThisError << "$";
       }
       inFile << " \\\\" << endl;
     }
     inFile << "\\hline" << endl;
     inFile << "SM Processes ";
     for (int ibin=1; ibin<=nbins; ibin++) {
       float ThisYield = y[ibin];
       float ThisError = _ErrorGr->GetErrorY(ibin);
       FormatTableYields(&ThisYield, &ThisError);
       inFile << " &  $" << ThisYield << " \\pm " << ThisError << "$";
     }
     inFile << " \\\\" << endl;
     inFile << " \\hline" << endl;
     inFile << " Data ";
     for (int ibin=1; ibin<=nbins; ibin++)
       if (_datahist)
	 inFile << " & $" << _datahist->GetBinContent(ibin) << "$";
       else
	 inFile << " & $" << "blind" << "$";
     inFile << " \\\\" << endl;
     inFile << " \\hline" << endl;
     for (int kproce=0; kproce<nsignals; kproce++) {
       TString ThisLabel = _signallabel[kproce].Data();
       ThisLabel.ReplaceAll("m_{#tilde{t}}=", "");
       ThisLabel.ReplaceAll("m_{#tilde{#chi}^{0}_{1}}=", "");   
       ThisLabel.ReplaceAll("#", "\\");
       inFile << "$" << ThisLabel << "$";  
       for (int ibin=1; ibin<=nbins; ibin++) {	
	 float ThisYield = yieldSign[kproce][ibin];
	 float ThisError = (sqrt(errSignUp[kproce][ibin])+sqrt(errSignDo[kproce][ibin]))/2.;
	 FormatTableYields(&ThisYield, &ThisError);
	 inFile << " & $" << ThisYield <<  " \\pm " << ThisError << "$";
       }
       inFile << " \\\\" << endl;
     }
     inFile << " \\hline\\hline" << endl;
     inFile << "\\end{tabular}" << endl;
     inFile << "\\end{center}" << endl;
     //inFile << "\\end{table}" << endl;
     
     inFile.close();
   }
   
}
