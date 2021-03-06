//------------------------------------------------------------------------------
//
// root -l nice2d.C
//
//------------------------------------------------------------------------------
void AxisFonts     (TAxis*  axis,
		    TString coordinate,
		    TString title);

void TH2FAxisFonts (TH2F*   h,
		    TString coordinate,
		    TString title);

void setupGreyScale();

void setupColors   ();


//------------------------------------------------------------------------------
//
// main
//
//------------------------------------------------------------------------------
void nice2d(Bool_t isBlackAndWhite = false)
{
  gInterpreter->ExecuteMacro("../test/PaperStyle.C");
  
  gStyle->SetOptStat(0);

  (isBlackAndWhite) ? setupGreyScale() : setupColors();


  // Fill a dummy 2D histogram
  //----------------------------------------------------------------------------
  TF2* f2 = new TF2("f2", "[0]*exp(-(x*y-[1])*(x*y-[1])/[2])+[3]*exp(-(y-[4])*(y-[4])/[5])", -10, 10, -10, 10);

  f2->SetParameters(11.2, 7.3, 1.4, 22.2, 3.9, 0.5);

  TH2F* h2 = new TH2F("h2", "", 500, -10, 10, 500, -10, 10);

  h2->FillRandom("f2", 50000);


  // Draw
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);

  canvas->SetRightMargin(3.5 * canvas->GetRightMargin());

  TH2FAxisFonts(h2, "x", "trees");
  TH2FAxisFonts(h2, "y", "flowers");

  h2->Draw("colz");

  canvas->Update();

  TPaletteAxis* palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");

  palette->SetLabelFont(42);

  //  TProfile* profile = h2->ProfileX();

  //  profile->Draw("same");

  canvas->Update();

  canvas->Modified();

  canvas->GetFrame()->DrawClone();
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
	       TString coordinate,
	       TString title)
{
  axis->SetLabelFont  (   42);
  axis->SetLabelOffset(0.015);
  axis->SetLabelSize  (0.050);
  axis->SetNdivisions (  505);
  axis->SetTitleFont  (   42);
  axis->SetTitleOffset(  1.5);
  axis->SetTitleSize  (0.050);

  if (coordinate == "y") axis->SetTitleOffset(1.6);

  axis->SetTitle(title);
}


//------------------------------------------------------------------------------
// TH2FAxisFonts
//------------------------------------------------------------------------------
void TH2FAxisFonts(TH2F*   h,
		   TString coordinate,
		   TString title)
{
  TAxis* axis = NULL;

  if (coordinate.Contains("x")) axis = h->GetXaxis();
  if (coordinate.Contains("y")) axis = h->GetYaxis();

  AxisFonts(axis, coordinate, title);
}


//------------------------------------------------------------------------------
// setupGreyScale
//------------------------------------------------------------------------------
void setupGreyScale() 
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red  [NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  double green[NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  double blue [NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}


//------------------------------------------------------------------------------
// setupColors
//------------------------------------------------------------------------------
void setupColors()
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = {0.00, 0.0625, 0.25, 0.5625, 1.00};
  Double_t red  [NRGBs] = {0.00, 0.00,   0.87, 1.00,   0.51};
  Double_t green[NRGBs] = {0.00, 0.81,   1.00, 0.20,   0.00};
  Double_t blue [NRGBs] = {0.51, 1.00,   0.12, 0.00,   0.00};
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}
