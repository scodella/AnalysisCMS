./make
mkdir -p Tables
rm Tables/*
rm -r figures
./runPlotter Stop/02_SR1_NoJet
./runPlotter Stop/02_SR2_NoJet
./runPlotter Stop/02_SR1_NoTag
./runPlotter Stop/02_SR2_NoTag
./runPlotter Stop/02_SR3_Veto
#./runPlotter Stop/02_SR1_Tag
#./runPlotter Stop/02_SR2_Tag
#./runPlotter Stop/02_SR3_Tag
rm -r ../Plots/PlotterSummer16TChiCWR
mv figures ../Plots/PlotterSummer16TChiCWR
#mkdir -p TablesCWR
#rm TablesCWR/*
#mv Tables/* TablesCWR/

