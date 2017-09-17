./make
mkdir -p Tables
rm Tables/*
rm -r figures
#./runPlotter Stop/02_SR1_NoJet
#./runPlotter Stop/02_SR2_NoJet
#./runPlotter Stop/02_SR1_NoTag
#./runPlotter Stop/02_SR2_NoTag
./runPlotter Stop/02_SR3_Veto
#./runPlotter Stop/02_SR1_Tag
#./runPlotter Stop/02_SR2_Tag
#./runPlotter Stop/02_SR3_Tag
#mv figures ../Plots/PlotterSummer16TChiV4
#mv Tables TablesTChiV4
mv figures/Stop/02_SR3_Veto/h_* ../Plots/PlotterSummer16TChiV4/Stop/02_SR3_Veto/
mv Tables/*Veto* TablesTChiV4/
