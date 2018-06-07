./make
mkdir -p Tables
rm Tables/*
rm -r figures
./runPlotter Stop/02_SR1_Veto
./runPlotter Stop/02_SR2_Veto
./runPlotter Stop/02_SR3_Veto
./runPlotter Stop/02_SR1_Tag
./runPlotter Stop/02_SR2_Tag
./runPlotter Stop/02_SR3_Tag
rm -r  ../Plots/PlotterSummer16T2ttCWR
mv figures ../Plots/PlotterSummer16T2ttCWR
#mkdir -p TablesCWR
#rm TablesCWR/*
#mv Tables/* TablesCWR/


