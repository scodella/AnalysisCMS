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
mv figures ../Plots/PlotterSummer16T2ttV4
mv Tables TablesT2ttV4
