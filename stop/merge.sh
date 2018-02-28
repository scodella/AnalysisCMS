#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "  "
    echo "  ./merge.sh ../rootfiles/<systematic>/<analysis>"
    echo "  ./merge.sh ../minitrees/<systematic>/<analysis>"
    echo "  "
    exit -1
fi

FOLDER="$1"

pushd $FOLDER

# Data
hadd -f -k 01_Data.root       *03Feb2017*
##hadd -f -k 01_DataBlind.root   *Run2016B-03Feb2017* *Run2016C-03Feb2017* *Run2016D-03Feb2017* 

# Top
hadd -f -k 04_TTTo2L2Nu.root   TTTo2L2Nu__part*.root
hadd -f -k 05_ST.root          ST_tW_antitop.root ST_tW_top.root

# ZJets

#hadd -f -k 07_ZJets.root      DYJetsToLL_M-10to50.root DYJetsToLL_M-50__part*.root # NLO
hadd -f -k 07_ZJetsHT.root     DYJetsToLL_M-10to50-LO.root DYJetsToLL_M-5to50_HT-*root DYJetsToLL_M-50-LO-ext1__part*.root DYJetsToLL_M-50_HT-70to100__part*.root DYJetsToLL_M-50_HT-100to200_ext1__part*.root DYJetsToLL_M-50_HT-200to400_ext1__part*.root DYJetsToLL_M-50_HT-400to600.root DYJetsToLL_M-50_HT-600to800__part*.root DYJetsToLL_M-50_HT-800to1200.root DYJetsToLL_M-50_HT-1200to2500.root DYJetsToLL_M-50_HT-2500toInf.root 

hadd -f -k 07_ZJetsHT_DYcorr.root   DYJetsToLL_M-10to50-LO_DYcorr.root DYJetsToLL_M-5to50_HT-*_DYcorr.root DYJetsToLL_M-50-LO-ext1__part*_DYcorr.root DYJetsToLL_M-50_HT-70to100__part*_DYcorr.root DYJetsToLL_M-50_HT-100to200_ext1__part*_DYcorr.root DYJetsToLL_M-50_HT-200to400_ext1__part*_DYcorr.root DYJetsToLL_M-50_HT-400to600_DYcorr.root DYJetsToLL_M-50_HT-600to800__part*_DYcorr.root DYJetsToLL_M-50_HT-800to1200_DYcorr.root DYJetsToLL_M-50_HT-1200to2500_DYcorr.root DYJetsToLL_M-50_HT-2500toInf_DYcorr.root 

# Dibosons

hadd -f -k 02_WZTo3LNu.root    WZTo3LNu.root
hadd -f -k 03_ZZ.root          ZZTo2L2Nu__part* ggZZTo2e2nu.root ggZZTo2mu2nu.root   
hadd -f -k 06_WW.root          WWTo2L2Nu.root GluGluWWTo2L2Nu_MCFM.root
hadd -f -k 14_ZZTo4L.root      ZZTo4L__part* ggHToZZTo4L.root ggZZTo2e2mu.root ggZZTo2e2nu.root ggZZTo2e2tau.root ggZZTo2mu2nu.root ggZZTo2mu2tau.root  ggZZTo4e.root ggZZTo4mu.root ggZZTo4tau.root qqHToZZTo4L.root
hadd -f -k 15_VZ.root          ZZTo2L2Q__part*.root  WZTo2L2Q__part*.root

# Tribosons
hadd -f -k 09_TTW.root         TTWJetsToLNu.root TTWJetsToQQ.root 
hadd -f -k 10_TTZ.root         TTZToQQ.root TTZToLLNuNu_M-10.root 
hadd -f -k 11_HWW.root         GluGluHToWWTo2L2NuAMCNLO_M125.root VBFHToWWTo2L2Nu_M125.root GluGluHToTauTau_M125.root VBFHToTauTau_M125.root HWminusJ_HToWW_M125.root HWplusJ_HToWW_M125.root
hadd -f -k 13_VVV.root         WWW.root WWZ.root WZZ.root #ZZZ.root
hadd -f -k 15_VZV3.root        15_VZ.root 13_VVV.root
popd
