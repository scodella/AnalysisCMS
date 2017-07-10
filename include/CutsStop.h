#if !defined(MYLIB_CUTSSTOP_H)
#define MYLIB_CUTSSTOP_H 1

#include "TString.h"


enum {
  Stop_00_Has2Leptons,
//  Stop_00_mll20,
//~~~~~~~~~~~~~~~~~~~~~~~ Z peak sigth ~~~~~~~~~~~~~~~~~~~
//    Look in the Z-peak
//
  Stop_00_VZ,
  Stop_00_ttZ_Tag,
  Stop_00_ttZ_2Jet,
  Stop_05_Zpeak,
//  Stop_05_NoTagZpeak,
//  Stop_05_TagZpeak,
//
//     DY estimation Rout/Rin
//
//  Stop_04_Routin,
//  Stop_04_Jet2Routin,
//  Stop_04_TagRoutin,
//  Stop_04_NoTagRoutin,
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Stop_00_SS_Tag,
  Stop_00_SS_2Jet,
  Stop_00_SS,
  Stop_00_Zveto,
 //~~~~~~~~~~~~~~~~~~~~~~ save minitree ~~~~~~~~~~~~~~~~~~
 //
 //~~~~~~~~~~~~~~~~~~~~~~ MC normalization study ~~~~~~~~~
  Stop_00_WWsel,
  Stop_00_WWselMET,
  Stop_00_TTsel,
  Stop_00_TTselMET,

 //~~~~~~~~~~~~~~~~~~~~~~ Analysis Cuts ~~~~~~~~~~~~~~~~~~
  Stop_01_Tag,
  Stop_01_NoTag,
  Stop_02_VR1_Tag,
  Stop_02_VR1_NoTag,
  Stop_02_VR1_Tag2Jet,
  Stop_02_VR1_NoJet,
  Stop_02_SR1_Tag,
  Stop_02_SR1_NoTag,
  Stop_02_SR2_Tag,
  Stop_02_SR2_NoTag,
  Stop_02_SR3_Tag,
  Stop_02_SR3_NoTag,
  Stop_02_SR1gen_Tag,
  Stop_02_SR1gen_NoTag,
  Stop_02_SR2gen_Tag,
  Stop_02_SR2gen_NoTag,
  Stop_02_SR3gen_Tag,
  Stop_02_SR3gen_NoTag,
  Stop_02_SRs,
  Stop_02_SRs_Tag,
  Stop_02_SRs_NoTag,

  ncut  // This line should be always last
};

const TString scut[ncut] = {
  "Stop/00_Has2Leptons",
// "Stop/00_mll20",
//~~~~~~~~~~~~~~~~~~~~~~~ Z peak sigth ~~~~~~~~~~~~~~~~~~~
//    Look in the Z-peak
  "Stop/00_VZ",
  "Stop/00_ttZ_Tag",
  "Stop/00_ttZ_2Jet",
  "Stop/05_Zpeak",
//  "Stop/05_NoTagZpeak",
//  "Stop/05_TagZpeak",
//
//     DY estimation Rout/Rin
//
//  "Stop/04_Routin",
//  "Stop/04_Jet2Routin",
//  "Stop/04_TagRoutin",
//  "Stop/04_NoTagRoutin",
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   "Stop/00_SS_Tag",
   "Stop/00_SS_2Jet",
   "Stop/00_SS",
   "Stop/00_Zveto",
 //~~~~~~~~~~~~~~~~~~~~~~ save minitree ~~~~~~~~~~~~~~~~~~
 //
 //~~~~~~~~~~~~~~~~~~~~~~ MC normalization study ~~~~~~~~~
   "Stop/00_WWsel",
   "Stop/00_WWselMET",
   "Stop/00_TTsel", 
   "Stop/00_TTselMET",

 //~~~~~~~~~~~~~~~~~~~~~~ Analysis Cuts ~~~~~~~~~~~~~~~~~~
  "Stop/01_Tag",
  "Stop/01_NoTag",
  "Stop/02_VR1_Tag",
  "Stop/02_VR1_NoTag",
  "Stop/02_VR1_Tag2Jet",
  "Stop/02_VR1_NoJet",
  "Stop/02_SR1_Tag",
  "Stop/02_SR1_NoTag",
  "Stop/02_SR2_Tag",
  "Stop/02_SR2_NoTag",
  "Stop/02_SR3_Tag",
  "Stop/02_SR3_NoTag",
  "Stop/02_SR1gen_Tag",
  "Stop/02_SR1gen_NoTag",
  "Stop/02_SR2gen_Tag",
  "Stop/02_SR2gen_NoTag",
  "Stop/02_SR3gen_Tag",
  "Stop/02_SR3gen_NoTag",
  "Stop/02_SRs",
  "Stop/02_SRs_Tag",
  "Stop/02_SRs_NoTag"

};

#endif
