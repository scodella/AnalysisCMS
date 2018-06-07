#include <map>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <fstream>
#include <iostream>


float TwoLFilterEff = 0.104970;
bool Apply2LFilterEff = false;

float susyMstop, susyMLSP, susyMChargino;
float nisrjet, eventW;
vector<float>   *std_vector_jet_pt;
vector<float>   *std_vector_jet_eta;
vector<float>   *std_vector_jet_phi;
vector<float>   *std_vector_jet_isFromISR;
vector<float>   *std_vector_partonGen_pt;
vector<float>   *std_vector_partonGen_eta;
vector<float>   *std_vector_partonGen_phi;
vector<float>   *std_vector_partonGen_pid;
vector<float>   *std_vector_partonGen_isHardProcess;
vector<float>   *std_vector_VBoson_pt;
vector<float>   *std_vector_VBoson_eta;
vector<float>   *std_vector_VBoson_phi;
vector<float>   *std_vector_VBoson_pid;
vector<float>   *std_vector_leptonGen_MotherPID;
vector<float>   *std_vector_lepton_isTightLepton;

TTree *GetTree(TFile *TreeFile, TString TreeType = "latino") {

  TTree *ThisTree = (TTree*) TreeFile->Get("latino");
  
  ThisTree->SetBranchAddress("susyMLSP",       &susyMLSP);
  ThisTree->SetBranchAddress("susyMstop",      &susyMstop);
  ThisTree->SetBranchAddress("susyMChargino",  &susyMChargino);

  if (TreeType=="latino") {

    ThisTree->SetBranchAddress("std_vector_jet_pt", &std_vector_jet_pt);
    ThisTree->SetBranchAddress("std_vector_jet_eta", &std_vector_jet_eta);
    ThisTree->SetBranchAddress("std_vector_jet_phi", &std_vector_jet_phi);
    ThisTree->SetBranchAddress("std_vector_jet_isFromISR", &std_vector_jet_isFromISR);
    ThisTree->SetBranchAddress("std_vector_partonGen_pt", &std_vector_partonGen_pt);
    ThisTree->SetBranchAddress("std_vector_partonGen_eta", &std_vector_partonGen_eta);
    ThisTree->SetBranchAddress("std_vector_partonGen_phi", &std_vector_partonGen_phi);
    ThisTree->SetBranchAddress("std_vector_partonGen_pid", &std_vector_partonGen_pid);
    ThisTree->SetBranchAddress("std_vector_partonGen_isHardProcess", &std_vector_partonGen_isHardProcess);
    ThisTree->SetBranchAddress("std_vector_VBoson_pt", &std_vector_VBoson_pt);
    ThisTree->SetBranchAddress("std_vector_VBoson_eta", &std_vector_VBoson_eta);
    ThisTree->SetBranchAddress("std_vector_VBoson_phi", &std_vector_VBoson_phi);
    ThisTree->SetBranchAddress("std_vector_VBoson_pid", &std_vector_VBoson_pid);
    ThisTree->SetBranchAddress("std_vector_leptonGen_MotherPID", &std_vector_leptonGen_MotherPID);
    ThisTree->SetBranchAddress("std_vector_lepton_isTightLepton", &std_vector_lepton_isTightLepton);
    
  } else {

    ThisTree->SetBranchAddress("nisrjet", &nisrjet);
    ThisTree->SetBranchAddress("eventW",  &eventW);
    
  } 

  return ThisTree;
  
}

float CrossSection, CrossSectionUncertainty;

void GetStopCrossSection(float StopMass) {

  int iStopMass = StopMass;
 
  if (StopMass<100. || StopMass>2000. || iStopMass%5!=0) {
    CrossSection = -1.; return;
  }

  if (iStopMass==100) { CrossSection = 	1521.11; CrossSectionUncertainty = CrossSection* 15.4038/100.; }
  if (iStopMass==105) { CrossSection = 	1233.18; CrossSectionUncertainty = CrossSection* 15.4059/100.; }
  if (iStopMass==110) { CrossSection = 	1013.76; CrossSectionUncertainty = CrossSection* 15.4088/100.; }
  if (iStopMass==115) { CrossSection = 	832.656; CrossSectionUncertainty = CrossSection* 15.1503/100.; }
  if (iStopMass==120) { CrossSection = 	689.799; CrossSectionUncertainty = CrossSection* 15.044/100.; }
  if (iStopMass==125) { CrossSection = 	574.981; CrossSectionUncertainty = CrossSection* 14.9895/100.; }
  if (iStopMass==130) { CrossSection = 	481.397; CrossSectionUncertainty = CrossSection* 14.8906/100.; }
  if (iStopMass==135) { CrossSection = 	405.159; CrossSectionUncertainty = CrossSection* 14.8952/100.; }
  if (iStopMass==140) { CrossSection = 	342.865; CrossSectionUncertainty = CrossSection* 14.9119/100.; }
  if (iStopMass==145) { CrossSection = 	291.752; CrossSectionUncertainty = CrossSection* 14.8022/100.; }
  if (iStopMass==150) { CrossSection = 	249.409; CrossSectionUncertainty = CrossSection* 14.7477/100.; }
  if (iStopMass==155) { CrossSection = 	214.221; CrossSectionUncertainty = CrossSection* 14.5928/100.; }
  if (iStopMass==160) { CrossSection = 	184.623; CrossSectionUncertainty = CrossSection* 14.5821/100.; }
  if (iStopMass==165) { CrossSection = 	159.614; CrossSectionUncertainty = CrossSection* 14.7859/100.; }
  if (iStopMass==170) { CrossSection = 	139.252; CrossSectionUncertainty = CrossSection* 14.547/100.; }
  if (iStopMass==175) { CrossSection = 	121.416; CrossSectionUncertainty = CrossSection* 14.6341/100.; }
  if (iStopMass==180) { CrossSection = 	106.194; CrossSectionUncertainty = CrossSection* 14.2033/100.; }
  if (iStopMass==185) { CrossSection = 	93.3347; CrossSectionUncertainty = CrossSection* 14.4893/100.; }
  if (iStopMass==190) { CrossSection = 	82.2541; CrossSectionUncertainty = CrossSection* 14.4677/100.; }
  if (iStopMass==195) { CrossSection = 	72.7397; CrossSectionUncertainty = CrossSection* 14.4452/100.; }
  if (iStopMass==200) { CrossSection = 	64.5085; CrossSectionUncertainty = CrossSection* 14.4098/100.; }
  if (iStopMass==205) { CrossSection = 	57.2279; CrossSectionUncertainty = CrossSection* 14.4191/100.; }
  if (iStopMass==210) { CrossSection = 	50.9226; CrossSectionUncertainty = CrossSection* 14.2457/100.; }
  if (iStopMass==215) { CrossSection = 	45.3761; CrossSectionUncertainty = CrossSection* 14.344/100.; }
  if (iStopMass==220) { CrossSection = 	40.5941; CrossSectionUncertainty = CrossSection* 14.2634/100.; }
  if (iStopMass==225) { CrossSection = 	36.3818; CrossSectionUncertainty = CrossSection* 14.2189/100.; }
  if (iStopMass==230) { CrossSection = 	32.6679; CrossSectionUncertainty = CrossSection* 14.1592/100.; }
  if (iStopMass==235) { CrossSection = 	29.3155; CrossSectionUncertainty = CrossSection* 14.2233/100.; }
  if (iStopMass==240) { CrossSection = 	26.4761; CrossSectionUncertainty = CrossSection* 14.1723/100.; }
  if (iStopMass==245) { CrossSection = 	23.8853; CrossSectionUncertainty = CrossSection* 13.9482/100.; }
  if (iStopMass==250) { CrossSection = 	21.5949; CrossSectionUncertainty = CrossSection* 14.0595/100.; }
  if (iStopMass==255) { CrossSection = 	19.5614; CrossSectionUncertainty = CrossSection* 13.8755/100.; }
  if (iStopMass==260) { CrossSection = 	17.6836; CrossSectionUncertainty = CrossSection* 13.9505/100.; }
  if (iStopMass==265) { CrossSection = 	16.112; CrossSectionUncertainty = CrossSection* 13.9531/100.; }
  if (iStopMass==270) { CrossSection = 	14.6459; CrossSectionUncertainty = CrossSection* 13.9278/100.; }
  if (iStopMass==275) { CrossSection = 	13.3231; CrossSectionUncertainty = CrossSection* 14.2549/100.; }
  if (iStopMass==280) { CrossSection = 	12.1575; CrossSectionUncertainty = CrossSection* 14.1584/100.; }
  if (iStopMass==285) { CrossSection = 	11.0925; CrossSectionUncertainty = CrossSection* 14.0904/100.; }
  if (iStopMass==290) { CrossSection = 	10.1363; CrossSectionUncertainty = CrossSection* 13.8967/100.; }
  if (iStopMass==295) { CrossSection = 	9.29002; CrossSectionUncertainty = CrossSection* 13.9107/100.; }
  if (iStopMass==300) { CrossSection = 	8.51615; CrossSectionUncertainty = CrossSection* 13.9223/100.; }
  if (iStopMass==305) { CrossSection = 	7.81428; CrossSectionUncertainty = CrossSection* 13.8996/100.; }
  if (iStopMass==310) { CrossSection = 	7.17876; CrossSectionUncertainty = CrossSection* 13.9357/100.; }
  if (iStopMass==315) { CrossSection = 	6.60266; CrossSectionUncertainty = CrossSection* 13.9256/100.; }
  if (iStopMass==320) { CrossSection = 	6.08444; CrossSectionUncertainty = CrossSection* 13.7957/100.; }
  if (iStopMass==325) { CrossSection = 	5.60471; CrossSectionUncertainty = CrossSection* 13.8144/100.; }
  if (iStopMass==330) { CrossSection = 	5.17188; CrossSectionUncertainty = CrossSection* 13.6954/100.; }
  if (iStopMass==335) { CrossSection = 	4.77871; CrossSectionUncertainty = CrossSection* 13.7554/100.; }
  if (iStopMass==340) { CrossSection = 	4.41629; CrossSectionUncertainty = CrossSection* 13.7945/100.; }
  if (iStopMass==345) { CrossSection = 	4.08881; CrossSectionUncertainty = CrossSection* 13.7075/100.; }
  if (iStopMass==350) { CrossSection = 	3.78661; CrossSectionUncertainty = CrossSection* 13.6877/100.; }
  if (iStopMass==355) { CrossSection = 	3.50911; CrossSectionUncertainty = CrossSection* 13.8089/100.; }
  if (iStopMass==360) { CrossSection = 	3.25619; CrossSectionUncertainty = CrossSection* 13.8002/100.; }
  if (iStopMass==365) { CrossSection = 	3.02472; CrossSectionUncertainty = CrossSection* 13.7093/100.; }
  if (iStopMass==370) { CrossSection = 	2.8077; CrossSectionUncertainty = CrossSection* 13.8064/100.; }
  if (iStopMass==375) { CrossSection = 	2.61162; CrossSectionUncertainty = CrossSection* 13.8477/100.; }
  if (iStopMass==380) { CrossSection = 	2.43031; CrossSectionUncertainty = CrossSection* 13.6999/100.; }
  if (iStopMass==385) { CrossSection = 	2.26365; CrossSectionUncertainty = CrossSection* 13.728/100.; }
  if (iStopMass==390) { CrossSection = 	2.10786; CrossSectionUncertainty = CrossSection* 13.732/100.; }
  if (iStopMass==395) { CrossSection = 	1.9665; CrossSectionUncertainty = CrossSection* 13.4737/100.; }
  if (iStopMass==400) { CrossSection = 	1.83537; CrossSectionUncertainty = CrossSection* 13.6985/100.; }
  if (iStopMass==405) { CrossSection = 	1.70927; CrossSectionUncertainty = CrossSection* 13.7114/100.; }
  if (iStopMass==410) { CrossSection = 	1.60378; CrossSectionUncertainty = CrossSection* 13.5468/100.; }
  if (iStopMass==415) { CrossSection = 	1.49798; CrossSectionUncertainty = CrossSection* 13.4453/100.; }
  if (iStopMass==420) { CrossSection = 	1.39688; CrossSectionUncertainty = CrossSection* 13.6719/100.; }
  if (iStopMass==425) { CrossSection = 	1.31169; CrossSectionUncertainty = CrossSection* 13.5013/100.; }
  if (iStopMass==430) { CrossSection = 	1.22589; CrossSectionUncertainty = CrossSection* 13.3237/100.; }
  if (iStopMass==435) { CrossSection = 	1.14553; CrossSectionUncertainty = CrossSection* 13.5478/100.; }
  if (iStopMass==440) { CrossSection = 	1.07484; CrossSectionUncertainty = CrossSection* 13.7238/100.; }
  if (iStopMass==445) { CrossSection = 	1.01019; CrossSectionUncertainty = CrossSection* 13.4187/100.; }
  if (iStopMass==450) { CrossSection = 	0.948333; CrossSectionUncertainty = CrossSection* 13.4559/100.; }
  if (iStopMass==455) { CrossSection = 	0.890847; CrossSectionUncertainty = CrossSection* 13.4587/100.; }
  if (iStopMass==460) { CrossSection = 	0.836762; CrossSectionUncertainty = CrossSection* 13.4468/100.; }
  if (iStopMass==465) { CrossSection = 	0.787221; CrossSectionUncertainty = CrossSection* 13.4149/100.; }
  if (iStopMass==470) { CrossSection = 	0.740549; CrossSectionUncertainty = CrossSection* 13.4127/100.; }
  if (iStopMass==475) { CrossSection = 	0.697075; CrossSectionUncertainty = CrossSection* 13.3926/100.; }
  if (iStopMass==480) { CrossSection = 	0.655954; CrossSectionUncertainty = CrossSection* 13.4392/100.; }
  if (iStopMass==485) { CrossSection = 	0.618562; CrossSectionUncertainty = CrossSection* 13.3705/100.; }
  if (iStopMass==490) { CrossSection = 	0.582467; CrossSectionUncertainty = CrossSection* 13.3914/100.; }
  if (iStopMass==495) { CrossSection = 	0.549524; CrossSectionUncertainty = CrossSection* 13.3691/100.; }
  if (iStopMass==500) { CrossSection = 	0.51848; CrossSectionUncertainty = CrossSection* 13.3797/100.; }
  if (iStopMass==505) { CrossSection = 	0.489324; CrossSectionUncertainty = CrossSection* 13.3608/100.; }
  if (iStopMass==510) { CrossSection = 	0.462439; CrossSectionUncertainty = CrossSection* 13.3046/100.; }
  if (iStopMass==515) { CrossSection = 	0.436832; CrossSectionUncertainty = CrossSection* 13.3703/100.; }
  if (iStopMass==520) { CrossSection = 	0.412828; CrossSectionUncertainty = CrossSection* 13.272/100.; }
  if (iStopMass==525) { CrossSection = 	0.390303; CrossSectionUncertainty = CrossSection* 13.3443/100.; }
  if (iStopMass==530) { CrossSection = 	0.368755; CrossSectionUncertainty = CrossSection* 13.3769/100.; }
  if (iStopMass==535) { CrossSection = 	0.348705; CrossSectionUncertainty = CrossSection* 13.2706/100.; }
  if (iStopMass==540) { CrossSection = 	0.330157; CrossSectionUncertainty = CrossSection* 13.2981/100.; }
  if (iStopMass==545) { CrossSection = 	0.312672; CrossSectionUncertainty = CrossSection* 13.277/100.; }
  if (iStopMass==550) { CrossSection = 	0.296128; CrossSectionUncertainty = CrossSection* 13.2687/100.; }
  if (iStopMass==555) { CrossSection = 	0.280734; CrossSectionUncertainty = CrossSection* 13.2363/100.; }
  if (iStopMass==560) { CrossSection = 	0.266138; CrossSectionUncertainty = CrossSection* 13.193/100.; }
  if (iStopMass==565) { CrossSection = 	0.251557; CrossSectionUncertainty = CrossSection* 13.1731/100.; }
  if (iStopMass==570) { CrossSection = 	0.238537; CrossSectionUncertainty = CrossSection* 13.3409/100.; }
  if (iStopMass==575) { CrossSection = 	0.226118; CrossSectionUncertainty = CrossSection* 13.2741/100.; }
  if (iStopMass==580) { CrossSection = 	0.214557; CrossSectionUncertainty = CrossSection* 13.1697/100.; }
  if (iStopMass==585) { CrossSection = 	0.203566; CrossSectionUncertainty = CrossSection* 13.3257/100.; }
  if (iStopMass==590) { CrossSection = 	0.193079; CrossSectionUncertainty = CrossSection* 13.2037/100.; }
  if (iStopMass==595) { CrossSection = 	0.183604; CrossSectionUncertainty = CrossSection* 13.0973/100.; }
  if (iStopMass==600) { CrossSection = 	0.174599; CrossSectionUncertainty = CrossSection* 13.2074/100.; }
  if (iStopMass==605) { CrossSection = 	0.166131; CrossSectionUncertainty = CrossSection* 13.0154/100.; }
  if (iStopMass==610) { CrossSection = 	0.158242; CrossSectionUncertainty = CrossSection* 13.142/100.; }
  if (iStopMass==615) { CrossSection = 	0.150275; CrossSectionUncertainty = CrossSection* 13.285/100.; }
  if (iStopMass==620) { CrossSection = 	0.142787; CrossSectionUncertainty = CrossSection* 13.0642/100.; }
  if (iStopMass==625) { CrossSection = 	0.136372; CrossSectionUncertainty = CrossSection* 12.7962/100.; }
  if (iStopMass==630) { CrossSection = 	0.129886; CrossSectionUncertainty = CrossSection* 13.2957/100.; }
  if (iStopMass==635) { CrossSection = 	0.123402; CrossSectionUncertainty = CrossSection* 13.016/100.; }
  if (iStopMass==640) { CrossSection = 	0.11795; CrossSectionUncertainty = CrossSection* 12.7132/100.; }
  if (iStopMass==645) { CrossSection = 	0.112008; CrossSectionUncertainty = CrossSection* 12.808/100.; }
  if (iStopMass==650) { CrossSection = 	0.107045; CrossSectionUncertainty = CrossSection* 12.9232/100.; }
  if (iStopMass==655) { CrossSection = 	0.102081; CrossSectionUncertainty = CrossSection* 13.0012/100.; }
  if (iStopMass==660) { CrossSection = 	0.09725; CrossSectionUncertainty = CrossSection* 12.9038/100.; }
  if (iStopMass==665) { CrossSection = 	0.0927515; CrossSectionUncertainty = CrossSection* 12.9548/100.; }
  if (iStopMass==670) { CrossSection = 	0.0885084; CrossSectionUncertainty = CrossSection* 13.0218/100.; }
  if (iStopMass==675) { CrossSection = 	0.0844877; CrossSectionUncertainty = CrossSection* 13.0703/100.; }
  if (iStopMass==680) { CrossSection = 	0.0806192; CrossSectionUncertainty = CrossSection* 13.1131/100.; }
  if (iStopMass==685) { CrossSection = 	0.0769099; CrossSectionUncertainty = CrossSection* 13.1517/100.; }
  if (iStopMass==690) { CrossSection = 	0.0734901; CrossSectionUncertainty = CrossSection* 13.2344/100.; }
  if (iStopMass==695) { CrossSection = 	0.0701805; CrossSectionUncertainty = CrossSection* 13.2716/100.; }
  if (iStopMass==700) { CrossSection = 	0.0670476; CrossSectionUncertainty = CrossSection* 13.3429/100.; }
  if (iStopMass==705) { CrossSection = 	0.0641426; CrossSectionUncertainty = CrossSection* 13.363/100.; }
  if (iStopMass==710) { CrossSection = 	0.0612942; CrossSectionUncertainty = CrossSection* 13.3941/100.; }
  if (iStopMass==715) { CrossSection = 	0.0585678; CrossSectionUncertainty = CrossSection* 13.4663/100.; }
  if (iStopMass==720) { CrossSection = 	0.0560753; CrossSectionUncertainty = CrossSection* 13.4984/100.; }
  if (iStopMass==725) { CrossSection = 	0.0536438; CrossSectionUncertainty = CrossSection* 13.5804/100.; }
  if (iStopMass==730) { CrossSection = 	0.0513219; CrossSectionUncertainty = CrossSection* 13.5682/100.; }
  if (iStopMass==735) { CrossSection = 	0.0491001; CrossSectionUncertainty = CrossSection* 13.6268/100.; }
  if (iStopMass==740) { CrossSection = 	0.0470801; CrossSectionUncertainty = CrossSection* 13.6895/100.; }
  if (iStopMass==745) { CrossSection = 	0.045061; CrossSectionUncertainty = CrossSection* 13.6816/100.; }
  if (iStopMass==750) { CrossSection = 	0.0431418; CrossSectionUncertainty = CrossSection* 13.7455/100.; }
  if (iStopMass==755) { CrossSection = 	0.0413447; CrossSectionUncertainty = CrossSection* 13.7833/100.; }
  if (iStopMass==760) { CrossSection = 	0.0396264; CrossSectionUncertainty = CrossSection* 13.8518/100.; }
  if (iStopMass==765) { CrossSection = 	0.0379036; CrossSectionUncertainty = CrossSection* 13.8537/100.; }
  if (iStopMass==770) { CrossSection = 	0.0363856; CrossSectionUncertainty = CrossSection* 13.9334/100.; }
  if (iStopMass==775) { CrossSection = 	0.0348796; CrossSectionUncertainty = CrossSection* 13.9597/100.; }
  if (iStopMass==780) { CrossSection = 	0.0334669; CrossSectionUncertainty = CrossSection* 14.0267/100.; }
  if (iStopMass==785) { CrossSection = 	0.0320548; CrossSectionUncertainty = CrossSection* 14.0406/100.; }
  if (iStopMass==790) { CrossSection = 	0.0307373; CrossSectionUncertainty = CrossSection* 14.115/100.; }
  if (iStopMass==795) { CrossSection = 	0.0295348; CrossSectionUncertainty = CrossSection* 14.1397/100.; }
  if (iStopMass==800) { CrossSection = 	0.0283338; CrossSectionUncertainty = CrossSection* 14.171/100.; }
  if (iStopMass==805) { CrossSection = 	0.0272206; CrossSectionUncertainty = CrossSection* 14.241/100.; }
  if (iStopMass==810) { CrossSection = 	0.0261233; CrossSectionUncertainty = CrossSection* 14.2891/100.; }
  if (iStopMass==815) { CrossSection = 	0.0251107; CrossSectionUncertainty = CrossSection* 14.3632/100.; }
  if (iStopMass==820) { CrossSection = 	0.0241099; CrossSectionUncertainty = CrossSection* 14.3805/100.; }
  if (iStopMass==825) { CrossSection = 	0.0230866; CrossSectionUncertainty = CrossSection* 14.4428/100.; }
  if (iStopMass==830) { CrossSection = 	0.0221834; CrossSectionUncertainty = CrossSection* 14.4791/100.; }
  if (iStopMass==835) { CrossSection = 	0.0213766; CrossSectionUncertainty = CrossSection* 14.5511/100.; }
  if (iStopMass==840) { CrossSection = 	0.0204715; CrossSectionUncertainty = CrossSection* 14.6131/100.; }
  if (iStopMass==845) { CrossSection = 	0.0197653; CrossSectionUncertainty = CrossSection* 14.6602/100.; }
  if (iStopMass==850) { CrossSection = 	0.0189612; CrossSectionUncertainty = CrossSection* 14.702/100.; }
  if (iStopMass==855) { CrossSection = 	0.0182516; CrossSectionUncertainty = CrossSection* 14.7648/100.; }
  if (iStopMass==860) { CrossSection = 	0.0175509; CrossSectionUncertainty = CrossSection* 14.7944/100.; }
  if (iStopMass==865) { CrossSection = 	0.0168336; CrossSectionUncertainty = CrossSection* 14.8528/100.; }
  if (iStopMass==870) { CrossSection = 	0.0162314; CrossSectionUncertainty = CrossSection* 14.8772/100.; }
  if (iStopMass==875) { CrossSection = 	0.015625; CrossSectionUncertainty = CrossSection* 14.9567/100.; }
  if (iStopMass==880) { CrossSection = 	0.0150143; CrossSectionUncertainty = CrossSection* 15.0389/100.; }
  if (iStopMass==885) { CrossSection = 	0.0144112; CrossSectionUncertainty = CrossSection* 15.0614/100.; }
  if (iStopMass==890) { CrossSection = 	0.0138979; CrossSectionUncertainty = CrossSection* 15.1/100.; }
  if (iStopMass==895) { CrossSection = 	0.0133962; CrossSectionUncertainty = CrossSection* 15.1325/100.; }
  if (iStopMass==900) { CrossSection = 	0.0128895; CrossSectionUncertainty = CrossSection* 15.2026/100.; }
  if (iStopMass==905) { CrossSection = 	0.0123843; CrossSectionUncertainty = CrossSection* 15.2968/100.; }
  if (iStopMass==910) { CrossSection = 	0.0119837; CrossSectionUncertainty = CrossSection* 15.3089/100.; }
  if (iStopMass==915) { CrossSection = 	0.0114713; CrossSectionUncertainty = CrossSection* 15.3678/100.; }
  if (iStopMass==920) { CrossSection = 	0.0110688; CrossSectionUncertainty = CrossSection* 15.4082/100.; }
  if (iStopMass==925) { CrossSection = 	0.0106631; CrossSectionUncertainty = CrossSection* 15.4806/100.; }
  if (iStopMass==930) { CrossSection = 	0.0102629; CrossSectionUncertainty = CrossSection* 15.5313/100.; }
  if (iStopMass==935) { CrossSection = 	0.0098874; CrossSectionUncertainty = CrossSection* 15.6066/100.; }
  if (iStopMass==940) { CrossSection = 	0.00952142; CrossSectionUncertainty = CrossSection* 15.6055/100.; }
  if (iStopMass==945) { CrossSection = 	0.00916636; CrossSectionUncertainty = CrossSection* 15.6849/100.; }
  if (iStopMass==950) { CrossSection = 	0.00883465; CrossSectionUncertainty = CrossSection* 15.7177/100.; }
  if (iStopMass==955) { CrossSection = 	0.00851073; CrossSectionUncertainty = CrossSection* 15.8094/100.; }
  if (iStopMass==960) { CrossSection = 	0.00820884; CrossSectionUncertainty = CrossSection* 15.844/100.; }
  if (iStopMass==965) { CrossSection = 	0.00791403; CrossSectionUncertainty = CrossSection* 15.9216/100.; }
  if (iStopMass==970) { CrossSection = 	0.00763112; CrossSectionUncertainty = CrossSection* 15.9742/100.; }
  if (iStopMass==975) { CrossSection = 	0.00735655; CrossSectionUncertainty = CrossSection* 16.0548/100.; }
  if (iStopMass==980) { CrossSection = 	0.00710317; CrossSectionUncertainty = CrossSection* 16.0626/100.; }
  if (iStopMass==985) { CrossSection = 	0.00684867; CrossSectionUncertainty = CrossSection* 16.144/100.; }
  if (iStopMass==990) { CrossSection = 	0.00660695; CrossSectionUncertainty = CrossSection* 16.1813/100.; }
  if (iStopMass==995) { CrossSection = 	0.00637546; CrossSectionUncertainty = CrossSection* 16.2158/100.; }
  if (iStopMass==1000) { CrossSection = 	0.00615134; CrossSectionUncertainty = CrossSection* 16.2953/100.; }
  if (iStopMass==1005) { CrossSection = 	0.00593765; CrossSectionUncertainty = CrossSection* 16.3716/100.; }
  if (iStopMass==1010) { CrossSection = 	0.00572452; CrossSectionUncertainty = CrossSection* 16.3857/100.; }
  if (iStopMass==1015) { CrossSection = 	0.00553094; CrossSectionUncertainty = CrossSection* 16.4628/100.; }
  if (iStopMass==1020) { CrossSection = 	0.00533968; CrossSectionUncertainty = CrossSection* 16.4963/100.; }
  if (iStopMass==1025) { CrossSection = 	0.00514619; CrossSectionUncertainty = CrossSection* 16.5762/100.; }
  if (iStopMass==1030) { CrossSection = 	0.00497235; CrossSectionUncertainty = CrossSection* 16.5838/100.; }
  if (iStopMass==1035) { CrossSection = 	0.00479906; CrossSectionUncertainty = CrossSection* 16.6646/100.; }
  if (iStopMass==1040) { CrossSection = 	0.00463806; CrossSectionUncertainty = CrossSection* 16.6947/100.; }
  if (iStopMass==1045) { CrossSection = 	0.00447537; CrossSectionUncertainty = CrossSection* 16.7071/100.; }
  if (iStopMass==1050) { CrossSection = 	0.00432261; CrossSectionUncertainty = CrossSection* 16.7859/100.; }
  if (iStopMass==1055) { CrossSection = 	0.00417983; CrossSectionUncertainty = CrossSection* 16.8637/100.; }
  if (iStopMass==1060) { CrossSection = 	0.00403886; CrossSectionUncertainty = CrossSection* 16.8981/100.; }
  if (iStopMass==1065) { CrossSection = 	0.0038962; CrossSectionUncertainty = CrossSection* 16.9794/100.; }
  if (iStopMass==1070) { CrossSection = 	0.00376343; CrossSectionUncertainty = CrossSection* 16.9764/100.; }
  if (iStopMass==1075) { CrossSection = 	0.00364174; CrossSectionUncertainty = CrossSection* 17.0634/100.; }
  if (iStopMass==1080) { CrossSection = 	0.00352093; CrossSectionUncertainty = CrossSection* 17.0908/100.; }
  if (iStopMass==1085) { CrossSection = 	0.00339813; CrossSectionUncertainty = CrossSection* 17.1929/100.; }
  if (iStopMass==1090) { CrossSection = 	0.00328695; CrossSectionUncertainty = CrossSection* 17.2274/100.; }
  if (iStopMass==1095) { CrossSection = 	0.00317628; CrossSectionUncertainty = CrossSection* 17.2617/100.; }
  if (iStopMass==1100) { CrossSection = 	0.00307413; CrossSectionUncertainty = CrossSection* 17.3377/100.; }
  if (iStopMass==1105) { CrossSection = 	0.00297377; CrossSectionUncertainty = CrossSection* 17.3822/100.; }
  if (iStopMass==1110) { CrossSection = 	0.00287148; CrossSectionUncertainty = CrossSection* 17.4725/100.; }
  if (iStopMass==1115) { CrossSection = 	0.00278078; CrossSectionUncertainty = CrossSection* 17.5091/100.; }
  if (iStopMass==1120) { CrossSection = 	0.00268873; CrossSectionUncertainty = CrossSection* 17.5883/100.; }
  if (iStopMass==1125) { CrossSection = 	0.00260821; CrossSectionUncertainty = CrossSection* 17.6126/100.; }
  if (iStopMass==1130) { CrossSection = 	0.00251529; CrossSectionUncertainty = CrossSection* 17.6836/100.; }
  if (iStopMass==1135) { CrossSection = 	0.00243484; CrossSectionUncertainty = CrossSection* 17.7128/100.; }
  if (iStopMass==1140) { CrossSection = 	0.00236295; CrossSectionUncertainty = CrossSection* 17.7977/100.; }
  if (iStopMass==1145) { CrossSection = 	0.00228192; CrossSectionUncertainty = CrossSection* 17.8507/100.; }
  if (iStopMass==1150) { CrossSection = 	0.00221047; CrossSectionUncertainty = CrossSection* 17.9259/100.; }
  if (iStopMass==1155) { CrossSection = 	0.00213907; CrossSectionUncertainty = CrossSection* 18.0255/100.; }
  if (iStopMass==1160) { CrossSection = 	0.00206845; CrossSectionUncertainty = CrossSection* 18.0518/100.; }
  if (iStopMass==1165) { CrossSection = 	0.0020063; CrossSectionUncertainty = CrossSection* 18.0954/100.; }
  if (iStopMass==1170) { CrossSection = 	0.00194569; CrossSectionUncertainty = CrossSection* 18.1194/100.; }
  if (iStopMass==1175) { CrossSection = 	0.0018741; CrossSectionUncertainty = CrossSection* 18.2145/100.; }
  if (iStopMass==1180) { CrossSection = 	0.00182266; CrossSectionUncertainty = CrossSection* 18.3074/100.; }
  if (iStopMass==1185) { CrossSection = 	0.00176211; CrossSectionUncertainty = CrossSection* 18.3375/100.; }
  if (iStopMass==1190) { CrossSection = 	0.00170006; CrossSectionUncertainty = CrossSection* 18.4075/100.; }
  if (iStopMass==1195) { CrossSection = 	0.00164968; CrossSectionUncertainty = CrossSection* 18.4438/100.; }
  if (iStopMass==1200) { CrossSection = 	0.00159844; CrossSectionUncertainty = CrossSection* 18.5209/100.; }
  if (iStopMass==1205) { CrossSection = 	0.0015472; CrossSectionUncertainty = CrossSection* 18.5977/100.; }
  if (iStopMass==1210) { CrossSection = 	0.00149657; CrossSectionUncertainty = CrossSection* 18.6485/100.; }
  if (iStopMass==1215) { CrossSection = 	0.00145544; CrossSectionUncertainty = CrossSection* 18.7347/100.; }
  if (iStopMass==1220) { CrossSection = 	0.00140288; CrossSectionUncertainty = CrossSection* 18.8774/100.; }
  if (iStopMass==1225) { CrossSection = 	0.00136155; CrossSectionUncertainty = CrossSection* 18.989/100.; }
  if (iStopMass==1230) { CrossSection = 	0.00131271; CrossSectionUncertainty = CrossSection* 18.8763/100.; }
  if (iStopMass==1235) { CrossSection = 	0.0012717; CrossSectionUncertainty = CrossSection* 18.9588/100.; }
  if (iStopMass==1240) { CrossSection = 	0.00123066; CrossSectionUncertainty = CrossSection* 19.049/100.; }
  if (iStopMass==1245) { CrossSection = 	0.00119994; CrossSectionUncertainty = CrossSection* 19.1442/100.; }
  if (iStopMass==1250) { CrossSection = 	0.0011583; CrossSectionUncertainty = CrossSection* 19.3006/100.; }
  if (iStopMass==1255) { CrossSection = 	0.00112694; CrossSectionUncertainty = CrossSection* 19.4441/100.; }
  if (iStopMass==1260) { CrossSection = 	0.00108716; CrossSectionUncertainty = CrossSection* 19.4141/100.; }
  if (iStopMass==1265) { CrossSection = 	0.00105517; CrossSectionUncertainty = CrossSection* 19.6361/100.; }
  if (iStopMass==1270) { CrossSection = 	0.00102241; CrossSectionUncertainty = CrossSection* 19.6297/100.; }
  if (iStopMass==1275) { CrossSection = 	0.000991293; CrossSectionUncertainty = CrossSection* 19.762/100.; }
  if (iStopMass==1280) { CrossSection = 	0.000961012; CrossSectionUncertainty = CrossSection* 19.7926/100.; }
  if (iStopMass==1285) { CrossSection = 	0.000932394; CrossSectionUncertainty = CrossSection* 19.8682/100.; }
  if (iStopMass==1290) { CrossSection = 	0.000903404; CrossSectionUncertainty = CrossSection* 19.9924/100.; }
  if (iStopMass==1295) { CrossSection = 	0.000876957; CrossSectionUncertainty = CrossSection* 20.0777/100.; }
  if (iStopMass==1300) { CrossSection = 	0.000850345; CrossSectionUncertainty = CrossSection* 20.1604/100.; }
  if (iStopMass==1305) { CrossSection = 	0.00082443; CrossSectionUncertainty = CrossSection* 20.2883/100.; }
  if (iStopMass==1310) { CrossSection = 	0.00079983; CrossSectionUncertainty = CrossSection* 20.373/100.; }
  if (iStopMass==1315) { CrossSection = 	0.000775222; CrossSectionUncertainty = CrossSection* 20.4622/100.; }
  if (iStopMass==1320) { CrossSection = 	0.000751372; CrossSectionUncertainty = CrossSection* 20.5919/100.; }
  if (iStopMass==1325) { CrossSection = 	0.000728912; CrossSectionUncertainty = CrossSection* 20.6884/100.; }
  if (iStopMass==1330) { CrossSection = 	0.000706867; CrossSectionUncertainty = CrossSection* 20.7763/100.; }
  if (iStopMass==1335) { CrossSection = 	0.000685372; CrossSectionUncertainty = CrossSection* 20.8587/100.; }
  if (iStopMass==1340) { CrossSection = 	0.000664649; CrossSectionUncertainty = CrossSection* 20.9879/100.; }
  if (iStopMass==1345) { CrossSection = 	0.000644804; CrossSectionUncertainty = CrossSection* 21.1487/100.; }
  if (iStopMass==1350) { CrossSection = 	0.000625155; CrossSectionUncertainty = CrossSection* 21.2761/100.; }
  if (iStopMass==1355) { CrossSection = 	0.000606802; CrossSectionUncertainty = CrossSection* 21.3529/100.; }
  if (iStopMass==1360) { CrossSection = 	0.000588512; CrossSectionUncertainty = CrossSection* 21.4428/100.; }
  if (iStopMass==1365) { CrossSection = 	0.000570506; CrossSectionUncertainty = CrossSection* 21.6584/100.; }
  if (iStopMass==1370) { CrossSection = 	0.000553379; CrossSectionUncertainty = CrossSection* 21.6036/100.; }
  if (iStopMass==1375) { CrossSection = 	0.000536646; CrossSectionUncertainty = CrossSection* 21.775/100.; }
  if (iStopMass==1380) { CrossSection = 	0.000521404; CrossSectionUncertainty = CrossSection* 21.8383/100.; }
  if (iStopMass==1385) { CrossSection = 	0.000505008; CrossSectionUncertainty = CrossSection* 21.9675/100.; }
  if (iStopMass==1390) { CrossSection = 	0.000490353; CrossSectionUncertainty = CrossSection* 22.1444/100.; }
  if (iStopMass==1395) { CrossSection = 	0.000476164; CrossSectionUncertainty = CrossSection* 22.2016/100.; }
  if (iStopMass==1400) { CrossSection = 	0.000461944; CrossSectionUncertainty = CrossSection* 22.2704/100.; }
  if (iStopMass==1405) { CrossSection = 	0.000448172; CrossSectionUncertainty = CrossSection* 22.4911/100.; }
  if (iStopMass==1410) { CrossSection = 	0.000435082; CrossSectionUncertainty = CrossSection* 22.5606/100.; }
  if (iStopMass==1415) { CrossSection = 	0.000422967; CrossSectionUncertainty = CrossSection* 22.6095/100.; }
  if (iStopMass==1420) { CrossSection = 	0.000410381; CrossSectionUncertainty = CrossSection* 22.797/100.; }
  if (iStopMass==1425) { CrossSection = 	0.000398106; CrossSectionUncertainty = CrossSection* 22.8949/100.; }
  if (iStopMass==1430) { CrossSection = 	0.000386792; CrossSectionUncertainty = CrossSection* 23.1319/100.; }
  if (iStopMass==1435) { CrossSection = 	0.000375724; CrossSectionUncertainty = CrossSection* 23.1724/100.; }
  if (iStopMass==1440) { CrossSection = 	0.000364616; CrossSectionUncertainty = CrossSection* 23.2234/100.; }
  if (iStopMass==1445) { CrossSection = 	0.000353965; CrossSectionUncertainty = CrossSection* 23.4637/100.; }
  if (iStopMass==1450) { CrossSection = 	0.000343923; CrossSectionUncertainty = CrossSection* 23.4948/100.; }
  if (iStopMass==1455) { CrossSection = 	0.000333885; CrossSectionUncertainty = CrossSection* 23.5468/100.; }
  if (iStopMass==1460) { CrossSection = 	0.000324344; CrossSectionUncertainty = CrossSection* 23.771/100.; }
  if (iStopMass==1465) { CrossSection = 	0.0003153; CrossSectionUncertainty = CrossSection* 23.8004/100.; }
  if (iStopMass==1470) { CrossSection = 	0.00030583; CrossSectionUncertainty = CrossSection* 24.0064/100.; }
  if (iStopMass==1475) { CrossSection = 	0.000296811; CrossSectionUncertainty = CrossSection* 24.0314/100.; }
  if (iStopMass==1480) { CrossSection = 	0.000288149; CrossSectionUncertainty = CrossSection* 23.9248/100.; }
  if (iStopMass==1485) { CrossSection = 	0.000279711; CrossSectionUncertainty = CrossSection* 24.1257/100.; }
  if (iStopMass==1490) { CrossSection = 	0.000271724; CrossSectionUncertainty = CrossSection* 24.1274/100.; }
  if (iStopMass==1495) { CrossSection = 	0.000264275; CrossSectionUncertainty = CrossSection* 24.3545/100.; }
  if (iStopMass==1500) { CrossSection = 	0.000256248; CrossSectionUncertainty = CrossSection* 24.372/100.; }
  if (iStopMass==1505) { CrossSection = 	0.000248853; CrossSectionUncertainty = CrossSection* 24.5827/100.; }
  if (iStopMass==1510) { CrossSection = 	0.000241844; CrossSectionUncertainty = CrossSection* 24.6187/100.; }
  if (iStopMass==1515) { CrossSection = 	0.000234438; CrossSectionUncertainty = CrossSection* 24.8442/100.; }
  if (iStopMass==1520) { CrossSection = 	0.000227374; CrossSectionUncertainty = CrossSection* 24.8909/100.; }
  if (iStopMass==1525) { CrossSection = 	0.000221045; CrossSectionUncertainty = CrossSection* 25.0895/100.; }
  if (iStopMass==1530) { CrossSection = 	0.000214431; CrossSectionUncertainty = CrossSection* 24.8728/100.; }
  if (iStopMass==1535) { CrossSection = 	0.000208092; CrossSectionUncertainty = CrossSection* 25.1043/100.; }
  if (iStopMass==1540) { CrossSection = 	0.000201748; CrossSectionUncertainty = CrossSection* 25.3207/100.; }
  if (iStopMass==1545) { CrossSection = 	0.000196399; CrossSectionUncertainty = CrossSection* 25.5641/100.; }
  if (iStopMass==1550) { CrossSection = 	0.000190474; CrossSectionUncertainty = CrossSection* 25.5213/100.; }
  if (iStopMass==1555) { CrossSection = 	0.000185188; CrossSectionUncertainty = CrossSection* 25.7329/100.; }
  if (iStopMass==1560) { CrossSection = 	0.000179263; CrossSectionUncertainty = CrossSection* 25.6931/100.; }
  if (iStopMass==1565) { CrossSection = 	0.000174021; CrossSectionUncertainty = CrossSection* 25.9111/100.; }
  if (iStopMass==1570) { CrossSection = 	0.000169176; CrossSectionUncertainty = CrossSection* 25.8106/100.; }
  if (iStopMass==1575) { CrossSection = 	0.000163861; CrossSectionUncertainty = CrossSection* 26.0597/100.; }
  if (iStopMass==1580) { CrossSection = 	0.000159583; CrossSectionUncertainty = CrossSection* 26.2958/100.; }
  if (iStopMass==1585) { CrossSection = 	0.000154719; CrossSectionUncertainty = CrossSection* 26.195/100.; }
  if (iStopMass==1590) { CrossSection = 	0.000150506; CrossSectionUncertainty = CrossSection* 26.4111/100.; }
  if (iStopMass==1595) { CrossSection = 	0.000145626; CrossSectionUncertainty = CrossSection* 26.3077/100.; }
  if (iStopMass==1600) { CrossSection = 	0.000141382; CrossSectionUncertainty = CrossSection* 26.5291/100.; }
  if (iStopMass==1605) { CrossSection = 	0.000137131; CrossSectionUncertainty = CrossSection* 26.7424/100.; }
  if (iStopMass==1610) { CrossSection = 	0.000132187; CrossSectionUncertainty = CrossSection* 26.668/100.; }
  if (iStopMass==1615) { CrossSection = 	0.000127929; CrossSectionUncertainty = CrossSection* 26.9117/100.; }
  if (iStopMass==1620) { CrossSection = 	0.000124086; CrossSectionUncertainty = CrossSection* 26.7738/100.; }
  if (iStopMass==1625) { CrossSection = 	0.00011982; CrossSectionUncertainty = CrossSection* 27.0483/100.; }
  if (iStopMass==1630) { CrossSection = 	0.000116042; CrossSectionUncertainty = CrossSection* 26.8071/100.; }
  if (iStopMass==1635) { CrossSection = 	0.000112767; CrossSectionUncertainty = CrossSection* 27.127/100.; }
  if (iStopMass==1640) { CrossSection = 	0.000108936; CrossSectionUncertainty = CrossSection* 26.9351/100.; }
  if (iStopMass==1645) { CrossSection = 	0.000105746; CrossSectionUncertainty = CrossSection* 27.1783/100.; }
  if (iStopMass==1650) { CrossSection = 	0.000102693; CrossSectionUncertainty = CrossSection* 27.292/100.; }
  if (iStopMass==1655) { CrossSection = 	0.000100112; CrossSectionUncertainty = CrossSection* 27.4445/100.; }
  if (iStopMass==1660) { CrossSection = 	9.75763e-05; CrossSectionUncertainty = CrossSection* 27.5431/100.; }
  if (iStopMass==1665) { CrossSection = 	9.52062e-05; CrossSectionUncertainty = CrossSection* 27.6946/100.; }
  if (iStopMass==1670) { CrossSection = 	9.29857e-05; CrossSectionUncertainty = CrossSection* 27.7869/100.; }
  if (iStopMass==1675) { CrossSection = 	9.08285e-05; CrossSectionUncertainty = CrossSection* 27.9347/100.; }
  if (iStopMass==1680) { CrossSection = 	8.87433e-05; CrossSectionUncertainty = CrossSection* 28.1539/100.; }
  if (iStopMass==1685) { CrossSection = 	8.66618e-05; CrossSectionUncertainty = CrossSection* 28.3509/100.; }
  if (iStopMass==1690) { CrossSection = 	8.46535e-05; CrossSectionUncertainty = CrossSection* 28.4432/100.; }
  if (iStopMass==1695) { CrossSection = 	8.27102e-05; CrossSectionUncertainty = CrossSection* 28.591/100.; }
  if (iStopMass==1700) { CrossSection = 	8.07774e-05; CrossSectionUncertainty = CrossSection* 28.7497/100.; }
  if (iStopMass==1705) { CrossSection = 	7.8666e-05; CrossSectionUncertainty = CrossSection* 28.8194/100.; }
  if (iStopMass==1710) { CrossSection = 	7.6572e-05; CrossSectionUncertainty = CrossSection* 29.0265/100.; }
  if (iStopMass==1715) { CrossSection = 	7.45994e-05; CrossSectionUncertainty = CrossSection* 29.1193/100.; }
  if (iStopMass==1720) { CrossSection = 	7.25199e-05; CrossSectionUncertainty = CrossSection* 29.3013/100.; }
  if (iStopMass==1725) { CrossSection = 	7.05189e-05; CrossSectionUncertainty = CrossSection* 29.3697/100.; }
  if (iStopMass==1730) { CrossSection = 	6.85712e-05; CrossSectionUncertainty = CrossSection* 29.4972/100.; }
  if (iStopMass==1735) { CrossSection = 	6.67296e-05; CrossSectionUncertainty = CrossSection* 29.6167/100.; }
  if (iStopMass==1740) { CrossSection = 	6.49184e-05; CrossSectionUncertainty = CrossSection* 29.7686/100.; }
  if (iStopMass==1745) { CrossSection = 	6.30949e-05; CrossSectionUncertainty = CrossSection* 29.8524/100.; }
  if (iStopMass==1750) { CrossSection = 	6.13637e-05; CrossSectionUncertainty = CrossSection* 29.9789/100.; }
  if (iStopMass==1755) { CrossSection = 	5.97301e-05; CrossSectionUncertainty = CrossSection* 30.0928/100.; }
  if (iStopMass==1760) { CrossSection = 	5.80751e-05; CrossSectionUncertainty = CrossSection* 30.2585/100.; }
  if (iStopMass==1765) { CrossSection = 	5.65479e-05; CrossSectionUncertainty = CrossSection* 30.366/100.; }
  if (iStopMass==1770) { CrossSection = 	5.49998e-05; CrossSectionUncertainty = CrossSection* 30.5241/100.; }
  if (iStopMass==1775) { CrossSection = 	5.35686e-05; CrossSectionUncertainty = CrossSection* 30.6718/100.; }
  if (iStopMass==1780) { CrossSection = 	5.20828e-05; CrossSectionUncertainty = CrossSection* 30.6799/100.; }
  if (iStopMass==1785) { CrossSection = 	5.07079e-05; CrossSectionUncertainty = CrossSection* 30.9201/100.; }
  if (iStopMass==1790) { CrossSection = 	4.93948e-05; CrossSectionUncertainty = CrossSection* 31.0043/100.; }
  if (iStopMass==1795) { CrossSection = 	4.80635e-05; CrossSectionUncertainty = CrossSection* 31.138/100.; }
  if (iStopMass==1800) { CrossSection = 	4.67492e-05; CrossSectionUncertainty = CrossSection* 31.2291/100.; }
  if (iStopMass==1805) { CrossSection = 	4.55055e-05; CrossSectionUncertainty = CrossSection* 31.4321/100.; }
  if (iStopMass==1810) { CrossSection = 	4.42835e-05; CrossSectionUncertainty = CrossSection* 31.5499/100.; }
  if (iStopMass==1815) { CrossSection = 	4.30744e-05; CrossSectionUncertainty = CrossSection* 31.6302/100.; }
  if (iStopMass==1820) { CrossSection = 	4.19954e-05; CrossSectionUncertainty = CrossSection* 31.7151/100.; }
  if (iStopMass==1825) { CrossSection = 	4.08527e-05; CrossSectionUncertainty = CrossSection* 31.9048/100.; }
  if (iStopMass==1830) { CrossSection = 	3.97561e-05; CrossSectionUncertainty = CrossSection* 31.9718/100.; }
  if (iStopMass==1835) { CrossSection = 	3.87041e-05; CrossSectionUncertainty = CrossSection* 32.2028/100.; }
  if (iStopMass==1840) { CrossSection = 	3.76008e-05; CrossSectionUncertainty = CrossSection* 32.268/100.; }
  if (iStopMass==1845) { CrossSection = 	3.66914e-05; CrossSectionUncertainty = CrossSection* 32.4529/100.; }
  if (iStopMass==1850) { CrossSection = 	3.56995e-05; CrossSectionUncertainty = CrossSection* 32.5039/100.; }
  if (iStopMass==1855) { CrossSection = 	3.47689e-05; CrossSectionUncertainty = CrossSection* 32.6767/100.; }
  if (iStopMass==1860) { CrossSection = 	3.38528e-05; CrossSectionUncertainty = CrossSection* 32.8878/100.; }
  if (iStopMass==1865) { CrossSection = 	3.29644e-05; CrossSectionUncertainty = CrossSection* 32.8975/100.; }
  if (iStopMass==1870) { CrossSection = 	3.20679e-05; CrossSectionUncertainty = CrossSection* 32.9608/100.; }
  if (iStopMass==1875) { CrossSection = 	3.12583e-05; CrossSectionUncertainty = CrossSection* 33.1541/100.; }
  if (iStopMass==1880) { CrossSection = 	3.04342e-05; CrossSectionUncertainty = CrossSection* 33.3117/100.; }
  if (iStopMass==1885) { CrossSection = 	2.96516e-05; CrossSectionUncertainty = CrossSection* 33.2866/100.; }
  if (iStopMass==1890) { CrossSection = 	2.88952e-05; CrossSectionUncertainty = CrossSection* 33.6279/100.; }
  if (iStopMass==1895) { CrossSection = 	2.81145e-05; CrossSectionUncertainty = CrossSection* 33.6845/100.; }
  if (iStopMass==1900) { CrossSection = 	2.73974e-05; CrossSectionUncertainty = CrossSection* 33.8247/100.; }
  if (iStopMass==1905) { CrossSection = 	2.66796e-05; CrossSectionUncertainty = CrossSection* 33.9708/100.; }
  if (iStopMass==1910) { CrossSection = 	2.59941e-05; CrossSectionUncertainty = CrossSection* 33.9526/100.; }
  if (iStopMass==1915) { CrossSection = 	2.52784e-05; CrossSectionUncertainty = CrossSection* 34.1137/100.; }
  if (iStopMass==1920) { CrossSection = 	2.46598e-05; CrossSectionUncertainty = CrossSection* 34.2714/100.; }
  if (iStopMass==1925) { CrossSection = 	2.39932e-05; CrossSectionUncertainty = CrossSection* 34.2328/100.; }
  if (iStopMass==1930) { CrossSection = 	2.33737e-05; CrossSectionUncertainty = CrossSection* 34.394/100.; }
  if (iStopMass==1935) { CrossSection = 	2.27623e-05; CrossSectionUncertainty = CrossSection* 34.5138/100.; }
  if (iStopMass==1940) { CrossSection = 	2.21454e-05; CrossSectionUncertainty = CrossSection* 34.6933/100.; }
  if (iStopMass==1945) { CrossSection = 	2.15924e-05; CrossSectionUncertainty = CrossSection* 35.0815/100.; }
  if (iStopMass==1950) { CrossSection = 	2.10232e-05; CrossSectionUncertainty = CrossSection* 34.9444/100.; }
  if (iStopMass==1955) { CrossSection = 	2.05211e-05; CrossSectionUncertainty = CrossSection* 35.0155/100.; }
  if (iStopMass==1960) { CrossSection = 	1.98996e-05; CrossSectionUncertainty = CrossSection* 35.2135/100.; }
  if (iStopMass==1965) { CrossSection = 	1.9408e-05; CrossSectionUncertainty = CrossSection* 35.3328/100.; }
  if (iStopMass==1970) { CrossSection = 	1.88974e-05; CrossSectionUncertainty = CrossSection* 35.4643/100.; }
  if (iStopMass==1975) { CrossSection = 	1.84612e-05; CrossSectionUncertainty = CrossSection* 35.7904/100.; }
  if (iStopMass==1980) { CrossSection = 	1.79562e-05; CrossSectionUncertainty = CrossSection* 35.8898/100.; }
  if (iStopMass==1985) { CrossSection = 	1.75673e-05; CrossSectionUncertainty = CrossSection* 35.989/100.; }
  if (iStopMass==1990) { CrossSection = 	1.70612e-05; CrossSectionUncertainty = CrossSection* 36.0953/100.; }
  if (iStopMass==1995) { CrossSection = 	1.66228e-05; CrossSectionUncertainty = CrossSection* 36.4709/100.; }
  if (iStopMass==2000) { CrossSection = 	1.62355e-05; CrossSectionUncertainty = CrossSection* 36.5277/100.; }
 
} 

// https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf

float ISRBinWeight[7] = {1., 0.920, 0.821, 0.715, 0.662, 0.561, 0.511};

bool IsISRJet(float pt, float eta, float phi)
{

  // https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295

  // https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L373-L395
  // https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/interface/jet_met_tools.hh#L34-L36
  if (pt<=30. || fabs(eta)>2.4) return false;
  
  TLorentzVector ThisJet; ThisJet.SetPtEtaPhiM(pt, eta, phi, 0.);

  for (uint i = 0; i<std_vector_partonGen_pt->size(); i++) {
    
    if (fabs(std_vector_partonGen_pid->at(i))!=5) continue;
    if (std_vector_partonGen_isHardProcess->at(i)!=1) continue;
    
    TLorentzVector bQuark; bQuark.SetPtEtaPhiM(std_vector_partonGen_pt ->at(i),
					       std_vector_partonGen_eta->at(i),
					       std_vector_partonGen_phi->at(i),
					       4.18);
    
    for (int wb = std_vector_VBoson_pt->size()-1; wb>=0; wb--) {
      
      if (std_vector_VBoson_pt->at(wb) <= 0.) continue;
      
      int Wid  = std_vector_VBoson_pid->at(wb);  // Wid  = -24 for W- and +24 for W+
    
      if ((Wid==+24 && std_vector_partonGen_pid->at(i)==+5) || 
	  (Wid==-24 && std_vector_partonGen_pid->at(i)==-5)) {
	cout << "got pair"<< endl;
	TLorentzVector WBoson; WBoson.SetPtEtaPhiM(std_vector_VBoson_pt ->at(wb),
						   std_vector_VBoson_eta->at(wb),
						   std_vector_VBoson_phi->at(wb),
						   80.385);

	TLorentzVector tCandidate = bQuark + WBoson;

	for (uint j = 0; j<std_vector_partonGen_pt->size(); j++) {

	  if (fabs(std_vector_partonGen_pid->at(j))!=6) continue;
	  if (std_vector_partonGen_isHardProcess->at(j)!=1) continue;
	  if (std_vector_partonGen_pid->at(i)*std_vector_partonGen_pid->at(j)<0) continue;
	  cout <<"got top"<<endl;
	  TLorentzVector tQuark; tQuark.SetPtEtaPhiM(std_vector_partonGen_pt ->at(j),
						     std_vector_partonGen_eta->at(j),
						     std_vector_partonGen_phi->at(j),
						     172.44);
	  cout << "   " << tQuark.DeltaR(tCandidate) << " " << ThisJet.DeltaR(bQuark) << endl;
	  if (tQuark.DeltaR(tCandidate)<0.1) 
	    if (ThisJet.DeltaR(bQuark)<0.3) {
	      cout << "Got b"<<endl;
	      return false;
	    }
	}

      }

    }

  }	  

  return true;

}

void CreateMap(bool IsTest = false, int DataRange = -1) {

  //TString TreeDirectory = "/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoSkims/Feb2017_Summer16_stop_ghent_others/MCl2stopnosel/";
  //TString TreeDirectory = "/eos/user/s/scodella/Stop/LatinoSkims/Feb2017_Summer16_stop_ghent_others_isr/MCl2stop/";
  //TString TreeDirectory = "/eos/user/s/scodella/Stop/LatinoSkims/Feb2017_Summer16_stop_ghent_others/MCl2stop__BlindSFWeights__hadd/";
  TString TreeDirectory = "../minitrees/nominal/Stop/";
  const int nMassScans = 8;
  /*TString MassScanName[nMassScans] = {"latino_T2tt_mStop-150to250__part",
				      "latino_T2tt_mStop-250to350__part",
				      "latino_T2tt_mStop-350to400__part",
				      "latino_T2tt_mStop-400to1200__part",
				      "latino_T2bW__part",
				      "latino_T2bt__part",
				      "latino_T2tt_dM10to80__part",
				      "latino_T2bW_dM10to80__part"};
  //int nMassScanFiles[nMassScans] = {45, 48, 43, 43,  0,  0,  0,  0};
  //int nMassScanFiles[nMassScans] = {45,  0,  0, 43,  0,  0,  0,  0};
  //int nMassScanFiles[nMassScans] = { 0,  0,  0,  0, 44,  0,  0,  0};
  //int nMassScanFiles[nMassScans] = { 0,  0,  0,  0,  0, 48,  0,  0};
  //int nMassScanFiles[nMassScans] = { 0,  0,  0,  0,  0,  0,  0,  0};
  //int nMassScanFiles[nMassScans] = { 0,  0,  0,  0,  0,  0,  9,  0};
  int nMassScanFiles[nMassScans] = { 0,  0,  0,  0,  0,  0,  0,  9};*/
  TString MassScanName[nMassScans] = {"T2tt_mStop-150to250",
				      "T2tt_mStop-250to350",
				      "T2tt_mStop-350to400",
				      "T2tt_mStop-400to1200",
				      "T2bW_isr",
				      "T2bt",
				      "T2tt_dM10to80",
				      "T2bW_dM10to80"};
  int nMassScanFiles[nMassScans] = { 0,  0,  0,  0,  0,  0,  0,  1};

  TString TreeType = MassScanName[0].Contains("latino") ? "latino" : "minitrees";

  TH2F *MassPlan = new TH2F("MassPlan", "", 1300, 0., 1300., 1300, 0., 1300.);

  int nMPfound = 0;
  int stopMP[1000], chi0MP[1000];
  TH1F *nISRJets[1000];
  for (int mp = 0; mp<1000; mp++) {
    
    TString MP = "MP"; MP += mp;
    stopMP[mp] = chi0MP[mp] = -1;
    nISRJets[mp] = new TH1F(MP, "", 7, 0, 7);
      
  }

  int FirstStopMass = -1, FirstNeutralinoMass = -1;
  
  if (DataRange>=-1) {

    for (int ms = 0; ms<nMassScans; ms++) 
      for (int tf = 0; tf<nMassScanFiles[ms]; tf++) {

	if (DataRange==-1 || tf==DataRange) {
	  
	  TString FileName = TreeDirectory + MassScanName[ms]; 
	  if (FileName.Contains("__part")) FileName += tf;
	  TFile *TreeFile = TFile::Open(FileName + ".root");
      
	  if (FileName.Contains("dM10to80")) Apply2LFilterEff = true;
	  
	  TTree *Tree = GetTree(TreeFile, TreeType);
	  
	  Int_t nentries = IsTest ? 100 : (Int_t) Tree->GetEntries();
	  
	  cout << "Looping on " << FileName << endl;
	  
	  bool IsDone[10] = {false};
	  
	  for (Int_t i = 0; i<nentries; i++) {
	    
	    if (100*(i+1)/nentries%10==0)
	      if (!IsDone[10*(i)/nentries]) {
		cout << "  Completed " << 100*i/nentries << "% of the events" << endl; 
		IsDone[10*(i)/nentries] = true;
	      }
	    
	    Tree->GetEntry(i);
	    
	    if (i==0) { FirstStopMass = susyMstop; FirstNeutralinoMass = susyMLSP; }

	    MassPlan->Fill(susyMstop, susyMLSP);

	    int ThisMP = -1;
	    for (int mp = 0; mp<nMPfound; mp++) 
	      if (stopMP[mp]==int(susyMstop) && chi0MP[mp]==int(susyMLSP)) 
		ThisMP = mp;
	    if (ThisMP==-1) {
	      ThisMP = nMPfound;
	      stopMP[nMPfound] = int(susyMstop);
	      chi0MP[nMPfound] = int(susyMLSP);
	      nMPfound++;
	    }
	    /*
	    int nLeptonFromW = 0;
	    for (uint ll = 0; ll<std_vector_leptonGen_MotherPID->size(); ll++)
	      if (fabs(std_vector_leptonGen_MotherPID->at(ll))==24) nLeptonFromW++;
	    if (nLeptonFromW<2) continue;
	    
	    int nTightLeptons = 0;
	    for (uint ll = 0; ll<std_vector_lepton_isTightLepton->size(); ll++)
	      if (std_vector_lepton_isTightLepton->at(ll)==1) nTightLeptons++;
	    if (nTightLeptons<2) continue;
	    */
	    int _nisrjets = 0;
	    float EventWeight = 1.;

	    if (TreeType=="minitrees") {

	      _nisrjets = nisrjet;
	      EventWeight = eventW;
	      
	    } else {
	      
	      for (uint j=0; j<std_vector_jet_pt->size(); j++) {
		
		float pt  = std_vector_jet_pt ->at(j); if (pt<0.) continue;
		float eta = std_vector_jet_eta->at(j);
		//float phi = std_vector_jet_phi->at(j);
		
		if (pt<=30. || fabs(eta)>2.4) continue;
		
		if (std_vector_jet_isFromISR->at(j)==1) _nisrjets++;
		//if (IsISRJet(pt, eta, phi)) _nisrjets++;
		
	      }

	    }
	      
	    if (_nisrjets>6) _nisrjets = 6;
	    
	    nISRJets[ThisMP]->Fill(_nisrjets, EventWeight);
	    
	  }

	}

      }

  }
  
  if (DataRange>=0) {
    
    TString HistoFileName = "./MassPlanHistograms/MassPlan_part"; HistoFileName += DataRange;
    TFile *HistoFile = new TFile(HistoFileName + ".root", "recreate");
   
    MassPlan->Write();

    for (int mp = 0; mp<nMPfound; mp++) {

      TString HistoName = "nISRJets_mStop-"; HistoName += stopMP[mp];
      HistoName += "_mLSP-"; HistoName += chi0MP[mp];
      nISRJets[mp]->SetName(HistoName);
      nISRJets[mp]->Write();

    }

    HistoFile->Close();
 

  } else {
    
    if (DataRange==-2) {

      /*      
      for (int ms = 0; ms<nMassScans; ms++) {
	for (int tf = 0; tf<nMassScanFiles[ms]; tf++) {
	  
	  TString FileName = "./MassPlanHistograms/MassPlan_T2tt__part"; FileName += tf; FileName += ".root";
	  TFile *HistoFile = TFile::Open(FileName);
    
	  TH2F *ThisMassPlan = (TH2F*) HistoFile->Get("MassPlan");
	  MassPlan->Add(ThisMassPlan);

	}

	cout << MassScanName[ms] << " " << MassPlan->GetEntries() << endl;

      }
      */

      TString FileName = "./MassPlanHistograms/MassPlan_T2bW_dM10to80.root";
      TFile *HistoFile = TFile::Open(FileName);
      
      if (FileName.Contains("dM10to80")) Apply2LFilterEff = true;

      TH2F *ThisMassPlan = (TH2F*) HistoFile->Get("MassPlan");
      MassPlan->Add(ThisMassPlan);

      for (int mp1 = 0; mp1<1300; mp1++) 
	for (int mp2 = 0; mp2<1300; mp2++) {

	  TString HistoName = "nISRJets_mStop-"; HistoName += mp1;
	  HistoName += "_mLSP-"; HistoName += mp2;
	  TH1F *ThisnISRJets = (TH1F*) HistoFile->Get(HistoName);

	  if (ThisnISRJets) {

	    stopMP[nMPfound] = mp1;
	    chi0MP[nMPfound] = mp2;
	    nISRJets[nMPfound]->Add(ThisnISRJets);
	    nMPfound++;
	    
	  }
	  
	}
      
    }

    ofstream OutputMap; OutputMap.open("./StopNeutralinoMap.txt");
    ofstream OutputISRMap; OutputISRMap.open("./StopNeutralinoISRMap.txt");
    
    typedef pair<int, int> MassPoint;
    typedef pair<float, float> StopCrossSection;
    typedef pair<StopCrossSection, int> MassPointParameters;
    
    typedef map<MassPoint, MassPointParameters> MassPointMap;
    MassPointMap StopNeutralinoMap;
    
    typedef map<MassPoint, float> MassPointISRMap;
    MassPointISRMap StopNeutralinoISRMap;
    
    int nMassPoints = 0, TotalEvents = 0;
    
    for (int iStopMass = 0; iStopMass<1300; iStopMass++) 
      for (int iNeutralinoMass = 0; iNeutralinoMass<1300; iNeutralinoMass++) {
       
	float StopMass = 1.*iStopMass;
	float NeutralinoMass  = 1.*iNeutralinoMass;
	
	int nEntries = MassPlan->GetBinContent(MassPlan->FindBin(StopMass, NeutralinoMass));
	//if (iStopMass<250 || (iStopMass==250 && iNeutralinoMass<=163)) nEntries *= 32283695./20333800.;
	TotalEvents += nEntries;
	
	if (nEntries==0) continue;
	
	CrossSection = -1.;
	
	if (iStopMass%5==0.) 
	  GetStopCrossSection(StopMass);
	else {
	  
	  float StopMass1 = 5.*int(StopMass/5);
	  GetStopCrossSection(StopMass1);
	  float CrossSection1 = CrossSection, CrossSection1Uncertainty = CrossSectionUncertainty;
	  
	  float StopMass2 = 5.*(int(StopMass/5)+1);
	  GetStopCrossSection(StopMass2);
	  float CrossSection2 = CrossSection, CrossSection2Uncertainty = CrossSectionUncertainty;
	  
	  CrossSection = CrossSection1 + (StopMass-StopMass1)*(CrossSection2 - CrossSection1)/5.;
	  CrossSectionUncertainty = CrossSection1Uncertainty + (StopMass-StopMass1)*(CrossSection2Uncertainty - CrossSection1Uncertainty)/5.;
	  
	}
	
	if (CrossSection==-1) 
	  cout << "Warning: cross section error for stop mass " << StopMass << endl;
	
	if (Apply2LFilterEff) {

	  CrossSection *= TwoLFilterEff;
	  CrossSectionUncertainty *= TwoLFilterEff;

	}

	MassPoint mp (iStopMass, iNeutralinoMass);
	StopCrossSection cs (CrossSection, CrossSectionUncertainty);
	MassPointParameters mpp (cs, nEntries);
	
	StopNeutralinoMap.insert(std::make_pair(mp, mpp));
	
	OutputMap << "    MassPoint mp" << nMassPoints << " (" << iStopMass << ", " << iNeutralinoMass << ");" << endl;
	OutputMap << "    StopCrossSection cs" << nMassPoints << " (" << CrossSection << ", " << CrossSectionUncertainty << ");" << endl;
	OutputMap << "    MassPointParameters mpp" << nMassPoints << " (cs" << nMassPoints << ", " << nEntries << ");" << endl;
	OutputMap << "    StopNeutralinoMap.insert(std::make_pair(mp" << nMassPoints << ", mpp" << nMassPoints << "));" << endl << endl;

	bool GotISRPoint = false;
	
	for (int imp = 0; imp<nMPfound; imp++) {
	  if (iStopMass==stopMP[imp] && iNeutralinoMass== chi0MP[imp]) {
	    
	    float Dvalue = nISRJets[imp]->Integral(0, 8);
	    for (int ib = 1; ib<=7; ib++) {
	      float ThisBinContent = nISRJets[imp]->GetBinContent(ib);
	      nISRJets[imp]->SetBinContent(ib, ThisBinContent*ISRBinWeight[ib-1]);
	    }

	    Dvalue /= nISRJets[imp]->Integral(0, 8);
	    StopNeutralinoISRMap.insert(std::make_pair(mp, Dvalue));

	    OutputISRMap << "    MassPoint mpisr" << nMassPoints << " (" << iStopMass << ", " << iNeutralinoMass << ");" << endl;
	    OutputISRMap << "    StopNeutralinoISRMap.insert(std::make_pair(mpisr" << nMassPoints << ", " <<  Dvalue << "));" << endl << endl;

	    GotISRPoint = true;

	  }
	}

	if (!GotISRPoint) cout << "  ISR iformation not found for mass point " << nMassPoints << endl;
	
	nMassPoints++;
	
      }

    if (IsTest) {

      MassPoint ThisMassPoint (FirstStopMass, FirstNeutralinoMass);
      MassPointParameters TheseMassPointParameters = StopNeutralinoMap.at(ThisMassPoint);
      StopCrossSection ThisStopCrossSection = TheseMassPointParameters.first;
      int SampleSize = TheseMassPointParameters.second;
      
      cout << "Stop mass: " << FirstStopMass << ", Neutralino mass: " << FirstNeutralinoMass 
	   << ", Cross section " << ThisStopCrossSection.first 
	   << ", Cross section uncertainty " << ThisStopCrossSection.second
	   << ", Sample size " << SampleSize;

      float ThisDvalue = StopNeutralinoISRMap.at(ThisMassPoint); 

      cout << ", D = " << ThisDvalue;

      cout << endl;

    }
    
    cout << "Total Events: " << TotalEvents << endl;

  }
  
}
