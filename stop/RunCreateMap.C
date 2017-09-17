R__LOAD_LIBRARY(CreateMap_C.so) 

void RunCreateMap() {

  cout << "RunCreateMap" << endl; 

  TString DataRangeString = gSystem->Getenv("DATARANGE");
  int DataRange = DataRangeString.Atoi() - 1; 

  CreateMap(false, DataRange);

}
