void rootlogon() 
{
  // change below to your preferred default value of omega
  Double_t defOmega = -1; // -1 as well as 0 means hard-sphere approximation (default)
  //
  const char *version = "v3.3.1";
  gSystem->Load("libMathMore");  
  TString libnucgenPath = Form("%s/trnucgen/libtrnucgen.so", gSystem->WorkingDirectory());
  if (gSystem->AccessPathName(libnucgenPath)==0) {
    gSystem->AddDynamicPath(Form("%s/trnucgen/", gSystem->WorkingDirectory()));
    //R__LOAD_LIBRARY(trnucgen/libtrnucgen.so)
    int status = gSystem->Load("libtrnucgen.so");
    if (status != 0) {
      cout << "Error loading libtrnucgen.so" << endl;
    } else {
      cout << "Setting up TGlauberMC " <<version<< " with libtrnucgen" << endl;
      gSystem->AddIncludePath("-DUSE_TRNUCGEN -I ./trnucgen");
    }
  } else {
    cout << "Setting up TGlauberMC " <<version<< " without libstrnucgen" << endl;
  }
  if (1) { // set random seed to 0
    delete gRandom;
    gRandom = new TRandom3(0);
  }
  // load the runglauber macro
  gROOT->LoadMacro("runglauber_v3.3.C+");
  if (defOmega>=0) {
    cout << "Setting default value of omega to " << defOmega << endl;
    gROOT->ProcessLine(Form("TGlauberMC::SetDefOmega(%f)",defOmega));
  }
}
