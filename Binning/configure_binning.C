#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>

#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TGaxis.h>

#include "Binning.h"
#include "Binning.cxx"
#endif

void configure_binning(){

  const int N_cost_bins_2pt6 = 22;
  int min_cost_bin_2pt6[N_cost_bins_2pt6] = {1,11,21,26,31,36,41,43,45,47,49,51,53,55,57,59,61,66,71,76,81,91};
  int max_cost_bin_2pt6[N_cost_bins_2pt6] = {10,20,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100};

  const int N_phi_bins_2pt6 = 10;
  int min_phi_bin_2pt6[N_phi_bins_2pt6] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_2pt6[N_phi_bins_2pt6] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_2pt6 = new Binning();
  binning_2pt6 -> ConfigureBinValues(N_cost_bins_2pt6,min_cost_bin_2pt6,max_cost_bin_2pt6,N_phi_bins_2pt6,min_phi_bin_2pt6,max_phi_bin_2pt6);

  TFile *file_2pt6 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_2pt6.root","RECREATE");
  binning_2pt6 -> Write();
  file_2pt6 -> Close();

  const int N_cost_bins_6pt12 = 17;
  int min_cost_bin_6pt12[N_cost_bins_6pt12] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_6pt12[N_cost_bins_6pt12] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_6pt12 = 10;
  int min_phi_bin_6pt12[N_phi_bins_6pt12] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_6pt12[N_phi_bins_6pt12] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_6pt12 = new Binning();
  binning_6pt12 -> ConfigureBinValues(N_cost_bins_6pt12,min_cost_bin_6pt12,max_cost_bin_6pt12,N_phi_bins_6pt12,min_phi_bin_6pt12,max_phi_bin_6pt12);

  TFile *file_6pt12 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_6pt12.root","RECREATE");
  binning_6pt12 -> Write();
  file_6pt12 -> Close();

  const int N_cost_bins_6pt10 = 17;
  int min_cost_bin_6pt10[N_cost_bins_6pt10] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_6pt10[N_cost_bins_6pt10] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_6pt10 = 10;
  int min_phi_bin_6pt10[N_phi_bins_6pt10] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_6pt10[N_phi_bins_6pt10] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_6pt10 = new Binning();
  binning_6pt10 -> ConfigureBinValues(N_cost_bins_6pt10,min_cost_bin_6pt10,max_cost_bin_6pt10,N_phi_bins_6pt10,min_phi_bin_6pt10,max_phi_bin_6pt10);

  TFile *file_6pt10 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_6pt10.root","RECREATE");
  binning_6pt10 -> Write();
  file_6pt10 -> Close();

  //============================================================================
  const int N_cost_bins_6pt10_test = 19;
  int min_cost_bin_6pt10_test[N_cost_bins_6pt10_test] = {1,11,16,21,26,31,36,41,45,49,53,57,61,66,71,76,81,86,91};
  int max_cost_bin_6pt10_test[N_cost_bins_6pt10_test] = {10,15,20,25,30,35,40,44,48,52,56,60,65,70,75,80,,85,90,100};

  const int N_phi_bins_6pt10_test = 10;
  int min_phi_bin_6pt10_test[N_phi_bins_6pt10_test] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_6pt10_test[N_phi_bins_6pt10_test] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_6pt10_test = new Binning();
  binning_6pt10_test -> ConfigureBinValues(N_cost_bins_6pt10_test,min_cost_bin_6pt10_test,max_cost_bin_6pt10_test,N_phi_bins_6pt10_test,min_phi_bin_6pt10_test,max_phi_bin_6pt10_test);

  TFile *file_6pt10_test = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_6pt10_test.root","RECREATE");
  binning_6pt10_test -> Write();
  file_6pt10_test -> Close();

  //============================================================================
  // Attempt to have 3 bins in pT

  const int N_cost_bins_2pt4 = 17;
  int min_cost_bin_2pt4[N_cost_bins_2pt4] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_2pt4[N_cost_bins_2pt4] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_2pt4 = 10;
  int min_phi_bin_2pt4[N_phi_bins_2pt4] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_2pt4[N_phi_bins_2pt4] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_2pt4 = new Binning();
  binning_2pt4 -> ConfigureBinValues(N_cost_bins_2pt4,min_cost_bin_2pt4,max_cost_bin_2pt4,N_phi_bins_2pt4,min_phi_bin_2pt4,max_phi_bin_2pt4);

  TFile *file_2pt4 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_2pt4.root","RECREATE");
  binning_2pt4 -> Write();
  file_2pt4 -> Close();

  //============================================================================
  const int N_cost_bins_2pt4_test = 19;
  int min_cost_bin_2pt4_test[N_cost_bins_2pt4_test] = {1,11,16,21,26,31,36,41,45,49,53,57,61,66,71,76,81,86,91};
  int max_cost_bin_2pt4_test[N_cost_bins_2pt4_test] = {10,15,20,25,30,35,40,44,48,52,56,60,65,70,75,80,,85,90,100};

  const int N_phi_bins_2pt4_test = 10;
  int min_phi_bin_2pt4_test[N_phi_bins_2pt4_test] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_2pt4_test[N_phi_bins_2pt4_test] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_2pt4_test = new Binning();
  binning_2pt4_test -> ConfigureBinValues(N_cost_bins_2pt4_test,min_cost_bin_2pt4_test,max_cost_bin_2pt4_test,N_phi_bins_2pt4_test,min_phi_bin_2pt4_test,max_phi_bin_2pt4_test);

  TFile *file_2pt4_test = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_2pt4_test.root","RECREATE");
  binning_2pt4_test -> Write();
  file_2pt4_test -> Close();

  //============================================================================

  const int N_cost_bins_4pt7 = 17;
  int min_cost_bin_4pt7[N_cost_bins_4pt7] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_4pt7[N_cost_bins_4pt7] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_4pt7 = 10;
  int min_phi_bin_4pt7[N_phi_bins_4pt7] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_4pt7[N_phi_bins_4pt7] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_4pt7 = new Binning();
  binning_4pt7 -> ConfigureBinValues(N_cost_bins_4pt7,min_cost_bin_4pt7,max_cost_bin_4pt7,N_phi_bins_4pt7,min_phi_bin_4pt7,max_phi_bin_4pt7);

  TFile *file_4pt7 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_4pt7.root","RECREATE");
  binning_4pt7 -> Write();
  file_4pt7 -> Close();

  const int N_cost_bins_7pt10 = 17;
  int min_cost_bin_7pt10[N_cost_bins_7pt10] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_7pt10[N_cost_bins_7pt10] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_7pt10 = 10;
  int min_phi_bin_7pt10[N_phi_bins_7pt10] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_7pt10[N_phi_bins_7pt10] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_7pt10 = new Binning();
  binning_7pt10 -> ConfigureBinValues(N_cost_bins_7pt10,min_cost_bin_7pt10,max_cost_bin_7pt10,N_phi_bins_7pt10,min_phi_bin_7pt10,max_phi_bin_7pt10);

  TFile *file_7pt10 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_7pt10.root","RECREATE");
  binning_7pt10 -> Write();
  file_7pt10 -> Close();

  const int N_cost_bins_4pt6 = 17;
  int min_cost_bin_4pt6[N_cost_bins_4pt6] = {1,11,21,26,31,36,41,45,49,53,57,61,66,71,76,81,91};
  int max_cost_bin_4pt6[N_cost_bins_4pt6] = {10,20,25,30,35,40,44,48,52,56,60,65,70,75,80,90,100};

  const int N_phi_bins_4pt6 = 10;
  int min_phi_bin_4pt6[N_phi_bins_4pt6] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_4pt6[N_phi_bins_4pt6] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_4pt6 = new Binning();
  binning_4pt6 -> ConfigureBinValues(N_cost_bins_4pt6,min_cost_bin_4pt6,max_cost_bin_4pt6,N_phi_bins_4pt6,min_phi_bin_4pt6,max_phi_bin_4pt6);

  TFile *file_4pt6 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_4pt6.root","RECREATE");
  binning_4pt6 -> Write();
  file_4pt6 -> Close();

  //============================================================================
  const int N_cost_bins_4pt6_test = 19;
  int min_cost_bin_4pt6_test[N_cost_bins_4pt6_test] = {1,11,16,21,26,31,36,41,45,49,53,57,61,66,71,76,81,86,91};
  int max_cost_bin_4pt6_test[N_cost_bins_4pt6_test] = {10,15,20,25,30,35,40,44,48,52,56,60,65,70,75,80,,85,90,100};

  const int N_phi_bins_4pt6_test = 10;
  int min_phi_bin_4pt6_test[N_phi_bins_4pt6_test] = {1,9,17,21,24,26,28,31,35,43};
  int max_phi_bin_4pt6_test[N_phi_bins_4pt6_test] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning_4pt6_test = new Binning();
  binning_4pt6_test -> ConfigureBinValues(N_cost_bins_4pt6_test,min_cost_bin_4pt6_test,max_cost_bin_4pt6_test,N_phi_bins_4pt6_test,min_phi_bin_4pt6_test,max_phi_bin_4pt6_test);

  TFile *file_4pt6_test = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_4pt6_test.root","RECREATE");
  binning_4pt6_test -> Write();
  file_4pt6_test -> Close();

  //============================================================================


}
