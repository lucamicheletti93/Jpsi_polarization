#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>

#include <TROOT.h>
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
#include <TCollection.h>
#include <TKey.h>
#include <TGaxis.h>

#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

Double_t Func_VWG(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_fix(Double_t *, Double_t *);
Double_t Func_Psi2s_CB2(Double_t *, Double_t *);
Double_t Func_Psi2s_CB2_fix(Double_t *, Double_t *);
Double_t Func_tot(Double_t *, Double_t *);

void single_histo_fit();
void loop_on_histos();
void fit_of_minv(TH1D *, int, int, string, string);

double PI = TMath::Pi();
Double_t scaling_factor = 1.05154; //factor introduced to pass from the sigma of Jpsi to the sigma of Psi(2S)
double nJpsi = 0, statJpsi = 0;
double sigmaJpsi = 0, errSigmaJpsi = 0;
double ChiSquare_NDF = 0;
char fit_status[10];
vector <double> CostValues;
vector <double> PhiValues;

//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void single_histo_fit(){

  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  string dataset = "6pt10";
  string fileInputName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/mass_histos_cost_phi_" + dataset + ".root";
  TFile *fileInput = new TFile(fileInputName.c_str());

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";
  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  //TH1D *histMinv = (TH1D*) fileInput -> Get("HE_11cost20_17phi20");
  TH1D *histMinv = (TH1D*) fileInput -> Get("HE_11cost20_35phi42");
  histMinv -> SetMarkerStyle(20);
  histMinv -> SetMarkerSize(0.7);
  string outputDir = "No outputDir";
  fit_of_minv(histMinv,1,1,dataset,outputDir);
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void loop_on_histos(){

  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  string dataset = "4pt6";
  //string plotsDirectory = "mkdir binned_1D_" + dataset;
  string plotsDirectory = "mkdir binned_1D_" + dataset + "_test";
  gSystem -> Exec(plotsDirectory.c_str());
  //string fileInputName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/mass_histos_cost_phi_" + dataset + ".root";
  //string fileInputName = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/" + dataset + ".root";
  string fileInputName = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/" + dataset + "_test.root";
  TFile *fileInput = new TFile(fileInputName.c_str());

  if(dataset == "0pt12"){
    string outputDir = "binned_1D_" + dataset;
    TH1D *histMinv = (TH1D*) fileInput -> Get("histIntegrated");
    fit_of_minv(histMinv,100,100,dataset,outputDir);
    //histMinv -> Draw();
    return;
  }

  //============================================================================
  printf("---> Setting histograms bins ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);

  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  //string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";
  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + "_test.root";
  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  binning -> PrintBinValues();
  fileInput -> ls();

  char histMinvName[100];
  vector <string> fit_errors;

  TH2D *histNJpsi = new TH2D("histNJpsi","histNJpsi",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *histSigmaJpsi = new TH2D("histSigmaJpsi","histSigmaJpsi",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  gStyle -> SetPaintTextFormat("2.0f");

  //============================================================================
  printf("Fitting along Cost slices ... \n");
  //============================================================================
  TH1D *histNJpsiCost = new TH1D("histNJpsiCost","",NCostBins,&CostValues[0]);
  TH1D *histSigmaJpsiCost = new TH1D("histSigmaJpsiCost","",NCostBins,&CostValues[0]);

  /*for(int i = 1;i < NCostBins-1;i++){
    sprintf(histMinvName,"HE_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[1],PhiBinsMax[1]);
    TH1D *histMinvSliceCost = (TH1D*) fileInput -> Get(histMinvName);
    for(int j = 2;j < NPhiBins-1;j++){
      sprintf(histMinvName,"HE_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      TH1D *tmpMinvSlice = (TH1D*) fileInput -> Get(histMinvName);
      histMinvSliceCost -> Add(tmpMinvSlice);
    }
    histMinvSliceCost -> SetName(Form("HE_%icost%i",CostBinsMin[i],CostBinsMax[i]));
    //string outputDir = "binned_1D_" + dataset;
    string outputDir = "binned_1D_" + dataset + "_test";
    fit_of_minv(histMinvSliceCost,i,100,dataset,outputDir);
    histMinvSliceCost -> Draw();
    histNJpsiCost -> SetBinContent(i+1,nJpsi);
    histNJpsiCost -> SetBinError(i+1,statJpsi);
    histSigmaJpsiCost -> SetBinContent(i+1,sigmaJpsi);
    histSigmaJpsiCost -> SetBinError(i+1,errSigmaJpsi);
  }*/

  for(int i = 1;i < NCostBins-1;i++){
    sprintf(histMinvName,"histCost_%i",i);
    TH1D *histMinvSliceCost = (TH1D*) fileInput -> Get(histMinvName);
    histMinvSliceCost -> SetName(Form("HE_%icost%i",CostBinsMin[i],CostBinsMax[i]));
    //string outputDir = "binned_1D_" + dataset;
    string outputDir = "binned_1D_" + dataset + "_test";
    fit_of_minv(histMinvSliceCost,i,100,dataset,outputDir);
    histMinvSliceCost -> Draw();
    histNJpsiCost -> SetBinContent(i+1,nJpsi);
    histNJpsiCost -> SetBinError(i+1,statJpsi);
    histSigmaJpsiCost -> SetBinContent(i+1,sigmaJpsi);
    histSigmaJpsiCost -> SetBinError(i+1,errSigmaJpsi);
  }

  //============================================================================
  printf("Fitting along Phi slices ... \n");
  //============================================================================
  TH1D *histNJpsiPhi = new TH1D("histNJpsiPhi","",NPhiBins,&PhiValues[0]);
  TH1D *histSigmaJpsiPhi = new TH1D("histSigmaJpsiPhi","",NPhiBins,&PhiValues[0]);

  /*for(int i = 1;i < NPhiBins-1;i++){
    sprintf(histMinvName,"HE_%icost%i_%iphi%i",CostBinsMin[1],CostBinsMax[1],PhiBinsMin[i],PhiBinsMax[i]);
    TH1D *histMinvSlicePhi = (TH1D*) fileInput -> Get(histMinvName);
    for(int j = 2;j < NCostBins-1;j++){
      sprintf(histMinvName,"HE_%icost%i_%iphi%i",CostBinsMin[j],CostBinsMax[j],PhiBinsMin[i],PhiBinsMax[i]);
      TH1D *tmpMinvSlice = (TH1D*) fileInput -> Get(histMinvName);
      histMinvSlicePhi -> Add(tmpMinvSlice);
      cout << histMinvName << endl;
    }
    histMinvSlicePhi -> SetName(Form("HE_%iphi%i",PhiBinsMin[i],PhiBinsMax[i]));
    string outputDir = "binned_1D_" + dataset;
    fit_of_minv(histMinvSlicePhi,100,i,dataset,outputDir);
    histNJpsiPhi -> SetBinContent(i+1,nJpsi);
    histNJpsiPhi -> SetBinError(i+1,statJpsi);
    histSigmaJpsiPhi -> SetBinContent(i+1,sigmaJpsi);
    histSigmaJpsiPhi -> SetBinError(i+1,errSigmaJpsi);
  }*/

  for(int i = 1;i < NPhiBins-1;i++){
    sprintf(histMinvName,"histPhi_%i",i);
    TH1D *histMinvSlicePhi = (TH1D*) fileInput -> Get(histMinvName);
    histMinvSlicePhi -> SetName(Form("HE_%iphi%i",PhiBinsMin[i],PhiBinsMax[i]));
    //string outputDir = "binned_1D_" + dataset;
    string outputDir = "binned_1D_" + dataset + "_test";
    fit_of_minv(histMinvSlicePhi,100,i,dataset,outputDir);
    histNJpsiPhi -> SetBinContent(i+1,nJpsi);
    histNJpsiPhi -> SetBinError(i+1,statJpsi);
    histSigmaJpsiPhi -> SetBinContent(i+1,sigmaJpsi);
    histSigmaJpsiPhi -> SetBinError(i+1,errSigmaJpsi);
  }

  //============================================================================
  printf("---> Saving results ... \n");
  //============================================================================
  //string fileSigmaJpsiName = "binned_1D_" + dataset + "/" + dataset + ".root";
  string fileSigmaJpsiName = "binned_1D_" + dataset + "_test/" + dataset + ".root";
  TFile *fileSigmaJpsi = new TFile(fileSigmaJpsiName.c_str(),"RECREATE");
  histNJpsiCost -> Write();
  histSigmaJpsiCost -> Write();
  histNJpsiPhi -> Write();
  histSigmaJpsiPhi -> Write();
  fileSigmaJpsi -> Close();
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void fit_of_minv(TH1D *histMinv, int counter_cost, int counter_phi, string dataset, string outputDir){

  string sig = "CB2";
  string bck = "VWG";
  string fit_range = "min";
  double min_fit_range = 2.1;
  double max_fit_range = 4.9;
  string tails = "pp8";
  string tails_fix = "yes";

  gStyle -> SetOptStat(0);
  TGaxis::SetMaxDigits(2);

  if(dataset.compare("2pt6") == 0){
    if(counter_cost == 1 || counter_cost == (int) CostValues.size() - 2){histMinv -> Rebin(2);}
  }
  if(dataset.compare("2pt4") == 0){
    //if(counter_cost == 1 || counter_cost == (int) CostValues.size() - 2){histMinv -> Rebin(2);}
  }
  if(dataset.compare("6pt12") == 0){histMinv -> Rebin(2);}
  //if(dataset.compare("6pt10") == 0){histMinv -> Rebin(2);}
  if(dataset.compare("4pt7") == 0){histMinv -> Rebin(2);}
  if(dataset.compare("7pt10") == 0){histMinv -> Rebin(2);}
  //if(dataset.compare("4pt6") == 0){histMinv -> Rebin(2);}

  //============================================================================
  // FIT STARTING PARAMETERS
  //============================================================================
  double par_bck[4], par_signal[8];
  if(dataset == "0pt12"){
    par_bck[0] = 500000.; par_bck[1] = 0.6; par_bck[2] = 0.2; par_bck[3] = 0.2;
    par_signal[0] = 500.; par_signal[1] = 3.096; par_signal[2] = 7.0e-02;
    par_signal[3] = 1.00011e+00; par_signal[4] = 3.70078e+00; par_signal[5] = 1.68359e+00; par_signal[6] = 3.63002e+00; par_signal[7] = 0.01;
    string histTmpName = histMinv -> GetName();
    //if(histTmpName == "HE_11cost20_9phi16"){histMinv -> Rebin(2); min_fit_range = 2.; max_fit_range = 5.; par_bck[2] = 1.;}
  }

  if(dataset == "2pt4"){
    par_bck[0] = 500000.; par_bck[1] = 0.6; par_bck[2] = 0.2; par_bck[3] = 0.2;
    par_signal[0] = 500.; par_signal[1] = 3.096; par_signal[2] = 7.0e-02;
    par_signal[3] = 1.00011e+00; par_signal[4] = 3.70078e+00; par_signal[5] = 1.68359e+00; par_signal[6] = 3.63002e+00; par_signal[7] = 0.01;
    string histTmpName = histMinv -> GetName();
    //if(histTmpName == "HE_11cost20_9phi16"){histMinv -> Rebin(2); min_fit_range = 2.; max_fit_range = 5.; par_bck[2] = 1.;}
  }

  if(dataset == "4pt6"){
    par_bck[0] = 500000.; par_bck[1] = 0.6; par_bck[2] = 0.2; par_bck[3] = 0.2;
    par_signal[0] = 500.; par_signal[1] = 3.096; par_signal[2] = 7.0e-02;
    par_signal[3] = 1.00011e+00; par_signal[4] = 3.70078e+00; par_signal[5] = 1.68359e+00; par_signal[6] = 3.63002e+00; par_signal[7] = 0.01;
    string histTmpName = histMinv -> GetName();
    //if(histTmpName == "HE_11cost20_9phi16"){histMinv -> Rebin(2); min_fit_range = 2.; max_fit_range = 5.; par_bck[2] = 1.;}
  }

  if(dataset == "6pt10"){
    par_bck[0] = 500000.; par_bck[1] = 0.6; par_bck[2] = 0.2; par_bck[3] = 0.2;
    par_signal[0] = 500.; par_signal[1] = 3.096; par_signal[2] = 7.0e-02;
    par_signal[3] = 1.00011e+00; par_signal[4] = 3.70078e+00; par_signal[5] = 1.68359e+00; par_signal[6] = 3.63002e+00; par_signal[7] = 0.01;
    string histTmpName = histMinv -> GetName();
    // rebin
    //if(histTmpName == "HE_11cost20_9phi16"){histMinv -> Rebin(2);}
    //if(histTmpName == "HE_11cost20_17phi20"){histMinv -> Rebin(4);}
    //if(histTmpName == "HE_11cost20_21phi23"){histMinv -> Rebin(2);}
    //if(histTmpName == "HE_11cost20_24phi25"){histMinv -> Rebin(2);}
    //if(histTmpName == "HE_11cost20_26phi27"){histMinv -> Rebin(2);}
    //if(histTmpName == "HE_11cost20_28phi30"){histMinv -> Rebin(2);}
  }

  //double par_bck[4] = {500000.,0.6,0.2,0.2};
  // Signal
  //double par_signal[8] = {50000.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01};
  //double par_signal[8] = {500.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01}; // fot 4 < pT < 7 GeV/c
  //double par_signal[8] = {50.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01}; // fot 7 < pT < 10 GeV/c
  //============================================================================
  // FIT OF THE BACKGROUND
  //============================================================================
  TH1D *histo_for_bck = (TH1D*) histMinv -> Clone("histo_for_bck");
  double m_width = histo_for_bck -> GetBinWidth(1);
  double m_min = histo_for_bck -> FindBin(2.9);
  double m_max = histo_for_bck -> FindBin(3.3);

  cout << "-----------------------------------" << endl;
  cout << "Bin width : " << m_width << endl;
  cout << "bin min : " << m_min << "; bin max : " << m_max << endl;
  cout << "-----------------------------------" << endl;

  for(int i = m_min; i < m_max;i++){
    histo_for_bck -> SetBinContent(i+1,0);
    histo_for_bck -> SetBinError(i+1,0);
  }

  double *par_bck_ptr = par_bck;

  TF1 *func_bck_VWG = new TF1("func_bck_VWG",Func_VWG,min_fit_range,max_fit_range,4);
  //TF1 *func_bck_VWG = new TF1("func_bck_VWG",Func_VWG,2.1,4.9,4);
  func_bck_VWG -> SetLineColor(kRed);
  func_bck_VWG -> SetParameters(par_bck);
  histo_for_bck -> Fit(func_bck_VWG,"R0");

  //============================================================================
  // TAILS {alpha1,n1,alpha2,n2}
  //============================================================================
  //double tails_par[4] = {1.012,2.874,2.232,2.924}; //p-A MC COLLISIONS
  //double tails_par[4] = {1.089,3.393,2.097,8.694}; //p-p 8 TeV
  double tails_par[4] = {0.98,6.97,1.86,14.99}; //p-p 13 TeV
  //double tails_par[4] = {1.06,3.23,2.55,1.56}; //GEANT 4
  //double tails_par[4] = {0.97,3.98,2.3,3.03}; //GEANT 3
  //============================================================================
  // FIT OF THE JPSI SIGNAL
  //============================================================================

  TF1 *func_Jpsi_CB2 = new TF1("func_Jpsi_CB2",Func_Jpsi_CB2,2.9,3.2,11);
  func_Jpsi_CB2 -> SetNpx(100000);
  func_Jpsi_CB2 -> SetParameter(4,par_signal[0]);
  func_Jpsi_CB2 -> SetParameter(5,par_signal[1]);
  func_Jpsi_CB2 -> FixParameter(6,par_signal[2]);
  func_Jpsi_CB2 -> FixParameter(7,tails_par[0]);
  func_Jpsi_CB2 -> FixParameter(8,tails_par[1]);
  func_Jpsi_CB2 -> FixParameter(9,tails_par[2]);
  func_Jpsi_CB2 -> FixParameter(10,tails_par[3]);
  histMinv -> Fit(func_Jpsi_CB2,"R0");

  //============================================================================
  // FIT OF THE TOTAL SPECTRUM
  //============================================================================
  //TF1 *func_tot = new TF1("func_tot",Func_tot,min_fit_range,max_fit_range,12);
  TF1 *func_tot = new TF1("func_tot",Func_tot,min_fit_range,max_fit_range,11);
  TFitResultPtr fit_ptr;
  for(int i = 0;i < 100;i++){
    if(i == 0){
      func_tot -> SetParameter(0,func_bck_VWG -> GetParameter(0));
      func_tot -> SetParameter(1,func_bck_VWG -> GetParameter(1));
      func_tot -> SetParameter(2,func_bck_VWG -> GetParameter(2));
      func_tot -> SetParameter(3,func_bck_VWG -> GetParameter(3));
      func_tot -> SetParameter(4,func_Jpsi_CB2 -> GetParameter(4));
      func_tot -> SetParLimits(4,0,10000000);
      func_tot -> SetParameter(5,3.096);
      func_tot -> SetParLimits(5,2.9,3.2);
      func_tot -> SetParameter(6,7.0e-02);
      func_tot -> SetParLimits(6,2.0e-02,2.0e-01);
      //if(counter_cost == 1){func_tot -> FixParameter(6,1.0e-01);}
      //else{
        //func_tot -> SetParameter(6,7.0e-02);
        //func_tot -> SetParLimits(6,3.0e-02,2.0e-01);
      //}


      //func_tot -> SetParLimits(6,5.0e-02,1.2e-01);
      //func_tot -> SetParLimits(6,5.0e-02,9.0e-02); // To be resetted
      //========================================================================
      // FIX THE WIDTH TO MC
      //char INPUT_FILE_NAME[300] = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/accxeff_new.root";
      //TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
      //TH2D *histo_sigmaJpsi = (TH2D*) input_file -> Get("histo_sigmaJpsi");
      //double sigmaJpsi_MC = histo_sigmaJpsi -> GetBinContent(counter_cost+1,counter_phi+1);
      //func_tot -> FixParameter(6,sigmaJpsi_MC);
      //========================================================================

      if(strcmp(tails_fix.c_str(),"yes")==0){
        func_tot -> FixParameter(7,func_Jpsi_CB2 -> GetParameter(7));
        func_tot -> FixParameter(8,func_Jpsi_CB2 -> GetParameter(8));
        func_tot -> FixParameter(9,func_Jpsi_CB2 -> GetParameter(9));
        func_tot -> FixParameter(10,func_Jpsi_CB2 -> GetParameter(10));
      }
      if(strcmp(tails_fix.c_str(),"no")==0){
        func_tot -> SetParameter(7,func_Jpsi_CB2 -> GetParameter(7));
        func_tot -> SetParameter(8,func_Jpsi_CB2 -> GetParameter(8));
        func_tot -> SetParameter(9,func_Jpsi_CB2 -> GetParameter(9));
        func_tot -> SetParameter(10,func_Jpsi_CB2 -> GetParameter(10));
      }
      //func_tot -> SetParameter(11,func_Jpsi_CB2 -> GetParameter(11));
    }
    else{func_tot -> SetParameters(func_tot -> GetParameters());}
    fit_ptr = (TFitResultPtr) histMinv -> Fit(func_tot,"RLS0");
    if(gMinuit->fCstatu.Contains("CONVERGED")) break;
  }

  //============================================================================
  // EXTRA CONDITION FOR BAD FITS
  //============================================================================
  /*double min_sigma_limit[1], max_sigma_limit[1];
  func_tot -> GetParLimits(6,min_sigma_limit[0],max_sigma_limit[0]);
  if(func_tot -> GetParameter(6) > max_sigma_limit[0] - 5.0e-03){
    char INPUT_FILE_NAME[300] = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/accxeff_new.root";
    TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
    TH2D *histo_sigmaJpsi = (TH2D*) input_file -> Get("histo_sigmaJpsi");
    double sigmaJpsi_MC = histo_sigmaJpsi -> GetBinContent(counter_cost,counter_phi);
    func_tot -> FixParameter(6,6.5e-02);
    fit_ptr = (TFitResultPtr) histMinv -> Fit(func_tot,"RLS0");
  }*/
  //============================================================================

  ChiSquare_NDF = func_tot -> GetChisquare()/func_tot -> GetNDF();

  //printf("\n\nfit status: %s \n\n",gMinuit.fCstatu.Data());
  if(gMinuit -> fCstatu.Contains("FAILED")){
    cout << "WARNING : FIT STATUS FAILED" << endl;
    sprintf(fit_status,"FAILED");
    return;
  }
  else{sprintf(fit_status,"SUCCESS");}

  TMatrixDSym cov = fit_ptr -> GetCovarianceMatrix();

  double *fullmat;
  fullmat = cov.GetMatrixArray();

  //PSI(2S) MATRIX
  double Jpsi_mat[49];
  for(Int_t i = 0;i < 7;i++){
    for(Int_t j = 0;j < 7;j++){
        Jpsi_mat[7*i+j] = fullmat[48+j+11*i];
    }
  }

  double Jpsi_par[7];
  for(int i = 0;i < 7;i++){
    Jpsi_par[i] = func_tot -> GetParameter(4+i);
  }

  //============================================================================
  // PLOT OF JPSI AND PSI(2S) SHAPES
  //============================================================================
  TF1 *func_bck_VWG_fix = new TF1("func_bck_VWG_fix",Func_VWG,2.,5.,4);
  func_bck_VWG_fix -> SetNpx(1000);
  func_bck_VWG_fix -> SetParameter(0,func_tot -> GetParameter(0));
  func_bck_VWG_fix -> SetParameter(1,func_tot -> GetParameter(1));
  func_bck_VWG_fix -> SetParameter(2,func_tot -> GetParameter(2));
  func_bck_VWG_fix -> SetParameter(3,func_tot -> GetParameter(3));
  func_bck_VWG_fix -> SetLineStyle(4);
  func_bck_VWG_fix -> SetLineColor(kBlue+1);
  func_bck_VWG_fix -> Draw("same");

  TF1 *func_Jpsi_CB2_fix = new TF1("func_Jpsi_CB2_fix",Func_Jpsi_CB2_fix,2.,5.,7);
  func_Jpsi_CB2_fix -> SetNpx(1000);
  func_Jpsi_CB2_fix -> SetParameter(0,func_tot -> GetParameter(4));
  func_Jpsi_CB2_fix -> SetParameter(1,func_tot -> GetParameter(5));
  func_Jpsi_CB2_fix -> SetParameter(2,func_tot -> GetParameter(6));
  func_Jpsi_CB2_fix -> SetParameter(3,func_tot -> GetParameter(7));
  func_Jpsi_CB2_fix -> SetParameter(4,func_tot -> GetParameter(8));
  func_Jpsi_CB2_fix -> SetParameter(5,func_tot -> GetParameter(9));
  func_Jpsi_CB2_fix -> SetParameter(6,func_tot -> GetParameter(10));
  func_Jpsi_CB2_fix -> SetLineStyle(2);
  func_Jpsi_CB2_fix -> Draw("same");

  double mass_Jpsi = func_tot -> GetParameter(5);
  sigmaJpsi = func_tot -> GetParameter(6);
  errSigmaJpsi = func_tot -> GetParError(6);
  double sigma_min_Jpsi = func_tot -> GetParameter(5) - 3*(func_tot -> GetParameter(6));
  double sigma_max_Jpsi = func_tot -> GetParameter(5) + 3*(func_tot -> GetParameter(6));
  double N_Jpsi_3sigma = func_Jpsi_CB2_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double N_bck_Jpsi_3sigma = func_bck_VWG_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double SB_Jpsi = N_Jpsi_3sigma/N_bck_Jpsi_3sigma;

  nJpsi = func_Jpsi_CB2_fix -> Integral(0,5)/m_width;
  statJpsi = func_Jpsi_CB2_fix -> IntegralError(0.,5.,Jpsi_par,Jpsi_mat)/m_width;

  //============================================================================
  // PLOT OF THE TOTAL SPECTRUM
  //============================================================================
  char title[100];
  double max_histo_value = histMinv -> GetMaximum();
  max_histo_value = max_histo_value + 0.6*max_histo_value;

  TH2D *h_spectrum = new TH2D("h_spectrum","",120,2,5,100,0,max_histo_value);
  h_spectrum -> GetXaxis() -> SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  h_spectrum -> GetXaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetXaxis() -> SetTitleOffset(0.95);
  h_spectrum -> GetXaxis() -> SetLabelSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitle(Form("Counts per %3.0f MeV/#it{c}^{2}",histMinv -> GetBinWidth(1)*1000));
  h_spectrum -> GetYaxis() -> CenterTitle(true);
  h_spectrum -> GetYaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitleOffset(1.2);
  h_spectrum -> GetYaxis() -> SetLabelSize(0.05);

  sprintf(title,"N_{J/#psi} = %2.0f #pm %2.0f",nJpsi,statJpsi);
  TLatex *lat2 = new TLatex(0.5,0.82,title);
  lat2 -> SetTextSize(0.04);
  lat2 -> SetNDC();
  lat2 -> SetTextFont(42);

  sprintf(title,"#it{m}_{J/#psi} = %4.3f GeV/#it{c}^{2}",mass_Jpsi);
  TLatex *lat3 = new TLatex(0.5,0.76,title);
  lat3 -> SetTextSize(0.04);
  lat3 -> SetNDC();
  lat3 -> SetTextFont(42);

  sprintf(title,"#it{#sigma}_{J/#psi} = %3.0f #pm %3.0f MeV/#it{c}^{2}",sigmaJpsi*1000,errSigmaJpsi*1000);
  TLatex *lat4 = new TLatex(0.5,0.70,title);
  lat4 -> SetTextSize(0.04);
  lat4 -> SetNDC();
  lat4 -> SetTextFont(42);

  if(counter_cost == 100) sprintf(title,"%2.1f < cos#it{#theta}^{HX} < %2.1f",-1.,1.);
  else{sprintf(title,"%3.2f < cos#it{#theta}^{HX} < %3.2f",CostValues[counter_cost],CostValues[counter_cost+1]);}
  TLatex *lat5 = new TLatex(0.55,0.62,title);
  lat5 -> SetTextSize(0.04);
  lat5 -> SetNDC();
  lat5 -> SetTextFont(42);

  if(counter_phi == 100) sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",0.,PI);
  else{sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",PhiValues[counter_phi],PhiValues[counter_phi+1]);}
  TLatex *lat6 = new TLatex(0.55,0.55,title);
  lat6 -> SetTextSize(0.04);
  lat6 -> SetNDC();
  lat6 -> SetTextFont(42);

  sprintf(title,"#chi^{2}/ndf = %3.1f",ChiSquare_NDF);
  TLatex *lat7 = new TLatex(0.55,0.48,title);
  lat7 -> SetTextSize(0.04);
  lat7 -> SetNDC();
  lat7 -> SetTextFont(42);

  //============================================================================
  // DRAVING THE TH2D
  //============================================================================
  TCanvas *c_spectrum = new TCanvas("c_spectrum","c_spectrum",65,73,900,806);
  c_spectrum -> Range(1.825,-6776.052,5.019444,37862.12);
  c_spectrum -> SetFillColor(0);
  c_spectrum -> SetBorderMode(0);
  c_spectrum -> SetBorderSize(0);
  c_spectrum -> SetTickx(1);
  c_spectrum -> SetTicky(1);
  c_spectrum -> SetLeftMargin(0.18);
  c_spectrum -> SetBottomMargin(0.1518219);
  c_spectrum -> SetFrameBorderMode(0);
  c_spectrum -> SetFrameBorderMode(0);

  h_spectrum -> Draw();
  histMinv -> Draw("sameE");
  func_bck_VWG_fix -> Draw("same");
  func_tot -> Draw("same");
  func_Jpsi_CB2_fix -> Draw("same");
  //lat0 -> Draw();
  //lat1 -> Draw();
  lat2 -> Draw();
  lat3 -> Draw();
  lat4 -> Draw();
  lat5 -> Draw();
  lat6 -> Draw();
  lat7 -> Draw();

  string hist_name = outputDir + "/" + histMinv -> GetName() + ".png";

  c_spectrum -> SaveAs(hist_name.c_str());
  c_spectrum -> Close();
  delete h_spectrum;
  delete c_spectrum;
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
Double_t Func_VWG(Double_t *x, Double_t *par){
  Double_t sigma = par[2] + par[3]*((x[0] - par[1])/par[1]);
  Double_t FitBck = par[0]*TMath::Exp(-(x[0] - par[1])*(x[0] - par[1])/(2.*sigma*sigma));
  return FitBck;
}
//==============================================================================
Double_t Func_Jpsi_CB2(Double_t *x, Double_t *par){

  Double_t t = (x[0] - par[5])/par[6];
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;

  return  par[4]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Jpsi_CB2_fix(Double_t *x, Double_t *par){
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Psi2s_CB2(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[5]+(3.686-3.097)))/(par[6]*scaling_factor);
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*par[11]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*par[11]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;
  return  par[4]*par[11]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_Psi2s_CB2_fix(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[1]+(3.686-3.097)))/(par[2]*scaling_factor);
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*par[7]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*par[7]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*par[7]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}
//==============================================================================
Double_t Func_tot(Double_t *x, Double_t *par){
  //return Func_VWG(x,par) + Func_Jpsi_CB2(x,par) + Func_Psi2s_CB2(x,par);
  return Func_VWG(x,par) + Func_Jpsi_CB2(x,par);
}
