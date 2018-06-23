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
void fit_of_minv(TH1D *, int, int, string);

double PI = TMath::Pi();
Double_t scaling_factor = 1.05154; //factor introduced to pass from the sigma of Jpsi to the sigma of Psi(2S)
double n_jpsi = 0, stat_jpsi = 0;
double ChiSquare_NDF = 0;
double sigma_Jpsi = 0;
char fit_status[10];
vector <double> CostValues;
vector <double> PhiValues;

//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void single_histo_fit(){

  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  string dataset = "4pt7";
  string filepathin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/mass_histos_cost_phi_";
  string filenamein = filepathin + dataset + ".root";
  TFile *file = new TFile(filenamein.c_str());

  string filebinningname = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";
  TFile *filebinning = new TFile(filebinningname.c_str(),"READ");
  Binning *binning = (Binning*) filebinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  TH1D *hist_minv = (TH1D*) file -> Get("HE_11cost20_21phi23");
  hist_minv -> SetMarkerStyle(20);
  hist_minv -> SetMarkerSize(0.7);
  fit_of_minv(hist_minv,1,1,dataset);
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
void loop_on_histos(){

  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  string dataset = "6pt10";
  string filepathin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/mass_histos_cost_phi_";
  string filenamein = filepathin + dataset + ".root";
  TFile *file = new TFile(filenamein.c_str());

  //============================================================================
  printf("---> Setting histograms bins ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);

  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  string filebinningname = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";
  TFile *filebinning = new TFile(filebinningname.c_str(),"READ");
  Binning *binning = (Binning*) filebinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  //============================================================================

  bool save_tree = kTRUE;
  int counter_cost = 0; // counter to fill the matrix
  int counter_phi = 0; // counter to fill the matrix

  //TH1D *hist_minv_integrated = new TH1D("hist_minv_integrated","hist_minv_integrated",120,2,5);

  //============================================================================
  char hist_mass_HE_name[100];
  vector <string> fit_errors;

  TH2D *hist_N_Jpsi_HE = new TH2D("hist_N_Jpsi_HE","hist_N_Jpsi_HE",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      sprintf(hist_mass_HE_name,"HE_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      TH1D *hist_minv = (TH1D*) file -> Get(hist_mass_HE_name);
      //hist_minv_integrated -> Add(hist_minv);
      if(i != 0 && i != NCostBins-1 &&  j != 0 && j != NPhiBins-1) fit_of_minv(hist_minv,i,j,dataset);
      else{n_jpsi = 0; stat_jpsi = 0;}
      hist_N_Jpsi_HE -> SetBinContent(i+1,j+1,(int) n_jpsi);
      hist_N_Jpsi_HE -> SetBinError(i+1,j+1,(int) stat_jpsi);

      if(ChiSquare_NDF > 3. || sigma_Jpsi > 0.09) fit_errors.push_back(hist_mass_HE_name);
    }
  }

  hist_N_Jpsi_HE -> Draw("COLZtexterror");

  //============================================================================
  printf("List of bad fits ... \n");
  //============================================================================
  for(int i = 0;i < (int) fit_errors.size();i++) cout << fit_errors[i] << endl;
  //============================================================================
  printf("Saving results ... \n");
  //============================================================================

  if(save_tree){
    string filepathout = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_";
    //string filenameout = filepathout + dataset + "_sigmaMC.root";
    string filenameout = filepathout + dataset + ".root";

    TFile *output_file = new TFile(filenameout.c_str(),"RECREATE");
    hist_N_Jpsi_HE -> Write();
    output_file -> Close();
    printf("The file is saved in %s \n",filenameout.c_str());
  }

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void fit_of_minv(TH1D *hist_minv, int counter_cost, int counter_phi, string dataset){

  string sig = "CB2";
  string bck = "VWG";
  string fit_range = "min";
  double min_fit_range = 2.;
  double max_fit_range = 5.;
  string tails = "pp8";
  string tails_fix = "yes";

  gStyle -> SetOptStat(0);
  TGaxis::SetMaxDigits(2);

  if(dataset.compare("2pt6") == 0){
    if(counter_cost == 1 || counter_cost == (int) CostValues.size() - 2){hist_minv -> Rebin(2);}
  }
  if(dataset.compare("6pt12") == 0){hist_minv -> Rebin(2);}
  if(dataset.compare("6pt10") == 0){hist_minv -> Rebin(2);}
  if(dataset.compare("4pt7") == 0){hist_minv -> Rebin(2);}
  if(dataset.compare("7pt10") == 0){hist_minv -> Rebin(2);}
  if(dataset.compare("4pt6") == 0){hist_minv -> Rebin(2);}

  //============================================================================
  // FIT OF THE BACKGROUND
  //============================================================================
  TH1D *histo_for_bck = (TH1D*) hist_minv -> Clone("histo_for_bck");
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

  //double par_bck[4] = {500000.,0.6,0.4,0.2};
  double par_bck[4] = {500000.,0.6,0.2,0.2};
  double *par_bck_ptr = par_bck;

  //TF1 *func_bck_VWG = new TF1("func_bck_VWG",Func_VWG,2.2,5.5,4);
  TF1 *func_bck_VWG = new TF1("func_bck_VWG",Func_VWG,2.1,4.9,4);
  func_bck_VWG -> SetLineColor(kRed);
  func_bck_VWG -> SetParameters(par_bck);
  histo_for_bck -> Fit(func_bck_VWG,"R0");

  //============================================================================
  // TAILS {alpha1,n1,alpha2,n2}
  //============================================================================
  //double tails_par[4] = {1.012,2.874,2.232,2.924}; //p-A MC COLLISIONS
  //double tails_par[4] = {1.089,3.393,2.097,8.694}; //p-p 8 TeV
  //double tails_par[4] = {0.98,6.97,1.86,14.99}; //p-p 13 TeV
  double tails_par[4] = {1.06,3.23,2.55,1.56}; //GEANT 4
  //double tails_par[4] = {0.97,3.98,2.3,3.03}; //GEANT 3
  //============================================================================
  // FIT OF THE JPSI SIGNAL
  //============================================================================
  //double par_signal[8] = {50000.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01};
  //double par_signal[8] = {500.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01}; // fot 4 < pT < 7 GeV/c
  double par_signal[8] = {50.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01}; // fot 7 < pT < 10 GeV/c

  TF1 *func_Jpsi_CB2 = new TF1("func_Jpsi_CB2",Func_Jpsi_CB2,2.9,3.2,11);
  func_Jpsi_CB2 -> SetNpx(100000);
  func_Jpsi_CB2 -> SetParameter(4,par_signal[0]);
  func_Jpsi_CB2 -> SetParameter(5,par_signal[1]);
  func_Jpsi_CB2 -> FixParameter(6,par_signal[2]);
  func_Jpsi_CB2 -> FixParameter(7,tails_par[0]);
  func_Jpsi_CB2 -> FixParameter(8,tails_par[1]);
  func_Jpsi_CB2 -> FixParameter(9,tails_par[2]);
  func_Jpsi_CB2 -> FixParameter(10,tails_par[3]);
  hist_minv -> Fit(func_Jpsi_CB2,"R0");

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
      func_tot -> SetParLimits(5,2.7,3.3);
      func_tot -> SetParameter(6,7.0e-02);
      //func_tot -> SetParLimits(6,5.0e-02,1.2e-01);
      func_tot -> SetParLimits(6,5.0e-02,9.0e-02); // To be resetted
      //========================================================================
      // FIX THE WIDTH TO MC
      //char INPUT_FILE_NAME[300] = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/accxeff_new.root";
      //TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
      //TH2D *histo_sigma_Jpsi = (TH2D*) input_file -> Get("histo_sigma_Jpsi");
      //double sigma_Jpsi_MC = histo_sigma_Jpsi -> GetBinContent(counter_cost+1,counter_phi+1);
      //func_tot -> FixParameter(6,sigma_Jpsi_MC);
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
    fit_ptr = (TFitResultPtr) hist_minv -> Fit(func_tot,"RLS0");
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
    TH2D *histo_sigma_Jpsi = (TH2D*) input_file -> Get("histo_sigma_Jpsi");
    double sigma_Jpsi_MC = histo_sigma_Jpsi -> GetBinContent(counter_cost,counter_phi);
    func_tot -> FixParameter(6,6.5e-02);
    fit_ptr = (TFitResultPtr) hist_minv -> Fit(func_tot,"RLS0");
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
  sigma_Jpsi = func_tot -> GetParameter(6);
  double sigma_min_Jpsi = func_tot -> GetParameter(5) - 3*(func_tot -> GetParameter(6));
  double sigma_max_Jpsi = func_tot -> GetParameter(5) + 3*(func_tot -> GetParameter(6));
  double N_Jpsi_3sigma = func_Jpsi_CB2_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double N_bck_Jpsi_3sigma = func_bck_VWG_fix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double SB_Jpsi = N_Jpsi_3sigma/N_bck_Jpsi_3sigma;
  n_jpsi = func_Jpsi_CB2_fix -> Integral(0,5)/m_width;
  stat_jpsi = func_Jpsi_CB2_fix -> IntegralError(0.,5.,Jpsi_par,Jpsi_mat)/m_width;

  //============================================================================
  // PLOT OF THE TOTAL SPECTRUM
  //============================================================================
  char title[100];
  double max_histo_value = hist_minv -> GetMaximum();
  max_histo_value = max_histo_value + 0.6*max_histo_value;

  TH2D *h_spectrum = new TH2D("h_spectrum","",120,2,5,100,0,max_histo_value);
  h_spectrum -> GetXaxis() -> SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  h_spectrum -> GetXaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetXaxis() -> SetTitleOffset(0.95);
  h_spectrum -> GetXaxis() -> SetLabelSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitle("Counts per 25 MeV/#it{c}^{2}");
  h_spectrum -> GetYaxis() -> CenterTitle(true);
  h_spectrum -> GetYaxis() -> SetTitleSize(0.05);
  h_spectrum -> GetYaxis() -> SetTitleOffset(1.2);
  h_spectrum -> GetYaxis() -> SetLabelSize(0.05);

  sprintf(title,"N_{J/#psi} = %2.0f #pm %2.0f",n_jpsi,stat_jpsi);
  TLatex *lat2 = new TLatex(0.55,0.82,title);
  lat2 -> SetTextSize(0.04);
  lat2 -> SetNDC();
  lat2 -> SetTextFont(42);

  sprintf(title,"#it{m}_{J/#psi} = %4.3f GeV/#it{c}^{2}",mass_Jpsi);
  TLatex *lat3 = new TLatex(0.55,0.76,title);
  lat3 -> SetTextSize(0.04);
  lat3 -> SetNDC();
  lat3 -> SetTextFont(42);

  sprintf(title,"#it{#sigma}_{J/#psi} = %3.2f MeV/#it{c}^{2}",sigma_Jpsi*1000);
  TLatex *lat4 = new TLatex(0.55,0.70,title);
  lat4 -> SetTextSize(0.04);
  lat4 -> SetNDC();
  lat4 -> SetTextFont(42);

  sprintf(title,"%3.2f < cos#it{#theta}^{HX} < %3.2f",CostValues[counter_cost],CostValues[counter_cost+1]);
  if(counter_cost == 100) sprintf(title,"%2.1f < cos#it{#theta}^{HX} < %2.1f",-1.,1.);
  TLatex *lat5 = new TLatex(0.55,0.62,title);
  lat5 -> SetTextSize(0.04);
  lat5 -> SetNDC();
  lat5 -> SetTextFont(42);

  sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",PhiValues[counter_phi],PhiValues[counter_phi+1]);
  if(counter_phi == 100) sprintf(title,"%3.2f < #it{#varphi}^{HX} < %3.2f rad",0.,PI);
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
  hist_minv -> Draw("sameE");
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

  string hist_name = "../fit_plots/" + dataset + "/" + hist_minv -> GetName() + ".png";
  c_spectrum -> SaveAs(hist_name.c_str());
  c_spectrum -> Close();
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
