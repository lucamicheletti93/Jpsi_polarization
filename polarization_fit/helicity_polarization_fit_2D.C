#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include <TROOT.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
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

#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

double Func_W(double *, double *);
double Func_cost(double *, double *);
double Func_phi(double *, double *);

void helicity_polarization_fit_2D(int ptMin, int ptMax){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  TGaxis::SetMaxDigits(2);
  double PI = TMath::Pi();

  ostringstream convertPtMin;
  convertPtMin << ptMin;
  string strPtMin =  convertPtMin.str();

  ostringstream convertPtMax;
  convertPtMax << ptMax;
  string strPtMax =  convertPtMax.str();

  string dataset = strPtMin + "pt" + strPtMax;
  string nameOption = "_test";

  double min_fit_range_Cost = -0.7;
  double max_fit_range_Cost = 0.7;
  double min_fit_range_Phi = 0.502655;
  double max_fit_range_Phi = 2.63894;

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + nameOption + ".root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma" + nameOption + ".root";
  string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";
  string fileAccxEffName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + nameOption + ".root";

  //string filepathin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/N_Jpsi_";
  //string fileNJpsiName = filepathin + dataset + "_new.root";
  //string fileNJpsiName = filepathin + dataset + "_sigmaMC.root"; // only for 2 < pT < 6 GeV/c

  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "free_sigma.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma_test.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "_test.root";
  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileNJpsi = new TFile(fileNJpsiName.c_str(),"READ");
  TFile *fileAccxEff = new TFile(fileAccxEffName.c_str(),"READ");

  //============================================================================
  printf("---> Setting histograms bins ... \n");
  //============================================================================
  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;
  vector < vector <double> > CellAreaMatrix;

  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  CellAreaMatrix = binning -> GetCellAreaMatrix();

  const int NCostLines = NCostBins - 1;
  TLine *line_cost[NCostLines];
  for(int i = 0;i < NCostLines;i++){line_cost[i] = new TLine(CostValues[i+1],0,CostValues[i+1],PI);}

  const int NPhiLines = NPhiBins - 1;
  TLine *line_phi[NPhiLines];
  for(int i = 0;i < NPhiLines;i++){line_phi[i] = new TLine(-1,PhiValues[i+1],1,PhiValues[i+1]);}

  TH2D *h_grid = new TH2D("h_grid","",100,-1,1,50,0,PI);
  h_grid -> GetXaxis() -> SetTitle("cos#it{#theta}^{HX}");
  h_grid -> GetYaxis() -> SetTitle("#it{#varphi}^{HX}");

  //============================================================================
  printf("---> Reading the the file.root in ... \n");
  //============================================================================
  TH2D *hist_N_Jpsi_HE = (TH2D*) fileNJpsi -> Get("histNJpsi");
  //TH1D *proj_Cost_N_Jpsi_HE = (TH1D*) hist_N_Jpsi_HE -> ProjectionX("proj_Cost_N_Jpsi_HE");
  //TH1D *proj_Phi_N_Jpsi_HE = (TH1D*) hist_N_Jpsi_HE -> ProjectionY("proj_Phi_N_Jpsi_HE");
  //printf("%3.2f +- %3.2f \n",hist_N_Jpsi_HE -> GetBinContent(4,4),hist_N_Jpsi_HE -> GetBinError(4,4));

  TCanvas *c_N_Jpsi_HE = new TCanvas("c_N_Jpsi_HE","c_N_Jpsi_HE",4,132,1024,768);
  TGaxis::SetMaxDigits(2);
  h_grid -> Draw();
  hist_N_Jpsi_HE -> Draw("COLZtextsame");
  for(int i = 0;i < NCostLines;i++) line_cost[i] -> Draw("same");
  for(int i = 0;i < NPhiLines;i++) line_phi[i] -> Draw("same");

  //============================================================================
  printf("---> Reading the the Acc X Eff and projecting ... \n");
  //============================================================================
  //string fileAccxEffName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + ".root";
  //string histoAccxEffName = "hist_accxeff_HE_" + dataset + "_rebin";
  string histoAccxEffName = "hist_accxeff_HE_rebin";
  TH2D *hist_accxeff_HE_NOpol = (TH2D*) fileAccxEff -> Get(histoAccxEffName.c_str());
  TCanvas *c_accxeff_HE = new TCanvas("c_accxeff_HE","c_accxeff_HE",4,132,1024,768);
  hist_accxeff_HE_NOpol -> Draw("COLZtext");

  //TH1D *proj_Cost_accxeff_HE_NOpol = (TH1D*) fileAccxEff -> Get("proj_Cost_HE_rebin");
  //TH1D *proj_Phi_accxeff_HE_NOpol = (TH1D*) fileAccxEff -> Get("proj_Phi_HE_rebin");
  //printf("%5.4f +- %5.4f \n",hist_accxeff_HE_NOpol -> GetBinContent(4,4),hist_accxeff_HE_NOpol -> GetBinError(4,4));


  /*string filePathInput = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_";
  string file_accxeff_name = filePathInput + dataset + ".root";
  TFile *fileAccxEff = new TFile(file_accxeff_name.c_str(),"READ");
  TH2D *hist_accxeff_HE_NOpol = (TH2D*) fileAccxEff -> Get("hist_accxeff_HE_rebin");*/

  TH2D *hist_accxeff_HE_NOpol_clone = new TH2D("hist_accxeff_HE_NOpol_clone","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  for(int i = 1;i < NCostBins-1;i++){
    for(int j = 1;j < NPhiBins-1;j++){
      hist_accxeff_HE_NOpol_clone -> SetBinContent(i+1,j+1,hist_accxeff_HE_NOpol -> GetBinContent(i+1,j+1));
    }
  }

  TCanvas *c_accxeff = new TCanvas("c_accxeff","c_accxeff",4,132,1024,768);
  TGaxis::SetMaxDigits(2);
  h_grid -> Draw();
  hist_accxeff_HE_NOpol_clone -> Draw("COLZtextsame");
  for(int i = 0;i < NCostLines;i++) line_cost[i] -> Draw("same");
  for(int i = 0;i < NPhiLines;i++) line_phi[i] -> Draw("same");

  //============================================================================
  printf("---> Normalizing the Rebin histo to bin area and projecting ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *hist_N_Jpsi_HE_AN = new TH2D("hist_N_Jpsi_HE_AN","hist_N_Jpsi_HE_AN",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  //TH1D *proj_Cost_N_Jpsi_HE_AN = new TH1D("proj_Cost_N_Jpsi_HE_AN","proj_Cost_N_Jpsi_HE_AN",N_cost_bins,value_cost);
  //TH1D *proj_Phi_N_Jpsi_HE_AN = new TH1D("proj_Phi_N_Jpsi_HE_AN","proj_Phi_N_Jpsi_HE_AN",N_phi_bins,value_phi);

  for(int i = 0;i < NCostBins;i++){
    //proj_Cost_N_Jpsi_HE_AN -> SetBinContent(i+1,(proj_Cost_N_Jpsi_HE -> GetBinContent(i+1))/width_cost[i]);
    //proj_Cost_N_Jpsi_HE_AN -> SetBinError(i+1,(proj_Cost_N_Jpsi_HE -> GetBinError(i+1))/width_cost[i]);

    for(int j = 0;j < NPhiBins;j++){
      //proj_Phi_N_Jpsi_HE_AN -> SetBinContent(j+1,(proj_Phi_N_Jpsi_HE -> GetBinContent(j+1))/width_phi[j]);
      //proj_Phi_N_Jpsi_HE_AN -> SetBinError(j+1,(proj_Phi_N_Jpsi_HE -> GetBinError(j+1))/width_phi[j]);

      hist_N_Jpsi_HE_AN -> SetBinContent(i+1,j+1,(hist_N_Jpsi_HE -> GetBinContent(i+1,j+1))/CellAreaMatrix[i][j]);
      hist_N_Jpsi_HE_AN -> SetBinError(i+1,j+1,(hist_N_Jpsi_HE -> GetBinError(i+1,j+1))/CellAreaMatrix[i][j]);
    }
  }
  //printf("%3.2f +- %3.2f -> [%f] \n",hist_N_Jpsi_HE_AN -> GetBinContent(4,4),hist_N_Jpsi_HE_AN -> GetBinError(4,4),CellAreaMatrix[3][3]);

  //============================================================================
  printf("---> Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected
  //============================================================================

  TH2D *hist_N_Jpsi_HE_ANAC = new TH2D("hist_N_Jpsi_HE_ANAC","hist_N_Jpsi_HE_ANAC",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  hist_N_Jpsi_HE_ANAC -> Sumw2();
  hist_N_Jpsi_HE_ANAC -> Divide(hist_N_Jpsi_HE_AN,hist_accxeff_HE_NOpol,1,1);
  //printf("%3.2f +- %3.2f \n",hist_N_Jpsi_HE_ANAC -> GetBinContent(4,4),hist_N_Jpsi_HE_ANAC -> GetBinError(4,4));

  TCanvas *c_N_Jpsi_HE_ANAC_HE = new TCanvas("c_N_Jpsi_HE_ANAC_HE","c_N_Jpsi_HE_ANAC_HE",4,132,1024,768);
  hist_N_Jpsi_HE_ANAC -> Draw("COLZtext");

  TH2D *hist_N_Jpsi_HE_ANACX2 = (TH2D*) hist_N_Jpsi_HE_ANAC -> Clone("hist_N_Jpsi_HE_ANACX2");
  //hist_N_Jpsi_HE_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) #chi^{2} fit");
  hist_N_Jpsi_HE_ANACX2 -> SetTitle("");
  TH2D *hist_N_Jpsi_HE_ANACLW = (TH2D*) hist_N_Jpsi_HE_ANAC -> Clone("hist_N_Jpsi_HE_ANACLW");
  //hist_N_Jpsi_HE_ANACLW -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) W.L. fit");
  hist_N_Jpsi_HE_ANACLW -> SetTitle("");

  TH1D *proj2D_Cost_N_Jpsi_HE_ANACX2 = (TH1D*) hist_N_Jpsi_HE_ANAC -> ProjectionX("proj2D_Cost_N_Jpsi_HE_ANACX2");
  //proj2D_Cost_N_Jpsi_HE_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta}) #chi^{2} fit");
  proj2D_Cost_N_Jpsi_HE_ANACX2 -> SetTitle("");
  TH1D *proj2D_Cost_N_Jpsi_HE_ANACLW = (TH1D*) hist_N_Jpsi_HE_ANAC -> ProjectionX("proj2D_Cost_N_Jpsi_HE_ANACLW");
  //proj2D_Cost_N_Jpsi_HE_ANACLW -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta}) W.L. fit");
  proj2D_Cost_N_Jpsi_HE_ANACLW -> SetTitle("");

  //TH1D *proj_Cost_N_Jpsi_HE_ANAC = new TH1D("proj_Cost_N_Jpsi_HE_ANAC","proj_Cost_N_Jpsi_HE_ANAC",N_cost_bins,value_cost);
  //proj_Cost_N_Jpsi_HE_ANAC -> Sumw2();
  //proj_Cost_N_Jpsi_HE_ANAC -> Divide(proj_Cost_N_Jpsi_HE_AN,proj_Cost_accxeff_HE_NOpol,1,1);

  //TH1D *proj_Cost_N_Jpsi_HE_ANACX2 = (TH1D*) proj_Cost_N_Jpsi_HE_ANAC -> Clone("proj_Cost_N_Jpsi_HE_ANACX2");
  //TH1D *proj_Cost_N_Jpsi_HE_ANACLW = (TH1D*) proj_Cost_N_Jpsi_HE_ANAC -> Clone("proj_Cost_N_Jpsi_HE_ANACLW");

  //TH1D *proj_Phi_N_Jpsi_HE_ANAC = new TH1D("proj_Phi_N_Jpsi_HE_ANAC","proj_Phi_N_Jpsi_HE_ANAC",N_phi_bins,value_phi);
  //proj_Phi_N_Jpsi_HE_ANAC -> Sumw2();
  //proj_Phi_N_Jpsi_HE_ANAC -> Divide(proj_Phi_N_Jpsi_HE_AN,proj_Phi_accxeff_HE_NOpol,1,1);

  //TH1D *proj_Phi_N_Jpsi_HE_ANACX2 = (TH1D*) proj_Phi_N_Jpsi_HE_ANAC -> Clone("proj_Phi_N_Jpsi_HE_ANACX2");
  //TH1D *proj_Phi_N_Jpsi_HE_ANACLW = (TH1D*) proj_Phi_N_Jpsi_HE_ANAC -> Clone("proj_Phi_N_Jpsi_HE_ANACLW");

  //============================================================================
  printf("---> Fitting ... \n");
  //============================================================================

  TF2 *func2D_W_HE_X2 = new TF2("func2D_W_HE_X2",Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
  func2D_W_HE_X2 -> SetParameter(0,1000);
  func2D_W_HE_X2 -> SetParName(0,"N");
  func2D_W_HE_X2 -> SetParameter(1,0);
  func2D_W_HE_X2 -> SetParName(1,"#lambda_{#theta}");
  func2D_W_HE_X2 -> SetParameter(2,0);
  func2D_W_HE_X2 -> SetParName(2,"#lambda_{#phi}");
  func2D_W_HE_X2 -> SetParameter(3,0);
  func2D_W_HE_X2 -> SetParName(3,"#lambda_{#theta#phi}");
  hist_N_Jpsi_HE_ANACX2 -> Fit(func2D_W_HE_X2,"RSI0");

  TCanvas *c_fit2D_HE_X2 = new TCanvas("c_fit2D_HE_X2","c_fit2D_HE_X2",4,132,1024,768);
  hist_N_Jpsi_HE_ANACX2 -> Draw("COLZsame");
  func2D_W_HE_X2 -> Draw("same");

  //============================================================================

  TF1 *func1D_Cost_W_HE_X2 = new TF1("func1D_Cost_W_HE_X2",Func_cost,min_fit_range_Cost,max_fit_range_Cost,2);
  func1D_Cost_W_HE_X2 -> SetParameter(0,1000);
  func1D_Cost_W_HE_X2 -> SetParName(0,"N");
  func1D_Cost_W_HE_X2 -> SetParameter(1,0);
  func1D_Cost_W_HE_X2 -> SetParName(1,"#lambda_{#theta}");
  proj2D_Cost_N_Jpsi_HE_ANACX2 -> Fit(func1D_Cost_W_HE_X2,"RSI0");

  TCanvas *c_fit1D_Cost_HE_X2 = new TCanvas("c_fit1D_Cost_HE_X2","c_fit1D_Cost_HE_X2",4,132,1024,768);
  proj2D_Cost_N_Jpsi_HE_ANACX2 -> Draw();
  func1D_Cost_W_HE_X2 -> Draw("same");

  /*TF1 *func1D_Phi_W_HE_X2 = new TF1("func1D_Phi_W_HE_X2",Func_cost,min_fit_range_Phi,max_fit_range_Phi,2);
  func1D_Phi_W_HE_X2 -> SetParameter(0,1000);
  func1D_Phi_W_HE_X2 -> SetParName(0,"N");
  func1D_Phi_W_HE_X2 -> SetParameter(1,0);
  func1D_Phi_W_HE_X2 -> SetParName(1,"#lambda_{#theta}");
  proj_Phi_N_Jpsi_HE_ANACX2 -> Fit(func1D_Phi_W_HE_X2,"RSI0");

  TCanvas *c_fit1D_HE_X2 = new TCanvas("c_fit1D_HE_X2","c_fit1D_HE_X2",4,132,1024,768);
  proj_Phi_N_Jpsi_HE_ANACX2 -> Draw("COLZ");
  func1D_Phi_W_HE_X2 -> Draw("same");*/

  //============================================================================

  /*TF2 *func2D_W_HE_LW = new TF2("func2D_W_HE_LW",Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
  func2D_W_HE_LW -> SetParameter(0,1000);
  func2D_W_HE_LW -> SetParName(0,"N");
  func2D_W_HE_LW -> SetParameter(1,0);
  func2D_W_HE_LW -> SetParName(1,"#lambda_{#theta}");
  func2D_W_HE_LW -> SetParameter(2,0);
  func2D_W_HE_LW -> SetParName(2,"#lambda_{#phi}");
  func2D_W_HE_LW -> SetParameter(3,0);
  func2D_W_HE_LW -> SetParName(3,"#lambda_{#theta#phi}");
  hist_N_Jpsi_HE_ANACLW -> Fit(func2D_W_HE_LW,"RSI0WL");

  TCanvas *c_fit2D_HE_LW = new TCanvas("c_fit2D_HE_LW","c_fit2D_HE_LW",4,132,1024,768);
  hist_N_Jpsi_HE_ANACLW -> Draw("COLZ");
  func2D_W_HE_LW -> Draw("same");*/

  /*double central_binx;
  double central_biny;

  TH2D *histo_diff_fit_LW = new TH2D("histo_diff_fit_LW","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      central_binx = hist_N_Jpsi_HE_ANACLW -> GetXaxis() -> GetBinCenter(i+1,j+1);
      central_biny = hist_N_Jpsi_HE_ANACLW -> GetYaxis() -> GetBinCenter(i+1,j+1);
      histo_diff_fit_LW -> SetBinContent(i+1,j+1,hist_N_Jpsi_HE_ANACLW -> GetBinContent(i+1,j+1)/func2D_W_HE_LW -> Eval(central_binx,central_biny));
    }
  }*/

  //============================================================================

  /*TF1 *func1D_Cost_W_HE_LW = new TF1("func1D_Cost_W_HE_LW",Func_cost,min_fit_range_Cost,max_fit_range_Cost,2);
  func1D_Cost_W_HE_LW -> SetParameter(0,1000);
  func1D_Cost_W_HE_LW -> SetParName(0,"N");
  func1D_Cost_W_HE_LW -> SetParameter(1,0);
  func1D_Cost_W_HE_LW -> SetParName(1,"#lambda_{#theta}");
  proj2D_Cost_N_Jpsi_HE_ANACLW -> Fit(func1D_Cost_W_HE_LW,"RSI0WL");

  TCanvas *c_fit1D_2pt6_HE_LW = new TCanvas("c_fit1D_2pt6_HE_LW","c_fit1D_2pt6_HE_LW",4,132,1024,768);
  //proj_Cost_N_Jpsi_2pt6_HE_ANACLW -> Draw();
  proj2D_Cost_N_Jpsi_HE_ANACLW -> Draw();
  func1D_Cost_W_HE_LW -> Draw("same");*/

}
////////////////////////////////////////////////////////////////////////////////
double Func_W(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];
  double L_thph = par[3];

  double costh = x[0];
  double phi = x[1];

  double cosph = TMath::Cos(phi);
  double cos2ph = TMath::Cos(2*phi);

  double W =  (N/(3 + L_th))*(1 + (L_th + L_ph*cos2ph)*costh*costh + 2*L_thph*costh*cosph*TMath::Sqrt(1 - costh*costh) + L_ph*cos2ph);
  return W;
}
//------------------------------------------------------------------------------
double Func_cost(double *x, double *par){

  double N = par[0];
  double L_th = par[1];

  double costh = x[0];

  double W =  (N/(3 + L_th))*(1 + L_th*costh*costh);
  return W;
}
//------------------------------------------------------------------------------
double Func_phi(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];

  double phi = x[0];

  double cos2ph = TMath::Cos(2*phi);

  double W =  N*(1 + ((2*L_ph)/(3 + L_th))*cos2ph);
  return W;
}
