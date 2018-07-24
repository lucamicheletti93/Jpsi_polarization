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

  double min_fit_range_Cost = -0.8;
  double max_fit_range_Cost = 0.8;
  double min_fit_range_Phi = 0.502655;
  double max_fit_range_Phi = 2.63894;

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + nameOption + ".root";
  string fileNJpsiSigmaFixedName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "_fixed_sigma" + nameOption + ".root";
  string fileNJpsiSigmaFreeName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";
  string fileAccxEffName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + nameOption + ".root";

  //string filepathin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/N_Jpsi_";
  //string fileNJpsiName = filepathin + dataset + "_new.root";
  //string fileNJpsiName = filepathin + dataset + "_sigmaMC.root"; // only for 2 < pT < 6 GeV/c

  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "free_sigma.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma_test.root";
  //string fileNJpsiName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "_test.root";
  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileNJpsiSigmaFixed = new TFile(fileNJpsiSigmaFixedName.c_str(),"READ");
  TFile *fileNJpsiSigmaFree = new TFile(fileNJpsiSigmaFreeName.c_str(),"READ");
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
  TH2D *histNJpsiSigmaFixed = (TH2D*) fileNJpsiSigmaFixed -> Get("histNJpsi");
  TH2D *histNJpsiSigmaFree = (TH2D*) fileNJpsiSigmaFree -> Get("histNJpsi");
  //TH1D *proj_Cost_N_Jpsi_HE = (TH1D*) histNJpsiSigmaFixed -> ProjectionX("proj_Cost_N_Jpsi_HE");
  //TH1D *proj_Phi_N_Jpsi_HE = (TH1D*) histNJpsiSigmaFixed -> ProjectionY("proj_Phi_N_Jpsi_HE");
  //printf("%3.2f +- %3.2f \n",histNJpsiSigmaFixed -> GetBinContent(4,4),histNJpsiSigmaFixed -> GetBinError(4,4));

  TCanvas *c_N_Jpsi_HE = new TCanvas("c_N_Jpsi_HE","c_N_Jpsi_HE",4,132,1024,768);
  TGaxis::SetMaxDigits(2);
  h_grid -> Draw();
  histNJpsiSigmaFixed -> Draw("COLZtextsame");
  for(int i = 0;i < NCostLines;i++) line_cost[i] -> Draw("same");
  for(int i = 0;i < NPhiLines;i++) line_phi[i] -> Draw("same");

  printf("---> Integral : %i \n",(int) histNJpsiSigmaFixed -> Integral());

  //============================================================================
  printf("---> Reading the the Acc X Eff and projecting ... \n");
  //============================================================================
  //string fileAccxEffName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + ".root";
  //string histAccxEffName = "hist_accxeff_HE_" + dataset + "_rebin";
  string histAccxEffName = "hist_accxeff_HE_rebin";
  TH2D *histAccxEff = (TH2D*) fileAccxEff -> Get(histAccxEffName.c_str());
  TCanvas *c_accxeff_HE = new TCanvas("c_accxeff_HE","c_accxeff_HE",4,132,1024,768);
  histAccxEff -> Draw("COLZtext");

  //TH1D *proj_Cost_accxeff_HE_NOpol = (TH1D*) fileAccxEff -> Get("proj_Cost_HE_rebin");
  //TH1D *proj_Phi_accxeff_HE_NOpol = (TH1D*) fileAccxEff -> Get("proj_Phi_HE_rebin");
  //printf("%5.4f +- %5.4f \n",histAccxEff -> GetBinContent(4,4),histAccxEff -> GetBinError(4,4));


  /*string filePathInput = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_";
  string file_accxeff_name = filePathInput + dataset + ".root";
  TFile *fileAccxEff = new TFile(file_accxeff_name.c_str(),"READ");
  TH2D *histAccxEff = (TH2D*) fileAccxEff -> Get("hist_accxeff_HE_rebin");*/

  TH2D *histAccxEffClone = new TH2D("histAccxEffClone","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  for(int i = 1;i < NCostBins-1;i++){
    for(int j = 1;j < NPhiBins-1;j++){
      histAccxEffClone -> SetBinContent(i+1,j+1,histAccxEff -> GetBinContent(i+1,j+1));
    }
  }

  TCanvas *c_accxeff = new TCanvas("c_accxeff","c_accxeff",4,132,1024,768);
  TGaxis::SetMaxDigits(2);
  h_grid -> Draw();
  histAccxEffClone -> Draw("COLZtextsame");
  for(int i = 0;i < NCostLines;i++) line_cost[i] -> Draw("same");
  for(int i = 0;i < NPhiLines;i++) line_phi[i] -> Draw("same");

  //============================================================================
  printf("---> Normalizing the Rebin histo to bin area and projecting ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *histNJpsiSigmaFixed_AN = new TH2D("histNJpsiSigmaFixed_AN","histNJpsiSigmaFixed_AN",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *histNJpsiSigmaFree_AN = new TH2D("histNJpsiSigmaFree_AN","histNJpsiSigmaFree_AN",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  //TH1D *proj_Cost_N_Jpsi_HE_AN = new TH1D("proj_Cost_N_Jpsi_HE_AN","proj_Cost_N_Jpsi_HE_AN",N_cost_bins,value_cost);
  //TH1D *proj_Phi_N_Jpsi_HE_AN = new TH1D("proj_Phi_N_Jpsi_HE_AN","proj_Phi_N_Jpsi_HE_AN",N_phi_bins,value_phi);

  for(int i = 0;i < NCostBins;i++){
    //proj_Cost_N_Jpsi_HE_AN -> SetBinContent(i+1,(proj_Cost_N_Jpsi_HE -> GetBinContent(i+1))/width_cost[i]);
    //proj_Cost_N_Jpsi_HE_AN -> SetBinError(i+1,(proj_Cost_N_Jpsi_HE -> GetBinError(i+1))/width_cost[i]);

    for(int j = 0;j < NPhiBins;j++){
      //proj_Phi_N_Jpsi_HE_AN -> SetBinContent(j+1,(proj_Phi_N_Jpsi_HE -> GetBinContent(j+1))/width_phi[j]);
      //proj_Phi_N_Jpsi_HE_AN -> SetBinError(j+1,(proj_Phi_N_Jpsi_HE -> GetBinError(j+1))/width_phi[j]);

      histNJpsiSigmaFixed_AN -> SetBinContent(i+1,j+1,(histNJpsiSigmaFixed -> GetBinContent(i+1,j+1))/CellAreaMatrix[i][j]);
      histNJpsiSigmaFixed_AN -> SetBinError(i+1,j+1,(histNJpsiSigmaFixed -> GetBinError(i+1,j+1))/CellAreaMatrix[i][j]);

      histNJpsiSigmaFree_AN -> SetBinContent(i+1,j+1,(histNJpsiSigmaFree -> GetBinContent(i+1,j+1))/CellAreaMatrix[i][j]);
      histNJpsiSigmaFree_AN -> SetBinError(i+1,j+1,(histNJpsiSigmaFree -> GetBinError(i+1,j+1))/CellAreaMatrix[i][j]);
    }
  }
  //printf("%3.2f +- %3.2f -> [%f] \n",histNJpsiSigmaFixed_AN -> GetBinContent(4,4),histNJpsiSigmaFixed_AN -> GetBinError(4,4),CellAreaMatrix[3][3]);

  //============================================================================
  printf("---> Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected
  //============================================================================

  TH2D *histNJpsiSigmaFixed_ANAC = new TH2D("histNJpsiSigmaFixed_ANAC","histNJpsiSigmaFixed_ANAC",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histNJpsiSigmaFixed_ANAC -> Sumw2();
  histNJpsiSigmaFixed_ANAC -> Divide(histNJpsiSigmaFixed_AN,histAccxEff,1,1);
  //printf("%3.2f +- %3.2f \n",histNJpsiSigmaFixed_ANAC -> GetBinContent(4,4),histNJpsiSigmaFixed_ANAC -> GetBinError(4,4));

  TH2D *histNJpsiSigmaFree_ANAC = new TH2D("histNJpsiSigmaFree_ANAC","histNJpsiSigmaFree_ANAC",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histNJpsiSigmaFree_ANAC -> Sumw2();
  histNJpsiSigmaFree_ANAC -> Divide(histNJpsiSigmaFree_AN,histAccxEff,1,1);

  TCanvas *c_N_Jpsi_HE_ANAC_HE = new TCanvas("c_N_Jpsi_HE_ANAC_HE","c_N_Jpsi_HE_ANAC_HE",4,132,1024,768);
  histNJpsiSigmaFixed_ANAC -> Draw("COLZtext");

  TH2D *histNJpsiSigmaFixed_ANACX2 = (TH2D*) histNJpsiSigmaFixed_ANAC -> Clone("histNJpsiSigmaFixed_ANACX2");
  //histNJpsiSigmaFixed_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) #chi^{2} fit");
  histNJpsiSigmaFixed_ANACX2 -> SetTitle("");
  TH2D *histNJpsiSigmaFixed_ANACLW = (TH2D*) histNJpsiSigmaFixed_ANAC -> Clone("histNJpsiSigmaFixed_ANACLW");
  //histNJpsiSigmaFixed_ANACLW -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) W.L. fit");
  histNJpsiSigmaFixed_ANACLW -> SetTitle("");

  TH2D *histNJpsiSigmaFree_ANACX2 = (TH2D*) histNJpsiSigmaFree_ANAC -> Clone("histNJpsiSigmaFree_ANACX2");
  //histNJpsiSigmaFree_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) #chi^{2} fit");
  histNJpsiSigmaFree_ANACX2 -> SetTitle("");
  TH2D *histNJpsiSigmaFree_ANACLW = (TH2D*) histNJpsiSigmaFree_ANAC -> Clone("histNJpsiSigmaFree_ANACLW");
  //histNJpsiSigmaFree_ANACLW -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta} #upoint #Delta#it{#varphi}) W.L. fit");
  histNJpsiSigmaFree_ANACLW -> SetTitle("");

  TH1D *histNJpsiCostSigmaFixed_ANACX2 = (TH1D*) histNJpsiSigmaFixed_ANAC -> ProjectionX("histNJpsiCostSigmaFixed_ANACX2");
  //histNJpsiCostSigmaFixed_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta}) #chi^{2} fit");
  histNJpsiCostSigmaFixed_ANACX2 -> SetTitle("");
  histNJpsiCostSigmaFixed_ANACX2 -> SetLineColor(kBlue);

  TH1D *histNJpsiCostSigmaFree_ANACX2 = (TH1D*) histNJpsiSigmaFree_ANAC -> ProjectionX("histNJpsiCostSigmaFree_ANACX2");
  //histNJpsiCostSigmaFree_ANACX2 -> SetTitle("N_{J/#psi}/(A#times#epsilon #upoint #Deltacos#it{#theta}) #chi^{2} fit");
  histNJpsiCostSigmaFree_ANACX2 -> SetTitle("");
  histNJpsiCostSigmaFree_ANACX2 -> SetLineColor(kRed);

  TH1D *proj2D_Cost_N_Jpsi_HE_ANACLW = (TH1D*) histNJpsiSigmaFixed_ANAC -> ProjectionX("proj2D_Cost_N_Jpsi_HE_ANACLW");
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

  TF2 *func2DSigmaFixedX2 = new TF2("func2DSigmaFixedX2",Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
  func2DSigmaFixedX2 -> SetParameter(0,1000);
  func2DSigmaFixedX2 -> SetParName(0,"N");
  func2DSigmaFixedX2 -> SetParameter(1,0);
  func2DSigmaFixedX2 -> SetParName(1,"#lambda_{#theta}");
  func2DSigmaFixedX2 -> SetParameter(2,0);
  func2DSigmaFixedX2 -> SetParName(2,"#lambda_{#phi}");
  func2DSigmaFixedX2 -> SetParameter(3,0);
  func2DSigmaFixedX2 -> SetParName(3,"#lambda_{#theta#phi}");
  histNJpsiSigmaFixed_ANACX2 -> Fit(func2DSigmaFixedX2,"RSI0");

  TCanvas *cFit2DSigmaFixedX2 = new TCanvas("cFit2DSigmaFixedX2","cFit2DSigmaFixedX2",4,132,1024,768);
  histNJpsiSigmaFixed_ANACX2 -> Draw("COLZsame");
  func2DSigmaFixedX2 -> Draw("same");

  TF2 *func2DSigmaFreeX2 = new TF2("func2DSigmaFreeX2",Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
  func2DSigmaFreeX2 -> SetParameter(0,1000);
  func2DSigmaFreeX2 -> SetParName(0,"N");
  func2DSigmaFreeX2 -> SetParameter(1,0);
  func2DSigmaFreeX2 -> SetParName(1,"#lambda_{#theta}");
  func2DSigmaFreeX2 -> SetParameter(2,0);
  func2DSigmaFreeX2 -> SetParName(2,"#lambda_{#phi}");
  func2DSigmaFreeX2 -> SetParameter(3,0);
  func2DSigmaFreeX2 -> SetParName(3,"#lambda_{#theta#phi}");
  histNJpsiSigmaFree_ANACX2 -> Fit(func2DSigmaFreeX2,"RSI0");

  TCanvas *cFit2DSigmaFreeX2 = new TCanvas("cFit2DSigmaFreeX2","cFit2DSigmaFreeX2",4,132,1024,768);
  histNJpsiSigmaFree_ANACX2 -> Draw("COLZsame");
  func2DSigmaFreeX2 -> Draw("same");

  //============================================================================

  TF2 *func2DCostSigmaFixed[3], *func2DCostSigmaFree[3];
  TF1 *func1DCostSigmaFixed[3], *func1DCostSigmaFree[3];
  double fitRangeMin[3] = {-0.6,-0.7,-0.8};
  double fitRangeMax[3] = {0.6,0.7,0.8};
  TLatex *latLambdaTheta[3];
  double LTh2DSigmaFixed, errLTh2DSigmaFixed, LTh2DSigmaFree, errLTh2DSigmaFree;
  double LTh1DSigmaFixed, errLTh1DSigmaFixed, LTh1DSigmaFree, errLTh1DSigmaFree;

  for(int i = 0;i < 3;i++){
    func2DCostSigmaFixed[i] = new TF2(Form("func2DCostSigmaFixed%i",i),Func_W,fitRangeMin[i],fitRangeMax[i],min_fit_range_Phi,max_fit_range_Phi,4);
    func2DCostSigmaFixed[i] -> SetParameter(0,1000);
    func2DCostSigmaFixed[i] -> SetParName(0,"N");
    func2DCostSigmaFixed[i] -> SetParameter(1,0);
    func2DCostSigmaFixed[i] -> SetParName(1,"#lambda_{#theta}");
    func2DCostSigmaFixed[i] -> SetParameter(2,0);
    func2DCostSigmaFixed[i] -> SetParName(2,"#lambda_{#phi}");
    func2DCostSigmaFixed[i] -> SetParameter(3,0);
    func2DCostSigmaFixed[i] -> SetParName(3,"#lambda_{#theta#phi}");
    histNJpsiSigmaFixed_ANACX2 -> Fit(func2DCostSigmaFixed[i],"RSI0");

    LTh2DSigmaFixed = func2DCostSigmaFixed[i] -> GetParameter(1);
    errLTh2DSigmaFixed = func2DCostSigmaFixed[i] -> GetParError(1);

    func1DCostSigmaFixed[i] = new TF1(Form("func1DCostSigmaFixed%i",i),Func_cost,fitRangeMin[i],fitRangeMax[i],2);
    func1DCostSigmaFixed[i] -> SetParameter(0,1000);
    func1DCostSigmaFixed[i] -> SetParName(0,"N");
    func1DCostSigmaFixed[i] -> SetParameter(1,0);
    func1DCostSigmaFixed[i] -> SetParName(1,"#lambda_{#theta}");
    func1DCostSigmaFixed[i] -> SetLineColor(kBlue);
    func1DCostSigmaFixed[i] -> SetLineStyle(i+1);
    histNJpsiCostSigmaFixed_ANACX2 -> Fit(func1DCostSigmaFixed[i],"RSI0");

    LTh1DSigmaFixed = func1DCostSigmaFixed[i] -> GetParameter(1);
    errLTh1DSigmaFixed = func1DCostSigmaFixed[i] -> GetParError(1);

    func2DCostSigmaFree[i] = new TF2(Form("func2DCostSigmaFree%i",i),Func_W,fitRangeMin[i],fitRangeMax[i],min_fit_range_Phi,max_fit_range_Phi,4);
    func2DCostSigmaFree[i] -> SetParameter(0,1000);
    func2DCostSigmaFree[i] -> SetParName(0,"N");
    func2DCostSigmaFree[i] -> SetParameter(1,0);
    func2DCostSigmaFree[i] -> SetParName(1,"#lambda_{#theta}");
    func2DCostSigmaFree[i] -> SetParameter(2,0);
    func2DCostSigmaFree[i] -> SetParName(2,"#lambda_{#phi}");
    func2DCostSigmaFree[i] -> SetParameter(3,0);
    func2DCostSigmaFree[i] -> SetParName(3,"#lambda_{#theta#phi}");
    histNJpsiSigmaFree_ANACX2 -> Fit(func2DCostSigmaFree[i],"RSI0");

    LTh2DSigmaFree = func2DCostSigmaFree[i] -> GetParameter(1);
    errLTh2DSigmaFree = func2DCostSigmaFree[i] -> GetParError(1);

    func1DCostSigmaFree[i] = new TF1(Form("func1DCostSigmaFree%i",i),Func_cost,fitRangeMin[i],fitRangeMax[i],2);
    func1DCostSigmaFree[i] -> SetParameter(0,1000);
    func1DCostSigmaFree[i] -> SetParName(0,"N");
    func1DCostSigmaFree[i] -> SetParameter(1,0);
    func1DCostSigmaFree[i] -> SetParName(1,"#lambda_{#theta}");
    func1DCostSigmaFree[i] -> SetLineColor(kRed);
    func1DCostSigmaFree[i] -> SetLineStyle(i+1);
    histNJpsiCostSigmaFree_ANACX2 -> Fit(func1DCostSigmaFree[i],"RSI0");

    LTh1DSigmaFree = func1DCostSigmaFree[i] -> GetParameter(1);
    errLTh1DSigmaFree = func1DCostSigmaFree[i] -> GetParError(1);

    latLambdaTheta[i] = new TLatex(0.18,0.15 + i*0.06,Form("[%2.1f<cos#theta<%2.1f] #color[4]{#lambda_{#theta} = %4.3f #pm %4.3f} (#lambda_{#theta}^{2D} = %4.3f #pm %4.3f) ; #color[2]{#lambda_{#theta} = %4.3f #pm %4.3f} (#lambda_{#theta}^{2D} = %4.3f #pm %4.3f)",fitRangeMin[i],fitRangeMax[i],LTh1DSigmaFixed,errLTh1DSigmaFixed,LTh2DSigmaFixed,errLTh2DSigmaFixed,LTh1DSigmaFree,errLTh1DSigmaFree,LTh2DSigmaFree,errLTh2DSigmaFree));
    latLambdaTheta[i] -> SetTextSize(0.02);
    latLambdaTheta[i] -> SetNDC();
    latLambdaTheta[i] -> SetTextFont(42);
  }

  TLegend *legNJpsiCost = new TLegend(0.42,0.7,0.65,0.8,"","brNDC");
  legNJpsiCost -> SetBorderSize(0);
  legNJpsiCost -> SetFillColor(10);
  legNJpsiCost -> SetFillStyle(1);
  legNJpsiCost -> SetLineStyle(0);
  legNJpsiCost -> SetLineColor(0);
  legNJpsiCost -> SetTextFont(42);
  legNJpsiCost -> SetTextSize(0.035);
  legNJpsiCost -> AddEntry(histNJpsiCostSigmaFree_ANACX2,"Free #sigma_{J/#psi}","L");
  legNJpsiCost -> AddEntry(histNJpsiCostSigmaFixed_ANACX2,"Fixed #sigma_{J/#psi}","L");

  double maximumNJpsiCost = histNJpsiCostSigmaFree_ANACX2 -> GetBinContent(histNJpsiCostSigmaFree_ANACX2 -> FindBin(0));

  TCanvas *cFit1DCost = new TCanvas("cFit1DCost","cFit1DCost",4,132,1024,768);
  TH2D *hNJpsiCost = new TH2D("hNJpsiCost","",100,-1,1,100,maximumNJpsiCost - 0.8*maximumNJpsiCost,maximumNJpsiCost + 0.8*maximumNJpsiCost);
  hNJpsiCost -> Draw();
  legNJpsiCost -> Draw("same");
  histNJpsiCostSigmaFixed_ANACX2 -> Draw("same");
  histNJpsiCostSigmaFree_ANACX2 -> Draw("same");
  for(int i = 0;i < 3;i++){
    func1DCostSigmaFixed[i] -> Draw("same");
    func1DCostSigmaFree[i] -> Draw("same");
    latLambdaTheta[i] -> Draw("same");
  }

  cFit1DCost -> SaveAs(Form("/home/luca/Scrivania/%s.png",dataset.c_str()));

  return;
  /*TCanvas *cFit1DCostSigmaFixedX2 = new TCanvas("cFit1DCostSigmaFixedX2","cFit1DCostSigmaFixedX2",4,132,1024,768);
  hNJpsiCost -> Draw();
  histNJpsiCostSigmaFixed_ANACX2 -> Draw("same");
  func1DCostSigmaFixedX2 -> Draw("same");
  //legNJpsiCost -> Draw("same");


  TCanvas *cFit1DCostSigmaFreeX2 = new TCanvas("cFit1DCostSigmaFreeX2","cFit1DCostSigmaFreeX2",4,132,1024,768);
  hNJpsiCost -> Draw();
  histNJpsiCostSigmaFree_ANACX2 -> Draw("same");
  func1DCostSigmaFreeX2 -> Draw("same");*/
  //legNJpsiCost -> Draw("same");

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
  histNJpsiSigmaFixed_ANACLW -> Fit(func2D_W_HE_LW,"RSI0WL");

  TCanvas *c_fit2D_HE_LW = new TCanvas("c_fit2D_HE_LW","c_fit2D_HE_LW",4,132,1024,768);
  histNJpsiSigmaFixed_ANACLW -> Draw("COLZ");
  func2D_W_HE_LW -> Draw("same");*/

  /*double central_binx;
  double central_biny;

  TH2D *histo_diff_fit_LW = new TH2D("histo_diff_fit_LW","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      central_binx = histNJpsiSigmaFixed_ANACLW -> GetXaxis() -> GetBinCenter(i+1,j+1);
      central_biny = histNJpsiSigmaFixed_ANACLW -> GetYaxis() -> GetBinCenter(i+1,j+1);
      histo_diff_fit_LW -> SetBinContent(i+1,j+1,histNJpsiSigmaFixed_ANACLW -> GetBinContent(i+1,j+1)/func2D_W_HE_LW -> Eval(central_binx,central_biny));
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
