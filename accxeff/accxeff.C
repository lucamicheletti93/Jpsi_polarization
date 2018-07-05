#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include <TMinuit.h>
#include <TDirectory.h>
#include <TSystem.h>
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

#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

Double_t Func_Jpsi_CB2(Double_t *, Double_t *);

void accxeff(int pt_min, int pt_max){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();


  TDirectory *dir = new TDirectory("directory","directory");
  return;

  ostringstream convert_pt_min;
  convert_pt_min << pt_min;
  string str_pt_min =  convert_pt_min.str();

  ostringstream convert_pt_max;
  convert_pt_max << pt_max;
  string str_pt_max =  convert_pt_max.str();

  string dataset = str_pt_min + "pt" + str_pt_max;
  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";

  //============================================================================
  printf("---> Reading the binning file ... \n");
  //============================================================================
  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  Binning *binning = (Binning*) fileBinning -> Get("Binning");

  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  string fileNameInput = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/HistosFromOfficialTree_Jpsi_PbPb_Nopol_TH3rec.root";
  printf("Opening %s ... \n",fileNameInput.c_str());
  TFile *fileInput = new TFile(fileNameInput.c_str(),"READ");

  //============================================================================
  printf("Defining the pT ranges ... \n");
  //============================================================================
  TH2D *hist_CostPhiHE_gen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",pt_min,pt_min+1));
  TH3D *hist_CostPhiMassHE_rec = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",pt_min,pt_min+1));
  cout << pt_min << " - " << pt_min+1 << endl;
  //cout << hist_CostPhiHE_gen -> GetBinContent(20,20) << " " << hist_CostPhiHE_gen -> GetBinError(20,20) << endl;

  for(int i = pt_min+1;i < pt_max;i++){
    cout << i << " - " << i+1 << endl;
    TH2D *tmpTH2HistoGen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",i,i+1));
    //cout << tmpTH2HistoGen -> GetBinContent(20,20) << " " << tmpTH2HistoGen -> GetBinError(20,20) << endl;
    hist_CostPhiHE_gen -> Add(tmpTH2HistoGen);
    TH3D *tmpTH3Histo = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",i,i+1));
    hist_CostPhiMassHE_rec -> Add(tmpTH3Histo);
  }
  TH2D *hist_CostPhiHE_rec = (TH2D*) hist_CostPhiMassHE_rec -> Project3D("yx");
  //cout << hist_CostPhiHE_gen -> GetBinContent(20,20) << " " << hist_CostPhiHE_gen -> GetBinError(20,20) << endl;

  //============================================================================
  printf("Creating the TH2D with the binning tuned on data ... \n");
  //============================================================================
  TH2D *hist_CostPhiHE_gen_rebin = new TH2D("hist_CostPhiHE_gen_rebin","hist_CostPhiHE_gen_rebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *hist_CostPhiHE_rec_rebin = new TH2D("hist_CostPhiHE_rec_rebin","hist_CostPhiHE_rec_rebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      hist_CostPhiHE_gen_rebin -> SetBinContent(i+1,j+1,hist_CostPhiHE_gen -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
      hist_CostPhiHE_rec_rebin -> SetBinContent(i+1,j+1,hist_CostPhiHE_rec -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
    }
  }

  //============================================================================
  printf("Fitting the mass spectra ... \n");
  //============================================================================
  TH1D *histMass;
  char histMassName[100];

  double parTails[4] = {1.06,3.23,2.55,1.56};
  double parSignal[3] = {20.,3.096,7.0e-02};

  TH2D *histoSigmaJpsi = new TH2D("histoSigmaJpsi","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  TF1 *funcJpsiCB2 = new TF1("funcJpsiCB2",Func_Jpsi_CB2,2.2,3.5,7);
  funcJpsiCB2 -> SetParameter(0,parSignal[0]);
  funcJpsiCB2 -> SetParameter(1,parSignal[1]);
  funcJpsiCB2 -> SetParameter(2,parSignal[2]);
  funcJpsiCB2 -> SetParLimits(2,0.,0.15);
  funcJpsiCB2 -> FixParameter(3,parTails[0]);
  funcJpsiCB2 -> FixParameter(4,parTails[1]);
  funcJpsiCB2 -> FixParameter(5,parTails[2]);
  funcJpsiCB2 -> FixParameter(6,parTails[3]);

  for(int i = 1;i < NCostBins-1;i++){
    for(int j = 1;j < NPhiBins-1;j++){
      sprintf(histMassName,"HE_2pt6_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMass = (TH1D*) hist_CostPhiMassHE_rec -> ProjectionZ(histMassName,CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      sprintf(histMassName,"sigma_plots/HE_2pt6_%icost%i_%iphi%i.png",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      TCanvas *c_histoSigmaJpsi = new TCanvas("c_histoSigmaJpsi","c_histoSigmaJpsi",20,20,600,600);
      histMass -> Fit(funcJpsiCB2,"RL");
      histoSigmaJpsi -> SetBinContent(i+1,j+1,funcJpsiCB2 -> GetParameter(2));
      histMass -> Draw();
      c_histoSigmaJpsi -> SaveAs(histMassName);
    }
  }
  histoSigmaJpsi -> Draw("COLZtext");


  return;

  //============================================================================
  printf("Computing Acc X Eff ... \n");
  //============================================================================
  hist_CostPhiHE_rec -> Sumw2();
  hist_CostPhiHE_gen -> Sumw2();
  hist_CostPhiHE_rec_rebin -> Sumw2();
  hist_CostPhiHE_gen_rebin -> Sumw2();

  TH2D *hist_accxeff_HE = new TH2D("hist_accxeff_HE","",100,-1,1,50,0,PI);
  hist_accxeff_HE -> Divide(hist_CostPhiHE_rec,hist_CostPhiHE_gen,1,1,"B");

  TH2D *hist_accxeff_HE_rebin = new TH2D("hist_accxeff_HE_rebin","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  hist_accxeff_HE_rebin -> Divide(hist_CostPhiHE_rec_rebin,hist_CostPhiHE_gen_rebin,1,1,"B");

  //============================================================================
  printf("Saving Acc X Eff ... \n");
  //============================================================================
  return;
  string filePathOutput = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_";
  string fileNameOutput = filePathOutput + dataset + ".root";
  printf("Opening %s ... \n",fileNameOutput.c_str());
  TFile *fileOutput = new TFile(fileNameOutput.c_str(),"RECREATE");
  hist_accxeff_HE -> Write();
  hist_accxeff_HE_rebin -> Write();
  fileOutput -> Close();

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Double_t Func_Jpsi_CB2(Double_t *x, Double_t *par){
  Double_t t = (x[0] - par[1])/par[2];
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
