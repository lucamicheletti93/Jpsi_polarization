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

void check_MC(int pt_min, int pt_max){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();

  ostringstream convert_pt_min;
  convert_pt_min << pt_min;
  string str_pt_min =  convert_pt_min.str();

  ostringstream convert_pt_max;
  convert_pt_max << pt_max;
  string str_pt_max =  convert_pt_max.str();

  string dataset = str_pt_min + "pt" + str_pt_max;
  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";

  //============================================================================
  printf("---> Defining the pT ranges ... \n");
  //============================================================================
  string fileNameInput = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/HistosFromOfficialTree_Jpsi_PbPb_Nopol_TH3rec.root";
  printf("Opening %s ... \n",fileNameInput.c_str());
  TFile *fileInput = new TFile(fileNameInput.c_str(),"READ");

  TH2D *histCostPhiGen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",pt_min,pt_min+1));
  TH3D *histCostPhiMassRec = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",pt_min,pt_min+1));
  cout << pt_min << " - " << pt_min+1 << endl;
  //cout << histCostPhiGen -> GetBinContent(20,20) << " " << histCostPhiGen -> GetBinError(20,20) << endl;

  for(int i = pt_min+1;i < pt_max;i++){
    cout << i << " - " << i+1 << endl;
    TH2D *tmpTH2HistoGen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",i,i+1));
    //cout << tmpTH2HistoGen -> GetBinContent(20,20) << " " << tmpTH2HistoGen -> GetBinError(20,20) << endl;
    histCostPhiGen -> Add(tmpTH2HistoGen);
    TH3D *tmpTH3Histo = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",i,i+1));
    histCostPhiMassRec -> Add(tmpTH3Histo);
  }
  TH1D *histMass = (TH1D*) histCostPhiMassRec -> ProjectionZ();

  //============================================================================
  printf("---> Fitting the mass spectra and getting the J/psi sigma ... \n");
  //============================================================================
  char histMassName[100];

  double parTails[4] = {1.06,3.23,2.55,1.56};
  double parSignal[3] = {20.,3.096,7.0e-02};

  TF1 *funcJpsiCB2 = new TF1("funcJpsiCB2",Func_Jpsi_CB2,2.2,3.5,7);
  funcJpsiCB2 -> SetParameter(0,parSignal[0]);
  funcJpsiCB2 -> SetParameter(1,parSignal[1]);
  funcJpsiCB2 -> SetParameter(2,parSignal[2]);
  funcJpsiCB2 -> SetParLimits(2,0.,0.15);
  funcJpsiCB2 -> FixParameter(3,parTails[0]);
  funcJpsiCB2 -> FixParameter(4,parTails[1]);
  funcJpsiCB2 -> FixParameter(5,parTails[2]);
  funcJpsiCB2 -> FixParameter(6,parTails[3]);

  histMass -> Fit(funcJpsiCB2,"RL");
  histMass -> Draw();

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
