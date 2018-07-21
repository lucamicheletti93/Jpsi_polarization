#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include <TROOT.h>
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
void LoadStyle();

void accxeff(int pt_min, int pt_max){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  //LoadStyle();
  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();

  ostringstream convert_pt_min;
  convert_pt_min << pt_min;
  string str_pt_min =  convert_pt_min.str();

  ostringstream convert_pt_max;
  convert_pt_max << pt_max;
  string str_pt_max =  convert_pt_max.str();

  string dataset = str_pt_min + "pt" + str_pt_max;
  string nameOption = "_test";

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + nameOption + ".root";
  string fileInputName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/HistosFromOfficialTree_Jpsi_PbPb_Nopol_TH3rec.root";
  string fileInputName1DBinned = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_" + dataset + nameOption + "/" + dataset + ".root";
  string fileInputName2D = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";
  string fileSigmaName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/sigma_Jpsi_" + dataset + nameOption + ".root";
  string fileOutputName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + nameOption + ".root";

  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileInput = new TFile(fileInputName.c_str(),"READ");

  double fitRangeMinCost = -0.7;
  double fitRangeMaxCost = 0.7;

  //============================================================================
  printf("---> Reading the binning file ... \n");
  //============================================================================
  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  printf("Opening %s ... \n",fileInputName.c_str());

  //============================================================================
  printf("---> Defining the pT ranges ... \n");
  //============================================================================
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
  TH2D *histCostPhiRec = (TH2D*) histCostPhiMassRec -> Project3D("yx");
  //cout << histCostPhiGen -> GetBinContent(20,20) << " " << histCostPhiGen -> GetBinError(20,20) << endl;

  //============================================================================
  printf("---> Creating 1D AccxEff along Cost and Phi ... \n");
  //============================================================================
  TH1D *histCostGen = new TH1D("histCostGen","histCostGen",100,-1,1);
  TH1D *histCostRec = new TH1D("histCostRec","histCostRec",100,-1,1);

  for(int i = 0;i < 100;i++){
      histCostGen -> SetBinContent(i+1,histCostPhiGen -> Integral(i+1,i+2,1,49));
      histCostRec -> SetBinContent(i+1,histCostPhiRec -> Integral(i+1,i+2,1,49));
  }

  TH1D *histPhiGen = new TH1D("histPhiGen","histPhiGen",50,0,PI);
  TH1D *histPhiRec = new TH1D("histPhiRec","histPhiRec",50,0,PI);

  for(int i = 0;i < 50;i++){
      histPhiGen -> SetBinContent(i+1,histCostPhiGen -> Integral(i+1,i+2,1,99));
      histPhiRec -> SetBinContent(i+1,histCostPhiRec -> Integral(i+1,i+2,1,99));
  }

  //============================================================================
  printf("---> Creating the TH2D with the binning tuned on data ... \n");
  //============================================================================
  TH1D *histCostGenRebin = new TH1D("histCostGenRebin","histCostGenRebin",NCostBins,&CostValues[0]);
  TH1D *histCostRecRebin = new TH1D("histCostRecRebin","histCostRecRebin",NCostBins,&CostValues[0]);
  TH2D *histCostPhiGenRebin = new TH2D("histCostPhiGenRebin","histCostPhiGenRebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *histCostPhiRecRebin = new TH2D("histCostPhiRecRebin","histCostPhiRecRebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  int rec = 0;
  int gen = 0;

  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      histCostPhiGenRebin -> SetBinContent(i+1,j+1,histCostPhiGen -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
      histCostPhiRecRebin -> SetBinContent(i+1,j+1,histCostPhiRec -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
      gen += histCostPhiGen -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      rec += histCostPhiRec -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
    }
    //histCostGenRebin -> SetBinContent(i+1,histCostPhiGen -> Integral(i+1,i+2,PhiBinsMin[2],PhiBinsMax[NPhiBins-2]));
    //histCostRecRebin -> SetBinContent(i+1,histCostPhiRec -> Integral(i+1,i+2,PhiBinsMin[2],PhiBinsMax[NPhiBins-2]));
    histCostGenRebin -> SetBinContent(i+1,gen);
    histCostRecRebin -> SetBinContent(i+1,rec);
    gen = 0;
    rec = 0;
  }

  //============================================================================
  printf("---> Fitting the mass spectra and getting the J/psi sigma ... \n");
  //============================================================================
  TH1D *histMass;
  char histMassName[100];

  double parTails[4] = {1.06,3.23,2.55,1.56};
  double parSignal[3] = {20.,3.096,7.0e-02};

  string titleHistSigmaJpsi = str_pt_min + " < #it{p}_{T} < " + str_pt_max + " GeV/#it{c}";
  TH1D *histSigmaJpsiCost = new TH1D("histSigmaJpsiCost","",NCostBins,&CostValues[0]);
  TH2D *histSigmaJpsi = new TH2D("histSigmaJpsi","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histSigmaJpsi -> SetTitle(titleHistSigmaJpsi.c_str());
  gStyle -> SetPaintTextFormat("2.0f");

  TF1 *funcJpsiCB2 = new TF1("funcJpsiCB2",Func_Jpsi_CB2,2.2,3.5,7);
  funcJpsiCB2 -> SetParameter(0,parSignal[0]);
  funcJpsiCB2 -> SetParameter(1,parSignal[1]);
  funcJpsiCB2 -> SetParameter(2,parSignal[2]);
  funcJpsiCB2 -> SetParLimits(2,0.,0.15);
  funcJpsiCB2 -> FixParameter(3,parTails[0]);
  funcJpsiCB2 -> FixParameter(4,parTails[1]);
  funcJpsiCB2 -> FixParameter(5,parTails[2]);
  funcJpsiCB2 -> FixParameter(6,parTails[3]);

  char title[100];
  double sigmaJpsi, errSigmaJpsi;

  //============================================================================
  printf("-> Fit along (Cost,Phi) ... \n");
  //============================================================================
  for(int i = 1;i < NCostBins-1;i++){
    TH1D *histMassIntegrated = new TH1D("histMassIntegrated","",120,2,5);
    for(int j = 1;j < NPhiBins-1;j++){
      sprintf(histMassName,"HE_2pt6_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMass = (TH1D*) histCostPhiMassRec -> ProjectionZ(histMassName,CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMassIntegrated -> Add(histMass);
      sprintf(histMassName,"sigma_plots/HE_2pt6_%icost%i_%iphi%i.png",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMass -> Fit(funcJpsiCB2,"RL0");

      sigmaJpsi = funcJpsiCB2 -> GetParameter(2);
      errSigmaJpsi = funcJpsiCB2 -> GetParError(2);

      histSigmaJpsi -> SetBinContent(i+1,j+1,sigmaJpsi*1000);
      histSigmaJpsi -> SetBinError(i+1,j+1,errSigmaJpsi*1000);

      //TCanvas *c_histSigmaJpsi = new TCanvas("c_histSigmaJpsi","c_histSigmaJpsi",20,20,600,600);
      //sprintf(title,"#it{#sigma}_{J/#psi} = (%3.2f +- %3.2f) MeV/#it{c}^{2}",sigmaJpsi*1000,errSigmaJpsi*1000);
      //TLatex *lat = new TLatex(0.55,0.70,title);
      //lat -> SetTextSize(0.04);
      //lat -> SetNDC();
      //lat -> SetTextFont(42);
      //lat -> Draw("same");

      //c_histSigmaJpsi -> SaveAs(histMassName);
      //delete c_histSigmaJpsi;
    }
  }

  //============================================================================
  printf("-> Fit along Phi ... \n");
  //============================================================================
  for(int i = 1;i < NCostBins-1;i++){
    TH1D *histMassIntegrated = new TH1D("histMassIntegrated","",120,2,5);
    for(int j = 1;j < NPhiBins-1;j++){
      sprintf(histMassName,"HE_2pt6_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMass = (TH1D*) histCostPhiMassRec -> ProjectionZ(histMassName,CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMassIntegrated -> Add(histMass);
    }
    histMassIntegrated -> Fit(funcJpsiCB2,"RL0");

    sigmaJpsi = funcJpsiCB2 -> GetParameter(2);
    errSigmaJpsi = funcJpsiCB2 -> GetParError(2);

    histSigmaJpsiCost -> SetBinContent(i+1,sigmaJpsi*1000);
    histSigmaJpsiCost -> SetBinError(i+1,errSigmaJpsi*1000);
  }

  //============================================================================
  printf("-> Fit along Phi ... \n");
  //============================================================================
  TH1D *histSigmaJpsiPhi = new TH1D("histSigmaJpsiPhi","",NPhiBins,&PhiValues[0]);
  for(int i = 1;i < NPhiBins-1;i++){
    TH1D *histMassIntegrated = new TH1D("histMassIntegrated","",120,2,5);
    for(int j = 1;j < NCostBins-1;j++){
      sprintf(histMassName,"HE_2pt6_%icost%i_%iphi%i",CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]);
      histMass = (TH1D*) histCostPhiMassRec -> ProjectionZ(histMassName,CostBinsMin[j],CostBinsMax[j],PhiBinsMin[i],PhiBinsMax[i]);
      histMassIntegrated -> Add(histMass);
    }
    histMassIntegrated -> Fit(funcJpsiCB2,"RL0");

    sigmaJpsi = funcJpsiCB2 -> GetParameter(2);
    errSigmaJpsi = funcJpsiCB2 -> GetParError(2);

    histSigmaJpsiPhi -> SetBinContent(i+1,sigmaJpsi*1000);
    histSigmaJpsiPhi -> SetBinError(i+1,errSigmaJpsi*1000);
  }

  //============================================================================
  printf("---> Drawing J/psi sigma ... \n");
  //============================================================================
  TFile *fileInput1DBinned = new TFile(fileInputName1DBinned.c_str(),"READ");
  TFile *fileInput2D = new TFile(fileInputName2D.c_str(),"READ");

  TCanvas *cSigmaJpsi = new TCanvas("cSigmaJpsi","cSigmaJpsi",20,20,600,600);
  cSigmaJpsi -> Divide(2,2);

  cSigmaJpsi -> cd(1);
  TH2D *hSigmaJpsi = new TH2D("hSigmaJpsi","#sigma_{J/#psi} vs (cos(#it{#theta}),#it{#varphi})",100,-1,1,50,0,PI);
  hSigmaJpsi -> GetXaxis() -> SetTitle("cos(#it{#theta})");
  hSigmaJpsi -> GetYaxis() -> SetTitle("#it{#varphi}");
  hSigmaJpsi -> Draw();
  histSigmaJpsi -> Draw("COLZtexterrsame");

  cSigmaJpsi -> cd(2);
  TH2D *hSigmaJpsiPhi = new TH2D("hSigmaJpsiPhi","#sigma_{J/#psi} vs #it{#varphi}",50,0,PI,100,20,100);
  hSigmaJpsiPhi -> GetXaxis() -> SetTitle("#it{#varphi}");
  hSigmaJpsiPhi -> GetYaxis() -> SetTitle("#sigma_{J/#psi}");
  hSigmaJpsiPhi -> GetYaxis() -> SetTitleOffset(0.95);
  hSigmaJpsiPhi -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisSigmaJpsiPhi = new TGaxis(0,25,0,95,25,95,510,"");
  axisSigmaJpsiPhi -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisSigmaJpsiPhi -> SetLabelSize(15);

  TH1D *histSigmaJpsiPhiBinned = (TH1D*) fileInput1DBinned -> Get("histSigmaJpsiPhi");
  histSigmaJpsiPhiBinned -> SetMarkerColor(kRed);
  histSigmaJpsiPhiBinned -> SetLineColor(kRed);
  histSigmaJpsiPhiBinned -> Scale(1000);

  TPad *pad1 = new TPad("pad1","pad1",0.,0.3,1.,1.);
  pad1 -> SetBottomMargin(0);
  pad1 -> Draw();
  pad1 -> cd();
  hSigmaJpsiPhi -> Draw();
  axisSigmaJpsiPhi -> Draw("same");
  histSigmaJpsiPhi -> Draw("same");
  histSigmaJpsiPhiBinned -> Draw("same");

  cSigmaJpsi -> cd(2);
  TPad *pad2 = new TPad("pad2","pad2",0.,0.05,1.,0.3);
  pad2 -> SetTopMargin(0);
  pad2 -> SetBottomMargin(0.2);
  pad2 -> Draw();
  pad2 -> cd();
  TH2D *hRatioSigmaDataMCPhi = new TH2D("hRatioSigmaDataMCPhi","",50,0.,PI,100,0.5,1.6);
  hRatioSigmaDataMCPhi -> GetXaxis() -> SetTitle("#it{#varphi}");
  hRatioSigmaDataMCPhi -> GetXaxis() -> SetTitleSize(0.1);
  hRatioSigmaDataMCPhi -> GetXaxis() -> SetLabelSize(0.1);
  hRatioSigmaDataMCPhi -> GetYaxis() -> SetTitle("#sigma_{Data}/#sigma_{MC}");
  hRatioSigmaDataMCPhi -> GetYaxis() -> SetTitleOffset(0.95);
  hRatioSigmaDataMCPhi -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisRatioSigmaDataMCPhi = new TGaxis(0,0.6,0,1.5,0.6,1.5,510,"");
  axisRatioSigmaDataMCPhi -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisRatioSigmaDataMCPhi -> SetLabelSize(10);

  TLine *lUnityPhi = new TLine(0,1,PI,1);
  lUnityPhi -> SetLineStyle(2);
  lUnityPhi -> SetLineWidth(2);

  TH1D *histRatioSigmaDataMCPhi = new TH1D("histRatioSigmaDataMCPhi","",NPhiBins,&PhiValues[0]);
  histRatioSigmaDataMCPhi -> Divide(histSigmaJpsiPhiBinned,histSigmaJpsiPhi,1,1);
  hRatioSigmaDataMCPhi -> Draw();
  axisRatioSigmaDataMCPhi -> Draw("same");
  histRatioSigmaDataMCPhi -> Draw("same");
  lUnityPhi -> Draw("same");

  cSigmaJpsi -> cd(3);
  TH2D *hSigmaJpsiCost = new TH2D("hSigmaJpsiCost","#sigma_{J/#psi} vs cos(#it{#theta})",100,-1,1,100,20,120);
  hSigmaJpsiCost -> GetXaxis() -> SetTitle("cos(#it{#theta})");
  hSigmaJpsiCost -> GetYaxis() -> SetTitle("#sigma_{J/#psi}");
  hSigmaJpsiCost -> GetYaxis() -> SetTitleOffset(0.95);
  hSigmaJpsiCost -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisSigmaJpsiCost = new TGaxis(-1.,25,-1.,115,25,115,510,"");
  axisSigmaJpsiCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisSigmaJpsiCost -> SetLabelSize(15);

  TH1D *histSigmaJpsiCostBinned = (TH1D*) fileInput1DBinned -> Get("histSigmaJpsiCost");
  histSigmaJpsiCostBinned -> SetMarkerColor(kRed);
  histSigmaJpsiCostBinned -> SetLineColor(kRed);
  histSigmaJpsiCostBinned -> Scale(1000);

  TF1 *funcSigmaJpsiCost = new TF1("funcSigmaJpsiCost","pol2",fitRangeMinCost,fitRangeMaxCost);
  funcSigmaJpsiCost -> SetLineColor(kBlue);
  funcSigmaJpsiCost -> SetLineStyle(kDashed);
  histSigmaJpsiCost -> Fit(funcSigmaJpsiCost,"0");

  TF1 *funcSigmaJpsiCostBinned = new TF1("funcSigmaJpsiCostBinned","pol2",fitRangeMinCost,fitRangeMaxCost);
  funcSigmaJpsiCostBinned -> SetLineColor(kRed);
  funcSigmaJpsiCostBinned -> SetLineStyle(kDashed);
  histSigmaJpsiCostBinned -> Fit(funcSigmaJpsiCostBinned,"0");

  TPad *pad3 = new TPad("pad3","pad3",0.,0.3,1.,1.);
  pad3 -> SetBottomMargin(0);
  pad3 -> Draw();
  pad3 -> cd();
  hSigmaJpsiCost -> Draw();
  axisSigmaJpsiCost -> Draw("same");
  histSigmaJpsiCost -> Draw("same");
  histSigmaJpsiCostBinned -> Draw("same");
  funcSigmaJpsiCost -> Draw("same");
  funcSigmaJpsiCostBinned -> Draw("same");

  cSigmaJpsi -> cd(3);
  gStyle -> SetPaintTextFormat("2.1f");
  TPad *pad4 = new TPad("pad4","pad4",0.,0.05,1.,0.3);
  pad4 -> SetTopMargin(0);
  pad4 -> SetBottomMargin(0.2);
  pad4 -> Draw();
  pad4 -> cd();
  TH2D *hRatioSigmaDataMCCost = new TH2D("hRatioSigmaDataMCCost","",100,-1.,1.,100,0.5,1.6);
  hRatioSigmaDataMCCost -> GetXaxis() -> SetTitle("cos(#it{#theta})");
  hRatioSigmaDataMCCost -> GetXaxis() -> SetTitleSize(0.1);
  hRatioSigmaDataMCCost -> GetXaxis() -> SetLabelSize(0.1);
  hRatioSigmaDataMCCost -> GetYaxis() -> SetTitle("#sigma_{Data}/#sigma_{MC}");
  hRatioSigmaDataMCCost -> GetYaxis() -> SetTitleOffset(0.95);
  hRatioSigmaDataMCCost -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisRatioSigmaDataMCCost = new TGaxis(-1.,0.6,-1.,1.5,0.6,1.5,510,"");
  axisRatioSigmaDataMCCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisRatioSigmaDataMCCost -> SetLabelSize(10);

  TLine *lUnityCost = new TLine(-1.,1,1.,1);
  lUnityCost -> SetLineStyle(2);
  lUnityCost -> SetLineWidth(2);

  TH1D *histRatioSigmaDataMCCost = new TH1D("histRatioSigmaDataMCCost","",NCostBins,&CostValues[0]);
  histRatioSigmaDataMCCost -> Divide(histSigmaJpsiCostBinned,histSigmaJpsiCost,1,1);
  hRatioSigmaDataMCCost -> Draw();
  axisRatioSigmaDataMCCost -> Draw("same");
  histRatioSigmaDataMCCost -> Draw("same");
  lUnityCost -> Draw("same");

  cSigmaJpsi -> cd(4);
  TH2D *histSigmaJpsiCostPhiData = (TH2D*) fileInput2D -> Get("histSigmaJpsi");
  TH2D *histRatioSigmaDataMCCostPhi = new TH2D("histRatioSigmaDataMCCostPhi","#sigma_{J/#psi}^{DATA}/#sigma_{J/#psi}^{MC}",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histRatioSigmaDataMCCostPhi -> Divide(histSigmaJpsiCostPhiData,histSigmaJpsi,1,1);
  hSigmaJpsi -> Draw();
  histRatioSigmaDataMCCostPhi -> Draw("COLZtextsame");

  //============================================================================
  printf("---> Saving sigma Jpsi ... \n");
  //============================================================================
  printf("Saving J/psi width in %s ... \n",fileSigmaName.c_str());
  TFile *fileSigma = new TFile(fileSigmaName.c_str(),"RECREATE");
  histSigmaJpsi -> Write();
  fileSigma -> Close();

  //============================================================================
  printf("---> Computing Acc X Eff ... \n");
  //============================================================================
  histCostRec -> Sumw2();
  histCostGen -> Sumw2();
  histCostGenRebin -> Sumw2();
  histCostRecRebin -> Sumw2();
  histCostPhiRec -> Sumw2();
  histCostPhiGen -> Sumw2();
  histCostPhiRecRebin -> Sumw2();
  histCostPhiGenRebin -> Sumw2();

  TH1D *hist_accxeff_HE_Cost = new TH1D("hist_accxeff_HE_Cost","",100,-1,1);
  hist_accxeff_HE_Cost -> Divide(histCostRec,histCostGen,1,1,"B");

  TH1D *hist_accxeff_HE_Cost_rebin = new TH1D("hist_accxeff_HE_Cost_rebin","",NCostBins,&CostValues[0]);
  hist_accxeff_HE_Cost_rebin -> Divide(histCostRecRebin,histCostGenRebin,1,1,"B");
  hist_accxeff_HE_Cost_rebin -> SetLineColor(kRed);

  TH2D *hist_accxeff_HE = new TH2D("hist_accxeff_HE","",100,-1,1,50,0,PI);
  hist_accxeff_HE -> Divide(histCostPhiRec,histCostPhiGen,1,1,"B");

  TH2D *hist_accxeff_HE_rebin = new TH2D("hist_accxeff_HE_rebin","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  hist_accxeff_HE_rebin -> Divide(histCostPhiRecRebin,histCostPhiGenRebin,1,1,"B");

  TCanvas *cAccxEff = new TCanvas("cAccxEff","cAccxEff",800,800);

  TPad *pad5 = new TPad("pad5","pad5",0,0.3,1,1.0);
  pad5 -> SetBottomMargin(0);
  pad5 -> Draw();
  pad5 -> cd();
  TH2D *hAccxEff = new TH2D("hAccxEff","",100,-1,1,50,0,PI);
  hSigmaJpsi -> Draw();
  hist_accxeff_HE -> Draw("COLZsame");

  cAccxEff -> cd();
  TPad *pad6 = new TPad("pad6","pad6",0,0.05,1,0.3);
  pad6 -> SetTopMargin(0);
  pad6 -> SetBottomMargin(0.2);
  pad6 -> Draw();
  pad6 -> cd();
  TH2D *hAccxEffCost = new TH2D("hAccxEffCost","",100,-1,1,100,0,1);
  hAccxEffCost -> Draw();
  hist_accxeff_HE_Cost -> Draw("same");
  hist_accxeff_HE_Cost_rebin -> Draw("same");

  //============================================================================
  printf("---> Saving Acc X Eff ... \n");
  //============================================================================
  printf("Saving Acc x Eff in %s ... \n",fileOutputName.c_str());
  TFile *fileOutput = new TFile(fileOutputName.c_str(),"RECREATE");
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
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
  int font = 42;
  gROOT -> SetStyle("Plain");
  gStyle -> SetPalette(1);
  gStyle -> SetFrameBorderMode(0);
  gStyle -> SetFrameFillColor(0);
  gStyle -> SetCanvasBorderMode(0);
  gStyle -> SetPadBorderMode(0);
  gStyle -> SetPadColor(10);
  gStyle -> SetCanvasColor(10);
  gStyle -> SetTitleFillColor(10);
  gStyle -> SetTitleBorderSize(1);
  gStyle -> SetStatColor(10);
  gStyle -> SetStatBorderSize(1);
  gStyle -> SetLegendBorderSize(1);
  gStyle -> SetDrawBorder(0);
  gStyle -> SetTextFont(font);
  gStyle -> SetStatFont(font);
  gStyle -> SetStatFontSize(0.05);
  gStyle -> SetStatX(0.97);
  gStyle -> SetStatY(0.98);
  gStyle -> SetStatH(0.03);
  gStyle -> SetStatW(0.3);
  gStyle -> SetTickLength(0.02,"y");
  gStyle -> SetEndErrorSize(3);
  gStyle -> SetLabelSize(0.04,"xyz");
  gStyle -> SetLabelFont(font,"xyz");
  gStyle -> SetLabelOffset(0.01,"xyz");
  gStyle -> SetTitleFont(font,"xyz");
  gStyle -> SetTitleOffset(0.9,"x");
  gStyle -> SetTitleOffset(0.7,"y");
  gStyle -> SetTitleSize(0.05,"xyz");
  gStyle -> SetMarkerSize(1.3);
  gStyle -> SetPalette(1,0);
  gROOT -> ForceStyle();
  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);
  gStyle -> SetEndErrorSize(0);
}
