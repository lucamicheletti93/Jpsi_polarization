#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TObjectTable.h>
#include <TGraphErrors.h>
#endif

void SumRunPbPbJPsi_AllRuns(const char *period="LHC15o", Bool_t PS=kTRUE){

 TH1D *hNEv = new TH1D("hNEv","hNEv",6,0.,6.);
 TString namelabel1[6]={"TotEv","CINT7","CMUL7","CMLL7","CMSL7","CMSH7"};
 for(int k=0;k<6;k++) hNEv->GetXaxis()->SetBinLabel(k+1,namelabel1[k]);
 hNEv->Sumw2();
 
 //
 // dimuon histos
 //
TH1D *hMassOS_notrig = new TH1D("hMassOS_notrig","hMassOS_notrig",600,0.,15.);
TH1D *hMassOS = new TH1D("hMassOS","hMassOS",600,0.,15.);
TH1D *hMassOS_CINT7 = new TH1D("hMassOS_CINT7","hMassOS_CINT7",600,0.,15.);
TH1D *hMassOS_CMUL7 = new TH1D("hMassOS_CMUL7","hMassOS_CMUL7",600,0.,15.);
TH1D *hMassOS_1m = new TH1D("hMassOS_1m","hMassOS_1m",600,0.,15.);
TH1D *hMassOS_CINT7_1m = new TH1D("hMassOS_CINT7_1m","hMassOS_CINT7_1m",600,0.,15.);
TH1D *hMassOS_CMUL7_1m = new TH1D("hMassOS_CMUL7_1m","hMassOS_CMUL7_1m",600,0.,15.);
TH1D *hMassOS_CMSL7_1m = new TH1D("hMassOS_CMSL7_1m","hMassOS_CMSL7_1m",600,0.,15.);
TH1D *hMassOS_2m = new TH1D("hMassOS_2m","hMassOS_2m",600,0.,15.);
TH1D *hMassOS_CINT7_2m = new TH1D("hMassOS_CINT7_2m","hMassOS_CINT7_2m",600,0.,15.);
TH1D *hMassOS_CMUL7_2m = new TH1D("hMassOS_CMUL7_2m","hMassOS_CMUL7_2m",600,0.,15.);
TH1D *hMassOS_CMSL7_2m = new TH1D("hMassOS_CMSL7_2m","hMassOS_CMSL7_2m",600,0.,15.);

TH1D *hMassPP_CMLL7_2m = new TH1D("hMassPP_CMLL7_2m","hMassPP_CMLL7_2m",600,0.,15.);
TH1D *hMassMM_CMLL7_2m = new TH1D("hMassMM_CMLL7_2m","hMassMM_CMLL7_2m",600,0.,15.);

TH1D *hPt_CMUL7_2m = new TH1D("hPt_CMUL7_2m","hPt_CMUL7_2m",150,0.,15.);
TH1D *hY_CMUL7_2m = new TH1D("hY_CMUL7_2m","hY_CMUL7_2m",100,-5,0);

TH1D *hMu0MatchTrigger = new TH1D("hMu0MatchTrigger","hMu0MatchTrigger",4,0.,4.);
TH1D *hNContrib = new TH1D("hNContrib","hNContrib",2000,0.,2000.);
TH1D *hZVertex = new TH1D("hZVertex","hZVertex",2000,-100.,100.);

hMassOS_notrig->Sumw2();
hMassOS->Sumw2(); 
hMassOS_CINT7->Sumw2(); 
hMassOS_CMUL7->Sumw2();  
hMassOS_1m->Sumw2(); 
hMassOS_CINT7_1m->Sumw2();
hMassOS_CMUL7_1m->Sumw2();
hMassOS_CMSL7_1m->Sumw2();
hMassOS_2m->Sumw2(); 
hMassOS_CINT7_2m->Sumw2();
hMassOS_CMUL7_2m->Sumw2();
hMassOS_CMSL7_2m->Sumw2();

hMassPP_CMLL7_2m->Sumw2();
hMassMM_CMLL7_2m->Sumw2();
hPt_CMUL7_2m->Sumw2();
hY_CMUL7_2m->Sumw2();

hMu0MatchTrigger->Sumw2();
hNContrib->Sumw2();  
hZVertex->Sumw2(); 
 
 // single muon histos
 //
TH1D *hNMu = new TH1D("hNMu","hNMu",50,0.,50.);
//TH1D *hNMu_Match = new TH1D("hNMu_Match","hNMu_Match",50,0.,50.);
TH1D *hNMuons_1m = new TH1D("hNMuons_1m","hNMuons_1m",50,0.,50.);
TH1D *hPtMu_Match = new TH1D("hPtMu_Match","hPtMu_Match",600,0.,15.);
hNMu->Sumw2();
//hNMu_Match->Sumw2();
hNMuons_1m->Sumw2();
hPtMu_Match->Sumw2();
 //
 // stat histos
 //
TH1D *hNDimuons = new TH1D("hNDimuons","hNDimuons",50,0.,50.); 
TH1D *hNDimuons_2m = new TH1D("hNDimuons_2m","hNDimuons_2m",50,0.,50.); 
hNDimuons->Sumw2();
hNDimuons_2m->Sumw2();
 //
 // centrality
 //
 TH1F *hPercentileV0M = new TH1F("hPercentileV0M","hPercentileV0M",201,-1,200);  //V0M
 TH1F *hPercentileV0M_CMUL7 = new TH1F("hPercentileV0M_CMUL7","hPercentileV0M_CMUL7",201,-1,200);
 TH1F *hPercentileV0M_CINT7 = new TH1F("hPercentileV0M_CINT7","hPercentileV0M_CINT7",201,-1,200);

 hPercentileV0M->Sumw2();
 hPercentileV0M_CMUL7->Sumw2();
 hPercentileV0M_CINT7->Sumw2();

TH1D *hAllpt_1M = new TH1D("hAllpt_1M","hAllpt_1M",500,0.,50.);      
TH1D *hLpt_1M = new TH1D("hLpt_1M","hLpt_1M",500,0.,50.);

TH1D *hAllpt_25y3 = new TH1D("hAllpt_25y3","hAllpt_25y3",500,0.,50.);   
TH1D *hLpt_25y3 = new TH1D("hLpt_25y3","hLpt_25y3",500,0.,50.);
TH1D *hAllpt_3y35 = new TH1D("hAllpt_3y35","hAllpt_3y35",500,0.,50.); 
TH1D *hLpt_3y35 = new TH1D("hLpt_3y35","hLpt_3y35",500,0.,50.);
TH1D *hAllpt_35y4 = new TH1D("hAllpt_35y4","hAllpt_35y4",500,0.,50.);
TH1D *hLpt_35y4 = new TH1D("hLpt_35y4","hLpt_35y4",500,0.,50.); 

TH1D *hAllpt_25y3_1M = new TH1D("hAllpt_25y3_1M","hAllpt_25y3_1M",500,0.,50.);   
TH1D *hLpt_25y3_1M = new TH1D("hLpt_25y3_1M","hLpt_25y3_1M",500,0.,50.);
TH1D *hAllpt_3y35_1M = new TH1D("hAllpt_3y35_1M","hAllpt_3y35_1M",500,0.,50.); 
TH1D *hLpt_3y35_1M = new TH1D("hLpt_3y35_1M","hLpt_3y35_1M",500,0.,50.);
TH1D *hAllpt_35y4_1M = new TH1D("hAllpt_35y4_1M","hAllpt_35y4_1M",500,0.,50.);
TH1D *hLpt_35y4_1M = new TH1D("hLpt_35y4_1M","hLpt_35y4_1M",500,0.,50.); 
  
TH1D *hAllpt_010_1M = new TH1D("hAllpt_010_1M","hAllpt_010_1M",500,0.,50.);
TH1D *hLpt_010_1M = new TH1D("hLpt_010_1M","hLpt_010_1M",500,0.,50.);
TH1D *hAllpt_1020_1M = new TH1D("hAllpt_1020_1M","hAllpt_1020_1M",500,0.,50.);
TH1D *hLpt_1020_1M = new TH1D("hLpt_1020_1M","hLpt_1020_1M",500,0.,50.);
TH1D *hAllpt_2030_1M = new TH1D("hAllpt_2030_1M","hAllpt_2030_1M",500,0.,50.);
TH1D *hLpt_2030_1M = new TH1D("hLpt_2030_1M","hLpt_2030_1M",500,0.,50.);
TH1D *hAllpt_3040_1M = new TH1D("hAllpt_3040_1M","hAllpt_3040_1M",500,0.,50.);
TH1D *hLpt_3040_1M = new TH1D("hLpt_3040_1M","hLpt_3040_1M",500,0.,50.);
TH1D *hAllpt_4050_1M = new TH1D("hAllpt_4050_1M","hAllpt_4050_1M",500,0.,50.);
TH1D *hLpt_4050_1M = new TH1D("hLpt_4050_1M","hLpt_4050_1M",500,0.,50.);
TH1D *hAllpt_5060_1M = new TH1D("hAllpt_5060_1M","hAllpt_5060_1M",500,0.,50.);
TH1D *hLpt_5060_1M = new TH1D("hLpt_5060_1M","hLpt_5060_1M",500,0.,50.);
TH1D *hAllpt_6070_1M = new TH1D("hAllpt_6070_1M","hAllpt_6070_1M",500,0.,50.);
TH1D *hLpt_6070_1M = new TH1D("hLpt_6070_1M","hLpt_6070_1M",500,0.,50.);
TH1D *hAllpt_7080_1M = new TH1D("hAllpt_7080_1M","hAllpt_7080_1M",500,0.,50.);
TH1D *hLpt_7080_1M = new TH1D("hLpt_7080_1M","hLpt_7080_1M",500,0.,50.);
TH1D *hAllpt_8090_1M = new TH1D("hAllpt_8090_1M","hAllpt_8090_1M",500,0.,50.);
TH1D *hLpt_8090_1M = new TH1D("hLpt_8090_1M","hLpt_8090_1M",500,0.,50.);


 hAllpt_1M->Sumw2();	 	
 hLpt_1M->Sumw2();

 hAllpt_25y3->Sumw2();
 hLpt_25y3->Sumw2(); 
 hAllpt_3y35->Sumw2(); 
 hLpt_3y35->Sumw2();
 hAllpt_35y4->Sumw2(); 
 hLpt_35y4->Sumw2();

 	 	
 hAllpt_25y3_1M->Sumw2();
 hLpt_25y3_1M->Sumw2();
 hAllpt_3y35_1M->Sumw2();
 hLpt_3y35_1M->Sumw2(); 	 	
 hAllpt_35y4_1M->Sumw2();	 	
 hLpt_35y4_1M->Sumw2(); 
 	 	
 hAllpt_010_1M->Sumw2();
 hLpt_010_1M->Sumw2();
 hAllpt_1020_1M->Sumw2();
 hLpt_1020_1M->Sumw2();
 hAllpt_2030_1M->Sumw2();
 hLpt_2030_1M->Sumw2();
 hAllpt_3040_1M->Sumw2();
 hLpt_3040_1M->Sumw2();
 hAllpt_4050_1M->Sumw2();
 hLpt_4050_1M->Sumw2();
 hAllpt_5060_1M->Sumw2();
 hLpt_5060_1M->Sumw2();
 hAllpt_6070_1M->Sumw2();
 hLpt_6070_1M->Sumw2();
 hAllpt_7080_1M->Sumw2();
 hLpt_7080_1M->Sumw2();
 hAllpt_8090_1M->Sumw2();
 hLpt_8090_1M->Sumw2();
//------------------------------------------------
// Loop on files
//------------------------------------------------

FILE *flist;
char runlist[200];

//sprintf(runlist,"dummy.dat"); 
sprintf(runlist,"RunList_%s.dat",period); 

if((flist = fopen(runlist,"r")) == NULL){
  printf("Cannot open file %s \n", runlist);
  return;
} else printf("Reading List %s\n",runlist);
Int_t Run=0;
Int_t NTotRun=1000;
Int_t ntotruns=0;
for(int irun=1; irun<=NTotRun; irun++){
  fscanf(flist,"%d",&Run);
  if (feof(flist)!=0) {
    break;
  }
  TFile *f;
  char file[300];
  printf("\nrun number = %d\n",irun);
    if(PS){
       sprintf(file,"PbPb_PS_%d.root",Run);
    } else {
      sprintf(file,"PbPb_noPS_%d.root",Run);
    }

  printf("file= %s\n",file);
  Long_t *dummy1 =0, *dummy2 =0, *dummy3 =0, *dummy4 =0;
  if (gSystem->GetPathInfo(file,dummy1,dummy2,dummy3, dummy4) == 0) {
    printf(" Opening %s\n",file);
    f= new TFile(file);
    ntotruns++;
  }  
 else continue;
       
 TList *list = (TList*) f->Get("chist0");
 printf("N CMUL Run %d = %f\n",Run,((TH1D*)list->FindObject("hNEv"))->GetBinContent(3));
 hNEv->Add(((TH1D*)list->FindObject("hNEv")));
 hMassOS_notrig->Add(((TH1D*)list->FindObject("hMassOS_notrig")));;
 hMassOS->Add(((TH1D*)list->FindObject("hMassOS")));; 
 hMassOS_CINT7->Add(((TH1D*)list->FindObject("hMassOS_CINT7")));; 
 hMassOS_CMUL7->Add(((TH1D*)list->FindObject("hMassOS_CMUL7")));;  
 hMassOS_1m->Add(((TH1D*)list->FindObject("hMassOS_1m")));; 
 hMassOS_CINT7_1m->Add(((TH1D*)list->FindObject("hMassOS_CINT7_1m")));;
 hMassOS_CMUL7_1m->Add(((TH1D*)list->FindObject("hMassOS_CMUL7_1m")));;
 hMassOS_CMSL7_1m->Add(((TH1D*)list->FindObject("hMassOS_CMSL7_1m")));;
 hMassOS_2m->Add(((TH1D*)list->FindObject("hMassOS_2m")));; 
 hMassOS_CINT7_2m->Add(((TH1D*)list->FindObject("hMassOS_CINT7_2m")));;
 hMassOS_CMUL7_2m->Add(((TH1D*)list->FindObject("hMassOS_CMUL7_2m")));;
 hMassOS_CMSL7_2m->Add(((TH1D*)list->FindObject("hMassOS_CMSL7_2m")));;

 hMassPP_CMLL7_2m->Add(((TH1D*)list->FindObject("hMassPP_CMLL7_2m")));;
 hMassMM_CMLL7_2m->Add(((TH1D*)list->FindObject("hMassMM_CMLL7_2m")));;

 hPt_CMUL7_2m->Add(((TH1D*)list->FindObject("hPt_CMUL7_2m")));;
 hY_CMUL7_2m->Add(((TH1D*)list->FindObject("hY_CMUL7_2m")));;

 hMu0MatchTrigger->Add(((TH1D*)list->FindObject("hMu0MatchTrigger")));;
 hNContrib->Add(((TH1D*)list->FindObject("hNContrib")));;  
 hZVertex->Add(((TH1D*)list->FindObject("hZVertex")));; 
 hNMu->Add(((TH1D*)list->FindObject("hNMu")));;
// hNMu_Match->Add(((TH1D*)list->FindObject("hNMu_Match")));;
 hNMuons_1m->Add(((TH1D*)list->FindObject("hNMuons_1m")));;
 hPtMu_Match->Add(((TH1D*)list->FindObject("hPtMu_Match")));;
 hNDimuons->Add(((TH1D*)list->FindObject("hNDimuons")));;
 hNDimuons_2m->Add(((TH1D*)list->FindObject("hNDimuons_2m")));;

 hPercentileV0M->Add(((TH1F*)list->FindObject("hPercentileV0M")));;
 hPercentileV0M_CMUL7->Add(((TH1F*)list->FindObject("hPercentileV0M_CMUL7")));;
 hPercentileV0M_CINT7->Add(((TH1F*)list->FindObject("hPercentileV0M_CINT7")));;

hAllpt_1M->Add(((TH1D *) list->FindObject("hAllpt_1M")));
hLpt_1M->Add(((TH1D *) list->FindObject("hLpt_1M")));

hAllpt_25y3->Add(((TH1D *) list->FindObject("hAllpt_25y3")));
hLpt_25y3->Add(((TH1D *) list->FindObject("hLpt_25y3")));
hAllpt_3y35->Add(((TH1D *) list->FindObject("hAllpt_3y35")));
hLpt_3y35->Add(((TH1D *) list->FindObject("hLpt_3y35")));
hAllpt_35y4->Add(((TH1D *) list->FindObject("hAllpt_35y4")));
hLpt_35y4->Add(((TH1D *) list->FindObject("hLpt_35y4")));

hAllpt_25y3_1M->Add(((TH1D *) list->FindObject("hAllpt_25y3_1M")));
hLpt_25y3_1M->Add(((TH1D *) list->FindObject("hLpt_25y3_1M")));
hAllpt_3y35_1M->Add(((TH1D *) list->FindObject("hAllpt_3y35_1M")));
hLpt_3y35_1M->Add(((TH1D *) list->FindObject("hLpt_3y35_1M")));
hAllpt_35y4_1M->Add(((TH1D *) list->FindObject("hAllpt_35y4_1M")));
hLpt_35y4_1M->Add(((TH1D *) list->FindObject("hLpt_35y4_1M")));

 hAllpt_010_1M->Add(((TH1D *) list->FindObject("hAllpt_010_1M")));
 hLpt_010_1M->Add(((TH1D *) list->FindObject("hLpt_010_1M")));
 hAllpt_1020_1M->Add(((TH1D *) list->FindObject("hAllpt_1020_1M")));
 hLpt_1020_1M->Add(((TH1D *) list->FindObject("hLpt_1020_1M")));
 hAllpt_2030_1M->Add(((TH1D *) list->FindObject("hAllpt_2030_1M")));
 hLpt_2030_1M->Add(((TH1D *) list->FindObject("hLpt_2030_1M")));
 hAllpt_3040_1M->Add(((TH1D *) list->FindObject("hAllpt_3040_1M")));
 hLpt_3040_1M->Add(((TH1D *) list->FindObject("hLpt_3040_1M")));
 hAllpt_4050_1M->Add(((TH1D *) list->FindObject("hAllpt_4050_1M")));
 hLpt_4050_1M->Add(((TH1D *) list->FindObject("hLpt_4050_1M")));
 hAllpt_5060_1M->Add(((TH1D *) list->FindObject("hAllpt_5060_1M")));
 hLpt_5060_1M->Add(((TH1D *) list->FindObject("hLpt_5060_1M")));
 hAllpt_6070_1M->Add(((TH1D *) list->FindObject("hAllpt_6070_1M")));
 hLpt_6070_1M->Add(((TH1D *) list->FindObject("hLpt_6070_1M")));
 hAllpt_7080_1M->Add(((TH1D *) list->FindObject("hAllpt_7080_1M")));
 hLpt_7080_1M->Add(((TH1D *) list->FindObject("hLpt_7080_1M")));
 hAllpt_8090_1M->Add(((TH1D *) list->FindObject("hAllpt_8090_1M")));
 hLpt_8090_1M->Add(((TH1D *) list->FindObject("hLpt_8090_1M")));
}

printf("Total number of analyzed runs = %d\n",ntotruns);     

//------------------------------------------------
// Output file
//------------------------------------------------
 char fnameout[300];
 if(PS) sprintf(fnameout,"Sum_%s_PS.root",period,period);
 else sprintf(fnameout,"Sum_%s_noPS.root",period,period);
 TFile *fout = new TFile(fnameout,"recreate");
 fout->cd();
 hNEv->Write();
 hMassOS_notrig->Write();
 hMassOS->Write(); 
 hMassOS_CINT7->Write(); 
 hMassOS_CMUL7->Write();  
 hMassOS_1m->Write(); 
 hMassOS_CINT7_1m->Write();
 hMassOS_CMUL7_1m->Write();
 hMassOS_CMSL7_1m->Write();
 hMassOS_2m->Write(); 
 hMassOS_CINT7_2m->Write();
 hMassOS_CMUL7_2m->Write();
 hMassOS_CMSL7_2m->Write();

 hMassPP_CMLL7_2m->Write();
 hMassMM_CMLL7_2m->Write();

 hPt_CMUL7_2m->Write();
 hY_CMUL7_2m->Write();

 hMu0MatchTrigger->Write();
 hNContrib->Write();  
 hZVertex->Write();

 hNMu->Write();
// hNMu_Match->Write();
 hNMuons_1m->Write();
 hPtMu_Match->Write();
 hNDimuons->Write();
 hNDimuons_2m->Write();

 hPercentileV0M->Write();
 hPercentileV0M_CMUL7->Write();
 hPercentileV0M_CINT7->Write();

 hAllpt_1M->Write();	 	
 hLpt_1M->Write();

 hAllpt_25y3->Write();
 hLpt_25y3->Write(); 
 hAllpt_3y35->Write(); 
 hLpt_3y35->Write();
 hAllpt_35y4->Write(); 
 hLpt_35y4->Write();

 	 	
 hAllpt_25y3_1M->Write();
 hLpt_25y3_1M->Write();
 hAllpt_3y35_1M->Write();
 hLpt_3y35_1M->Write(); 	 	
 hAllpt_35y4_1M->Write();	 	
 hLpt_35y4_1M->Write(); 
 	 	
 hAllpt_010_1M->Write();
 hLpt_010_1M->Write();
 hAllpt_1020_1M->Write();
 hLpt_1020_1M->Write();
 hAllpt_2030_1M->Write();
 hLpt_2030_1M->Write();
 hAllpt_3040_1M->Write();
 hLpt_3040_1M->Write();
 hAllpt_4050_1M->Write();
 hLpt_4050_1M->Write();
 hAllpt_5060_1M->Write();
 hLpt_5060_1M->Write();
 hAllpt_6070_1M->Write();
 hLpt_6070_1M->Write();
 hAllpt_7080_1M->Write();
 hLpt_7080_1M->Write();
 hAllpt_8090_1M->Write();
 hLpt_8090_1M->Write();

 fout->Close();	
 printf("Writing file %s\n",fnameout);
 
}
