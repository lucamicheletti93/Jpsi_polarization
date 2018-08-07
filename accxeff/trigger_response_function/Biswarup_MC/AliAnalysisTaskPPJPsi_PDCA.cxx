/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskPPJPsi_PDCA.cxx 48309 2011-03-10 21:32:21Z martinez $ */

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskPPJPsi_PDCA_CXX
//#define AliAnalysisTaskPPJPsi_PDCA_CXX

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include "AliInputEventHandler.h"
#include "AliCentrality.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "/home/biswarup/alice/aliphysics/master/inst/include/AliMuonTrackCuts.h"
#include "/home/biswarup/alice/aliphysics/master/inst/include/AliMuonPairCuts.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliAnalysisTaskPPJPsi_PDCA.h"

ClassImp(AliAnalysisTaskPPJPsi_PDCA)
//__________________________________________________________________________
AliAnalysisTaskPPJPsi_PDCA::AliAnalysisTaskPPJPsi_PDCA():
  AliAnalysisTaskSE(),
  //fMuonPairCuts(),  
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fCountDimu(0x0),
  fCountDimu_CMUU(0x0),
  fCountMB(0x0),
  fCountDimu_CMUL(0x0),
  fCountTotEv(0x0),
  fCountDimu_LL(0x0),
  fCountDimu_ULLL(0x0),
  fCountCPBI1MSL(0x0),
  fCountCPBI2(0x0),
  fCountCVHN(0x0),
  fCountCVLN(0x0),
  fAODEvent(0x0),
  fOutput(0x0)
  
{
  /// Default ctor.
}

//__________________________________________________________________________
AliAnalysisTaskPPJPsi_PDCA::AliAnalysisTaskPPJPsi_PDCA(const char *name) :
  AliAnalysisTaskSE(name),
  //fMuonPairCuts(),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fCountDimu(0x0),
  fCountDimu_CMUU(0x0),
  fCountMB(0x0),
  fCountDimu_CMUL(0x0),
  fCountTotEv(0x0),
  fCountDimu_LL(0x0),
  fCountDimu_ULLL(0x0),
  fCountCPBI1MSL(0x0),
  fCountCPBI2(0x0),
  fCountCVHN(0x0),
  fCountCVLN(0x0),
  fAODEvent(0x0),
  fOutput(0x0)
  
{
 //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskPPJPsi_PDCA","Calling Constructor");
  
  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskPPJPsi_PDCA& AliAnalysisTaskPPJPsi_PDCA::operator=(const AliAnalysisTaskPPJPsi_PDCA& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskPPJPsi_PDCA::AliAnalysisTaskPPJPsi_PDCA(const AliAnalysisTaskPPJPsi_PDCA& c) :
  AliAnalysisTaskSE(c),
  //fMuonPairCuts(c.fMuonPairCuts),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fCountDimu(c.fCountDimu),
  fCountDimu_CMUU(c.fCountDimu_CMUU),
  fCountMB(c.fCountMB),
  fCountDimu_CMUL(c.fCountDimu_CMUL),
  fCountTotEv(c.fCountTotEv),
  fCountDimu_LL(c.fCountDimu_LL),
  fCountDimu_ULLL(c.fCountDimu_ULLL),
  fCountCPBI1MSL(c.fCountCPBI1MSL),
  fCountCPBI2(c.fCountCPBI2),
  fCountCVHN(c.fCountCVHN),
  fCountCVLN(c.fCountCVLN),
  fAODEvent(c.fAODEvent),
  fOutput(c.fOutput)
 {
  //
  // Copy Constructor										
  //
}

//___________________________________________________________________________
AliAnalysisTaskPPJPsi_PDCA::~AliAnalysisTaskPPJPsi_PDCA() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskPPJPsi_PDCA","Calling Destructor");
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {

  if (fOutput){
    delete fOutput;
    fOutput = 0;
   }
  }
}

//___________________________________________________________________________
void AliAnalysisTaskPPJPsi_PDCA::NotifyRun() {
  //if ( fMuonPairCuts.GetFilterMask() ) fMuonPairCuts.SetRun(fCurrentRunNumber); // Run number needed to get parameters from OADB
  //  if ( fMuonPairCuts.GetFilterMask() ) fMuonPairCuts.SetPassNumber(2);
  //if(fMuonPairCuts.GetFilterMask() ) fMuonPairCuts.SetRun(fInputHandler);
  //fMuonPairCuts.SetRun(fInputHandler);
    //if ( fMuonPairCuts.GetFilterMask() ) fMuonPairCuts.SetRun(fInputEvent);
}
//___________________________________________________________________________
void AliAnalysisTaskPPJPsi_PDCA::UserCreateOutputObjects(){
 //
 // output objects creation
 //	 
 fOutput = new TList();
 fOutput->SetOwner(); 



/*TH1D *hNEv = new TH1D("hNEv","hNEv",4,0.,4.);
TString namelabel[4]={"Dimu_CMUU","MB","Dimu_CMUL","Dimu_CMUUCMUL"};
for(int kk=0;kk<4;kk++) hNEv->GetXaxis()->SetBinLabel(kk+1,namelabel[kk]);*/


TH1D *hAllpt_25y3 = new TH1D("hAllpt_25y3","hAllpt_25y3",500,0.,50.);     
TH1D *hAllpt_25y3_1M = new TH1D("hAllpt_25y3_1M","hAllpt_25y3_1M",500,0.,50.);

TH1D *hLpt_25y3 = new TH1D("hLpt_25y3","hLpt_25y3",500,0.,50.);      
TH1D *hLpt_25y3_1M = new TH1D("hLpt_25y3_1M","hLpt_25y3_1M",500,0.,50.);

TH1D *hAllpt_3y35 = new TH1D("hAllpt_3y35","hAllpt_3y35",500,0.,50.);     
TH1D *hAllpt_3y35_1M = new TH1D("hAllpt_3y35_1M","hAllpt_3y35_1M",500,0.,50.);

TH1D *hLpt_3y35 = new TH1D("hLpt_3y35","hLpt_3y35",500,0.,50.);      
TH1D *hLpt_3y35_1M = new TH1D("hLpt_3y35_1M","hLpt_3y35_1M",500,0.,50.);

TH1D *hAllpt_35y4 = new TH1D("hAllpt_35y4","hAllpt_35y4",500,0.,50.);     
TH1D *hAllpt_35y4_1M = new TH1D("hAllpt_35y4_1M","hAllpt_35y4_1M",500,0.,50.);

TH1D *hLpt_35y4 = new TH1D("hLpt_35y4","hLpt_35y4",500,0.,50.);      
TH1D *hLpt_35y4_1M = new TH1D("hLpt_35y4_1M","hLpt_35y4_1M",500,0.,50.);
 
 // stat histos
 TH1D *hNDimuons = new TH1D("hNDimuons","hNDimuons",2000,0.,2000.); 

	
 fOutput->Add(hAllpt_25y3); 	 	
 fOutput->Add(hAllpt_25y3_1M);

 fOutput->Add(hLpt_25y3); 	 	
 fOutput->Add(hLpt_25y3_1M);

 fOutput->Add(hAllpt_3y35); 	 	
 fOutput->Add(hAllpt_3y35_1M);

 fOutput->Add(hLpt_3y35); 	 	
 fOutput->Add(hLpt_3y35_1M);

 fOutput->Add(hAllpt_35y4); 	 	
 fOutput->Add(hAllpt_35y4_1M);

 fOutput->Add(hLpt_35y4); 	 	
 fOutput->Add(hLpt_35y4_1M);
 
 fOutput->Add(hNDimuons);	
 fOutput->ls(); 

 PostData(1,fOutput);
} 

//___________________________________________________________
void AliAnalysisTaskPPJPsi_PDCA::UserExec(Option_t *)
{
//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }
 // Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
 // if(isSelected){  
   // to be modified
   // 1 GeV bin 
  /*const int nbinsPt = 20;    // 20 bins
  Double_t PtMin[nbinsPt]={0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.};  //test
  Double_t PtMax[nbinsPt]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.};
  const int nbinsY = 6;
  Double_t YMin[nbinsY]={-4.,-3.75,-3.5,-3.25,-3.,-2.75};
  Double_t YMax[nbinsY]={-3.75,-3.5,-3.25,-3.,-2.75,-2.5};*/
  
  //char hname[200];
  // char hnamec[200];
 
/*   fCountTotEv++;
  AliAODHeader *aodheader =  fAODEvent->GetHeader();
  TString firedtrigger = aodheader->GetFiredTriggerClasses();
  
  Bool_t TriggerSelected=kFALSE;
  Bool_t TriggerSelected_CMUU=kFALSE;
  Bool_t TriggerSelected_MB=kFALSE;
  Bool_t TriggerSelected_CMUL=kFALSE; 

  if(firedtrigger.Contains("CMUU7-B-NOPF-ALLNOTRD") || firedtrigger.Contains("CMUU7-B-NOPF-MUON")||firedtrigger.Contains("CMUL7-B-NOPF-MUON")) TriggerSelected = kTRUE; 
  else TriggerSelected = kFALSE; 

  if(firedtrigger.Contains("CMUU7-B-NOPF-ALLNOTRD") || firedtrigger.Contains("CMUU7-B-NOPF-MUON")) TriggerSelected_CMUU = kTRUE; 
  else TriggerSelected_CMUU = kFALSE; 

if(firedtrigger.Contains("CINT1-B-NOPF-ALLNOTRD") || firedtrigger.Contains("CINT7-B-NOPF-ALLNOTRD")||firedtrigger.Contains("CINT7-I-NOPF-ALLNOTRD")) TriggerSelected_MB = kTRUE;
  else TriggerSelected_MB = kFALSE;

 if(firedtrigger.Contains("CMUL7-B-NOPF-MUON")) TriggerSelected_CMUL = kTRUE; 
  else TriggerSelected_CMUL = kFALSE; 

  fCountTotEv++;

  if(TriggerSelected) fCountDimu++;
  if(TriggerSelected_CMUU) fCountDimu_CMUU++;
  if(TriggerSelected_MB) fCountMB++;
  if(TriggerSelected_CMUL) fCountDimu_CMUL++;*/

  //if(TriggerSelected){
     ((TH1D*)(fOutput->FindObject("hNDimuons")))->Fill(fAODEvent->GetNumberOfDimuons());
    
    for(int nd=0;nd<fAODEvent->GetNumberOfDimuons();nd++){
       AliAODDimuon *dimu = dynamic_cast<AliAODDimuon*>(fAODEvent->GetDimuon(nd));
       AliAODTrack *mu0 = dimu->GetMu(0); 
       AliAODTrack *mu1 = dimu->GetMu(1); 
       /*TList *ldimu = new TList();
       AliVParticle *vmu0 = 0x0;
       AliVParticle *vmu1 = 0x0;
       vmu0 = dimu->GetMu(0); 
       vmu1 = dimu->GetMu(1); 
       ldimu->Add(vmu0);
       ldimu->Add(vmu1);*/
      
       Double_t DimuMass=999;
       Double_t DimuPt=-999;
       Double_t DimuY=-999;
       Double_t Match_Mu0=-999;
       Double_t Match_Mu1=-999;
       Double_t Pt_Mu0=-999;
       Double_t Pt_Mu1=-999;
       Double_t Eta_Mu0=-999;
       Double_t Eta_Mu1=-999;
       Double_t RAbs_Mu0=-999;
       Double_t RAbs_Mu1=-999;	     
       
       Match_Mu0=mu0->GetMatchTrigger();
       Match_Mu1=mu1->GetMatchTrigger();
       Eta_Mu0=mu0->Eta();
       Eta_Mu1=mu1->Eta();
       RAbs_Mu0=mu0->GetRAtAbsorberEnd();
       RAbs_Mu1=mu1->GetRAtAbsorberEnd();	
       Pt_Mu0=mu0->Pt();
       Pt_Mu1=mu1->Pt();

         DimuMass=dimu->Mass();
//	 if (DimuMass<1.8) continue;
           if((Eta_Mu0>-4 && Eta_Mu0<-2.5) && (Eta_Mu1>-4 && Eta_Mu1<-2.5)){
             if((RAbs_Mu0>17.6 && RAbs_Mu0<89.5) && (RAbs_Mu1>17.6 && RAbs_Mu1<89.5)){
               DimuY=dimu->Y();
//               if(DimuY>-4 && DimuY<-2.5){
                  if(dimu->Charge()==0) {
		    DimuPt=dimu->Pt();
		    if(DimuPt>0 && DimuPt<50){ 	

               if(DimuY>-3 && DimuY<-2.5){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_25y3_1M")))->Fill(Pt_Mu0);
                                        }
		       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_25y3_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0 >=1 && Match_Mu1>=1){  
			    ((TH1D*)(fOutput->FindObject("hAllpt_25y3")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hAllpt_25y3")))->Fill(Pt_Mu1);
                                     } 
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_25y3_1M")))->Fill(Pt_Mu0);
                                     }
		       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_25y3_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0>=2 && Match_Mu1>=2){  
			    ((TH1D*)(fOutput->FindObject("hLpt_25y3")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hLpt_25y3")))->Fill(Pt_Mu1);
                                     } 
                                   }//2.5y3

               if(DimuY>-3.5 && DimuY<-3){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_3y35_1M")))->Fill(Pt_Mu0);
                                        }
		       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_3y35_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0 >=1 && Match_Mu1>=1){  
			    ((TH1D*)(fOutput->FindObject("hAllpt_3y35")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hAllpt_3y35")))->Fill(Pt_Mu1);
                                     } 
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_3y35_1M")))->Fill(Pt_Mu0);
                                     }
		       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_3y35_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0>=2 && Match_Mu1>=2){  
			    ((TH1D*)(fOutput->FindObject("hLpt_3y35")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hLpt_3y35")))->Fill(Pt_Mu1);
                                     } 
                                   }//3y3.5

               if(DimuY>-4 && DimuY<-3.5){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_35y4_1M")))->Fill(Pt_Mu0);
                                        }
		       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_35y4_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0 >=1 && Match_Mu1>=1){  
			    ((TH1D*)(fOutput->FindObject("hAllpt_35y4")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hAllpt_35y4")))->Fill(Pt_Mu1);
                                     } 
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_35y4_1M")))->Fill(Pt_Mu0);
                                     }
		       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_35y4_1M")))->Fill(Pt_Mu1);                                               
                                  }
                        if(Match_Mu0>=2 && Match_Mu1>=2){  
			    ((TH1D*)(fOutput->FindObject("hLpt_35y4")))->Fill(Pt_Mu0);
			    ((TH1D*)(fOutput->FindObject("hLpt_35y4")))->Fill(Pt_Mu1);
                                     } 
                                   }//3.5y4


                         }//Dimuon Pt loop 
                      }//Dimuon charge 
//                   }//Y loop
                 }//muon Rabs 
                }//muon eta
              }//dimuon loop
            
          //} 
 

 /*  ((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(1,(Double_t)fCountDimu_CMUU);
   ((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(2,(Double_t)fCountMB);
   ((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(3,(Double_t)fCountDimu_CMUL);
   ((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(4,(Double_t)fCountDimu);*/
//}
PostData(1,fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskPPJPsi_PDCA::Terminate(Option_t *) 
{
//
// Draw histos
//
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  //Int_t xmin=20; 
  // Int_t ymin=20;
  
  printf("Using beam Energy=%f \n",fBeamEnergy);
  
 }
