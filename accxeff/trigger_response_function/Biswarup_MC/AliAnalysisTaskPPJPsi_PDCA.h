#ifndef ALIANALYSISTASKPPJPSI_PDCA_H
#define ALIANALYSISTASKPPJPSI_PDCA_H
#endif
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "/home/biswarup/alice/aliphysics/master/inst/include/AliMuonTrackCuts.h"
#include "/home/biswarup/alice/aliphysics/master/inst/include/AliMuonPairCuts.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
 
class AliAnalysisTaskPPJPsi_PDCA : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskPPJPsi_PDCA();
  AliAnalysisTaskPPJPsi_PDCA(const char *name);
  virtual ~AliAnalysisTaskPPJPsi_PDCA();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun(); 
  
  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
    
 private:
  AliAnalysisTaskPPJPsi_PDCA(const AliAnalysisTaskPPJPsi_PDCA&);
  AliAnalysisTaskPPJPsi_PDCA& operator=(const AliAnalysisTaskPPJPsi_PDCA&);
     
 protected:
 
  //AliMuonPairCuts fMuonPairCuts;   ///< Muon pair track cuts 
   
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)
    
  const char* fkAnalysisType; //ESD or AOD based analysis
  Int_t fCountDimu;
  Int_t fCountDimu_CMUU;    //!< counter
  Int_t fCountMB;
  Int_t fCountDimu_CMUL;
  Int_t fCountTotEv;   //!< counter
  Int_t fCountDimu_LL; //!< counter
  Int_t fCountDimu_ULLL; //!< counter
  Int_t fCountCPBI1MSL;//!< counter
  Int_t fCountCPBI2;   //!< counter
  Int_t fCountCVHN;    //!< counter
  Int_t fCountCVLN;    //!< counter
  AliAODEvent* fAODEvent;      //!< AOD event
  TList *fOutput;  //!< List of histograms for data
 
 ClassDef(AliAnalysisTaskPPJPsi_PDCA,1);
};

