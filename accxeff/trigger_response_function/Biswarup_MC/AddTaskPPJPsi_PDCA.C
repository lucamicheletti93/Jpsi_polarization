AliAnalysisTaskPPJPsi_PDCA *AddTaskPPJPsi_PDCA(const char *kAnalysisType, TString run, Bool_t usePhysicsSelection){

//****************************************************************************************
// Add task class.
// The attached class prepares and draws some kinematical distributions of muons/dimuons
// Roberta
//****************************************************************************************

   printf("Creating Task for Muon/Dimuon Histos\n");
    
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPbPbJPsi", "No analysis manager to connect to.");
      return NULL;
   }   
   TString fnameout;
   fnameout = "PbPb2015_"+run+".root";
   //else if(!usePhysicsSelection)fnameout = "PP_noPS_"+run+"_PDCA.root";
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0",TList::Class(),AliAnalysisManager::kOutputContainer,fnameout);

  /* AliMuonPairCuts* muonPairCuts = new AliMuonPairCuts("StandardMuonPairCuts", "StandardMuonPairCuts");
   muonPairCuts->SetFilterMask(AliMuonPairCuts::kBothMuPdca);
   AliAnalysisTaskPPJPsi_PDCA *PPJPsiTask = new AliAnalysisTaskPPJPsi_PDCA("AliAnalysisTaskPPJPsi_PDCA",*muonPairCuts);
   PPJPsiTask->SetAnalysisType(kAnalysisType);*/

   AliAnalysisTaskPPJPsi_PDCA *PPJPsiTask = new AliAnalysisTaskPPJPsi_PDCA("AliAnalysisTaskPPJPsi_PDCA");
   PPJPsiTask->SetAnalysisType(kAnalysisType);
   //
   // define by hand the beam energy
   //
   PPJPsiTask->SetBeamEnergy(5020.);
   mgr->AddTask(PPJPsiTask);

   mgr->ConnectInput(PPJPsiTask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(PPJPsiTask,1,coutput1);
  
   return PPJPsiTask;
}
