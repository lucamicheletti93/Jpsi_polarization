void RunLocalPP2015(TString runnumber = "999", 
		      Bool_t usePhysicsSelection = kTRUE,
                      TString datadir = "/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/With_Y_Charge_cuts", 
                      TString workdir = "/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/With_Y_Charge_cuts", 
		      Int_t nev=999999999,
		      Int_t firstev=0)

//========================================================================
{
    printf("\n#################################################\n"); 
    printf("Running local train with following wagons:\n\n"); 
    if(usePhysicsSelection) printf("\n Using Physics Selection\n");
    else printf("\n Not using Physics Selection\n");
    printf("#################################################\n\n"); 
    
    
    //=====================================================================
    // Global configuration flags 
    //=====================================================================
    Bool_t debug         = kFALSE;
    Bool_t useMUONlib    = kFALSE;

    //=====================================================================
    // Load common and specific libraries 
    //=====================================================================
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");

    // Common packages
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");

    //gSystem->Load("libPWG3base");
    gSystem->Load("libPWGmuon");
    
    //gROOT->ProcessLine(".L $ALICE_ROOT/PWG3/muon/AliMuonPairCuts.cxx+g");

    //=====================================================================
    // Create the chain. 
    //=====================================================================
    Int_t firstchunk=1;
    Int_t nchunk=100;
    
    TChain *chain = new TChain("aodTree");
    char file[200];
    for(int jc=firstchunk;jc<=nchunk;jc++){

if(jc<10)
sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%s/00%d/AliAOD.Muons.root",runnumber.Data(),jc);
     else if(jc>=10 && jc<100)
sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%s/0%d/AliAOD.Muons.root",runnumber.Data(),jc);
     else
 if(jc>100)
sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%s/%d/AliAOD.Muons.root",runnumber.Data(),jc);

//    sprintf(file,"/Data/Official_AOD_2011/LHC11c/%s/000%d/AliAOD.root",runnumber.Data(),jc);
      Long_t *dummy1 =0, *dummy2 =0, *dummy3 =0, *dummy4 =0;
      if(gSystem->GetPathInfo(file,dummy1,dummy2,dummy3, dummy4) == 0) {
        printf(" Opening %s\n",file);
        chain->Add(file);
      }  
    }   
    printf("entries %d\n",chain->GetEntries());  
    //=====================================================================
    // Create the train and set-up the handlers
    //=====================================================================
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis Train");

    // to save plots produced in terminate
    mgr->SetSaveCanvases(kFALSE);
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
       
    // Debugging if requested
    //if (debug) mgr->SetDebugLevel(3);
    
    // load the physics selection
    /*gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask; 
    if(usePhysicsSelection) physSelTask = AddTaskPhysicsSelection(); */
     
    // load the analysis wagon 
    gROOT->ProcessLine(".L AliAnalysisTaskPPJPsi_PDCA.cxx++");
    gROOT->LoadMacro("AddTaskPPJPsi_PDCA.C");
    AliAnalysisTaskPPJPsi_PDCA *PPJPsiTask = AddTaskPPJPsi_PDCA("AOD",runnumber,usePhysicsSelection);
    if(usePhysicsSelection) {};
//PPJPsiTask->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB | AliVEvent::kINT7);
//PPJPsiTask->SelectCollisionCandidates(AliVEvent::kAny); 
 
   if (!mgr->InitAnalysis()) return;
   
    mgr->StartAnalysis("local",chain,nev,firstev);

}  
