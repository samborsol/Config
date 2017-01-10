
#void forest2diJetSkim(
#          TString fname = "/home/kbg777/CMSwork/4_UPCTriggers_Beomgon.root",
#	  TString outputFname = "upcDiJetSkim",
#	  TString trig = "",  # if trig=="" then the trigger cut is not functioning
#	  TString jetCollection = "ak5PFJetAnalyzer", or  "akPu5PFJetAnalyzer",
#	  float minjPt = 30,
#	  int nevt=-1  ) 

# pp di-jet data 
#Input file is at MIT server 
inputFile="/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root"
jetCollection="ak5PFJetAnalyzer"
trig=""

for  minjPt in 20 30 
do
    root -l -b -q 'forest2diJetSkim.C+("'$inputFile'","upcDijetSkim","'$trig'","'$jetCollection'",'$minjPt',-1)'
done

# pp PYTHIA MC :
# input file is at XXX

# PbPb UPC triggered event 

# PbPb high pt jet triggered event 

# STARLIGHT DPMJET event


