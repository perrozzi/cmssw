####################
# #
#Summer 15 #
# #
####################

#datasets = [
#"/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#"/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
##"/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-25nsNoPURaw_magnetOn_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
##"/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
##"/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_magnetOn_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
##"/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM",
#]

#datasets = [
#"/ZeroBias/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias1/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias2/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias3/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias4/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias5/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias6/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias7/Run2015D-16Dec2015-v1/MINIAOD",
#"/ZeroBias8/Run2015D-16Dec2015-v1/MINIAOD",
#]


datasets = [
"/JetHT/Run2015D-16Dec2015-v1/MINIAOD",
]

replacePatterns=[
("pythia","Py"),
("76X_mcRun2_asymptotic","76r2as"),
("RunIIFall15MiniAODv2","fall15MAv2"),
("PU25nsData2015","pu25ns15")
]

if __name__ == '__main__':
 from CRABAPI.RawCommand import crabCommand

def submit(config):
 res = crabCommand('submit', config = config)

import re
import copy
import heppy_crab_config
config0 = heppy_crab_config.config
for dataset  in datasets:
     config=copy.deepcopy(config0)
     m=re.match("\/(.*)\/(.*)\/(.*)",dataset)
     if not m : 
        print "NOT A GOOD DATASET", dataset
        continue
     sample=m.group(1)+"__"+m.group(2)
     for (s,r) in replacePatterns :
         sample=re.sub(s,r,sample)
     config.General.requestName+= "_"+sample
     if len(config.General.requestName) > 100 :
         config.General.requestName=config.General.requestName[:90]+config.General.requestName[-10:] 
     config.Data.inputDataset = dataset
     config.Data.outputDatasetTag += "_"+sample
     
     print "======== SUBMIT " , config.General.requestName , "=========="
     print "outdir",config.Data.outLFNDirBase
     try:
       submit(config)
     except:
        print "FAILED",config.General.requestName
