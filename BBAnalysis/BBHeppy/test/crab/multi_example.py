####################
# #
#Summer 15 #
# #
####################

dataset = {
 'QCD15': '/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD30': '/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD50': '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD80': '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD120': '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW', 
 'QCD170': '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD300': '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCD470': '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',

# 'QCD600': '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD800': '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD1000': '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD1400': '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD1800': '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD2400': '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',
# 'QCD3200': '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_20_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',


 'QCDEM15': '/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDEM20': '/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDEM30': '/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDEM50': '/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDEM80': '/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDEM120': '/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',

'QCDMu15': '/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDMu20': '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDMu30': '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDMu50': '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDMu80': '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'QCDMu120': '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',

'DYToLL': '/DYToLL_M_1_TuneCUETP8M1_13TeV_pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
 'WJets': '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW'
}

nfiles = {
 'QCD15': -1,
 'QCD30': -1,
 'QCD50': -1,
 'QCD80': -1,
 'QCD120': -1,
 'QCD170': -1,
 'QCD300': -1,
 'QCD470': -1,
 'QCD600': -1,
 'QCD800': -1,
 'QCD1000': -1,
 'QCD1400': -1,
 'QCD1800': -1,
 'QCD2400': -1,
 'QCD3200': -1,

'QCDEM15': -1,
 'QCDEM20': -1,
 'QCDEM30': -1,
 'QCDEM50': -1,
 'QCDEM80': -1,
 'QCDEM120': -1,

'QCDMu15': -1,
 'QCDMu20': -1,
 'QCDMu30': -1,
 'QCDMu50': -1,
 'QCDMu80': -1,
 'QCDMu120': -1,
 'QCDMu120': -1,
 
 'DYToLL': -1,
 'WJets': -1

}

filesPerJob = {
 'QCD15': 75,
 'QCD30': 140,
 'QCD50': 130,
 'QCD80': 45,
 'QCD120': 19,
 'QCD170': 10,
 'QCD300': 7,
 'QCD470': 7,
 'QCD600': 7,
 'QCD800': 7,
 'QCD1000': 6,
 'QCD1400': 5,
 'QCD1800': 3,
 'QCD2400': 3,
 'QCD3200': 3,

'QCDEM15': 100,
 'QCDEM20': 140,
 'QCDEM30': 140,
 'QCDEM50': 10,
 'QCDEM80': 10,
 'QCDEM120': 5,

'QCDMu15': 45,
 'QCDMu20': 45,
 'QCDMu30': 45,
 'QCDMu50': 10,
 'QCDMu80': 5,
 'QCDMu120': 5,

'DYToLL': 11,
 'WJets': 11
 
}

if __name__ == '__main__':
 from CRABAPI.RawCommand import crabCommand

def submit(config):
 res = crabCommand('submit', config = config)

from CRABClient.UserUtilities import config
 config = config()
 name = 'HLTRates_2e33_25ns_V4p4_V1'
 config.General.workArea = 'crab_'+name
 config.General.transferLogs = True
# config.General.transferOutputs = True
 config.JobType.pluginName = 'Analysis'
 config.JobType.psetName = 'hlt.py'
 config.Data.inputDBS = 'global'
 config.Data.splitting = 'FileBased'
# config.Data.publication = True
 config.Data.publication = False
 config.JobType.outputFiles = ['hltbits.root'] #,'DQMIO.root']
 config.Site.storageSite = 'T2_CH_CERN'
 
# listOfSamples = ['QCDEM15','QCDEM20','QCDEM30','QCDEM50','QCDEM80','QCDEM120','QCDMu15','QCDMu20','QCDMu30','QCDMu50','QCDMu80','QCDMu120','QCD15','QCD30','QCD50','QCD80','QCD120','QCD170','QCD300','QCD470','QCD600','QCD800','QCD1000','QCD1400','QCD1800','QCD2400','QCD3200','DYToLL','WJets']
 listOfSamples = ['QCDEM15','QCDEM20','QCDEM30','QCDEM50','QCDEM80','QCDEM120','QCDMu15','QCDMu20','QCDMu30','QCDMu50','QCDMu80','QCDMu120','QCD15','QCD30','QCD50','QCD80','QCD120','QCD170','QCD300','QCD470','DYToLL','WJets']
 listOfSamples.reverse()
 for sample in listOfSamples:
 config.General.requestName = sample
 config.Data.inputDataset = dataset[sample]
 config.Data.unitsPerJob = filesPerJob[sample]
 config.Data.totalUnits = nfiles[sample]
 config.Data.outputDatasetTag = sample
 config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/Spring15/' + name # or '/store/group/<subdir>'
 submit(config)
