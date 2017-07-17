from WMCore.Configuration import Configuration
version='V95JETHT'
config = Configuration()

config.section_("General")
config.General.requestName = 'BBHeppy_'+version
config.General.workArea = 'crab_projects_'+version+'_001'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")

config.JobType.inputFiles = ['heppy_config.py',  #do not remove
                             'heppy_crab_script.py', #do not remove
                             'python.tar.gz',
                             '../btod.py']
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/berger_p2/gbb/ntuples/BBHeppy'+version+'/'
config.Data.publication = False #True
config.Data.outputDatasetTag = 'BBHeppy_'+version
config.Data.allowNonValidInputDataset = True

config.section_("Site")
#config.Site.storageSite = "T2_IT_Legnaro"
config.Site.storageSite = "T2_IT_Pisa"
#config.Data.ignoreLocality = True
