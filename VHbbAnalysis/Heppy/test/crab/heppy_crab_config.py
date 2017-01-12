from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'VHBB_A25_004'
config.General.workArea = 'crab_projects_A25_004'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
#os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")

os.system("tar czf python.tar.gz --directory $CMSSW_BASE python `find $CMSSW_BASE/src -name python | perl -pe s#$CMSSW_BASE/## `")
#onfig.JobType.sendPythonFolder = True
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
                             'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                              '../vhbb_combined.py',
                             'TMVAClassification_BDT.weights.xml',
                             'puData.root',
                             'puDataMinus.root',
                             'puDataPlus.root',
                             '../triggerEmulation.root',
                             'puMC.root',
                              'json.txt',
                              #"../Zll-spring15.weights.xml",
                              #"../Wln-spring15.weights.xml",
                              #"../Znn-spring15.weights.xml",
                              #"../VBF-spring15.weights.xml",
                              #"../ttbar-fall15_TargetGenOverPt_GenPtCut0.weights.xml",
			      #'../ttbar-spring16-80X.weights.xml',
			      '../ttbar-spring16-500k-13d-300t.weights.xml',	
			      '../TMVA_blikelihood_vbf_cmssw76_h21trained.weights.xml'
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/ZH_HToBB_ZToLL_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 20
config.Data.allowNonValidInputDataset = True # to run on datasets in PRODUCTION
config.Data.outLFNDirBase = '/store/user/arizzi/VHBBHeppyA25/'
config.Data.publication = True
config.Data.outputDatasetTag = 'VHBB_HEPPY_A25'

config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"
#config.Site.storageSite = "T3_CH_PSI"

#config.Data.ignoreLocality = True
