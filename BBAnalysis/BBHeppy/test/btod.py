import FWCore.ParameterSet.Config as cms

process = cms.Process("BTOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(       
        #'file:QCD_Pt-120to170.root',
        #'file:QCD_Pt-80to120.root',
        'file:QCD_Pt-1800to2400.root'
    )#,
    #skipEvents=cms.untracked.uint32(70)
)

process.ProducerBToD = cms.EDProducer('ProducerBToD',
)

process.OUT = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('test.root'),
        outputCommands = cms.untracked.vstring(['drop *','keep *_*_*_BTOD'])
    )


process.p = cms.Path(process.ProducerBToD)
process.endpath= cms.EndPath(process.OUT)

