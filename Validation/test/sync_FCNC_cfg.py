import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
#    '/store/user/jhgoh/CATTools/sync/v7-6-5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',
    'file:/state/partition1/store/user/jhgoh/FCNC/Synchronization/201610/v801/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',
]

#process.load("CATTools.CatAnalyzer.filters_cff")
#process.load("CATTools.Validation.ttllEventSelector_cfi")
#process.load("CATTools.Validation.validation_cff")
process.load("CATTools.Validation.eventsTopFCNC_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(process.eventsTopFCNC)

