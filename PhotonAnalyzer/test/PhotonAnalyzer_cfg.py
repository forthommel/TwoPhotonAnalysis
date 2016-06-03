import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
	'/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/000FAE50-82A6-E511-BC87-00261894397F.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('76X_dataRun2_16Dec2015_v0')
#process.GlobalTag.globaltag = cms.string('')


#from TwoPhoton.PhotonAnalyzer.PhotonAnalyzer_cfi import lightPhotonAnalyzer
process.load("TwoPhoton.PhotonAnalyzer.PhotonAnalyzer_cfi")

process.p = cms.Path(process.lightPhotonAnalyzer)
