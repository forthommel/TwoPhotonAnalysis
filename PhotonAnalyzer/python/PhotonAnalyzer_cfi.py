import FWCore.ParameterSet.Config as cms

lightPhotonAnalyzer = cms.EDAnalyzer('PhotonAnalyzer',
	photonTag = cms.InputTag('slimmedPhotons'),
	beamspotTag = cms.InputTag('offlineBeamSpot'),
)
