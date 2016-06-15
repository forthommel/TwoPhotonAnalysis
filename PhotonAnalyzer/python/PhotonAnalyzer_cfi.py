import FWCore.ParameterSet.Config as cms

lightPhotonAnalyzerNoPhId = cms.EDAnalyzer('PhotonAnalyzer',
	photonTag = cms.InputTag('slimmedPhotons'),
	beamspotTag = cms.InputTag('offlineBeamSpot'),
	conversionTag = cms.InputTag('reducedEgamma', 'reducedConversions'),
        pfCandidateTag = cms.InputTag('packedPFCandidates'),
	vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
	#photonMedIdBoolMapTag = cms.untracked.InputTag('egmPhotonIDs:mvaPhoID-Spring15-50ns-nonTrig-V2-wp90'),
	#photonMedIdInfoMapTag = cms.untracked.InputTag('egmPhotonIDs:mvaPhoID-Spring15-50ns-nonTrig-V2-wp90'),
	#photonMedIdBoolMapTag = cms.untracked.InputTag('egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90'),
	#photonMedIdInfoMapTag = cms.untracked.InputTag('egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90'),
)

from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import photonIDValueMapProducer
from RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi import photonMVAValueMapProducer

lightPhotonAnalyzer = cms.Sequence(photonIDValueMapProducer * photonMVAValueMapProducer * lightPhotonAnalyzerNoPhId)
