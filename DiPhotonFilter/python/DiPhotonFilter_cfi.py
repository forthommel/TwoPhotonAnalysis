import FWCore.ParameterSet.Config as cms

diphotonCandidates = cms.EDFilter(
    "DiPhotonFilter",
    #cut = cms.string("isGlobalMuon = 1 & pt > 5.0"),
    #src = cms.InputTag("selectedLayer1Muons"),
    #filter = cms.bool(True)
    photonTag = cms.InputTag('slimmedPhotons'),
    vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    minimalDiPhotonMass = cms.double(-1),
    maximalDiPhotonMass = cms.double(-1),
    maximalVertexDistance = cms.double(2.),
)
