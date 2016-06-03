// -*- C++ -*-
//
// Package:    TwoPhoton/PhotonAnalyzer
// Class:      PhotonAnalyzer
// 
/**\class PhotonAnalyzer PhotonAnalyzer.cc TwoPhoton/PhotonAnalyzer/plugins/PhotonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Sun, 29 May 2016 19:05:16 GMT
//
//

#include "PhotonAnalyzer.h"

//
// constructors and destructor
//
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig) :
  fPhotonToken(consumes< edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"))),
  //fBeamSpotToken(consumes< edm::View<reco::BeamSpot> >(iConfig.getParameter<edm::InputTag>("beamspotTag")))
  fTree(0)
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TFileDirectory dir = fs->mkdir("DiPhotonAnalyzer");

  fTree = fs->make<TTree>("ntp1", "diphoton candidates tree");
}


PhotonAnalyzer::~PhotonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called for each event  ------------
void
PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByToken(fPhotonToken, photons);

  /*edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(fVertexToken, vertices);*/
 
  unsigned int i=0, j=0;
  for (edm::View<pat::Photon>::const_iterator photon=photons->begin(); photon!=photons->end(); photon++) {
    //if (!passElectronId())
    aPhotonPt[i] = photon->pt();
    aPhotonEta[i] = photon->eta();
    aPhotonPhi[i] = photon->phi();
    aPhotonVtxX[i] = photon->vx();
    aPhotonVtxY[i] = photon->vy();
    aPhotonVtxZ[i] = photon->vz();

    for (edm::View<pat::Photon>::const_iterator photon2=photon+1; photon2!=photons->end(); photon2++) {
      const reco::Candidate::LorentzVector ppair = photon->p4()+photon2->p4();
      aPhotonPairMass[j] = ppair.mass();
      aPhotonPairPt[j] = ppair.pt();
      j++;
    }

    i++;
  }
  aNumPhotons = i;
  aNumPhotonPairs = j;

  /*for (edm::View<reco::Vertex>::const_iterator vtx=vertices->begin(); vtx!=vertices->end(); vtx++) {
    if (passVertexCriteria(*vtx, *photon1, *photon2)) matched_vertex = true;
  }*/

  fTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonAnalyzer::beginJob()
{
  fTree->Branch("num_photons", &aNumPhotons, "num_photons/I");
  fTree->Branch("photon_pt", aPhotonPt, "photon_pt[num_photons]/D");
  fTree->Branch("photon_eta", aPhotonEta, "photon_eta[num_photons]/D");
  fTree->Branch("photon_phi", aPhotonPhi, "photon_phi[num_photons]/D");
  fTree->Branch("photon_vtxx", aPhotonVtxX, "photon_vtxx[num_photons]/D");
  fTree->Branch("photon_vtxy", aPhotonVtxY, "photon_vtxy[num_photons]/D");
  fTree->Branch("photon_vtxz", aPhotonVtxZ, "photon_vtxz[num_photons]/D");
  fTree->Branch("num_photonpairs", &aNumPhotonPairs, "num_photonpairs/I");
  fTree->Branch("photonpair_mass", aPhotonPairMass, "photonpair_mass[num_photonpairs]/D");
  fTree->Branch("photonpair_pt", aPhotonPairPt, "photonpair_pt[num_photonpairs]/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
//DEFINE_FWK_MODULE(PhotonAnalyzer);
