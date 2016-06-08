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
//         Created:  Sun, 22 May 2016 19:05:16 GMT
//
//

#include "PhotonAnalyzer.h"

//
// constructors and destructor
//
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig) :
  fPhotonToken(consumes< edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"))),
  fConversionToken(consumes< edm::View<reco::Conversion> >(iConfig.getParameter<edm::InputTag("conversionTag")>)),
  //fBeamSpotToken(consumes< edm::View<reco::BeamSpot> >(iConfig.getParameter<edm::InputTag>("beamspotTag")))
  fTree(0)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TFileDirectory dir = fs->mkdir("DiPhotonAnalyzer");

  fTree = fs->make<TTree>("ntp1", "diphoton candidates tree");
}


PhotonAnalyzer::~PhotonAnalyzer()
{}


// ------------ method called for each event  ------------
void
PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByToken(fPhotonToken, photons);

  edm::Handle< edm::View<reco::Conversion> > conversions;
  iEvent.getByToken(fConversionToken, conversions);

  /*edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(fVertexToken, vertices);*/

  TLorentzVector ph1, ph2;
 
  unsigned int i=0, j=0;
  for (edm::View<pat::Photon>::const_iterator photon=photons->begin(); photon!=photons->end(); photon++) {
    if (!passPhotonId(*photon)) continue;

    aPhotonPt[i] = photon->pt();
    aPhotonEta[i] = photon->eta();
    aPhotonPhi[i] = photon->phi();
    aPhotonVtxX[i] = photon->vx();
    aPhotonVtxY[i] = photon->vy();
    aPhotonVtxZ[i] = photon->vz();

    ph1.SetPxPyPzE(photon->p4().px(), photon->p4().py(), photon->p4().pz(), photon->p4().E());

    for (edm::View<pat::Photon>::const_iterator photon2=photon+1; photon2!=photons->end(); photon2++) {
      if (!passPhotonId(*photon2)) continue;

      ph2.SetPxPyPzE(photon2->p4().px(), photon2->p4().py(), photon2->p4().pz(), photon2->p4().E());
      const TLorentzVector photon_pair = ph1+ph2;
      double dphi = fabs(ph1.Phi()-ph2.Phi());

      //aPhotonPairVertexDist[j] = std::sqrt((photon->vertex()-photon2->vertex()).mag2());
      aPhotonPairVertexDist[j] = std::sqrt(pow(photon->vx()-photon2->vx(), 2)+
                                           pow(photon->vy()-photon2->vy(), 2)+
                                           pow(photon->vz()-photon2->vz(), 2));
      aPhotonPairDphi[j] = (dphi<TMath::Pi()) ? dphi : 2.*TMath::Pi()-dphi; // dphi lies in [-pi, pi]
      aPhotonPairDpt[j] = fabs(ph1.Pt()-ph2.Pt());
      aPhotonPairMass[j] = photon_pair.M();
      aPhotonPairPt[j] = photon_pair.Pt();
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

bool
PhotonAnalyzer::passPhotonId(const pat::Photon& photon) const
{
  return true;
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
  fTree->Branch("photonpair_vertex_dist", aPhotonPairVertexDist, "photonpair_vertex_dist[num_photonpairs]/D");
  fTree->Branch("photonpair_mass", aPhotonPairMass, "photonpair_mass[num_photonpairs]/D");
  fTree->Branch("photonpair_pt", aPhotonPairPt, "photonpair_pt[num_photonpairs]/D");
  fTree->Branch("photonpair_dphi", aPhotonPairDphi, "photonpair_dphi[num_photonpairs]/D");
  fTree->Branch("photonpair_dpt", aPhotonPairDpt, "photonpair_dpt[num_photonpairs]/D");
  
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
