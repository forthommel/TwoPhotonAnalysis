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


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#define MAX_PHOTONS 100
#define MAX_PHOTON_PAIRS 100

class PhotonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PhotonAnalyzer(const edm::ParameterSet&);
      ~PhotonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT< edm::View<pat::Photon> > fPhotonToken;
      //edm::EDGetTokenT< reco::BeamSpot > fBeamSpotToken;
      //edm::EDGetTokenT< edm::View<reco::Vertex> > fVertexToken;
      
      TTree* fTree;

      int aNumPhotons;
      double aPhotonPt[MAX_PHOTONS];
      double aPhotonEta[MAX_PHOTONS];
      double aPhotonPhi[MAX_PHOTONS];
      double aPhotonVtxX[MAX_PHOTONS];
      double aPhotonVtxY[MAX_PHOTONS];
      double aPhotonVtxZ[MAX_PHOTONS];

      int aNumPhotonPairs;
      double aPhotonPairMass[MAX_PHOTON_PAIRS];
      double aPhotonPairPt[MAX_PHOTON_PAIRS];
};

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonAnalyzer);
