// -*- C++ -*-
//
// Package:    TwoPhoton/DiPhotonFilter
// Class:      DiPhotonFilter
// 
/**\class DiPhotonFilter DiPhotonFilter.cc TwoPhoton/DiPhotonFilter/plugins/DiPhotonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Sun, 22 May 2016 15:30:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//#include "VertexSelectorBase.h"

//
// class declaration
//

class DiPhotonFilter : public edm::stream::EDFilter<> {
   public:
      explicit DiPhotonFilter(const edm::ParameterSet&);
      ~DiPhotonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      bool passSinglePhotonCriteria(const pat::Photon& photon) const;
      bool passVertexCriteria(const reco::Vertex& vertex, const pat::Photon& photon1, const pat::Photon& photon2) const;
      bool passDoublePhotonCriteria(const pat::Photon& photon1, const pat::Photon& photon2) const;
      
      edm::EDGetTokenT< edm::View<pat::Photon> > fPhotonToken;
      edm::EDGetTokenT< edm::View<reco::Vertex> > fVertexToken;

      double fMinMpair, fMaxMpair;
      double fMaxVertexDistance;

      //std::unique_ptr<flashgg::VertexSelectorBase> fVertexSelector;


};

//define this as a plug-in
DEFINE_FWK_MODULE(DiPhotonFilter);
