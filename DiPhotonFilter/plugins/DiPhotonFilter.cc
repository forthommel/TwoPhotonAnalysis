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

#include "DiPhotonFilter.h"

//
// constructors and destructor
//
DiPhotonFilter::DiPhotonFilter(const edm::ParameterSet& iConfig) :
  fPhotonToken(consumes< edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"))),
  fVertexToken(consumes< edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexTag"))),
  fMinMpair(iConfig.getParameter<double>("minimalDiPhotonMass")),
  fMaxMpair(iConfig.getParameter<double>("maximalDiPhotonMass")),
  fMaxVertexDistance(iConfig.getParameter<double>("maximalVertexDistance"))
{
}


DiPhotonFilter::~DiPhotonFilter()
{
}


//
// member functions
//

bool
DiPhotonFilter::passSinglePhotonCriteria(const pat::Photon& photon) const
{
  return true;
}

bool
DiPhotonFilter::passVertexCriteria(const reco::Vertex& vertex, const pat::Photon& photon1, const pat::Photon& photon2) const
{
  const reco::Candidate::Point vp1 = photon1.vertex(), vp2 = photon2.vertex(), vv = vertex.position();
  
  // FIXME overly simplified!!
  if (std::sqrt((vp1-vv).mag2())>fMaxVertexDistance) return false;
  if (std::sqrt((vp2-vv).mag2())>fMaxVertexDistance) return false;

  return true;
}

bool
DiPhotonFilter::passDoublePhotonCriteria(const pat::Photon& photon1, const pat::Photon& photon2) const
{
  /*std::cout << "-----> Photons with pt=" << photon1.pt() << "///" << photon2.pt() << std::endl;
  std::cout << "  ---> vertex: " << photon1.vertex() << "///" << photon2.vertex() << std::endl;*/
  const reco::Candidate::LorentzVector photon_pair = photon1.p4()+photon2.p4();
  const double mpair = photon_pair.M();

  if (fMinMpair>0. and mpair<fMinMpair) return false;
  if (fMaxMpair>0. and mpair>fMaxMpair) return false;

  return true;
}

// ------------ method called on each new Event  ------------
bool
DiPhotonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByToken(fPhotonToken, photons);
  if (!photons.isValid()) {
    edm::LogError("DiPhotonFilter") << "No photon collection retrieved!";
    return false;
  }
  if (!photons->size()) {
    //edm::LogError("DiPhotonFilter") << "Photon collection is empty!";
    return false;
  }

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(fVertexToken, vertices);

  /*edm::Handle<VertexCandidateMap> vertexCandidateMap;
        evt.getByToken( vertexCandidateMapToken_, vertexCandidateMap );*/

  bool has_candidate = false;

  for (edm::View<pat::Photon>::const_iterator photon1=photons->begin(); photon1!=photons->end(); photon1++) {
    if (!passSinglePhotonCriteria(*photon1)) continue;
    for (edm::View<pat::Photon>::const_iterator photon2=photon1+1; photon2!=photons->end(); photon2++) {
      if (!passSinglePhotonCriteria(*photon2)) continue;

      /*edm::Ptr<reco::Vertex> pvx = fVertexSelector->select(photon1, photon2, vertices->ptrs(), *vertexCandidateMap, conversions->ptrs(), conversionsSingleLeg->ptrs(),
                                        vertexPoint, useSingleLeg_ );*/

      bool matched_vertex = false;
      for (edm::View<reco::Vertex>::const_iterator vtx=vertices->begin(); vtx!=vertices->end(); vtx++) {
        if (passVertexCriteria(*vtx, *photon1, *photon2)) matched_vertex = true;
      }
      if (!matched_vertex) continue;
      if (!passDoublePhotonCriteria(*photon1, *photon2)) continue;
      has_candidate = true; // if all criteria are passed, flag the event as a candidate

      if (has_candidate) {
        /*const reco::Candidate::LorentzVector photon_pair = photon1->p4()+photon2->p4();
        const double mpair = photon_pair.M(),
                     ptpair = photon_pair.Pt();*/
      }
    }
  }

  //if(!iEvent.isRealData()) { }

  return has_candidate;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DiPhotonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DiPhotonFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DiPhotonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DiPhotonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DiPhotonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DiPhotonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiPhotonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
