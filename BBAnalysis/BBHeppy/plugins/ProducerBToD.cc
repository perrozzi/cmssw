// -*- C++ -*-
//
// Package:    cartellaProducer/ProducerBToD
// Class:      ProducerBToD
// 
/**\class ProducerBToD ProducerBToD.cc cartellaProducer/ProducerBToD/plugins/ProducerBToD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giulio MANDORLI
//         Created:  Fri, 05 Feb 2016 15:36:04 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// includes from clusteringAnalyzerHeader
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <DataFormats/Candidate/interface/CompositeRefCandidateT.h>
#include "DataFormats/Common/interface/RefVectorIterator.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/Common/interface/Ptr.h>
#include <DataFormats/Common/interface/AssociativeIterator.h>
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/HandleBase.h"


// tutti gli include che trovo
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector2.h"
#include <TVector3.h>

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"
#include <RecoVertex/VertexPrimitives/interface/BasicVertexState.h>

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include <stack>

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/BTauReco/interface/ParticleMasses.h"




//#include "FWCore/ParameterSet/interface/InputTag.h"
//
// class declaration
//

class ProducerBToD : public edm::stream::EDProducer<> {
   public:
      explicit ProducerBToD(const edm::ParameterSet&);
      ~ProducerBToD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;



       struct VertexProxy{
         reco::VertexCompositePtrCandidate vert;
         double significance3D;
         bool itIsMerged;
       };


       // comparison operator for VertexProxy, used in sorting
       friend bool operator<(VertexProxy v1, VertexProxy v2){
         return (v1.significance3D<v2.significance3D);
       }

      bool PassAFisrtSelection(reco::VertexCompositePtrCandidate secVert);
      bool isSelected(VertexProxy secVertProxy);
      double GetSignificance(reco::VertexCompositePtrCandidate secVert);
      double GetSignificanceXY(reco::VertexCompositePtrCandidate secVert);
      double deltaR(TVector3 v1, TVector3 v2);

      void resolveBtoDchain(std::vector<VertexProxy> & coll, unsigned int k);


      typedef reco::TemplatedSecondaryVertex<reco::VertexCompositePtrCandidate> SecondaryVertex;




      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  //Where to get the informations

      edm::InputTag m_slimmedSecondaryVertices;
      edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > slimmedSecondaryVerticesToken;

      edm::InputTag m_offlineSlimmedPrimaryVertices;
      edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken;



// variables
//      const reco::Vertex & pv;
      reco::Vertex   pv;


      double maxDRForUnique;
      double maxPtreltomerge;
      double minCosPAtomerge;
      double maxvecSumIMCUTForUnique;
      double maxTOTALmassForUnique;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
ProducerBToD::ProducerBToD(const edm::ParameterSet& iConfig)
{

  //get the informations 
  slimmedSecondaryVerticesToken = consumes<std::vector<reco::VertexCompositePtrCandidate> >(edm::InputTag("slimmedSecondaryVertices"));
  offlineSlimmedPrimaryVerticesToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));


   produces<std::vector<reco::VertexCompositePtrCandidate> >(); 
   produces<std::vector<bool> >("ifMerged"); 

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


ProducerBToD::~ProducerBToD()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProducerBToD::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

  maxDRForUnique = 0.3;
  maxPtreltomerge = 6;
  minCosPAtomerge = 0.8;    // 36 deg
//  maxvecSumIMCUTForUnique = 4.5;
  maxTOTALmassForUnique = 6.5;

//get the informations in slimmedSecondaryVertices
  Handle<std::vector<reco::VertexCompositePtrCandidate> > SecondaryVerticesCollection;
  iEvent.getByToken(slimmedSecondaryVerticesToken, SecondaryVerticesCollection);
  std::vector<reco::VertexCompositePtrCandidate>  SecondaryVertices = *SecondaryVerticesCollection.product();

//get the informations in Primary Vertex
  Handle<std::vector<reco::Vertex> > vertices; 
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken, vertices);
  pv = (*vertices)[0];



    std::auto_ptr<std::vector<reco::VertexCompositePtrCandidate> >  output(new std::vector<reco::VertexCompositePtrCandidate>);
    std::auto_ptr<std::vector<bool> >  ifMerged(new std::vector<bool>);

// make the first selected vertices
  std::vector<VertexProxy>  SecondaryVerticesProxy;
  for(std::vector<reco::VertexCompositePtrCandidate>::const_iterator itSV = SecondaryVertices.begin(); itSV!=SecondaryVertices.end(); itSV++) {
    if(PassAFisrtSelection(*itSV)) {
      VertexProxy secondaryAUX;
      secondaryAUX.vert = *itSV;
      secondaryAUX.significance3D = GetSignificance(*itSV);
      secondaryAUX.itIsMerged = false;
      SecondaryVerticesProxy.push_back(secondaryAUX);
    }
  }

  sort( SecondaryVerticesProxy.begin(), SecondaryVerticesProxy.end()); // SecondaryVerticesProxy[0] is the one with less significance


//  bool b1= false;
    for(unsigned int kVtx=SecondaryVerticesProxy.size(); kVtx>0 && SecondaryVerticesProxy.size()>1; --kVtx){
//      int tempSize = SecondaryVerticesProxy.size();
      // remove D vertices from the collection and add the tracks to the original one
      resolveBtoDchain(SecondaryVerticesProxy, kVtx-1);

      // check how many times it merges tre SV all together
//      if(tempSize - SecondaryVerticesProxy.size() > 0) {
//        if(b1)
//          cout << "tre SV di fila" << endl;
//        else b1 = true;
//      }
//      else b1 = false;
    }

  // print the significance of the vertices in order
//  for (unsigned int i=0; i < SecondaryVerticesProxy.size(); i++)
//    cout << SecondaryVerticesProxy[i].significance3D << "\t";
//  cout << endl;


//  if(numberOfSteps>1)
//    cout << "numberOfSteps:   " << numberOfSteps << endl;



  for(std::vector<VertexProxy>::const_iterator itSV = SecondaryVerticesProxy.begin(); itSV!=SecondaryVerticesProxy.end(); itSV++) {
    //if(isSelected(*itSV)) 
      output->push_back(itSV->vert);
      ifMerged->push_back(itSV->itIsMerged);
  }


   iEvent.put(output );
   iEvent.put(ifMerged, "ifMerged");


}



// MERGING THE SVs
void ProducerBToD::resolveBtoDchain(std::vector<VertexProxy> & coll,  unsigned int k){

  TVector3 ppv(pv.position().x(),pv.position().y(),pv.position().z());
  TVector3 SecondaryVertexPosition1(coll[k].vert.position().x(), coll[k].vert.position().y(), coll[k].vert.position().z());
  TVector3 pvToNear;
  TVector3 pvToFar ;

  reco::Candidate::LorentzVector p4Near;
  reco::Candidate::LorentzVector p4Far;


  bool found = false;
  unsigned int nearIdxTemp = 0;
  unsigned int farIdxTemp = 0;
  unsigned int nearIdx = 0;
  unsigned int farIdx = 0;
  double ptRelMin =100000;

  for(unsigned int I = 0; I < coll.size(); I++) {
    if(I!=k) {
      TVector3 SecondaryVertexPosition2(coll[I].vert.position().x(), coll[I].vert.position().y(), coll[I].vert.position().z());
      pvToFar = SecondaryVertexPosition2 -ppv;
      pvToNear= SecondaryVertexPosition1 -ppv;
      if(deltaR(pvToFar, pvToNear) < maxDRForUnique) {
        TVector3 nearToFar = pvToFar - pvToNear;
        //swap if k is the farther
        if(pvToNear.Mag() > pvToFar.Mag()) {
          nearIdxTemp = k;
          farIdxTemp = I;
          nearToFar = pvToFar - pvToNear;
        }
        else {
          nearIdxTemp = I;
          farIdxTemp = k;
          nearToFar = pvToNear - pvToFar;
          pvToNear = pvToFar;
        }

        p4Near = coll[nearIdxTemp].vert.p4();
        p4Far = coll[farIdxTemp].vert.p4();
        TVector3 momentumFar(p4Far.X(), p4Far.Y(), p4Far.Z());
        double cosPA =  nearToFar.Dot(momentumFar) / momentumFar.Mag()/ nearToFar.Mag();
        double cosa  =  pvToNear. Dot(momentumFar) / pvToNear.Mag()   / momentumFar.Mag();
        double ptRel = sqrt(1.0 - cosa*cosa)* momentumFar.Mag();

        // Qui stanno tutte le condizioni per unire i due vertici, tranne una. La massa invariante sta dopo
        if((cosPA > minCosPAtomerge) && (ptRel < maxPtreltomerge) && (ptRel < ptRelMin) && (p4Far.mass()<2.0) && (p4Near.mass() > 1.2 )) {
          farIdx=farIdxTemp;
          nearIdx=nearIdxTemp;
          found=true;
          ptRelMin = ptRel;
        }
      }
    }
  }


  unsigned int index_moreSignVertex = nearIdx;
  unsigned int index_lessSignVertex = farIdx;
  if(coll[nearIdx].significance3D < coll[farIdx].significance3D) {
    index_moreSignVertex = farIdx;
    index_lessSignVertex = nearIdx;
  }

  if(found) {

    reco::Candidate::LorentzVector p4ToThrow = coll[index_lessSignVertex].vert.p4();
    GlobalVector momentumToThrow(p4ToThrow.X(), p4ToThrow.Y(), p4ToThrow.Z());
    SecondaryVertex sToThrow(pv, coll[index_lessSignVertex].vert, momentumToThrow, true);

    reco::Candidate::LorentzVector p4ToKeep = coll[index_moreSignVertex].vert.p4();
    GlobalVector momentumToKeep(p4ToKeep.X(), p4ToKeep.Y(), p4ToKeep.Z());
    SecondaryVertex sToKeep(pv, coll[index_moreSignVertex].vert, momentumToKeep, true);

  // create a set of all tracks from both vertices, avoid double counting by using a std::set<>
    std::set<reco::CandidatePtr> trackrefs;
  // first vertex
    for(size_t j=0; j < sToThrow.numberOfSourceCandidatePtrs(); ++j)
      trackrefs.insert(sToThrow.daughterPtr(j));
  // second vertex
    for(size_t j=0; j < sToKeep.numberOfSourceCandidatePtrs(); ++j)
      trackrefs.insert(sToKeep.daughterPtr(j));

  // now calculate one LorentzVector from the track momenta
    reco::Candidate::LorentzVector mother;
    for(std::set<reco::CandidatePtr>::const_iterator it = trackrefs.begin(); it!= trackrefs.end(); ++it)
      mother += (*it)->p4();


    if(mother.mass() < maxTOTALmassForUnique) {   //questa è l'ultima condizione

      const std::vector<reco::CandidatePtr> & tracks1 = sToKeep.daughterPtrVector();
      const std::vector<reco::CandidatePtr> & tracks2 = sToThrow.daughterPtrVector();
      for(std::vector<reco::CandidatePtr>::const_iterator ti = tracks2.begin(); ti!=tracks2.end(); ++ti) {
        std::vector<reco::CandidatePtr>::const_iterator it = find(tracks1.begin(), tracks1.end(), *ti);
        if (it==tracks1.end()) {
          coll[index_moreSignVertex].vert.addDaughter( *ti );
          coll[index_moreSignVertex].vert.setP4( (*ti)->p4() + coll[index_moreSignVertex].vert.p4() );
        }
      }
    coll[index_moreSignVertex].itIsMerged = true;
    coll.erase( coll.begin() + index_lessSignVertex  );
    }
  }

}  //parentesi di fine resolveBtoDchain

//SECONDARY VERTEX SELECTION
//bool ProducerBToD::isSelected(reco::VertexCompositePtrCandidate secVert) {

//  double significance = GetSignificance(secVert);
//  double significanceXY = GetSignificanceXY(secVert);

//  TVector3 ppv(pv.position().x(),pv.position().y(),pv.position().z());
//  TVector3 SecondaryVertexPosition(secVert.position().x(), secVert.position().y(), secVert.position().z());
//  TVector3 VertexDirection = SecondaryVertexPosition -ppv;

//  reco::Candidate::LorentzVector momentum = secVert.p4();
//  TVector3 pDirection(momentum.X(), momentum.Y(), momentum.Z());


//  double distances_SV_p = pDirection.Dot(VertexDirection); 
//  distances_SV_p = distances_SV_p/pDirection.Mag()/VertexDirection.Mag(); 

//  if((momentum.mass()>1.5) && (momentum.mass()<6.5) && (secVert.numberOfDaughters()>2) && (abs(VertexDirection.PseudoRapidity())<2.5) && (momentum.pt()>8.) && (abs(VertexDirection.Perp())<2.) && (significanceXY > 3.) && (significance > 5.) && (distances_SV_p>0.95))
//    return true;
//  else
//    return false;
//}

//SECONDARY VERTEX FIRST SELECTION
bool ProducerBToD::PassAFisrtSelection(reco::VertexCompositePtrCandidate secVert) {

//  double significance = GetSignificance(secVert);
//  double significanceXY = GetSignificanceXY(secVert);

//  TVector3 ppv(pv.position().x(),pv.position().y(),pv.position().z());
//  TVector3 SecondaryVertexPosition(secVert.position().x(), secVert.position().y(), secVert.position().z());
//  TVector3 VertexDirection = SecondaryVertexPosition -ppv;

//  TVector3 pDirection(secVert.p4().Px(), secVert.p4().Py(), secVert.p4().Pz());

//  double distances_SV_p = pDirection.Dot(VertexDirection); 
//  distances_SV_p = distances_SV_p/pDirection.Mag()/VertexDirection.Mag(); 

//  if((secVert.mass()>1.) && (secVert.mass()<10.) && (abs(VertexDirection.PseudoRapidity())<3.) && (abs(VertexDirection.Perp())<2.) && (significanceXY > 2.) && (significance > 3.) && (distances_SV_p>0.80))
//    return true;
//  else
//    return false;

  return true;
}
//SIGNIFICANCES
double ProducerBToD::GetSignificanceXY(reco::VertexCompositePtrCandidate secVert) {

  VertexState s2(RecoVertex::convertPos(secVert.position()),RecoVertex::convertError(secVert.error()));
  VertexState s1(RecoVertex::convertPos(pv.position()),RecoVertex::convertError(pv.error()));

  VertexDistanceXY dist;
  return dist.distance(s1,s2).significance();

}

double ProducerBToD::GetSignificance(reco::VertexCompositePtrCandidate secVert) {

  VertexState s2(RecoVertex::convertPos(secVert.position()),RecoVertex::convertError(secVert.error()));
  VertexState s1(RecoVertex::convertPos(pv.position()),RecoVertex::convertError(pv.error()));

  VertexDistance3D dist;
  return dist.distance(s1,s2).significance();

}

double ProducerBToD::deltaR(TVector3 v1, TVector3 v2) {

  double Phi = TVector2::Phi_mpi_pi(v1.Phi() - v2.Phi());
  double Eta = v1.PseudoRapidity() - v2.PseudoRapidity();
  double R = std::sqrt(Phi*Phi+Eta*Eta);

  return R;
}




// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ProducerBToD::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ProducerBToD::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ProducerBToD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ProducerBToD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ProducerBToD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ProducerBToD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProducerBToD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerBToD);
