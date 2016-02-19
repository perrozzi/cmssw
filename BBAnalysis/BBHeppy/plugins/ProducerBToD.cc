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


//  P4 WITH PION MASS HYPOTESIS
reco::Candidate::LorentzVector vtxP4(const reco::VertexCompositePtrCandidate & vtx) {
  reco::Candidate::LorentzVector sum;
  const std::vector<reco::CandidatePtr> & tracks = vtx.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;
    vec.SetPx((*track)->px());
    vec.SetPy((*track)->py());
    vec.SetPz((*track)->pz());
    vec.SetM(reco::ParticleMasses::piPlus );
    sum += vec;
  }
  return sum;
}


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
         
         bool itIsMerged;
       };

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
  maxvecSumIMCUTForUnique = 5;
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
      secondaryAUX.itIsMerged = false;
      SecondaryVerticesProxy.push_back(secondaryAUX);
    }
  }

//  int numberOf_SecondaryVerticesProxy_BeforMerging = SecondaryVerticesProxy.size();
//  int mergedVertices=0;

  // loop forward over all vertices and
  // check all vertices against each other for B->D chain
  int numberOfSteps = 0;
  unsigned int numberOf_SecondaryVerticesProxy_BeforMerging_lastStep;
  do {
    numberOf_SecondaryVerticesProxy_BeforMerging_lastStep = SecondaryVerticesProxy.size();
    for(unsigned int kVtx=SecondaryVerticesProxy.size(); kVtx>0 && SecondaryVerticesProxy.size()>1; --kVtx){

      // remove D vertices from the collection and add the tracks to the original one
      resolveBtoDchain(SecondaryVerticesProxy, kVtx-1);
    }

  numberOfSteps++;
  } while (numberOf_SecondaryVerticesProxy_BeforMerging_lastStep != SecondaryVerticesProxy.size());

//  if(numberOfSteps>1)
//    cout << "numberOfSteps:   " << numberOfSteps << endl;



  for(std::vector<VertexProxy>::const_iterator itSV = SecondaryVerticesProxy.begin(); itSV!=SecondaryVerticesProxy.end(); itSV++) {
    //if(isSelected(*itSV)) 
      output->push_back(itSV->vert);
      ifMerged->push_back(itSV->itIsMerged);
//      if(itSV->itIsMerged)
//        mergedVertices++;
  }

//int size = output->size();
////if(size + mergedVertices - numberOf_SecondaryVerticesProxy_BeforMerging!=0)
//if( mergedVertices !=0)
//  cout << " TRE SV DI FILA:    ci sono    " << output->size() << "   "<< mergedVertices << "   " << numberOf_SecondaryVerticesProxy_BeforMerging << 
//          " la differenza è:   "  << size + mergedVertices - numberOf_SecondaryVerticesProxy_BeforMerging << "   fine" << endl;

   iEvent.put(output );
   iEvent.put(ifMerged, "ifMerged");

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}



// MERGING THE SVs
void ProducerBToD::resolveBtoDchain(std::vector<VertexProxy> & coll,  unsigned int k){

  TVector3 ppv(pv.position().x(),pv.position().y(),pv.position().z());
  TVector3 SecondaryVertexPosition(coll[k].vert.position().x(), coll[k].vert.position().y(), coll[k].vert.position().z());
  TVector3 pvToNear= SecondaryVertexPosition -ppv;
  TVector3 pvToFar ;

  reco::Candidate::LorentzVector p4Near = vtxP4(coll[k].vert);
  GlobalVector momentumNear(p4Near.X(), p4Near.Y(), p4Near.Z());
  SecondaryVertex sNear(pv, coll[k].vert, momentumNear, true);
  reco::Candidate::LorentzVector p4Far;

  bool found = false;
  unsigned int i = 0;
  double ptRelMin =100000;

  for(unsigned int I = 0; I < coll.size(); I++) {
    if(I!=k) {
      TVector3 SecondaryVertexPositionFar(coll[I].vert.position().x(), coll[I].vert.position().y(), coll[I].vert.position().z());
      pvToFar = SecondaryVertexPositionFar -ppv;
      if((pvToFar.Mag() > pvToNear.Mag()) && (deltaR(pvToFar, pvToNear) < maxDRForUnique)) {
        TVector3 nearToFar = pvToFar - pvToNear;
        p4Far = vtxP4(coll[I].vert);
        TVector3 momentumFar(p4Far.X(), p4Far.Y(), p4Far.Z());
        double cosPA =  nearToFar.Dot(momentumFar) / momentumFar.Mag()/ nearToFar.Mag();
        double cosa  =  pvToNear. Dot(momentumFar) / pvToNear.Mag()   / momentumFar.Mag();
        double ptRel = sqrt(1.0 - cosa*cosa)* momentumFar.Mag();

        // Qui stanno tutte le condizioni per unire i due vertici, tranne una. La massa invariante sta dopo
        if((cosPA > minCosPAtomerge) && (ptRel < maxPtreltomerge) && (ptRel < ptRelMin) && (p4Far.mass() + p4Near.mass() < maxvecSumIMCUTForUnique )) {
          i=I;
          found=true;
          ptRelMin = ptRel;
        }
      }
    }
  }



  if(found) {

    GlobalVector momentumFar(p4Far.X(), p4Far.Y(), p4Far.Z());
    SecondaryVertex sFar(pv, coll[i].vert, momentumFar, true);

  // create a set of all tracks from both vertices, avoid double counting by using a std::set<>
    std::set<reco::CandidatePtr> trackrefs;
  // first vertex
    for(size_t j=0; j < sNear.numberOfSourceCandidatePtrs(); ++j)
      trackrefs.insert(sNear.daughterPtr(j));
  // second vertex
    for(size_t j=0; j < sFar.numberOfSourceCandidatePtrs(); ++j)
      trackrefs.insert(sFar.daughterPtr(j));

  // now calculate one LorentzVector from the track momenta
    reco::Candidate::LorentzVector mother;
    for(std::set<reco::CandidatePtr>::const_iterator it = trackrefs.begin(); it!= trackrefs.end(); ++it){
      reco::Candidate::LorentzVector temp ( (*it)->px(),(*it)->py(),(*it)->pz(), reco::ParticleMasses::piPlus );
      mother += temp;
    }

    if(mother.mass() < maxTOTALmassForUnique) {   //questa è l'ultima condizione

      const std::vector<reco::CandidatePtr> & tracks1 = sNear.daughterPtrVector();
      const std::vector<reco::CandidatePtr> & tracks2 = sFar.daughterPtrVector();
      for(std::vector<reco::CandidatePtr>::const_iterator ti = tracks2.begin(); ti!=tracks2.end(); ++ti) {
        std::vector<reco::CandidatePtr>::const_iterator it = find(tracks1.begin(), tracks1.end(), *ti);
        if (it==tracks1.end()) {
          coll[k].vert.addDaughter( *ti );
          coll[k].vert.setP4( (*ti)->p4() + coll[k].vert.p4() );
        }
      }
    coll[k].itIsMerged = true;
    coll.erase( coll.begin() + i  );
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

//  reco::Candidate::LorentzVector momentum = vtxP4(secVert);
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
