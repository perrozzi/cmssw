//  \class MuScleFitPlotter
//  Plotter for simulated,generated and reco info of muons
//
//  $Date: 2009/10/14 12:58:35 $
//  $Revision: 1.10 $
//  \author  C.Mariotti, S.Bolognesi - INFN Torino / T.Dorigo, M.De Mattia - INFN Padova
//
// ----------------------------------------------------------------------------------

#include "MuonAnalysis/MomentumScaleCalibration/interface/MuScleFitPlotter.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/Histograms.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/MuScleFitUtils.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include <CLHEP/Vector/LorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TMinuit.h"
#include <vector>

using namespace std;
using namespace edm;
using namespace reco; // For AODSIM MC objects

// Constructor
// ----------
MuScleFitPlotter::MuScleFitPlotter(string theGenInfoRootFileName){
  outputFile = new TFile(theGenInfoRootFileName.c_str(),"RECREATE");
  fillHistoMap();
}

MuScleFitPlotter::~MuScleFitPlotter(){
  outputFile->cd();
  writeHistoMap();
  outputFile->Close();
}

// Find and store in histograms the generated resonance and muons
// --------------------------------------------------------------
void MuScleFitPlotter::fillGen1(Handle<GenParticleCollection> genParticles)
{
  bool prova = false;
  //Loop on generated particles
  pair<reco::Particle::LorentzVector,reco::Particle::LorentzVector> muFromRes;
  reco::Particle::LorentzVector genRes;

  int mothersFound[] = {0, 0, 0, 0, 0, 0};

  for( GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter!=genParticles->end(); ++mcIter ) {
    int status = mcIter->status();
    int pdgId = abs(mcIter->pdgId());
    //Check if it's a resonance
    if( status == 2 &&
        ( pdgId==23  || pdgId==443    || pdgId==100443 ||
          pdgId==553 || pdgId==100553 || pdgId==200553 ) ) {
      genRes = mcIter->p4();
      if( pdgId == 23 ) mapHisto["hGenResZ"]->Fill(genRes);
      else if( pdgId == 443 ) mapHisto["hGenResJPsi"]->Fill(genRes);
      else if( pdgId == 553 ) mapHisto["hGenResUpsilon1S"]->Fill(genRes);
    }
    //Check if it's a muon from a resonance
    if( status==1 && pdgId==13 ) {
      int momPdgId = abs(mcIter->mother()->pdgId());
      if( momPdgId==23  || momPdgId==443    || momPdgId==100443 || 
          momPdgId==553 || momPdgId==100553 || momPdgId==200553 ) {
        if( momPdgId == 23 ) mothersFound[0] = 1;
        if( momPdgId == 443 ) mothersFound[5] = 1;
        if( momPdgId == 553 ) mothersFound[3] = 1;
	mapHisto["hGenMu"]->Fill(mcIter->p4());
	cout<<"genmu "<<mcIter->p4()<<endl;
	if(mcIter->charge()>0){
	  muFromRes.first = mcIter->p4();
	  // prova = true;
	}
	else muFromRes.second = mcIter->p4();
      }
    }
  }
  //   if(!prova)
  //     cout<<"hgenmumu not found"<<endl;

  if( mothersFound[0] == 1 ) {
    mapHisto["hGenMuMuZ"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuZ"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuZ"]->Fill( muFromRes.second,genRes, -1 );
  }
  if( mothersFound[3] == 1 ) {
    mapHisto["hGenMuMuUpsilon1S"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuUpsilon1S"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuUpsilon1S"]->Fill( muFromRes.second,genRes, -1 );
  }
  if( mothersFound[5] == 1 ) {
    mapHisto["hGenMuMuJPsi"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuJPsi"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuJPsi"]->Fill( muFromRes.second,genRes, -1 );
  }

  mapHisto["hGenResVsSelf"]->Fill( genRes, genRes, 1 );
}

// Find and store in histograms the generated resonance and muons
// --------------------------------------------------------------
void MuScleFitPlotter::fillGen2(Handle<HepMCProduct> evtMC)
{
  //Loop on generated particles
  const HepMC::GenEvent* Evt = evtMC->GetEvent();
  pair<reco::Particle::LorentzVector,reco::Particle::LorentzVector> muFromRes; 
  reco::Particle::LorentzVector genRes;

  int mothersFound[] = {0, 0, 0, 0, 0, 0};

  for (HepMC::GenEvent::particle_const_iterator part=Evt->particles_begin(); 
       part!=Evt->particles_end(); part++) {
    int status = (*part)->status();
    int pdgId = abs((*part)->pdg_id());
    //cout<<"PDG ID "<< (*part)->pdg_id() <<"    status "<< (*part)->status()
    //<<"   pt "<<(*part)->momentum().perp()<< "     eta  "<<(*part)->momentum().eta()<<endl    ;
     //Check if it's a resonance	
    if( status==2 && 
        ( pdgId==23  || pdgId==443    || pdgId==100443 ||
          pdgId==553 || pdgId==100553 || pdgId==200553 ) ) {
      genRes = reco::Particle::LorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
                                             (*part)->momentum().pz(),(*part)->momentum().e());
      if( pdgId == 23 ) mapHisto["hGenResZ"]->Fill(genRes);
      if( pdgId == 443 ) mapHisto["hGenResJPsi"]->Fill(genRes);
      if( pdgId == 553 ) {
        // cout << "genRes mass = " << CLHEP::HepLorentzVector(genRes.x(),genRes.y(),genRes.z(),genRes.t()).m() << endl;
        mapHisto["hGenResUpsilon1S"]->Fill(genRes);
      }
    }
    //Check if it's a muon from a resonance
    if (pdgId==13 && status==1) {      
      bool fromRes=false;
      for (HepMC::GenVertex::particle_iterator mother = 
	     (*part)->production_vertex()->particles_begin(HepMC::ancestors);
	   mother != (*part)->production_vertex()->particles_end(HepMC::ancestors); ++mother) {
        int motherPdgId = (*mother)->pdg_id();
	if (motherPdgId==23  || motherPdgId==443    || motherPdgId==100443 || 
	    motherPdgId==553 || motherPdgId==100553 || motherPdgId==200553) {
	  fromRes=true;
          if( motherPdgId == 23 ) mothersFound[0] = 1;
          if( motherPdgId == 443 ) mothersFound[3] = 1;
          if( motherPdgId == 553 ) mothersFound[5] = 1;
	}
      }

      if(fromRes) {	
	mapHisto["hGenMu"]->Fill(reco::Particle::LorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
							       (*part)->momentum().pz(),(*part)->momentum().e()));
	mapHisto["hGenMuVSEta"]->Fill(reco::Particle::LorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
								    (*part)->momentum().pz(),(*part)->momentum().e()));
	if((*part)->pdg_id()==-13)
	  muFromRes.first = (reco::Particle::LorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
							   (*part)->momentum().pz(),(*part)->momentum().e()));
	else
	  muFromRes.second = (reco::Particle::LorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
							    (*part)->momentum().pz(),(*part)->momentum().e()));
      }
    }
  }
  if( mothersFound[0] == 1 ) {
    mapHisto["hGenMuMuZ"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuZ"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuZ"]->Fill( muFromRes.second,genRes, -1 );
  }
  if( mothersFound[3] == 1 ) {
    mapHisto["hGenMuMuUpsilon1S"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuUpsilon1S"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuUpsilon1S"]->Fill( muFromRes.second,genRes, -1 );
  }
  if( mothersFound[5] == 1 ) {
    mapHisto["hGenMuMuJPsi"]->Fill(muFromRes.first+muFromRes.second);
    mapHisto["hGenResVSMuJPsi"]->Fill( muFromRes.first, genRes, 1 );
    mapHisto["hGenResVSMuJPsi"]->Fill( muFromRes.second,genRes, -1 );
  }
  mapHisto["hGenResVsSelf"]->Fill( genRes, genRes, 1 );
}

// Find and store in histograms the simulated resonance and muons
// --------------------------------------------------------------
 void MuScleFitPlotter::fillSim(Handle<SimTrackContainer> simTracks){

   vector<SimTrack> simMuons;

   //Loop on simulated tracks
   for (SimTrackContainer::const_iterator simTrack=simTracks->begin(); simTrack!=simTracks->end(); ++simTrack) {
     // Select the muons from all the simulated tracks
     if (fabs((*simTrack).type())==13) {
       simMuons.push_back(*simTrack);	  
       mapHisto["hSimMu"]->Fill((*simTrack).momentum());
     }
   }
   mapHisto["hSimMu"]->Fill(simMuons.size());

   // Recombine all the possible Z from simulated muons
   if (simMuons.size()>=2) {
     for (vector<SimTrack>::const_iterator  imu=simMuons.begin(); imu != simMuons.end(); ++imu) {   
       for (vector<SimTrack>::const_iterator imu2=imu+1; imu2!=simMuons.end(); ++imu2) {
	 if (imu==imu2) continue;
	    
	 // Try all the pairs with opposite charge
	 if (((*imu).charge()*(*imu2).charge())<0) {
	   reco::Particle::LorentzVector Z = (*imu).momentum()+(*imu2).momentum();
	   mapHisto["hSimMuPMuM"]->Fill(Z); 
	 }
       }
     }
   
     // Plots for the best possible simulated resonance
     pair<SimTrack,SimTrack> simMuFromBestRes = MuScleFitUtils::findBestSimuRes(simMuons);
     reco::Particle::LorentzVector bestSimZ = (simMuFromBestRes.first).momentum()+(simMuFromBestRes.second).momentum();
     mapHisto["hSimBestRes"]->Fill(bestSimZ);
     if (fabs(simMuFromBestRes.first.momentum().eta())<2.5 && fabs(simMuFromBestRes.second.momentum().eta())<2.5 &&
	 simMuFromBestRes.first.momentum().pt()>2.5 && simMuFromBestRes.second.momentum().pt()>2.5) {
       mapHisto["hSimBestResVSMu"]->Fill (simMuFromBestRes.first.momentum(), bestSimZ, int(simMuFromBestRes.first.charge()));
       mapHisto["hSimBestResVSMu"]->Fill (simMuFromBestRes.second.momentum(),bestSimZ, int(simMuFromBestRes.second.charge()));
    }
   }  
 }

// Find and store in histograms the RIGHT simulated resonance and muons
// --------------------------------------------------------------
 void MuScleFitPlotter::fillGenSim(Handle<HepMCProduct> evtMC, Handle<SimTrackContainer> simTracks){
   pair <reco::Particle::LorentzVector, reco::Particle::LorentzVector> simMuFromRes = 
     MuScleFitUtils::findSimMuFromRes(evtMC,simTracks);
   //Fill resonance info
   reco::Particle::LorentzVector rightSimRes = (simMuFromRes.first)+(simMuFromRes.second);
   mapHisto["hSimRightRes"]->Fill(rightSimRes);
   /*if ((fabs(simMuFromRes.first.Eta())<2.5 && fabs(simMuFromRes.second.Eta())<2.5) 
       && simMuFromRes.first.Pt()>2.5 && simMuFromRes.second.Pt()>2.5) {
       }*/
 }


// Find and store in histograms the reconstructed resonance and muons
// --------------------------------------------------------------
 void MuScleFitPlotter::fillRec(vector<reco::LeafCandidate>& muons){
   for(vector<reco::LeafCandidate>::const_iterator mu1 = muons.begin(); mu1!=muons.end(); mu1++){
     mapHisto["hRecMu"]->Fill(mu1->p4());
     mapHisto["hRecMuVSEta"]->Fill(mu1->p4());
     for(vector<reco::LeafCandidate>::const_iterator mu2 = muons.begin(); mu2!=muons.end(); mu2++){  
       if (mu1==mu2) continue;
       reco::Particle::LorentzVector Res (mu1->p4()+mu2->p4());
        mapHisto["hRecMuPMuM"]->Fill(Res);	  
     } 
   }
   mapHisto["hRecMu"]->Fill(muons.size());
 }


// Histogram booking
// -----------------
void MuScleFitPlotter::fillHistoMap() {

  // Generated Z and muons
  // ---------------------
  mapHisto["hGenResJPsi"]      = new HParticle   ("hGenResJPsi", 3.09685, 3.09695);
  mapHisto["hGenResUpsilon1S"] = new HParticle   ("hGenResUpsilon1S", 9., 11.);
  mapHisto["hGenResZ"]         = new HParticle   ("hGenResZ", 60., 120.);
  mapHisto["hGenMu"]      = new HParticle  ("hGenMu");
  mapHisto["hGenMuVSEta"] = new HPartVSEta ("hGenMuVSEta");

  mapHisto["hGenMuMuJPsi"]      = new HParticle   ("hGenMuMuJPsi",3.09685, 3.09695 );
  mapHisto["hGenResVSMuJPsi"]   = new HMassVSPart ("hGenResVSMuJPsi",3.09685, 3.09695);
  mapHisto["hGenMuMuUpsilon1S"]      = new HParticle   ("hGenMuMuUpsilon1S", 9., 11.);
  mapHisto["hGenResVSMuUpsilon1S"]   = new HMassVSPart ("hGenResVSMuUpsilon1S", 9., 11.);
  mapHisto["hGenMuMuZ"]      = new HParticle   ("hGenMuMuZ", 60., 120.);
  mapHisto["hGenResVSMuZ"]   = new HMassVSPart ("hGenResVSMuZ", 60., 120.);

  mapHisto["hGenResVsSelf"] = new HMassVSPart ("hGenResVsSelf");

  // Simulated resonance and muons
  // -----------------------------
  mapHisto["hSimMu"]      = new HParticle ("hSimMu");

  mapHisto["hSimMuPMuM"]      = new HParticle ("hSimMuPMuM");      
                                                                 
  mapHisto["hSimBestMu"]      = new HParticle ("hSimBestMu");
  mapHisto["hSimBestRes"]         = new HParticle  ("hSimBestRes");
  mapHisto["hSimBestResVSMu"]  = new HMassVSPart ("hSimBestResVSMu");
    
  mapHisto["hSimRightRes"]         = new HParticle  ("hSimRightZ");
 
  // Reconstructed resonance and muons
  // -----------------------------  
  mapHisto["hRecMu"]      = new HParticle ("hRecMu");
  mapHisto["hRecMuVSEta"]      = new HPartVSEta ("hRecMuVSEta");
  mapHisto["hRecMuPMuM"]         = new HParticle  ("hRecMuPMuM");
}  


// Histogram saving
// -----------------
void MuScleFitPlotter::writeHistoMap() {
  outputFile->cd();
  for (map<string, Histograms*>::const_iterator histo=mapHisto.begin(); 
       histo!=mapHisto.end(); histo++) {
    (*histo).second->Write();
  }
}

