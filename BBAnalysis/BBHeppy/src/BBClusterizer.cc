//#include "PhysicsTools/Heppy/interface/ReclusterJets.h"
#include "BBAnalysis/BBHeppy/interface/BBClusterizer.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "fastjet/tools/Pruner.hh"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
using namespace std;

//using namespace std;
using namespace fastjet;

namespace heppy{

BBClusterizer::BBClusterizer(const std::vector<LorentzVector> & candidates, const LorentzVector & b1, const LorentzVector & b2, double ktpower, double rparam) :
   b1_(b1),b2_(b2), ktpower_(ktpower), rparam_(rparam)
{
  //  preparing fastjets
      fjInputs.clear();
      unsigned int index=0;
      for (index=0; index<candidates.size(); index++) {
        fastjet::PseudoJet j(candidates[index].px(), candidates[index].py(), candidates[index].pz(), candidates[index].energy());
        j.set_user_index(index);  // in case we want to know which piece ended where
        fjInputs.push_back(j);
      }

      fastjet::PseudoJet j1_(b1.px(), b1.py(), b1.pz(), b1.energy());
      j1_.set_user_index(-1);
      fjInputs.push_back(j1_);
      fastjet::PseudoJet j2_(b2.px(), b2.py(), b2.pz(), b2.energy());
      j2_.set_user_index(-2);
      fjInputs.push_back(j2_);
 

  // choose a jet definition
  fastjet::JetDefinition jet_def;

  // prepare jet def 
  if (ktpower_ == 1.0) {
    jet_def = JetDefinition(kt_algorithm, rparam_);
  }  else if (ktpower_ == 0.0) {
    jet_def = JetDefinition(cambridge_algorithm, rparam_);
  }  else if (ktpower_ == -1.0) {
    jet_def = JetDefinition(antikt_algorithm, rparam_);
  }  else {
    throw cms::Exception("InvalidArgument", "Unsupported ktpower value");
  }
  
  // print out some infos
  //  cout << "Clustering with " << jet_def.description() << endl;
  ///
  // define jet clustering sequence
  fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, jet_def)); 
}


// ast jet with 4B
BBClusterizer::BBClusterizer(const std::vector<LorentzVector> & candidates, const LorentzVector & b1, const LorentzVector & b2, const LorentzVector & b3, const LorentzVector & b4, double ktpower, double rparam) :
   b1_(b1),b2_(b2), ktpower_(ktpower), rparam_(rparam)
{
  //  preparing fastjets
      fjInputs.clear();
      unsigned int index=0;
      for (index=0; index<candidates.size(); index++) {
        fastjet::PseudoJet j(candidates[index].px(), candidates[index].py(), candidates[index].pz(), candidates[index].energy());
        j.set_user_index(index);  // in case we want to know which piece ended where
        fjInputs.push_back(j);
      }

      fastjet::PseudoJet j1_(b1.px(), b1.py(), b1.pz(), b1.energy());
      j1_.set_user_index(-1);
      fjInputs.push_back(j1_);
      fastjet::PseudoJet j2_(b2.px(), b2.py(), b2.pz(), b2.energy());
      j2_.set_user_index(-2);
      fjInputs.push_back(j2_);
      fastjet::PseudoJet j3_(b3.px(), b3.py(), b3.pz(), b3.energy());
      j3_.set_user_index(-3);
      fjInputs.push_back(j3_);
      fastjet::PseudoJet j4_(b4.px(), b4.py(), b4.pz(), b4.energy());
      j4_.set_user_index(-4);
      fjInputs.push_back(j4_);


  // choose a jet definition
  fastjet::JetDefinition jet_def;

  // prepare jet def 
  if (ktpower_ == 1.0) {
    jet_def = JetDefinition(kt_algorithm, rparam_);
  }  else if (ktpower_ == 0.0) {
    jet_def = JetDefinition(cambridge_algorithm, rparam_);
  }  else if (ktpower_ == -1.0) {
    jet_def = JetDefinition(antikt_algorithm, rparam_);
  }  else {
    throw cms::Exception("InvalidArgument", "Unsupported ktpower value");
  }
  
  // print out some infos
  //  cout << "Clustering with " << jet_def.description() << endl;
  ///
  // define jet clustering sequence
  fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, jet_def)); 
}



std::vector<math::XYZTLorentzVector> BBClusterizer::makeP4s(const std::vector<fastjet::PseudoJet> &jets) {
  std::vector<math::XYZTLorentzVector> JetObjectsAll;
  for (const fastjet::PseudoJet & pj : jets) {
/*    std::vector<fastjet::PseudoJet> constituents = pj.constituents();
    std::cout << "Constituents for " << pj.pt() << " " ;
    for (unsigned j = 0; j < constituents.size(); j++) {
		std::cout << constituents[j].user_index() << " ";
    }
    std::cout << std::endl;*/
    JetObjectsAll.push_back( LorentzVector( pj.px(), pj.py(), pj.pz(), pj.e() ) );
  }
  return JetObjectsAll;
}



std::vector<math::XYZTLorentzVector> BBClusterizer::getBJets() {
  // recluster jet
  vector<PseudoJet> jets = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets());
  vector<PseudoJet> JetsWithB;


  for(vector<PseudoJet>::const_iterator it = jets.begin(); it != jets.end(); it++) 
    if(SVIsInTheJet(*it) > 0) JetsWithB.push_back(*it);


  PseudoJet parent1, parent2;
  bool had_parents = JetsWithB[0].has_parents(parent1,parent2);

//  std::cout << "Prima di while:  SV in 1   " << SVIsInTheJet(parent1) << "    SV in 2   " << SVIsInTheJet(parent2) << std::endl;

  while ((JetsWithB.size()==1) && (had_parents)) {
      int parent1HasSV = SVIsInTheJet(parent1);
      int parent2HasSV = SVIsInTheJet(parent2);

    if ((parent1HasSV == 1) && (parent2HasSV == 1)) {
      JetsWithB[0] = parent1;
      JetsWithB.push_back(parent2);
    }
    if ((parent1HasSV == 2) && (parent2HasSV < 2)) {
      had_parents = parent1.has_parents(parent1,parent2);
    }
    if ((parent1HasSV < 2) && (parent2HasSV == 2)) {
      had_parents = parent2.has_parents(parent1,parent2);
    }

  }

  bool toSwap=false;
  std::vector<fastjet::PseudoJet> constituents1 = JetsWithB[0].constituents();
  for (unsigned j = 0; j < constituents1.size(); j++)
        if(constituents1[j].user_index()==-2)
              toSwap = true;
 
  if(toSwap) {
        auto tmp = JetsWithB[0];
        JetsWithB[0]=JetsWithB[1];
        JetsWithB[1]=tmp;
  } 
  // return
  return makeP4s(JetsWithB);
}


//get jet with 4B


std::vector<math::XYZTLorentzVector> BBClusterizer::getBJets4B() {
  // recluster jet
  vector<PseudoJet> jets = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets());
  vector<PseudoJet> JetsWithB;

  for(vector<PseudoJet>::const_iterator it = jets.begin(); it != jets.end(); it++) 
    if(SVIsInTheJet(*it) > 0) JetsWithB.push_back(*it);

  bool beginAgainTheLoop = false;
  for(unsigned int n = 0; n < JetsWithB.size(); n++) {
    if (beginAgainTheLoop) n = 0;
    beginAgainTheLoop = false;
    if (SVIsInTheJet(JetsWithB[n]) > 1) {
        vector<PseudoJet> twoBJetVector = DivideTheJet(JetsWithB[n]);
        JetsWithB[n] = twoBJetVector[0];
        JetsWithB.push_back(twoBJetVector[1]);
        beginAgainTheLoop = true;
//        std::cout << "Siamo qui" << std::endl;
    }
  }

  if (JetsWithB.size() != 4) 
    std::cout << "Error: not 4 jets" << std::endl;

  vector<PseudoJet> JetsWithOrderedB;
  vector<int> order;

  for(unsigned int n = 0; n < JetsWithB.size(); n++) {
    std::vector<fastjet::PseudoJet> constituents = JetsWithB[n].constituents();
//        std::cout << "jet indexes:  " << std::endl;
        for (unsigned j = 0; j < constituents.size(); j++) {
//            std::cout << constituents[j].user_index() << "  \t ";
            if(constituents[j].user_index() < 0)
                order.push_back((-1)*constituents[j].user_index());
        }
//        std::cout << std::endl;
  }

  for(unsigned int n = 0; n < JetsWithB.size(); n++) 
    for(int j = 0; j < 4; j++) 
        if (j == order[n]-1)
            JetsWithOrderedB.push_back(JetsWithB[j]);


//    std::cout << "JetsWithOrderedB index are " << JetsWithOrderedB.size() << ":  " << order[0] << "    " << order[1] << "    " << order[2] << "    " << order[3] << "    " << std::endl;

  return makeP4s(JetsWithOrderedB);
//  return makeP4s(JetsWithB);
}

 

std::vector<math::XYZTLorentzVector> BBClusterizer::GetLeadingJets() {

  vector<PseudoJet> jets = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets());
  std::vector<math::XYZTLorentzVector> LeadingJets;
//  math::XYZTLorentzVector nullLorentzVector = LorentzVector(0, 0, 0, 0);
  math::XYZTLorentzVector tempLorentzVector = LorentzVector(0, 0, 0, 0);

  for(unsigned int i = 0; i<10; i++) {

    if(jets.size() > i) {
      tempLorentzVector = LorentzVector( jets[i].px(), jets[i].py(), jets[i].pz(), jets[i].e() );
      LeadingJets.push_back(tempLorentzVector);
    }
//    else
//      LeadingJets.push_back(nullLorentzVector);
  }

  return LeadingJets;
}

int BBClusterizer::SVIsInTheJet(const fastjet::PseudoJet & pj) {
  int numberOfSV = 0;
  std::vector<fastjet::PseudoJet> constituents = pj.constituents();
  for (unsigned j = 0; j < constituents.size(); j++)
    if(constituents[j].user_index()<0)
      numberOfSV = numberOfSV + 1;
  return numberOfSV;
}


std::vector<PseudoJet>  BBClusterizer::DivideTheJet(const fastjet::PseudoJet & pj) {
    std::vector<PseudoJet> splittedBjets;
    PseudoJet parent1, parent2;
    bool had_parents = pj.has_parents(parent1,parent2);

    if(had_parents) {
        int parent1HasSV = SVIsInTheJet(parent1);
        int parent2HasSV = SVIsInTheJet(parent2);

        if ((parent1HasSV > 0) && (parent2HasSV > 0)) {
          splittedBjets.push_back(parent1);
          splittedBjets.push_back(parent2);
        }
        if (parent1HasSV == 0) splittedBjets = DivideTheJet(parent2);

        if (parent2HasSV == 0) splittedBjets = DivideTheJet(parent1);
    }

    return splittedBjets;
}

}




