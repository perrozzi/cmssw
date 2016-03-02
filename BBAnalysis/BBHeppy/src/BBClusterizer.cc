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
    if(SVIsInTheJet(*it)) JetsWithB.push_back(*it);


  PseudoJet parent1, parent2;
  bool had_parents = JetsWithB[0].has_parents(parent1,parent2);
 
  while ((JetsWithB.size()==1) && (had_parents)) {
      bool parent1HasSV = SVIsInTheJet(parent1);
      bool parent2HasSV = SVIsInTheJet(parent2);

    if ((parent1HasSV) && (parent2HasSV)) {
      JetsWithB[0] = parent1;
      JetsWithB.push_back(parent2);
    }
    if ((parent1HasSV) && (!parent2HasSV)) {
      had_parents = parent1.has_parents(parent1,parent2);
    }
    if ((!parent1HasSV) && (parent2HasSV)) {
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

math::XYZTLorentzVector BBClusterizer::GetLeadingJet() {

  vector<PseudoJet> jets = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets());
  return LorentzVector( jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].e() );

}

bool BBClusterizer::SVIsInTheJet(const fastjet::PseudoJet & pj) {
  bool pjHasSV = false;
  std::vector<fastjet::PseudoJet> constituents = pj.constituents();
  for (unsigned j = 0; j < constituents.size(); j++)
    if(constituents[j].user_index()<0)
      pjHasSV = true;
  return pjHasSV;
}

}

