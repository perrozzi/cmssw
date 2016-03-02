#ifndef PhysicsTools_Heppy_BBClusterizer_h
#define PhysicsTools_Heppy_BBClusterizer_h

#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include "DataFormats/Math/interface/LorentzVector.h"

#include <boost/shared_ptr.hpp>
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

namespace heppy{
class BBClusterizer {
    
 public:
  typedef math::XYZTLorentzVector LorentzVector;

  BBClusterizer(const std::vector<LorentzVector> & ,
		   const LorentzVector &,const LorentzVector &,
                   double ktpower, double rparam);

  /// get grouping (inclusive jets)
  std::vector<LorentzVector> getBJets();
  LorentzVector GetLeadingJet();

 private:
  // pack the returns in a fwlite-friendly way
  std::vector<LorentzVector> makeP4s(const std::vector<fastjet::PseudoJet> &jets) ;

  // look if there is the SV in the jet
  bool SVIsInTheJet(const fastjet::PseudoJet & pj);

  // used to handle the inputs
  std::vector<fastjet::PseudoJet> fjInputs;        // fastjet inputs

  LorentzVector b1_;
  LorentzVector b2_;
  double ktpower_;
  double rparam_;
 
  /// fastjet outputs
  typedef boost::shared_ptr<fastjet::ClusterSequence>  ClusterSequencePtr;
  ClusterSequencePtr fjClusterSeq_;    
};
}
#endif   
 
