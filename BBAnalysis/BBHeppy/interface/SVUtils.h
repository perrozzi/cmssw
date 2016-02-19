#ifndef PhysicsTools_Heppy_SVUtils_h
#define PhysicsTools_Heppy_SVUtils_h

#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include "DataFormats/Math/interface/LorentzVector.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"
#include <RecoVertex/VertexPrimitives/interface/BasicVertexState.h>


namespace heppy{
class SVUtils {
    
 public:
  SVUtils(const reco::Vertex & pv, const  reco::VertexCompositePtrCandidate & sv): 
  pv_(RecoVertex::convertPos(pv.position()),RecoVertex::convertError(pv.error())),
  sv_(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error()))
  {}
  Measurement1D dist3D();
  Measurement1D dist2D();

 private:
    VertexState pv_;
    VertexState sv_;
};
}
#endif   
 
