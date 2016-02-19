//#include "PhysicsTools/Heppy/interface/ReclusterJets.h"
#include "BBAnalysis/BBHeppy/interface/SVUtils.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "fastjet/tools/Pruner.hh"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
using namespace std;

//using namespace std;
using namespace fastjet;

namespace heppy{

Measurement1D SVUtils::dist3D()
{
  VertexDistance3D dist;
  return dist.distance(sv_,pv_);
}
Measurement1D SVUtils::dist2D()
{
  VertexDistanceXY dist;
  return dist.distance(sv_,pv_);
}

}
