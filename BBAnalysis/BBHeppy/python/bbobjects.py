#!/bin/env python
from math import *
import ROOT
#from CMGTools.TTHAnalysis.signedSip import *
from PhysicsTools.Heppy.analyzers.objects.autophobj import *
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi


import copy, os

##------------------------------------------  
## JET
##------------------------------------------  

jetTypeBB = NTupleObjectType("jet",  baseObjectTypes = [ jetType ], variables = [
    NTupleVariable("mcIdx",   lambda x : x.mcJet.index if hasattr(x,"mcJet") and x.mcJet is not None else -1, int, mcOnly=False,help="index of the matching gen jet"),
 ])

genBbPairType = NTupleObjectType("genBbPair",  baseObjectTypes = [ fourVectorType ], variables = [
#     NTupleVariable("numberOfSVinThisEvent",  lambda x : x.numberOfSVinThisEvent, help="number of SV in this event"),
#     NTupleVariable("numberOfSelectedSVinThisEvent",  lambda x : x.numberOfSelectedSVinThisEvent, help="number of Selected SV in this event"),
     NTupleVariable("deltaRhad",  lambda x : x.deltaRHad, help="deltaR between B"),
     NTupleVariable("deltaRjet",  lambda x : x.deltaRJet, help="deltaR between the jets"),
     NTupleVariable("deltaRLastb",  lambda x : x.deltaRLastb, help="deltaR between last bs"),
     NTupleVariable("deltaRFirstb",  lambda x : x.deltaRFirstb, help="deltaR between the jets"),
     NTupleVariable("hadPt",  lambda x : x.hadronPair.pt(), help="pt of the hadron pair"),
     NTupleVariable("hadEta",  lambda x : x.hadronPair.eta(), help="eta of the hadron pair"),
     NTupleVariable("hadPhi",  lambda x : x.hadronPair.phi(), help="phi of the hadron pair"),
     NTupleVariable("hadMass",  lambda x : x.hadronPair.M(), help="mass of the hadron pair"),
     NTupleVariable("deltaRForMatching0",  lambda x : x.deltaRForMatching0, help="deltaR of matched SV with B0 "),
     NTupleVariable("deltaRForMatching1",  lambda x : x.deltaRForMatching1, help="deltaR of matched SV with B1 "),

])

bbPairType = NTupleObjectType("bbPair",  baseObjectTypes = [ fourVectorType ], variables = [
     NTupleVariable("numberOfBinThisEvent",  lambda x : x.numberOfBinThisEvent, mcOnly=True, help="number of B in this event"),
#     NTupleVariable("numberOfSVinThisEvent",  lambda x : x.numberOfSVinThisEvent, help="number of SV in this event"),
     NTupleVariable("B0idx",  lambda x : x.B0, help="index of B0 vertex"),
     NTupleVariable("B1idx",  lambda x : x.B1, help="index of B1 vertex"),
     NTupleVariable("deltaR",  lambda x : x.deltaR, help="deltaR between B with SVs"),
     NTupleVariable("deltaRpp",  lambda x : x.deltaRpp, help="deltaR between B with momenta"),
     NTupleVariable("deltaRjet",  lambda x : x.deltaRjet, help="deltaR between the jets"),
     NTupleVariable("mcDeltaR",  lambda x : x.mcDeltaR, mcOnly=True, help="deltaR between Hadrons"),
     NTupleVariable("nSharedTracks",  lambda x : x.numberOfSharedTracks, help="number of shared tracks between the SVs"),
     NTupleVariable("maxSIPShared",  lambda x : x.SIPsOfTheShareds[0], help="max 3d impact parameter significance of the the shared tracks"),
     NTupleVariable("secSIPShared",  lambda x : x.SIPsOfTheShareds[1], help="second max impact parameter significance of the the shared tracks"),
     NTupleVariable("deltaRForBMatch0",  lambda x : x.deltaRForBMatch0, mcOnly=True, help="deltaR of matched B with SV0"),
     NTupleVariable("deltaRForBMatch1",  lambda x : x.deltaRForBMatch1, mcOnly=True, help="deltaR of matched B with SV1"),
     NTupleVariable("deltaRForDMatch0",  lambda x : x.deltaRForDMatch0, mcOnly=True, help="deltaR of matched D with SV0"),
     NTupleVariable("deltaRForDMatch1",  lambda x : x.deltaRForDMatch1, mcOnly=True, help="deltaR of matched D with SV1"),
     NTupleVariable("jetsPt",  lambda x : x.jetsPt, help="Pt of the jet Pair"),
     NTupleVariable("jetsMass",  lambda x : x.jetsMass, help="invariant mass of the jet Pair"),
     NTupleVariable("jetsPtByMass",  lambda x : x.jetsPtByMass, help="Pt divided by mass of the jet Pair"),
])

bQuarkType = NTupleObjectType("bQuarks",  baseObjectTypes = [ fourVectorType ], variables = [
     NTupleVariable("deltaRwithB",  lambda x : x.deltaRwithB, help="delta R between this b and the B hadron"),
])

vertexPairType = NTupleObjectType("vertexPair",  baseObjectTypes = [ fourVectorType ], variables = [
     NTupleVariable("Bidx",  lambda x : x.B, help="index of B vertex"),
     NTupleVariable("Didx",  lambda x : x.D, help="index of D vertex"),
     NTupleVariable("ptRel",  lambda x : x.ptRel, help="ptRel of D wrt B-PV direction"),
     NTupleVariable("deltaR",  lambda x : x.deltaR, help="deltaR between B and D"),
     NTupleVariable("mcDeltaR",  lambda x : x.mcDeltaR, help="deltaR between Hadrons"),

])

vertexType = NTupleObjectType("vertex",  baseObjectTypes = [ fourVectorType ], variables = [
#    NTupleVariable("x",  lambda x : x.x(), help="x coordinate"),
#    NTupleVariable("y",  lambda x : x.y(), help="y coordinate"),
#    NTupleVariable("z",  lambda x : x.z(), help="z coordinate"),
])

GenDfromBHadType = NTupleObjectType("genBToDHadrons",  baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("PtRel",  lambda x : x.PtRel, help="PtRel of D respecto to B flight direction"),
    NTupleVariable("BToDdeltaR2",  lambda x : x.BtoDdeltaR, help="delta R squared between B and D"),
#    NTupleVariable("y",  lambda x : x.y(), help="y coordinate"),
#    NTupleVariable("z",  lambda x : x.z(), help="z coordinate"),
])
 

heavyFlavourHadronType = NTupleObjectType("heavyFlavourHadron", baseObjectTypes = [ genParticleType ], variables = [
    NTupleVariable("flav", lambda x : x.flav, int, mcOnly=True, help="Flavour"),
    NTupleVariable("sourceId", lambda x : x.sourceId, int, mcOnly=True, help="pdgId of heaviest mother particle (stopping at the first one heaviest than 175 GeV)"),
    NTupleVariable("bIndex",  lambda x : x.bIndex, help="index of b in bQuarks"),
    NTupleVariable("svMass",   lambda x : x.sv.mass() if x.sv else 0, help="SV: mass"),
    NTupleVariable("svPt",   lambda x : x.sv.pt() if x.sv else 0, help="SV: pt"),
    NTupleVariable("svCharge",   lambda x : x.sv.charge() if x.sv else -99., int, help="SV: charge"),
    NTupleVariable("svNtracks", lambda x : x.sv.numberOfDaughters() if x.sv else 0, int, help="SV: Number of tracks (with weight > 0.5)"),
    NTupleVariable("svChi2", lambda x : x.sv.vertexChi2() if x.sv else -99., help="SV: Chi2 of the vertex fit"),
    NTupleVariable("svNdof", lambda x : x.sv.vertexNdof() if x.sv else -99., help="SV: Degrees of freedom of the fit, ndof = (2*ntracks - 3)" ),
    NTupleVariable("svDxy",  lambda x : x.sv.dxy.value() if x.sv else -99., help="SV: Transverse distance from the PV [cm]"),
    NTupleVariable("svEdxy", lambda x : x.sv.dxy.error() if x.sv else -99., help="SV: Uncertainty on the transverse distance from the PV [cm]"),
    NTupleVariable("svIp3d",  lambda x : x.sv.d3d.value() if x.sv else -99., help="SV: 3D distance from the PV [cm]"),
    NTupleVariable("svEip3d", lambda x : x.sv.d3d.error() if x.sv else -99., help="SV: Uncertainty on the 3D distance from the PV [cm]"),
    NTupleVariable("svSip3d", lambda x : x.sv.d3d.significance() if x.sv else -99., help="SV: S_{ip3d} with respect to PV (absolute value)"),
    NTupleVariable("svCosTheta", lambda x : x.sv.cosTheta if x.sv else -99., help="SV: Cosine of the angle between the 3D displacement and the momentum"),
#    NTupleVariable("jetPt",  lambda x : x.jet.pt() if x.jet != None else 0, help="Jet: pT"),
#    NTupleVariable("jetBTag",  lambda x : x.jet.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags') if x.jet != None else -99, help="CSV b-tag of associated jet"),
])

primaryVertexType = NTupleObjectType("primaryVertex", variables = [
    NTupleVariable("x",    lambda x : x.x()),
    NTupleVariable("y",   lambda x : x.y()),
    NTupleVariable("z",   lambda x : x.z()),
    NTupleVariable("isFake",   lambda x : x.isFake()),
    NTupleVariable("ndof",   lambda x : x.ndof()),
    NTupleVariable("Rho",   lambda x : x.position().Rho()),
    NTupleVariable("score",  lambda x : x.score),
])

svType = NTupleObjectType("sv", baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
#<<<<<<< HEAD
##    NTupleVariable("dirEta",   lambda x : x.direction.eta()),
##    NTupleVariable("dirPhi",   lambda x : x.direction.phi()),
#=======
    NTupleVariable("x",   lambda x : x.CMSCoordinates.x(), help="x coordinate in CMS"),
    NTupleVariable("y",   lambda x : x.CMSCoordinates.y(), help="y coordinate in CMS"),
    NTupleVariable("z",   lambda x : x.CMSCoordinates.z(), help="z coordinate in CMS"),
    NTupleVariable("CMSeta",   lambda x : x.CMSCoordinates.eta(), help="eta coordinate in CMS"),
    NTupleVariable("CMSphi",   lambda x : x.CMSCoordinates.phi(), help="phi coordinate in CMS"),
    NTupleVariable("dirEta",   lambda x : x.direction.eta()),
    NTupleVariable("dirPhi",   lambda x : x.direction.phi()),
    NTupleVariable("PtRel",   lambda x : x.PtRel, help="PtRel between momentum and flight direction"),
#>>>>>>> nuovoHbb
    NTupleVariable("ntracks", lambda x : x.numberOfDaughters(), int, help="Number of tracks (with weight > 0.5)"),
    NTupleVariable("isMerged", lambda x : x.ism, help="number of times this vertex merged"),
    NTupleVariable("isBSelected", lambda x : x.isBSelected, help="boolean true if it is one of the two ivf selected as B"),
    NTupleVariable("chi2", lambda x : x.vertexChi2(), help="Chi2 of the vertex fit"),
    NTupleVariable("ndof", lambda x : x.vertexNdof(), help="Degrees of freedom of the fit, ndof = (2*ntracks - 3)" ),
    NTupleVariable("dxy",  lambda x : x.dxy.value(), help="Transverse distance from the PV [cm]"),
    NTupleVariable("sdxy",  lambda x : x.dxy.significance(), help="Significance Transverse distance from the PV "),
    NTupleVariable("edxy", lambda x : x.dxy.error(), help="Uncertainty on the transverse distance from the PV [cm]"),
    NTupleVariable("ip3d",  lambda x : x.d3d.value(), help="3D distance from the PV [cm]"),
    NTupleVariable("eip3d", lambda x : x.d3d.error(), help="Uncertainty on the 3D distance from the PV [cm]"),
    NTupleVariable("sip3d", lambda x : x.d3d.significance(), help="S_{ip3d} with respect to PV (absolute value)"),
    NTupleVariable("cosTheta", lambda x : x.cosTheta, help="Cosine of the angle between the 3D displacement and the momentum"),
#    NTupleVariable("mva", lambda x : x.mva, help="MVA discriminator"),
 #   NTupleVariable("jetPt",  lambda x : x.jet.pt() if x.jet != None else 0, help="pT of associated jet"),
 #   NTupleVariable("jetBTagCSV",   lambda x : x.jet.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags') if x.jet != None else -99, help="CSV b-tag of associated jet"),
 #   NTupleVariable("jetBTagCMVA",  lambda x : x.jet.btag('pfCombinedMVABJetTags') if x.jet != None else -99, help="CMVA b-tag of associated jet"),
    NTupleVariable("mcMatchNTracks", lambda x : getattr(x, 'mcMatchNTracks', -1), int, mcOnly=True, help="Number of mc-matched tracks in SV"),
    NTupleVariable("mcMatchNTracksHF", lambda x : getattr(x, 'mcMatchNTracksHF', -1), int, mcOnly=True, help="Number of mc-matched tracks from b/c in SV"),
    NTupleVariable("mcMatchFraction", lambda x : getattr(x, 'mcMatchFraction', -1), mcOnly=True, help="Fraction of mc-matched tracks from b/c matched to a single hadron (or -1 if mcMatchNTracksHF < 2)"),
    NTupleVariable("mcFlavFirst", lambda x : getattr(x,'mcFlavFirst', -1), int, mcOnly=True, help="Flavour of last ancestor with maximum number of matched daughters"),
    NTupleVariable("mcFlavHeaviest", lambda x : getattr(x,'mcFlavHeaviest', -1), int, mcOnly=True, help="Flavour of heaviest hadron with maximum number of matched daughters"),
    NTupleVariable("mcNDifferentBHads", lambda x : getattr(x,'mcNDifferentBHads', -1), int, mcOnly=True, help="Number of different B-hadrons when no B parents are available"),
    NTupleVariable("maxDxyTracks", lambda x : x.maxDxyTracks, help="highest |dxy| of vertex tracks"),
    NTupleVariable("secDxyTracks", lambda x : x.secDxyTracks, help="second highest |dxy| of vertex tracks"),
    NTupleVariable("maxD3dTracks", lambda x : x.maxD3dTracks, help="highest |ip3D| of vertex tracks"),
    NTupleVariable("secD3dTracks", lambda x : x.secD3dTracks, help="second highest |ip3D| of vertex tracks"),

])



triggerObjectsType = NTupleObjectType("triggerObjects",  baseObjectTypes = [ fourVectorType ], variables = [
])

triggerObjectsOnlyPtType = NTupleObjectType("triggerObjects",  baseObjectTypes = [ ], variables = [
    NTupleVariable("pt", lambda x : x.pt(), float, mcOnly=False, help="trigger object pt"),
])

triggerObjectsNothingType = NTupleObjectType("triggerObjects",  baseObjectTypes = [ ], variables = [
])



#JetsType = NTupleObjectType("Jets",  baseObjectTypes = [ fourVectorType ], variables = [
#])



