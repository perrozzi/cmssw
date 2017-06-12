from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR,deltaPhi
from copy import deepcopy
from math import *
import itertools
import ROOT
from sets import Set

def ptRel(p4,axis):
        a=ROOT.TVector3(axis.X(),axis.Y(),axis.Z())
        o=ROOT.TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E())
        return o.Perp(a)

class BBHeppy( Analyzer ):
    '''Analyze BB events
    '''

    def declareHandles(self):
        super(BBHeppy, self).declareHandles()
        self.handles['cands'] =  AutoHandle( 'packedPFCandidates','std::vector<pat::PackedCandidate>')
#        self.handles['SV'] =  AutoHandle( 'slimmedSecondaryVertices','std::vector<reco::VertexCompositePtrCandidate>')
        self.mchandles['packedGen'] = AutoHandle( 'packedGenParticles', 'std::vector<pat::PackedGenParticle>' )
        self.handles['PV'] =  AutoHandle( 'offlineSlimmedPrimaryVertices','std::vector<reco::Vertex>' )
    def beginLoop(self,setup):
        super(BBHeppy,self).beginLoop(setup)
        if "outputfile" in setup.services :
            setup.services["outputfile"].file.cd()
            self.inputCounter = ROOT.TH1F("Count","Count",1,0,2)
#            self.inputCounterWeighted = ROOT.TH1F("CountWeighted","Count with gen weight and pu weight",1,0,2)
#            self.inputCounterPosWeight = ROOT.TH1F("CountPosWeight","Count genWeight>0",1,0,2)
#            self.inputCounterNegWeight = ROOT.TH1F("CountNegWeight","Count genWeight<0",1,0,2)
#            for LHE_scale in range(6):
#               setattr(self, "inputCounterWeightedLHEWeightScale_"+str(LHE_scale), ROOT.TH1F("CountWeightedLHEWeightScale_"+str(LHE_scale),"Count with gen weight x LHE_weights_scale["+str(LHE_scale)+"] and pu weight",1,0,2))
#            for LHE_pdf in range(2):
#               setattr(self, "inputCounterWeightedLHEWeightPdf_"+str(LHE_pdf), ROOT.TH1F("CountWeightedLHEWeightPdf_"+str(LHE_pdf),"Count with gen weight x LHE_weights_pdf["+str(LHE_pdf)+"] and pu weight",1,0,2))
    def notDaughtersOf(self,anchestors,packedGens):
      result=[]
      for p in packedGens :
           mom = p.mother(0)
           found = False
           while mom:
               if mom in anchestors :
                  found=True
                  break
               mom = mom.mother(0) if mom.numberOfMothers() > 0 else None
           if not found :   
              result.append(p)
      return result

    def clusterizeGenParticles(self,event,b1,b2) :
         packedGens=list(self.mchandles['packedGen'].product())
         goodCands=self.notDaughtersOf([b1,b2],packedGens)
   #      goodCands = [x for x  in packedGens if x not in excludedCands ]
         lorentzVectorForFJ=ROOT.std.vector(ROOT.reco.Particle.LorentzVector)()
         map(lambda x:lorentzVectorForFJ.push_back(x.p4()), goodCands)
         clusterizer = ROOT.heppy.BBClusterizer(lorentzVectorForFJ,b1.p4(),b2.p4(),0,0.4);
         outJets=clusterizer.getBJets()
         return outJets

    def clusterize(self,event,b1,b2,alldaughters) :
         pfCands=list(self.handles['cands'].product())
         excludedCands = [x.get() for x in alldaughters] 
         goodPFCands = [x for x  in pfCands if x not in excludedCands ]

#         print len(pfCands),len(excludedCands), len(goodPFCands)
 
         lorentzVectorForFJ=ROOT.std.vector(ROOT.reco.Particle.LorentzVector)()
         map(lambda x:lorentzVectorForFJ.push_back(x.p4()), goodPFCands)
         clusterizer = ROOT.heppy.BBClusterizer(lorentzVectorForFJ,b1,b2,0,0.4);
         outJets=clusterizer.getBJets()
         return outJets

    def studyMergeBToD(self,event) :
        for svpair in itertools.combinations(event.ivf,2) :
              svs=sorted(svpair,key=lambda x: x.d3d.value())
              if deltaR(svs[0].direction,svs[1].direction) < 0.3 and svs[0].p4().M() > 1.2 and svs[1].p4().M()<2.0  :
                   daughters = Set()
                   map(daughters.add,svs[0].daughterPtrVector())
                   map(daughters.add,svs[1].daughterPtrVector())
                   thisPair=sum(daughters, ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.))
                   thisPair.B=event.ivf.index(svs[0]) 
                   thisPair.D=event.ivf.index(svs[1])
                   thisPair.ptRel=ptRel(svs[1].p4(),(svs[0].vertex()-event.PV.position()))
                   thisPair.deltaR=deltaR(svs[0].direction,svs[1].direction) 
                   thisPair.mcDeltaR=deltaR(svs[0].mcHadron.p4(),svs[1].mcHadron.p4()) if svs[0].mcHadron is not None and  svs[1].mcHadron is not None else -1
                   if thisPair.M() < 6 :
                       event.mergeablePairs.append(thisPair)  

    def numberOfSharedTracks(self,event,svs1,svs2) :
        tracksVector1 = svs1.daughterPtrVector()
        tracksVector2 = svs2.daughterPtrVector()
        shared = sum(1  for t in tracksVector1 if t in tracksVector2 )
        return shared

    def infForMatching(self,event, Hadrons) :
        final = [-1, -1, -1, -1, 20., 20.]
#    calculte all the distances between Bs and SVs
        allDistances = []
        HLen=len(Hadrons)
        svLen=len(event.selectedSVs)
        for H in Hadrons :
            distancesOfOneH = []
            svIdx=0
            for sv in event.selectedSVs :
                distance = deltaR(sv.direction, H.p4())
                distancesOfOneH.append(distance)
            allDistances.append(distancesOfOneH)
        minDist = 20
        secMinDist = 20
        firstHidx = 0
        secondHidx = 0
        firstSVidx = 0
        secondSVidx = 0
#       Loop fo matching
        for iteration in range(1,min(svLen,HLen,2)+1) :
            for k in range(0,HLen) :
                for l in range(0,svLen) :
                    if minDist>allDistances[k][l] :
                        secMinDist  = minDist
                        minDist     = allDistances[k][l]
                        firstHidx = k
                        secondHidx = firstHidx
                        firstSVidx = l
                        secondSVidx = firstSVidx
                    elif secMinDist>allDistances[k][l] :
                        secMinDist  = allDistances[k][l]
                        secondHidx = k
                        secondSVidx = l
#           record in final
            if iteration==1 :
                final[0] = firstHidx
                final[1] = firstSVidx
                final[4] = allDistances[firstHidx][firstSVidx]
                if firstHidx != secondHidx and firstSVidx != secondSVidx :
                    final[2] = secondHidx
                    final[3] = secondSVidx
                    final[5] = allDistances[secondHidx][secondSVidx]
                    return final
                if min(svLen,HLen)>1 :
                    minDist = 20
                    secMinDist = 20
#       "delete" the SV and H I have already matched
                    temporaneyList = []
                    for l in range(0,svLen) :
                        temporaneyList.append(21)
                        allDistances[firstHidx][l] = 21
                    allDistances[firstHidx] = temporaneyList
            if iteration==2 :
                final[2] = firstHidx
                final[3] = firstSVidx
                final[5] = allDistances[firstHidx][firstSVidx]
        return final


    def process(self, event):
	#print "Event number",event.iEv
        self.readCollections( event.input )
        self.inputCounter.Fill(1)
        if self.cfg_comp.isMC and False: #AR: disable for now
            genWeight = self.handles['GenInfo'].product().weight()
            self.inputCounterWeighted.Fill(1,copysign(1.0,genWeight)*event.puWeight)
            for LHE_scale in range(min(len(event.LHE_weights_scale),6)): 
               getattr(self, "inputCounterWeightedLHEWeightScale_"+str(LHE_scale)).Fill(1,copysign(1.0, genWeight)*event.puWeight*(event.LHE_weights_scale[LHE_scale]).wgt) 
            for LHE_pdf in range(min(len(event.LHE_weights_pdf),2)): 
               getattr(self, "inputCounterWeightedLHEWeightPdf_"+str(LHE_pdf)).Fill(1,copysign(1.0, genWeight)*event.puWeight*(event.LHE_weights_pdf[LHE_pdf]).wgt) 
            if genWeight > 0:
                self.inputCounterPosWeight.Fill(1)
            elif genWeight < 0:
                self.inputCounterNegWeight.Fill(1)
#        event.SVs=list(self.handles['SV'].product())
        event.PV=self.handles['PV'].product()[0]

        event.mergeablePairs = [] 
        event.bbPairSystem = [] 
        event.genBbPairSystem = [] 
        event.bjets = [] 
        event.genBjets = [] 
#       self.studyMergeBToD(event)

        for sv in event.ivf :
            sv.direction=(sv.vertex()-event.PV.position())
#            sv.directionUnit=sv.direction/sv.direction.mag()
        event.selectedSVs = [sv for sv in event.ivf if sv.p4().M() > 1.5 and sv.p4().M() <6.5 and sv.numberOfDaughters()>2 and abs(sv.direction.eta()) < 2 and sv.p4().pt() > 8.
          and sv.direction.perp2() < 4. and sv.dxy.significance() > 3 and sv.d3d.significance() > 5  and sv.cosTheta > 0.95 ] 
#        print len(event.selectedSVs)       , len(event.ivf) 
        if len(event.selectedSVs) != 2 and len(event.genBHadrons) < 2 :
          return False
        infBMatch = self.infForMatching(event, event.genBHadrons)
        infDMatch = self.infForMatching(event, event.genDHadrons)
        if len(event.selectedSVs) == 2  :
              svs=event.selectedSVs
              daughters = Set()
              map(daughters.add,svs[0].daughterPtrVector())
              map(daughters.add,svs[1].daughterPtrVector())
              thisPair=sum([x.p4() for x in daughters], ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.))
              thisPair.numberOfBinThisEvent=len(event.genBHadrons)
              thisPair.B0=event.ivf.index(svs[0])
              thisPair.B1=event.ivf.index(svs[1])
              thisPair.deltaR=deltaR(svs[0].direction,svs[1].direction)
              thisPair.deltaRpp=deltaR(svs[0].p4(),svs[1].p4())
              thisPair.mcDeltaR=deltaR(svs[0].mcHadron.p4(),svs[1].mcHadron.p4()) if svs[0].mcHadron is not None and  svs[1].mcHadron is not None else -1
              event.bjets=self.clusterize(event,svs[0].p4(),svs[1].p4(),daughters) 
              thisPair.deltaRjet=deltaR(event.bjets[0],event.bjets[1])
              thisPair.numberOfSharedTracks=self.numberOfSharedTracks(event,svs[0],svs[1])
#              thisPair.mcMatchFraction0=svs[0].mcMatchFraction
#              thisPair.mcMatchFraction1=svs[1].mcMatchFraction
              if infBMatch[1]==0 :
                thisPair.deltaRForBMatch0 = infBMatch[4]
                thisPair.deltaRForBMatch1 = infBMatch[5]
              else :
                thisPair.deltaRForBMatch0 = infBMatch[5]
                thisPair.deltaRForBMatch1 = infBMatch[4]
              if infDMatch[1]==0 :
                thisPair.deltaRForDMatch0 = infDMatch[4]
                thisPair.deltaRForDMatch1 = infDMatch[5]
              else :
                thisPair.deltaRForDMatch0 = infDMatch[5]
                thisPair.deltaRForDMatch1 = infDMatch[4]
              event.bbPairSystem.append(thisPair)

        if len(event.genBHadrons) == 2 :
              event.genBjets=self.clusterizeGenParticles(event,event.genBHadrons[0],event.genBHadrons[1])
              thisGenPair=event.genBjets[0]+event.genBjets[1]
              thisGenPair.numberOfSVinThisEvent=len(event.selectedSVs)
              thisGenPair.deltaRHad=deltaR(event.genBHadrons[0],event.genBHadrons[1])
              thisGenPair.deltaRJet=deltaR(event.genBjets[0],event.genBjets[1])
              thisGenPair.hadronPair=event.genBHadrons[0].p4()+event.genBHadrons[1].p4()
              if infBMatch[0]==0 :
                thisGenPair.deltaRForMatching0 = infBMatch[4]
                thisGenPair.deltaRForMatching1 = infBMatch[5]
              else :
                thisGenPair.deltaRForMatching0 = infBMatch[5]
                thisGenPair.deltaRForMatching1 = infBMatch[4]
              event.genBbPairSystem.append(thisGenPair)


        return True


