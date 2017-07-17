from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import NTupleVariable
from PhysicsTools.HeppyCore.utils.deltar import matchObjectCollection, matchObjectCollection3
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.utils.deltar import deltaR,deltaPhi
from BBAnalysis.BBHeppy.signedSip import SignedImpactParameterComputer
from copy import deepcopy
from math import *
import itertools
import ROOT
import math
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
        self.handles['slimmedJets'] = AutoHandle( 'slimmedJets', 'std::vector<pat::Jet>' )
        self.mchandles['packedGen'] = AutoHandle( 'packedGenParticles', 'std::vector<pat::PackedGenParticle>' )
        self.handles['PV'] =  AutoHandle( 'offlineSlimmedPrimaryVertices','std::vector<reco::Vertex>' )
        self.mchandles['pileUp_source'] =  AutoHandle( 'slimmedAddPileupInfo','vector<PileupSummaryInfo>' )
        self.mchandles['generator_source'] =  AutoHandle( 'generator','GenEventInfoProduct' )

#        triggerObjectsCfgs = getattr(self.cfg_ana,"triggerObjectsCfgs",[])
#        self.triggerObjectInputTag = getattr(self.cfg_ana,"triggerObjectInputTag",("","",""))
#        self.handles['TriggerBits']     = AutoHandle( self.triggerBitsInputTag, 'edm::TriggerResults' )
#        self.handles['TriggerObjects']  = AutoHandle( self.triggerObjectInputTag, 'std::vector<pat::TriggerObjectStandAlone>' )

    def beginLoop(self,setup):
        super(BBHeppy,self).beginLoop(setup)
        if "outputfile" in setup.services :
            setup.services["outputfile"].file.cd()
        elif 'PhysicsTools.HeppyCore.framework.services.tfile.TFileService_outputfile' in setup.services :
            setup.services['PhysicsTools.HeppyCore.framework.services.tfile.TFileService_outputfile'].file.cd()
        else:
            print "no outputfile in services:"
            print setup.services
        
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

    def myDeltaR(self, pvToDdirection, pvToBdirection) :
        Phi = ROOT.TVector2.Phi_mpi_pi(pvToDdirection.Phi()-pvToBdirection.Phi())
        Eta = pvToDdirection.PseudoRapidity()-pvToBdirection.PseudoRapidity()
        R = Phi*Phi + Eta*Eta
        return R

    def clusterizeGenParticles(self,event,b1,b2) :
         packedGens=list(self.mchandles['packedGen'].product())
         goodCands=self.notDaughtersOf([b1,b2],packedGens)
   #      goodCands = [x for x  in packedGens if x not in excludedCands ]
         lorentzVectorForFJ=ROOT.std.vector(ROOT.reco.Particle.LorentzVector)()
         map(lambda x:lorentzVectorForFJ.push_back(x.p4()), goodCands)
         clusterizer = ROOT.heppy.BBClusterizer(lorentzVectorForFJ,b1.p4(),b2.p4(),0,0.4);
         outJets=clusterizer.getBJets()
         return outJets

    def clusterizeGenParticles4B(self,event,b1,b2,b3,b4) :
         packedGens=list(self.mchandles['packedGen'].product())
         goodCands=self.notDaughtersOf([b1,b2,b3,b4],packedGens)
         lorentzVectorForFJ=ROOT.std.vector(ROOT.reco.Particle.LorentzVector)()
         map(lambda x:lorentzVectorForFJ.push_back(x.p4()), goodCands)
         clusterizer = ROOT.heppy.BBClusterizer(lorentzVectorForFJ,b1.p4(),b2.p4(),b3.p4(),b4.p4(),0,0.4);
#         clusterizer = ROOT.heppy.BBClusterizer(lorentzVectorForFJ,b1.p4(),b2.p4(),0,0.4);
         outJets=clusterizer.getBJets4B()
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
         event.CAJets = clusterizer.GetLeadingJets()
         return outJets

    def  svIsSelected(self, sv) :
        isSelected = False
        if sv.p4().M() > 1. and sv.p4().M() <5. and sv.numberOfDaughters()>1 and abs(sv.direction.eta()) < 2.0 and sv.p4().pt() > 8. :
         if sv.d3d.significance() > 5  and sv.dxy.significance() > 3 and sv.cosTheta > 0.95  :
          if  sv.direction.perp2()*sv.p4().M()*sv.p4().M()/sv.p4().pt()/sv.p4().pt() > 0.00025 or sv.p4().M() > 2.:
           if sv.direction.perp2()*sv.p4().M()*sv.p4().M()/sv.p4().pt()/sv.p4().pt() < 0.15 and  sv.p4().pt() < 400 :
            if log10(sv.direction.perp2()*sv.p4().M()*sv.p4().M()) - 0.4*log10(sv.p4().pt()) < 1.9  :#  and sv.direction.perp2() < 4. 
             if sv.secD3dTracks > 4. and sv.PtRel < 2.5 and sv.PtRel > 0 :
              isSelected = True
        return isSelected

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

    def selectVertices(self,selectedSVs) :
        selectedSVs.sort(key = lambda sv : abs(sv.p4().M()), reverse = True)
        newSV = [sv for sv in selectedSVs if sv.numberOfDaughters()> 2]
        if len(newSV) == 2 :
            if newSV[0].p4().M() > 2. and newSV[0].p4().M() + newSV[1].p4().M() < 9. :
                return newSV
        elif len(selectedSVs) == 2 :
            if selectedSVs[0].p4().M() > 2. and selectedSVs[0].p4().M() + selectedSVs[1].p4().M() < 9. :
                if selectedSVs[0].numberOfDaughters()> 2 or selectedSVs[1].numberOfDaughters()> 2 :
                    return selectedSVs
        return []

    def generatedFromTheSameGluonSplitting(self, BHadrons, svs) :
        BHadronIdx1 = 0
        svsIdx = 0
        minDeltaRdistance = 100
        for n in range(0, len(BHadrons)) :
          for m in range(0, len(svs)) :
            distanceInR = deltaR(BHadrons[n].p4(), svs[m].p4())
            if distanceInR < minDeltaRdistance :
                BHadronIdx1 = n
                svsIdx = m
                minDeltaRdistance = distanceInR


        b1 = BHadrons[BHadronIdx1].firstb
        BHadronIdx2 = 0
        minDeltaRdistance = 100
        for n in range(0, len(BHadrons)) :
            if n != BHadronIdx1 :
                for m in range(0, len(svs)) :
                    if m != svsIdx :
                        distanceInR = deltaR(BHadrons[n].p4(), svs[m].p4())
                        if distanceInR < minDeltaRdistance :
                            BHadronIdx2 = n
                            minDeltaRdistance = distanceInR

        b2 = BHadrons[BHadronIdx2].firstb
        mom1 = b1.motherRefVector()[0] 
        mom2 = b2.motherRefVector()[0] 
        if mom1 == None or mom1.isNull() or not mom1.isAvailable(): 
            print "ERROR1"
        if mom2 == None or mom2.isNull() or not mom2.isAvailable(): 
            print "ERROR2"
        if mom1 == mom2 :
            return True
        else :
            return False


    def numberOfSharedTracks(self,event,svs1,svs2) :
        tracksVector1 = svs1.daughterPtrVector()
        tracksVector2 = svs2.daughterPtrVector()
        shared = sum(1  for t in tracksVector1 if t in tracksVector2 )
        return shared

    def SIPsOfTheShareds(self,event,svs1,svs2) :
        SIPs = [0, 0]
        AllSIPs = []
        for dau in svs1.daughterPtrVector() :
            if dau in svs2.daughterPtrVector() :
                sharedTrack = dau.get()
                AllSIPs.append(SignedImpactParameterComputer.signedIP3D(sharedTrack.pseudoTrack(), event.PV, svs1.momentum()).significance())
        AllSIPs.sort(key = lambda s : s, reverse = True)
        if len(AllSIPs) > 0 :
            SIPs[0] = AllSIPs[0]
        if len(AllSIPs) > 1 :
            SIPs[1] = AllSIPs[1]
        return SIPs

    def infForMatching(self,event, Hadrons) :
      final = [-1, -1, -1, -1, 20., 20.]
      HLen=len(Hadrons)
      svLen=len(event.selectedSVsSelected)
      if HLen>0 and svLen>0 :
#    calculte all the distances between Bs and SVs
        allDistances = []
        for H in Hadrons :
            distancesOfOneH = []
            for sv in event.selectedSVsSelected :
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
        for iteration in range(0,min(svLen,HLen,2)) :
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
            if iteration==0 :
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
                    for k in range(0,HLen) :
                        allDistances[k][firstSVidx] = 21
                    for l in range(0,svLen) :
                        temporaneyList.append(21)
                    allDistances[firstHidx] = temporaneyList
            if iteration==1 :
                final[2] = firstHidx
                final[3] = firstSVidx
                final[5] = allDistances[firstHidx][firstSVidx]
      return final


    def process(self, event):

	#print "Event number",event.iEv
        self.readCollections( event.input )
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
#        allTriggerObjects = self.handles['TriggerObjects'].product()
#        triggerBits = self.handles['TriggerBits'].product()

#       keep the event only if PV is the vertex with bigger pt_hat
        if self.cfg_comp.isMC :
            event.generatorSource = self.mchandles['generator_source']
            event.pileUpSource = self.mchandles['pileUp_source']
            event.ptHat = event.generatorSource.product().qScale()
            maxPUptHat = -1
            for PUInteraction in range(event.pileUpSource.product().size()) :
                if event.pileUpSource.product().at(PUInteraction).getBunchCrossing() == 0 :
                    for PUpyHadIterator in event.pileUpSource.product().at(PUInteraction).getPU_pT_hats() :
                        maxPUptHat = max(maxPUptHat, PUpyHadIterator) 
    #        print event.ptHat, maxPUptHat
            event.maxPUptHat = maxPUptHat
            if event.ptHat < maxPUptHat :
                return False
    #            print "hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
            self.inputCounter.Fill(1)
        else :
            self.inputCounter.Fill(1)


        event.isSignal = False
        event.mergeablePairs = [] 
        event.bbPairSystem = [] 
        event.genBbPairSystem = [] 
        event.bjets = [] 
        event.genBjets = [] 
        event.CAJets=[]
        event.AK4Jets=[]
#       self.studyMergeBToD(event)

        if self.cfg_comp.isMC :
          for D in event.genBToDHadrons :
            D.PtRel = -1
            D.BtoDdeltaR = -1
            if D.BDecayPoint != None :
                pvToB=(D.BDecayPoint-event.PV.position())
                pvToBdirection = ROOT.TVector3(pvToB.x(), pvToB.y(), pvToB.z())
                momentum = ROOT.TVector3(D.X(), D.Y(), D.Z())
                cos = momentum.Dot(pvToBdirection)/momentum.Mag()/pvToBdirection.Mag()
                sin = math.sqrt(1-cos*cos)
                D.PtRel = sin*momentum.Mag()
#                if D.DDecayPoint != None :
#                    pvToD=(D.DDecayPoint-event.PV.position())
#                    pvToDdirection = ROOT.TVector3(pvToD.x(), pvToD.y(), pvToD.z())
#                    D.BtoDdeltaR = self.myDeltaR(pvToDdirection, pvToBdirection)


        for sv in event.ivf :
            sv.isBSelected = False
            sv.direction=(sv.vertex()-event.PV.position())
#<<<<<<< HEAD
#            sv.directionUnit=sv.direction/sv.direction.mag()
#        event.selectedSVs = [sv for sv in event.ivf if sv.p4().M() > 1.5 and sv.p4().M() <6.5 and sv.numberOfDaughters()>2 and abs(sv.direction.eta()) < 2 and sv.p4().pt() > 8.
#          and sv.direction.perp2() < 4. and sv.dxy.significance() > 3 and sv.d3d.significance() > 5  and sv.cosTheta > 0.95 ] 

        #print len(event.selectedSVs)       , len(event.ivf) 
        
#        if self.cfg_comp.isMC:
#            if len(event.selectedSVs) != 2 and len(event.genBHadrons) < 2 :
#              return False
#            infBMatch = self.infForMatching(event, event.genBHadrons)
#            infDMatch = self.infForMatching(event, event.genDHadrons)
#
#        if len(event.selectedSVs) == 2  :
#              svs=event.selectedSVs
#=======
            sv.CMSCoordinates = sv.vertex()
            SVDirection = ROOT.TVector3(sv.direction.x(), sv.direction.y(), sv.direction.z())
            momentum = ROOT.TVector3(sv.p4().X(), sv.p4().Y(), sv.p4().Z())
            sv.PtRel = -1
            if momentum.Mag() > 0.1 :
                cosTheta = momentum.Dot(SVDirection)/SVDirection.Mag()/momentum.Mag()
                sinTheta = math.sqrt(1-cosTheta*cosTheta)
                sv.PtRel = sinTheta*momentum.Mag()
        event.selectedSVs = [sv for sv in event.ivf if self.svIsSelected(sv)]
#        print len(event.selectedSVs)       , len(event.ivf) 
        event.selectedSVsSelected = self.selectVertices(event.selectedSVs)
        if len(event.ivf)  < 1 :
            return False
        infBMatch = 0
        infDMatch = 0
        if self.cfg_comp.isMC :
            infBMatch = self.infForMatching(event, event.genBHadrons)
            infDMatch = self.infForMatching(event, event.genDHadrons)
        if len(event.selectedSVsSelected) == 2  :
          SIPsOfTheShareds_temp_for_check=self.SIPsOfTheShareds(event,event.selectedSVsSelected[0],event.selectedSVsSelected[1])
          if SIPsOfTheShareds_temp_for_check[0] < 6 :
              svs=event.selectedSVsSelected
#>>>>>>> nuovoHbb
              daughters = Set()
              map(daughters.add,svs[0].daughterPtrVector())
              map(daughters.add,svs[1].daughterPtrVector())
              thisPair=sum([x.p4() for x in daughters], ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.))
#<<<<<<< HEAD
#              if self.cfg_comp.isMC: 
#                  thisPair.numberOfBinThisEvent=len(event.genBHadrons)
#              else:
#                  thisPair.numberOfBinThisEvent = -1
#
#=======
#              thisPair.numberOfSVinThisEvent=len(event.selectedSVs)
              if self.cfg_comp.isMC :
                thisPair.numberOfBinThisEvent=len(event.genBHadrons)
                thisPair.mcDeltaR=deltaR(svs[0].mcHadron.p4(),svs[1].mcHadron.p4()) if svs[0].mcHadron is not None and  svs[1].mcHadron is not None else -1
#>>>>>>> nuovoHbb
              thisPair.B0=event.ivf.index(svs[0])
              thisPair.B1=event.ivf.index(svs[1])
              event.ivf[event.ivf.index(svs[0])].isBSelected = True
              event.ivf[event.ivf.index(svs[1])].isBSelected = True
              thisPair.deltaR=deltaR(svs[0].direction,svs[1].direction)
              thisPair.deltaRpp=deltaR(svs[0].p4(),svs[1].p4())
#<<<<<<< HEAD
#              if self.cfg_comp.isMC: 
#                  thisPair.mcDeltaR=deltaR(svs[0].mcHadron.p4(),svs[1].mcHadron.p4()) if svs[0].mcHadron is not None and  svs[1].mcHadron is not None else -1
#              else:
#                  thisPair.mcDeltaR=-1
#=======
#>>>>>>> nuovoHbb
              event.bjets=self.clusterize(event,svs[0].p4(),svs[1].p4(),daughters) 
              #AK4=self.handles['slimmedJets'].product()
              AK4=list(self.handles['slimmedJets'].product())
              event.AK4Jets = AK4
              thisPair.deltaRjet=deltaR(event.bjets[0],event.bjets[1])
              thisPair.numberOfSharedTracks=self.numberOfSharedTracks(event,svs[0],svs[1])

#<<<<<<< HEAD
##              thisPair.mcMatchFraction0=svs[0].mcMatchFraction
##              thisPair.mcMatchFraction1=svs[1].mcMatchFraction
#              if self.cfg_comp.isMC:
#                  if infBMatch[1]==0 :
#                    thisPair.deltaRForBMatch0 = infBMatch[4]
#                    thisPair.deltaRForBMatch1 = infBMatch[5]
#                  else :
#                    thisPair.deltaRForBMatch0 = infBMatch[5]
#                    thisPair.deltaRForBMatch1 = infBMatch[4]
#                  if infDMatch[1]==0 :
#                    thisPair.deltaRForDMatch0 = infDMatch[4]
#                    thisPair.deltaRForDMatch1 = infDMatch[5]
#                  else :
#                    thisPair.deltaRForDMatch0 = infDMatch[5]
#                    thisPair.deltaRForDMatch1 = infDMatch[4]
#              else:
#                    thisPair.deltaRForBMatch0 = -1
#                    thisPair.deltaRForBMatch1 = -1
#                    thisPair.deltaRForDMatch0 = -1
#                    thisPair.deltaRForDMatch1 = -1
#              event.bbPairSystem.append(thisPair)
#        
#        if self.cfg_comp.isMC:
#            if len(event.genBHadrons) == 2 :
#                  event.genBjets=self.clusterizeGenParticles(event,event.genBHadrons[0],event.genBHadrons[1])
#                  thisGenPair=event.genBjets[0]+event.genBjets[1]
#                  thisGenPair.numberOfSVinThisEvent=len(event.selectedSVs)
#                  thisGenPair.deltaRHad=deltaR(event.genBHadrons[0],event.genBHadrons[1])
#                  thisGenPair.deltaRJet=deltaR(event.genBjets[0],event.genBjets[1])
#                  thisGenPair.hadronPair=event.genBHadrons[0].p4()+event.genBHadrons[1].p4()
#                  if infBMatch[0]==0 :
#                    thisGenPair.deltaRForMatching0 = infBMatch[4]
#                    thisGenPair.deltaRForMatching1 = infBMatch[5]
#                  else :
#                    thisGenPair.deltaRForMatching0 = infBMatch[5]
#                    thisGenPair.deltaRForMatching1 = infBMatch[4]
#                  event.genBbPairSystem.append(thisGenPair)
#=======
              thisPair.SIPsOfTheShareds=SIPsOfTheShareds_temp_for_check
              jetsMomentum = event.bjets[0] + event.bjets[1]
              thisPair.jetsPt = jetsMomentum.pt()
              thisPair.jetsMass = jetsMomentum.mass()
              thisPair.jetsPtByMass = jetsMomentum.pt()/jetsMomentum.mass()
              if self.cfg_comp.isMC :
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

              if thisPair.SIPsOfTheShareds[0] < 4 :
                event.bbPairSystem.append(thisPair)


        if self.cfg_comp.isMC :
          if len(event.genBHadrons) == 2 :
              event.isSignal = True
              event.genBjets=self.clusterizeGenParticles(event,event.genBHadrons[0],event.genBHadrons[1])
              thisGenPair=event.genBjets[0]+event.genBjets[1]
#              thisGenPair.numberOfSVinThisEvent=len(event.selectedSVs)
#              thisGenPair.numberOfSelectedSVinThisEvent=len(event.selectedSVsSelected)
              thisGenPair.deltaRHad=deltaR(event.genBHadrons[0],event.genBHadrons[1])
              thisGenPair.deltaRJet=deltaR(event.genBjets[0],event.genBjets[1])
              thisGenPair.deltaRLastb =deltaR(event.genBHadrons[0].lastb,event.genBHadrons[1].lastb)
              thisGenPair.deltaRFirstb =deltaR(event.genBHadrons[0].firstb.p4(),event.genBHadrons[1].firstb.p4())
              thisGenPair.hadronPair=event.genBHadrons[0].p4()+event.genBHadrons[1].p4()
              if infBMatch[0]==0 :
                thisGenPair.deltaRForMatching0 = infBMatch[4]
                thisGenPair.deltaRForMatching1 = infBMatch[5]
              else :
                thisGenPair.deltaRForMatching0 = infBMatch[5]
                thisGenPair.deltaRForMatching1 = infBMatch[4]
              event.genBbPairSystem.append(thisGenPair)

          elif len(event.genBPair) == 2 and len(event.genBHadrons) == 4:
            if len(event.selectedSVsSelected) == 2  :
                if self.generatedFromTheSameGluonSplitting(event.genBHadrons, event.selectedSVsSelected) :
                    event.isSignal = True
                event.genBjets=self.clusterizeGenParticles4B(event,event.genBPair[0][0],event.genBPair[0][1],event.genBPair[1][0],event.genBPair[1][1])
                for i in range(0, len(event.genBPair)) :
#                    i = 0
                    thisGenPair=event.genBjets[2*i]+event.genBjets[1+2*i]
                    thisGenPair.deltaRHad=deltaR(event.genBPair[i][0],event.genBPair[i][1])
                    thisGenPair.deltaRJet=deltaR(event.genBjets[2*i],event.genBjets[1+2*i])
                    thisGenPair.deltaRLastb =deltaR(event.genBPair[i][0].lastb,event.genBPair[i][1].lastb)
                    thisGenPair.deltaRFirstb =deltaR(event.genBPair[i][0].firstb.p4(),event.genBPair[i][1].firstb.p4())
                    thisGenPair.hadronPair=event.genBPair[i][0].p4()+event.genBPair[i][1].p4()
                    thisGenPair.deltaRForMatching0 = 20
                    thisGenPair.deltaRForMatching1 = 20
                    event.genBbPairSystem.append(thisGenPair)

#>>>>>>> nuovoHbb

        return True





