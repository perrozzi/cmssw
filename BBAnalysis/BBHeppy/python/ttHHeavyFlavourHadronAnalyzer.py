from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from BBAnalysis.BBHeppy.ttHSVAnalyzer import matchToGenHadron
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class ttHHeavyFlavourHadronAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(ttHHeavyFlavourHadronAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

    def declareHandles(self):
        super(ttHHeavyFlavourHadronAnalyzer, self).declareHandles()
        self.handles['genParticles'] =  AutoHandle( 'prunedGenParticles','std::vector<reco::GenParticle>')

    def beginLoop(self, setup):
        super(ttHHeavyFlavourHadronAnalyzer,self).beginLoop(setup)

    def makeBToDHadrons(self, B, D, event) :
        Dp4 = D.p4()
        Dp4.BDecayPoint = None
        Dp4.DDecayPoint = None
        if B.numberOfDaughters() > 0 :
            Dp4.BDecayPoint = B.daughter(0).vertex()
        if D.numberOfDaughters() > 0 :
            Dp4.DDecayPoint = D.daughter(0).vertex()
        event.genBToDHadrons.append(Dp4)


    def funcionToCheck(self, event) :
        for particle in event.genParticles :
            if abs(particle.pdgId()) != 5 :
                bDaughterNumber = 0
                for dau in particle.daughterRefVector() :
                    if abs(particle.pdgId()) == 5 :
                        bDaughterNumber = bDaughterNumber + 1

                if bDaughterNumber != 0 and bDaughterNumber != 2 :
                    print "particella con non due b"
                    for dau in particle.daughterRefVector() :
                        print particle.pdgId(), "  \t"



    def makeFirstAndLastb(self, event) :

        self.funcionToCheck(event)

        print "Event -------------------------------------------------------------------------------"
        for n in range(0,len(event.genBHadrons)):
            event.genBHadrons[n].bIndex = n
            chainIt_beforSearchForMother = None
            chainIt = event.genBHadrons[n]
            while chainIt != chainIt_beforSearchForMother :
                chainIt_beforSearchForMother = chainIt
                for mom in chainIt.motherRefVector() :
                    if max((abs(mom.pdgId())/1000) % 10, (abs(mom.pdgId())/100) % 10) == 5 :
                        chainIt = mom

            BFromHadronization = True
            for mom in chainIt.motherRefVector() :
                if mom == None or mom.isNull() or not mom.isAvailable(): 
                    print "ERROR: mother of B not available"
                else :
                    if abs(mom.get().pdgId()) == 5 : #and abs(mom.eta()) < 15000 :
                        chainIt = mom
                        BFromHadronization = False
            if BFromHadronization :
                print "ERROR: B is not from b hadronization"
            if chainIt == chainIt_beforSearchForMother :
                print "ERROR: B has not b mother"

            lastb = chainIt

            BDaughterNumber = 0
            for dau in lastb.daughterRefVector() :
                if max((abs(dau.pdgId())/1000) % 10, (abs(dau.pdgId())/100) % 10) == 5 :
                    BDaughterNumber = BDaughterNumber + 1
            if BDaughterNumber == 3 :
                print "---Prima figlia del b con 3 B:                     ", lastb.daughterRefVector()[0].pdgId()


            while chainIt != chainIt_beforSearchForMother :
                chainIt_beforSearchForMother = chainIt
                for mom in chainIt.motherRefVector() :
                    if mom == None or mom.isNull() or not mom.isAvailable(): 
                        print "ERROR: mother of b not available"
                    else :
                        if abs(mom.get().pdgId()) == 5 : #and abs(mom.eta()) < 15000 :
                            chainIt = mom

            firstb = chainIt
            event.genBHadrons[n].firstb = firstb
            firstbP4 = firstb.p4()

            event.genBHadrons[n].gluon = firstb
            gluonOfb = chainIt.motherRefVector()[0] 
            if gluonOfb == None or gluonOfb.isNull() or not gluonOfb.isAvailable(): 
                print "ERROR: gluon not available"
            else :
                event.genBHadrons[n].gluon = gluonOfb
            if len(event.genBHadrons) != 2 :
                print "Number Of daughter of the gluon:     ", len(gluonOfb.daughterRefVector())

#            try :
#            lastb = lastb.get()


#               lastb is a RefVector and firstb in a genParticle------------
            event.genBHadrons[n].lastb = lastb
            lastbP4 = lastb.p4()

            firstbP4.deltaRwithB = deltaR(firstbP4, event.genBHadrons[n].p4())
            lastbP4.deltaRwithB = deltaR(lastbP4, event.genBHadrons[n].p4())
            event.genFirstb.append(firstbP4)
            event.genLastb.append(lastbP4)

#            except :
#                pass


    def makeGenBPair(self, event) :

        gluonList = []
        for Bhad in event.genBHadrons : 
            gluonAlreadyPresent = False
            for g in gluonList : 
                if g == Bhad.gluon :
                    gluonAlreadyPresent = True
            if not gluonAlreadyPresent :
#                if abs(Bhad.firstb.eta()) < 15000 :
                    gluonList.append(Bhad.gluon)

        if len(event.genBHadrons) != 2*len(gluonList) :
            print "---B number:  ", len(event.genBHadrons), "gluon number:  ", len(gluonList)
            if len(event.genBHadrons) == 2 and len(gluonList) == 2 :
                pass
            else :
                print "ERRORE ------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#            print "gluons: "
#            for g in gluonList : 
#                print g


        for g in gluonList : 
            booleanCheck = True
            pairOfB = []
            for Bhad in event.genBHadrons : 
                if g == Bhad.gluon :
                    pairOfB.append(Bhad)
                    booleanCheck = False
            if len(pairOfB) != 2 :
                print "ERROR: not 2 bottom from a gluon. b from the same gluon:  ", len(pairOfB), "  \t number of B:   ", len(event.genBHadrons)
                for Bparticle in pairOfB :
                    print "pdg:  \t quark pdgId: ", Bparticle.firstb.pdgId(), "  \t hadron pdgId: ",  Bparticle.pdgId()
            else :
#                print " 2 bottom from a gluon:  ", len(pairOfB)
                event.genBPair.append(pairOfB)

            if booleanCheck :
                print "Non arriva alla fine!!!!!------------------------------------------------------------------"



    def process(self, event):
        self.readCollections( event.input )
        if not self.cfg_comp.isMC: return True
        event.genParticles = list(self.handles['genParticles'].product() )
       
        def ref2id(ref):
            return (ref.id().processIndex(), ref.id().productIndex(), ref.key()) if ref else (0,0,0)
        def flav(gp):
            id = abs(gp.pdgId())
            return max((id/1000) % 10, (id/100) % 10)
        def same(gp1,gp2):
            return gp1.pdgId() == gp2.pdgId() and gp1.status() == gp2.status() and abs(gp1.pt()-gp2.pt()) < 1e-4  and abs(gp1.eta()-gp2.eta()) < 1e-4  and abs(gp1.phi()-gp2.phi()) < 1e-4
        
        def descendent(child, bhadron):
            mom = child.mother(0) if child.numberOfMothers() > 0 else None
            if mom.status() != 2 or abs(mom.pdgId()) < 100 and abs(mom.pdgId()) != 15:
                return False
            elif same(bhadron,mom):
                return True
            elif mom == None:
                return False
            else:
                return descendent(mom, bhadron)

        # has a memory leak ("mom = mom.motherRef()" ...) , todo: remove
        #def descendent_leaky(child, bhadron):
        #    mom = child.mother(0) if child.numberOfMothers() > 0 else None
        #    while mom != None:
        #        if mom.status() != 2 or abs(mom.pdgId()) < 100 and abs(mom.pdgId()) != 15: break
        #        if same(bhadron,mom):
        #            return True
        #        mom = mom.motherRef() if mom.numberOfMothers() > 0 else None
        #        if mom == None or mom.isNull() or not mom.isAvailable(): break
        #    return False

        event.genLastb = [] 
        event.genFirstb = [] 
        heavyHadrons = []
        event.genAllBHadrons = []
        event.genBToDHadrons = []
        event.genBPair = []

        for g in event.genParticles:
            if g.status() != 2 or abs(g.pdgId()) < 100: continue
            myflav = flav(g)
            if myflav not in [4,5]: continue
            lastInChain = True
            for idau in xrange(g.numberOfDaughters()):
                if flav(g.daughter(idau)) == myflav:
                    lastInChain = False
                    break
            if not lastInChain: continue
            if myflav == 4:
                heaviestInChain = True
                mom = g.motherRef() if g.numberOfMothers() > 0 else None
                while mom != None and mom.isNonnull() and mom.isAvailable():
                    if mom.status() != 2 or abs(mom.pdgId()) < 100: break
                    if flav(mom) == 5:
                        heaviestInChain = False
                        self.makeBToDHadrons(mom, g, event)
                        break
                    mom = mom.motherRef() if mom.numberOfMothers() > 0 else None
                if not heaviestInChain: continue
            # OK, here we are
            g.flav = myflav
            g.firstb = g
            g.lastb = g
            g.bIndex = -1
            heavyHadrons.append(g)
        
        # if none is found, give up here without going through the rest, so we avoid e.g. mc matching for jets
        if len(heavyHadrons) == 0:
            event.genHeavyHadrons = heavyHadrons
            event.genBHadrons = [ h for h in heavyHadrons if h.flav == 5 ]
            event.genDHadrons = [ h for h in heavyHadrons if h.flav == 4 ]
            return True

        # match with IVF 
        had_ivf_pairs = []
        #print "\nNew event"
        for ihad, had in enumerate(heavyHadrons):
            #print "HAD %2d with flav %d,  %d daughters, mass %5.2f, pt %5.2f, eta %+4.2f, phi %+4.2f: " % (ihad, had.flav, had.numberOfDaughters(), had.mass(), had.pt(), had.eta(), had.phi())
            had.sv = None
            for isv,s in enumerate(event.ivf):
                #print "   SV %2d with %d mc tracks, mass %5.2f, pt %5.2f, eta %+4.2f, phi %+4.2f: %d mc-matched tracks" % (isv, s.numberOfDaughters(), s.mass(), s.pt(), s.eta(), s.phi(),len(s.mctracks))
                shared_n, shared_pt = 0, 0 
                for mct in s.mctracks:
                    isDescendent = descendent(mct,had)

                    #DEBUG: test if new descendent function gives same result, todo: remove later
                    #isDescendentLeaky = descendent_leaky(mct,had)
                    #if isDescendent != isDescendentLeaky:
                    #    print "\x1b[31mERROR: no match",mct, " ", had , "\x1b[0m"

                    if isDescendent:
                        shared_n += 1; shared_pt += mct.pt()
                if shared_n:
                    #print "       matched %d tracks (total pt: %.2f) " % (shared_n, shared_pt)
                    had_ivf_pairs.append( (ihad, isv, shared_n, shared_pt) )
        had_ivf_pairs.sort(key = lambda (i1,i2,n,pt) : n + 0.0001*pt, reverse=True)

        for ihad,isv,n,pt in had_ivf_pairs:
            had = heavyHadrons[ihad]
            #print "( had %d, sv %d ): shared %d tracks, %.2f pt ==> %s" % (ihad, isv, n, pt, had.sv)
            if had.sv == None:
                had.sv = event.ivf[isv]
                #print " had %d --> sv %d " % (ihad, isv)
            #else:
            #    print " had %d is already matched " % (ihad,) 
        # match with jets:
        had_jet_pairs = []
        # first loop on jets, get and match daughters
        cccc='''        jetsWithMatchedDaughters = [] 
        for j in event.jetsIdOnly:
            dausWithMatch = []
            for idau in xrange(j.numberOfDaughters()):
                dau = j.daughter(idau)
                if dau.charge() == 0 or abs(dau.eta()) > 2.5: continue
                mct, dr, dpt =  matchToGenHadron(dau, event, minDR=0.05, minDpt=0.1)
                if mct == None: continue
                dausWithMatch.append((dau,mct))
            jetsWithMatchedDaughters.append((j,dausWithMatch))
        for ihad, had in enumerate(heavyHadrons):
            had.jet = None
            for ij,(j,dausWithMatch) in enumerate(jetsWithMatchedDaughters):
                shared_n, shared_pt = 0, 0 
                for dau,mct in dausWithMatch:
                   if descendent(mct,had):
                        shared_n += 1; shared_pt += mct.pt()
                if shared_n:
                    had_jet_pairs.append( (ihad, ij, shared_n, shared_pt) )
        had_jet_pairs.sort(key = lambda (i1,i2,n,pt) : n + 0.0001*pt, reverse=True)
        for ihad,ij,n,pt in had_jet_pairs:
            had = heavyHadrons[ihad]
            if had.jet == None:
                had.jet = event.jetsIdOnly[ij]
'''
        # match with hard scattering
        for had in heavyHadrons:
            had.sourceId = 0
            srcmass = 0
            mom = had.motherRef() if had.numberOfMothers() > 0 else None
            while mom != None and mom.isNonnull() and mom.isAvailable():
                if mom.status() > 2: 
                    if mom.mass() > srcmass:
                        srcmass = mom.mass()
                        had.sourceId = mom.pdgId() 
                    if srcmass > 175:
                        break
                mom = mom.motherRef() if mom.numberOfMothers() > 0 else None
        # sort and save
        heavyHadrons.sort(key = lambda h : h.pt(), reverse=True)
        event.genHeavyHadrons = heavyHadrons
#        event.genAllBHadrons = [ h for h in heavyHadrons if h.flav == 5 ]
#        event.genBHadrons = [ h for h in heavyHadrons if h.flav == 5 and h.pt() > 15.]
        event.genBHadrons = [ h for h in heavyHadrons if h.flav == 5 ]
        event.genAllBHadrons = [ h for h in heavyHadrons if h.flav == 5 and h.pt() > 15.]
        event.genDHadrons = [ h for h in heavyHadrons if h.flav == 4 ]
        if len(event.genBHadrons) > 1 :
            self.makeFirstAndLastb(event)
            self.makeGenBPair(event)
        #print "Summary: "
        #for had in event.genBHadrons:
        #    print "    HAD with %d daughters, mass %5.2f, pt %5.2f, eta %+4.2f, phi %+4.2f: sv %s, jet %s" % (had.numberOfDaughters(), had.mass(), had.pt(), had.eta(), had.phi(), had.sv != None, had.jet != None)
        return True
