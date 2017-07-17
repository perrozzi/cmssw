from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
import PhysicsTools.HeppyCore.framework.config as cfg
from copy import deepcopy
import ROOT
from sets import Set
import os
import resource
import gc
import pprint
import sys

class MemoryLeakAnalyzer( Analyzer ):
   
    def getMemory(self, alternativeMethod = False):
        # in kB
        mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return mem

    def showMemory(self, verbose = False):
        print 'memory:', self.getMemory()/1024, "MB"
        n = gc.collect()
        if verbose:
            print 'Unreachable objects:', n
            newGarbage = gc.garbage
            print 'additional garbage:',len(newGarbage)
            newGarbage.sort(key=lambda x: sys.getsizeof(x))
            totalSize = 0
            for x in newGarbage:
                objSize = sys.getsizeof(x)
                totalSize += objSize
                if objSize > 1024:
                    print ' --> ', objSize/1024, '%r'%x
            print 'total:', totalSize/1024.0, ' kB'
            print 'gc -> memory:', self.getMemory()/1024, "MB"


    def process(self, event):
        try:
            if event.iEv % 100 == 0:
                print "Event number",event.iEv
                self.showMemory() 
        except:
            pass
       
        # measure leak every 1000 events
        try:
            if event.iEv == 0:
                self.lastMemoryUsage = self.getMemory()
            elif event.iEv % 1000 == 0:
                currentMemoryUsage = self.getMemory()
                memoryLeaked = currentMemoryUsage - self.lastMemoryUsage
                # don't printout after first 1000 events since it is unreliable 
                if memoryLeaked > 10 and event.iEv > 1000:
                    print "*"*80
                    print "\x1b[31mMemory leak detected: ", memoryLeaked*0.001, " kB/event!\x1b[0m"
                    print "*"*80
                self.lastMemoryUsage = currentMemoryUsage
        except:
            self.lastMemoryUsage = 0
