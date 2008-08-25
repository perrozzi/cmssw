#!/usr/bin/env python

import time, os, sys
import optparse as opt
#from ROOT import gROOT, TCanvas, TF1
import ROOT

def getParameters():
    parser = opt.OptionParser()
    #
    # Options
    #
    parser.add_option('-n',
                      type="int",
                      help='Number of secs per bin. Default is 1.' ,
                      default=1,
                      dest='startevt')     
    (options,args) = parser.parse_args()
    if not len(args) == 2:
        print "ERROR: Not enough arguments"
        sys.exit()

    path1 = os.path.abspath(args[0])
    path2 = os.path.abspath(args[1])    
    if os.path.exists(path1) and os.path.exists(path2):
        return (path1, path2, options.startevt)
    else:
        print "Error: one of the paths does not exist"
        sys.exit()

def get_max(data,index=1):
    max_time=-1
    for el in data:
        sec=el[index]
        if max_time<sec:
            max_time=sec
    return max_time

def get_min(data,index=1):
    min_time=1e20
    for el in data:
        sec=el[index]
        if min_time>sec:
            min_time=sec
    return min_time  

def createROOT(outdir,filename):
    __argv=sys.argv # trick for a strange behaviour of the TApp..
    sys.argv=sys.argv[:1]
    ROOT.gROOT.SetStyle("Plain") # style paranoia
    sys.argv=__argv
    #Cannot use this option when the logfile includes
    #a large number of events... PyRoot seg-faults.
    #Set ROOT in batch mode to avoid canvases popping up!
    #ROOT.gROOT.SetBatch(1)

    # Save in file
    rootfilename='%s/%s' %(outdir,filename)
    myfile=ROOT.TFile(rootfilename,'RECREATE')
    return myfile

def getDataFromTimingLog(logfile_name):
    data=[]
    
    # open file and read it and fill the structure!
    logfile=open(logfile_name,'r')
    logfile_lines=logfile.readlines()
    logfile.close()

    # we get the info we need!
    i=0
    while i < len(logfile_lines):
        line=logfile_lines[i]
        if 'TimeEvent>' in line:
            line=line.strip()
            line_content_list = line.split(' ')[0:]
            event_number = int(line_content_list[1])
            seconds = float(line_content_list[3])
            data.append((event_number,seconds))
        i+=1
        
    return data

def newGraphAndHisto(npoints,nbins,min_val,max_val,data,graph_num):
    colors = [2,4]

    histo=ROOT.TH1F('Seconds per event','Seconds per event',nbins,min_val,max_val)

    graph=ROOT.TGraph(npoints)
    evt_counter=0
    total=0
    for evt_num,secs in data:
        graph.SetPoint(evt_counter,evt_num,secs)
        histo.Fill(secs)
        total+=secs
        evt_counter+=1
        
    print 'Total Time=', total
    
    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(.7)
    graph.SetMarkerColor(1)
    graph.SetLineWidth(3)
    graph.SetLineColor(colors[graph_num]) # for each iterate through colors
    histo.SetLineColor(colors[graph_num])

    return (graph,histo)

def getLimits(data,secsperbin):
    min_val=get_min(data,1)
    max_val=get_max(data,1)
    interval=int(max_val-min_val)
    
    min_val=min_val-interval*0.2
    max_val=max_val+interval*0.2
    interval=int(max_val-min_val)
    
    nbins=int(interval/secsperbin)

    npoints=len(data)

    last_event=data[-1][0]
    print 'last event =',last_event    

    return (min_val,max_val,interval,npoints,last_event)

def setupSuperimpose(graph1,graph2,last_event,max_val):
    graph1.SetTitle('Seconds per event')
    graph1.SetName('SecondsPerEvent')    
    graph1.GetXaxis().SetTitle("Event")
    graph1.GetYaxis().SetTitleOffset(1.3)
    graph1.GetYaxis().SetTitle("s")
    graph1.GetXaxis().SetLimits(0,last_event)
    graph1.GetYaxis().SetRangeUser(0,max_val)
    # Do I need to set limits on graph 2? I don't know
    # I'm doing it anyway, can't hurt.
    graph2.GetXaxis().SetLimits(0,last_event)    
    graph2.GetYaxis().SetRangeUser(0,max_val)    

def getMeanLines(histo,last_event,graph_num):
    colors = [2,4]
    avg=histo.GetMean()
    avg_line=ROOT.TLine(1,avg,last_event,avg)
    avg_line.SetLineColor(colors[graph_num])
    avg_line.SetLineWidth(2)

    return avg_line

def getDiff(data1,data2,npoints,last_event,orig_max_val):
    data3 = []
    for x in range(len(data2)):
        try:
            #avgEventNum = (data2[x][0] + data1[x][0]) / 2
            if data2[x][0] == data1[x][0]:
                avgEventNum = data2[x][0]
                diffSecs    = data2[x][1] - data1[x][1]                
                data3.append((avgEventNum,diffSecs))                
        except IndexError:
            pass
        except ValueError:
            pass

    graph=ROOT.TGraph(npoints)

    evt_counter=0
    total=0
    for evt_num,secs in data3:
        graph.SetPoint(evt_counter,evt_num,secs)
        total+=secs
        evt_counter+=1
        
    min_val=get_min(data3,1)
    max_val=get_max(data3,1)
    interval=int(max_val-min_val)
    
    min_val=min_val-interval*0.2
    max_val=max_val+interval*0.2
    interval=int(max_val-min_val)
    print min_val
    # Determine the max value to be something that makes the scale similar to what
    # the original graph had. Unless we can't seem the maximum value.

    new_max = min_val + orig_max_val
    print "new max", new_max, "min_val", min_val, "max_val", max_val
    if new_max < max_val:
        pass
    else :
        max_val = new_max
    
    graph.SetTitle('Change in Seconds per event')
    graph.SetName('SecondsPerEvent')    
    graph.GetXaxis().SetTitle("Change in time per Event between revs")
    graph.GetYaxis().SetTitleOffset(1.3)
    graph.GetYaxis().SetTitle("s")
    graph.GetXaxis().SetLimits(0,last_event)
    graph.GetYaxis().SetRangeUser(min_val,max_val)

    return graph

def drawGraphs(graph1,graph2,avg1,avg2):
    graph_canvas = ROOT.TCanvas("graph_canvas")
    graph_canvas.cd()
    graph1.Draw("ALP")
    graph2.Draw("LP")
    avg1.Draw("Same")
    avg2.Draw("Same")
    return graph_canvas

def drawHistos(histo1,histo2):
    histo_canvas = ROOT.TCanvas("histo_canvas")
    histo_canvas.cd()
    histo1.GetXaxis().SetTitle("s")        
    histo1.Draw("")
    histo2.Draw("")
    return histo_canvas

def drawChanges(graph):
    graph_canvas = ROOT.TCanvas("change_canvas")
    graph_canvas.cd()
    graph.Draw("ALP")
    return graph_canvas

def main():
    rootfilename = "regression.root"
    outdir = os.getcwd()
    
    (file1,file2,secsperbin)  = getParameters()

    data1 = getDataFromTimingLog(file1)
    data2 = getDataFromTimingLog(file2)

    newrootfile = createROOT(outdir,rootfilename)

    (min_val1,max_val1,nbins1,npoints1,last_event1) = getLimits(data1,secsperbin)
    (min_val2,max_val2,nbins2,npoints2,last_event2) = getLimits(data2,secsperbin)

    (graph1,histo1) = newGraphAndHisto(npoints1,nbins1,min_val1,max_val1,data1,0)
    (graph2,histo2) = newGraphAndHisto(npoints2,nbins2,min_val2,max_val2,data2,1)

    biggestLastEvt = last_event1
    biggestMaxval  = max_val1
    
    if last_event2 > biggestLastEvt:
        biggestLastEvt = last_event2
    if max_val2 > biggestMaxval:
        biggestMaxval = max_val2
        
    changegraph = getDiff(data1,data2,npoints2,biggestLastEvt,biggestMaxval)
    setupSuperimpose(graph1,graph2,biggestLastEvt,biggestMaxval)
    avg_line1 = getMeanLines(histo1,last_event1,0)
    avg_line2 = getMeanLines(histo2,last_event2,1)

    graph_canvas   = drawGraphs(graph1,graph2,avg_line1,avg_line2)
    changes_canvas = drawChanges(changegraph)
    histo_canvas   = drawHistos(histo1,histo2)

    
    #
    # Create a one dimensional function and draw it
    #

    while graph_canvas or histo_canvas or changes_canvas:
        time.sleep(2.5)

if __name__ == "__main__":
    main()

