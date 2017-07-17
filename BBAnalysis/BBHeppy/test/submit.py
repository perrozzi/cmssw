#!/usr/bin/env python
import os
import PhysicsTools.HeppyCore.framework.config as cfg
#cfg.Analyzer.nosubdir=True

import sys
import re
import imp
import hashlib
import subprocess
import argparse

# optional modules
try:
    import pwd
except:
    pass

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='input',
                    help='input root file')
parser.add_argument('-o', action='store', dest='output',
                    help='output folder')
parser.add_argument('-s', action='store', dest='subfolder',
                    help='subfolder')
parser.add_argument('-n', action='store', dest='number',
                    help='tree number')
args = parser.parse_args()

inputFileName = args.input #'/store/mc/RunIIFall15MiniAODv2/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0023A3AF-8FB8-E511-85EF-0025905AC99A.root'
ntuplesFolder = args.output #'root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/berger_p2/gbb/ntuples/2015/'
ntuplesSubfolder = args.subfolder #'170612_141100/0000/'
treeNumber = int(args.number)
userName = None # None = automatic

# try to get username if not specified
if not userName:
    try:
        userName = pwd.getpwuid(os.getuid()).pw_name
    except:
        print "userName could not be determined automatically, it has to be set manually in the file"
        exit(-1)

def getXRDname(fileName):
    #return 'root://cms-xrd-global.cern.ch/' + fileName
    return 'root://xrootd-cms.infn.it/' + fileName

def shellRun(command):
    print 'RUN COMMAND: \x1b[32m%s\x1b[0m'%command
    ret = subprocess.call([command], shell=True)
    print 'RETURN: ',ret

def mkdirRemote(folder):
    print "folder =", folder
    folderPnfsPath = '/pnfs/' + folder.strip('/').split('/pnfs/')[1]
    print "folderPnfsPath = ", folderPnfsPath
    xfPrefix = 'xrdfs t3dcachedb03.psi.ch '
    folder3 = '/'.join(folderPnfsPath.split('/')[0:-3])
    if not os.path.exists(folder3):
        command = xfPrefix + 'mkdir ' + folder3
        shellRun(command)
    folder3 = '/'.join(folderPnfsPath.split('/')[0:-2])
    if not os.path.exists(folder3):
        command = xfPrefix + 'mkdir ' + folder3
        shellRun(command)
    folder3 = '/'.join(folderPnfsPath.split('/')[0:-1])
    if not os.path.exists(folder3):
        command = xfPrefix + 'mkdir ' + folder3
        shellRun(command)
    if not os.path.exists(folderPnfsPath):
        command = xfPrefix + 'mkdir ' + folderPnfsPath
        shellRun(command)

def ntuplesExist(path):
    path2 = path
    if '/pnfs/' in path2:
        path2 = '/pnfs/' + path2.strip('/').split('/pnfs/')[1]
    return os.path.isfile(path2)

# prepare output folders
ntuplesFullPath = ntuplesFolder + ntuplesSubfolder 
mkdirRemote(ntuplesFullPath)
outputFileNameTree = ntuplesFullPath + '/tree_%d.root'%treeNumber

# check if file exists
if ntuplesExist(outputFileNameTree):
    print 'already existing! => SKIP'
    print outputFileNameTree
    exit()

# load gBB config
handle = open("bb.py", 'r')
cfo = imp.load_source("bb", "bb.py", handle)
config = cfo.config
handle.close()

# input files
config.components[0].files = [getXRDname(inputFileName)]
config.components[0].isEmbed = False

# prepare tmp/output folder
replacePatterns=[
("pythia","Py"),
("76X_mcRun2_asymptotic","76r2as"),
("RunIIFall15MiniAODv2","fall15MAv2"),
("PU25nsData2015","pu25ns15")
]
outputName = inputFileName.strip('/').split('/')[3] + '_' + inputFileName.strip('/').split('/')[5]
for replacePattern in replacePatterns:
    outputName = outputName.replace(replacePattern[0],replacePattern[1])
print outputName
tmpFolder = '/scratch/%s/heppy/%s/%d/output'%(userName, outputName, treeNumber)
try:
    os.makedirs(tmpFolder)
except:
    pass

# run loop on events
from PhysicsTools.HeppyCore.framework.looper import Looper
looper = Looper( tmpFolder, config, nPrint = 1)
looper.loop()
looper.write()

print "looper output:", looper.name

# copy to final storage
command = 'xrdcp -d 1 '+ looper.name +'/tree.root '+ outputFileNameTree
shellRun(command)

command = 'rm ' + looper.name + '/tree.root'
shellRun(command)
command = 'rm ' + looper.name + '/cmsswPreProcessing.root'
shellRun(command)
command = 'rm ' + looper.name + '/*.log'
shellRun(command)

deleteFolder = looper.name.strip()
if deleteFolder.startswith('/scratch/%s/'%userName) and ' ' not in deleteFolder:
    command = 'rm -rf ' + deleteFolder
    shellRun(command)




