#!/usr/bin/env python
import os
import sys
import argparse
import datetime
import subprocess
import glob

def shellRun(command, interactive = False):
    print 'RUN COMMAND: \x1b[32m%s\x1b[0m'%command
    if interactive:
        print '(press any key to run)'
        raw_input()
    ret = subprocess.call([command], shell=True)
    print ' --> RETURN: ',ret
    if interactive:
        print '(press any key to continue with next job)'
        raw_input()

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-S', action='store', dest='sample', help='sample textfile(s) with dataset names, can contain * as wildcard')
parser.add_argument('-q', action='store', dest='queue', help='queue to submit, defaults to all.q', default='all.q')
parser.add_argument('-i', action='store_true', dest='interactive', help='interactive mode')
args = parser.parse_args()
outputFolder = 'root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/berger_p2/gbb/ntuples/2016/'
ntuplesSubfolder = datetime.datetime.now().strftime("%y%m%d_%H%M%S") + '/0000/'

# PREPARE list of jobs
jobs = []
sampleTextFileNames = glob.glob(args.sample)
for sampleTextFileName in sampleTextFileNames:
    with open(sampleTextFileName, 'r') as sampleTextFile:
        sampleFiles = sampleTextFile.readlines()
        sampleFiles.sort()
        for i, file in enumerate(sampleFiles, 1):
            job = {'file': file.strip(), 'number': i}
            jobs.append(job)

# show number of jobs and files
print '-'*80
print ' found %d .txt files'%len(sampleTextFileNames)
print ' -> %d .root files to process'%len(jobs)
print '-'*80
print ' output: ' + outputFolder
print ' subfolder: ' + ntuplesSubfolder
if args.interactive:
    print '(press any key to start submission of jobs)'
    raw_input()

# log path
startOfSubmissionTimestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
logPath = 'log/' + datetime.datetime.now().strftime("%y%m%d") + '/'
try:
    os.makedirs(logPath)
except:
    pass

# SUBMIT jobs
for job in jobs:
    sampleName = job['file'].split('/')[3] + '_' + job['file'].split('/')[4]
    outputFolderSample = outputFolder + sampleName + '/'
    command = 'submit.py -i {input} -o {output} -s {subfolder} -n{number}'.format(input=job['file'], output=outputFolderSample, subfolder=ntuplesSubfolder, number=job['number'])
    print 'COMMAND: \x1b[33m' + command + '\x1b[0m'
    jobName = 'gbbheppy_' + job['file'].split('/')[4] + '_%d'%job['number']
    qsubCommand = 'qsub -V -cwd -q {queue} -N {name} -j y -o {logpath}/{name}_{timestamp}.log {cmd}'.format(queue=args.queue, name=jobName, logpath=logPath, timestamp=startOfSubmissionTimestamp, cmd=command)
    shellRun(qsubCommand, args.interactive)
