from __future__ import division
import os, sys, random, copy, bisect, time
import pygame as pg
import numpy as np
import scipy as sp
import shutil

# For importing from shared
sharepath = os.path.join(os.path.dirname(__file__),'..','SharedCode')
predpath = os.path.join(os.path.dirname(__file__),'..','PredictionModel')
trpath = os.path.join(os.path.dirname(__file__),'..','OccTrials')
sys.path.insert(1,sharepath)
sys.path.insert(1,predpath)
from dynamicTable import *
from dynamicTrial import *
from sampleTrial import sampletrial
from uncertaintyParameters import *
from noisyTable import *
from rectangles import *
from pyText import *
from particleModel import *
from EasyMultithread import *

def makedec(x):
    if x == RED: return 'R'
    elif x == GREEN: return 'G'
    else: return 'U'

# Simulate a single set of particles
def singleSim(tr, pathid, numpart, pcut, declim = 10., temp = 1., ofl = None, kapv = kapv_def, kapb = kapb_def, kapm = kapm_def, perr = perr_def, timeperstep = .1):
    
    partTab = ParticleTable(tr.makeTable(), nparticles = numpart, timeperstep = timeperstep, newp = pcut, tmp = temp)
    ps = partTab.getPartPs()
    decs = []
    t = 0.
    
    if ofl:
        poss = partTab.getPartPos()
        ds = map(lambda p: makedec(p.getdecision(partTab.t, declim)),partTab.particles)
        ofl.write(tr.name + ',' + str(pathid) + ',' + str(t) + ',' + partTab.getDecision(declim,ps))
        for pos, p, d in zip(poss, ps, ds):
            ofl.write(',' + str(pos[0]) + ',' + str(pos[1]) + ',' + str(p) + ',' + d)
        ofl.write('\n')
    
    
    
    running = True
    while running:
        s = partTab.step()
        t += timeperstep
        if s is not True: running = False
        ps = partTab.getPartPs()
        dec = partTab.getDecision(declim, ps)
        decs.append(dec)
        if ofl:
            poss = partTab.getPartPos()
            ds = map(lambda p: makedec(p.getdecision(partTab.t, declim)),partTab.particles)
            ofl.write(tr.name + ',' + str(pathid) + ',' + str(t) + ',' + dec)
            for pos, p, d in zip(poss, ps, ds):
                ofl.write(',' + str(pos[0]) + ',' + str(pos[1]) + ',' + str(p) + ',' + d)
            ofl.write('\n')
    
    return decs

def trialSim(tr, numpart, pcut, declim = 10., temp =1., nsims = 50, path = '.', kapv = kapv_def, kapb = kapb_def, kapm = kapm_def, perr = perr_def, timeperstep = .1):
    
    #tr = loadTrial(os.path.join('..','OccTrials',trnm+'.ptr'))
    trnm = tr.name
    
    switms = dict()
    
    flnm = 'Part{:02d}_Prob{:1.1e}_TLim{:02d}_Temp{:0.2f}_{}.csv'.format(numpart,pcut,int(declim),temp,tr.name)
    fl = open(os.path.join(path,flnm),'w')
    fl.write('Trial,SimNo,Time,Decision')
    for i in range(numpart):
        fl.write(',Part{0}_X,Part{0}_Y,Part{0}_p,Part{0}_choice'.format(i))
    fl.write('\n')
    
    firsttime = True
    
    for i in range(nsims):
        try:
            ds = singleSim(tr,i,numpart,pcut,declim,temp,fl,kapv,kapb,kapm,perr,timeperstep)
        except:
            print "Error in", trnm
            raise
        nsw = 1
        if firsttime:
            reds = [1*(x=='R') for x in ds]
            greens = [1*(x=='G') for x in ds]
            switch = [0]
            switch.extend([1*(ds[k] != ds[k-1]) for k in range(1,len(ds))])
            #print switch
            firsttime = False
            for j in range(1,len(ds)):
                if ds[j] != ds[j-1]:
                    try:
                        switms[nsw].append((j+1) * timeperstep)
                    except:
                        switms[nsw] = [(j+1)*timeperstep]
                    nsw += 1
                
        else:
            for j in range(len(ds)):
                if ds[j] == 'R': reds[j] += 1
                if ds[j] == 'G': greens[j] += 1
                if j != 0:
                    if ds[j] != ds[j-1]:
                        switch[j] += 1
                        try:
                            switms[nsw].append((j+1) * timeperstep)
                        except:
                            switms[nsw] = [(j+1)*timeperstep]
                        nsw += 1
                    
        #print i
    fl.close()
    reds = map(lambda x: x / float(nsims), reds)
    greens = map(lambda x: x / float(nsims), greens)

    print "Done:",trnm
    return [trnm,reds,greens, switch, switms]
    
def simulateThemAll(numpart, pcut, declim, temp=1., nsims = 50, kapv = kapv_def, kapb = kapb_def, kapm = kapm_def, perr = perr_def, timeperstep = .1):
    
    # Make a new folder for these (first check if one exists though)
    folnm = os.path.join('.','ParticleOutput',"OccPaths_Part{:02d}_Prob{:1.1e}_TLim{:02d}_Temp{:0.2f}".format(numpart,pcut,int(declim),temp))
   
    try:
        os.mkdir(folnm)
    except OSError:
        #print 'Directory already exists'
        asking = True
        while asking:
            ans = raw_input('Directory already exists - delete? (y/n): ')
            if ans == 'n': return None
            if ans == 'y': asking = False
        shutil.rmtree(folnm)
        os.mkdir(folnm)
    
    trials = []
    trpath = os.path.join(os.path.dirname(__file__),'..','OccTrials')
    for f in os.listdir(trpath):
        if f[-4:] == '.ptr': trials.append(loadTrial(os.path.join(trpath,f)))
    
    # For testing - comment out when done
    #trials = [loadTrial(os.path.join(trpath,'OccRandTrial_'+str(t)+'.ptr')) for t in range(4)]
    
    # Run each trial
    outs = multimap(lambda t: trialSim(t,numpart,pcut,declim,temp,nsims,folnm, kapv,kapb,kapm,perr, timeperstep), trials)
    
    # Summarize output
    sumfl = open(os.path.join(folnm,'Summary_Part{:02d}_Prob{:1.1e}_TLim{:02d}_Temp{:0.2f}.csv'.format(numpart,pcut,int(declim),temp)),'w')
    sumfl.write('Trial,T,Red,Green,Switches')
    for i in range(1,9):
        sumfl.write(',SwitchT' + str(i))
    sumfl.write('\n')
    for o in outs:
        tnm = o[0]
        reds = o[1]
        greens = o[2]
        switch = o[3]
        switms = o[4]
        clsw = map(lambda x: map(lambda y: round(y,1),switms[x]),range(1,min(len(switms)+1,9)))
        tm = 0.
        for i in range(len(reds)):
            tm += .1
            sumfl.write(tnm+','+str(tm)+','+str(reds[i])+','+str(greens[i])+','+str(switch[i]))
            for k in range(0,8):
                if k >= len(clsw): sumfl.write(',0')
                else:
                    thisct = sum(1 for x in clsw[k] if x == round(tm,1))
                    sumfl.write(',' + str(thisct))
            sumfl.write('\n')
        print 'Wrote',tnm
    sumfl.close()
    return True
    
    
if __name__ == '__main__':
    
    npart = int(sys.argv[1])
    pcut = float(sys.argv[2])
    tlim = float(sys.argv[3])
    temp = float(sys.argv[4])
    
    simulateThemAll(npart,pcut,tlim,temp)
    
        