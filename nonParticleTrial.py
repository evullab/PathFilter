# Table to extend PhysicsTable that adds noise on movement and bounces

from __future__ import division
import os, sys, random, time
import pygame as pg
import numpy as np

# For importing from shared
sharepath = os.path.join(os.path.dirname(__file__),'..','SharedCode')
predpath = os.path.join(os.path.dirname(__file__),'..','PredictionModel')
sys.path.insert(1,sharepath)
sys.path.insert(1,predpath)
from dynamicTable import *
from dynamicTrial import *
from sampleTrial import sampletrial
from uncertaintyParameters import *
from noisyTable import *
from EasyMultithread import *

outpath = os.path.join(os.path.dirname(__file__),'NonPartTrials')
trialpath = os.path.join(os.path.dirname(__file__),'..','OccTrials')

random.seed(1010101)

def simulateTrial(trial, oflnm, samplerate = 10., kapv = kapv_def, kapb = kapb_def, kapm = kapm_def, perr = perr_def, nsims = 500, maxsteps = 10000, notify = False):
    
    timepersample = 1. / samplerate
    ofl = open(os.path.join(outpath,oflnm), 'w')
    ofl.write('Trial,Time,Px,Py,Vx,Vy,Red,Green,Uncertain,KapV,KapB,KapM,Perr\n')
    
    ender = ',' + str(kapv) + ',' + str(kapb) + ',' + str(kapm) + ',' + str(perr) + '\n'
    
    nm = trial.name
    i = 0
    
    table = trial.makeTable()
    sims = predict(table,kapv,kapb,kapm,perr,nsims,maxsteps,return_path = False)
    time = i * timepersample
    oldsims = sims
    
    ofl.write(nm + ',' + str(time) + ',' + str(table.ball.x) + ',' + str(table.ball.y) + ',' + str(table.ball.v[0]) + ',')
    ofl.write(str(table.ball.v[1]) + ',' + str(sims['Red']) + ',' + str(sims['Green']) + ',' + str(sims['Uncertain']) + ender)
    
    running = True
    
    while running:
        i += 1
        e = table.step(t = timepersample)
        if notify: print i
        
        if table.fullyOcc(): sims = oldsims
        else: sims = predict(table,kapv,kapb,kapm,perr,nsims,maxsteps,return_path = False); oldsims = sims
        
        time = i * timepersample
    
        ofl.write(nm + ',' + str(time) + ',' + str(table.ball.x) + ',' + str(table.ball.y) + ',' + str(table.ball.v[0]) + ',')
        ofl.write(str(table.ball.v[1]) + ',' + str(sims['Red']) + ',' + str(sims['Green']) + ',' + str(sims['Uncertain']) + ender)
        
        if e is not True: running = False

    ofl.close()


def simulateByName(trname, yell = False):
    
    trpath = os.path.join(trialpath,trname+'.ptr')
    outname = 'Model_' + trname + '.csv'
    
    trial = loadTrial(trpath)
    
    simulateTrial(trial,outname, notify = yell)
    
    print "Done with trial:", trname
    
if __name__ == '__main__':
    #stime = time.time()
    #simulateByName('RandTrial_0')
    #print 'Done; t =', int(time.time()-stime)
    
    trialdir = os.listdir(trialpath)
    trs = []
    
    for t in trialdir:
        splt = t.split('.')
        if len(splt) == 2:
            if splt[1] == 'ptr': trs.append(splt[0])
    
    multimap(simulateByName,trs)
    