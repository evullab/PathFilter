from __future__ import division
import os, sys, random, copy, bisect
import pygame as pg
import numpy as np
import scipy as sp

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
from rectangles import *
from pyText import *
from mvncdf import mvnormcdf

# Likelihood of needing a new particle (will need to play around with this...)
NEWP = 10e-6
TEMP = 1

def euclidist(p1, p2): return np.sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 )

def selectReplace(items, weights, n):
    if len(items) != len(weights): raise(Exception("Needs equal number of items, weights"))
    cumw = [sum(weights[:(x+1)]) for x in range(len(weights))]
    return [items[bisect.bisect(cumw,random.random()*max(cumw))] for x in range(n)] 


def normpdf(x, mu=0, sigma=1):
    u = (x-mu)/abs(sigma)
    y = (1/(np.sqrt(2*np.pi)*abs(sigma)))*np.exp(-u*u/2)
    return y

class Particle(object):
    def __init__(self, parenttable, kv, kb, km, pe, timeperstep = .1, inittime = 0., spos = None, toff = 0):
        if spos is None: self.pos = parenttable.ball.getpos()
        else: self.pos = spos
        self.v = parenttable.ball.v
        self.tab = makeNoisy(parenttable, kv,kb,km,pe)
        self.inittime = inittime
        self.tps = timeperstep
        self.weight = 1
        
        self.setpath(self.tab)
    
    def setpath(self, noisytable, limitsteps = 50000):
        ntab = copy.deepcopy(noisytable)
        # Re-add the ball with noise
        b = ntab.addBall(self.pos,self.v)
        
        capturesteps = int(1000. * self.tps)
        
        self.path = [b.getpos()]
        
        for i in range(1,limitsteps):
            r = ntab.physicsstep()
            if i % capturesteps: self.path.append(b.getpos())
            if r == GREEN: self.end = GREEN; self.timedec = self.inittime + (i / 1000.); return self.path
            if r == RED: self.end = RED; self.timedec = self.inittime + (i / 1000.); return self.path
        
        self.end = -1
        self.timedec = self.inittime + (limitsteps / 1000.)
        return self.path
    
    def getpos(self, t):
        if t < self.inittime: raise Exception("Cannot call position from before simulation starts")
        idx = int( (t - self.inittime) * 1000)
        if idx > len(self.path): return self.path[-1]
        return self.path[idx]
    
    def getdecision(self, curtime, tlimit = 10.):
        if self.timedec > (curtime + tlimit): return -1
        return self.end
    
    def getPath(self, t, tlimit = 10.):
        begidx = int( (t - self.inittime) * 1000)
        endidx = int( (t - self.inittime + tlimit) * 1000)
        if endidx > len(self.path): endidx = len(self.path)
        return self.path[begidx:endidx]
    
    def getP(self, parenttable):
        pos = self.getpos(parenttable.t)
        if not parenttable.table.fullyOcc():
            dist = euclidist(pos, parenttable.table.ball.getpos())
            #print dist, parenttable.perr, normpdf(dist,0,parenttable.perr)
            return normpdf(dist,0,parenttable.perr)
        else:
            # If the ball is covered, the particle should test vs probability that it should be covered
            # NOTE: MUST BE A FASTER WAY OF DOING THIS - use pg.Rect.clip() to get intersection
            covmat = np.array([[parenttable.perr,0],[0,parenttable.perr]])
            brad = int(parenttable.table.rad)
            blocked = []
            for w in parenttable.table.walls:
                 blocked.append(w.r.inflate(2*brad,2*brad))
                 
            p = 0
            ps = []
            
            # Break up occluders into non-overlapping, then run
            uoccs = uniqueOccs(map(lambda x: x.r, parenttable.table.occludes), blocked)
            #print map(lambda x: x.r, parenttable.table.occludes)
            #print blocked
            #print uoccs
            for o in uoccs:
                lefto = o.left
                righto = o.right
                topo = o.top
                boto = o.bottom
                ''' OLD WAY: INDIVIDUAL PARAMETER CALC
                for x in range(lefto, righto):
                    for y in range(topo, boto):
                        pt = (x,y)
                        good = True
                        for b in blocked:
                            if b.r.collidepoint(pt): good = False; continue
                        if good:
                            dist = euclidist(pos, pt)
                            p += normpdf(dist,0,parenttable.perr)
                blocked.append(o)
                '''
                tst = mvnormcdf([lefto,topo],[righto,boto],pos,covmat)
                #if tst < 0: print o, pos, tst
                p += mvnormcdf([lefto,topo],[righto,boto],pos,covmat)

                ps.append(p)
                
            #if p < 0: print p, ps, uoccs
            return p
                
        
        

class ParticleTable(object):
    
    # Calls basic __init__ but also jitters position
    def __init__(self, phystable, kapv = kapv_def, kapb = kapb_def, kapm = kapm_def, perr = perr_def, nparticles = 5, timeperstep = .1):
        self.table = phystable
        self.kapv = kapv
        self.kapb = kapb
        self.kapm = kapm
        self.perr = perr
        self.npart = nparticles
        self.tps = timeperstep
        self.particles = [Particle(self.table, kapv, kapb, kapm, perr, timeperstep) for i in range(nparticles)]
        self.t = 0
        self.lastseen = phystable.ball.getpos()
        self.lastseent = 0
        self.getPartPs()
        
    def getPartPs(self):
        # Calculate P(obs | particlepos) given Gaussian error with sigma = perr
        return map(lambda p: p.getP(self), self.particles)
        
    def step(self):
        r = self.table.step(self.tps)
        self.t += self.tps
        if not self.table.fullyOcc():
            self.lastseent = self.t
            self.lastseen = map(int,self.table.ball.getpos())
        print self.lastseen, self.lastseent
        
        # Update particles
        weights = [p.weight for p in self.particles]
        ps = self.getPartPs()
        newws = [w*p for w,p in zip(weights, ps)]
        newws.append(NEWP)
        newws = map(lambda x: np.power(x,TEMP), newws)
        totw = sum(newws)
        newws = map(lambda x: x / totw, newws)
        #seff = sum(map(lambda w: 1 / (w*w), newws))
        newparts = copy.copy(self.particles); newparts.append("Empty")
        newps = selectReplace(newparts,newws,len(self.particles))
        for i in range(len(newps)):
            if newps[i] == "Empty": newps[i] = Particle(self.table,self.kapv,self.kapb,self.kapm,self.perr,self.tps,self.lastseent, self.lastseen)
            #else: newps[i] = copy.deepcopy(newps[i])
        for p in newps: p.weight = 1
        self.particles = newps
        return r
        
    def draw(self, drawlines = False, stillshow = False):
        sc = self.table.draw(stillshow = stillshow)
        
        # Draw on the particles
        for part, p in zip(self.particles, self.getPartPs()):
            # Draw particles proportional to confidence
            col = part.getdecision(self.t)
            if col == -1: col = GREY
            try:
                rad = 15 + int(np.log(p))
            except:
                print self.table.ball.getpos(), self.getPartPs()
                sys.exit(0)
            if rad < 2: rad = 2
            pg.draw.circle(sc, col, map(int,part.getpos(self.t)),rad)
            if drawlines:
                for pt in part.getPath(self.t):
                    sc.set_at(map(int,pt), col)
        return sc


p = os.path.join('..','OccTrials','Tr1_Tube.ptr')
sampocc = loadTrial(p)
pt = ParticleTable(sampocc.makeTable())
blocked = []
for w in pt.table.walls:
    blocked.append(w.r.inflate(40,40))
occs = map(lambda x: x.r, pt.table.occludes)

def simpleDraw(sc, occs, walls):
    sc.fill(WHITE)
    for o in occs:
        pg.draw.rect(sc,GREY,o)
    for w in walls:
        pg.draw.rect(sc,BLACK,w)
    pg.display.flip()
                
if __name__ == "__main__":
    pg.init()
    clock = pg.time.Clock()
    sc = pg.display.set_mode((1200,900))
    
    if len(sys.argv) == 1:
        sampocc = copy.deepcopy(sampletrial)
        sampocc.addObsc(pg.Rect(400,0,200,600))
    else:
        trnm = sys.argv[1]
        try:
            p = os.path.join('..','OccTrials',trnm+'.ptr')
            sampocc = loadTrial(p) 
        except:
            print "Trial name not found"
            sys.exit(0)
        
        
    pt = ParticleTable(sampocc.makeTable())
    tscreen = pt.draw(True)
    sc.blit(tscreen,(0,0))
    pg.display.flip()
    pause(waittime = .1)
    running = True
    while running:
        clock.tick(1/pt.tps)
        s = pt.step()
        if s is not True: running = False
        tscreen = pt.draw(True, stillshow = True)
        sc.blit(tscreen,(0,0))
        pg.display.flip()
        #print pt.table.ball.getpos(), pt.table.fullyOcc()
        #print [p.pos for p in pt.particles]
        #print pt.getPartPs()
        #print ''
        pause(waittime = .1)
    pause()
    pg.quit()
