from __future__ import division
import os,sys

outpath = os.path.join(os.path.dirname(__file__),'NonPartTrials')
catfl = open(os.path.join('..','OccOutput','NonPartModeledTrials.csv'),'w')
wrthdr = True

for fl in os.listdir(outpath):
    if fl[-4:] != '.csv': continue
    ofl = open(os.path.join(outpath,fl),'rU')
    
    if wrthdr:
        catfl.write(ofl.next())
        wrthdr = False
    else: ofl.next()
    
    for ln in ofl: catfl.write(ln)
    ofl.close()

catfl.close()