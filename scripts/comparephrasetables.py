#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
from collections import defaultdict
import itertools

try:
    phrasetables = sys.argv[1:]
except:
    print >>sys.stderr, "Usage: comparephrasetables.py phrasetable1 phrasetable2"

pairs = []
source = []
target = []
ptsize = []

for phrasetable in phrasetables:
    pairs.append( defaultdict(int) )
    source.append( defaultdict(int) )
    target.append( defaultdict(int) )

    print >>sys.stderr, "Loading " + phrasetable
    for line in codecs.open(phrasetable,'r','utf-8'):
        fields = line.split('|||')
        s = fields[0].strip()
        t = fields[1].strip()
        pairs[-1][(s,t)] += 1
        source[-1][s] += 1
        target[-1][t] += 1
        

    
pairoverlap = 0
pairtotal = 0
sourceoverlap = 0
sourcetotal = 0
targetoverlap = 0
targettotal = 0
    

print >>sys.stderr, "Computing pair overlap"
overlapmem = {}    
for x in itertools.chain( *( x.keys() for x in pairs ) ):    
    if not x in overlapmem:
        overlapmem[x] = True        
        pairtotal += 1
        overlap = True
        for i, phrasetable in enumerate(phrasetables): 
            if not (x in pairs[i]):
                overlap = False
                break
        if overlap:      
            s =  "PAIR MATCH: " +  x[0]  +  " <---> " +  x[1]
            print s.encode('utf-8')
            pairoverlap += 1
            
            
print >>sys.stderr, "Computing source overlap"                    
overlapmem = {}
for x in itertools.chain( *( x.keys() for x in source ) ):    
    if not x in overlapmem:
        overlapmem[x] = True        
        sourcetotal += 1
        overlap = True
        for i, phrasetable in enumerate(phrasetables): 
            if not (x in source[i]):
                overlap = False
                break
        if overlap:                
            s =  "SOURCE MATCH: " +  x
            print s.encode('utf-8')
            sourceoverlap += 1
            
            
                    

print >>sys.stderr, "Computing target overlap"
overlapmem = {}
for x in itertools.chain( *( x.keys() for x in target ) ):    
    if not x in overlapmem:
        overlapmem[x] = True        
        targettotal += 1
        overlap = True
        for i, phrasetable in enumerate(phrasetables): 
            if not (x in target[i]):
                overlap = False
                break
        if overlap:      
            s  = "TARGET MATCH: " +  x
            print s.encode('utf-8')
            targetoverlap += 1
            

print >>sys.stderr, "Total combined pairs:\t" , pairtotal
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "Total pairs in " + phrasetable + ":\t", len(pairs[i]), len(pairs[i]) / float(pairtotal)
print >>sys.stderr, "Overlapping pairs:\t" , pairoverlap,  pairoverlap / float(pairtotal)
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "   ...as fraction of " + phrasetable + ":\t",pairoverlap / float(len(pairs[i]))

print >>sys.stderr, "Total combined source patterns: " , sourcetotal
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "Total source patterns in " + phrasetable + ": ", len(source[i]), len(source[i]) / float(sourcetotal)
print >>sys.stderr, "Overlapping source patterns: " , sourceoverlap,  sourceoverlap / float(sourcetotal)
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "   ...as fraction of " + phrasetable + ":\t",sourceoverlap / float(len(source[i]))

print >>sys.stderr, "Total combined target patterns: " , targettotal
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "Total target patterns in " + phrasetable + ": ", len(target[i]), len(target[i]) / float(targettotal)
print >>sys.stderr, "Overlapping target patterns: " , targetoverlap,  targetoverlap / float(targettotal)
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "   ...as fraction of " + phrasetable + ":\t",targetoverlap / float(len(target[i]))
    
    
