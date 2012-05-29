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
        

# def shorten(item, items, pivot = 0):    
    # if isinstance(item, tuple):
        #TODO: all permutations
    # else:
        # words = item.split(' ')
        # l = 
        # while l >= 1:
            # words
        # yie
        
        

def computeoverlap(items, label= 'PAIR', pivot = 0):    
    total = 0
    overlap = 0
    overlapmem = {}    
    for item in itertools.chain( *( x.keys() for x in items ) ):    
        if not item in overlapmem:
            total += 1
            overlap = True
            for i, _ in enumerate(items): 
                if not (item in items[i]):
                    overlap = False
                    break
            if overlap:      
                s =  label + ' MATCH:' + '\t' + items[0][key]  +  "\t<--->\t" +  items[1][key]
                print s.encode('utf-8')
                overlap += 1
            overlapmem[item] = overlap 
    return total, overlap

    

print >>sys.stderr, "Computing source overlap"     
sourcetotal, sourceoverlap = computeoverlap(source, 'SOURCE')
               
print >>sys.stderr, "Computing target overlap"     
targettotal, targetoverlap = computeoverlap(target, 'TARGET')
                                
print >>sys.stderr, "Computing pair overlap"     
pairtotal, pairoverlap = computeoverlap(pairs, 'PAIR')


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
    
    
