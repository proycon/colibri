#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
from collections import defaultdict
import itertools

try:
    phrasetables = sys.argv[1:3]
    scoreside = int(sys.argv[3])
except:
    print >>sys.stderr, "Usage: intersectphrasetables.py phrasetable1 phrasetable2 [scores=0/1]"

pairs = []
scores = {}

for i, phrasetable in enumerate(phrasetables):
    pairs.append( defaultdict(int) )
    
    print >>sys.stderr, "Loading " + phrasetable
    for line in codecs.open(phrasetable,'r','utf-8'):
        fields = line.split('|||')
        s = fields[0].strip()
        t = fields[1].strip()
        if i == scoreside:
            sc = fields[2].strip().split()
            if len(sc) == 5: #strip lexical weighting:
                sc = (sc[0],sc[2])
            scores[(s,t)] = sc        
        pairs[-1][(s,t)] += 1
                
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
            overlaps = True
            for i, _ in enumerate(items): 
                if not (item in items[i]):
                    overlaps = False
                    break
    
            if overlaps:      
                s = item[0] + ' ||| ' + item[1] + ' ||| ' + ' '.join(scores[item]) 
                print s.encode('utf-8')
                overlap += 1
            overlapmem[item] = overlap 
    return total, overlap

    

print >>sys.stderr, "Computing pair overlap"     
pairtotal, pairoverlap = computeoverlap(pairs, 'PAIR')


print >>sys.stderr, "Total combined pairs:\t" , pairtotal
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "Total pairs in " + phrasetable + ":\t", len(pairs[i]), len(pairs[i]) / float(pairtotal)
print >>sys.stderr, "Overlapping pairs:\t" , pairoverlap,  pairoverlap / float(pairtotal)
for i, phrasetable in enumerate(phrasetables):
    print >>sys.stderr, "   ...as fraction of " + phrasetable + ":\t",pairoverlap / float(len(pairs[i]))
