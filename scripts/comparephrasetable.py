#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
from collections import defaultdict

try:
    phrasetables = sys.argv[1:]
except:
    print >>sys.stderr, "Usage: comparephrasetable.py phrasetable1 phrasetable2"

pairs = []
source = []
target = []

for phrasetable in phrasetables:
    pairs.append( defaultdict(int) )
    source.append( defaultdict(int) )
    target.append( defaultdict(int) )
        
    for line in codecs.open(phrasetable,'r','utf-8'):
        fields = line.split('|||')
        s = fields[0].strip()
        t = fields[1].strip()
        pairs[-1][(s,t)] += 1
        source[-1][s] += 1
        target[-1][t] += 1
    

overlapmem = {}    
for pair in ( x.keys() for x in pairs ):    
    overlap = True
    for i, phrasetable in enumerate(phrasetables): 
        if not (pair in pairs[i]):
            overlap = False
            break
    if overlap:
        if not pair in overlapmem:
            print "PAIR MATCH:", pair
        overlapmem[pair] = True
        
overlapmem = {}        
for s in ( x.keys() for x in source ):    
    overlap = True
    for i, phrasetable in enumerate(phrasetables): 
        if not (s in source[i]):
            overlap = False
            break
    if overlap:
        if not s in overlapmem:
            print "SOURCE MATCH:", s
        overlapmem[s] = True

overlapmem = {}
for t in ( x.keys() for x in target ):    
    overlap = True
    for i, phrasetable in enumerate(phrasetables): 
        if not (t in target[i]):
            overlap = False
            break
    if overlap:
        if not t in overlapmem:
            print "TARGET MATCH:",t
        overlapmem[t] = True


