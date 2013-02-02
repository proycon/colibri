#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import os


try:
    scorecount = int(sys.argv[1])
except:
    scorecount = 5;

try:    
    tweight = float(sys.argv[2])
except:
    tweight = 0.20


if not os.path.exists("model/moses.ini"):
    print >>sys.stderr, "No moses ini found"
    sys.exit(2)
    
intweight = False
for line in open("model/moses.ini"):
    line = line.strip()
    if not line: intweight = False
    if line[-16:] == '/phrase-table.gz':
        d = os.path.dirname(line.split(' ')[-1])
        if d[-5:] == "model": d = d[:-5]
        print "0 0 0 " + str(scorecount) + " " + d + "/tmp.phrasetable"
    elif line == '[weight-t]':
        intweight = True
        print line
        for i in range(0,scorecount):
            print str(tweight)           
    elif not intweight: 
        print line

