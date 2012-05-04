#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
from collections import defaultdict

try:
    phrasetable = sys.argv[1]
except:
    print >>sys.stderr, "Usage: flattenphrasetablescores.py phrasetable"
    
source = defaultdict(int)
target = defaultdict(int)

print >>sys.stderr, "Loading " + phrasetable
for line in codecs.open(phrasetable,'r','utf-8'):
    fields = line.split('|||')
    s = fields[0].strip()
    t = fields[1].strip()
    source[s] += 1
    target[t] += 1
    
print >>sys.stderr, "Loading " + phrasetable
for line in codecs.open(phrasetable,'r','utf-8'):
    fields = line.split('|||')
    s = fields[0].strip()
    t = fields[1].strip()
    out = s + ' ||| ' + t + ' ||| ' + str(1/float(source[s])) + ' ' + str(1/float(target[t]))
    print out.encode('utf-8')
