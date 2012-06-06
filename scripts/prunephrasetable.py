#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs


try:
    phrasetable = sys.argv[1]
    sourcepatternsfile = sys.argv[2]
    targetpatternsfile = sys.argv[3]
except:
    print >>sys.stderr, "Usage: prunephrasetable.py phrasetable sourcepatterns targetpatterns"
    
print >>sys.stderr, "Loading " + sourcepatternsfile
sourcepatterns = {}
for line in codecs.open(sourcepatternsfile,'r','utf-8'):
    fields = line.split('\t')
    if len(fields) == 1:
        sourcepatterns[fields[0].strip()] = True
    else:
        sourcepatterns[fields[1].strip()] = True

print >>sys.stderr, "Loading " + targetpatternsfile
targetpatterns = {}
for line in codecs.open(targetpatternsfile,'r','utf-8'):
    fields = line.split('\t')
    if len(fields) == 1:
        targetpatterns[fields[0].strip()] = True
    else:
        targetpatterns[fields[1].strip()] = True
        

print >>sys.stderr, "Loading " + phrasetable
for line in codecs.open(phrasetable,'r','utf-8'):    
    fields = line.split('|||')
    s = fields[0].strip()
    t = fields[1].strip()
    if s in sourcepatterns and t in targetpatterns:
        print line.strip().encode('utf-8')
