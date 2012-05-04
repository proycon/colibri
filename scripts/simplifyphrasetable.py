#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs


try:
    phrasetable = sys.argv[1]
except:
    print >>sys.stderr, "Usage: simplifyphrasetable.py phrasetable"

print >>sys.stderr, "Loading " + phrasetable
for line in codecs.open(phrasetable,'r','utf-8'):
    fields = line.split('|||')
    s = fields[0].strip()
    t = fields[1].strip()
    scores = fields[2].strip().split(' ')
    out =  s + ' ||| ' + t + ' ||| ' + str(scores[0]) + ' ' + str(scores[2])
    print out.encode('utf-8')
