#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import codecs
import itertools

try:
    phrasetable = sys.argv[1]
    column = sys.argv[2]
    patternmodel = sys.argv[3]
except:
    print >>sys.stderr, "Usage: comparephrasetables.py phrasetable column patternmodel"

phrases = ({},{})



print >>sys.stderr, "Loading " + phrasetable
for line in codecs.open(phrasetable,'r','utf-8'):
    fields = line.split('|||')
    phrase = fields[column].strip()
    phrases[0][phrase] += 1

print >>sys.stderr, "Loading " + patternmodel
for line in codecs.open(patternmodel,'r','utf-8'):
    fields = line.split('\t')
    phrase = fields[1].strip()
    phrases[1][phrase] += 1

    


# def shorten(item, items, pivot = 0):    
    # if isinstance(item, tuple):
        #TODO: all permutations
    # else:
        # words = item.split(' ')
        # l = 
        # while l >= 1:
            # words
        # yie
        
        

def computeoverlap(items, pivot = 0):    
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
                s =  'MATCH:' + '\t' + item[0]  +  "\t<--->\t" +  item[1]
                print s.encode('utf-8')
                overlap += 1
            overlapmem[item] = overlap 
    return total, overlap

    

print >>sys.stderr, "Computing overlap"     
total, overlap = computeoverlap(phrases, 'SOURCE')
               


print >>sys.stderr, "Total combined:\t" , total
for i in range(2):
    if i == 0: 
        label = phrasetable
    else:
        label = patternmodel    
    print >>sys.stderr, "Total in " + label + ":\t", len(phrases[i]), len(phrases[i]) / float(total)
print >>sys.stderr, "Overlapping:\t" , overlap, overlap / float(total)
for i in range(2):
    if i == 0: 
        label = phrasetable
    else:
        label = patternmodel
    print >>sys.stderr, "   ...as fraction of " + label + ":\t", overlap / float(len(phrases[i]))


    
