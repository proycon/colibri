#! /usr/bin/env python
#-*- coding: utf8 -*-

#Given a symmetric alignment produced by moses's symal, estimate a maximum likelihood lexical translation table. We estimate the w(e|f) as well as the inverse w(f|e) word translation table. This is a rewrite of a portion of Phillip Koehn's code in train-model.pl 

import sys
from collections import defaultdict

def get_lexical(sourcecorpusfilename, targetcorpusfilename, symalfilename, outputprefix):    
    sourcefile = open(sourcecorpusfilename)
    targetfile = open(targetcorpusfilename)
    symalfile = open(symalfilename)
    
    total_source = defaultdict(int)
    total_target = defaultdict(int)
    word_translation = {None: defaultdict(int) }
    
    while True:
        symalline = symalfile.read()
        if not symalline: 
            break
        try:
            sourceline = targetfile.read()
        except:
            print >>sys.stderr,"ERROR: Premature end of " + targetcorpusfilename
            break
        try:
            targetline = targetfile.read()
        except:
            print >>sys.stderr,"ERROR: Premature end of " + targetcorpusfilename
            break
    
        sourcewords = sourceline.strip().split(' ')
        targetwords = targetline.strip().split(' ')
        alignmentpoints = symalline.strip().split(' ')
    
        source_aligned = defaultdict(int)
        target_aligned = defaultdict(int)    
    
        for alignmentpoint in alignmentpoints:
            source, target = ( int(x) for x in  alignmentpoint.split('-') )
            if source >= len(sourcewords) or target >= len(targetwords):
                print >>sys.stderr,"WARNING: Alignment point " + str(source) + '-' + str(target) +  " is out of range, ignoring"
            else: 
                source_aligned[source] += 1
                target_aligned[target] += 1
                
                total_source[sourcewords[source]] += 1
                total_target[targetwords[target]] += 1
                
                if not source in word_translation: word_translation[source] = defaultdict(int)
                word_translation[sourcewords[source]][targetwords[target]] += 1
    
    
        #unaligned words
        for i, sourceword in enumerate(sourcewords):
            if i in source_aligned:
                continue
            else:
                word_translation[sourceword][None] += 1            
                total_source[sourceword] += 1
                total_target[None] += 1
        for i, targetword in enumerate(targetwords):
            if i in target_aligned:
                continue
            else:
                word_translation[None][targetword] += 1            
                total_source[None] += 1
                total_target[targetword] += 1                
      
    sourcefile.close()
    targetfile.close()
    symalfile.close()

    try:         
        s2tfile = open(outputprefix + '.s2t', 'w')
        t2sfile = open(outputprefix + '.t2s', 'w')
    except:
        print >>sys.stderr, "Error opening output files"
        sys.exit(2)
        
    for source, data in word_translation.items():
        for target, count in data.items():
            s2tfile.write("%s %s %.7f\n"%(target,source,count/float(total_source[source])))
            t2sfile.write("%s %s %.7f\n" %(source,target,count/float(total_target[target])))
      
    print >>sys.stderr, "Done";

if __name__ == "__main__":    
    try:
        sourcecorpusfilename = sys.argv[1]
        targetcorpusfilename = sys.argv[2]
        symalfilename = sys.argv[3]
        outputprefix = sys.argv[4]
    except:
        print >>sys.stderr,"Usage: moses-lexicaltranslationtable.py sourcecorpusfilename targetcorpusfilename symalfilename outputprefix"
        sys.exit(1)
    get_lexical(sourcecorpusfilename,targetcorpusfilename,symalfilename,outputprefix)
    
    
