#! /usr/bin/env python
# -*- coding: utf8 -*-

# MT Translation Wrapper
#   by Maarten van Gompel
#   Radboud University Nijmegen
#   Licensed under GPL v3

import sys
import os
import subprocess
import getopt
import codecs
import glob
from pynlpl.evaluation import filesampler

def bold(s):
   CSI="\x1B["
   return CSI+"1m" + s + CSI + "0m"
   
def white(s):
   CSI="\x1B["
   return CSI+"37m" + s + CSI + "0m"   


def red(s):
   CSI="\x1B["
   return CSI+"31m" + s + CSI + "0m"
   
def green(s):
   CSI="\x1B["
   return CSI+"32m" + s + CSI + "0m"   


def yellow(s):
   CSI="\x1B["
   return CSI+"33m" + s + CSI + "0m"   

   
def blue(s):
   CSI="\x1B["
   return CSI+"34m" + s + CSI + "0m"   
   

def magenta(s):
   CSI="\x1B["
   return CSI+"35m" + s + CSI + "0m"   

class MTWrapper(object):
    defaults = [
            ('WORKDIR','','Full path to the working directory that holds all data for the system'),
            ('CORPUSNAME', '','The name of the corpus (without language codes)'),
            ('EXPERIMENTNAME', '','A unique name for this experiment (optional)'),
            ('SOURCELANG', '','A language code identifying the source language'),
            ('TARGETLANG', '','A language code identifying the target language'),
            ('TRAINSOURCECORPUS', '','The file containing to the source-language part of the parallel corpus used for training, one sentence per line'),
            ('TRAINTARGETCORPUS', '','The file containing to the target-language part of the parallel corpus used for training, one sentence per line'),
            ('TESTSOURCECORPUS', '','The file containing to the source-language part of the parallel corpus used for testing, one sentence per line'),
            ('TESTTARGETCORPUS', '','The file containing to the target-language part of the parallel corpus used as reference for testing.'),
            ('DEVSOURCECORPUS', '','The file containing to the source-language part of the parallel corpus used for parameter tuning, one sentence per line'),
            ('DEVTARGETCORPUS', '','The file containing to the target-language part of the parallel corpus used for parameter tuning, one sentence per line. Used as reference in parameter tuning.'),           
            ('TOKENIZE_SOURCECORPUS', False,''),
            ('TOKENIZE_TARGETCORPUS', False,''),
            ('BUILD_SRILM_SOURCEMODEL',False,'Build a source-language model'),
            ('BUILD_SRILM_TARGETMODEL',False,'Build a target-language model'),
            ('BUILD_GIZA_WORDALIGNMENT',False,'Build GIZA++ Word Alignments'),
            ('BUILD_GIZA_WORDALIGNMENT_REV',False,'Build GIZA++ Reverse Word Alignment (target to source)'),
            ('BUILD_GIZA_WORDALIGNMENT_COOC',False,'Output extra co-occurence data'),
            ('BUILD_MOSES_SYMAL', False,'Symmetrise word alignments'),
            ('BUILD_MOSES_WORDTRANSTABLE', False,'Build lexical translation table'),
            ('BUILD_MOSES_PHRASEEXTRACT', False,'Extract phrases'),            
            ('BUILD_MOSES_PHRASETRANSTABLE', False,'Build phrase translation table'),
            ('BUILD_MOSES_MEMSCORE', False,'Use memscore to score phrases rather than the default phrase-extract scorer'),
            ('BUILD_MOSES', False,'Build moses configuration, necessary for decoding using moses'),
            ('BUILD_MOSES_MERT', False,'Do Minimum Error Rate Training for Moses (on development set)'),                       
            ('BUILD_PBMBMT', False, 'Build model for Phrase-Based Memory-based Machine Translation'),   
            ('BUILD_PBMBMT_PARAMSEARCH', False, 'Do parameter optimisation for PBMBMT using wrapped progressive sampling'),
            ('BUILD_COLIBRI_ALIGNMENT', False,'Create an alignment using colibri'),
            ('PATH_MOSES', '','Base directory where Moses is installed'),
            ('PATH_SRILM', '','Base directory where SRILM is installed'),
            ('PATH_GIZA', '','Base directory where GIZA++ is installed'),
            ('PATH_COLIBRI', '','Base directory where COLIBRI is installed'),
            ('PATH_MATREX','','Base directory for Matrex evaluation scripts'),
            ('PATH_PBMBMT','','Base directory to PBMBMT'),
            ('EXEC_UCTO', 'ucto','Path to ucto binary'),
            ('EXEC_SRILM', 'ngram-count','Path to ngram-count (SRILM)'),
            ('EXEC_TIMBL', 'timbl','Path to timbl binary'),
            ('EXEC_GIZA_MKCLS', 'mkcls','Path to mkcls (part of GIZA++)'),
            ('EXEC_GIZA', 'GIZA++','Path to GIZA++ binary'),
            ('EXEC_GIZA_PLAIN2SNT', 'plain2snt.out','Path to plain2snt.out (part of GIZA++)'),
            ('EXEC_GIZA_SNT2COOC', 'snt2cooc.out','Path to snt2cooc.out (part of GIZA++)'),
            ('EXEC_PERL', 'perl','Path to perl binary'),
            ('EXEC_JAVA', 'java','Path to java binary'),
            ('EXEC_MOSES', 'moses','Path to Moses binary'),
            ('EXEC_MOSES_GIZA2BAL', 'scripts/training/symal/giza2bal.pl', ''),
            ('EXEC_MOSES_SYMAL', 'scripts/training/symal/symal', ''),
            ('EXEC_MOSES_WORDTRANSTABLE','scripts/moses-lexicaltranslationtable.py',''),
            ('EXEC_MOSES_PHRASEEXTRACT','scripts/training/phrase-extract/extract',''),
            ('EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE','scripts/training/phrase-extract/consolidate',''),
            ('EXEC_MOSES_PHRASEEXTRACT_SCORE','scripts/training/phrase-extract/score',''),
            ('EXEC_MOSES_MEMSCORE','scripts/training/memscore/memscore',''),
            ('EXEC_MOSES_MERT','scripts/training/mert-moses.pl',''),
            ('EXEC_MATREX_WER','eval/WER_v01.pl',''),
            ('EXEC_MATREX_PER','eval/PER_v01.pl',''),
            ('EXEC_MATREX_BLEU','eval/bleu-1.04.pl',''),
            ('EXEC_MATREX_METEOR','meteor-0.6/meteor.pl',''),
            ('EXEC_MATREX_MTEVAL','mteval-v11b.pl','NIST and BLEU'),
            ('EXEC_MATREX_TER','tercom.jar',''),
            ('EXEC_PBMBMT_DECODER','pbmbmt-decode',''),
            ('EXEC_PBMBMT_INSTANCEGENERATOR','instancegenerator.py',''),
            ('EXEC_COLIBRI_CLASSENCODE','classencode',''),
            ('EXEC_COLIBRI_PATTERNFINDER','patternfinder',''),
            ('EXEC_COLIBRI_GRAPHMODEL','graphmodel',''),
            ('EXEC_COLIBRI_ALIGNER','aligner',''),
            ('MKCLS_OPTIONS','-m2 -c50',''),
            ('GIZA_OPTIONS','-p0 0.98 -m1 5 -m2 0 -m3 3 -m4 3 -nsmooth 4 -model4smoothfactor 0.4',''),
            ('SRILM_ORDER',3,'N-gram size for language model'),
            ('SRILM_OPTIONS','-interpolate -kndiscount','Further SRILM options (do not use -order here, use SRILM_ORDER instead)'),
            ('UCTO_OPTIONS','-m -n',''),
            ('SYMAL_OPTIONS','-alignment=grow -diagonal=yes -final=yes -both=no',''), #-hmmiterations 5 -hmmdumpfrequency -5'
            ('MOSES_MERT_OPTIONS','','See http://www.statmt.org/moses/?n=FactoredTraining.Tuning'),
            ('MOSES_MEMSCORE_METHOD','ml','Memscore scoring method:  ml: maximum likelihood, wittenbell: Witten-Bell smoothing, absdiscount: Absolute discounting'),
            ('PHRASEEXTRACT_MAX_PHRASE_LENGTH',7,''),
            ('PHRASEEXTRACT_REORDERING_FLAGS','',''), #" --model wbe-mslr --model phrase-mslr --model hier-mslr" #Maximum lexical reordering 
            ('PHRASESCORE_OPTIONS', '',''), #--Hierarchical --WordAlignment (--Inverse)
            ('PBMBMT_PHRASETABLE','','Use the following pre-existing phrase-table (rather than depending on Moses to create one from scratch)'),
            ('PBMBMT_GIZAALIGNMENT','','Use the following pre-existing GIZA Word Alignment (rather than depending on GIZA++ to create one from scratch)'),
            ('PBMBMT_MAXPHRASELENGTH',6,''),
            ('PBMBMT_LEFTCONTEXTSIZE',1,''),
            ('PBMBMT_RIGHTCONTEXTSIZE',1,''),
            ('PBMBMT_DECODER_OPTIONS','','Options for PBMBMT Decoder (do not include --srilm=, will be added automatically if BUILD_SRILM_TARGETMODEL is enabled)'),
            ('PBMBMT_TIMBL_OPTIONS','-k 1 -a4','Timbl options (+v+db+di is added automatically). See Timbl -h'),
            ('COLIBRI_GRAPHMODEL_OPTIONS','','Options for the Graphmodel, if empty, no graph model will be constructed for the aligner, see graphmodel -h'),
            ('COLIBRI_PATTERNFINDER_OPTIONS','-t 10 -s -B -E', 'Options for the pattern finder, see patternfinder -h'),
            ('COLIBRI_ALIGNER_OPTIONS','-E -I 100','Options for the colibri aligner, see aligner -h'),  
    ]

    
    def parsekwargs(self, key, default, **kwargs):
        if key in kwargs:
            return kwargs[key]
            del kwargs[key]
        else:
            return default
        

    def __init__(self, **kwargs):                
        
        for key, default, help in MTWrapper.defaults:
            if key in kwargs:
                setattr(self,key,kwargs[key])
                del kwargs[key]
            else:
                setattr(self,key,default)
        
        self.EXEC_PERL = self.findpath(self.EXEC_PERL)
        self.EXEC_JAVA = self.findpath(self.EXEC_JAVA)    
        
        self.EXEC_TIMBL = self.findpath(self.EXEC_TIMBL)
        self.EXEC_SRILM = self.findpath(self.EXEC_SRILM,self.PATH_SRILM)
        self.EXEC_UCTO = self.findpath(self.EXEC_UCTO)
        
        self.EXEC_GIZA = self.findpath(self.EXEC_GIZA,self.PATH_GIZA)
        self.EXEC_GIZA_PLAIN2SNT = self.findpath(self.EXEC_GIZA_PLAIN2SNT,self.PATH_GIZA)
        self.EXEC_GIZA_SNT2COOC = self.findpath(self.EXEC_GIZA_SNT2COOC, self.PATH_GIZA)
        self.EXEC_GIZA_MKCLS = self.findpath(self.EXEC_GIZA_MKCLS, self.PATH_GIZA)
                
        self.EXEC_MOSES = self.findpath(self.EXEC_MOSES,self.PATH_MOSES)        
        self.EXEC_MOSES_GIZA2BAL = self.findpath(self.EXEC_MOSES_GIZA2BAL,self.PATH_MOSES)
        self.EXEC_MOSES_SYMAL = self.findpath(self.EXEC_MOSES_SYMAL,self.PATH_MOSES)
        self.EXEC_MOSES_WORDTRANSTABLE = self.findpath(self.EXEC_MOSES_WORDTRANSTABLE,self.PATH_COLIBRI)
        self.EXEC_MOSES_PHRASEEXTRACT = self.findpath(self.EXEC_MOSES_PHRASEEXTRACT,self.PATH_MOSES)
        self.EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE = self.findpath(self.EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE,self.PATH_MOSES)
        self.EXEC_MOSES_PHRASEEXTRACT_SCORE = self.findpath(self.EXEC_MOSES_PHRASEEXTRACT_SCORE,self.PATH_MOSES)
        self.EXEC_MOSES_MERT  = self.findpath(self.EXEC_MOSES_MERT , self.PATH_MOSES)
        
        self.EXEC_MATREX_WER = self.findpath(self.EXEC_MATREX_WER, self.PATH_MATREX)
        self.EXEC_MATREX_PER = self.findpath(self.EXEC_MATREX_PER, self.PATH_MATREX)
        self.EXEC_MATREX_BLEU = self.findpath(self.EXEC_MATREX_BLEU, self.PATH_MATREX)
        self.EXEC_MATREX_METEOR = self.findpath(self.EXEC_MATREX_METEOR, self.PATH_MATREX)
        self.EXEC_MATREX_MTEVAL = self.findpath(self.EXEC_MATREX_MTEVAL, self.PATH_MATREX)
        self.EXEC_MATREX_TER = self.findpath(self.EXEC_MATREX_TER, self.PATH_MATREX)
                
        self.EXEC_COLIBRI_CLASSENCODE = self.findpath(self.EXEC_COLIBRI_CLASSENCODE , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_PATTERNFINDER = self.findpath(self.EXEC_COLIBRI_PATTERNFINDER , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_GRAPHMODEL = self.findpath(self.EXEC_COLIBRI_GRAPHMODEL , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_ALIGNER = self.findpath(self.EXEC_COLIBRI_ALIGNER , self.PATH_COLIBRI)
        
        self.EXEC_PBMBMT_DECODER = self.findpath(self.EXEC_PBMBMT_DECODER, self.PATH_PBMBMT)
        self.EXEC_PBMBMT_INSTANCEGENERATOR = self.findpath(self.EXEC_PBMBMT_INSTANCEGENERATOR, self.PATH_PBMBMT)
        
        if self.PATH_MOSES:
            self.PATH_MOSES_MERT = self.PATH_MOSES + '/mert'
        else:        
            self.PATH_MOSES_MERT = ''            
        
        for key in kwargs:
            print >>sys.stderr, "Unknown configuration directive: " + key
            sys.exit(2)
            
        
        

    def findpath(self, name, basepath = ''):                        
        for path in os.environ['PATH'].split(':'):
            if os.path.exists(path + '/' + name) and not os.path.isdir(path + '/' + name):
                print >>sys.stderr, green("Found " + name + " in " + path)
                return path + '/' + name
        if basepath and os.path.exists(basepath + '/' + name): 
            print >>sys.stderr, green("Found " + name + " in " + basepath)
            return basepath + '/' + name        
        print >>sys.stderr, yellow("Warning: Did not find " + name)
        return ""

    def check_common(self):
        sane = True
        if self.WORKDIR[-1] != '/':     
            self.WORKDIR += '/'
            if not os.path.isdir(self.WORKDIR):
                print >>sys.stderr,yellow("Work directory does not exist, creating " + self.WORKDIR)
                try:
                    os.mkdir(self.WORKDIR)
                except:
                    print >>sys.stderr, red("Configuration error: Unable to create work directory " + self.WORKDIR)
                    sane = False
        
        sane = True
        if not self.CORPUSNAME:
            print >>sys.stderr,red("Configuration error: CORPUSNAME not specified!")
            sane = False
        if not self.SOURCELANG:
            print >>sys.stderr,red("Configuration error: SOURCELANG not specified!")
            sane = False
        if not self.TARGETLANG:
            print >>sys.stderr,red("Configuration error: TARGETLANG not specified!")
            sane = False
        return sane


    def check_run(self):        
        if self.BUILD_MOSES:
            if not self.EXEC_MOSES or not os.path.isfile(self.EXEC_MOSES):
                print >>sys.stderr,bold(red("Error: Moses not found!"))
                return False
            elif not os.path.exists(self.WORKDIR + '/moses.ini'):
                print >>sys.stderr,bold(red("Error: No Moses configuration found. Did you forget to train the system first?"))
                return False
            elif not os.path.exists(self.gets2tfilename('phrasetable')):
                print >>sys.stderr,bold(red("Error: No Moses phrasetable found ("+ self.gets2tfilename('phrasetable')+") . Did you forget to train the system first?"))
                return False
        elif self.BUILD_PBMBMT:
            #TODO: implement
            return False
        else:
            print >>sys.stderr,bold(red("Error: System is not runnable, no MT decoder enabled"))
            return False
        return True
    
    def check_test(self):
        if not (self.EXEC_MATREX_WER or self.EXEC_MATREX_PER or self.EXEC_MATREX_BLEU or self.EXEC_MATREX_MTEVAL or self.EXEC_MATREX_METEOR or self.EXEC_MATREX_TER):
            print >>sys.stderr,bold(red("Error: No evaluation scripts found, set at least one of EXEC_MATREX_* and PATH_MATREX"))
            return False            
        return True

    def check_train(self):
        sane = True                            
        if not self.TRAINSOURCECORPUS:
            print >>sys.stderr,red("Configuration error: TRAINSOURCECORPUS not specified!")
            sane = False
        if not self.TRAINTARGETCORPUS:
            print >>sys.stderr,red("Configuration error: TRAINTARGETCORPUS not specified!")
            sane = False
            
        if (self.TOKENIZE_SOURCECORPUS or self.TOKENIZE_TARGETCORPUS) and (not self.EXEC_UCTO or not os.path.isfile(self.EXEC_UCTO)):
            print >>sys.stderr,red("Dependency error: ucto not found (EXEC_UCTO=" + self.EXEC_UCTO + ")")
            sane = False
            
        if self.BUILD_COLIBRI_ALIGNMENT:
            if not self.EXEC_COLIBRI_PATTERNFINDER or not os.path.exists(self.EXEC_COLIBRI_PATTERNFINDER):
                sane = False
                print >>sys.stderr,red("Configuration error: EXEC_COLIBRI_PATTERNFINDER not found ! Required for BUILD_COLIBRI_ALIGNMENT !")            
            if not self.EXEC_COLIBRI_ALIGNER or not os.path.exists(self.EXEC_COLIBRI_ALIGNER):
                sane = False
                print >>sys.stderr,red("Configuration error: EXEC_COLIBRI_ALIGNER not found ! Required for BUILD_COLIBRI_ALIGNMENT !")
            if not self.EXEC_COLIBRI_GRAPHMODEL or not os.path.exists(self.EXEC_COLIBRI_GRAPHMODEL):
                sane = False
                print >>sys.stderr,red("Configuration error: EXEC_COLIBRI_GRAPHMODEL not found ! Required for BUILD_COLIBRI_ALIGNMENT !")                
            if not self.EXEC_COLIBRI_CLASSENCODE or not os.path.exists(self.EXEC_COLIBRI_CLASSENCODE):
                sane = False
                print >>sys.stderr,red("Configuration error: EXEC_COLIBRI_CLASSENCODE not found ! Required for BUILD_COLIBRI_ALIGNMENT !")
                
            

        if self.BUILD_MOSES_MERT:
            if not self.BUILD_MOSES: 
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES automatically enabled because BUILD_MOSES_MERT is too")
                self.BUILD_MOSES = True
            if not self.DEVSOURCECORPUS or not os.path.exists(self.DEVSOURCECORPUS):
                sane = False
                print >>sys.stderr,red("Configuration error: DEVSOURCECORPUS not found! Required for BUILD_MOSES_MERT !")
            if not self.DEVTARGETCORPUS or not os.path.exists(self.DEVTARGETCORPUS):
                sane = False
                print >>sys.stderr,red("Configuration error: DEVTARGETCORPUS not found! Required for BUILD_MOSES_MERT !")
            if not self.PATH_MOSES_MERT or not os.path.isdir(self.PATH_MOSES_MERT):
                sane = False
                print >>sys.stderr,red("PATH_MOSES_MERT not found, please set PATH_MOSES !")

        if self.BUILD_MOSES:
            if self.BUILD_PBMBMT:
                print >>sys.stderr,red("Configuration error: Ambiguous selection of MT system: Select only one of BUILD_MOSES or BUILD_PBMBMT")
                sane = False
            if not self.BUILD_MOSES_PHRASETRANSTABLE:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_MOSES is too")
                self.BUILD_MOSES_PHRASETRANSTABLE = True
            if not self.BUILD_SRILM_TARGETMODEL:                 
                print >>sys.stderr,yellow("Configuration update: BUILD_SRILM_TARGETMODEL automatically enabled because BUILD_MOSES is too")
                self.BUILD_SRILM_TARGETMODEL = True
            if not self.EXEC_MOSES or not os.path.isfile(self.EXEC_MOSES):
                sane = False
                print >>sys.stderr,red("Moses not found! Set EXEC_MOSES or PATH_MOSES !")                
                
                
                
        if self.BUILD_PBMBMT:
            if not self.PBMBMT_PHRASETABLE and not self.BUILD_MOSES_PHRASETABLE:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_PBMBMT is enabled and no pre-existing phrasetable is set (PBMBMT_PHRASETABLE)")
            if not self.PBMBMT_GIZAALIGNMENT and not self.BUILD_GIZA_WORDALIGNMENT:
                print >>sys.stderr,yellow("Configuration update: BUILD_GIZA_WORDALIGNMENT automatically enabled because BUILD_PBMBMT is enabled and no pre-existing word alignment file is set (PBMBMT_GIZAALIGNMENT)")               
            if not self.BUILD_SRILM_TARGETMODEL:
                print >>sys.stderr,yellow("Configuration update: BUILD_SRILM_TARGETMODEL automatically enabled because BUILD_PBMBMT is too")
                self.BUILD_SRILM_TARGETMODEL = True         
            if self.PBMBMT_PHRASETABLE:
                if not os.path.isfile(self.PBMBMT_PHRASETABLE):
                    print >>sys.stderr,yellow("Configuration error: PBMBMT_PHRASETABLE does not exist!")
                    sane = False
                else:
                    os.symlink(self.PBMBMT_PHRASETABLE, self.gets2tfilename('phrasetable'))
            if self.PBMBMT_GIZAALIGNMENT: 
                if not os.path.isfile(self.PBMBMT_GIZAALIGNMENT):
                    print >>sys.stderr,yellow("Configuration error: PBMBMT_GIZAALIGNMENT does not exist!")
                    sane = False
                else:
                    os.symlink(self.PBMBMT_GIZAALIGNMENT, self.gets2tfilename('A3.final'))
            if not self.EXEC_PBMBMT_DECODER or not os.path.isfile(self.EXEC_PBMBMT_DECODER):
                sane = False
                print >>sys.stderr,red("PBMBMT decoder not found! Set EXEC_PBMBMT_DECODER or PATH_PBMBMT !")                
            if not self.EXEC_PBMBMT_INSTANCEGENERATOR or not os.path.isfile(self.EXEC_PBMBMT_INSTANCEGENERATOR):
                sane = False
                print >>sys.stderr,red("PBMBMT instance generator not found! Set EXEC_PBMBMT_DECODER or PATH_PBMBMT !")
            if not self.EXEC_TIMBL or not os.path.isfile(self.EXEC_TIMBL):
                sane = False
                print >>sys.stderr,red("TiMBL was not found, but is required for PBMBMT! Set EXEC_TIMBL or PATH_TIMBL !")



        
        if self.BUILD_MOSES_MEMSCORE:
            if not self.BUILD_MOSES_PHRASETRANSTABLE:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_MOSES_MEMSCORE is too")
                self.BUILD_MOSES_PHRASETRANSTABLE = True  
            
        if self.BUILD_MOSES_PHRASETRANSTABLE:
            if not self.BUILD_MOSES_PHRASEEXTRACT:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_PHRASEEXTRACT automatically enabled because BUILD_MOSES_PHRASECORE is too")
                self.BUILD_MOSES_PHRASEEXTRACT = True        
                
        
        if self.BUILD_MOSES_PHRASEEXTRACT:
            if not self.BUILD_MOSES_SYMAL:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_WORDTRANSTABLE automatically enabled because BUILD_MOSES_PHRASEEXTRACT is too")
                self.BUILD_MOSES_WORDTRANSTABLE = True
                
        
        if self.BUILD_MOSES_WORDTRANSTABLE:
            if not self.BUILD_MOSES_SYMAL:
                print >>sys.stderr,yellow("Configuration update: BUILD_MOSES_SYMAL automatically enabled because BUILD_MOSES_WORDTRANSTABLE is too")
                self.BUILD_MOSES_SYMAL = True
        
                
        if self.BUILD_MOSES_SYMAL:            
            if not self.BUILD_GIZA_WORDALIGNMENT:
                print >>sys.stderr,yellow("Configuration update: BUILD_GIZA_WORDALIGNMENT automatically enabled because BUILD_MOSES_SYMAL is too")
                self.BUILD_GIZA_WORDALIGNMENT = True
            if not self.BUILD_GIZA_WORDALIGNMENT_REV:
                print >>sys.stderr,yellow("Configuration update: BUILD_GIZA_WORDALIGNMENT_REV automatically enabled because BUILD_MOSES_SYMAL is too")
                self.BUILD_GIZA_WORDALIGNMENT_REV = True                            

            if not self.EXEC_MOSES_GIZA2BAL or not os.path.isfile(self.EXEC_MOSES_GIZA2BAL):
                sane = False
                print >>sys.stderr,red("Dependency error: giza2bal.pl (provided by Moses) not found (EXEC_MOSES_GIZA2BAL=" + self.EXEC_MOSES_GIZA2BAL + ")")
               
            
            if not self.EXEC_MOSES_SYMAL or not os.path.isfile(self.EXEC_MOSES_SYMAL):
                sane = False
                print >>sys.stderr,red("Dependency error: symal (provided by Moses) not found (EXEC_MOSES_SYMAL=" + self.EXEC_MOSES_SYMAL + ")")
               
            
                
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.EXEC_GIZA or not os.path.isfile(self.EXEC_GIZA)): 
            print >>sys.stderr,red("Dependency error: GIZA++ not found (EXEC_GIZA=" + self.EXEC_GIZA + ")")
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.EXEC_GIZA_PLAIN2SNT or not os.path.isfile(self.EXEC_GIZA_PLAIN2SNT)): 
            print >>sys.stderr,red("Dependency error: plain2snt.out (provided by GIZA++) not found (EXEC_GIZA_PLAIN2SNT=" + self.EXEC_GIZA_PLAIN2SNT + ")")            
            sane = False
        if self.BUILD_GIZA_WORDALIGNMENT_COOC and (not self.EXEC_GIZA_SNT2COOC or not os.path.isfile(self.EXEC_GIZA_SNT2COOC)): 
            print >>sys.stderr,red("Dependency error: snt2cooc.out (provided by GIZA++) not found (EXEC_GIZA_SNT2COOC=" + self.EXEC_GIZA_SNT2COOC + ")")            
            sane = False                        
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.EXEC_GIZA_MKCLS or not os.path.isfile(self.EXEC_GIZA_MKCLS)): 
            print >>sys.stderr,red("Dependency error: mkcls (provided by GIZA++) not found (EXEC_GIZA_MKCLS=" + self.EXEC_GIZA_MKCLS + ")")            
            sane = False                            
        if (self.BUILD_SRILM_TARGETMODEL or self.BUILD_SRILM_SOURCEMODEL) and (not self.EXEC_SRILM or not os.path.isfile(self.EXEC_SRILM)):
            print >>sys.stderr,red("Dependency error: ngram-count (provided by SRILM) not found (EXEC_SRILM=" + self.EXEC_SRILM + ")")
            sane = False
        return sane
    
    

    def getsourcefilename(self, extension):
        return self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '.' + extension

    def gettargetfilename(self, extension):
        return self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.TARGETLANG + '.' + extension
    
    def gets2tfilename(self, extension='', longform =False):
        s = self.WORKDIR + self.CORPUSNAME + '-' + self.SOURCELANG 
        if longform: s += '_' + self.CORPUSNAME 
        s += '-' + self.TARGETLANG 
        if extension: s += '.' + extension
        return s
    
    def gett2sfilename(self, extension='',longform=False):
        s = self.WORKDIR + self.CORPUSNAME + '-' + self.TARGETLANG 
        if longform: s += '_' + self.CORPUSNAME 
        s += '-' + self.SOURCELANG 
        if extension: s += '.' + extension
        return s    
            


    def usage(self):
        print >>sys.stderr,"Usage: " + os.path.basename(sys.argv[0]) + ' [command]'
        print >>sys.stderr,"Commands:"
        print >>sys.stderr,"\ttrain                              Train the MT system"
        print >>sys.stderr,"\trun <inputfile> [options]          Run the MT system on the specified input file"
        print >>sys.stderr,"\t\t-t                               Tokenise the input file"
        print >>sys.stderr,"\t\t-o <outputfile>                  Output file (default: stdout)"        
        print >>sys.stderr,"\ttest <inputfile> <referencefile>   Run and evaluate the MT system on the specified input file and reference file (one sentence per line). If <inputfile> and <referencefile> are not given, the default test files from the system configuration are used."                               
        print >>sys.stderr,"\tscore <inputfile> <referencefile>  Like test, but work on a pre-run system, does  not run the translation again."                   
        print >>sys.stderr,"\tclean [all|giza|moses|colibri|score]    Clean generated files"
        print >>sys.stderr,"\tbranch <new-directory>             Create a new branch based on this project (files are symlinked instead of copied)"

    def start(self):        
        try:
            cmd = sys.argv[1]
        except:
            self.usage()
            sys.exit(2)
        if cmd == 'train':
            if not self.starttrain():
                sys.exit(1)                
        elif cmd == 'clean':
            targets = sys.argv[2:]
            if not self.clean(targets):
                sys.exit(1)
        elif cmd == 'run':
            try:
                inputfile = sys.argv[2]
            except:
                print >>sys.stderr, "ERROR: Expected: run <inputfile>"                
                self.usage()
                sys.exit(2)
            try:
                opts, args = getopt.getopt(sys.argv[3:], "to:")
            except getopt.GetoptError, err:
                # print help information and exit:
                print str(err) # will print something like "option -a not recognized"
                self.usage()
                sys.exit(2)
                
            outputfile = 'output.txt'
            tokenise = False             
            for o, a in opts:
                if o == "-t":
                    tokenise = True
                elif o in ("-o", "--output"):
                    outputfile = a
                else:
                    assert False, "unhandled option"
            
            if self.run(inputfile, outputfile, tokenise):
                os.system('cat ' + outputfile)                
            else:
                print >>sys.stderr, "An error occurred whilst trying to run the system" 
                sys.exit(1)
    
        elif cmd == 'branch':          
            
            #TODO: IMplement
            raise NotImplemented
            
        elif cmd == 'test':                        
            if len(sys.argv) == 4:
                inputfile = sys.argv[2]
                referencefile = sys.argv[3]
            elif len(sys.argv) == 2:
                if not self.TESTSOURCECORPUS or not os.path.exists(self.TESTSOURCECORPUS):
                    print >>sys.stderr, bold(red("No predefined default test corpus set for input (set TESTSOURCECORPUS) or specify on command line"))
                    sys.exit(2)
                if not self.TESTTARGETCORPUS or not os.path.exists(self.TESTTARGETCORPUS):
                    print >>sys.stderr, bold(red("No predefined default test corpus set for reference (set TESTSOURCECORPUS and TESTTARGETCORPUS) or specify on command line "))
                    sys.exit(2)                    
                inputfile = self.TESTSOURCECORPUS                
                referencefile = self.TESTTARGETCORPUS                
            else:
                self.usage()
                sys.exit(2)
            
            if not self.test(inputfile, referencefile): 
                sys.exit(1)

        elif cmd == 'score':                        
            if len(sys.argv) == 4:
                inputfile = sys.argv[2]
                referencefile = sys.argv[3]
            elif len(sys.argv) == 2:
                if not self.TESTSOURCECORPUS or not os.path.exists(self.TESTSOURCECORPUS):
                    print >>sys.stderr, bold(red("No predefined default test corpus set for input (set TESTSOURCECORPUS) or specify on command line"))
                    sys.exit(2)
                if not self.TESTTARGETCORPUS or not os.path.exists(self.TESTTARGETCORPUS):
                    print >>sys.stderr, bold(red("No predefined default test corpus set for reference (set TESTSOURCECORPUS and TESTTARGETCORPUS) or specify on command line "))
                    sys.exit(2)                    
                inputfile = self.TESTSOURCECORPUS                
                referencefile = self.TESTTARGETCORPUS                
            else:
                self.usage()
                sys.exit(2)
            
            if not self.score(inputfile, referencefile,'output.txt'): 
                sys.exit(1)
                
        elif cmd == 'help' or cmd == '-h':
            self.usage()
        else:
            print >>sys.stderr,"Error, no such command: " + cmd
            self.usage()
            sys.exit(2)
        
        sys.exit(0)
            
    def clean(self, targets):            
        if not targets:
            print >>sys.stderr,"Nothing to clean, please specify one or more targets: all, giza, moses, colibri, srilm"
            sys.exit(2)        
        
        if 'giza' in targets or 'all' in targets:
            self.cleanfiles('*.final', '*.vcb','*.snt','*.classes','*.classes.cats','*.gizacfg','*.Decoder.config','*.perp','*.cooc')
        if 'moses' in targets or 'all' in targets:
            self.cleanfiles('*.bal', '*.symal','*.s2t','*.s2t.sorted','*.t2s','*.t2s','*.sorted','*.phrasetable', '*.phraseextract', '*.phraseextract.inv','*.half','moses.ini')
        if 'srilm' in targets or 'all' in targets:
            self.cleanfiles('*.srilm')
        if 'colibri' in targets or 'all' in targets:
            self.cleanfiles('*.colibri')
        if 'test' in targets or 'all' in targets:
            self.cleanfiles('output.txt','*.score')
        if 'score' in targets or 'all' in targets:
            self.cleanfiles('*.score')
        return True
            
    def cleanfiles(self, *args):
        ok = True
        for mask in args:
            for filename in glob.glob(self.WORKDIR + '/' + mask):
                try:
                    os.unlink(filename)
                    print >>sys.stderr, green("Removed " + filename)
                except:
                    ok = False
                    print >>sys.stderr, bold(red("Unable to remove " + filename))
        return ok
            
    def starttrain(self):                
        self.init()
        if not self.check_common(): return False
        if not self.check_train(): return False
        
        if self.TOKENIZE_SOURCECORPUS and not self.tokenize_sourcecorpus(): return False
        if self.TOKENIZE_TARGETCORPUS and not self.tokenize_targetcorpus(): return False

        if self.BUILD_COLIBRI_ALIGNMENT and not self.build_colibri_alignment(): return False

        if self.BUILD_SRILM_TARGETMODEL and not self.build_srilm_targetmodel(): return False    
        if self.BUILD_SRILM_SOURCEMODEL and not self.build_srilm_sourcemodel(): return False       
                
        if self.BUILD_GIZA_WORDALIGNMENT and not self.build_giza_wordalignment(): return False
        if self.BUILD_GIZA_WORDALIGNMENT_REV and not self.build_giza_wordalignment_rev(): return False    
        
        if self.BUILD_MOSES_SYMAL and not self.build_moses_symal(): return False
        if self.BUILD_MOSES_WORDTRANSTABLE and not self.build_moses_wordtranstable(): return False
        if self.BUILD_MOSES_PHRASEEXTRACT and not self.build_moses_phraseextract(): return False
        if self.BUILD_MOSES_PHRASETRANSTABLE and not self.build_moses_phrasescore(): return False
        
        #TODO: Moses reordering model and generation model
        
        if self.BUILD_MOSES and not self.build_moses(): return False
        if self.BUILD_MOSES_MERT and not self.build_moses_mert(): return False
        
        if self.BUILD_PBMBMT and not self.build_pbmbmt(): return False
        
        return True    


    def header(self,name,*outputfiles, **kwargs):
        print >>sys.stderr, "----------------------------------------------------"
        if outputfiles:
            skip = True
            for outputfile in outputfiles:
                if not os.path.exists(outputfile):                
                    skip = False
                    break                                
            if skip:
                print >>sys.stderr, bold(yellow("Skipping " + name))  + " (output files already present)"
                return False        
        if 'cmd' in kwargs:
            print >>sys.stderr, bold(white("Calling " + name)) + ": " + kwargs['cmd']            
        else:
            print >>sys.stderr, bold(white("Calling " + name))
        return True
            
            
    def footer(self, name, r, *outputfiles, **kwargs):
        if 'successcodes' in kwargs:
            successcodes = kwargs['successcodes']
        else:
            successcodes = [0]
        if r in successcodes:
           print >>sys.stderr, bold(green("Finished " + name))
        else:
           print >>sys.stderr, bold(red("Runtime error from " + name + '(return code ' + str(r) + ')'))
           return False
        if outputfiles:
            error = False
            for outputfile in outputfiles:
                if os.path.exists(outputfile):                
                    print >>sys.stderr, green("Produced output file " + outputfile)
                else:
                    print >>sys.stderr, bold(red("Expected output file " + outputfile+ ", not produced!"))
                    error = True
            if error: 
                return False    
        return True            
        
    def runcmd(self, cmd, name, *outputfiles, **kwargs):        
        if not self.header(name,*outputfiles, cmd=cmd): return True  
        r = subprocess.call(cmd, shell=True)
        return self.footer(name, r, *outputfiles,**kwargs)
        
    def init(self):
        if not os.path.exists(self.getsourcefilename('txt')):
            try:
                os.symlink(self.TRAINSOURCECORPUS, self.getsourcefilename('txt') )
            except:
                pass
        if not os.path.exists(self.gettargetfilename('txt')):
            try:
                os.symlink(self.TRAINTARGETCORPUS, self.gettargetfilename('txt') )
            except:
                pass
        return True        
        
        
    #---------------------------------- Methods for building sub-parts ----------------------------
    
    def build_colibri_alignment(self):
        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.getsourcefilename('txt'), "Encoding source corpus for Colibri",self.getsourcefilename('cls'), self.getsourcefilename('clsenc') ): return False
         
        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.gettargetfilename('txt'), "Encoding target corpus for Colibri",self.gettargetfilename('cls'), self.gettargetfilename('clsenc') ): return False
        
        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.getsourcefilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS, "Building source-side pattern model",self.getsourcefilename('indexedpatternmodel.colibri') ): return False
        
        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.gettargetfilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS, "Building target-side pattern model",self.gettargetfilename('indexedpatternmodel.colibri') ): return False
        
        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.getsourcefilename('indexedpatternmodel.colibri') + ' ' + self.COLIBRI_GRAPHMODEL_OPTIONS, "Building source-side graph model",self.getsourcefilename('graphmodel.colibri') ): return False

        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.gettargetfilename('indexedpatternmodel.colibri') + ' ' + self.COLIBRI_GRAPHMODEL_OPTIONS, "Building target-side graph model",self.gettargetfilename('graphmodel.colibri') ): return False
        
        if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -s ' + self.getsourcefilename('graphmodel.colibri') + ' -t ' + self.gettargetfilename('graphmodel.colibri') + ' '+ self.COLIBRI_ALIGNER_OPTIONS, "Building alignment model",self.gettargetfilename('alignmodel.colibri') ): return False
        
        return True
        
    def build_giza_wordalignment(self):
        if not self.runcmd(self.EXEC_GIZA_PLAIN2SNT + ' ' + self.getsourcefilename('txt') + ' ' + self.gettargetfilename('txt'),'GIZA++ Input Preparation', self.getsourcefilename('vcb'), self.gettargetfilename('vcb'), self.gets2tfilename('snt',longform=True) ): return False
        if self.BUILD_GIZA_WORDALIGNMENT_COOC and not self.runcmd(self.EXEC_GIZA_SNT2COOC + ' ' + self.gettargetfilename('vcb') + ' ' + self.getsourcefilename('vcb') + ' ' + self.getsourcefilename('txt') + ' > ' + self.gets2tfilename('cooc'), 'GIZA++ Co-occurrence output',  self.gets2tfilename('cooc')): return False
        if not self.runcmd(self.EXEC_GIZA_MKCLS + ' -p' + self.getsourcefilename('txt') + ' ' + self.MKCLS_OPTIONS + ' -V' + self.getsourcefilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for source corpus', self.getsourcefilename('vcb.classes')): return False
        if not self.runcmd(self.EXEC_GIZA_MKCLS + ' -p' + self.gettargetfilename('txt') +  ' ' + self.MKCLS_OPTIONS  + ' -V' + self.gettargetfilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for target corpus', self.gettargetfilename('vcb.classes')): return False       
        if not self.runcmd(self.EXEC_GIZA + ' -S ' + self.getsourcefilename('vcb') + ' -T ' + self.gettargetfilename('vcb') + ' -C ' +  self.gets2tfilename('snt',longform=True) + ' ' + self.GIZA_OPTIONS + ' -o ' + self.gets2tfilename(),'GIZA++ Word Alignment',  self.gets2tfilename('A3.final') ): return False
        return True        

    def build_giza_wordalignment_rev(self):
        if not self.runcmd(self.EXEC_GIZA_PLAIN2SNT + ' ' + self.gettargetfilename('txt') + ' ' + self.getsourcefilename('txt'),'GIZA++ Input Preparation (reversed))', self.getsourcefilename('vcb'), self.gettargetfilename('vcb'), self.gett2sfilename('snt',longform=True) ): return False
        if self.BUILD_GIZA_WORDALIGNMENT_COOC and not self.runcmd(self.EXEC_GIZA_SNT2COOC + ' ' + self.getsourcefilename('vcb') + ' ' + self.gettargetfilename('vcb') + ' ' + self.gettargetfilename('txt') + ' > ' + self.gett2sfilename('cooc'), 'GIZA++ Co-occurrence output (reversed)',  self.gett2sfilename('cooc')): return False        
        if not self.runcmd(self.EXEC_GIZA_MKCLS + ' -p' + self.getsourcefilename('txt') + ' ' + self.MKCLS_OPTIONS + ' -V' + self.getsourcefilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for source corpus (reversed)', self.getsourcefilename('vcb.classes')): return False
        if not self.runcmd(self.EXEC_GIZA_MKCLS + ' -p' + self.gettargetfilename('txt') +  ' ' + self.MKCLS_OPTIONS  + ' -V' + self.gettargetfilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for target corpus (reversed)', self.gettargetfilename('vcb.classes')): return False       
        if not self.runcmd(self.EXEC_GIZA + ' -S ' + self.gettargetfilename('vcb') + ' -T ' + self.getsourcefilename('vcb') + ' -C ' +  self.gett2sfilename('snt',longform=True) + ' ' + self.GIZA_OPTIONS + ' -o ' + self.gett2sfilename(),'GIZA++ Word Alignment (reversed)',  self.gett2sfilename('A3.final') ): return False
        return True        

    def build_srilm_targetmodel(self):
        if not self.runcmd(self.EXEC_SRILM + ' -order ' + str(self.SRILM_ORDER) + ' ' + self.SRILM_OPTIONS + ' -text ' + self.gettargetfilename('txt') + ' -lm ' + self.gettargetfilename('srilm'),'SRILM Target-language Model', self.gettargetfilename('srilm')): return False
        return True
        
    def build_srilm_sourcemodel(self):
        if not self.runcmd(self.EXEC_SRILM +' -order ' + str(self.SRILM_ORDER) + ' ' + self.SRILM_OPTIONS + ' -text ' + self.getsourcefilename('txt') + ' -lm ' + self.getsourcefilename('srilm'),'SRILM Source-language Model', self.getsourcefilename('srilm')): return False
        return True        

    def tokenize_sourcecorpus(self):
        if not os.path.exists(self.getsourcefilename('notok')):
            os.rename( os.path.exists(self.getsourcefilename('txt')), os.path.exists(self.getsourcefilename('notok') ) )
        if not self.runcmd(self.EXEC_UCTO + ' ' + self.UCTO_OPTIONS +  ' -L' + self.SOURCELANG +  ' ' + self.getsourcefilename('notok') + ' ' + self.getsourcefilename('tok'),'Tokenisation Source Corpus', self.getsourcefilename('tok')): return False    
        if os.path.exists(self.getsourcefilename('tok')):
            try:
                os.unlink(self.getsourcefilename('txt'))
            except:
                pass
            os.symlink( self.getsourcefilename('tok'), self.getsourcefilename('txt') )
        return True
            
    def tokenize_targetcorpus(self):
        if not os.path.exists(self.gettargetfilename('notok')):
            os.rename( os.path.exists(self.gettargetfilename('txt')), os.path.exists(self.gettargetfilename('notok') ) )
        if not self.runcmd(self.EXEC_UCTO + ' ' + self.UCTO_OPTIONS +  ' -L' + self.SOURCELANG +  ' ' + self.gettargetfilename('notok') + ' ' + self.gettargetfilename('tok'),'Tokenisation Target Corpus', self.gettargetfilename('tok')): return False    
        if os.path.exists(self.gettargetfilename('tok')):
            try:
                os.unlink(self.gettargetfilename('txt'))
            except:
                pass
            os.symlink( self.gettargetfilename('tok'), self.gettargetfilename('txt') )
        return True
        

    def build_moses_symal(self):
        if not self.runcmd(self.EXEC_MOSES_GIZA2BAL + ' -d ' + self.gets2tfilename('A3.final')  + ' -i ' + self.gett2sfilename('A3.final') + ' > ' + self.gets2tfilename('bal'),'Data preparation for Symmetric Aligner', self.gets2tfilename('bal')): return False
        if not self.runcmd(self.EXEC_MOSES_SYMAL + ' ' + self.SYMAL_OPTIONS + ' < ' + self.gets2tfilename('bal') + ' > '  + self.gets2tfilename('symal'), 'Moses Symmetric Alignment',self.gets2tfilename('symal')): return False 
        return True        
    
    def build_moses_wordtranstable(self):
        if not self.runcmd(self.EXEC_MOSES_WORDTRANSTABLE + ' ' + self.getsourcefilename('txt') + ' ' + self.gettargetfilename('txt') + ' ' + self.gets2tfilename('symal')  + ' ' + self.gets2tfilename(),'Build of Lexical Translation Table', self.gets2tfilename('s2t'), self.gets2tfilename('t2s')): return False 
        return True
    
    def build_moses_phraseextract(self):        
        #TODO: Support for EPPEX and Hierchical rule extraction and reordering models
        cmd = self.EXEC_MOSES_PHRASEEXTRACT + ' ' + self.gettargetfilename('txt') + ' ' + self.getsourcefilename('txt') + ' ' + self.gets2tfilename('symal') + ' ' + self.gets2tfilename('phraseextract') + ' ' + str(self.PHRASEEXTRACT_MAX_PHRASE_LENGTH)
        if self.PHRASEEXTRACT_REORDERING_FLAGS: cmd += ' orientation ' + self.PHRASEEXTRACT_REORDERING_FLAGS
        if not self.runcmd(cmd ,'Phrase Extraction', self.gets2tfilename('phraseextract')): return False 
        return True        
    
    def build_moses_phrasescore(self):
        
        
        #if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('s2t') + ' > ' +  self.gets2tfilename('s2t.sorted'),'Sorting Lexical Translation Table (source->target)',  self.gets2tfilename('s2t.sorted') ): return False
        #if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('t2s') + ' > ' +  self.gets2tfilename('t2s.sorted'),'Sorting Lexical Translation Table (target->source)',  self.gets2tfilename('t2s.sorted') ): return False
        
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('phraseextract') + ' > ' +  self.gets2tfilename('phraseextract.sorted'),'Sorting Extracted Phrases (source->target)',  self.gets2tfilename('phraseextract.sorted') ): return False
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('phraseextract.inv') + ' > ' +  self.gets2tfilename('phraseextract.inv.sorted'),'Sorting Extracted Phrases Table (target->source)',  self.gets2tfilename('phraseextract.inv.sorted') ): return False
                

        if self.BUILD_MOSES_MEMSCORE:
            if not self.runcmd(self.EXEC_MOSES_MEMSCORE + ' -s ' + self.MOSES_MEMSCORE_METHOD + ' -s lexweights ' + self.gets2tfilename('t2s') + ' -r ' +  self.MOSES_MEMSCORE_METHOD + ' -r lexweights ' + self.gets2tfilename('s2t') + ' -s const 2.718 < ' + self.gets2tfilename('phraseextract') +  ' > ' + self.gets2tfilename('phrasetable')): return False            
        else:
            if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gets2tfilename('phraseextract.sorted') + ' ' +  self.gets2tfilename('s2t') + ' ' + self.gets2tfilename('half.s2t') + ' ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (source->target)', self.gets2tfilename('half.s2t') ): return False        
            if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gets2tfilename('phraseextract.inv.sorted') + ' '  + self.gets2tfilename('t2s') + ' ' + self.gets2tfilename('half.t2s') + ' --Inverse ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (target->source)', self.gets2tfilename('half.t2s') ): return False        
            
            if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('half.t2s') + ' > ' +  self.gets2tfilename('half.t2s.sorted'),'Sorting Inverse Table',  self.gets2tfilename('half.t2s.sorted') ): return False
                    
            if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE + ' ' + self.gets2tfilename('half.s2t') + ' ' + self.gets2tfilename('half.t2s.sorted') + ' ' + self.gets2tfilename('phrasetable') + ' ' + self.PHRASESCORE_OPTIONS, 'Consolidating two phrase table halves', self.gets2tfilename('phrasetable') ): return False
        
        return True
        
    def build_moses_reorderingmodel(self):
        #TODO, skipped for now
        return False
    
    def build_moses_generationmodel(self):
        #TODO, skipped for now
        return False
        
    def build_moses(self):
        outputfiles = ['moses.ini']
        if not self.header('Build Moses Configuration',*outputfiles): return True
        f = open('moses.ini','w')
        f.write('#Moses INI, produced by mtwrapper.py\n')
        f.write('[input-factors]\n')
        f.write('0\n\n')
        f.write('[mapping]\n')
        f.write('T 0\n\n') 
        f.write('# translation tables: source-factors, target-factors, number of scores, file\n')
        f.write('[ttable-file]\n')
        if self.BUILD_MOSES_PHRASETRANSTABLE:
            f.write('0 0 0 5 ' + self.gets2tfilename('phrasetable') + '\n\n')
        f.write('[lmodel-file]\n')
        if self.BUILD_SRILM_TARGETMODEL:
            f.write('0 0 ' + str(self.SRILM_ORDER) + ' ' + self.gettargetfilename('srilm') + '\n\n')
        f.write('[ttable-limit]\n20\n\n')
        f.write('[weight-d]\n1\n\n')
        f.write('[weight-l]\n1\n\n')
        f.write('[weight-t]\n1\n1\n1\n1\n1\n\n')
        f.write('[weight-w]\n0\n\n')        
        f.close()
        return self.footer('Build Moses Configuration', 0, *outputfiles)
    
    def build_moses_mert(self):            
        if not self.runcmd(self.EXEC_MOSES_MERT + ' --mertdir=' + self.PATH_MOSES_MERT + ' ' + self.MOSES_MERT_OPTIONS + ' ' + self.DEVSOURCECORPUS + ' ' + self.DEVTARGETCORPUS + ' ' + self.EXEC_MOSES  + ' moses.ini', 'Parameter tuning for Moses using MERT'): return False         
        return True
    
    def build_pbmbmt(self):
        outputfiles = []
        for i in range(1,self.PBMBMT_MAXPHRASELENGTH+1):
            outputfiles.append( self.gets2tfilename('train.' + str(self.PBMBMT_LEFTCONTEXTSIZE) + str(i) + str(self.PBMBMT_RIGHTCONTEXTSIZE) + '.0x0.inst') )
        
        if not self.runcmd(self.EXEC_PBMBMT_INSTANCEGENERATOR + ' --train=' +  self.gets2tfilename('A3.final') + ' -p ' + self.gets2tfilename('phrasetable') + ' --nfeatleft=' + str(self.PBMBMT_LEFTCONTEXTSIZE) + ' --nfeatright='+str(self.PBMBMT_RIGHTCONTEXTSIZE) + ' --phraselength='+str(self.PBMBMT_MAXPHRASELENGTH),'Extracting Training Instances for PBMBMT', *outputfiles ): return False
        
                
        outputfiles_ibase = [ x  + '.ibase' for x in outputfiles ]                
        for trainfile in glob.glob(self.WORKDIR + '/*.train.*.inst'):
            if not self.runcmd(self.EXEC_TIMBL + self.PBMBMT_TIMBL_OPTIONS + ' +v+db+di +D -s -f ' + trainfile + ' -I ' + trainfile + '.ibase', 'Building classifier ' + os.path.basename(trainfile), *outputfiles_ibase): return False
        
         
        return True
    
    def run(self, inputfile, outputfile='output.txt', tokenise=False):        
        if tokenise and (not self.EXEC_UCTO or not os.path.isfile(self.EXEC_UCTO)):
            print >>sys.stderr,red("Error: Ucto not found! Unable to tokenise!" )
            return False
        
        if not self.check_common(): return False
        if not self.check_run(): return False
        
        if not os.path.isfile(inputfile):
            print >>sys.stderr,red("Error: Input file " + inputfile + " not found!" )
            return False

        if os.path.exists( outputfile):
            os.unlink( outputfile)

        if tokenise:        
            if not self.runcmd(self.EXEC_UCTO + ' -n -L' + self.SOURCELANG +  ' ' + inputfile + ' ' + 'input.txt','Tokenisation of Input File'): return False                                        
        else:
            if os.path.exists( self.WORKDIR + '/input.txt'):
                os.unlink( self.WORKDIR + '/input.txt' )
            os.symlink(inputfile, self.WORKDIR + '/input.txt' )
        
        if self.BUILD_MOSES and not self.run_moses(): return False
        if self.BUILD_PBMBMT and not self.run_pbmbmt(): return False
        
        
        
        os.rename('output.txt',outputfile)        
        return True


    
    def run_moses(self):
        if not self.runcmd(self.EXEC_MOSES + ' -f moses.ini < input.txt > output.txt','Moses Decoder'): return False
        return True 
    
    def run_pbmbmt(self):
        for trainfile in glob.glob(self.WORKDIR + '/*.train.*.inst'):
            basename = '.'.join(os.path.basename(trainfile).split('.')[:-4])
        
        os.symlink(self.WORKDIR + '/input.txt', self.WORKDIR + '/' + basename + '.txt')  

        if not self.runcmd(self.EXEC_PBMBMT_INSTANCEGENERATOR + ' --test='+self.WORKDIR + '/' +basename + '.txt -p ' + self.gets2tfilename('phrasetable') + ' --nfeatleft=' + str(self.PBMBMT_LEFTCONTEXTSIZE) + ' --nfeatright='+str(self.PBMBMT_RIGHTCONTEXTSIZE) + ' --phraselength='+str(self.PBMBMT_MAXPHRASELENGTH),'Extracting Test Instances for PBMBMT'): return False
                    
        for testfile in glob.glob(self.WORKDIR + '/*.test.*.inst'):
            if not self.runcmd(self.EXEC_TIMBL + self.PBMBMT_TIMBL_OPTIONS + ' +v+db+di +D -s -f ' + testfile + ' -i ' + testfile.replace('test','train') + '.ibase', 'Running classifier ' + os.path.basename(testfile)): return False
        
        if not self.runcmd(self.EXEC_PBMBMT_DECODER + ' -t ' + basename + '.txt --srilm=' + self.gets2tfilename('srilm') + ' ' + self.PBMBMT_DECODER_OPTIONS + ' > ' + self.WORKDIR + '/output.txt','Running PBMBMT Decoder'): return False        
        
        return True
            

    
    def score(self, sourcefile, reffile, targetfile):
        if not os.path.isfile(targetfile):
            print >>sys.stderr,red("Error: Output file " + targetfile + " not found!" )
            return False    
     
        self.header('Converting source to XML for evaluation')
        r = self.xmlize(sourcefile,'src')
        sourcexml = sourcefile + '.xml'
        sourcexml = sourcexml.replace('//','/')
        if not self.footer('Converting source to XML for evaluation',int(not r), sourcexml): return False
        self.header('Converting reference to XML for evaluation')
        r = self.xmlize(reffile,'ref')
        refxml = reffile + '.xml'
        refxml = refxml.replace('//','/')
        if not self.footer('Converting reference to XML for evaluation',int(not r),refxml): return False
        self.header('Converting output to XML for evaluation')
        r = self.xmlize(targetfile,'tst')
        targetxml = targetfile + '.xml'
        targetxml = targetxml.replace('//','/')
        if not self.footer('Converting output to XML for evaluation',int(not r),targetxml): return False        
        
        
        self.per = 0
        self.wer = 0
        self.bleu = 0
        self.meteor = 0
        self.nist = 0
        self.ter = 0        
        
        
        errors = False
        if self.EXEC_MATREX_BLEU and os.path.exists(self.EXEC_MATREX_BLEU) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_BLEU + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml + ' -ci > ' + 'bleu.score',  'Computing BLEU score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/bleu.score')
                    for line in f:
                        if line[0:9] == "BLEUr1n4,":
                             self.bleu = float(line[10:].strip())
                             print >>sys.stderr,"BLEU score: ", self.bleu
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading bleu.score")
                    errors = True            
        else:
            print >>sys.stderr, yellow("Skipping BLEU (no script found ["+self.EXEC_MATREX_BLEU+"])")
            
        if self.EXEC_MATREX_WER and os.path.exists(self.EXEC_MATREX_WER) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_WER + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml + '  > ' + 'wer.score', 'Computing WER score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/wer.score')
                    for line in f:
                        if line[0:11] == "WER score =":
                             self.wer = float(line[12:20].strip())
                             print >>sys.stderr,"WER score: ", self.wer
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading wer.score")
                    errors = True     
        else:
            print >>sys.stderr, yellow("Skipping WER (no script found ["+self.EXEC_MATREX_WER+"]) ")
     
        if self.EXEC_MATREX_PER and os.path.exists(self.EXEC_MATREX_PER) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_PER + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml + '  > ' + 'per.score',  'Computing PER score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/per.score')
                    for line in f:
                        if line[0:11] == "PER score =":
                             self.per = float(line[12:20].strip())
                             print >>sys.stderr,"PER score: ", self.per
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading per.score")
                    errors = True                     
        else:
            print >>sys.stderr, yellow("Skipping PER (no script found ["+self.EXEC_MATREX_PER+"])")
        
        if self.EXEC_MATREX_METEOR and os.path.exists(self.EXEC_MATREX_METEOR) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' -I ' + os.path.dirname(self.EXEC_MATREX_METEOR) + ' ' + self.EXEC_MATREX_METEOR + " -s " + self.CORPUSNAME + " -r " + refxml + ' -t ' + targetxml + ' --modules "exact"  > ' + 'meteor.score',  'Computing METEOR score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/meteor.score')
                    for line in f:
                        if line[0:6] == "Score:":
                             self.meteor = float(line[7:].strip())
                             print >>sys.stderr,"METEOR score: ", self.meteor
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading meteor.score")
                    errors = True                      
        else:
            print >>sys.stderr, yellow("Skipping METEOR (no script found ["+self.EXEC_MATREX_METEOR+"])")

        if self.EXEC_MATREX_MTEVAL and os.path.exists(self.EXEC_MATREX_MTEVAL) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_MTEVAL + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml +  '  > ' + 'mteval.score',  'Computing NIST & BLEU scores'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/mteval.score')
                    for line in f:
                        if line[0:12] == "NIST score =":
                            self.nist = float(line[13:21].strip())
                            print >>sys.stderr,"NIST score: ", self.blist
                        if line[21:33] == "BLEU score =":
                            if self.bleu > 0:
                                self.bleu = float(line[34:40].strip())
                                print >>sys.stderr,"BLEU score: ", self.bleu
                            else:
                                print >>sys.stderr,"BLEU score (not stored): ", float(line[10:].strip())
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading mteval.score")
                    errors = True                   
        else:
            print >>sys.stderr, yellow("Skipping MTEVAL (BLEU & NIST) (no script found)")
     
        if self.EXEC_MATREX_TER and os.path.exists(self.EXEC_MATREX_TER) and self.EXEC_JAVA and os.path.exists(self.EXEC_JAVA):
            if not self.runcmd(self.EXEC_JAVA + ' -jar ' + self.EXEC_MATREX_TER + " -r " + refxml + ' -h ' + targetxml + '  > ' + 'ter.score',  'Computing TER score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/ter.score')
                    for line in f:
                        if line[0:10] == "Total TER:":
                             self.ter = float(line[11:].strip().split(' ')[0])
                             print >>sys.stderr,"TER score: ", self.ter
                    f.close()
                except:                
                    print >>sys.stderr, red("Error reading ter.score")
        else:
            print >>sys.stderr, yellow("Skipping TER (no script found)")     
    
        if not errors:
            print >>sys.stderr,"SCORE SUMMARY\n===================\n"
            f = open(self.WORKDIR + '/summary.score','w')
            s = "BLEU METEOR NIST TER WER PER"
            f.write(s+ "\n")
            print >>sys.stderr, s            
            s = str(round(self.bleu,4)) + " " + str(round(self.meteor,4)) + " " + str(round(self.nist,4))  + " " + str(round(self.ter,2)) + " " + str(round(self.wer,2))  + " " + str(round(self.per,2))
            f.write(s + "\n")
            print >>sys.stderr, s
            f.close()
                        
        elif os.path.exists(self.WORKDIR + '/summary.score'):
            os.unlink(self.WORKDIR + '/summary.score')
                             
        return not errors
    
    def test(self, sourcefile, reffile):
        if not self.check_common(): return False
        if not self.check_run(): return False
        if not self.check_test(): return False
        
        if not os.path.isfile(sourcefile):
            print >>sys.stderr,red("Error: Source file " + sourcefile + " not found!" )
            return False        
             
        if not os.path.isfile(reffile):
            print >>sys.stderr,red("Error: Reference file " + reffile + " not found!" )
            return False        
        
        if not self.run(sourcefile):
            return False
    
        targetfile = 'output.txt'
        if not self.score(sourcefile,reffile,targetfile):
            return False
        
        return True
    
    
    def xmlize(self, inputfile, type):
        assert type in ('tst','ref','src')
        try:        
            fin = open(inputfile,'r')
            fout = open(inputfile + '.xml','w')
            fout.write( "<" + type + "set setid=\"mteval\" srclang=\"" + self.SOURCELANG + "\" trglang=\"" + self.TARGETLANG + "\">\n")
            fout.write("<DOC docid=\"" + self.CORPUSNAME + "\" sysid=\"" + self.CORPUSNAME + "\">\n")
            for linenum, line in enumerate(fin):
                fout.write("<seg id=\"" + str(linenum+1) + "\">\n")
                fout.write(line.strip())
                fout.write("</seg>\n")    
            fout.write("</DOC>\n</" + type + "set>")       
            fout.close()
            fin.close()
        except IOError:
            return False
        return True

   



    
def usage():
    print >>sys.stderr,"mtwrapper.py -- MT wrapper - Outputs a MT wrapper script (python)"
    print >>sys.stderr,"Mandatory Input:"
    print >>sys.stderr,"\t-n <name>         Name of the corpus [MANDATORY!]"
    print >>sys.stderr,"\t-s <file>         Corpus in source language (for training)"
    print >>sys.stderr,"\t-t <file>         Corpus in target language (for training)"
    print >>sys.stderr,"\t-S <code>         Source language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-T <code>         Target language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-w <dir>          Work directory (by default: current dir)"
    print >>sys.stderr,"Optional Input:"
    print >>sys.stderr,"\t--testset=n          Extract a random sample of n lines as test set, and exclude from training"
    print >>sys.stderr,"\t--devset=n           Extract a random sample of n lines as development set, and exclude from training"
    print >>sys.stderr,"\t--trainset=n         Restrict the training set to a random sample of n lines"
    print >>sys.stderr,"\t-i <dirs>         Colon-separated directories where python can find non-standard modules"

if __name__ == "__main__":        
    
    
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:t:S:T:n:x:w:d:i:", ['testset=','devset=','trainset='])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    sourcecorpusfile = targetcorpusfile = sourcelang = targetlang = corpusname = expname = workdir = ""
    devsourcecorpusfile = devtargetcorpusfile = testsourcecorpusfile = testtargetcorpusfile = ""
    trainset = testset = devset = 0
    includedirs = []

    for o, a in opts:
        if o == '-s':
            sourcecorpusfile = a
        elif o == '-t':
            targetcorpusfile = a
        elif o == '-S':
            sourcelang = a        
        elif o == '-T':
            targetlang = a        
        elif o == '-n':
            corpusname = a
        elif o == '-x':
            expname = a
        elif o == '-w':
            workdir = a
        elif o == '-d':
            parentdir = a
        elif o == '-i':
            includedirs = a.split(':')
        elif o == '--testset':
            testset = int(a)
        elif o == '--devset':
            devset = int(a)            
        elif o == '--trainset':
            trainset = int(a)
        elif o == '-h':            
            usage()
            sys.exit(0)
        else:
            print >>sys.stderr,"Error, invalid option: ", o
            usage()
            sys.exit(2)
    
    if not corpusname or not sourcelang or not targetlang:
        print >>sys.stderr,"Specify at least -n -S and -T"
        usage()
        sys.exit(2)
        
    if workdir and not os.path.isdir(workdir):            
        print>>sys.stderr, "Creating work directory " + workdir
        os.mkdir(workdir)
    elif workdir:
        print>>sys.stderr, yellow("WARNING: work directory " +  workdir + " already exists! Press ENTER to continue or ctrl-C to abort")
        raw_input()
    elif parentdir:
        if not os.path.isdir(parentdir):
            print>>sys.stderr, "Creating parent directory " + parentdir
            os.mkdir(parentdir)
        workdir = parentdir + '/' + corpusname + '-' + sourcelang + '-' + targetlang
        if expname: workdir += '-' + expname
        if workdir and not os.path.isdir(workdir):            
            print>>sys.stderr, "Creating work directory " + workdir
            os.mkdir(workdir)
        elif workdir:
            print>>sys.stderr, yellow("WARNING: work directory " +  workdir + " already exists! Press ENTER to continue or ctrl-C to abort")
            raw_input()
    else:
        workdir = os.getcwd()        
        
    if testset or devset:
        if not sourcecorpusfile or not targetcorpusfile:
            print>>sys.stderr, "Error: You need to specify -s and -t on the command line!"
        filesampler([sourcecorpusfile, targetcorpusfile],  testset, devset, trainset, workdir )
        
        #rename files
        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.dev' 
        if os.path.exists(oldfile): 
            os.rename(oldfile, workdir + '/' + corpusname + '-' + sourcelang + '-dev.txt')
            devsourcecorpusfile = workdir + '/' + corpusname + '-' + sourcelang + '-dev.txt'
        
        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.dev' 
        if os.path.exists(oldfile): 
            os.rename(oldfile, workdir + '/'+ corpusname + '-' + targetlang + '-dev.txt')
            devtargetcorpusfile = workdir + '/' +corpusname + '-' + targetlang + '-dev.txt'

        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.test' 
        if os.path.exists(oldfile): 
            os.rename(oldfile, workdir + '/' +corpusname + '-' + sourcelang + '-test.txt')
            testsourcecorpusfile = workdir + '/' +corpusname + '-' + sourcelang + '-test.txt'
        
        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.test' 
        if os.path.exists(oldfile): 
            os.rename(oldfile, workdir + '/' + corpusname + '-' + targetlang + '-test.txt')
            testtargetcorpusfile = workdir + '/' + corpusname + '-' + targetlang + '-test.txt'
        
        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.train' 
        if os.path.exists(oldfile): 
            os.rename(oldfile, workdir + '/' + corpusname + '-' + sourcelang + '-train.txt')
            sourcecorpusfile = workdir + '/' + corpusname + '-' + sourcelang + '-train.txt'
        
        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.train' 
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' + corpusname + '-' + targetlang + '-train.txt')
            targetcorpusfile = workdir + '/' + corpusname + '-' + targetlang + '-train.txt'
                
    if expname:
        settingsfile = workdir + '/mt-' + corpusname + '-' + sourcelang + '-' + targetlang + '-' + expname + '.py'
    else:
        settingsfile = workdir + '/mt-' + corpusname + '-' + sourcelang + '-' + targetlang + '.py'
    f = codecs.open(settingsfile,'w','utf-8')
    f.write("#! /usr/bin/env python\n# -*- coding: utf8 -*-#\n\n")
    f.write("#Generated with: " + ' '.join(sys.argv) + "\n\n")
    if includedirs:         
        f.write('import sys\n')
        for dir in includedirs:
            f.write("sys.path.append('" + dir + "')\n")
    f.write("\n")
    f.write("from mtwrapper import MTWrapper\n")    
    f.write("mtwrapper = MTWrapper(\n")
    for key, default, help in MTWrapper.defaults:            
        if key == 'CORPUSNAME': 
            default = corpusname
        if key == 'EXPERIMENTNAME': 
            default = expname
        elif key == 'WORKDIR':
            default = workdir
        elif key == 'TRAINSOURCECORPUS':
            default = sourcecorpusfile
        elif key == 'TRAINTARGETCORPUS':
            default = targetcorpusfile
        elif key == 'SOURCELANG':
            default = sourcelang
        elif key == 'TARGETLANG':
            default = targetlang
        elif key == 'TESTSOURCECORPUS':
            default = testsourcecorpusfile
        elif key == 'TESTTARGETCORPUS':
            default = testtargetcorpusfile            
        elif key == 'DEVSOURCECORPUS':
            default = devsourcecorpusfile
        elif key == 'DEVTARGETCORPUS':
            default = devtargetcorpusfile
            
        
        if isinstance(default, str) or isinstance(default,  unicode):            
            f.write("    " + key + "=\"" + default + "\"")
        else:
            f.write("    " + key + "=" + str(default))
            
        if help:
            f.write(", #" + help + "\n")
        else:
            f.write(",\n")
    f.write(")\nmtwrapper.start()\n")
    f.close()
    os.chmod(settingsfile, 0754)
    os.system('vim ' + settingsfile)    




