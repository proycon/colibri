#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import os
import subprocess
import getopt
import codecs


class MTWrapper(object):
    

    def __init__(self):        
        self.CORPUSNAME = ""
        self.WORKDIR = "./" 

        self.TRAINSOURCECORPUS = ""
        self.TRAINTARGETCORPUS = ""

        self.SOURCELANG = "" #language code
        self.TARGETLANG = "" #language code

        self.TOKENIZE_SOURCECORPUS = False #Set to true if the corpus is not tokenised yet
        self.TOKENIZE_TARGETCORPUS = False #Set to true if the corpus is not tokenised yet
        self.BUILD_SRILM_SOURCEMODEL = False
        self.BUILD_SRILM_TARGETMODEL = False
        self.BUILD_GIZA_WORDALIGNMENT = False
        self.BUILD_GIZA_WORDALIGNMENT_REV = False #Reverse word-alignment
        self.BUILD_GIZA_WORDALIGNMENT_COOC = False
        self.BUILD_MOSES_SYMAL = True #Symmetrise word alignments
        self.BUILD_MOSES_WORDTRANSTABLE = True #Extract word translation table
        self.BUILD_MOSES_PHRASEEXTRACT = True  #Phrase extraction
        self.BUILD_MOSES_PHRASETRANSTABLE = True  #Phrase scoring (this will make a phrase table!)
        #self.BUILD_MOSES_PHRASETRANSTABLE = False
        
        #defaults
        self.PATH_UCTO = self.findpath('ucto')      
        self.PATH_TIMBL = self.findpath('timbl')
        self.PATH_GIZA_MKCLS = self.findpath('mkcls')
        self.PATH_GIZA = self.findpath('GIZA++')
        self.PATH_GIZA_PLAIN2SNT = self.findpath('plain2snt.out')
        self.PATH_GIZA_SNT2COOC = self.findpath('snt2cooc.out')                
        self.PATH_MOSES = self.findpath('moses')
        self.PATH_MOSES_GIZA2BAL = self.findpath('scripts/giza2bal.pl')
        self.PATH_MOSES_SYMAL = self.findpath('symal')
        self.PATH_MOSES_WORDTRANSTABLE = self.findpath('scripts/moses-lexicalextractiontable.py')
        self.PATH_MOSES_PHRASEEXTRACT = self.findpath('scripts/training/phrase-extract/extract')
        self.PATH_MOSES_PHRASEEXTRACT_CONSOLIDATE = self.findpath('scripts/training/phrase-extract/consolidate')
        self.PATH_MOSES_PHRASEEXTRACT_SCORE = self.findpath('scripts/training/phrase-extract/score')        
        self.PATH_SRILM = self.findpath('ngram-count')   
        
        #default options
        self.MKCLS_OPTIONS = '-m2 -C50' 
        self.GIZA_OPTIONS = '-p0 0.98 -m1 5 -m2 0 -m3 3 -m4 3 -nsmooth 4 -model4smoothfactor 0.4'
        #self.GIZA_OPTIONS = '-p0 0.98 -m1 5 -m2 0 -m3 0 -m4 0 -nsmooth 4 -hmmiterations 5 -hmmdumpfrequency -5'
        self.SRILM_OPTIONS = '-order 3 -interpolate -kndiscount' #trigram model
        self.UCTO_OPTIONS = '-m -n'
        self.SYMAL_OPTIONS = "-alignment=grow -diagonal=yes -final=yes -both=no" #alignment: union/intersect/grow/srctotgt/tgttosrc ,  diagonal: yes|no , -final: yes|no , -both: yes|no 
        self.PHRASEEXTRACT_MAX_PHRASE_LENGTH = 7
        self.PHRASEEXTRACT_REORDERING_FLAGS = "" #default distance-based reordering
        #" --model wbe-mslr --model phrase-mslr --model hier-mslr" #Maximum lexical reordering

    def findpath(self, name):
        for path in os.environ['PATH'].split(':'):
            if os.path.exists(path + '/' + name) and not os.path.isdir(path + '/' + name):
                print >>sys.stderr, green("Found " + name + " in " + path)
                return path + '/' + name
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

    def check_train(self):
        sane = True                            
        if not self.TRAINSOURCECORPUS:
            print >>sys.stderr,red("Configuration error: TRAINSOURCECORPUS not specified!")
            sane = False
        if not self.TRAINTARGETCORPUS:
            print >>sys.stderr,red("Configuration error: TRAINTARGETCORPUS not specified!")
            sane = False
            
        if (self.TOKENIZE_SOURCECORPUS or self.TOKENIZE_TARGETCORPUS) and (not self.PATH_UCTO or not os.path.isfile(self.PATH_UCTO)):
            print >>sys.stderr,red("Dependency error: ucto not found (PATH_UCTO=" + self.PATH_UCTO + ")")
            sane = False
            
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

            if not self.PATH_MOSES_GIZA2BAL or not os.path.isfile(self.PATH_MOSES_GIZA2BAL):
                sane = False
                print >>sys.stderr,red("Dependency error: giza2bal.pl (provided by Moses) not found (PATH_MOSES_GIZA2BAL=" + self.PATH_MOSES_GIZA2BAL + ")")
               
            
            if not self.PATH_MOSES_SYMAL or not os.path.isfile(self.PATH_MOSES_SYMAL):
                sane = False
                print >>sys.stderr,red("Dependency error: symal (provided by Moses) not found (PATH_MOSES_SYMAL=" + self.PATH_MOSES_SYMAL + ")")
               
            if not self.PATH_MOSES or not os.path.isfile(self.PATH_MOSES):
                sane = False
                print >>sys.stderr,red("Dependency error: Moses not found (PATH_MOSES=" + self.PATH_MOSES + ")")
                
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.PATH_GIZA or not os.path.isfile(self.PATH_GIZA)): 
            print >>sys.stderr,red("Dependency error: GIZA++ not found (PATH_GIZA=" + self.PATH_GIZA + ")")
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.PATH_GIZA_PLAIN2SNT or not os.path.isfile(self.PATH_GIZA_PLAIN2SNT)): 
            print >>sys.stderr,red("Dependency error: plain2snt.out (provided by GIZA++) not found (PATH_GIZA_PLAIN2SNT=" + self.PATH_GIZA_PLAIN2SNT + ")")            
            sane = False
        if self.BUILD_GIZA_WORDALIGNMENT_COOC and (not self.PATH_GIZA_SNT2COOC or not os.path.isfile(self.PATH_GIZA_SNT2COOC)): 
            print >>sys.stderr,red("Dependency error: snt2cooc.out (provided by GIZA++) not found (PATH_GIZA_SNT2COOC=" + self.PATH_GIZA_SNT2COOC + ")")            
            sane = False                        
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.PATH_GIZA_MKCLS or not os.path.isfile(self.PATH_GIZA_MKCLS)): 
            print >>sys.stderr,red("Dependency error: mkcls (provided by GIZA++) not found (PATH_GIZA_MKCLS=" + self.PATH_GIZA_MKCLS + ")")            
            sane = False                            
        if (self.BUILD_SRILM_TARGETMODEL or self.BUILD_SRILM_SOURCEMODEL) and (not self.PATH_SRILM or not os.path.isfile(self.PATH_SRILM)):
            print >>sys.stderr,red("Dependency error: ngram-count (provided by SRILM) not found (PATH_SRILM=" + self.PATH_SRILM + ")")
            sane = False
        return sane
    
    

    def getsourcefilename(self, extension):
        return self.WORKDIR + self.CORPUSNAME + '-' + self.SOURCELANG + '.' + extension

    def gettargetfilename(self, extension):
        return self.WORKDIR + self.CORPUSNAME + '-' + self.TARGETLANG + '.' + extension
    
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
        print >>sys.stderr,"\ttrain                           Train the MT system"      
        #print >>sys.stderr,"\tclean [all|moses|giza|srilm]    Clean data"


    def start(self):        
        try:
            cmd = sys.argv[1]
        except:
            self.usage()
            sys.exit(2)
        if cmd == 'train':
            self.starttrain()
        elif cmd == 'help' or cmd == '-h':
            self.usage()
        else:
            print >>sys.stderr,"Error, no such command: " + cmd
            self.usage()
            sys.exit(2)
            
    def starttrain(self):                
        self.init()
        if not self.check_common(): return False
        if not self.check_train(): return False
        
        if self.TOKENIZE_SOURCECORPUS and not self.tokenize_sourcecorpus(): return False
        if self.TOKENIZE_TARGETCORPUS and not self.tokenize_targetcorpus(): return False


        if self.BUILD_SRILM_TARGETMODEL and not self.build_srilm_targetmodel(): return False    
        if self.BUILD_SRILM_SOURCEMODEL and not self.build_srilm_sourcemodel(): return False       
        
        if self.BUILD_GIZA_WORDALIGNMENT and not self.build_giza_wordalignment(): return False
        if self.BUILD_GIZA_WORDALIGNMENT_REV and not self.build_giza_wordalignment_rev(): return False    
        
        if self.BUILD_MOSES_SYMAL and not self.build_moses_symal(): return False
        if self.BUILD_MOSES_WORDTRANSTABLE and not self.build_moses_wordtranstable(): return False
        if self.BUILD_MOSES_PHRASEEXTRACT and not self.build_moses_phraseextract(): return False
        if self.BUILD_MOSES_PHRASETRANSTABLE and not self.build_moses_phrasescore(): return False
        return True    

        
    def runcmd(self, cmd, name, *outputfiles, **kwargs):
        if 'successcodes' in kwargs:
            successcodes = kwargs['successcodes']
        else:
            successcodes = [0]
        print >>sys.stderr, "----------------------------------------------------"
        if outputfiles:
            skip = True
            for outputfile in outputfiles:
                if not os.path.exists(outputfile):                
                    skip = False
                    break                                
            if skip:
                print >>sys.stderr, bold(yellow("Skipping " + name))  + " (output files already present)"
                return True        
        print >>sys.stderr, bold(white("Calling " + name)) + ": " + cmd        
        r = subprocess.call(cmd, shell=True)
        if r in successcodes:
           print >>sys.stderr, bold(green("Finished " + name))
        else:
           print >>sys.stderr, bold(red("Runtime error from " + name + '(return code ' + str(r) + ')'))
           return False
        return True
        
    def init(self):
        if not os.path.exists(self.getsourcefilename('txt')):
            os.symlink(self.TRAINSOURCECORPUS, self.getsourcefilename('txt') )
        if not os.path.exists(self.gettargetfilename('txt')):
            os.symlink(self.TRAINTARGETCORPUS, self.gettargetfilename('txt') )
        return True        
        
        
    #---------------------------------- Methods for building sub-parts ----------------------------
        
    def build_giza_wordalignment(self):
        if not self.runcmd(self.PATH_GIZA_PLAIN2SNT + ' ' + self.getsourcefilename('txt') + ' ' + self.gettargetfilename('txt'),'GIZA++ Input Preparation', self.getsourcefilename('vcb'), self.gettargetfilename('vcb'), self.gets2tfilename('snt',longform=True) ): return False
        if not self.runcmd(self.PATH_GIZA_SNT2COOC + ' ' + self.gettargetfilename('vcb') + ' ' + self.getsourcefilename('vcb') + ' ' + self.getsourcefilename('txt') + ' > ' + self.gets2tfilename('cooc'), 'GIZA++ Co-occurrence output',  self.gets2tfilename('cooc')): return False
        if not self.runcmd(self.PATH_GIZA_MKCLS + ' -p' + self.getsourcefilename('txt') + ' ' + self.MKCLS_OPTIONS + ' -V' + self.getsourcefilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for source corpus', self.getsourcefilename('vcb.classes')): return False
        if not self.runcmd(self.PATH_GIZA_MKCLS + ' -p' + self.gettargetfilename('txt') +  ' ' + self.MKCLS_OPTIONS  + ' -V' + self.gettargetfilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for target corpus', self.gettargetfilename('vcb.classes')): return False       
        if not self.runcmd(self.PATH_GIZA + ' -S ' + self.getsourcefilename('vcb') + ' -T ' + self.gettargetfilename('vcb') + ' -C ' +  self.gets2tfilename('snt',longform=True) + ' ' + self.GIZA_OPTIONS + ' -o ' + self.gets2tfilename(),'GIZA++ Word Alignment',  self.gets2tfilename('A3.final') ): return False
        return True        

    def build_giza_wordalignment_rev(self):
        if not self.runcmd(self.PATH_GIZA_PLAIN2SNT + ' ' + self.gettargetfilename('txt') + ' ' + self.getsourcefilename('txt'),'GIZA++ Input Preparation (reversed))', self.getsourcefilename('vcb'), self.gettargetfilename('vcb'), self.gett2sfilename('snt',longform=True) ): return False
        if not self.runcmd(self.PATH_GIZA_SNT2COOC + ' ' + self.getsourcefilename('vcb') + ' ' + self.gettargetfilename('vcb') + ' ' + self.gettargetfilename('txt') + ' > ' + self.gett2sfilename('cooc'), 'GIZA++ Co-occurrence output (reversed)',  self.gett2sfilename('cooc')): return False        
        if not self.runcmd(self.PATH_GIZA_MKCLS + ' -p' + self.getsourcefilename('txt') + ' ' + self.MKCLS_OPTIONS + ' -V' + self.getsourcefilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for source corpus (reversed)', self.getsourcefilename('vcb.classes')): return False
        if not self.runcmd(self.PATH_GIZA_MKCLS + ' -p' + self.gettargetfilename('txt') +  ' ' + self.MKCLS_OPTIONS  + ' -V' + self.gettargetfilename('vcb.classes') + ' opt','GIZA++ Word Categoriser for target corpus (reversed)', self.gettargetfilename('vcb.classes')): return False       
        if not self.runcmd(self.PATH_GIZA + ' -S ' + self.gettargetfilename('vcb') + ' -T ' + self.getsourcefilename('vcb') + ' -C ' +  self.gett2sfilename('snt',longform=True) + ' ' + self.GIZA_OPTIONS + ' -o ' + self.gett2sfilename(),'GIZA++ Word Alignment (reversed)',  self.gett2sfilename('A3.final') ): return False
        return True        

    def build_srilm_targetmodel(self):
        if not self.runcmd(self.PATH_SRILM + ' ' + self.SRILM_OPTIONS + ' -text ' + self.gettargetfilename('txt') + ' -lm ' + self.gettargetfilename('lm'),'SRILM Target-language Model', self.gettargetfilename('lm')): return False
        
    def build_srilm_sourcemodel(self):
        if not self.runcmd(self.PATH_SRILM + ' ' + self.SRILM_OPTIONS + ' -text ' + self.getsourcefilename('txt') + ' -lm ' + self.getsourcefilename('lm'),'SRILM Source-language Model', self.getsourcefilename('lm')): return False        

    def tokenize_sourcecorpus(self):
        if not os.path.exists(self.getsourcefilename('notok')):
            os.rename( os.path.exists(self.getsourcefilename('txt')), os.path.exists(self.getsourcefilename('notok') ) )
        if not self.runcmd(self.PATH_UCTO + ' ' + self.UCTO_OPTIONS +  ' -L' + self.SOURCELANG +  ' ' + self.getsourcefilename('notok') + ' ' + self.getsourcefilename('tok'),'Tokenisation Source Corpus', self.getsourcefilename('tok')): return False    
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
        if not self.runcmd(self.PATH_UCTO + ' ' + self.UCTO_OPTIONS +  ' -L' + self.SOURCELANG +  ' ' + self.gettargetfilename('notok') + ' ' + self.gettargetfilename('tok'),'Tokenisation Target Corpus', self.gettargetfilename('tok')): return False    
        if os.path.exists(self.gettargetfilename('tok')):
            try:
                os.unlink(self.gettargetfilename('txt'))
            except:
                pass
            os.symlink( self.gettargetfilename('tok'), self.gettargetfilename('txt') )
        return True
        

    def build_moses_symal(self):
        if not self.runcmd(self.PATH_MOSES_GIZA2BAL + ' -d ' + self.gett2sfilename('A3.final') + ' -i ' + self.gets2tfilename('A3.final') + ' > ' + self.gets2tfilename('bal'),'Data preparation for Symmetric Aligner', self.gets2tfilename('bal')): return False
        if not self.runcmd(self.PATH_MOSES_SYMAL + ' ' + self.SYMAL_OPTIONS + ' < ' + self.gets2tfilename('bal') + ' > '  + self.gets2tfilename('symal'), 'Moses Symmetric Alignment',self.gets2tfilename('symal')): return False 
        return True        
    
    def build_moses_wordtranstable(self):
        if not self.runcmd(self.PATH_MOSES_WORDTRANSTABLE + ' ' + self.getsourcefilename('txt') + ' ' + self.gettargetfilename('txt') + ' ' + self.gets2tfilename('bal')  + ' ' + self.gets2tfilename(),'Build of Lexical Translation Table', self.gets2tfilename('s2t'), self.gets2tfilename('t2s')): return False 
        return True
    
    def build_moses_phraseextract(self):        
        #TODO: Support for EPPEX and Hierchical rule extraction and reordering models     
        if not self.runcmd(self.PATH_MOSES_PHRASEEXTRACT + ' ' + self.getsourcefilename('txt') + ' ' + self.gettargetfilename('txt') + ' ' + self.gets2tfilename('phraseextract') + ' ' + str(self.PHRASEEXTRACT_MAX_PHRASE_LENGTH) + ' orientation ' + self.PHRASEEXTRACT_REORDERING_FLAGS ,'Phrase Extraction', self.gets2tfilename('phraseextract')): return False 
        return True        
    
    def build_moses_phrasescore(self):
        #TODO: IMplement memscore alternative?
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('s2t') + ' > ' +  self.gets2tfilename('s2t.sorted'),'Sorting Lexical Translation Table (source->target)',  self.gets2tfilename('s2t.sorted') ): return False
        if not self.runcmd('LC_ALL=C sort ' + self.gett2sfilename('t2s') + ' > ' +  self.gett2sfilename('t2s.sorted'),'Sorting Lexical Translation Table (target->source)',  self.gets2tfilename('t2s.sorted') ): return False
        
        self.PHRASESCORE_OPTIONS= '' #--Hierarchical --WordAlignment (--Inverse)
        
        if not self.runcmd(self.PATH_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gets2tfilename('phraseextract') + ' ' +  self.gets2tfilename('s2t.sorted') + ' ' + self.gets2tfilename('.halfphrasescore') + ' ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (source->target)', self.gets2tfilename('.halfphrasescore') ): return False        
        if not self.runcmd(self.PATH_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gett2sfilename('phraseextract') + ' '  + self.gets2tfilename('t2s.sorted') + ' ' + self.gett2sfilename('.halfphrasescore') + ' --Inverse ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (target->source)', self.gett2sfilename('.halfphrasescore') ): return False        
        if not self.runcmd('LC_ALL=C sort ' + self.gett2sfilename('halfphrasescore') + ' > ' +  self.gett2sfilename('halfphrasescore.sorted'),'Sorting Inverse Table',  self.gets2tfilename('halfphrasescore.sorted') ): return False        
        if not self.runcmd(self.PATH_MOSES_PHRASEEXTRACT_CONSOLIDATE + ' ' + self.gets2tfilename('.halfphrasescore') + ' ' + self.gett2sfilename('.halfphrasescore.sorted') + ' ' + self.gets2tfilename('.phrasetable') + ' ' + self.PHRASESCORE_OPTIONS, 'Consolidating two phrase table halves', self.gets2tfilename('.halfphrasescore') ): return False
        return True
        
    def build_moses_reorderingmodel(self):
        #TODO, skipped for now
        return False
    
    def build_moses_generationmodel(self):
        #TODO, skipped for now
        return False
        
    
def usage():
    print >>sys.stderr,"mtwrapper.py -- MT wrapper - Outputs a MT wrapper script (python)"
    print >>sys.stderr,"Mandatory Input:"
    print >>sys.stderr,"\t-n <name>         Name of the corpus [MANDATORY!]"
    print >>sys.stderr,"\t-s <file>         Corpus in source language"
    print >>sys.stderr,"\t-t <file>         Corpus in target language"
    print >>sys.stderr,"\t-S <code>         Source language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-T <code>         Target language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-w <dir>          Work directory (by default: current dir)"
    
    
if __name__ == "__main__":        
    
    workdir = os.getcwd()
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:t:S:T:n:w:")
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    sourcecorpusfile = targetcorpusfile = sourcelang = targetlang = corpusname = workdir = ""


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
        elif o == '-w':
            workdir = a
        else:
            usage()
            sys.exit(0)
    
    if not corpusname:
        usage()
        sys.exit(2)
        
    settingsfile = workdir + '/mt-' + corpusname + '.py'
    f = codecs.open(settingsfile,'w','utf-8')
    f.write("""#! /usr/bin/env python
# -*- coding: utf8 -*-#     
    
from mtwrapper import MTWrapper

mtwrapper = MTWrapper()
    
#Primary configuration
    
mtwrapper.CORPUSNAME = "%s"
mtwrapper.WORKDIR = "%s"

mtwrapper.TRAINSOURCECORPUS = "%s"
mtwrapper.TRAINTARGETCORPUS = "%s"

mtwrapper.SOURCELANG = "%s" #language code
mtwrapper.TARGETLANG = "%s" #language code

mtwrapper.TOKENIZE_SOURCECORPUS = False #Set to true if the corpus is not tokenised yet
mtwrapper.TOKENIZE_TARGETCORPUS = False #Set to true if the corpus is not tokenised yet
mtwrapper.BUILD_SRILM_SOURCEMODEL = False
mtwrapper.BUILD_SRILM_TARGETMODEL = True
mtwrapper.BUILD_GIZA_WORDALIGNMENT = True
mtwrapper.BUILD_GIZA_WORDALIGNMENT_REV = False #Reverse word-alignment
mtwrapper.BUILD_GIZA_WORDALIGNMENT_COOC = False #Produce co-occurrence files?
mtwrapper.BUILD_MOSES_SYMAL = True #Symmetrize word alignments
mtwrapper.BUILD_MOSES_PHRASETRANSTABLE = True


#Paths
#mtwrapper.PATH_UCTO = ""      
#mtwrapper.PATH_TIMBL = ""
#mtwrapper.PATH_MKCLS = ""
#mtwrapper.PATH_GIZA = ""
#mtwrapper.PATH_GIZA_PLAIN2SNT = ""                
#mtwrapper.PATH_MOSES = ""
#mtwrapper.PATH_SRILM = "" #path to ngram-count from SRILM

#Options for building word alignments
#mtwrapper.MKCLS_OPTIONS = '-m2 -C50' 
mtwrapper.GIZA_OPTIONS = '-p0 0.98 -m1 5 -m2 0 -m3 3 -m4 3 -nsmooth 4 -model4smoothfactor 0.4' #Using IBM Model 1,3,4
#mtwrapper.GIZA_OPTIONS = '-p0 0.98 -m1 5 -m2 0 -m3 0 -m4 0 -nsmooth 4 -hmmiterations 5 -hmmdumpfrequency -5' #Using Hidden-Markov Models

#Options for building language models with SRILM
#mtwrapper.SRILM_OPTIONS = '-order 3 -interpolate -kndiscount'

mtwrapper.start()
    """ % (corpusname, workdir, sourcecorpusfile, targetcorpusfile, sourcelang, targetlang))
    f.close()
    os.chown(settingsfile, 0754)
    os.system('vim ' + settingsfile)    



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

