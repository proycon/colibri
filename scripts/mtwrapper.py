#! /usr/bin/env python
# -*- coding: utf8 -*-

import sys
import os
import subprocess
import getopt
import codecs


class MTWrapper(object):
    defaults = [
            ('WORKDIR','','Full path to the working directory that holds all data for the system'),
            ('CORPUSNAME', '','The name of the corpus (without language codes)'),
            ('TRAINSOURCECORPUS', '','The file containing to the source-language part of the parallel corpus used for training, one sentence per line'),
            ('TRAINTARGETCORPUS', '','The file containing to the taret-language part of the parallel corpus used for training, one sentence per line'),
            ('SOURCELANG', '','A language code identifying the source language'),
            ('TARGETLANG', '','A language code identifying the target language'),
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
            ('PATH_MOSES', '','Base directory where Moses is installed'),
            ('PATH_SRILM', '','Base directory where SRILM is installed'),
            ('PATH_GIZA', '','Base directory where GIZA++ is installed'),
            ('PATH_COLIBRI', '','Base directory where COLIBRI is installed'),
            ('EXEC_UCTO', 'ucto','Path to ucto binary'),
            ('EXEC_SRILM', 'ngram-count','Path to ngram-count (SRILM)'),
            ('EXEC_TIMBL', 'timbl','Path to timbl binary'),
            ('EXEC_GIZA_MKCLS', 'mkcls','Path to mkcls (part of GIZA++)'),
            ('EXEC_GIZA', 'GIZA++','Path to GIZA++ binary'),
            ('EXEC_GIZA_PLAIN2SNT', 'plain2snt.out','Path to plain2snt.out (part of GIZA++)'),
            ('EXEC_GIZA_SNT2COOC', 'snt2cooc.out','Path to snt2cooc.out (part of GIZA++)'),
            ('EXEC_MOSES', 'moses','Path to Moses binary'),
            ('EXEC_MOSES_GIZA2BAL', 'scripts/training/symal/giza2bal.pl', ''),
            ('EXEC_MOSES_SYMAL', 'scripts/training/symal/symal', ''),
            ('EXEC_MOSES_WORDTRANSTABLE','scripts/moses-lexicaltranslationtable.py',''),
            ('EXEC_MOSES_PHRASEEXTRACT','scripts/training/phrase-extract/extract',''),
            ('EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE','scripts/training/phrase-extract/consolidate',''),
            ('EXEC_MOSES_PHRASEEXTRACT_SCORE','scripts/training/phrase-extract/score',''),
            ('MKCLS_OPTIONS','-m2 -c50',''),
            ('GIZA_OPTIONS','-p0 0.98 -m1 5 -m2 0 -m3 3 -m4 3 -nsmooth 4 -model4smoothfactor 0.4',''),
            ('SRILM_OPTIONS','-order 3 -interpolate -kndiscount',''),
            ('UCTO_OPTIONS','-m -n',''),
            ('SYMAL_OPTIONS','-alignment=grow -diagonal=yes -final=yes -both=no',''), #-hmmiterations 5 -hmmdumpfrequency -5'
            ('PHRASEEXTRACT_MAX_PHRASE_LENGTH',7,''),
            ('PHRASEEXTRACT_REORDERING_FLAGS','',''), #" --model wbe-mslr --model phrase-mslr --model hier-mslr" #Maximum lexical reordering 
            ('PHRASESCORE_OPTIONS', '',''), #--Hierarchical --WordAlignment (--Inverse)        
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
        
        
        for key in kwargs:
            print >>sys.stderr, "Unknown configuration directive: " + key
            sys.exit(2)
            
        
        

    def findpath(self, name, basepath = ''):                        
        for path in os.environ['PATH'].split(':'):
            if os.path.exists(path + '/' + name) and not os.path.isdir(path + '/' + name):
                print >>sys.stderr, green("Found " + name + " in " + path)
                return path + '/' + name
        if basepath: 
            return basepath + '/' + name    
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
            
        if (self.TOKENIZE_SOURCECORPUS or self.TOKENIZE_TARGETCORPUS) and (not self.EXEC_UCTO or not os.path.isfile(self.EXEC_UCTO)):
            print >>sys.stderr,red("Dependency error: ucto not found (EXEC_UCTO=" + self.EXEC_UCTO + ")")
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

            if not self.EXEC_MOSES_GIZA2BAL or not os.path.isfile(self.EXEC_MOSES_GIZA2BAL):
                sane = False
                print >>sys.stderr,red("Dependency error: giza2bal.pl (provided by Moses) not found (EXEC_MOSES_GIZA2BAL=" + self.EXEC_MOSES_GIZA2BAL + ")")
               
            
            if not self.EXEC_MOSES_SYMAL or not os.path.isfile(self.EXEC_MOSES_SYMAL):
                sane = False
                print >>sys.stderr,red("Dependency error: symal (provided by Moses) not found (EXEC_MOSES_SYMAL=" + self.EXEC_MOSES_SYMAL + ")")
               
            if not self.EXEC_MOSES or not os.path.isfile(self.EXEC_MOSES):
                sane = False
                print >>sys.stderr,red("Dependency error: Moses not found (EXEC_MOSES=" + self.EXEC_MOSES + ")")
                
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
        if not self.runcmd(self.EXEC_SRILM + ' ' + self.SRILM_OPTIONS + ' -text ' + self.gettargetfilename('txt') + ' -lm ' + self.gettargetfilename('lm'),'SRILM Target-language Model', self.gettargetfilename('lm')): return False
        
    def build_srilm_sourcemodel(self):
        if not self.runcmd(self.EXEC_SRILM + ' ' + self.SRILM_OPTIONS + ' -text ' + self.getsourcefilename('txt') + ' -lm ' + self.getsourcefilename('lm'),'SRILM Source-language Model', self.getsourcefilename('lm')): return False        

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
        if not self.runcmd(self.EXEC_MOSES_GIZA2BAL + ' -d ' + self.gett2sfilename('A3.final') + ' -i ' + self.gets2tfilename('A3.final') + ' > ' + self.gets2tfilename('bal'),'Data preparation for Symmetric Aligner', self.gets2tfilename('bal')): return False
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
        #TODO: IMplement memscore alternative?
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('s2t') + ' > ' +  self.gets2tfilename('s2t.sorted'),'Sorting Lexical Translation Table (source->target)',  self.gets2tfilename('s2t.sorted') ): return False
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('t2s') + ' > ' +  self.gett2sfilename('t2s.sorted'),'Sorting Lexical Translation Table (target->source)',  self.gets2tfilename('t2s.sorted') ): return False
        
        
        
        if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gets2tfilename('phraseextract') + ' ' +  self.gets2tfilename('s2t.sorted') + ' ' + self.gets2tfilename('.sourcehalf') + ' ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (source->target)', self.gets2tfilename('.sourcehalf') ): return False        
        if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_SCORE + ' ' + self.gets2tfilename('phraseextract') + ' '  + self.gets2tfilename('t2s.sorted') + ' ' + self.gets2tfilename('.targethalf') + ' --Inverse ' + self.PHRASESCORE_OPTIONS, 'Scoring phrases (target->source)', self.gett2sfilename('.targethalf') ): return False        
        if not self.runcmd('LC_ALL=C sort ' + self.gets2tfilename('targethalf') + ' > ' +  self.gett2sfilename('targethalf.sorted'),'Sorting Inverse Table',  self.gets2tfilename('targethalf.sorted') ): return False
                
        if not self.runcmd(self.EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE + ' ' + self.gets2tfilename('.sourcehalf') + ' ' + self.gets2tfilename('.targethalf.sorted') + ' ' + self.gets2tfilename('.phrasetable') + ' ' + self.PHRASESCORE_OPTIONS, 'Consolidating two phrase table halves', self.gets2tfilename('.phrasetable') ): return False
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
    f.write("#! /usr/bin/env python\n# -*- coding: utf8 -*-#\n\nfrom mtwrapper import MTWrapper\n")
    f.write("mtwrapper = MTWrapper(\n")
    for key, default, help in MTWrapper.defaults:            
        if key == 'CORPUSNAME': 
            default = corpusname
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

