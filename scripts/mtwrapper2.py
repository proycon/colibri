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
import datetime
import time
import shutil
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from pynlpl.evaluation import filesampler, ExperimentPool, AbstractExperiment
from pynlpl.net import GenericWrapperServer


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




def serveroutputproc(line):
    #print >>sys.stderr, "SERVEROUTPUTPROC: ", line
    prefix = 'BEST TRANSLATION: '
    if line[:len(prefix)] == prefix:
        outputwords = []
        words = line[18:].split(' ')
        for word in words:
            if word:
                if word[0] == '[' and word[-1] == ']':
                    ones = True
                    for c in word[1:-1]:
                        if c != '1':
                            ones = False
                            break
                    if ones:
                        return " ".join(outputwords) + "<br />"
                else:
                    factors = word.split('|')
                    if len(factors) == 1:
                        outputwords.append(word)
                    elif factors[1] == 'UNK':
                        outputwords.append("<em>" + factors[0] + "</em> ")
                    else:
                        outputwords.append(factors[0])
        return "<pre>ERROR! UNABLE TO FIND END OF TRANSLATION OUTPUT! DEBUG OUTPUT: ["+line+"]</pre>"
    else:
        return ""

class BatchExperiment(AbstractExperiment):

    def start(self):
        mtwrapper, batch, conf, train, test = self.inputdata

        mtwrapper.log("Branching for batch " + batch,white,True)
        mtwrapper.branch(batch, conf, False, True, False)


        batchdir = mtwrapper.WORKDIR + '/' + mtwrapper.CORPUSNAME + '-' + mtwrapper.SOURCELANG + '-' + mtwrapper.TARGETLANG + '-' + batch
        if train and test:
            mtwrapper.log("Starting batch " + batch + " " + mtwrapper.timestamp(),white,True)
            mtwrapper.log("Logs will be in " + batchdir + "/train.log and " + batchdir + "/test.log")
            self.startcommand(batchdir + '/mt-' +  mtwrapper.CORPUSNAME + '-' + mtwrapper.SOURCELANG + '-' + mtwrapper.TARGETLANG + '-' + batch + '.py',batchdir,sys.stdout, sys.stderr,'start')
        elif train:
            mtwrapper.log("Starting training batch " + batch + " " + mtwrapper.timestamp(),white,True)
            mtwrapper.log(" Log will be in " + batchdir + '/train.log')
            self.startcommand(batchdir + '/mt-' +  mtwrapper.CORPUSNAME + '-' + mtwrapper.SOURCELANG + '-' + mtwrapper.TARGETLANG + '-' + batch + '.py',batchdir,sys.stdout, sys.stderr,'train')
        elif test:
            mtwrapper.log("Starting testing batch " + batch + " " + mtwrapper.timestamp(),white,True)
            mtwrapper.log(" Log will be in " + batchdir + '/test.log')
            self.startcommand(batchdir + '/mt-' +  mtwrapper.CORPUSNAME + '-' + mtwrapper.SOURCELANG + '-' + mtwrapper.TARGETLANG + '-' + batch + '.py',batchdir,sys.stdout, sys.stderr,'test')
        else:
            raise Exception("Don't know what to do")


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
            ('TRAINLIMIT','','If set, only use the first x sentences for training'),
            ('TOKENIZE_SOURCECORPUS', False,''),
            ('TOKENIZE_TARGETCORPUS', False,''),
            ('BUILD_SRILM_SOURCEMODEL',False,'Build a source-language model'),
            ('BUILD_SRILM_TARGETMODEL',False,'Build a target-language model'),
            ('BUILD_GIZA_WORDALIGNMENT',False,'Build GIZA++ Word Alignments'),
            ('BUILD_GIZA_WORDALIGNMENT_REV',False,'Build GIZA++ Reverse Word Alignment (target to source)'),
            ('BUILD_GIZA_WORDALIGNMENT_COOC',False,'Output extra co-occurence data'),
            #('BUILD_MOSES_SYMAL', False,'Symmetrise word alignments'),
            #('BUILD_MOSES_WORDTRANSTABLE', False,'Build lexical translation table'),
            #('BUILD_MOSES_PHRASEEXTRACT', False,'Extract phrases'),
            ('BUILD_MOSES_PHRASETRANSTABLE', False,'Build phrase translation table'),
            ('BUILD_MOSES_CLASSIFIERS', False, 'Use context-aware classifiers. Set COLIBRI_CLASSIFIER_OPTIONS to select a classifier method.'),
            ('BUILD_MOSES_MEMSCORE', False,'Use memscore to score phrases rather than the default phrase-extract scorer'),
            ('BUILD_MOSES', False,'Build moses configuration, necessary for decoding using moses'),
            ('BUILD_MOSES_MERT', False,'Do Minimum Error Rate Training for Moses (on development set)'),
            ('BUILD_PBMBMT', False, 'Build model for Phrase-Based Memory-based Machine Translation'),
            ('BUILD_PBMBMT_PARAMSEARCH', False, 'Do parameter optimisation for PBMBMT using wrapped progressive sampling'),
            ('BUILD_COLIBRI_ALIGNMENT', False,'Create an alignment using colibri'),
            ('BUILD_COLIBRI_GIZA', False,'Base aligner on word-alignments using giza (do not manually specify -W -s -t in COLIBRI_ALIGNER_OPTIONS)'),
            ('BUILD_COLIBRI_TRANSTABLE', False,'Build a translation table using colibri'),
            ('BUILD_COLIBRI_MOSESPHRASETABLE', False,'Build a Moses-style phrasetable using colibri'),
            ('BUILD_COLIBRI_CLASSIFIERS', False,'Build and use classifiers. Also set MOSES_CLASSIFIER_OPTIONS to select a classifier method'),
            ('BUILD_COLIBRI_SKIPGRAMS', False,'Include support for skipgrams (automatically adds -s -B -E to patternfinder options, creates proper graph models if enabled, extracts skipgrams in alignment). If alignment is enabled, skipgrams will be extracted from the aligned models, do not use if skipgrams are already to be extracted during EM or Jaccard alignment.'),
            ('BUILD_COLIBRI_GRAPHMODEL', False,'Build a graph model using colibri'),
            ('BUILD_COLIBRI', False,'Build for colibri decoder'),
            ('BUILD_PHRASAL', False,'Build phrasal configuration, necessary for decoding using phrasal'),
            ('BUILD_PHRASAL_MERT', False,'Do Minumum Error Rate Training for Phrasal (on development set)'),
            ('BUILD_PHRASAL_WORDALIGN', False,'Align words using berkeley aligner (supplied with phrasal)'),
            ('BUILD_PHRASAL_PHRASEEXTRACT', False,'Extract phrases'),
            ('PATH_MOSES', '','Base directory where Moses is installed'),
            ('PATH_SRILM', '','Base directory where SRILM is installed'),
            ('PATH_GIZA', '','Base directory where GIZA++ is installed'),
            ('PATH_COLIBRI', '','Base directory where COLIBRI is installed'),
            ('PATH_MATREX','','Base directory for Matrex evaluation scripts'),
            ('PATH_PBMBMT','','Base directory to PBMBMT'),
            ('PATH_CORENLP','','Base directory to Stanford Core NLP'),
            ('PATH_PHRASAL','','Base directory to Stanford Phrasal'),
            ('PATH_MOSES_EXTERNALBIN','/vol/customopt/machine-translation/bin',''),
            ('JAR_BERKELEYALIGNER','berkeleyaligner.jar','Path to berkeleyaligner.jar (v2.1+) (can be automatically determined if PATH_PHRASAL is set)'),
            ('JAR_FASTUTIL','fastutil.jar','Path to fastutil.jar (can be automatically determined if PATH_PHRASAL is set)'),
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
            ('EXEC_MOSES_GIZA2BAL', 'scripts/training/giza2bal.pl', ''),
            ('EXEC_MOSES_SYMAL', 'symal', ''),
            ('EXEC_MOSES_WORDTRANSTABLE','scripts/moses-lexicaltranslationtable.py',''),
            ('EXEC_MOSES_PHRASEEXTRACT','scripts/training/phrase-extract/extract',''),
            ('EXEC_MOSES_PHRASEEXTRACT_CONSOLIDATE','scripts/training/phrase-extract/consolidate',''),
            ('EXEC_MOSES_PHRASEEXTRACT_SCORE','scripts/training/phrase-extract/score',''),
            ('EXEC_MOSES_MEMSCORE','scripts/training/memscore/memscore',''),
            ('EXEC_MOSES_MERT','scripts/training/mert-moses.pl',''),
            ('EXEC_MOSES_TRAINMODEL','scripts/training/train-model.perl',''),
            ('EXEC_PHRASAL','scripts/decode',''),
            ('EXEC_PHRASAL_MERT','scripts/phrasal-mert.pl',''),
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
            ('EXEC_COLIBRI_CONTEXTMOSES','contextmoses',''),
            ('EXEC_COLIBRI_GRAPHER','grapher',''),
            ('EXEC_COLIBRI_ALIGNER','aligner',''),
            ('EXEC_COLIBRI_DECODER','decoder',''),
            ('EXEC_COLIBRI_TRAINCLASSIFIERS','trainclassifiers',''),
            ('MKCLS_OPTIONS','-m2 -c50',''),
            ('GIZA_OPTIONS','-p0 0.98 -m1 5 -m2 0 -m3 3 -m4 3 -nsmooth 4 -model4smoothfactor 0.4',''),
            ('SRILM_ORDER',3,'N-gram size for language model'),
            ('SRILM_OPTIONS','-interpolate -kndiscount -unk','Further SRILM options (do not use -order here, use SRILM_ORDER instead)'),
            ('UCTO_OPTIONS','-m -n',''),
            ('SYMAL_OPTIONS','-alignment=grow -diagonal=yes -final=yes -both=no',''), #-hmmiterations 5 -hmmdumpfrequency -5'
            ('MOSES_MERT_OPTIONS','','See http://www.statmt.org/moses/?n=FactoredTraining.Tuning'),
            ('MOSES_MERT_RUNS',1,'Number of MERT runs to perform (results will be averaged)'),
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
            ('COLIBRI_TIMBL_OPTIONS', '-k 1 -a4','Timbl options. See Timbl -h'),
            ('COLIBRI_GRAPHER_OPTIONS','-P -X -r','Options for the Graphmodel, if empty, no graph model will be constructed for the aligner, see graphmodel -h'),
            ('COLIBRI_PATTERNFINDER_OPTIONS','-t 10', 'Options for the pattern finder, see patternfinder -h'),
            ('COLIBRI_ALIGNER_OPTIONS','-J -p 0.1','Options for the colibri aligner, see aligner -h'),
            ('COLIBRI_DECODER_OPTIONS','','Options for the colibri decoder, see decoder -h'),
            ('COLIBRI_LEFTCONTEXTSIZE',1,'For use with BUILD_COLIBRI_CLASSIFIERS=True'),
            ('COLIBRI_RIGHTCONTEXTSIZE',1,'For use with BUILD_COLIBRI_CLASSIFIERS=True'),
            ('COLIBRI_CLASSIFIER_OPTIONS','-N','For use with BUILD_COLIBRI_CLASSIFIERS=True. Make sure to select at least a classifier here (-N, -X,-M). See trainclassifiers -h for all options'),
            ('COLIBRI_GLOBALKEYWORDS',False,'Enable global context keywords (also affects BUILD_MOSES_CLASSIFIERS)'),
            ('COLIBRI_GLOBALKEYWORDS_OPTIONS','1,3,20,0.0000000009','Four comma-separated values for computation of global context keywords: include_threshold,absolute_threshold,filter_threshold,probability_threshold'),
            ('MOSES_LEFTCONTEXTSIZE',1,'For use with BUILD_MOSES_CLASSIFIERS=True'),
            ('MOSES_RIGHTCONTEXTSIZE',1,'For use with BUILD_MOSES_CLASSIFIERS=True'),
            ('MOSES_CLASSIFIER_OPTIONS','-N','For use with BUILD_MOSES_CLASSIFIERS=True. Make sure to select at least a classifier (-N, -X,-M) here. See contextmoses -h for all options.'),
            ('PHRASAL_MAXMEM', '4g', 'Memory allocated for word alignment, phrase extraction and decoding using phrasal (java)'),
            ('PHRASAL_WITHGAPS', True, 'Consider gaps if using Phrasal?'),
            ('PHRASAL_MAXSOURCEPHRASESPAN', 15, 'Maximum span for a source-side phrase with gaps (phrasal)'),
            ('PHRASAL_MAXTARGETPHRASESPAN', 7, 'Maximum span for a target-side phrase with gap (phrasal)'),
            ('PHRASAL_PHRASEEXTRACT_OPTIONS','-symmetrization grow-diag-final-and','Options for Phrase Extraction using Phrasal')
    ]

    def initlog(self, logfile):
        try:
            self.logfile = open(self.WORKDIR + '/' + logfile + '.log','w')
        except:
            print >>sys.stderr,"Unable to write to logfile: " + self.WORKDIR + '/' + logfile + '.log'


    def log(self, msg, color=None, dobold = False):
        if color:
            msg = color(msg)
        if dobold:
            msg = bold(msg)
        if self.logfile:
            self.logfile.write(msg+"\n")
        print >>sys.stderr, msg


    def parsekwargs(self, key, default, **kwargs):
        if key in kwargs:
            return kwargs[key]
            del kwargs[key]
        else:
            return default


    def __init__(self, *args, **kwargs):
        self.batches = []
        self.logfile = None

        self.confdata = {}
        self.sources = args
        if args:
            for arg in self.sources:
                exec 'from ' + arg + ' import confdata'
                for key, value in confdata.items():
                    self.confdata[key] = value

        for key, value in kwargs.items():
            self.confdata[key] = value
        self.setargs(**self.confdata)



    def setargs(self, **kwargs):
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
        self.EXEC_MOSES_MEMSCORE  = self.findpath(self.EXEC_MOSES_MEMSCORE , self.PATH_MOSES)

        self.EXEC_PHRASAL = self.findpath(self.EXEC_PHRASAL,self.PATH_PHRASAL)
        self.JAR_BERKELEYALIGNER = self.findpath(self.JAR_BERKELEYALIGNER, self.PATH_PHRASAL + '/lib/')
        self.JAR_FASTUTIL = self.findpath(self.JAR_FASTUTIL, self.PATH_PHRASAL + '/lib/')
        self.EXEC_PHRASAL_MERT = self.findpath(self.EXEC_PHRASAL_MERT,self.PATH_PHRASAL)

        self.EXEC_MATREX_WER = self.findpath(self.EXEC_MATREX_WER, self.PATH_MATREX)
        self.EXEC_MATREX_PER = self.findpath(self.EXEC_MATREX_PER, self.PATH_MATREX)
        self.EXEC_MATREX_BLEU = self.findpath(self.EXEC_MATREX_BLEU, self.PATH_MATREX)
        self.EXEC_MATREX_METEOR = self.findpath(self.EXEC_MATREX_METEOR, self.PATH_MATREX)
        self.EXEC_MATREX_MTEVAL = self.findpath(self.EXEC_MATREX_MTEVAL, self.PATH_MATREX)
        self.EXEC_MATREX_TER = self.findpath(self.EXEC_MATREX_TER, self.PATH_MATREX)

        self.EXEC_COLIBRI_CLASSENCODE = self.findpath(self.EXEC_COLIBRI_CLASSENCODE , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_PATTERNFINDER = self.findpath(self.EXEC_COLIBRI_PATTERNFINDER , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_GRAPHER = self.findpath(self.EXEC_COLIBRI_GRAPHER , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_ALIGNER = self.findpath(self.EXEC_COLIBRI_ALIGNER , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_DECODER = self.findpath(self.EXEC_COLIBRI_DECODER , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_TRAINCLASSIFIERS = self.findpath(self.EXEC_COLIBRI_TRAINCLASSIFIERS , self.PATH_COLIBRI)
        self.EXEC_COLIBRI_CONTEXTMOSES = self.findpath(self.EXEC_COLIBRI_CONTEXTMOSES , self.PATH_COLIBRI)

        self.EXEC_PBMBMT_DECODER = self.findpath(self.EXEC_PBMBMT_DECODER, self.PATH_PBMBMT)
        self.EXEC_PBMBMT_INSTANCEGENERATOR = self.findpath(self.EXEC_PBMBMT_INSTANCEGENERATOR, self.PATH_PBMBMT)

        if self.PATH_MOSES:
            self.PATH_MOSES_MERT = self.PATH_MOSES + '/mert'
        else:
            self.PATH_MOSES_MERT = ''

        for key in kwargs:
            self.log("Unknown configuration directive: " + key,red)
            sys.exit(2)


    def addbatch(self, batchname, **kwargs):
        for key in kwargs:
            try:
                getattr(self,key)
            except AttributeError:
                self.log("Unknown configuration directive in batch " + batchname + ":" + key,red)
        self.batches.append( (batchname, kwargs) )

    def findpath(self, name, basepath = ''):
        if not name.strip():
            return ""
        if os.path.exists(name):
            return name
        for path in os.environ['PATH'].split(':'):
            if os.path.exists(path + '/' + name) and not os.path.isdir(path + '/' + name):
                self.log("Found " + name + " in " + path,green)
                return path + '/' + name
            elif os.path.exists(path + '/' + os.path.basename(name) ) and not os.path.isdir(path + '/' + os.path.basename(name) ):
                self.log("Found " + os.path.basename(name) + " in " + path,green)
                return path + '/' + os.path.basename(name)
        if basepath and os.path.exists(basepath + '/' + name):
            self.log("Found " + name + " in " + basepath,green)
            return basepath + '/' + name
        print >>sys.stderr, yellow("Warning: Did not find " + name + ' ('+ basepath+')' )
        return ""

    def check_common(self):
        sane = True

        if not os.path.isdir(self.WORKDIR):
            self.log("Work directory does not exist, creating " + self.WORKDIR,yellow)
            try:
                os.mkdir(self.WORKDIR)
            except:
                self.log("Configuration error: Unable to create work directory " + self.WORKDIR,red)
                sane = False

        if self.WORKDIR[-1] != '/':
            self.WORKDIR += '/'

        sane = True
        if not self.CORPUSNAME:
            self.log("Configuration error: CORPUSNAME not specified!",red)
            sane = False
        if not self.SOURCELANG:
            self.log("Configuration error: SOURCELANG not specified!",red)
            sane = False
        if not self.TARGETLANG:
            self.log("Configuration error: TARGETLANG not specified!",red)
            sane = False
        return sane


    def check_run(self):
        if self.BUILD_MOSES:
            if not self.EXEC_MOSES or not os.path.isfile(self.EXEC_MOSES):
                self.log("Error: Moses not found! (" + self.EXEC_MOSES+")",red,True)
                return False
            elif not os.path.exists(self.WORKDIR + '/model/moses.ini'):
                self.log("Error: No Moses configuration found. Did you forget to train the system first?",red,True)
                return False
            elif not os.path.exists(self.gets2tfilename('phrasetable')) and not not os.path.exists('model/' + self.gets2tfilename('phrase-table')):
                self.log("Error: No Moses phrasetable found ("+ self.gets2tfilename('phrasetable')+") . Did you forget to train the system first?",red,True)
                return False
        elif self.BUILD_PBMBMT:
            #TODO: implement
            self.log("Error: PBMBMT not implemented yet",red,True)
            return False
        elif self.BUILD_PHRASAL:
            if not self.EXEC_JAVA or not os.path.isfile(self.EXEC_JAVA):
                self.log("Error: Java not found! Required for Phrasal",red,True)
                return False
            if not self.PATH_PHRASAL or not os.path.isdir(self.PATH_PHRASAL):
                self.log("Error: Phrasal not found! (" + self.PATH_PHRASAL+")",red,True)
                return False
            elif not os.path.exists(self.WORKDIR + '/phrasal.conf'):
                self.log("Error: No Phrasal configuration found. Did you forget to train the system first?",red,True)
                return False
            elif not os.path.exists(self.gets2tfilename('phrasetable')):
                self.log("Error: No phrasetable found ("+ self.gets2tfilename('phrasetable')+") . Did you forget to train the system first?",red,True)
                return False
        elif self.BUILD_COLIBRI:
            if not self.EXEC_COLIBRI_DECODER or not os.path.isfile(self.EXEC_COLIBRI_DECODER):
                self.log("Error: Colibri decoder not found! (" + self.EXEC_COLIBRI_DECODER+")",red,True)
                return False
            elif (not os.path.exists(self.gets2tfilename('translationtable.colibri'))) and (not os.path.exists(self.gets2tfilename('alignmodel.colibri'))) and (not os.path.exists(self.gets2tfilename('alignmodelS.colibri')))  and (not os.path.exists(self.gets2tfilename('phrasetable'))):
                self.log("Error: No translation table found. Did you forget to train the system first?",red,True)
                return False
        else:
            self.log("Error: System is not runnable, no MT decoder enabled",red,True)
            return False
        return True

    def check_test(self):
        if not (self.EXEC_MATREX_WER or self.EXEC_MATREX_PER or self.EXEC_MATREX_BLEU or self.EXEC_MATREX_MTEVAL or self.EXEC_MATREX_METEOR or self.EXEC_MATREX_TER):
            self.log("Error: No evaluation scripts found, set at least one of EXEC_MATREX_* and PATH_MATREX",red,True)
            return False
        return True

    def check_train(self):
        sane = True
        if not self.TRAINSOURCECORPUS:
            self.log("Configuration error: TRAINSOURCECORPUS not specified!",red,True)
            sane = False
        if not self.TRAINTARGETCORPUS:
            self.log("Configuration error: TRAINTARGETCORPUS not specified!",red,True)
            sane = False

        if (self.TOKENIZE_SOURCECORPUS or self.TOKENIZE_TARGETCORPUS) and (not self.EXEC_UCTO or not os.path.isfile(self.EXEC_UCTO)):
            self.log("Dependency error: ucto not found (EXEC_UCTO=" + self.EXEC_UCTO + ")",red)
            sane = False


        if self.BUILD_COLIBRI_GIZA:
            if not self.BUILD_COLIBRI_ALIGNMENT:
                self.log("Configuration update: BUILD_COLIBRI_ALIGNMENT automatically enabled because BUILD_COLIBRI_GIZA is too",yellow)
                self.BUILD_COLIBRI_ALIGNMENT = True

        if self.BUILD_COLIBRI_MOSESPHRASETABLE:
            if not self.BUILD_COLIBRI_TRANSTABLE:
                self.log("Configuration update: BUILD_COLIBRI_TRANSTABLE automatically enabled because BUILD_COLIBRI_MOSESPHRASETABLE is too",yellow)
                self.BUILD_COLIBRI_TRANSTABLE = True

        if self.BUILD_COLIBRI_TRANSTABLE:
            if not self.BUILD_COLIBRI_ALIGNMENT:
                self.log("Configuration update: BUILD_COLIBRI_ALIGNMENT automatically enabled because BUILD_COLIBRI_TRANSTABLE is too",yellow)
                self.BUILD_COLIBRI_ALIGNMENT = True


        if self.BUILD_COLIBRI_ALIGNMENT:
            if not self.EXEC_COLIBRI_PATTERNFINDER or not os.path.exists(self.EXEC_COLIBRI_PATTERNFINDER):
                sane = False
                self.log("Configuration error: EXEC_COLIBRI_PATTERNFINDER not found ! Required for BUILD_COLIBRI_ALIGNMENT !",red)
            if not self.EXEC_COLIBRI_ALIGNER or not os.path.exists(self.EXEC_COLIBRI_ALIGNER):
                sane = False
                self.log("Configuration error: EXEC_COLIBRI_ALIGNER not found ! Required for BUILD_COLIBRI_ALIGNMENT !",red)
            if not self.EXEC_COLIBRI_GRAPHER or not os.path.exists(self.EXEC_COLIBRI_GRAPHER):
                sane = False
                self.log("Configuration error: EXEC_COLIBRI_GRAPHER not found ! Required for BUILD_COLIBRI_ALIGNMENT !",red)
            if not self.EXEC_COLIBRI_CLASSENCODE or not os.path.exists(self.EXEC_COLIBRI_CLASSENCODE):
                sane = False
                self.log("Configuration error: EXEC_COLIBRI_CLASSENCODE not found ! Required for BUILD_COLIBRI_ALIGNMENT !",red)

            if self.BUILD_COLIBRI_SKIPGRAMS:
                self.log("Configuration update: BUILD_COLIBRI_GRAPHMODEL automatically enabled because BUILD_COLIBRI_SKIPGRAMS and BUILD_COLIBRI_ALIGNMENT is too",yellow)
                self.BUILD_COLIBRI_GRAPHMODEL = True

        if self.BUILD_PHRASAL_MERT:
            if not self.BUILD_PHRASAL:
                self.log("Configuration update: BUILD_PHRASAL automatically enabled because BUILD_PHRASAL_MERT is too",yellow)
                self.BUILD_PHRASAL = True
            if not self.DEVSOURCECORPUS or not os.path.exists(self.DEVSOURCECORPUS):
                sane = False
                self.log("Configuration error: DEVSOURCECORPUS not found! Required for BUILD_PHRASAL_MERT !",red)
            if not self.DEVTARGETCORPUS or not os.path.exists(self.DEVTARGETCORPUS):
                sane = False
                self.log("Configuration error: DEVTARGETCORPUS not found! Required for BUILD_PHRASAL_MERT !",red)
            if not self.PATH_PHRASAL_MERT or not os.path.isdir(self.PATH_PHRASAL_MERT):
                sane = False
                self.log("PATH_PHRASAL_MERT not found, please set PATH_MOSES !",red)

        if self.BUILD_PHRASAL:
            if self.BUILD_MOSES or self.BUILD_PBMBMT or self.BUILD_COLIBRI:
                self.log("Configuration error: Ambiguous selection of MT system: Select only one of BUILD_MOSES, BUILD_PBMBMT, BUILD_PHRASAL, BUILD_COLIBRI",red)
                sane = False
            if not self.BUILD_PHRASAL_PHRASEEXTRACT and not self.BUILD_COLIBRI_MOSESPHRASETABLE:
                self.log("Configuration update: BUILD_PHRASAL_EXTRACTPHRASES automatically enabled because BUILD_PHRASAL is too",yellow)
                self.BUILD_PHRASAL_PHRASEEXTRACT = True
            if not self.BUILD_SRILM_TARGETMODEL:
                self.log("Configuration update: BUILD_SRILM_TARGETMODEL automatically enabled because BUILD_PHRASAL is too",yellow)
                self.BUILD_SRILM_TARGETMODEL = True
            if not self.EXEC_PHRASAL or not os.path.isfile(self.EXEC_PHRASAL):
                sane = False
                self.log("Phrasal not found! Set EXEC_PHRASAL or PATH_PHRASAL !",red)



        if self.BUILD_PHRASAL_PHRASEEXTRACT:
            if not self.BUILD_PHRASAL_WORDALIGN and not (self.BUILD_GIZA_WORDALIGNMENT and self.BUILD_GIZA_WORDALIGNMENT_REV):
                self.log("Configuration update: BUILD_PHRASAL_WORDALIGN automatically enabled because BUILD_PHRASAL_PHRASEEXTRACT is too",yellow)
                self.BUILD_PHRASAL_WORDALIGN = True
            if not self.DEVSOURCECORPUS:
                sane = False
                self.log("A development corpus for the source-language is needed as a filter for phrase extraction for phrasal. Set DEVSOURCECORPUS",red)

        if self.BUILD_PHRASAL_WORDALIGN and not self.JAR_BERKELEYALIGNER:
            sane = False
            self.log("Berkeley Aligner (v2.1+) not found! Set JAR_BERKELEYALIGNER or PATH_PHRASAL !",red)


        if self.BUILD_MOSES_CLASSIFIERS:
            if not self.BUILD_MOSES:
                self.log("Configuration update: BUILD_MOSES automatically enabled because BUILD_MOSES_CLASSIFIERS is too",yellow)
                self.BUILD_MOSES = True

        if self.BUILD_MOSES_MERT:
            if not self.BUILD_MOSES:
                self.log("Configuration update: BUILD_MOSES automatically enabled because BUILD_MOSES_MERT is too",yellow)
                self.BUILD_MOSES = True
            if not self.DEVSOURCECORPUS or not os.path.exists(self.DEVSOURCECORPUS):
                sane = False
                self.log("Configuration error: DEVSOURCECORPUS not found! Required for BUILD_MOSES_MERT !",red)
            if not self.DEVTARGETCORPUS or not os.path.exists(self.DEVTARGETCORPUS):
                sane = False
                self.log("Configuration error: DEVTARGETCORPUS not found! Required for BUILD_MOSES_MERT !",red)
            if not self.PATH_MOSES_MERT or not os.path.isdir(self.PATH_MOSES_MERT):
                sane = False
                self.log("PATH_MOSES_MERT not found, please set PATH_MOSES !",red)

        if self.BUILD_MOSES:
            if self.BUILD_PBMBMT or self.BUILD_PHRASAL or self.BUILD_COLIBRI:
                self.log("Configuration error: Ambiguous selection of MT system: Select only one of BUILD_MOSES, BUILD_PBMBMT, BUILD_PHRASAL, BUILD_COLIBRI",red)
                sane = False
            if not self.BUILD_MOSES_PHRASETRANSTABLE and not self.BUILD_COLIBRI_MOSESPHRASETABLE:
                self.log("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_MOSES is too",yellow)
                self.BUILD_MOSES_PHRASETRANSTABLE = True
            if not self.BUILD_SRILM_TARGETMODEL:
                self.log("Configuration update: BUILD_SRILM_TARGETMODEL automatically enabled because BUILD_MOSES is too",yellow)
                self.BUILD_SRILM_TARGETMODEL = True
            if not self.EXEC_MOSES or not os.path.isfile(self.EXEC_MOSES):
                sane = False
                self.log("Moses not found! Set EXEC_MOSES or PATH_MOSES !",red)



        if self.BUILD_PBMBMT:
            if not self.PBMBMT_PHRASETABLE and not self.BUILD_MOSES_PHRASETABLE:
                self.log("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_PBMBMT is enabled and no pre-existing phrasetable is set (PBMBMT_PHRASETABLE)",yellow)
            if not self.PBMBMT_GIZAALIGNMENT and not self.BUILD_GIZA_WORDALIGNMENT:
                self.log("Configuration update: BUILD_GIZA_WORDALIGNMENT automatically enabled because BUILD_PBMBMT is enabled and no pre-existing word alignment file is set (PBMBMT_GIZAALIGNMENT)",yellow)
            if not self.BUILD_SRILM_TARGETMODEL:
                self.log("Configuration update: BUILD_SRILM_TARGETMODEL automatically enabled because BUILD_PBMBMT is too",yellow)
                self.BUILD_SRILM_TARGETMODEL = True
            if self.PBMBMT_PHRASETABLE:
                if not os.path.isfile(self.PBMBMT_PHRASETABLE):
                    self.log("Configuration error: PBMBMT_PHRASETABLE does not exist!",red)
                    sane = False
                else:
                    os.symlink(self.PBMBMT_PHRASETABLE, self.gets2tfilename('phrasetable'))
            if self.PBMBMT_GIZAALIGNMENT:
                if not os.path.isfile(self.PBMBMT_GIZAALIGNMENT):
                    self.log("Configuration error: PBMBMT_GIZAALIGNMENT does not exist!",red)
                    sane = False
                else:
                    os.symlink(self.PBMBMT_GIZAALIGNMENT, self.gets2tfilename('A3.final'))
            if not self.EXEC_PBMBMT_DECODER or not os.path.isfile(self.EXEC_PBMBMT_DECODER):
                sane = False
                self.log("PBMBMT decoder not found! Set EXEC_PBMBMT_DECODER or PATH_PBMBMT !",red)
            if not self.EXEC_PBMBMT_INSTANCEGENERATOR or not os.path.isfile(self.EXEC_PBMBMT_INSTANCEGENERATOR):
                sane = False
                self.log("PBMBMT instance generator not found! Set EXEC_PBMBMT_DECODER or PATH_PBMBMT !",red)
            if not self.EXEC_TIMBL or not os.path.isfile(self.EXEC_TIMBL):
                sane = False
                self.log("TiMBL was not found, but is required for PBMBMT! Set EXEC_TIMBL or PATH_TIMBL !",red)




        if self.BUILD_MOSES_MEMSCORE:
            if not self.BUILD_MOSES_PHRASETRANSTABLE:
                self.log("Configuration update: BUILD_MOSES_PHRASETRANSTABLE automatically enabled because BUILD_MOSES_MEMSCORE is too",yellow)
                self.BUILD_MOSES_PHRASETRANSTABLE = True


        if self.BUILD_COLIBRI_GIZA:
            if not self.BUILD_GIZA_WORDALIGNMENT:
                self.log("Configuration update: BUILD_GIZA_WORDALIGNMENT automatically enabled because BUILD_COLIBRI_GIZA is too",yellow)
                self.BUILD_GIZA_WORDALIGNMENT = True
            if not self.BUILD_GIZA_WORDALIGNMENT_REV:
                self.log("Configuration update: BUILD_GIZA_WORDALIGNMENT_REV automatically enabled because BUILD_COLIBRI_GIZA is too",yellow)
                self.BUILD_GIZA_WORDALIGNMENT_REV = True


        if self.BUILD_GIZA_WORDALIGNMENT and (not self.EXEC_GIZA or not os.path.isfile(self.EXEC_GIZA)):
            self.log("Dependency error: GIZA++ not found (EXEC_GIZA=" + self.EXEC_GIZA + ")",red)
        if self.BUILD_GIZA_WORDALIGNMENT and (not self.EXEC_GIZA_PLAIN2SNT or not os.path.isfile(self.EXEC_GIZA_PLAIN2SNT)):
            self.log("Dependency error: plain2snt.out (provided by GIZA++) not found (EXEC_GIZA_PLAIN2SNT=" + self.EXEC_GIZA_PLAIN2SNT + ")",red)
            sane = False

        if (self.BUILD_SRILM_TARGETMODEL or self.BUILD_SRILM_SOURCEMODEL) and (not self.EXEC_SRILM or not os.path.isfile(self.EXEC_SRILM)):
            self.log("Dependency error: ngram-count (provided by SRILM) not found (EXEC_SRILM=" + self.EXEC_SRILM + ")",red)
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
        print >>sys.stderr,"\tstart                              Train and test the MT system (on testset in configuration)"
        print >>sys.stderr,"\ttrain                              Train the MT system"
        print >>sys.stderr,"\trun <inputfile> [options]          Run the MT system on the specified input file"
        print >>sys.stderr,"\t\t-t                               Tokenise the input file"
        print >>sys.stderr,"\t\t-o <outputfile>                  Output file (default: stdout)"
        print >>sys.stderr,"\tstartserver <port> <html>          Start the MT system as a server on the specified port. Second argument may be set to 1 (default 0) if you want HTML output, which highlights unknown words."
        print >>sys.stderr,"\ttest <inputfile> <referencefile>   Run and evaluate the MT system on the specified input file and reference file (one sentence per line). If <inputfile> and <referencefile> are not given, the default test files from the system configuration are used."
        print >>sys.stderr,"\tscore <inputfile> <referencefile>  Like test, but work on a pre-run system, does not run the translation again."
        print >>sys.stderr,"\tclean [all|giza|moses|colibri|score|batch]    Clean generated files"
        print >>sys.stderr,"\tbranch <expname> [VARIABLE value]       Create a new branch based on this project (files are symlinked instead of copied)"
        print >>sys.stderr,"\tconf VARIABLE value [VARIABLE2 value2]  Change configuration"
        print >>sys.stderr,"\tls    List all batches"
        print >>sys.stderr,"\tstartbatch [batchname] [batchname2]     Start batch (if none are specified, all specified batches will be started sequentially)"
        print >>sys.stderr,"\tbatchreport [batchname] [batchname2]    Write a batch report for the specified batched if none are specified, all specified batches will be included)"
        print >>sys.stderr,"\tbatchtest [batchname] [batchname2]      (Re-)score batches (if none are specified, all specified batches will be tested sequentially)"
        print >>sys.stderr,"\tbatchscore [batchname] [batchname2]     (Re-)score batches (if none are specified, all specified batches will be scored sequentially)"
        print >>sys.stderr,"\rresource [sourcecorpus] [targetcorpus] [trainsize] [testsize] [devsize]    Regenerate training/test/dev sets from specified sources, trainsize may be '*' for all remaining data. Will automatically issue a 'clean all' as well"
        print >>sys.stderr,"\tbatchscore [batchname] [batchname2]     (Re-)score batches (if none are specified, all specified batches will be scored sequentially)"
        print >>sys.stderr,"Note: The batch commands (except batchreport), can run batches in parallel, append the number of threads you want directly to the command (no space)"

    def start(self):
        print >>sys.stderr, repr(sys.argv)
        try:
            os.chdir(self.WORKDIR)
            self.WORKDIR = os.getcwd()
        except:
            if not os.path.exists(sys.argv[0]):
                #not yet in working directory
                print >>sys.stderr, "Unable to find working directory ", self.WORKDIR, ". Path does not exist or not an absolute path?"
        try:
            cmd = sys.argv[1]
        except:
            print >>sys.stderr, "Please specify a command..."
            self.usage()
            sys.exit(2)
        if cmd == 'start':
            self.initlog('train')
            if os.path.isfile(self.WORKDIR + '/.frozen'):
                print >>sys.stderr, "Courageously refusing to train system because it is frozen"
                sys.exit(2)
            if not self.starttrain():
                sys.exit(1)

            self.initlog('test')
            if len(sys.argv) == 4:
                inputfile = sys.argv[2]
                referencefile = sys.argv[3]
            elif len(sys.argv) == 2:
                if not self.TESTSOURCECORPUS or not os.path.exists(self.TESTSOURCECORPUS):
                    self.log("No predefined default test corpus set for input (set TESTSOURCECORPUS) or specify on command line",red,True)
                    sys.exit(2)
                if not self.TESTTARGETCORPUS or not os.path.exists(self.TESTTARGETCORPUS):
                    self.log("No predefined default test corpus set for reference (set TESTSOURCECORPUS and TESTTARGETCORPUS) or specify on command line ",red,True)
                    sys.exit(2)
                inputfile = self.TESTSOURCECORPUS
                referencefile = self.TESTTARGETCORPUS
            else:
                self.usage()
                sys.exit(2)

            if not self.test(inputfile, referencefile):
                sys.exit(1)

        elif cmd == 'train':
            self.initlog('train')
            if os.path.isfile(self.WORKDIR + '/.frozen'):
                print >>sys.stderr, "Courageously refusing to train system because it is frozen"
                sys.exit(2)
            if not self.starttrain():
                sys.exit(1)
        elif cmd == 'clean':
            self.initlog('clean')
            if os.path.isfile(self.WORKDIR + '/.frozen'):
                print >>sys.stderr, "Courageously refusing to clean system because it is frozen"
                sys.exit(2)
            targets = sys.argv[2:]
            if not self.clean(targets):
                sys.exit(1)
        elif cmd == 'run':
            self.initlog('run')
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

        elif cmd == 'startserver':
            self.initlog('server')
            port = int(sys.argv[2])
            if len(sys.argv) > 3 and sys.argv[3]:
                html = True
            else:
                html = False
            self.server(port, html)

        elif cmd == 'freeze':
            self.initlog('freeze')
            if os.path.isfile(self.WORKDIR + '/.frozen'):
                print >>sys.stderr, "System is already frozen"
                sys.exit(0)


            open(self.WORKDIR + '/.frozen','w').close()
            print >>sys.stderr, "System frozen"


        elif cmd == 'unfreeze' or cmd=='thaw':
            self.initlog('unfreeze')
            if os.path.isfile(self.WORKDIR + '/.frozen'):
                os.unlink(self.WORKDIR + '/.frozen')
                print >>sys.stderr, "Thawing system"
                sys.exit(0)
            else:
                print >>sys.stderr, "System is not frozen"
                sys.exit(0)

        elif cmd == 'branch':
            self.initlog('branch')
            expname = sys.argv[2]

            branchdir, branchsettings = self.branch(expname, sys.argv[3:])

            if 'EDITOR' in os.environ:
                editor = os.environ['EDITOR']
            else:
                editor = 'vim'
            os.system(editor + ' ' + branchsettings)

            open(self.WORKDIR + '/.frozen','w').close()
            print >>sys.stderr, "Current system automatically frozen after branching"
            print >>sys.stderr,"All done, to go to the newly branched system: $ cd " + branchdir

        elif cmd == 'test':
            self.initlog('test')
            if len(sys.argv) == 4:
                inputfile = sys.argv[2]
                referencefile = sys.argv[3]
            elif len(sys.argv) == 2:
                if not self.TESTSOURCECORPUS or not os.path.exists(self.TESTSOURCECORPUS):
                    self.log("No predefined default test corpus set for input (set TESTSOURCECORPUS) or specify on command line",red,True)
                    sys.exit(2)
                if not self.TESTTARGETCORPUS or not os.path.exists(self.TESTTARGETCORPUS):
                    self.log("No predefined default test corpus set for reference (set TESTSOURCECORPUS and TESTTARGETCORPUS) or specify on command line ",red,True)
                    sys.exit(2)
                inputfile = self.TESTSOURCECORPUS
                referencefile = self.TESTTARGETCORPUS
            else:
                self.usage()
                sys.exit(2)

            if not self.test(inputfile, referencefile):
                sys.exit(1)

        elif cmd == 'score':
            self.initlog('score')
            if len(sys.argv) == 4:
                inputfile = sys.argv[2]
                referencefile = sys.argv[3]
            elif len(sys.argv) == 2:
                if not self.TESTSOURCECORPUS or not os.path.exists(self.TESTSOURCECORPUS):
                    self.log("No predefined default test corpus set for input (set TESTSOURCECORPUS) or specify on command line",red,True)
                    sys.exit(2)
                if not self.TESTTARGETCORPUS or not os.path.exists(self.TESTTARGETCORPUS):
                    self.log("No predefined default test corpus set for reference (set TESTSOURCECORPUS and TESTTARGETCORPUS) or specify on command line ",red,True)
                    sys.exit(2)
                inputfile = self.TESTSOURCECORPUS
                referencefile = self.TESTTARGETCORPUS
            else:
                self.usage()
                sys.exit(2)

            if not self.score(inputfile, referencefile,'output.txt'):
                sys.exit(1)
        elif cmd == 'conf':
            self.initlog('conf')
            if self.parseconf(sys.argv[2:]):
                self.log("Writing new configuration...")
                self.writesettings()
        elif cmd[:10] == 'startbatch' or cmd[:10] == 'batchstart':
            if len(cmd) > 10 and cmd[10] != ' ':
                space = cmd.find(' ',10)
                if space > -1:
                    threads = int(cmd[10:space])
                else:
                    threads = int(cmd[10:])
            else:
                threads = 1

            self.initlog('batch')
            if not self.batches:
                self.log("No batch jobs in configuration...",red)
                sys.exit(2)

            self.log("Threads: " + str(threads))

            if not os.path.isfile(self.WORKDIR + '/.frozen'):
                self.log("Refusing to start batches from a non-frozen system, please explicitly freeze the system first",red)
                sys.exit(2)

            if sys.argv[2:]:
                selectedbatches = sys.argv[2:]
                for batch in selectedbatches:
                    if not batch in [ key for (key,conf) in self.batches ]:
                        self.log( "No such batch: " + batch,red,True)
            else:
                selectedbatches= None


            xpool = ExperimentPool(threads)
            for batch, conf in self.batches:
                if not selectedbatches or batch in selectedbatches:
                    xpool.append(BatchExperiment( (self, batch, conf, True, True ) ) )
            for x in xpool.run(False):
                pass
            self.log("Done")
        elif cmd == 'ls':
            for batch in self.batches:
                print batch[0] #, "\t", batch[1]
        elif cmd == 'batchconf':


            if sys.argv[2:]:
                selectedbatches = []
                confargs = []
                for x in sys.argv[2:]:
                    if x in self.batches:
                        selectedbatches.append(x)
                    else:
                        confargs.append(x)
            else:
                selectedbatches= None
                confargs = []


            self.initlog('conf')
            if self.parseconf(sys.argv[2:]):
                self.log("Writing new configuration...")
                self.writesettings()
            for batch, conf in self.batches:
                if not selectedbatches or batch in selectedbatches:
                    batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                    r = os.system(batchdir + '/mt-' +  self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch + '.py conf ' + ' '.join(confargs))
                    if r == 0:
                        self.log("Configuring batch " + batch + " finished succesfully " + self.timestamp(),green,True)
                    else:
                        self.log("Configuring batch " + batch + " finished with error code " + str(r) + " " + self.timestamp(),red,True)
        elif cmd == 'batchscore':
            self.initlog('batchscore')
            if not self.batches:
                self.log("No batch jobs in configuration...",red)
                sys.exit(2)

            if not os.path.isfile(self.WORKDIR + '/.frozen'):
                self.log("Refusing to start batches from a non-frozen system, please explicitly freeze the system first",red)
                sys.exit(2)

            if sys.argv[2:]:
                selectedbatches = sys.argv[2:]
                for batch in selectedbatches:
                    if not batch in [ key for (key,conf) in self.batches ]:
                        self.log( "No such batch: " + batch,red,True)
            else:
                selectedbatches= None

            for batch, conf in self.batches:
                if not selectedbatches or batch in selectedbatches:
                    batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                    if os.path.exists(batchdir):
                        self.log("Starting scoring batch " + batch + " " + self.timestamp(),white,True)
                        os.system(batchdir + '/mt-' +  self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch + '.py score')
                    else:
                        self.log("Batch " + batch + " has not been trained or tested yet.. skipping",yellow,True)
        elif cmd[:9] == 'batchtest':
            if len(cmd) > 9 and cmd[9] != ' ':
                space = cmd.find(' ',9)
                if space > -1:
                    threads = int(cmd[9:space])
                else:
                    threads = int(cmd[9:])
            else:
                threads = 1

            self.initlog('batchtest')
            if not self.batches:
                self.log("No batch jobs in configuration...",red)
                sys.exit(2)

            self.log("Threads: " + str(threads))

            if not os.path.isfile(self.WORKDIR + '/.frozen'):
                self.log("Refusing to start batches from a non-frozen system, please explicitly freeze the system first",red)
                sys.exit(2)

            if sys.argv[2:]:
                selectedbatches = sys.argv[2:]
                for batch in selectedbatches:
                    if not batch in [ key for (key,conf) in self.batches ]:
                        self.log( "No such batch: " + batch,red,True)
            else:
                selectedbatches= None

            # for batch, conf in self.batches:
                # if not selectedbatches or batch in selectedbatches:
                    # batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                    # if os.path.exists(batchdir):
                        # self.log("Starting scoring batch " + batch + " " + self.timestamp(),white,True)
                        # rtrain = os.system(batchdir + '/mt-' +  self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch + '.py test')
                    # else:
                        # self.log("Batch " + batch + " has not been trained yet.. skipping",yellow,True)

            xpool = ExperimentPool(threads)
            for batch, conf in self.batches:
                if not selectedbatches or batch in selectedbatches:
                     batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                     if os.path.exists(batchdir):
                        xpool.append(BatchExperiment( (self, batch, conf, False, True ) ) )
                     else:
                        self.log("Batch " + batch + " has not been trained yet.. skipping",yellow,True)
            for x in xpool.run(False):
                pass
            self.log("Done")

        elif cmd == 'batchreport':
            self.initlog('batchreport')
            if not self.batches:
                self.log("No batch jobs in configuration...",red)
                sys.exit(2)

            if not os.path.isfile(self.WORKDIR + '/.frozen'):
                self.log("Refusing to start batches from a non-frozen system, please explicitly freeze the system first",red)
                sys.exit(2)

            if sys.argv[2:]:
                selectedbatches = sys.argv[2:]
                for batch in selectedbatches:
                    if not batch in [ key for (key,conf) in self.batches ]:
                        self.log( "No such batch: " + batch,red,True)
            else:
                selectedbatches= None

            self.batchreport(selectedbatches)
            self.log("Done")

        elif cmd == 'resource' or cmd == 'resource':
            try:
                sourcecorpusfile = sys.argv[2]
                targetcorpusfile = sys.argv[3]
                trainsize = int(sys.argv[4])
                testsize = int(sys.argv[5])
                devsize = int(sys.argv[6])
            except:
                self.log("Invalid arguments for generate",red)
                sys.exit(2)

            if not self.clean('all'):
                sys.exit(1)

            sourcecorpusfile, targetcorpusfile, testsourcecorpusfile, testtargetcorpusfile, devsourcecorpusfile, devtargetcorpusfile = resource(sourcecorpusfile, targetcorpusfile, testsize, devsize, trainsize, self.WORKDIR, self.SOURCELANG, self.TARGETLANG, self.CORPUSNAME)
            if sourcecorpusfile != None:  self.TRAINSOURCECORPUS = sourcecorpusfile
            if targetcorpusfile != None:  self.TRAINTARGETCORPUS = targetcorpusfile
            if testsourcecorpusfile != None:  self.TESTSOURCECORPUS = testsourcecorpusfile
            if testtargetcorpusfile != None:  self.TESTTARGETCORPUS = testtargetcorpusfile
            if devsourcecorpusfile != None:  self.DEVSOURCECORPUS = devsourcecorpusfile
            if devtargetcorpusfile != None:  self.DEVTARGETCORPUS = devtargetcorpusfile

            self.writesettings()
            self.log("Done",green)

        elif cmd == 'help' or cmd == '-h':
            self.usage()
        else:
            self.log("Error, no such command: " + cmd,red)
            self.usage()
            sys.exit(2)

        sys.exit(0)

    def batchreport(self, selectedbatches):
        if self.EXPERIMENTNAME:
            title = self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + self.EXPERIMENTNAME
        else:
            title = self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG

        scores = []
        names = []
        phrasetablesize = []
        names_phrasetable = []
        for batch, conf in self.batches:
            if not selectedbatches or batch in selectedbatches:
                self.log("Gathering scores from " + batch)
                batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                if os.path.isfile(batchdir + '/summary.score'):
                    f = open(batchdir + '/summary.score','r')
                    f.readline() #skip table header
                    blue, meteor, nist, ter, wer, per = [ float(x) for x in f.readline().split() ]
                    scores.append( ( blue, meteor, nist, ter, wer, per) )
                    names.append(batch)
                    f.close()
                if os.path.isfile(batchdir + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG  + '.alignmodel.colibri'):
                    os.system(self.EXEC_COLIBRI_ALIGNER + " -d " + batchdir + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG  + '.alignmodel.colibri' + " --stats > " + batchdir + "/alignmodel.stats")
                    f = open(batchdir + "/alignmodel.stats")
                    count = uniquecount = 0
                    for line in f:
                        if line[:17] == "sources aligned: ":
                            uniquecount = int(line[17:])
                        elif line[:18] == "total alignments: ":
                            count = int(line[18:])
                    f.close()
                    phrasetablesize.append( (count, uniquecount) )
                    names_phrasetable.append(batch)
                elif os.path.isfile(batchdir + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG  + '.phrasetable'):
                    count = 0
                    uniquecount = 0
                    f = open(batchdir + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG  + '.phrasetable')
                    prev = ''
                    for line in f:
                        sourcephrase = line.split('|||')[0].strip()
                        if sourcephrase != prev:
                            uniquecount += 1
                        prev = sourcephrase
                        count += 1
                    f.close()
                    phrasetablesize.append( (count, uniquecount) )
                    names_phrasetable.append(batch)

        if not scores:
            self.log("Error, no scores found! Did you forget to train/test/score the batches first?", red)
            sys.exit(2)

        f = open(self.WORKDIR + '/' + title + '.batchreport.tex','w')
        f.write("\\documentclass[a4paper,10pt]{article}\n\\usepackage{graphicx}\n")
        f.write("\\begin{document}\n\n")
        f.write('\\section*{Results for ' + title + '}\n\n')
        f.write("\\begin{table}[h]\n")
        f.write("\\begin{tabular}{|l|c|c|c|c|c|c|}\n")
        f.write("\\hline \n Name & BLEU & METEOR & NIST & TER & WER & PER \\\\ \n \\hline \n")
        bestbleu = max( [ x[0] for x in scores ] )
        bestmeteor = max( [ x[1] for x in scores ] )
        bestnist = max( [ x[2] for x in scores ] )
        bestter = min( [ x[3] for x in scores ] )
        bestwer = min( [ x[4] for x in scores ] )
        bestper = min( [ x[5] for x in scores ] )
        for name, ( bleu, meteor, nist, ter, wer, per)  in zip(names,scores):
               f.write('\\textbf{' + name + '} & ')
               if bleu == bestbleu:
                    f.write( '\\textbf{' + '%.3f' % bleu + '} &')
               else:
                    f.write( '%.3f' % bleu + ' &')
               if meteor == bestmeteor:
                    f.write( '\\textbf{' + '%.3f' % meteor + '} &')
               else:
                    f.write( '%.3f' % meteor + ' &')
               if nist == bestnist:
                    f.write( '\\textbf{' + '%.3f' % nist + '} &')
               else:
                    f.write( '%.3f' % nist + ' &')
               if ter == bestter:
                    f.write( '\\textbf{' + '%.2f' % ter + '} &')
               else:
                    f.write( '%.3f' % ter + ' &')
               if wer == bestwer:
                    f.write( '\\textbf{' + '%.2f' % wer + '} &')
               else:
                    f.write( '%.3f' % wer + ' &')
               if per == bestper:
                    f.write( '\\textbf{' + '%.2f' % per + '} \\\\ \n')
               else:
                    f.write( '%.3f' % per + ' \\\\ \n')
               f.write('\\hline')
        f.write("\\end{tabular}\n")
        f.write("\\caption{Evaluation results for batches in " + title + "}\n")
        f.write("\\end{table}\n")
        f.write("\\begin{figure}\n")
        f.write("\\includegraphics[width=18cm]{batchreport-bleu.png}\n")
        f.write("\\caption{BLEU scores for " + title + "}\n")
        f.write("\\end{figure}\n")
        f.write("\\begin{figure}\n")
        f.write("\\includegraphics[width=18cm]{batchreport-meteor.png}\n")
        f.write("\\caption{METEOR scores for " + title + "}\n")
        f.write("\\end{figure}\n")
        f.write("\\begin{figure}\n")
        f.write("\\includegraphics[width=18cm]{batchreport-nist.png}\n")
        f.write("\\caption{NIST scores for " + title + "}\n")
        f.write("\\end{figure}\n")
        f.write("\\begin{figure}\n")
        f.write("\\includegraphics[width=18cm]{batchreport-er.png}\n")
        f.write("\\caption{TER/WER/PER scores for " + title + "}\n")
        f.write("\\end{figure}\n")
        f.write("\\begin{figure}\n")
        f.write("\\includegraphics[width=18cm]{batchreport-phrasetable.png}\n")
        f.write("\\caption{Phrasetable size for " + title + "}\n")
        f.write("\\end{figure}\n")
        f.write("\\end{document}\n")
        f.close()


        def autolabel(rects):
            # attach some text labels
            for rect in rects:
                height = rect.get_height()
                if height == round(height):
                    matplotlib.pyplot.text(rect.get_x()+rect.get_width()/2., height+0.001, '%d'%height,ha='center', va='bottom')
                else:
                    matplotlib.pyplot.text(rect.get_x()+rect.get_width()/2., height+0.001, '%.3f'%height,ha='center', va='bottom')


        def autolabelh(rects):
            # attach some text labels
            for rect in rects:
                width = rect.get_width()
                if width == round(width):
                    matplotlib.pyplot.text(width+0.01, rect.get_y()+rect.get_height()/2., '%d'%width,ha='center', va='center')
                else:
                    matplotlib.pyplot.text(width+0.01, rect.get_y()+rect.get_height()/2., '%.3f'%width,ha='center', va='center')




        hbarheight = 0.2

        locations = numpy.arange(len(scores))    # the x locations for the groups

        fig = matplotlib.pyplot.figure(figsize=(15,10))
        matplotlib.pyplot.gcf().subplots_adjust(bottom=0.2, left=0.4)
        matplotlib.pyplot.grid(True)
        p_bleu = matplotlib.pyplot.barh(locations ,  [x[0] for x in scores], align='center', color='b')
        matplotlib.pyplot.ylabel('BLEU score')
        matplotlib.pyplot.title('BLEU scores for ' + title)
        matplotlib.pyplot.yticks(locations+hbarheight/2., names)# size='small')
        matplotlib.pyplot.xticks(numpy.arange(0,max( (x[0] for x in scores)),0.01))
        autolabelh(p_bleu)
        fig.savefig(self.WORKDIR + '/batchreport-bleu.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        matplotlib.pyplot.clf()

        fig = matplotlib.pyplot.figure(figsize=(15,10))
        matplotlib.pyplot.gcf().subplots_adjust(bottom=0.2, left=0.4)
        matplotlib.pyplot.grid(True)
        p_meteor = matplotlib.pyplot.barh(locations ,  [x[1] for x in scores], align='center', color='y')
        matplotlib.pyplot.ylabel('METEOR score')
        matplotlib.pyplot.title('METEOR scores for ' + title)
        matplotlib.pyplot.yticks(locations+hbarheight/2., names)# size='small')
        matplotlib.pyplot.xticks(numpy.arange(0,max( (x[1] for x in scores)),0.05))
        autolabelh(p_meteor)
        fig.savefig(self.WORKDIR + '/batchreport-meteor.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        matplotlib.pyplot.clf()

        fig = matplotlib.pyplot.figure(figsize=(15,10))
        matplotlib.pyplot.gcf().subplots_adjust(bottom=0.2, left=0.4)
        matplotlib.pyplot.grid(True)
        p_nist = matplotlib.pyplot.barh(locations ,  [x[2] for x in scores], align='center', color='y')
        matplotlib.pyplot.ylabel('NIST score')
        matplotlib.pyplot.title('NIST scores for ' + title)
        matplotlib.pyplot.yticks(locations+hbarheight/2., names)# size='small')
        matplotlib.pyplot.xticks(numpy.arange(0,max( (x[2] for x in scores)),0.1))
        autolabelh(p_nist)
        fig.savefig(self.WORKDIR + '/batchreport-nist.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        matplotlib.pyplot.clf()



        fig = matplotlib.pyplot.figure(figsize=(15,10))
        matplotlib.pyplot.gcf().subplots_adjust(bottom=0.2, left=0.4)
        matplotlib.pyplot.grid(True)
        p_ter = matplotlib.pyplot.barh(locations ,  [x[3] for x in scores], align='center', color='g')
        p_wer = matplotlib.pyplot.barh(locations ,  [x[4] for x in scores], align='center', color='y')
        p_per = matplotlib.pyplot.barh(locations ,  [x[5] for x in scores], align='center', color='m')
        matplotlib.pyplot.ylabel('Error rate')
        matplotlib.pyplot.title('Error rates for ' + title)
        matplotlib.pyplot.yticks(locations+hbarheight/2., names)# size='small')
        matplotlib.pyplot.xticks(numpy.arange(0, max( ( max(x[3],x[4],x[5]) for x in scores ) ) ,5) )

        matplotlib.pyplot.legend( (p_ter[0],p_wer[0],p_per[0]), ('TER', 'WER','PER') )
        fig.savefig(self.WORKDIR + '/batchreport-er.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        matplotlib.pyplot.clf()


        locations2 = numpy.arange(len(phrasetablesize))    # the x locations for the groups

        fig = matplotlib.pyplot.figure(figsize=(15,10))
        matplotlib.pyplot.gcf().subplots_adjust(bottom=0.2, left=0.4)
        matplotlib.pyplot.grid(True)
        try:
            p_count = matplotlib.pyplot.barh(locations2 ,  [x[0] for x in phrasetablesize], align='center', color='b')
            p_uniquecount = matplotlib.pyplot.barh(locations2 ,  [x[1] for x in phrasetablesize], align='center', color='m')
            matplotlib.pyplot.ylabel('Number of phrase pairs in phrase table')
            matplotlib.pyplot.title('Phrase table size for ' + title)
            matplotlib.pyplot.legend( (p_count[0],p_uniquecount[0]), ('Count', 'Unique') )

            matplotlib.pyplot.yticks(locations2+hbarheight/2., names_phrasetable )# size='small')
            matplotlib.pyplot.xticks(numpy.arange(0,max( (x[0] for x in phrasetablesize) ),100000))
            autolabelh(p_count)
            fig.savefig(self.WORKDIR + '/batchreport-phrasetable.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        except:
            self.log("Error during plotting of phrasetable size",red)
        matplotlib.pyplot.clf()



        # width = 0.9       # the width of the bars: can also be len(x) sequence
        # fig = matplotlib.pyplot.figure(figsize=(15,10))
        # matplotlib.pyplot.gcf().subplots_adjust(bottom=0.4, left=0.2)
        # xlocations = numpy.arange(len(scores))    # the x locations for the groups


        # matplotlib.pyplot.grid(True)
        # p_bleu = matplotlib.pyplot.bar(xlocations, [ x[0] for x in scores] ,  width, color='b')
        # matplotlib.pyplot.ylabel('BLEU score')
        # matplotlib.pyplot.title('BLEU scores for ' + title)
        # matplotlib.pyplot.xticks(xlocations+width/2., names)# size='small')
        # fig.autofmt_xdate()
        # matplotlib.pyplot.yticks(numpy.arange(0,max( (x[0] for x in scores)),0.01))
        # autolabel(p_bleu)
        # fig.savefig(self.WORKDIR + '/batchreport-bleu.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        # matplotlib.pyplot.clf()

        # fig = matplotlib.pyplot.figure(figsize=(15,10))
        # matplotlib.pyplot.grid(True)
        # p_bleu = matplotlib.pyplot.bar(xlocations, [ x[1] for x in scores] ,  width, color='y')
        # matplotlib.pyplot.ylabel('METEOR score')
        # matplotlib.pyplot.title('METEOR scores for ' + title )
        # matplotlib.pyplot.xticks(xlocations+width/2., names)# size='small')
        # fig.autofmt_xdate()
        # matplotlib.pyplot.yticks(numpy.arange(0,max( (x[1] for x in scores)),0.05))
        # autolabel(p_bleu)
        # fig.savefig(self.WORKDIR + '/batchreport-meteor.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        # matplotlib.pyplot.clf()

        # fig = matplotlib.pyplot.figure(figsize=(15,10))
        # matplotlib.pyplot.grid(True)
        # p_bleu = matplotlib.pyplot.bar(xlocations, [ x[2] for x in scores] ,  width, color='r')
        # matplotlib.pyplot.ylabel('NIST score')
        # matplotlib.pyplot.title('NIST scores for ' + title )
        # matplotlib.pyplot.xticks(xlocations+width/2., names)# size='small')
        # fig.autofmt_xdate()
        # matplotlib.pyplot.yticks(numpy.arange(0,max( (x[2] for x in scores)),0.1))
        # autolabel(p_bleu)
        # fig.savefig(self.WORKDIR + '/batchreport-nist.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        # matplotlib.pyplot.clf()


        # fig = matplotlib.pyplot.figure(figsize=(15,10))
        # matplotlib.pyplot.gcf().subplots_adjust(bottom=0.4, left=0.2)
        # p_ter = matplotlib.pyplot.bar(xlocations, [x[3] for x in scores] ,  width, color='g')
        # p_wer = matplotlib.pyplot.bar(xlocations, [x[4] for x in scores] ,  width, color='y')
        # p_per = matplotlib.pyplot.bar(xlocations, [x[5] for x in scores] ,  width, color='m')
        # matplotlib.pyplot.ylabel('Score')
        # matplotlib.pyplot.title('TER/WER/PER scores for ' + title)
        # matplotlib.pyplot.xticks(xlocations+width/2., names)#, names , rotation='vertical')
        # matplotlib.pyplot.yticks(numpy.arange(0,max( (max(x[3],x[4],x[5]) for x in scores)),5))
        # matplotlib.pyplot.legend( (p_ter[0],p_wer[0],p_per[0]), ('TER', 'WER','PER') )
        # fig.autofmt_xdate()
        # fig.savefig(self.WORKDIR + '/batchreport-er.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)

        # fig = matplotlib.pyplot.figure(figsize=(15,10))
        # xlocations2 = numpy.arange(len(phrasetablesize))    # the x locations for the groups
        # matplotlib.pyplot.grid(True)
        # p_count = matplotlib.pyplot.bar(xlocations2, [ x[0] for x in phrasetablesize ] ,  width, color='b')
        # p_uniquecount = matplotlib.pyplot.bar(xlocations2, [ x[1] for x in phrasetablesize ] ,  width, color='m')
        # matplotlib.pyplot.ylabel('Number of phrase pairs in phrase table')
        # matplotlib.pyplot.title('Phrase table size for ' + title)
        # matplotlib.pyplot.xticks(xlocations2+width/2., names_phrasetable)# size='small')
        # fig.autofmt_xdate()
        # try:
            # matplotlib.pyplot.yticks(numpy.arange(0,max( (x[0] for x in phrasetablesize)),100000))
            # matplotlib.pyplot.legend( (p_count[0],p_uniquecount[0]), ('Count', 'Unique') )
            # autolabel(p_count)
            # fig.savefig(self.WORKDIR + '/batchreport-phrasetable.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.3)
        # except:
            # self.log("Error during plotting of phrasetable size",red)
        # matplotlib.pyplot.clf()





    def clean(self, targets):
        if not targets:
            self.log("Nothing to clean, please specify one or more targets: all, giza, berkeleyaligner, moses, colibri, classifiers, mert, phrasal, srilm, test, score, batches",red)
            sys.exit(2)

        if 'giza' in targets or 'all' in targets:
            self.cleanfiles('corpus/*','*.final', '*.vcb','*.snt','*.classes','*.classes.cats','*.gizacfg','*.Decoder.config','*.perp','*.cooc')
        if 'berkeleyaligner' in targets or 'all' in targets:
            self.cleanfiles('*.final','aligner.conf')
            if os.path.isdir(self.WORKDIR + '/alignerinput'):
                shutil.rmtree(self.WORKDIR + '/alignerinput')
                self.log("Removed alignerinput",green)
            if os.path.isdir(self.WORKDIR + '/aligneroutput'):
                shutil.rmtree(self.WORKDIR + '/aligneroutput')
                self.log("Removed aligneroutput",green)
        if 'moses' in targets or 'all' in targets:
            self.cleanfiles('*.bal', '*.symal','*.s2t','*.s2t.sorted','*.t2s','*.t2s','*.sorted','*.phrasetable', '*.phraseextract', '*.phraseextract.inv','*.half','moses.ini','model/*')
        if 'mert' in targets or 'all' in targets:
            self.cleanfiles('mert-work/*')
        if 'classifier' in targets or 'classifiers' in targets or 'timbl' in targets  or 'all' in targets:
            self.cleanfiles('classifier.*')
        if 'phrasal' in targets or 'all' in targets:
            self.cleanfiles('phrases*gz','*.phrasetable','phrasal.conf')
        if 'srilm' in targets or 'all' in targets:
            self.cleanfiles('*.srilm')
        if 'colibri' in targets or 'all' in targets:
            self.cleanfiles('*.colibri','*.cls','*.clsenc')
        if 'test' in targets or 'all' in targets:
            self.cleanfiles('output.txt','*.score')
        if 'score' in targets or 'all' in targets:
            self.cleanfiles('*.score')
        if 'batches' in targets:
            if not self.batches:
                self.log("No batch jobs in configuration...",red)
                sys.exit(2)

            if sys.argv[3:]:
                selectedbatches = sys.argv[2:]
                for batch in selectedbatches:
                    if not batch in [ key for (key,conf) in self.batches ]:
                        self.log( "No such batch: " + batch,red,True)


            for batch, conf in self.batches:
                if not selectedbatches or batch in selectedbatches:
                    self.log("Removing batch " + batch,white,True)
                    batchdir = self.WORKDIR + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + batch
                    if os.path.isdir(batchdir):
                        shutil.rmtree(batchdir)

        return True

    def cleanfiles(self, *args):
        ok = True
        for mask in args:
            for filename in glob.glob(self.WORKDIR + '/' + mask):
                try:
                    os.unlink(filename)
                    self.log("Removed " + filename,green)
                except:
                    ok = False
                    self.log("Unable to remove " + filename,red,True)
        return ok

    def starttrain(self):
        self.init()
        if not self.check_common(): return False
        if not self.check_train(): return False

        if self.TOKENIZE_SOURCECORPUS and not self.tokenize_sourcecorpus(): return False
        if self.TOKENIZE_TARGETCORPUS and not self.tokenize_targetcorpus(): return False

        if self.TRAINLIMIT and not self.trainlimit(): return False

        if self.BUILD_SRILM_TARGETMODEL and not self.build_srilm_targetmodel(): return False
        if self.BUILD_SRILM_SOURCEMODEL and not self.build_srilm_sourcemodel(): return False

        if self.BUILD_GIZA_WORDALIGNMENT and not self.build_giza_wordalignment(): return False
        if self.BUILD_GIZA_WORDALIGNMENT_REV and not self.build_giza_wordalignment_rev(): return False


        if self.BUILD_PHRASAL_WORDALIGN and not self.build_phrasal_wordalign(): return False

        if self.BUILD_COLIBRI_ALIGNMENT and not self.build_colibri_alignment(): return False


        if self.BUILD_PHRASAL_PHRASEEXTRACT and not self.build_phrasal_phraseextract(): return False

        #if self.BUILD_MOSES_SYMAL and not self.build_moses_symal(): return False
        #if self.BUILD_MOSES_WORDTRANSTABLE and not self.build_moses_wordtranstable(): return False
        #if self.BUILD_MOSES_PHRASEEXTRACT and not self.build_moses_phraseextract(): return False
        #if self.BUILD_MOSES_PHRASETRANSTABLE and not self.build_moses_phrasescore(): return False

        #TODO: Moses reordering model and generation model

        if self.BUILD_PHRASAL and not self.build_phrasal(): return False
        if self.BUILD_PHRASAL_MERT and not self.build_phrasal_mert(): return False

        if self.BUILD_MOSES and not self.build_moses(): return False
        if self.BUILD_MOSES_CLASSIFIERS and not self.build_moses_classifiers(): return False
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
                self.log("Skipping " + name, yellow, True)
                return False
        if 'cmd' in kwargs:
            self.log("Calling " + name + " " + self.timestamp() ,white, True)
            self.log("Command "+ ": " + kwargs['cmd'])
        else:
            self.log("Calling " + name + " " + self.timestamp(),white, True)
        return True

    def timestamp(self):
        return "\t" + magenta("@" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    def footer(self, name, r, *outputfiles, **kwargs):
        if 'successcodes' in kwargs:
            successcodes = kwargs['successcodes']
        else:
            successcodes = [0]
        if r in successcodes:
            self.log("Finished " + name + " " + self.timestamp(),green,True)
        else:
            self.log("Runtime error from " + name + ' (return code ' + str(r) + ') ' + self.timestamp(),red,True)
            return False
        if outputfiles:
            error = False
            for outputfile in outputfiles:
                if os.path.exists(outputfile):
                    self.log("Produced output file " + outputfile,green)
                else:
                    self.log("Expected output file " + outputfile+ ", not produced!",red)
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

    def parseconf(self, conf):
        changed = False
        if isinstance(conf,str) or isinstance(conf,unicode):
            conf = conf.split(' ')
        if isinstance(conf,list) or isinstance(conf, tuple):
            for i in range(0,len(conf) - 1):
                key = conf[i]
                value = conf[i+1]
                if value == "True": value = True
                if value == "False": value = False
                if value == "None": value = None
                try:
                    getattr(self, key)
                    setattr(self,key,value)
                    changed = True
                    self.log("Set variable " + key + " to \"" + str(value) +"\"")
                except AttributeError:
                    self.log("ERROR: No such variable exists: " + key,red)
                    sys.exit(2)
        elif isinstance(conf,dict):
            for key, value in conf.items():
                try:
                    getattr(self, key)
                    setattr(self,key,value)
                    changed = True
                    self.log("Set variable " + key + " to \"" + str(value) +"\"")
                except AttributeError:
                    self.log("ERROR: No such variable exists: " + key,red)
                    sys.exit(2)
        if changed:
            self.log("Configuration updated")
        else:
            self.log("Configuration unchanged")
        return changed


    #---------------------------------- Methods for building sub-parts ----------------------------

    def build_colibri_alignment(self):
        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.getsourcefilename('txt'), "Encoding source corpus for Colibri",self.getsourcefilename('cls'), self.getsourcefilename('clsenc') ): return False

        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.gettargetfilename('txt'), "Encoding target corpus for Colibri",self.gettargetfilename('cls'), self.gettargetfilename('clsenc') ): return False

        patternfinder_extraoptions = ''
        if self.BUILD_COLIBRI_SKIPGRAMS:
            if self.COLIBRI_PATTERNFINDER_OPTIONS.find('-s') == -1:
                patternfinder_extraoptions += ' -s'
            if self.COLIBRI_PATTERNFINDER_OPTIONS.find('-B') == -1:
                patternfinder_extraoptions += ' -B'
            if self.COLIBRI_PATTERNFINDER_OPTIONS.find('-E') == -1:
                patternfinder_extraoptions += ' -E'

        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.getsourcefilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS + patternfinder_extraoptions, "Building source-side pattern model",self.getsourcefilename('indexedpatternmodel.colibri') ): return False

        if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.gettargetfilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS + patternfinder_extraoptions, "Building target-side pattern model",self.gettargetfilename('indexedpatternmodel.colibri') ): return False


        grapher_extraoptions = ''
        if self.BUILD_COLIBRI_SKIPGRAMS:
            if self.COLIBRI_PATTERNFINDER_OPTIONS.find('-a') == -1 and self.COLIBRI_PATTERNFINDER_OPTIONS.find('-T') == -1:
                grapher_extraoptions += ' -T'
            if self.COLIBRI_PATTERNFINDER_OPTIONS.find('-a') == -1 and self.COLIBRI_PATTERNFINDER_OPTIONS.find('-S') == -1:
                grapher_extraoptions += ' -S'

        alignsource = 'indexedpatternmodel.colibri'
        if self.BUILD_COLIBRI_GRAPHMODEL:
            alignsource = 'graphpatternmodel.colibri'
            if not self.runcmd(self.EXEC_COLIBRI_GRAPHER + ' -f ' + self.getsourcefilename('indexedpatternmodel.colibri') + ' ' + self.COLIBRI_GRAPHER_OPTIONS + grapher_extraoptions, "Building source-side graph model",self.getsourcefilename('graphpatternmodel.colibri') ): return False

            if not self.runcmd(self.EXEC_COLIBRI_GRAPHER + ' -f ' + self.gettargetfilename('indexedpatternmodel.colibri') + ' ' + self.COLIBRI_GRAPHER_OPTIONS + grapher_extraoptions, "Building target-side graph model",self.gettargetfilename('graphpatternmodel.colibri') ): return False

        aligner_extraoptions = ''
        if self.BUILD_COLIBRI_CLASSIFIERS:
            aligner_extraoptions = ' -l ' + str(self.COLIBRI_LEFTCONTEXTSIZE) + ' -r ' + str(self.COLIBRI_RIGHTCONTEXTSIZE)

        if self.BUILD_COLIBRI_GIZA:
            if '-H' in self.COLIBRI_ALIGNER_OPTIONS:
                #unsupervised
                if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -o ' + self.gets2tfilename('alignmodel.colibri') + ' -I 2 -N ' + self.COLIBRI_ALIGNER_OPTIONS + ' ' + aligner_extraoptions + ' -W ' + self.gets2tfilename('A3.final') + ':' + self.gett2sfilename('A3.final') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls'), "Building alignment model",self.gets2tfilename('alignmodel.colibri') ): return False
            else:
                #semi-supervised
                if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -s ' + self.getsourcefilename(alignsource) + ' -t ' + self.gettargetfilename(alignsource) + ' -o ' + self.gets2tfilename('alignmodel.colibri') + ' -I 2 -N ' + self.COLIBRI_ALIGNER_OPTIONS + ' ' + aligner_extraoptions + ' -W ' + self.gets2tfilename('A3.final') + ':' + self.gett2sfilename('A3.final') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls'), "Building alignment model",self.gets2tfilename('alignmodel.colibri') ): return False
            #if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -t ' + self.getsourcefilename('graphpatternmodel.colibri') + ' -s ' + self.gettargetfilename('graphpatternmodel.colibri') + ' -o ' + self.gett2sfilename('alignmodel.colibri') + ' -N ' + self.COLIBRI_ALIGNER_OPTIONS + ' -W ' + self.gett2sfilename('A3.final') + ':' + self.gets2tfilename('A3.final') + ' -T ' + self.getsourcefilename('cls') + ' -S ' + self.gettargetfilename('cls'), "Building alignment model",self.gett2sfilename('alignmodel.colibri') ): return False #TODO: This step can be computed from the previous rather than from scratch as done here
        else:
            if self.BUILD_COLIBRI_CLASSIFIERS:
                self.log("WARNING: No classifiers will be built! Only possible with Giza extraction method!",red)
            if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -s ' + self.getsourcefilename(alignsource) + ' -t ' + self.gettargetfilename(alignsource) + ' -o ' + self.gets2tfilename('alignmodel.colibri') + ' -I 2 ' + self.COLIBRI_ALIGNER_OPTIONS, "Building alignment model",self.gets2tfilename('alignmodel.colibri') ): return False
            #if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -t ' + self.getsourcefilename('graphpatternmodel.colibri') + ' -s ' + self.gettargetfilename('graphpatternmodel.colibri') + ' -o ' + self.gett2sfilename('alignmodel.colibri') + ' -N ' + self.COLIBRI_ALIGNER_OPTIONS, "Building reverse alignment model",self.gett2sfilename('alignmodel.colibri') ): return False

        if self.BUILD_COLIBRI_SKIPGRAMS:
            if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -U -d ' + self.gets2tfilename('alignmodel.colibri') + ' -s ' + self.getsourcefilename(alignsource) + ' -t ' + self.gettargetfilename(alignsource) + ' -o ' + self.gets2tfilename('alignmodelS.colibri') + ' ' +  self.COLIBRI_ALIGNER_OPTIONS , "Extracting skipgrams",self.gets2tfilename('alignmodelS.colibri') ): return False
            #if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -U -d ' + self.gett2sfilename('alignmodel.colibri') + ' -t ' + self.getsourcefilename('graphpatternmodel.colibri') + ' -s ' + self.gettargetfilename('graphpatternmodel.colibri') + ' -o ' + self.gett2sfilename('alignmodelS.colibri') + ' ' +  self.COLIBRI_ALIGNER_OPTIONS, "Extracting skipgrams (for reverse model)",self.gett2sfilename('alignmodelS.colibri') ): return False

        if self.BUILD_COLIBRI_MOSESPHRASETABLE:
            if self.BUILD_COLIBRI_SKIPGRAMS:
                if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -d ' + self.gets2tfilename('alignmodelS.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' --moses > ' + self.gets2tfilename('phrasetable'), "Outputting moses-style phrasetable from colibri",   self.gets2tfilename('phrasetable') ): return False
            else:
                if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -d ' + self.gets2tfilename('alignmodel.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' --moses > ' + self.gets2tfilename('phrasetable'), "Outputting moses-style phrasetable from colibri",   self.gets2tfilename('phrasetable') ): return False


        if self.BUILD_COLIBRI_CLASSIFIERS and self.BUILD_COLIBRI_GIZA:
            if not self.runcmd(self.EXEC_COLIBRI_TRAINCLASSIFIERS + ' ' + self.COLIBRI_CLASSIFIER_OPTIONS + ' -C timbl -d ' + self.gets2tfilename('alignmodel.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls')  + ' -O "' + self.COLIBRI_TIMBL_OPTIONS + '"', "Building and training classifiers"): return False


        #if self.BUILD_COLIBRI_TRANSTABLE:
        #    extra = ''
        #    if self.BUILD_COLIBRI_MOSESPHRASETABLE:
        #        extra = '-S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' --moses'
        #    else:
        #        extra = ''
        #    if self.BUILD_COLIBRI_SKIPGRAMS:
        #        if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -d ' + self.gets2tfilename('alignmodelS.colibri') + ' -i ' + self.gett2sfilename('alignmodelS.colibri') + ' -o ' + self.gets2tfilename('translationtable.colibri') + ' -I 2 ' + extra + ' > ' + self.gets2tfilename('phrasetable') , "Building Translation Table", self.gets2tfilename('translationtable.colibri') ): return False
        #    else:
        #        if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -d ' + self.gets2tfilename('alignmodel.colibri') + ' -i ' + self.gett2sfilename('alignmodel.colibri') + ' -o ' + self.gets2tfilename('translationtable.colibri') + ' -I 2 ' + extra + ' > ' + self.gets2tfilename('phrasetable') , "Building Translation Table", self.gets2tfilename('translationtable.colibri') ): return False

        return True


    def build_srilm_targetmodel(self):
        if not self.runcmd(self.EXEC_SRILM + ' -order ' + str(self.SRILM_ORDER) + ' ' + self.SRILM_OPTIONS + ' -text ' + self.gettargetfilename('txt') + ' -lm ' + self.gettargetfilename('srilm'),'SRILM Target-language Model', self.gettargetfilename('srilm')): return False
        return True

    def build_srilm_sourcemodel(self):
        if not self.runcmd(self.EXEC_SRILM +' -order ' + str(self.SRILM_ORDER) + ' ' + self.SRILM_OPTIONS + ' -text ' + self.getsourcefilename('txt') + ' -lm ' + self.getsourcefilename('srilm'),'SRILM Source-language Model', self.getsourcefilename('srilm')): return False
        return True


    def trainlimit(self):
        if not os.path.exists(self.getsourcefilename(str(self.TRAINLIMIT) + '.txt')):
            if not self.runcmd('head -n ' + str(self.TRAINLIMIT) + ' ' + self.getsourcefilename('txt') + ' > ' + self.getsourcefilename(str(self.TRAINLIMIT) + '.txt'), 'Downsizing training corpus (source)', self.getsourcefilename(str(self.TRAINLIMIT) + '.txt')): return False
            os.unlink(self.getsourcefilename('txt'))
            os.symlink( self.getsourcefilename(str(self.TRAINLIMIT) + '.txt'), self.getsourcefilename('txt') )
        if not os.path.exists( self.gettargetfilename(str(self.TRAINLIMIT) + '.txt')):
            if not self.runcmd('head -n ' + str(self.TRAINLIMIT) + ' ' + self.gettargetfilename('txt') + ' > ' + self.gettargetfilename(str(self.TRAINLIMIT) + '.txt'), 'Downsizing training corpus (target)', self.gettargetfilename(str(self.TRAINLIMIT) + '.txt')): return False
            os.unlink(self.gettargetfilename('txt'))
            os.symlink( self.gettargetfilename(str(self.TRAINLIMIT) + '.txt'), self.gettargetfilename('txt') )
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





    def build_moses(self):

        steps = []
        if self.BUILD_MOSES_PHRASETRANSTABLE:
            steps.append(8)
        if self.BUILD_MOSES:
            steps.append(9)

        firststep = 1
        laststep = max(steps)

        try:
            os.symlink(self.getsourcefilename('txt') ,  "train." + self.SOURCELANG)
        except:
            pass
        try:
            os.symlink(self.gettargetfilename('txt') ,  "train." + self.TARGETLANG)
        except:
            pass


        if not os.path.exists("model/phrase-table"):
            if not self.runcmd(self.EXEC_MOSES_TRAINMODEL + ' -external-bin-dir ' + self.PATH_MOSES_EXTERNALBIN + " -root-dir . --corpus train --f " + self.SOURCELANG + " --e " + self.TARGETLANG + " --first-step " + str(firststep) + " --last-step " + str(laststep) + " --lm 0:3:" + self.gettargetfilename('srilm')  + ' >&2 2> train-model.log',"Training model (moses) (logged in train-model.log)", "model/phrase-table.gz", "model/moses.ini"): return False
            os.system("gunzip -f model/phrase-table.gz")
            os.system("sed -i s/phrase-table\.gz/phrase-table/ model/moses.ini")
        else:
            print >>sys.stderr,bold(yellow("Skipping training model (moses), phrasetable already exists"))



        try:
            os.symlink("giza." + self.SOURCELANG + "-" + self.TARGETLANG + "/" + self.SOURCELANG + "-" + self.TARGETLANG + ".A3.final" ,self.gets2tfilename('A3.final'))
        except:
            pass

        try:
            os.symlink("giza." + self.TARGETLANG + "-" + self.SOURCELANG + "/" + self.TARGETLANG + "-" + self.SOURCELANG + ".A3.final" ,self.gett2sfilename('A3.final'))
        except:
            pass


        return True

    def build_moses_classifiers(self):
        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.getsourcefilename('txt'), "Encoding source corpus for Colibri",self.getsourcefilename('cls'), self.getsourcefilename('clsenc') ): return False

        if not self.runcmd(self.EXEC_COLIBRI_CLASSENCODE + ' -f ' + self.gettargetfilename('txt'), "Encoding target corpus for Colibri",self.gettargetfilename('cls'), self.gettargetfilename('clsenc') ): return False

        if self.COLIBRI_GLOBALKEYWORDS:
            if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.getsourcefilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS , "Building source-side pattern model",self.getsourcefilename('indexedpatternmodel.colibri') ): return False
            if not self.runcmd(self.EXEC_COLIBRI_PATTERNFINDER + ' -f ' + self.gettargetfilename('clsenc') + ' ' + self.COLIBRI_PATTERNFINDER_OPTIONS , "Building target-side pattern model",self.gettargetfilename('indexedpatternmodel.colibri') ): return False


            if not self.runcmd(self.EXEC_COLIBRI_ALIGNER + ' -m model/phrase-table' + ' -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') +' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' -k -K ' + str(self.COLIBRI_GLOBALKEYWORDS_OPTIONS) + ' -o ' + self.gets2tfilename('withkeywords.alignmodel.colibri') + ' -s ' + self.getsourcefilename('indexedpatternmodel.colibri') + ' -t ' + self.gettargetfilename('indexedpatternmodel.colibri') + ' 2> contextmoses-globalkeywords.log', "Extracting global keywords (logged in contextmoses-globalkeywords.log)",self.gets2tfilename("withkeywords.alignmodel.colibri")): return False
            if not ('-I' in self.MOSES_CLASSIFIER_OPTIONS):
                if not self.runcmd(self.EXEC_COLIBRI_CONTEXTMOSES + ' -f ' + self.getsourcefilename('clsenc') + ' -g ' + self.gettargetfilename('clsenc') + ' -d ' +  self.gets2tfilename('withkeywords.alignmodel.colibri') + ' -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' -s ' + self.getsourcefilename('indexedpatternmodel.colibri') + ' -E ' + self.gettargetfilename('indexedpatternmodel.colibri') + ' -k ' + self.MOSES_CLASSIFIER_OPTIONS + ' 2> contextmoses-train.log', "Training classifiers for context-aware moses (logged in contextmoses-train.log)"): return False
            else:
                print >>sys.stderr, bold(yellow("Not training classifiers because -I (ignore) is set in options"))
        else:
            if not ('-I' in self.MOSES_CLASSIFIER_OPTIONS):
                if not self.runcmd(self.EXEC_COLIBRI_CONTEXTMOSES + ' -f ' + self.getsourcefilename('clsenc') + ' -g ' + self.gettargetfilename('clsenc') + ' -m ' +  'model/phrase-table' + ' -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' ' + self.MOSES_CLASSIFIER_OPTIONS + ' 2> contextmoses-train.log', "Training classifiers for context-aware moses (logged in contextmoses-train.log)"): return False
            else:
                print >>sys.stderr, bold(yellow("Not training classifiers because -I (ignore) is set in options"))



        return True


    def build_moses_mert(self):
        if self.BUILD_MOSES_CLASSIFIERS:

            if not self.runcmd(self.EXEC_COLIBRI_CONTEXTMOSES + ' -F ' + self.DEVSOURCECORPUS + ' -m model/phrase-table -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' ' + self.MOSES_CLASSIFIER_OPTIONS + ' -q -o devtmp 2> contextmoses-dev.log', "Preparing context-aware moses on development set (logged in contextmoses-dev.log)"): return False

            if self.MOSES_MERT_RUNS == 1:
                if not self.runcmd(self.EXEC_MOSES_MERT + ' --mertdir=' + self.PATH_MOSES_MERT + ' ' + self.MOSES_MERT_OPTIONS + ' ' + self.WORKDIR + '/devtmp.txt ' + self.DEVTARGETCORPUS + ' ' + self.EXEC_MOSES  + ' ' + self.WORKDIR + '/model/contextmoses.devtmp.ini 2> mert.log', 'Parameter tuning for Moses (+context) using MERT (logged in mert.log)'): return False
                shutil.copyfile(self.WORKDIR + '/mert-work/moses.ini', self.WORKDIR + '/model/contextmoses.tmp.ini')
                os.system("sed -i s/devtmp\.phrasetable/tmp.phrasetable/ model/contextmoses.tmp.ini")
            else:
                for i in range(1, self.MOSES_MERT_RUNS+1):
                    if not self.runcmd(self.EXEC_MOSES_MERT + ' --mertdir=' + self.PATH_MOSES_MERT + ' --working-dir=mert-work' + str(i) + ' ' + self.MOSES_MERT_OPTIONS + ' ' + self.WORKDIR + '/devtmp.txt ' + self.DEVTARGETCORPUS + ' ' + self.EXEC_MOSES  + ' ' + self.WORKDIR + '/model/contextmoses.devtmp.ini 2> mert.log', 'Parameter tuning for Moses (+context) using MERT (logged in mert.log)'): return False

                shutil.copyfile(self.WORKDIR + '/mert-work' + str(i) + '/moses.ini', self.WORKDIR + '/mert-work' + str(i) + '/contextmoses.tmp.ini')
                os.system("sed -i s/devtmp\.phrasetable/tmp.phrasetable/ mert-work" + str(i) + "/contextmoses.tmp.ini")


        else:
            if self.MOSES_MERT_RUNS == 1:
                if not self.runcmd(self.EXEC_MOSES_MERT + ' --mertdir=' + self.PATH_MOSES_MERT + ' ' + self.MOSES_MERT_OPTIONS + ' ' + self.DEVSOURCECORPUS + ' ' + self.DEVTARGETCORPUS + ' ' + self.EXEC_MOSES  + ' ' + self.WORKDIR + '/model/moses.ini 2> mert.log', 'Parameter tuning for Moses using MERT (logged in mert.log)'): return False
            else:
                for i in range(1, self.MOSES_MERT_RUNS+1):
                    if not self.runcmd(self.EXEC_MOSES_MERT + ' --mertdir=' + self.PATH_MOSES_MERT + ' --working-dir=' + self.WORKDIR + '/mert-work' + str(i) + ' ' + self.MOSES_MERT_OPTIONS + ' ' + self.DEVSOURCECORPUS + ' ' + self.DEVTARGETCORPUS + ' ' + self.EXEC_MOSES  + ' ' + self.WORKDIR + '/model/moses.ini 2> mert.log', 'Parameter tuning for Moses using MERT (logged in mert.log)'): return False
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


    def get_phrasal_classpath(self):
        classpath = []
        for filename in glob.glob(self.PATH_PHRASAL+'/*.jar'): classpath.append(filename)
        for filename in glob.glob(self.PATH_CORENLP+'/*.jar'): classpath.append(filename)
        if self.JAR_FASTUTIL: classpath.append(self.JAR_FASTUTIL)
        return ":".join(classpath)

    def build_phrasal_wordalign(self):
        if os.path.exists(self.gets2tfilename('A3.final')) and os.path.exists(self.gett2sfilename('A3.final')):
            self.log("Skipping Berkeley word aligner (output already exists)",yellow, True)
            return True


        try:
            os.mkdir(self.WORKDIR + '/alignerinput')
        except:
            pass

        try:
            os.symlink(self.getsourcefilename('txt'), self.WORKDIR + '/alignerinput/' + self.CORPUSNAME + '.' + self.SOURCELANG)
        except:
            pass

        try:
            os.symlink(self.gettargetfilename('txt'), self.WORKDIR + '/alignerinput/' + self.CORPUSNAME + '.' + self.TARGETLANG)
        except:
            pass



        f = open(self.WORKDIR + '/aligner.conf','w')
        f.write("""
# aligner.conf
##########################################
# Training: Defines the training regimen #
##########################################

forwardModels HMM
reverseModels HMM
mode JOINT
iters 2

###############################################
# Execution: Controls output and program flow #
###############################################

execDir aligneroutput
create
saveParams true
numThreads 4
msPerLine 10000
alignTraining
leaveTrainingOnDisk
safeConcurrency true
overwrite

#################
# Language/Data #
#################

foreignSuffix %s
englishSuffix %s
lowercase

# Choose the training sources, which can either be directories or files that list files/directories
trainSources alignerinput
sentences MAX

# The test sources must have hand alignments for all sentence pairs
testSources
maxTestSentences MAX
offsetTestSentences 0

##############
# Evaluation #
##############

competitiveThresholding
writeGIZA""" % (self.TARGETLANG, self.SOURCELANG) )
        f.close()
        classpath = self.get_phrasal_classpath()
        JAVA_OPTS="-XX:+UseCompressedOops -Xmx" + str(self.PHRASAL_MAXMEM) + ' -Xms' +  str(self.PHRASAL_MAXMEM)
        if self.runcmd('CLASSPATH=' + classpath + '  ' +  self.EXEC_JAVA + ' ' + JAVA_OPTS + ' -jar ' + self.JAR_BERKELEYALIGNER + ' ++aligner.conf', 'Word alignment'):
            shutil.copyfile(self.WORKDIR + '/aligneroutput/training.' + self.SOURCELANG + '-' + self.TARGETLANG + '.A3', self.gets2tfilename('A3.final'))
            shutil.copyfile(self.WORKDIR + '/aligneroutput/training.' + self.TARGETLANG + '-' + self.SOURCELANG + '.A3',  self.gett2sfilename('A3.final'))
        else:
            return False
        return True

    def build_phrasal_phraseextract(self):
        if os.path.exists(self.gets2tfilename('phrasetable')):
            self.log("Skipping Phrasal Phrase Extraction (output already exists)",yellow, True)
        else:
            classpath = self.get_phrasal_classpath()
            JAVA_OPTS="-XX:+UseCompressedOops -Xmx" + str(self.PHRASAL_MAXMEM) + ' -Xms' +  str(self.PHRASAL_MAXMEM)
            #EXTRACT_OPTS="-inputDir aligneroutput -outputFile phrases"
            EXTRACT_OPTS="-fCorpus " + self.getsourcefilename('txt') + ' -eCorpus ' + self.gettargetfilename('txt') + ' -feAlign ' + self.gets2tfilename('A3.final') + ' -efAlign ' + self.gett2sfilename('A3.final') + " -outputFile " + self.gets2tfilename('phrasetable') + ' ' + self.PHRASAL_PHRASEEXTRACT_OPTIONS
            cmd = 'CLASSPATH=' + classpath + ' ' + self.EXEC_JAVA + ' ' + JAVA_OPTS + ' edu.stanford.nlp.mt.train.PhraseExtract ' + EXTRACT_OPTS
            if self.DEVSOURCECORPUS: #should always exist, has been checked in checking stage
                cmd += ' -fFilterCorpus ' + self.DEVSOURCECORPUS
            if self.PHRASAL_WITHGAPS:
                cmd += ' -withGaps true'
            if not self.runcmd(cmd , 'Phrase extraction',  self.gets2tfilename('phrasetable') ): return False
        if not self.runcmd(self.EXEC_PERL + ' ' + self.PATH_PHRASAL + '/scripts/split-table ' + self.WORKDIR + '/phrases-tm.gz ' + self.WORKDIR + '/phrases-om.gz < ' + self.gets2tfilename('phrasetable') , 'Phrase splitting',  self.WORKDIR + '/phrases-tm.gz',  self.WORKDIR + '/phrases-om.gz'): return False
        return True

    def build_phrasal(self):
        self.log("Writing phrasal configuration",green, True)
        f = open(self.WORKDIR + '/phrasal.conf','w')
        f.write("""# filename: phrasal.conf

# translation table
[ttable-file]
phrases-tm.gz

# language model
[lmodel-file]
%s

# number of translation options for each phrase in f
[ttable-limit]
20

[additional-featurizers]
edu.stanford.nlp.mt.decoder.feat.HierarchicalReorderingFeaturizer(phrases-om.gz,msd2-bidirectional-fe,LexR,hierarchical,hierarchical,bin)

# maximum gap between covered spans
[distortion-limit]
6

[weights-file]
phrasal.wts

# detect processors present, and use them all
[localprocs]
0\n\n""" % (self.gettargetfilename('srilm'),) )
        if self.PHRASAL_WITHGAPS:
            f.write("[gaps]\n")
            f.write(str(self.PHRASAL_MAXSOURCEPHRASESPAN) + "\n")
            f.write(str(self.PHRASAL_MAXTARGETPHRASESPAN) + "\n")
        f.close()

        f = open(self.WORKDIR + '/phrasal.wts','w')
        f.write("""LM 2
LinearDistortion 1
TM:FPT.0 0.5
TM:FPT.1 0.5
TM:FPT.2 0.5
TM:FPT.3 0.5
TM:FPT.4 -1
WordPenalty: -0.5\n""")
        f.close()

        return True

    def build_phrasal_mert(self):
        cmd = self.EXEC_PERL + ' ' + self.EXEC_PHRASAL_MERT + ' ' + self.DEVSOURCECORPUS + ' ' + self.DEVTARGETCORPUS
        if self.PHRASAL_WITHGAPS:
            cmd += " --phrasal_flags=\"--gaps " + str(self.PHRASAL_MAXSOURCEPHRASESPAN) + ' ' + str(self.PHRASAL_MAXTARGETPHRASESPAN) + "\""
        if not self.runcmd(cmd):
            return False
        else:
            #TODO: copy MERT output
            pass
        return True



    def server(self, port, html=False):
        if not self.check_common(): return False
        if not self.check_run(): return False

        if self.BUILD_MOSES: self.server_moses(port, html)



    def run(self, inputfile, outputfile='output.txt', tokenise=False):
        if tokenise and (not self.EXEC_UCTO or not os.path.isfile(self.EXEC_UCTO)):
            self.log("Error: Ucto not found! Unable to tokenise!" ,red)
            return False

        if not self.check_common(): return False
        if not self.check_run(): return False

        if not os.path.isfile(inputfile):
            self.log("Error: Input file " + inputfile + " not found!" ,red)
            return False

        if os.path.exists( outputfile):
            os.unlink( outputfile)

        if tokenise:
            if not self.runcmd(self.EXEC_UCTO + ' -n -L' + self.SOURCELANG +  ' ' + inputfile + ' ' + 'input.txt','Tokenisation of Input File'): return False
        else:
            if os.path.exists( self.WORKDIR + '/input.txt'):
                os.unlink( self.WORKDIR + '/input.txt' )
            os.symlink(inputfile, self.WORKDIR + '/input.txt' )

        if self.BUILD_MOSES_CLASSIFIERS and not self.run_moses_classifiers(): return False
        if self.BUILD_MOSES and not self.BUILD_MOSES_CLASSIFIERS and not self.run_moses(): return False
        if self.BUILD_PBMBMT and not self.run_pbmbmt(): return False
        if self.BUILD_PHRASAL and not self.run_phrasal(): return False
        if self.BUILD_COLIBRI and not self.run_colibri(): return False


        os.rename('output.txt',outputfile)
        return True


    def mert_computeaveragescore(self):
        bleu = []
        meteor = []
        nist = []
        ter = []
        wer = []
        per = []
        for i in range(1,self.MOSES_MERT_RUNS+1):
            f = open('summary-mert' + str(i) +'.score')
            f.readline()
            scores = ( int(x) for x in f.readline().split() )
            for i, x in enumerate((bleu,meteor,nist,ter,wer,per)):
                x.append(scores[i])
            f.close()
        f = open('summary.score','w')
        f.write("BLEU METEOR NIST TER WER PER\n")
        f.write(str(sum(bleu)/float(len(bleu))) + " ")
        f.write(str(sum(meteor)/float(len(meteor))) + " ")
        f.write(str(sum(nist)/float(len(nist))) + " ")
        f.write(str(sum(ter)/float(len(ter))) + " ")
        f.write(str(sum(wer)/float(len(wer))) + " ")
        f.write(str(sum(per)/float(len(per))) + " ")
        f.close()

    def run_moses(self):
        if self.BUILD_MOSES_MERT and self.MOSES_MERT_RUNS > 1:
            for i in range(1,self.MOSES_MERT_RUNS+1):
                mosesini = self.WORKDIR + '/mert-work' + str(i) + '/moses.ini'
                if not self.runcmd(self.EXEC_MOSES + ' -f ' + mosesini + ' < input.txt > output.txt 2> moses.log','Moses Decoder (logged in moses.log)'): return False
                os.rename('summary.score','summary-mert'+str(i)+'.score')
            self.mert_computeaveragescore()
            return True
        else:
            if self.BUILD_MOSES_MERT:
                mosesini = self.WORKDIR + '/mert-work/moses.ini'
            else:
                mosesini = self.WORKDIR + '/model/moses.ini'
            if not self.runcmd(self.EXEC_MOSES + ' -f ' + mosesini + ' < input.txt > output.txt 2> moses.log','Moses Decoder (logged in moses.log)'): return False
            return True

    def run_moses_classifiers(self):
        #should work with MERT as well
        if not os.path.exists('tmp.srilm'):
            os.symlink(self.gettargetfilename('srilm'), 'tmp.srilm')
        if self.COLIBRI_GLOBALKEYWORDS:
            extra = ""
            if self.COLIBRI_GLOBALKEYWORDS_OPTIONS:
                extra = " -K " + self.COLIBRI_GLOBALKEYWORDS_OPTIONS.split(',')[-1] #keyword probability
            if not self.runcmd(self.EXEC_COLIBRI_CONTEXTMOSES + ' -k ' + extra + ' -F input.txt -d ' + self.gets2tfilename('withkeywords.alignmodel.colibri') + ' -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' ' + self.MOSES_CLASSIFIER_OPTIONS + ' > output.txt 2> contextmoses-test.log', "Testing classifiers and running context-aware moses decoder (logged in contextmoses-test.log)"): return False
        else:
            if not self.runcmd(self.EXEC_COLIBRI_CONTEXTMOSES + ' -F input.txt -m model/phrase-table -S ' +  self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' -l ' + str(self.MOSES_LEFTCONTEXTSIZE) + ' -r ' + str(self.MOSES_RIGHTCONTEXTSIZE) + ' ' + self.MOSES_CLASSIFIER_OPTIONS + ' > output.txt 2> contextmoses-test.log', "Testing classifiers and running context-aware moses decoder (logged in contextmoses-test.log)"): return False
        return True

    def run_phrasal(self):
        classpath = self.get_phrasal_classpath()
        JAVA_OPTS="-XX:+UseCompressedOops -Xmx" + str(self.PHRASAL_MAXMEM) + ' -Xms' +  str(self.PHRASAL_MAXMEM)
        cmd = 'CLASSPATH=' + classpath + ' ' + self.EXEC_JAVA + ' ' + JAVA_OPTS + ' edu.stanford.nlp.mt.Phrasal -config-file phrasal.conf < input.txt > output.txt'
        if not self.runcmd(cmd,"Phrasal Decoder"):
            return False
        return True

    def run_colibri(self):
        decoder_extraoptions = ''
        if self.BUILD_COLIBRI_CLASSIFIERS:
            decoder_extraoptions = '-C timbl'
            if self.COLIBRI_TIMBL_OPTIONS:
                decoder_extraoptions += ' -O "' + self.COLIBRI_TIMBL_OPTIONS + '"'

        if os.path.exists(self.gets2tfilename('transtable.colibri')): #backward compatibility
            if not self.runcmd(self.EXEC_COLIBRI_DECODER + ' -l ' + self.gettargetfilename('srilm') + ' -t ' + self.gets2tfilename('transtable.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' ' + self.COLIBRI_DECODER_OPTIONS +  ' ' + decoder_extraoptions + '  < input.txt > output.txt 2> decoder.log','Colibri Decoder (logged in decoder.log)'): return False
        elif os.path.exists(self.gets2tfilename('alignmodelS.colibri')):
            if not self.runcmd(self.EXEC_COLIBRI_DECODER + ' -l ' + self.gettargetfilename('srilm') + ' -t ' + self.gets2tfilename('alignmodelS.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' ' + self.COLIBRI_DECODER_OPTIONS +  ' ' + decoder_extraoptions + ' < input.txt > output.txt 2> decoder.log','Colibri Decoder (logged in decoder.log)'): return False
        elif os.path.exists(self.gets2tfilename('alignmodel.colibri')):
            if not self.runcmd(self.EXEC_COLIBRI_DECODER + ' -l ' + self.gettargetfilename('srilm') + ' -t ' + self.gets2tfilename('alignmodel.colibri') + ' -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' ' + self.COLIBRI_DECODER_OPTIONS + ' ' + decoder_extraoptions + ' < input.txt > output.txt 2> decoder.log','Colibri Decoder (logged in decoder.log)'): return False
        elif os.path.exists(self.gets2tfilename('phrasetable')):
            if self.BUILD_COLIBRI_CLASSIFIERS:
                self.log("WARNING: No classifiers will be used! Not possible with Moses phrasetable!",red)
            #moses-style phrase-table
            if not self.runcmd(self.EXEC_COLIBRI_DECODER + ' -l ' + self.gettargetfilename('srilm') + ' -t ' + self.gets2tfilename('phrasetable') + ' --moses -S ' + self.getsourcefilename('cls') + ' -T ' + self.gettargetfilename('cls') + ' ' + self.COLIBRI_DECODER_OPTIONS + ' < input.txt > output.txt 2> decoder.log','Colibri Decoder (logged in decoder.log)'): return False
        else:
            self.log("Error: No phrasetable found! Did you forget to train the system?" ,red)
            return False
        return True


    def server_moses(self, port, html):
        if self.BUILD_MOSES_MERT:
            mosesini = self.WORKDIR + '/mert-work/moses.ini'
        else:
            mosesini = self.WORKDIR + '/model/moses.ini'
        while True:
            if not html:
                GenericWrapperServer(self.EXEC_MOSES + ' -f ' + mosesini + ' 2> ' + self.WORKDIR + '/moses.log', port, True, False) #print stdout, not send stderr
            else:
                GenericWrapperServer(self.EXEC_MOSES + ' -f ' + mosesini + ' 2> ' + self.WORKDIR + '/moses.log', port, False, True, lambda x: None, serveroutputproc) #print stderr, not send stdout, but filter first
            print >>sys.stderr, "Server process failed? Restarting..."
            #server down? restart
            time.sleep(10)


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
            self.log("Error: Output file " + targetfile + " not found! Did you forget to test the system?" ,red)
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
                    self.log("Error reading bleu.score",red)
                    errors = True
        else:
            self.log("Skipping BLEU (no script found ["+self.EXEC_MATREX_BLEU+"])",yellow)

        if self.EXEC_MATREX_WER and os.path.exists(self.EXEC_MATREX_WER) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_WER + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml + '  > ' + 'wer.score', 'Computing WER score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/wer.score')
                    for line in f:
                        if line[0:11] == "WER score =":
                            self.wer = float(line[12:20].strip())
                            self.log("WER score: " + str(self.wer), white)
                    f.close()
                except:
                    self.log("Error reading wer.score",red)
                    errors = True
        else:
            self.log("Skipping WER (no script found ["+self.EXEC_MATREX_WER+"]) ",yellow)

        if self.EXEC_MATREX_PER and os.path.exists(self.EXEC_MATREX_PER) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_PER + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml + '  > ' + 'per.score',  'Computing PER score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/per.score')
                    for line in f:
                        if line[0:11] == "PER score =":
                            self.per = float(line[12:20].strip())
                            self.log("PER score: " + str(self.per), white)
                    f.close()
                except:
                    self.log("Error reading per.score",red)
                    errors = True
        else:
            self.log("Skipping PER (no script found ["+self.EXEC_MATREX_PER+"])",yellow)

        if self.EXEC_MATREX_METEOR and os.path.exists(self.EXEC_MATREX_METEOR) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' -I ' + os.path.dirname(self.EXEC_MATREX_METEOR) + ' ' + self.EXEC_MATREX_METEOR + " -s " + self.CORPUSNAME + " -r " + refxml + ' -t ' + targetxml + ' --modules "exact"  > ' + 'meteor.score',  'Computing METEOR score'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/meteor.score')
                    for line in f:
                        if line[0:6] == "Score:":
                            self.meteor = float(line[7:].strip())
                            self.log("METEOR score: " + str(self.meteor), white)
                    f.close()
                except:
                    self.log("Error reading meteor.score",red)
                    errors = True
        else:
            self.log("Skipping METEOR (no script found ["+self.EXEC_MATREX_METEOR+"])",yellow)

        if self.EXEC_MATREX_MTEVAL and os.path.exists(self.EXEC_MATREX_MTEVAL) and self.EXEC_PERL and os.path.exists(self.EXEC_PERL):
            if not self.runcmd(self.EXEC_PERL + ' ' + self.EXEC_MATREX_MTEVAL + " -r " + refxml + ' -t ' + targetxml + ' -s ' + sourcexml +  '  > ' + 'mteval.score',  'Computing NIST & BLEU scores'): errors = True
            if not errors:
                try:
                    f = open(self.WORKDIR + '/mteval.score')
                    for line in f:
                        if line[0:12] == "NIST score =":
                            self.nist = float(line[13:21].strip())
                            print >>sys.stderr,"NIST score: ", self.nist
                        if line[21:33] == "BLEU score =":
                            try:
                                bleu2 = float(line[34:40].strip())
                                if self.bleu == 0:
                                    self.bleu = bleu2
                                    self.log("BLEU score: " + str(self.bleu), white)
                                elif abs(self.bleu - bleu2) > 0.01:
                                    self.log("blue score from MTEVAL scripts differs too much: " + str(self.bleu) + " vs " + str(bleu2) +  ", choosing highest score")
                                    if bleu2 > self.bleu:
                                        self.bleu = bleu2
                                else:
                                    self.log("BLEU score (not stored): " + str(float(line[34:40].strip())))
                            except:
                                raise
                    f.close()
                except:
                    self.log("Error reading mteval.score",red)
                    errors = True
        else:
            self.log("Skipping MTEVAL (BLEU & NIST) (no script found)", yellow)

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
                    self.log("Error reading ter.score",red)
        else:
            self.log("Skipping TER (no script found)",yellow)


        print >>sys.stderr,"SCORE SUMMARY\n===================\n"
        f = open(self.WORKDIR + '/summary.score','w')
        s = "BLEU METEOR NIST TER WER PER"
        f.write(s+ "\n")
        print >>sys.stderr, s
        s = str(round(self.bleu,4)) + " " + str(round(self.meteor,4)) + " " + str(round(self.nist,4))  + " " + str(round(self.ter,2)) + " " + str(round(self.wer,2))  + " " + str(round(self.per,2))
        f.write(s + "\n")
        print >>sys.stderr, s
        f.close()


        return not errors

    def test(self, sourcefile, reffile):
        if not self.check_common(): return False
        if not self.check_run(): return False
        if not self.check_test(): return False

        if not os.path.isfile(sourcefile):
            self.log("Error: Source file " + sourcefile + " not found!" ,red)
            return False

        if not os.path.isfile(reffile):
            self.log("Error: Reference file " + reffile + " not found!" ,red)
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


    def downsize(self, sourcefile, targetfile, lines):
        fin = open(sourcefile,'r')
        fout = open(targetfile,'w')
        for i, line in enumerate(fin):
            if i < lines:
                fout.write(fout)
        fin.close()
        fout.close()
        self.log("Branched file " + sourcefile + " (downsized to " + str(lines)+")",green)

    def copychangepaths(self,sourcefile, targetfile, sourcepath, targetpath):
        f_in = codecs.open(sourcefile,'r','utf-8')
        f_out = codecs.open(targetfile,'w','utf-8')
        for line in f_in:
            line = line.replace(sourcepath, targetpath)
            f_out.write(line)
        f_in.close()
        f_out.close()


    def branch(self,expname, conf=None, useparentdir=True, quiet = False, writebatches=True):
        parentdir = self.WORKDIR
        if useparentdir:
            if parentdir[-1] == '/':
                parentdir = parentdir[:-1]
            parentdir = os.path.dirname(parentdir)

        if conf:
            conf = self.parseconf(conf)

        workdir = parentdir + '/' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + expname
        if workdir and not os.path.isdir(workdir):
            self.log("Creating branched work directory (as sibling): " + workdir,white)
            os.mkdir(workdir)
            for filename in glob.glob(self.WORKDIR + '/*'):
               basefilename = os.path.basename(filename)

               if basefilename[-4:] == '.ini' or basefilename[-4:] == '.cfg' or basefilename[-5:] == '.conf':
                    self.copychangepaths(filename, workdir + '/' + basefilename, self.WORKDIR, workdir)
               elif (basefilename == 'model' or basefilename == 'mert-work' or basefilename == 'corpus') and os.path.isdir(filename):
                    try:
                        os.mkdir(workdir + '/' + basefilename)
                    except:
                        pass
                    for filename2 in glob.glob(self.WORKDIR + '/' + basefilename + '/*'):
                        basefilename2 = os.path.basename(filename2)
                        if basefilename2[-4:] == '.ini' or basefilename2[-4:] == '.cfg' or basefilename2[-5:] == '.conf':
                            self.copychangepaths(filename2, workdir + '/' + basefilename + '/' + basefilename2, self.WORKDIR, workdir)
                        else:
                            try:
                                os.symlink(filename2, workdir + '/' + basefilename + '/' + basefilename2)
                                self.log("Branched file " + basefilename + "/" + basefilename2 + " (symlink)",green)
                            except:
                                self.log("Error making symlink for " + basefilename + "/" + basefilename2,red)
               elif (basefilename[-3:] == '.py' and basefilename[0:3] == 'mt-') or basefilename[0] == '.' or os.path.isdir(filename) or basefilename[-4:] in ['.log','.tex','.png','.jpg','.aux','.pdf'] or basefilename[-1] == '~':
                    continue
                #elif 'TRAINSOURCECORPUS' in conf and basefilename == os.path.basename(conf['TRAINSOURCECORPUS']) and 'TRAINSIZE' in conf:
                #    self.downsize(filename, workdir + '/' + basefilename, int(conf['TRAINSIZE']))
                #elif 'TRAINTARGETCORPUS' in conf and basefilename == os.path.basename(conf['TRAINTARGETCORPUS']) and 'TRAINSIZE' in conf:
                #    self.downsize(filename, workdir + '/' + basefilename, int(conf['TRAINSIZE']))
               else:
                    try:
                        os.symlink(filename, workdir + '/' + basefilename)
                        self.log("Branched file " + basefilename + " (symlink)",green)
                    except:
                        self.log("Error making symlink for " + basefilename,red)
        elif workdir and not quiet:
            self.log("WARNING: work directory " +  workdir + " already exists! Press ENTER to continue or ctrl-C to abort",white)
            raw_input()



        settingsfile = self.writesettings(expname, workdir, writebatches)

        self.setargs(**self.confdata) #reset

        return (workdir, settingsfile)


    def writesettings(self,expname=None, workdir = None, writebatches=True):
        if not workdir:
            workdir = self.WORKDIR
        if not expname:
            expname = self.EXPERIMENTNAME
        if expname:
            settingsfile = workdir + '/mt-' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '-' + expname + '.py'
        else:
            settingsfile = workdir + '/mt-' + self.CORPUSNAME + '-' + self.SOURCELANG + '-' + self.TARGETLANG + '.py'


        f = codecs.open(settingsfile,'w','utf-8')
        f.write("#! /usr/bin/env python\n# -*- coding: utf8 -*-#\n\n")
        f.write("#Generated as branch of " + self.WORKDIR + "\n\n")
        f.write('import sys\n')
        #for dir in sys.path:
        #    if dir:
        #        f.write("sys.path.append('" + dir + "')\n")
        f.write("\n")
        f.write("from mtwrapper2 import MTWrapper\n")
        f.write("mtwrapper = MTWrapper(")
        skip = {}
        for source in self.sources:
            f.write("'" + source + "',")
            exec 'from ' + source + ' import confdata'
            for key in confdata:
                skip[key] = True
        f.write("\n")
        for key, default, help in MTWrapper.defaults:
            if not key in skip:
                if key == 'EXPERIMENTNAME':
                    value = expname
                elif key == 'WORKDIR':
                    value = workdir
                else:
                    value = self.__getattribute__(key)

                if key[:5] == 'PATH_' and (value == default or value == ""):
                    continue
                if key[:5] == 'EXEC_' and (value == default or value == ""):
                    continue

                if isinstance(value, str) or isinstance(value,  unicode):
                    if '"' in value: value = value.replace('"','\\' + '"')
                    f.write("    " + key + "=\"" + value + "\"")
                else:
                    f.write("    " + key + "=" + str(value))

                if help:
                    f.write(", #" + help + "\n")
                else:
                    f.write(",\n")
        f.write(")\n")
        if writebatches:
            for batch, conf in self.batches:
                f.write("mtwrapper.addbatch('" + batch + "', \n**"+repr(conf)+")\n\n")
        f.write("mtwrapper.start()\n")
        f.close()
        os.chmod(settingsfile, 0754)

        return settingsfile


def resource(sourcecorpusfile, targetcorpusfile, testset, devset, trainset, workdir, sourcelang, targetlang, corpusname):
        """re-source: draw new samples from sources"""
        filesampler([sourcecorpusfile, targetcorpusfile],  testset, devset, trainset, workdir )


        #rename files
        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.dev'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' + corpusname + '-' + sourcelang + '-dev.txt')
            devsourcecorpusfile = workdir + '/' + corpusname + '-' + sourcelang + '-dev.txt'
        else:
            devsourcecorpusfile = None

        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.dev'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/'+ corpusname + '-' + targetlang + '-dev.txt')
            devtargetcorpusfile = workdir + '/' +corpusname + '-' + targetlang + '-dev.txt'
        else:
            devtargetcorpusfile = None

        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.test'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' +corpusname + '-' + sourcelang + '-test.txt')
            testsourcecorpusfile = workdir + '/' +corpusname + '-' + sourcelang + '-test.txt'
        else:
            testsourcecorpusfile = None

        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.test'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' + corpusname + '-' + targetlang + '-test.txt')
            testtargetcorpusfile = workdir + '/' + corpusname + '-' + targetlang + '-test.txt'
        else:
            testtargetcorpusfile = None

        oldfile = workdir + '/' + os.path.basename(sourcecorpusfile) + '.train'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' + corpusname + '-' + sourcelang + '-train.txt')
            sourcecorpusfile = workdir + '/' + corpusname + '-' + sourcelang + '-train.txt'
        else:
            sourcecorpusfile = None

        oldfile = workdir + '/' + os.path.basename(targetcorpusfile) + '.train'
        if os.path.exists(oldfile):
            os.rename(oldfile, workdir + '/' + corpusname + '-' + targetlang + '-train.txt')
            targetcorpusfile = workdir + '/' + corpusname + '-' + targetlang + '-train.txt'
        else:
            targetcorpusfile = None

        return (sourcecorpusfile, targetcorpusfile, testsourcecorpusfile, testtargetcorpusfile, devsourcecorpusfile, devtargetcorpusfile)

def usage():
    print >>sys.stderr,"mtwrapper.py -- MT wrapper - Outputs a MT wrapper script (python)"
    print >>sys.stderr,"Mandatory Input:"
    print >>sys.stderr,"\t-n <name>         Name of the corpus [MANDATORY!]"
    print >>sys.stderr,"\t-s <file>         Corpus in source language (for training)"
    print >>sys.stderr,"\t-t <file>         Corpus in target language (for training)"
    print >>sys.stderr,"\t-S <code>         Source language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-T <code>         Target language (iso-639-1 or 3)"
    print >>sys.stderr,"\t-w <dir>          Work directory (by default: current dir)"
    print >>sys.stderr,"\t-I <module>       Include settings module"
    print >>sys.stderr,"Optional Input:"
    print >>sys.stderr,"\t--testset=n          Extract a random sample of n lines as test set, and exclude from training"
    print >>sys.stderr,"\t--devset=n           Extract a random sample of n lines as development set, and exclude from training"
    print >>sys.stderr,"\t--trainset=n         Restrict the training set to a random sample of n lines"
    print >>sys.stderr,"\t--srcdev=filename    Development corpus in source language (can not be used with --testset, --devset,--trainset)"
    print >>sys.stderr,"\t--tgtdev=filename    Development corpus in target language  (can not be used with --testset, --devset,--trainset)"
    print >>sys.stderr,"\t--srctst=filename    Test corpus in source language  (can not be used with --testset, --devset,--trainset)"
    print >>sys.stderr,"\t--tgttst=filename    Target corpus in target language  (can not be used with --testset, --devset,--trainset)"
    print >>sys.stderr,"\t-i <dirs>         Colon-separated directories where python can find non-standard modules"




if __name__ == "__main__":



    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:t:S:T:n:x:w:d:i:I:", ['testset=','devset=','trainset=','srcdev=','tgtdev=','srctst=','tgttst='])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    parentdir = None
    includemod = None
    sourcecorpusfile = targetcorpusfile = sourcelang = targetlang = corpusname = expname = workdir = ""
    devsourcecorpusfile = devtargetcorpusfile = testsourcecorpusfile = testtargetcorpusfile = ""
    srcdev = tgtdev = srctst = tgttst = None
    trainset = testset = devset = 0
    includedirs = []
    confdata = {}

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
        elif o == '-I':
            includemod = a
            try:
                exec 'from ' + includemod + ' import confdata'
            except:
                print >>sys.stderr,"Unable to import " + includemod
                sys.exit(2)
        elif o == '-i':
            includedirs = a.split(':')
        elif o == '--testset':
            testset = int(a)
        elif o == '--devset':
            devset = int(a)
        elif o == '--trainset':
            trainset = int(a)
        elif o == '--srcdev':
            devsourcecorpusfile = a
        elif o == '--tgtdev':
            devtargetcorpusfile = a
        elif o == '--srctst':
            testsourcecorpusfile = a
        elif o == '--tgttst':
            testtargetcorpusfile = a
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
        workdir = os.getcwd() + '/' + workdir
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
            sys.exit(2)

        sourcecorpusfile, targetcorpusfile, testsourcecorpusfile, testtargetcorpusfile, devsourcecorpusfile, devtargetcorpusfile = resource( sourcecorpusfile, targetcorpusfile, testset, devset, trainset, workdir, sourcelang, targetlang, corpusname)



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
    f.write("from mtwrapper2 import MTWrapper\n")
    f.write("mtwrapper = MTWrapper(")
    if includemod:
        f.write("'" + includemod + "',")
    f.write("\n")
    for key, default, help in MTWrapper.defaults:
        if not key in confdata:
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
    if 'EDITOR' in os.environ:
        editor = os.environ['EDITOR']
    else:
        editor = 'vim'
    os.system(editor + ' ' + settingsfile)




