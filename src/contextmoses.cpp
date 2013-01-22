
#include <getopt.h>
#include "alignmodel.h"
#include "classifiers.h"

using namespace std;

int addunknownwords(ClassEncoder * sourceclassencoder, ClassDecoder * sourceclassdecoder,  ClassEncoder * targetclassencoder, ClassDecoder * targetclassdecoder) {
    //Sync sourceclassdecoder with sourceclassencoder, which may have added with unknown words 
    int added = 0;
    if (sourceclassencoder->gethighestclass() > sourceclassdecoder->gethighestclass()) {
        for (unsigned int i = sourceclassdecoder->gethighestclass() + 1; i <= sourceclassencoder->gethighestclass(); i++) {
            added++;
            unsigned int cls = i;
            const string word = sourceclassencoder->added[cls];
            sourceclassdecoder->add(cls, word);            
            cerr << "NOTICE: Unknown word in input: '" << word << "'" << endl; //, assigning classes (" << cls << ", " << targetcls << ")" << endl;                        
        }
    }
    sourceclassencoder->added.clear();
    return added;    
}

void usage() {
    cerr << "Training usage: contextmoses -f source-traindatafile [-N|-X|-M] [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Training usage: contextmoses -F testdatafile [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Classifier types: (pick one)" << endl;
    cerr << " -N           N-Classifier Array, one classifier per pattern size group" << endl;
    cerr << " -X           Construction experts, one classifier per construction" << endl;
    cerr << " -M           Monolithic joined classifier, focus words are joined (-1)" << endl;
    cerr << "Scorehandling:" << endl;
    cerr << " -H " << endl;
	cerr << "       weighed            Apply classifier score as weight to original scores (default)" << endl;
    cerr << "       append             Append classifier score to translation score vector (make sure to specify an extra weight using -W)" << endl;
    cerr << "       replace            Use only classifier score, replacing translation table scores (make to specify only one weight using -W)" << endl;
    cerr << "       ignore             Ignore, do not use classifiers" << endl;
    cerr << "Options:" << endl;
    cerr << " -l [size]    Left context size" << endl;
    cerr << " -r [size]    Right context size" << endl;
    cerr << " -1           Represent the focus feature as a single entity, rather than individual tokens" << endl;
    cerr << " -C [id]      Classifier prefix." << endl;
    cerr << " -c [int]     Context threshold. Only create a classifier when at least this many different contexts exist. Defaults to 1." << endl;
    cerr << " -t [int]     Target threshold. Only create a classifier when at least this many different target options exist. Defaults to 1." << endl;
    cerr << " -a [float]   Accuracy threshold for Construction experts (-X), only experts with a leave-one-out accuracy higher than specified will be included. Value between 0 and 1. Defaults to 0 (no threshold)." << endl;
    cerr << " -x           disable exemplar weighting" << endl;
    cerr << " -O [options] Timbl options" << endl;
    cerr << " -1           Represent the focus feature as a single entity, rather than individual tokens" << endl;
    cerr << " -D           Enable debug" << endl;
    
    //cerr << "\t-C number                 Classifier mode" << endl;
    //cerr << "\t   1 - Local context with Classifier Array" << endl;
}

int main( int argc, char *argv[] ) {
    string sourceclassfile = "";
    string targetclassfile = "";    
    string classifierid="classifier";
    bool exemplarweights = true;
    
    static struct option long_options[] = {                        
       {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    string timbloptions;
    if (exemplarweights) {
        timbloptions = "-a 1 -s -F Tabbed";
    } else {
        timbloptions = "-a 1 -F Tabbed";
    }
    
    
    ClassifierType mode = CLASSIFIERTYPE_NONE;
    ScoreHandling scorehandling;
    
    bool TRAIN = false;
    bool TEST = false;
    bool debug = false;
    
    string trainfile = "";
    string testfile = "";
    string mosesphrasetable = "";
    string alignmodelfile = "";

    int leftcontextsize = 1;
    int rightcontextsize = 1;

    int contextthreshold = 1;
    int targetthreshold = 1;
    bool singlefocusfeature = false;
    double accuracythreshold = 0;
    
    
    char c;    
    string s;
    while ((c = getopt_long(argc, argv, "hd:S:T:C:xO:XNc:t:M1a:f:t:l:r:F:DH:m:",long_options,&option_index)) != -1) {
        switch (c) {
        case 0:
            if (long_options[option_index].flag != 0)
               break;
        case 'd':
        	alignmodelfile = optarg;
        	break;
        case 'h':
        	usage();
        	exit(0);
        case 'S':
            sourceclassfile = optarg;
            break;
        case 'T':
            targetclassfile = optarg;
            break;    
        case 'C':
            classifierid = optarg;
            break;
        case 'X':
            mode = CLASSIFIERTYPE_CONSTRUCTIONEXPERTS;
            break;
        case 'N':
            mode = CLASSIFIERTYPE_NARRAY;
            break;
        case 'M':
            mode = CLASSIFIERTYPE_MONO;
            break;
        case 'O':
            if (exemplarweights) {
                timbloptions = "-s -F Tabbed " + std::string(optarg);
            } else {
                timbloptions = "-F Tabbed " + std::string(optarg);
            }
            break;
        case 'c':
            contextthreshold = atoi(optarg);
            break;
        case 't':
            targetthreshold = atoi(optarg);
            break;
        case 'a':
            accuracythreshold = atof(optarg);
            break; 
        case 'x':
            exemplarweights = false;
            break;
        case '1':
            singlefocusfeature = true;
            break;
        case 'm':
            mosesphrasetable = optarg;
            break; 	
        case 'f':
            TRAIN = true;
            trainfile = optarg;
            break;
        case 'F':
            TEST = true;
            testfile = optarg;
            break;            
        case 'l':
            leftcontextsize = atoi(optarg);
            break;
        case 'r':
            rightcontextsize = atoi(optarg);
            break;
        case 'D':
            debug = true;
            break;  
        case 'H':
            s = optarg;
            if (s == "weighed") {
                scorehandling = SCOREHANDLING_WEIGHED;
            } else if (s == "append") {
                scorehandling = SCOREHANDLING_APPEND;
            } else if (s == "replace") {
                scorehandling = SCOREHANDLING_REPLACE;
            } else if (s == "ignore") {
                scorehandling = SCOREHANDLING_IGNORE;                
            } else {
                cerr << "Invalid value for -x: '" << s << "'" << endl;
                exit(2);    
            }
            break;   
        }
    }
    
    if ((!TRAIN && !TEST)) {
        usage();
        exit(2);
    }
    
    if ( (sourceclassfile.empty() || targetclassfile.empty()) ) {
        cerr << "ERROR: Specify class files (-S -T)" << endl;
        usage();
        exit(2);
    }
    
    if (mode == CLASSIFIERTYPE_NONE) {
        cerr << "ERROR: Choose a classifier type" << endl;
        usage();
        exit(2);
    }
    
    ClassDecoder * sourceclassdecoder = NULL;
    ClassDecoder * targetclassdecoder = NULL;
    ClassEncoder * sourceclassencoder = NULL;
    ClassEncoder * targetclassencoder = NULL;
    AlignmentModel * alignmodel = NULL;
    

    cerr << "Loading source class decoder " << sourceclassfile << endl;
	sourceclassdecoder = new ClassDecoder(sourceclassfile);

	cerr << "Loading target class decoder " << targetclassfile << endl;
	targetclassdecoder = new ClassDecoder(targetclassfile);   


		
    
    if (!alignmodelfile.empty()) {
        cerr << "Loading alignment model " << alignmodelfile << endl;
        alignmodel = new AlignmentModel(alignmodelfile,false,true,0, false);
    } else if (!mosesphrasetable.empty()) {
	    cerr << "Loading target class encoder " << targetclassfile << endl;
	    targetclassencoder = new ClassEncoder(targetclassfile);  

        cerr << "Loading source class encoder " << sourceclassfile << endl;
        sourceclassencoder = new ClassEncoder(sourceclassfile);
    
        cerr << "Loading moses phrasetable " << mosesphrasetable << endl;
        alignmodel = new AlignmentModel(mosesphrasetable, sourceclassencoder, targetclassencoder);    
    } else {
        cerr << "ERROR: No moses phrasetable (-m) or colibri alignment model (-d) specified!" << endl;
        exit(2);
    }
    
    int maxn = 0;
    for (t_alignmatrix::iterator iter  = alignmodel->alignmatrix.begin(); iter != alignmodel->alignmatrix.end(); iter++) {
        const int n = iter->first->n();
        if (n > maxn) maxn = n;
    }
    
    ClassifierInterface * classifiers = NULL;
    
    if ((TRAIN) && (!trainfile.empty())) {
        /*
        train) 
	    - read moses phrasetable or colibri alignment model
	    - read source-side training data
	    - match with phrasetable
		     - extract context and features
			    - add to classifier training data
	    - train classifiers	
		*/
		
		AlignmentModel * contextalignmodel = new AlignmentModel((unsigned char) leftcontextsize, (unsigned char) rightcontextsize);
		
		if (mode == CLASSIFIERTYPE_NARRAY) {
    
            cerr << "Initialising N-Array classifiers" << endl;
            classifiers = new NClassifierArray(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
            
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
    
            cerr << "Initialising construction expert classifiers" << endl;
            classifiers = new ConstructionExperts(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);    
		
		} else if (mode == CLASSIFIERTYPE_MONO) {
		
            if (!singlefocusfeature) {
                    cerr << "ERROR: Monolithic classifier only supported with single focus feature" << endl;
                    exit(2);
            }
            
            cerr << "Initialising monolithic classifier" << endl;
            classifiers = new MonoClassifier(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);  
            
        }
		
		cerr << "Extracting context on training set" << endl;
		
		const int BUFFERSIZE = 65536;
        unsigned char linebuffer[BUFFERSIZE];
        unsigned char tmpbuffer[BUFFERSIZE];
        
        ifstream *IN =  new ifstream( trainfile.c_str() );
        if (!IN->good()) {
        	cerr << "ERROR: Unable to open file " << trainfile << endl;
        	exit(5);
        }        
        vector<unsigned int> words;
        int sentence = 0;
        while (IN->good()) {
            const int linesize = readline(IN, linebuffer, BUFFERSIZE );            
            if (!IN->good()) break;
            
            sentence++;

            //if ((sentence % 1000 == 0) || (debug))  { 
            cerr << "\t@" << sentence << " ";
            //}
                            
            
            const int l = countwords(linebuffer, linesize);            
            if (l >= 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping!" << endl;
                continue;
            } else if (l == 0) {
            	cerr << "WARNING: Sentence " << sentence << " contains no words, skipping!" << endl;
                continue;
            }
            int foundcount = 0;    
                                    
            if (linesize > 0) {
                EncData line = EncData(linebuffer, linesize);                        
                for (unsigned char i = 0; ((i < l) && (i < 256)); i++) {
                    bool found;
                    unsigned char n = 1;
                    do {
                        found = false;
                        EncNGram * ngram = line.slice(i,n);    
                        const EncAnyGram * key = alignmodel->getsourcekey((const EncAnyGram *) ngram);
                        if (key != NULL) {
                            foundcount++;
                            found = true;
                            if (debug) cerr << "found match" << endl;
                            //match found!
                            const EncAnyGram * incontext = alignmodel->addcontext(&line, (const EncAnyGram * ) ngram, (int) i, leftcontextsize, rightcontextsize);
                            //see if this one already exists:
                            const EncAnyGram * contextkey = contextalignmodel->getsourcekey(incontext, false);
                            
                            
                            //add to context-aware alignment model (classifier training data will be constructed on the basis of this)
                            for (t_aligntargets::iterator iter = alignmodel->alignmatrix[key].begin(); iter !=  alignmodel->alignmatrix[key].end(); iter++) {
                                const EncAnyGram * targetgram = iter->first;
                                const double score = (iter->second[0] < 0) ? pow(exp(1), iter->second[0]) : iter->second[0]; //no logprob
                                contextalignmodel->addextractedpattern(key, targetgram, score, 1, (contextkey != NULL) ? contextkey : incontext );
                            }

                            //add to classifier 
                            /*
                            if (mode == CLASSIFIERTYPE_NARRAY) {
                                if (debug) cerr << "adding to n-array classifier" << endl;                  
                                ((NClassifierArray *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[key], leftcontextsize, rightcontextsize, sourceclassdecoder, targetclassdecoder);                                                
                            } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
                                if (debug) cerr << "adding to expert classifier" << endl;
                                ((ConstructionExperts *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[key], leftcontextsize, rightcontextsize, sourceclassdecoder, targetclassdecoder);                        
                            } else if (mode == CLASSIFIERTYPE_MONO) {
                                if (debug) cerr << "adding to monolithic classifier" << endl;
                                ((MonoClassifier *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[key], leftcontextsize, rightcontextsize, sourceclassdecoder, targetclassdecoder);
                            }*/
                            
                            if (contextkey != NULL) { 
                                delete incontext; 
                            }
                        } 
                        delete ngram;                  
                        n++;
                    } while ((found) && (i+n <= l) && (n <= maxn));  
                }
            }
            
            cerr << foundcount << endl;
            
            
        }
		
		
		cerr << "Building classifiers" << endl;
		
		
        
        if (mode == CLASSIFIERTYPE_NARRAY) {
            if (debug) cerr << "Building n-array classifier" << endl;                  
            ((NClassifierArray *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder);                                                            
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
            if (debug) cerr << "Building expert classifier" << endl;
            ((ConstructionExperts *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder);
        } else if (mode == CLASSIFIERTYPE_MONO) {
            if (debug) cerr << "Building monolithic classifier" << endl;
            ((MonoClassifier *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder);
        }
		
		
		cerr << "Training classifiers" << endl;
        cerr << "   Timbl options: " << timbloptions << endl; 
        
        if (mode == CLASSIFIERTYPE_NARRAY) {

            ((NClassifierArray *) classifiers)->train(timbloptions);
        
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
        
            
            ((ConstructionExperts *) classifiers)->accuracythreshold = accuracythreshold;
            ((ConstructionExperts *) classifiers)->train(timbloptions);
        
        } else if (mode == CLASSIFIERTYPE_MONO) {
            
            ((MonoClassifier *) classifiers)->train(timbloptions);
            
        }
        
        writeclassifierconf(classifierid, mode, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
        cerr << "Training all done" << endl;
    }
    
    if (TEST) {   
        if (!classifierid.empty()) {
            
            //Load classifiers
            int scorecount = 0;
            //cerr << "Computing reverse index for translation table" << endl;
            //transtable->computereverse(); //not necessary 
            cerr << "Loading classifiers" << endl;
            cerr << "   ID: " << classifierid << endl;
            cerr << "   Timbl options: " << timbloptions << endl;
            cerr << "   Score handling: ";
            if (scorehandling == SCOREHANDLING_WEIGHED) {
                cerr << "weighed" << endl;
            } else if (scorehandling == SCOREHANDLING_APPEND) {
                cerr << "append" << endl;
            } else if (scorehandling == SCOREHANDLING_REPLACE) {
                cerr << "replace" << endl;
            }
            /*
            int contextthreshold; //will be set by getclassifiertype
            int targetthreshold; //will be set by getclassifiertype
            bool exemplarweights; //will be set by getclassifiertype
            bool singlefocusfeature; //will be set by getclassifiertype
            */        
            ClassifierType classifiertype = getclassifierconf(classifierid, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
            if (classifiertype == CLASSIFIERTYPE_NARRAY) {        
                classifiers = (ClassifierInterface*) new NClassifierArray(classifierid, leftcontextsize,rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
                classifiers->load(timbloptions, sourceclassdecoder, targetclassencoder, debug);
            } else if (classifiertype == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
                classifiers = (ClassifierInterface*) new ConstructionExperts(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
                 classifiers->load(timbloptions, sourceclassdecoder, targetclassencoder, debug);                    
            } else if (classifiertype == CLASSIFIERTYPE_MONO) {
                if (!singlefocusfeature) {
                    cerr << "ERROR: Monolithic classifier only supported with single focus feature" << endl;
                    throw InternalError();
                }
                classifiers = (ClassifierInterface*) new MonoClassifier(classifierid, leftcontextsize,rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
                classifiers->load(timbloptions, sourceclassdecoder, targetclassencoder, debug);            
            } else {
                cerr << "ERROR: Undefined classifier type:" << classifiertype << endl;
                throw InternalError();            
            }
            if (debug) classifiers->enabledebug(debug,sourceclassdecoder, targetclassdecoder); 
            
            /*
            - read moses phrasetable or colibri alignment model
	        - read test data
	        - match with phrasetable
		        - extract context and features
			        - classify
			        - replace source words with unique IDs
			        - write intermediate phrasetable with IDs instead of source words, store map IDs->words
	        - run decoder with intermediate phrasetable
	        - read decoding results
	        - lookup 
	        */
            
            cerr << "Reading test file" << endl;
    
            ifstream *IN =  new ifstream( testfile.c_str() );
            if (!IN->good()) {
            	cerr << "ERROR: Unable to open file " << trainfile << endl;
            	exit(5);
            }        
    
            ofstream *TMPTEST = new ofstream( "tmp.txt" ); //intermediate test file (IDs instead of words)
            ofstream *TMPTABLE = new ofstream( "tmp.phrasetable" ); //intermediate phrase table

            string input;
            unsigned char buffer[8192]; 
            int size;
            int sentence = 0;
            while (getline(*IN, input)) {        
                if (input.length() > 0) {
                    sentence++;                    
                    cerr << "INPUT: " << input << endl;
                    if (debug) cerr << "Processing input" ;        
                    size = sourceclassencoder->encodestring(input, buffer, true, true) - 1; //weird patch: - 1  to get n() right later
                    
                    if (debug) cerr << " (" << size << ") " << endl;
                    const EncData * line = new EncData(buffer,size);
                    if (debug) cerr << "Processing unknown words" << endl; 
                    addunknownwords(sourceclassencoder, sourceclassdecoder, targetclassencoder, targetclassdecoder);
                    //if (debug >= 1) cerr << "Setting up decoder (" << inputdata->length() << " stacks)" << endl;
                    
                    const int l = line->length();
                    
                    for (int i = 0; i < l; i++) {
                        if (i > 0) *TMPTEST << " ";
                        *TMPTEST << sentence << "_" << i;
                    }
                    *TMPTEST << endl;
                                            
                    for (unsigned char i = 0; ((i < l) && (i < 256)); i++) {
                        bool found;
                        unsigned char n = 1;
                        do {
                            found = false;
                            EncNGram * ngram = line->slice(i,n);
                            
                            stringstream ss;
                            for (int j = i; j < j+n; j++) {
                                if (j > i) ss << " ";
                                ss << sentence << "_" << j; 
                            } 
                            const string encodedngram = ss.str();
                            
                                  
                            const EncAnyGram * key = alignmodel->getsourcekey((const EncAnyGram *) ngram);
                            if (key != NULL) {
                                //match found!
                                const EncAnyGram * incontext = alignmodel->addcontext(line, (const EncAnyGram * ) ngram, (int) i, leftcontextsize, rightcontextsize);                                
                                alignmodel->alignmatrix[key];
                                
                                t_aligntargets * reftranslationoptions = &(alignmodel->alignmatrix[key]);
                                t_aligntargets translationoptions;
                                
                                //are there enough targets for this source to warrant a classifier?
                                if (alignmodel->alignmatrix[key].size() > scorecount) scorecount = alignmodel->alignmatrix[key].size(); 
                                if (alignmodel->alignmatrix[key].size() >= targetthreshold) {
                                    translationoptions = classifiers->classifyfragment(key, incontext, *reftranslationoptions, scorehandling, leftcontextsize, rightcontextsize);
                                } else {
                                    translationoptions = *reftranslationoptions;
                                }
                                
                                //write intermediate phrasetable
                                for (t_aligntargets::iterator iter = translationoptions.begin(); iter != translationoptions.end(); iter++) {
                                    *TMPTABLE << encodedngram << " ||| " << iter->first->decode(*targetclassdecoder) << " ||| ";
                                    for (vector<double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                                        if (iter2 != iter->second.begin()) *TMPTABLE << " ";
                                        *TMPTABLE << *iter2;
                                    }
                                    *TMPTABLE << endl;
                                }
                                
                                delete incontext;
                            }  
                            delete ngram;                  
                            n++;
                        } while (found);  
                        
                        
                        
                    }    
                                    
                    
                }
            }      

            TMPTABLE->close();
            TMPTEST->close(); 

            /*cerr << "Updating moses configuration..." << endl;
            ifstream *MOSESINI =  new ifstream( "moses.ini" );
            if (!IN->good()) {
            	cerr << "ERROR: Unable to open moses.ini" << endl;
            	exit(5);
            } */       

            cerr << "Invoking moses: ";
            stringstream ss;
            for (int i = 0; i < scorecount; i++) {
                if (i > 0) ss << " ";
                ss << 1;
            }
            stringstream cmd;
            cmd << "moses -config moses.ini -ttable-file \"0 0 0 " << scorecount << " tmp.phrasetable\" -weight-t \"" << ss.str() << "\" < tmp.txt";          
            cerr << cmd << endl;  
            system(cmd.str().c_str());

        }  
    }
}
