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
    cerr << "Training usage: contextmoses -f source-traindatafile -g target-traindatafile [-N|-X|-M] [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Test usage: contextmoses -F testdatafile [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Training usage with keywords: contextmoses -f source-traindatafile -g target-traindatafile -s sourcepatternmodel -E targetpatternmodel -k -X -d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Classifier types: (pick one, for training only)" << endl;
    cerr << " -N           N-Classifier Array, one classifier per pattern size group" << endl;
    cerr << " -X           Construction experts, one classifier per construction" << endl;
    cerr << " -M           Monolithic joined classifier, focus words are joined (-1)" << endl;
    cerr << " -I           Ignore classifier (this one can be specified when testing only, in which case the classifier will be ignored)" << endl;
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
    cerr << " -i [int]     Instance threshold for Construction experts (-X), prune all classifiers with less instances than this threshhold." << endl;
    cerr << " -k           enable global context keywords (only for -X, required an alignment model loaded with keywords (-d), and pattern models (-s -E)" << endl;
    cerr << " -K           probability threshold p(keyword|source,target) (default: 1e-99)" << endl;
    cerr << " -L [int]     best n keywords (default: 25)" << endl;
    cerr << " -s [file]    Source-side pattern model (needed for -k)" << endl;
    cerr << " -E [file]    Target-side pattern model (needed for -k)" << endl;
    cerr << " -x           disable exemplar weighting" << endl;
    cerr << " -O [options] Timbl options" << endl;
    cerr << " -1           Represent the focus feature as a single entity, rather than individual tokens" << endl;
    cerr << " -D           Enable debug" << endl;
    cerr << " -p [int]     Field of the score vector in which the forward probability p(t|s) is stored (default: 3)" << endl;
    cerr << " -e [float]   Small epsilon value used as score for unencountered options when score handling is set to append mode (default:  0.000001) " << endl;
    cerr << " -q           Skip decoder" << endl;
    cerr << " -o           Output prefix (default: tmp)" << endl;
    
    
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
    
    string traintimbloptions;
    string testtimbloptions = "-a 0 -F Tabbed";
    if (exemplarweights) {
        traintimbloptions = "-a 0 -s -F Tabbed";
    } else {
        traintimbloptions = "-a 0 -F Tabbed";
    }
    
    
    ClassifierType mode = CLASSIFIERTYPE_NONE;
    ScoreHandling scorehandling = SCOREHANDLING_WEIGHED;
    
    bool TRAIN = false;
    bool TEST = false;
    bool debug = false;
    
    bool DOKEYWORDS = false;
    
    string trainfile = "";
    string targettrainfile = "";
    string testfile = "";
    string mosesphrasetable = "";
    string alignmodelfile = "";
    string sourcepatternmodelfile = "";
    string targetpatternmodelfile = "";

    int leftcontextsize = 1;
    int rightcontextsize = 1;

    int contextthreshold = 1;
    int targetthreshold = 1;
    bool singlefocusfeature = false;
    double accuracythreshold = 0;
    int instancethreshold = 0;
    double keywordprobthreshold = 1e-99;
    int bestnkeywords = 25;
    
    bool timbloptionsset = false;
    double appendepsilon = 0.000001;
    
    int ptsfield = 3; //1-indexed
    
    bool skipdecoder = false;
    
    string outputprefix = "tmp";
    
    char c;    
    string s;
    while ((c = getopt_long(argc, argv, "hd:S:T:C:xO:XNc:t:M1a:f:g:t:l:r:F:DH:m:Ip:e:qo:i:kK:s:E:L:",long_options,&option_index)) != -1) {
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
        case 'I':
            mode = CLASSIFIERTYPE_IGNORE;
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
            timbloptionsset = true;
            if (exemplarweights) {
                traintimbloptions = "-s -F Tabbed " + std::string(optarg);
                testtimbloptions = "-s -F Tabbed " + std::string(optarg);
            } else {
                traintimbloptions = "-F Tabbed " + std::string(optarg);
                testtimbloptions = "-F Tabbed " + std::string(optarg);
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
        case 'i':
            instancethreshold = atoi(optarg);
            break;
        case 'x':
            if (timbloptionsset) {
                cerr << "ERROR: Only specify -x before -O, not after" << endl;
                exit(2);
            }
            traintimbloptions = "-a 0 -F Tabbed";
            exemplarweights = false;
            break;
        case '1':
            singlefocusfeature = true;
            break;
        case 'm':
            mosesphrasetable = optarg;
            break; 	
        case 's':
            sourcepatternmodelfile = optarg;
            break;
        case 'E':
            targetpatternmodelfile = optarg;
            break;
        case 'e':
            appendepsilon = atof(optarg);
            break;
        case 'f':
            TRAIN = true;
            trainfile = optarg;
            break;
        case 'g':
            TRAIN = true;
            targettrainfile = optarg;
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
        case 'p':
            ptsfield = atoi(optarg);
            break;
        case 'D':
            debug = true;
            break;  
        case 'q':
            skipdecoder = true;
            break;
        case 'o':
            outputprefix = optarg;
            break;
        case 'k':
            DOKEYWORDS = true;
            break;
        case 'K':
            keywordprobthreshold = atof(optarg);
            break;
        case 'L':
            bestnkeywords = atoi(optarg);
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
    
    if ((DOKEYWORDS) && (TRAIN) && ((sourcepatternmodelfile.empty()) || (targetpatternmodelfile.empty()))) {
        cerr << "ERROR: source and target pattern model must be specified when using -k";
        usage();
        exit(2);
    }
    if (DOKEYWORDS) mode = CLASSIFIERTYPE_CONSTRUCTIONEXPERTS;


    const string tmptestfile = outputprefix + ".txt";
    const string tmptablefile = outputprefix + ".phrasetable";

    ClassDecoder * sourceclassdecoder = NULL;
    ClassDecoder * targetclassdecoder = NULL;
    ClassEncoder * sourceclassencoder = NULL;
    ClassEncoder * targetclassencoder = NULL;
    AlignmentModel * alignmodel = NULL;
    SelectivePatternModel * sourcepatternmodel = NULL;
    SelectivePatternModel * targetpatternmodel = NULL;
    
    bool testexists = false;
    if (TEST) {
        ifstream TEST1(tmptestfile );
        ifstream TEST2( tmptablefile );
        testexists = (TEST1 && TEST2);
        if (testexists) {                       
            TEST1.close();
            TEST2.close();
        }
    }
    

    cerr << "Loading source class decoder " << sourceclassfile << endl;
	sourceclassdecoder = new ClassDecoder(sourceclassfile);

	cerr << "Loading target class decoder " << targetclassfile << endl;
	targetclassdecoder = new ClassDecoder(targetclassfile);  

    if ((DOKEYWORDS) && (TRAIN)) {
        GraphFilter graphfilter;
        sourcepatternmodel = new SelectivePatternModel(sourcepatternmodelfile, graphfilter, true,true);
        targetpatternmodel = new SelectivePatternModel(targetpatternmodelfile, graphfilter, true,true);
    }

    int maxn = 0;
	
    if (!alignmodelfile.empty()) {

            cerr << "Loading alignment model " << alignmodelfile << endl;
            alignmodel = new AlignmentModel(alignmodelfile,false,ptsfield, true,0, false); 
            cerr << "\tLoaded " << alignmodel->size() << " source patterns";
            if (alignmodel->keywords.size() > 0) {
                cerr << ", with keywords for " << alignmodel->keywords.size() << " of them";
            } else if (DOKEYWORDS) {
                cerr << "WARNING: Keywords are enabled but alignmodel model " << alignmodelfile << " contains no keywords!!!" << endl;
            }
            cerr << "." << endl;
    }


    if (!  ((!TRAIN) && (TEST) && (testexists))) {

        cerr << "Loading target class encoder " << targetclassfile << endl;
        targetclassencoder = new ClassEncoder(targetclassfile);  

        cerr << "Loading source class encoder " << sourceclassfile << endl;
        sourceclassencoder = new ClassEncoder(sourceclassfile);

        if (!mosesphrasetable.empty()) {
            if (DOKEYWORDS) {
                cerr << "ERROR: Global context features are enabled, need a colibri alignment model with keywords instead of a moses phrasetables" << endl;
                exit(2);
            }
        

            cerr << "Loading moses phrasetable " << mosesphrasetable << endl;
            alignmodel = new AlignmentModel(mosesphrasetable, sourceclassencoder, targetclassencoder, true, ptsfield);
            
            /*if (debug) {
                cerr << "Outputting phrasetable"<< endl;
                for (t_alignmatrix::iterator iter = alignmodel->alignmatrix.begin(); iter != alignmodel->alignmatrix.end(); iter++) {
                    const EncAnyGram * source = iter->first; 
                    for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                        const EncAnyGram * target = iter2->first;
                        cerr << source->decode(*sourceclassdecoder) << " ||| " << target->decode(*targetclassdecoder) << " ||| ";
                        for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                            cerr << *iter3 << " ";
                        }
                        cerr << endl;
                    }
                }
            }*/ 
            
        } else if (alignmodelfile.empty()) {
            cerr << "ERROR: No moses phrasetable (-m) or colibri alignment model (-d) specified!" << endl;
            exit(2);
        }
        
        for (t_alignmatrix::iterator iter  = alignmodel->alignmatrix.begin(); iter != alignmodel->alignmatrix.end(); iter++) {
            const int n = iter->first->n();
            if (n > maxn) maxn = n;
        }        
    }
    
    

    ClassifierInterface * classifiers = NULL;
    
    if ((TRAIN) && (!trainfile.empty()) && (!targettrainfile.empty())) {
    
        /*
        train) 
	    - read moses phrasetable or colibri alignment model
	    - read source-side training data
	    - match with phrasetable
		     - extract context and features
			    - add to classifier training data
	    - train classifiers	
		*/
		
		if (mode == CLASSIFIERTYPE_NONE) {
            cerr << "ERROR: Choose a classifier type" << endl;
            usage();
            exit(2);
        }
    
	
        //new alignment model to be build
		AlignmentModel * contextalignmodel = new AlignmentModel((unsigned char) leftcontextsize, (unsigned char) rightcontextsize, ptsfield);
	

		if (mode == CLASSIFIERTYPE_NARRAY) {
    
            cerr << "Initialising N-Array classifiers" << endl;
            classifiers = new NClassifierArray(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature);
            
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
    
            cerr << "Initialising construction expert classifiers ";
            if (DOKEYWORDS) {
                cerr << "with keywords" << endl;
            } else {
                cerr << "without keywords" << endl;
            }
            classifiers = new ConstructionExperts(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature, DOKEYWORDS, keywordprobthreshold, bestnkeywords); 
		
		} else if (mode == CLASSIFIERTYPE_MONO) {
		
            if (!singlefocusfeature) {
                    cerr << "ERROR: Monolithic classifier only supported with single focus feature" << endl;
                    exit(2);
            }
            
            cerr << "Initialising monolithic classifier" << endl;
            classifiers = new MonoClassifier(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature);
            
        } else if (mode == CLASSIFIERTYPE_IGNORE) {
            cerr << "ERROR: Can't ignore classifier when training!" << endl;
            exit(2);
        }
		
		cerr << "Extracting context on training set" << endl;
		
		const int BUFFERSIZE = 65536;
        unsigned char linebuffer[BUFFERSIZE];
        unsigned char targetlinebuffer[BUFFERSIZE];
        unsigned char tmpbuffer[BUFFERSIZE];
        
        ifstream *INSOURCE =  new ifstream( trainfile.c_str() );
        if (!INSOURCE->good()) {
        	cerr << "ERROR: Unable to open file " << trainfile << endl;
        	exit(5);
        }        
        ifstream *INTARGET =  new ifstream( targettrainfile.c_str() );
        if (!INTARGET->good()) {
        	cerr << "ERROR: Unable to open file " << trainfile << endl;
        	exit(5);
        }                
       
        t_keywordflags flaggedkeywords;

        unsigned int totalkwcount = 0;

        vector<unsigned int> words;
        int sentence = 0;
        while (INSOURCE->good()) {
            sentence++;
            int linesize = readline(INSOURCE, linebuffer, BUFFERSIZE );
            int targetlinesize = readline(INTARGET, targetlinebuffer, BUFFERSIZE );    
                    

            if (!INSOURCE->good()) linesize--; //silly fix, don't know why, but works


            if ((sentence % 1000 == 0) || (debug))  { 
                cerr << "\t@" << sentence << endl;
            }
                            
            
            const int l = countwords(linebuffer, linesize);            
            if (l == 0) {
            	cerr << "WARNING: Sentence " << sentence << " contains no words, skipping!" << endl;
                continue;
            }
            int foundcount = 0;    

            const int ltarget = countwords(targetlinebuffer,targetlinesize);

            if (targetpatternmodel->reverseindex.empty()) {
                cerr << "No reverse index loaded in target-side patternmodel" << endl;
                throw InternalError();
            }


                                    
            if (linesize > 0) {
                EncData line = EncData(linebuffer, linesize);
                EncData targetline = EncData(targetlinebuffer, targetlinesize);                        
                for (int i = 0; i < l; i++) {
                    bool found;
                    unsigned char n = 1;
                    do {
                        found = false;
                        EncNGram * ngram = line.slice(i,n);    
                        const EncAnyGram * key = alignmodel->getsourcekey((const EncAnyGram *) ngram);
                        if (key != NULL) {
                            foundcount++;
                            found = true;
                            if (debug) cerr << "found match @" << sentence << " " << i << ":" << (int) n << endl;
                            //match found!
                            const EncAnyGram * incontext = key; //no context
                            if ((leftcontextsize > 0) || (rightcontextsize > 0)) {
                                incontext = alignmodel->addcontext(&line, (const EncAnyGram * ) ngram, (int) i, leftcontextsize, rightcontextsize);
                            }
                            //see if this one already exists:
                            const EncAnyGram * contextkey = contextalignmodel->getsourcekey(incontext, false);

                            //const EncAnyGram * focuskey = alignmodel->getkeywordkey( (const EncAnyGram *) ngram);
                            //if (focuskey == NULL) focuskey = key;
                            
                            //see what targets in the target sentence match with the translation options (only one alignment is right, but if there is ambiguity we add them all as we don't have alignment data at this point).. 
                            bool targetfound = false;
                            for (t_aligntargets::iterator iter = alignmodel->alignmatrix[key].begin(); iter !=  alignmodel->alignmatrix[key].end(); iter++) {                            

                                if (debug) cerr << "    processing target";
                                const EncAnyGram * targetgram = iter->first;

                                const EncAnyGram * targetkey = targetpatternmodel->getkey(targetgram);
                                if (targetkey == NULL) {
                                    if (debug) cerr << " ... not found in target pattern model" << endl;
                                    continue;
                                } else {                                    
                                    if (debug) {
                                        cerr << "... target found.. ";
                                        if (DOKEYWORDS){
                                            if (key == NULL) { //can't happen
                                                cerr << "no keywords found for source " << ngram->decode(*sourceclassdecoder);
                                            } else {
                                                if (alignmodel->keywords[key].count(targetgram)) {
                                                    cerr << alignmodel->keywords[key][targetgram].size() << " possible keywords found";
                                                } else {
                                                    cerr << "no keywords found for target";
                                                }
                                            }
                                        }
                                        cerr << endl;
                                    } 
                                }
                                

                                if ((DOKEYWORDS) && (targetpatternmodel->reverseindex[sentence].count(targetkey))) {  //line.contains((const EncNGram *) targetgram)) { //use targetpatternmodel and reverse index!!
                                    if (debug) cerr << "\t\tCounting keywords" << endl;
                                    //loop over global context keywords and flag presence, store in separate datastructure: flaggedkeywords
                                    if ((alignmodel->keywords.count(key)) && (alignmodel->keywords[key].count(targetgram))) {
                                        unordered_set<const EncAnyGram *> keywords;
                                        //reverse: loop over patterns in
                                        //sentence and match each with
                                        //keywords
                                        if (debug) cerr << "\t\tKeywords for " << key->decode(*sourceclassdecoder) << " -> " << targetgram->decode(*targetclassdecoder) << " = " << alignmodel->keywords[key][targetgram].size();
                                        for (unordered_set<const EncAnyGram *>::iterator kwiter = sourcepatternmodel->reverseindex[sentence].begin(); kwiter != sourcepatternmodel->reverseindex[sentence].end(); kwiter++) { 
                                            const EncAnyGram * keyword = *kwiter; //candidate keyword
                                            if (alignmodel->keywords[key][targetgram].count(keyword)) { //check if this is a keyword
                                                if (alignmodel->keywords[key][targetgram][keyword] >= keywordprobthreshold) {
                                                    if (debug) cerr << " | " << keyword->decode(*sourceclassdecoder);
                                                    totalkwcount++;
                                                    keywords.insert(keyword);
                                                }
                                            }
                                        }
                                        if (debug) cerr << endl;


                                        //for (unordered_map<const EncAnyGram *, double>::iterator kwiter = alignmodel->keywords[key][targetgram].begin(); kwiter != alignmodel->keywords[key][targetgram].end(); kwiter++) { //problem: far too many keywords!!
                                        //    const EncAnyGram * keyword = kwiter->first;
                                        //    if (sourcepatternmodel->reverseindex[sentence].count(keyword)) {  //(line.contains((const EncNGram *) keyword)) { //TODO: use sourcepatternmodel and reverse index!!
                                        //        if (kwiter->second >= keywordprobthreshold) {
                                        //            keywords.insert(keyword);
                                        //        }
                                        //    }
                                        //}
                                        if (debug) cerr << "    found  " << keywords.size() << " of " << alignmodel->keywords[key][targetgram].size() << " keywords" << endl;                                            
                                        contextalignmodel->keywords[key][targetgram] = alignmodel->keywords[key][targetgram]; //simply copy keywords to contextalignmodel
                                        flaggedkeywords[(contextkey != NULL) ? contextkey : incontext][targetgram].push_back(keywords);
                                    }

                                    //add to context-aware alignment model (classifier training data will be constructed on the basis of this)
                                    targetfound = true;
                                    const double score = (exemplarweights) ?  (  (iter->second[0] < 0) ? pow(exp(1), iter->second[0]) : iter->second[0] ) : 1; //no logprob
                                    contextalignmodel->addextractedpattern(key, targetgram, score, 1, (contextkey != NULL) ? contextkey : incontext );
                                }
                                
                            }
                            
                            if (targetfound) {
                                contextalignmodel->sourcecontexts[key].insert((contextkey != NULL) ? contextkey : incontext);
                            } 
                            
                            
                            if ((contextkey != NULL) && (contextkey != key)) { 
                                delete incontext; 
                            }
                        } 
                        delete ngram;                  
                        n++;
                    } while ((found) && (i+n <= l) && (n <= maxn));  
                }
            }
            

            
            
        }
	

        if (DOKEYWORDS) {
            cerr << "Keywords found for " << flaggedkeywords.size() << " different contexts. Total keyword tokens counted: " << totalkwcount << endl;
        }

        
        if (mode == CLASSIFIERTYPE_NARRAY) {
            cerr << "Building n-array classifier" << endl;                  
            ((NClassifierArray *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder);                                                            
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
            cerr << "Building expert classifier" << endl;
            ((ConstructionExperts *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder, &flaggedkeywords);
        } else if (mode == CLASSIFIERTYPE_MONO) {
            cerr << "Building monolithic classifier" << endl;
            ((MonoClassifier *) classifiers)->build(contextalignmodel, sourceclassdecoder, targetclassdecoder);
        }
		
		
		cerr << "Training classifiers" << endl;
        cerr << "   Timbl options: " << traintimbloptions << endl; 
        
        if (mode == CLASSIFIERTYPE_NARRAY) {

            ((NClassifierArray *) classifiers)->train(traintimbloptions);
        
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
        
            
            ((ConstructionExperts *) classifiers)->accuracythreshold = accuracythreshold;
            ((ConstructionExperts *) classifiers)->instancethreshold = instancethreshold;
            ((ConstructionExperts *) classifiers)->train(traintimbloptions);
        
        } else if (mode == CLASSIFIERTYPE_MONO) {
            
            ((MonoClassifier *) classifiers)->train(traintimbloptions);
            
        }
        
        writeclassifierconf(classifierid, mode, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
        cerr << "Training all done" << endl;
    }
    
    if (TEST) {   
        if (!classifierid.empty()) {
            
            //Debug: sanity check  //TODO: REMOVE
            cerr << "sanity check" << endl;
            cerr << alignmodel->keywords.size() << endl;
            for (t_keywords::iterator kwiter0 = alignmodel->keywords.begin(); kwiter0 != alignmodel->keywords.end(); kwiter0++) {
                for (t_keywords_source::iterator kwiter = kwiter0->second.begin(); kwiter != kwiter0->second.end(); kwiter++) {
                    for (unordered_map<const EncAnyGram *, double>::iterator kwiter2 = kwiter->second.begin(); kwiter2 != kwiter->second.end(); kwiter2++) {
                        const EncAnyGram * keyword = kwiter2->first;
                        keyword->decode(*sourceclassdecoder);
                    }
                }
            }
            
            int scorecount = 0;
            cerr << "Score handling: ";
            if (scorehandling == SCOREHANDLING_WEIGHED) {
                cerr << "weighed" << endl;
                scorecount = 5;
            } else if (scorehandling == SCOREHANDLING_APPEND) {
                cerr << "append" << endl;
                scorecount = 6;
            } else if (scorehandling == SCOREHANDLING_REPLACE) {
                cerr << "replace" << endl;
                scorecount = 1;
            } else if (scorehandling == SCOREHANDLING_IGNORE) {
                cerr << "ignore" << endl;
                scorecount = 5;                
            }
           
            if (!testexists) {
            
                //Load classifiers
               
                //cerr << "Computing reverse index for translation table" << endl;
                //transtable->computereverse(); //not necessary 
                cerr << "Loading classifiers" << endl;
                cerr << "   ID: " << classifierid << endl;
                cerr << "   Timbl options: " << testtimbloptions << endl;

                /*
                int contextthreshold; //will be set by getclassifiertype
                int targetthreshold; //will be set by getclassifiertype
                bool exemplarweights; //will be set by getclassifiertype
                bool singlefocusfeature; //will be set by getclassifiertype
                */        
                ClassifierType classifiertype;
                if (mode != CLASSIFIERTYPE_IGNORE) {
                    classifiertype = getclassifierconf(classifierid, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
                    if (classifiertype == CLASSIFIERTYPE_NARRAY) {        
                        classifiers = (ClassifierInterface*) new NClassifierArray(classifierid, leftcontextsize,rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature);
                        classifiers->load(testtimbloptions, sourceclassdecoder, targetclassencoder, debug);
                    } else if (classifiertype == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
                        classifiers = (ClassifierInterface*) new ConstructionExperts(classifierid, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature, DOKEYWORDS, keywordprobthreshold, bestnkeywords);
                         classifiers->load(testtimbloptions, sourceclassdecoder, targetclassencoder, debug);                    
                    } else if (classifiertype == CLASSIFIERTYPE_MONO) {
                        if (!singlefocusfeature) {
                            cerr << "ERROR: Monolithic classifier only supported with single focus feature" << endl;
                            throw InternalError();
                        }
                        classifiers = (ClassifierInterface*) new MonoClassifier(classifierid, leftcontextsize,rightcontextsize, contextthreshold, targetthreshold, ptsfield, appendepsilon, exemplarweights, singlefocusfeature);
                        classifiers->load(testtimbloptions, sourceclassdecoder, targetclassencoder, debug);
                    } else if (classifiertype == CLASSIFIERTYPE_IGNORE) {
                        cerr << "Ignoring classifiers" << endl;                 
                    } else {
                        cerr << "ERROR: Undefined classifier type:" << classifiertype << endl;
                        throw InternalError();            
                    }
                    
                    if (debug && classifiers != NULL) classifiers->enabledebug(2,sourceclassdecoder, targetclassdecoder);
                } 
                
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
            

                ofstream *TMPTEST = new ofstream( tmptestfile ); //intermediate test file (IDs instead of words)
                ofstream *TMPTABLE = new ofstream( tmptablefile ); //intermediate phrase table

                string input;
                unsigned char buffer[8192]; 
                int size;
                int sentence = 0;
                
                unsigned int sourcefragmentcount = 0;
                int unknowncount = 0;
                int changedcount = 0;
                
                int maxn = 0;
                for (t_alignmatrix::iterator iter = alignmodel->alignmatrix.begin(); iter != alignmodel->alignmatrix.end(); iter++) {
                    const int n = iter->first->n();
                    if (n > maxn) maxn = n; 
                }
                
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
                                                
                        for (int i = 0; ((i < l) && (i < 256)); i++) {
                            bool found;
                            int n = 1;
                            while ((n <= maxn) && (i + n <= l))  {
                                found = false;
                                
                                const EncNGram * ngram = line->slice(i,n);
                                                                
                                stringstream ss;
                                for (int j = i; j < i+n; j++) {
                                    if (j > i) ss << " ";
                                    ss << sentence << "_" << j; 
                                } 
                                const string encodedngram = ss.str();
                                
                                      
                                const EncAnyGram * key = alignmodel->getsourcekey((const EncAnyGram *) ngram);
                                if ((n == 1) && (key == NULL)) {
                                    if (debug) cerr << "found unknown word '" << ngram->decode(*sourceclassdecoder) << "' (s=" << sentence << ",i=" << i << ",n=" << (int) n << "): " << ngram->decode(*sourceclassdecoder) << endl;
                                    //unknown word!! Add to phrasetable
                                    unknowncount += 1;
                                    
                                    *TMPTABLE << encodedngram << " ||| " << ngram->decode(*sourceclassdecoder) << " ||| ";
                                    for (int j = 0; j < scorecount; j++) {
                                        if (j > 0) *TMPTABLE << " ";
                                        if (j == 4) {
                                            *TMPTABLE << "2.718";
                                        } else {
                                            *TMPTABLE << pow(exp(1),-100); //unknown words : -100 (default weight used by moses
                                        }
                                    }                 
                                    *TMPTABLE << endl;                                     
                                } else if (key != NULL) {
                                    //match found!
                                    sourcefragmentcount += 1;
                                    
                                    found = true;
                                    if (debug) cerr << "found match (s=" << sentence << ",i=" << i << ",n=" << (int) n << "): " << ngram->decode(*sourceclassdecoder) << endl;
                                    
                                    const EncAnyGram * incontext = key;
                                    if ((leftcontextsize > 0) || (rightcontextsize > 0)) {
                                        incontext = alignmodel->addcontext(line, (const EncAnyGram * ) ngram, (int) i, leftcontextsize, rightcontextsize);                                
                                    }

                                    
                                    t_aligntargets * reftranslationoptions = &(alignmodel->alignmatrix[key]);
                                    if (debug) {
                                        for (t_aligntargets::iterator iter = reftranslationoptions->begin(); iter != reftranslationoptions->end(); iter++) {
                                            cerr << "BEFORE CLASSIFICATION: " << encodedngram << " ||| " << iter->first->decode(*targetclassdecoder) << " ||| ";
                                            for (vector<double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                                                if (iter2 != iter->second.begin()) cerr << " ";
                                               cerr << pow(exp(1), *iter2);
                                            }
                                            cerr << endl;
                                        }
                                    }
                                    t_aligntargets translationoptions;
                                    
                                    
                                    if (classifiers == NULL) {
                                        //ignore classifiers
                                        translationoptions = *reftranslationoptions;
                                    } else {
                                        //are there enough targets for this source to warrant a classifier?
                                        if (alignmodel->alignmatrix[key].size() >= targetthreshold) {
                                            if (debug) cerr << "classifying" << endl;
                                            vector<string> * extrafeatures = classifiers->computeextrafeatures(*line, alignmodel, scorehandling,  key, incontext, *reftranslationoptions, leftcontextsize, rightcontextsize);  
                                            translationoptions = classifiers->classifyfragment(key, incontext, *reftranslationoptions, scorehandling, leftcontextsize, rightcontextsize, changedcount, extrafeatures);
                                            if (extrafeatures != NULL) delete extrafeatures;
                                        } else {
                                            if (debug) cerr << "bypassing classifier, targetthreshold too low" << endl;

                                            for (t_aligntargets::iterator iter = reftranslationoptions->begin(); iter != reftranslationoptions->end(); iter++) {
                                                const EncAnyGram * target = iter->first;
                                                
                                                if (scorehandling == SCOREHANDLING_REPLACE) {
                                                    translationoptions[target].push_back((*reftranslationoptions)[target][0]); //fall back to first statistical value    
                                                } else if (scorehandling == SCOREHANDLING_APPEND) {
                                                    translationoptions[target] = (*reftranslationoptions)[target];
                                                    translationoptions[target].push_back((*reftranslationoptions)[target][0]); //fall back to first statistical value
                                                } else {
                                                    translationoptions[target] = (*reftranslationoptions)[target];
                                                }
                                             }

                                            
                                        }
                                    }
                                    
                                    //write intermediate phrasetable
                                    for (t_aligntargets::iterator iter = translationoptions.begin(); iter != translationoptions.end(); iter++) {
                                        if (iter->second.size() > scorecount) scorecount = iter->second.size();
                                        *TMPTABLE << encodedngram << " ||| " << iter->first->decode(*targetclassdecoder) << " ||| ";
                                        if (debug) cerr << "AFTER CLASSIFICATION: " << encodedngram << " ||| " << iter->first->decode(*targetclassdecoder) << " ||| ";
                                        for (vector<double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                                            if (iter2 != iter->second.begin()) {
                                                *TMPTABLE << " ";
                                                if (debug) cerr << " ";
                                            }
                                            *TMPTABLE << pow(exp(1), *iter2);
                                            if (debug) cerr << pow(exp(1), *iter2);
                                        }
                                        *TMPTABLE << endl;
                                        if (debug) cerr << endl;
                                    }
                                    
                                    //if ((incontext != NULL) && (incontext != key)) delete incontext; //TODO: reenable? segfault
                                }  
                                delete ngram;                  
                                n++;
                            } 
                            
                            
                            
                            
                        }    
                                        
                        
                    }
                }      

                TMPTABLE->close();
                TMPTEST->close(); 
                
                const double changedratio = changedcount / (sentence + 1);
                cerr << "Statistics:" << endl;
                cerr << "\tSentences: " << sentence << endl;
                cerr << "\tSource fragments: " << sourcefragmentcount +unknowncount << endl;
                cerr << "\tUnknown word instances: " << unknowncount << " " << ( (float) changedcount / (sourcefragmentcount+unknowncount)) * 100 << '%' << endl;
                cerr << "\tSource fragments affected by classifier outcome: " << changedcount << " " << ( (float) changedcount / sourcefragmentcount) * 100 << '%' << endl;
                
                stringstream cmd;
                cmd << "makecontextmosesini.py " << outputprefix << " " << scorecount << " > model/contextmoses." << outputprefix << ".ini";          
                cerr << cmd.str() << endl; 
                system(cmd.str().c_str());
            } else {
                cerr << "Classifier already tested (" << outputprefix <<".phrasetable and " << outputprefix << ".txt exist), not overwriting, proceeding with decoding..." << endl;
            }

            /*cerr << "Updating moses configuration..." << endl;
            ifstream *MOSESINI =  new ifstream( "moses.ini" );
            if (!IN->good()) {
            	cerr << "ERROR: Unable to open moses.ini" << endl;
            	exit(5);
            } */       
            if (!skipdecoder) {
                stringstream cmd;
                cmd << "moses -config model/contextmoses." << outputprefix << ".ini < " << outputprefix << ".txt";
                cerr << cmd.str() << endl;  
                system(cmd.str().c_str());
            }

        }  
    }
}
