
#include <getopt.h>
#include "alignmodel.h"
#include "classifiers.h"

using namespace std;

void usage() {
    cerr << "Training usage: contextmoses -f source-traindatafile [-N|-X|-M] [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Training usage: contextmoses -F testdatafile [-m mosesphrasetable|-d alignmentmodel -S source-class-file -T target-class-file]" << endl;
    cerr << "Classifier types: (pick one)" << endl;
    cerr << " -N           N-Classifier Array, one classifier per pattern size group" << endl;
    cerr << " -X           Construction experts, one classifier per construction" << endl;
    cerr << " -M           Monolithic joined classifier, focus words are joined (-1)" << endl;
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
    
    //cerr << "\t-C number                 Classifier mode" << endl;
    //cerr << "\t   1 - Local context with Classifier Array" << endl;
}

int main( int argc, char *argv[] ) {
    string sourceclassfile = "";
    string targetclassfile = "";    
    string outputprefix="classifier";
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
    
    bool TRAIN = false;
    bool TEST = false;
    
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
    while ((c = getopt_long(argc, argv, "hd:S:T:C:xO:XNc:t:M1a:f:t:l:r:F:",long_options,&option_index)) != -1) {
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
            outputprefix = optarg;
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
        }
    }
    
    if ((!TRAIN && !TEST)) {
        usage();
        exit(2);
    }
    
    if ( (!sourceclassfile.empty()) && (!targetclassfile.empty()) ) {
        cerr << "ERROR: Specify class files (-S -T)" << endl;
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


		
    cerr << "Loading alignment model " << alignmodelfile << endl;
    if (!alignmodelfile.empty()) {
        alignmodel = new AlignmentModel(alignmodelfile,false,true,0, false);
    } else if (!mosesphrasetable.empty()) {
	    cerr << "Loading target class encoder " << targetclassfile << endl;
	    targetclassencoder = new ClassEncoder(targetclassfile);  

        cerr << "Loading source class encoder " << sourceclassfile << endl;
        sourceclassencoder = new ClassEncoder(sourceclassfile);
    
        alignmodel = new AlignmentModel(alignmodelfile, sourceclassencoder, targetclassencoder);    
    } else {
        cerr << "ERROR: No moses phrasetable (-t) or colibri alignment model (-d) specified!" << endl;
        exit(2);
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
		
		
		
		if (mode == CLASSIFIERTYPE_NARRAY) {
    
            cerr << "Building N-Array classifiers" << endl;
            classifiers = new NClassifierArray(outputprefix, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
            
        } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
    
            cerr << "Building construction expert classifiers" << endl;
            classifiers = new ConstructionExperts(outputprefix, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);    
		
		} else if (mode == CLASSIFIERTYPE_MONO) {
		
            if (!singlefocusfeature) {
                    cerr << "ERROR: Monolithic classifier only supported with single focus feature" << endl;
                    exit(2);
            }
            
            cerr << "Building monolithic classifier" << endl;
            classifiers = new MonoClassifier(outputprefix, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);  
            
        }
		
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
                    
            sentence++;

            if (sentence % 10000 == 0) {
                cerr << "\t@" << sentence << endl;
            }
                            
            
            const int l = countwords(linebuffer, linesize);            
            if (l >= 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping!" << endl;
                continue;
            } else if (l == 0) {
            	cerr << "WARNING: Sentence " << sentence << " contains no words, skipping!" << endl;
                continue;
            }
  
                                    
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
                            //match found!
                            const EncAnyGram * incontext = alignmodel->addcontext(&line, (const EncAnyGram * ) ngram, (int) i, leftcontextsize, rightcontextsize);
                            
                            //add to classifier 
                            if (mode == CLASSIFIERTYPE_NARRAY) {                  
                                ((NClassifierArray *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[((const EncAnyGram *) ngram)], leftcontextsize, rightcontextsize);                                                
                            } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
                                ((ConstructionExperts *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[((const EncAnyGram *) ngram)], leftcontextsize, rightcontextsize);                        
                            } else if (mode == CLASSIFIERTYPE_MONO) {
                                ((MonoClassifier *) classifiers)->add((const EncAnyGram*) ngram, incontext, alignmodel->alignmatrix[((const EncAnyGram *) ngram)], leftcontextsize, rightcontextsize);
                            }
                            delete incontext;
                        }  
                        delete ngram;                  
                        n++;
                    } while (found);  
                }
            }
            
            
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
        
        writeclassifierconf(outputprefix, mode, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature);
        cerr << "Training all done" << endl;
    }
    
   
}
