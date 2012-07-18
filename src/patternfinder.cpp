#include <patternmodel.h>
#include <alignmodel.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithms.h>
#include <common.h>

using namespace std;



void usage() {
    cerr << "Syntax: patternfinder -f encoded-corpus" << endl;
    cerr << "Descriptions: Reads an encoded corpus and extracts patterns above a certain threshold. The patterns consist out of n-grams and skip-grams up to a certain maximum value of n. The resulting model will be built in-memory (very memory intestive) and subsequently stored in a binary representation on disk. If a class-file is specified, the output will immediately afterwards be decoded for you as well. Or use the modeldecode program and pass it the binary model." << endl;
    cerr << "Options:" << endl;
    cerr << "\t-c classfile     The classfile to use for decoding. If specified, decoded output will be produced (except in query mode)" << endl;
    cerr << "\t-d modelfile     Load and decode this patternmodel (instead of -f)" << endl;
    cerr << "\t-J modelfile2    Joint decoding of a second (test) model, (using the same classes), use with -d and -c" << endl;
    cerr << "\t-C               Generate a coverage report" << endl;   
    cerr << "\t-t <number>      Token threshold: n-grams and skipgrams occuring less than this will be pruned (default: 2)" << endl;
    cerr << "\t-l <number>      Maximum n-gram/skipgram length (in words, default: 9)" << endl;
    cerr << "\t-s               Compute skip-grams (costs extra memory and time)" << endl;    
    cerr << "\t-T <number>      Skip threshold: only skip content that occurs at least x times will be considered (default: same as -t). Value can never be lower than value for -t" << endl;
    cerr << "\t-u				Create an unindexed model instead of an indexed model. Significantly reduces memory requirements" << endl;
    //cerr << "\t-L               Compute and maintain content of skipgrams (costs extra memory)" << endl;
    cerr << "\t-S <number>      Skip type threshold: only skipgrams with at least x possible types for the skip will be considered, otherwise the skipgram will be pruned  (default: 2, this value is unchangable and fixed to 2 when -u is set)" << endl;
    cerr << "\t-B               Do NOT consider skipgrams that begin with a skip and have no further skips" << endl;
    cerr << "\t-E               Do NOT consider skipgrams that end in a skip and have no further skips" << endl;
    cerr << "\t-Q               Start query mode, allows for pattern lookup against the loaded model. Use with -c and -d to specify the class file and model to load" << endl; 
    cerr << "\t-A <alignmodel>  Prune existing model using source-side of an alignment model (use with -d and -P)" << endl;
    cerr << "\t-P <number>      Prune threshold (use with -A, -d)" << endl;
    cerr << "\t-o <string>      Output prefix" << endl;        
}


void decode(IndexedPatternModel & model, string classfile) {
    cerr << "Loading class decoder " << classfile << endl;
    ClassDecoder classdecoder = ClassDecoder(classfile);
    
    /*const string ngramoutputfile = outputprefix + ".ngrams";
    ofstream *NGRAMSOUT =  new ofstream( ngramoutputfile.c_str() );      
    const string skipgramoutputfile = outputprefix + ".skipgrams";
    ofstream *SKIPGRAMSOUT = NULL;*/
    //if (DOSKIPGRAMS) SKIPGRAMSOUT = new ofstream( skipgramoutputfile.c_str() );      
    cerr << "Decoding" << endl;
    model.decode(classdecoder, (ostream*) &stdout);   
}


int main( int argc, char *argv[] ) {
    
    string classfile = "";
    string corpusfile = "";
    string outputprefix = "";
    string modelfile = "";
    string modelfile2 = "";
    string covviewfile = "";
    string alignmodelfile = "";
    
    
    
    int MINTOKENS = 2;
    int MINSKIPTOKENS = 2;
    unsigned int MINSKIPTYPES = 2;
    int MAXLENGTH = 8;
    bool DOSKIPGRAMS = false;
    bool DOINDEX = true;
    bool DOINITIALONLYSKIP = true;
    bool DOFINALONLYSKIP = true;
    bool DOQUERIER = false;
    bool DOCOVERAGE = false;
    bool DOCOVVIEW = false;
    //bool DOCOMPOSITIONALITY = false;
    bool DEBUG = false;
    double alignthreshold = 0.0;
    char c;    
    while ((c = getopt(argc, argv, "c:f:d:t:T:S:l:o:suLhnBEQDJ:CVA:P:")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
            break;
        case 'd':
            modelfile = optarg;
            break;
        case 'J':
            modelfile2 = optarg;
            break;
        case 'D':
        	DEBUG = true;
        	break;
        case 'C':
            DOCOVERAGE = true;
            break;
        case 'V':
            DOCOVVIEW = true;
            break;
        case 'f':
            corpusfile = optarg;
            break;        
        case 't':
            MINTOKENS = atoi(optarg);
            break;
        case 'T':
            MINSKIPTOKENS = atoi(optarg);            
            break;
        case 'S':
            MINSKIPTYPES = atoi(optarg);            
            break;
        case 'l':
            MAXLENGTH = atoi(optarg);            
            break;
        case 's':
            DOSKIPGRAMS = true;
            break;
        case 'o': 
            outputprefix = optarg;
            break;
        case 'B':
            DOINITIALONLYSKIP = false;
            break;
        case 'E':
            DOFINALONLYSKIP = false;    
            break;
		case 'u':
			DOINDEX = false;    		
			break;
		case 'Q':
			DOQUERIER = true;
			break;
	    case 'A':
	        alignmodelfile = optarg;
	        break;
        case 'P':
            alignthreshold = atof(optarg);
            break;	        
        case 'h':
            usage();
            exit(0);
        case '?':
            if (optopt == 'c') {
                cerr <<  "Option -" << optopt << " requires an argument." << endl;
            } else {
                cerr << "Unknown option: -" <<  optopt << endl;
            }
            
            return 1;
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
    
    if (DOQUERIER && classfile.empty()) {
            cerr << "ERROR: To use the query mode (-Q), specify a -c classfile and an existing model (-d)" << endl;
            usage();
            exit(2);    	
    }
    
    if (DOCOVVIEW && (classfile.empty() || modelfile.empty() || corpusfile.empty() )) {
        cerr << "ERROR: When generating a coverage view, you need to specify a class file, model file and corpus file (-c, -d, -f)" << endl;
        usage();
        exit(2);
    }
    /*if (corpusfile.empty() ) {
        if (modelfile.empty() || classfile.empty()) {
            cerr << "ERROR: Need to specify -f corpusfile to compute pattern, or -d modelfile -c classfile to decode an existing model" << endl;
            usage();
            exit(2);
        }
    }*/
    
    if (outputprefix.empty()) {
        outputprefix = corpusfile; //TODO: strip .clsenc. .bin?
        strip_extension(outputprefix,"colibri");
        strip_extension(outputprefix,"bin");
        strip_extension(outputprefix,"clsenc");
        strip_extension(outputprefix,"txt");    
    }

    
    if ((!corpusfile.empty()) && (!modelfile.empty())) {
    
        if (DOCOVVIEW) {
            IndexedPatternModel model = IndexedPatternModel(modelfile, DEBUG);
            ClassDecoder classdecoder = ClassDecoder(classfile);
            model.coveragereport((ostream*) &cerr, corpusfile, (ostream*) &cout, &classdecoder);      
        } else {

            if (DOINDEX) {
                //not implemented yet
                cerr << "Loading reference model" << endl;
		        IndexedPatternModel refmodel = IndexedPatternModel(modelfile, DEBUG);
		        
		        cerr << "Computing model on " << corpusfile << endl;
		        IndexedPatternModel model = IndexedPatternModel(corpusfile, refmodel, MAXLENGTH, MINTOKENS, DOSKIPGRAMS, MINSKIPTOKENS, MINSKIPTYPES, DOINITIALONLYSKIP,DOFINALONLYSKIP);

		        cerr << "Saving " << outputprefix << ".indexedpatternmodel.colibri"  << endl;
		        const string outputfile = outputprefix + ".indexedpatternmodel.colibri";
		        model.save(outputfile);     
            } else {
                cerr << "Loading reference model" << endl;
		        UnindexedPatternModel refmodel = UnindexedPatternModel(modelfile, DEBUG);
		        
		        cerr << "Computing model on " << corpusfile << endl;
		        UnindexedPatternModel model = UnindexedPatternModel(corpusfile, refmodel, MAXLENGTH, MINTOKENS, DOSKIPGRAMS, MINSKIPTOKENS, MINSKIPTYPES, DOINITIALONLYSKIP,DOFINALONLYSKIP);

		        cerr << "Saving " << outputprefix << ".unindexedpatternmodel.colibri"  << endl;
		        const string outputfile = outputprefix + ".unindexedpatternmodel.colibri";
		        model.save(outputfile);       		  
            }
        }
        
    } else if (!corpusfile.empty()) {
    	if (DOINDEX) {
    
		    cerr << "Computing model on " << corpusfile << endl;
		    IndexedPatternModel model = IndexedPatternModel(corpusfile, MAXLENGTH, MINTOKENS, DOSKIPGRAMS, MINSKIPTOKENS, MINSKIPTYPES, DOINITIALONLYSKIP,DOFINALONLYSKIP);
		        
		    cerr << "Saving " << outputprefix << ".indexedpatternmodel.colibri"  << endl;
		    const string outputfile = outputprefix + ".indexedpatternmodel.colibri";
		    model.save(outputfile);            
		    
		    
		    if (!classfile.empty()) {
		        cerr << "Loading class decoder " << classfile << endl;
		        ClassDecoder classdecoder = ClassDecoder(classfile);
		        if (DOQUERIER) {
		        	cerr << "Loading class encoder " << classfile << endl;
		        	ClassEncoder classencoder = ClassEncoder(classfile);
		        	cerr << "Starting query mode:" << endl;
		        	model.querier(classencoder, classdecoder);
		        } else {
		        	cerr << "Decoding" << endl;
		        	model.decode(classdecoder, (ostream*) &cout);
		        }   
		    }
		    
        } else {    
		    cerr << "Computing model on " << corpusfile << endl;
		    UnindexedPatternModel model = UnindexedPatternModel(corpusfile, MAXLENGTH, MINTOKENS, DOSKIPGRAMS, MINSKIPTOKENS ,DOINITIALONLYSKIP,DOFINALONLYSKIP);
		        
		    cerr << "Saving " << outputprefix << ".unindexedpatternmodel.colibri"  << endl;
		    const string outputfile = outputprefix + ".unindexedpatternmodel.colibri";
		    model.save(outputfile);            
		    
		    
		    if (!classfile.empty()) {
		        cerr << "Loading class decoder " << classfile << endl;
		        ClassDecoder classdecoder = ClassDecoder(classfile);
		        if (DOQUERIER) {
		        	cerr << "Loading class encoder " << classfile << endl;
		        	ClassEncoder classencoder = ClassEncoder(classfile);
		        	cerr << "Starting query mode:" << endl;
		        	model.querier(classencoder, classdecoder);
		        } else {
		        	cerr << "Decoding" << endl;
			        model.decode(classdecoder, (ostream*) &cout);
			    }   
		    }
		            	
        }
    } else if ( (!modelfile.empty()) && (!alignmodelfile.empty()) && (alignthreshold > 0)  ) {
        if (DOINDEX) {
    	    cerr << "Loading model" << endl;
		    IndexedPatternModel model = IndexedPatternModel(modelfile, DEBUG);
		    cerr << "Loading Alignment Model" << endl;
		    AlignmentModel alignmodel = AlignmentModel(alignmodelfile);
		    unsigned int pruned = model.prunebyalignment(alignmodel.alignmatrix, alignthreshold);		    
            cerr << "pruned " << pruned << endl;
		    
		    if (!outputprefix.empty()) {
		        const string outputfile = outputprefix + ".indexedpatternmodel.colibri";
		        cerr << "Saving " << outputprefix << ".indexedpatternmodel.colibri"  << endl;
		        model.save(outputfile);
		    } else {
		        cerr << "Saving " << modelfile << endl;
		        model.save(modelfile);		    
		    }                 
        } else {
            cerr << "Not implemented yet for unindexed models" << endl;
            exit(2);
        }       
    } else if ( (!modelfile.empty()) && ((!classfile.empty()) || DOCOVERAGE  ) ) {
    	if (DOINDEX) {
    	    cerr << "Loading model" << endl;
		    IndexedPatternModel model = IndexedPatternModel(modelfile, DEBUG);
		    if (!classfile.empty()) {
		        cerr << "Loading class decoder " << classfile << endl;
		        ClassDecoder classdecoder = ClassDecoder(classfile);
		        if (DOQUERIER) {
		        	cerr << "Loading class encoder " << classfile << endl;
		        	ClassEncoder classencoder = ClassEncoder(classfile);
		        	cerr << "Starting query mode:" << endl;
		        	model.querier(classencoder, classdecoder);
		        } else {
		            if (modelfile2.empty()) {		        
    		        	cerr << "Decoding" << endl;
		        	    if (covviewfile.empty()) model.decode(classdecoder, (ostream*) &cout);
    			     } else {    			
    			        cerr << "Loading test model" << endl;     
    			        IndexedPatternModel testmodel = IndexedPatternModel(modelfile2, DEBUG);
    			        cerr << "Joint decoding" << endl;
    			        if (covviewfile.empty()) model.decode(testmodel, classdecoder, (ostream*) &cout);
    			     }		        	    
		        }
		    }
	       	if (DOCOVERAGE) {
   		        model.coveragereport((ostream*) &cout);		        
	        }

		    
		} else {
		    cerr << "Loading model" << endl;
		    UnindexedPatternModel model = UnindexedPatternModel(modelfile, DEBUG);
		    if (!classfile.empty()) {
		        cerr << "Loading class decoder " << classfile << endl;
		        ClassDecoder classdecoder = ClassDecoder(classfile);
				if (DOQUERIER) {
		        	cerr << "Loading class encoder " << classfile << endl;
		        	ClassEncoder classencoder = ClassEncoder(classfile);
		        	cerr << "Starting query mode:" << endl;
		        	model.querier(classencoder, classdecoder);
		        } else {
		            if (modelfile2.empty()) {
			            cerr << "Decoding" << endl;
    			        model.decode(classdecoder, (ostream*) &cout);
    			     } else {    			
    			        cerr << "Loading test model" << endl;     
    			        UnindexedPatternModel testmodel = UnindexedPatternModel(modelfile2, DEBUG);
    			        cerr << "Joint decoding" << endl;
    			        model.decode(testmodel, classdecoder, (ostream*) &cout);
    			     }
			    }   
		    }			
		}
        
    } else {
        cerr << "Nothing to do?" << endl;
    }


}
