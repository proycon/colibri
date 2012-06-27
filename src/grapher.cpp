#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>
#include <patternmodel.h>
#include <common.h>


using namespace std;


void usage() {
    cerr << "Syntax: graphmodel " << endl;        
    cerr << "Constructs a graph model" << endl;
    cerr << "\t-f filename.indexedpatternmodel.colibri		Indexed pattern model to load (required for building a graph model)" << endl;      
    cerr << "\t-P               Compute/load subsumption relations from children to parents (reverse of -C)" << endl;
    cerr << "\t-C               Compute/load subsumption relations from parents to children (reverse of -P)" << endl;
    cerr << "\t-S               Compute/load subsumption skipgram to skipcontent relations" << endl;
    cerr << "\t-s               Compute/load subsumption skip-content to skipgram relations (reverse of -S)" << endl;
    cerr << "\t-L               Compute/load subsumption predecessor relations (constructions to the left)" << endl;
    cerr << "\t-R               Compute/load subsumption sucessor relations (constructions to the right)" << endl;
    cerr << "\t-T               Compute/load subsumption template relations" << endl;
    cerr << "\t-I               Compute/load subsumption instance relations (reverse of -T)" << endl;
    cerr << "\t-a               Compute/load subsumption all relations" << endl;      
    cerr << "\t-X               Compute/load exclusive count" << endl;
    cerr << "\t------------------------------------------------------------------------------" << endl;
    cerr << "\t-r               Keep only transitive reduction (sizes down the model)" << endl;
    cerr << "\t-d filename.graphpatternmodel.colibri		Graph pattern model to load (for decoding an existing model, use with -c)" << endl;
    cerr << "\t-c classfile     The classfile to use for decoding. If specified, decoded output will be produced (use with -d)" << endl;
	cerr << "\t-g               Output relations" << endl;    
    cerr << "\t-G               Output graphviz graph for visualisation" << endl;
    cerr << "\t-q word          Query word (use with -G to output a selected graph)" << endl;
}

int main( int argc, char *argv[] ) {
    string classfile = "";
    string patternmodelfile = "";
    string modelfile = "";
    string outputprefix = "";
    string querystring = "";
    
    bool DOPARENTS = false;
    bool DOCHILDREN = false;
    bool DOXCOUNT = false;
    bool DOSUCCESSORS = false;
    bool DOPREDECESSORS = false;
    bool DOSKIPCONTENT = false;
    bool DOSKIPUSAGE = false;
    bool DOOUTPUTRELATIONS = false;
    bool TRANSITIVEREDUCTION = false;
    
    bool DOTEMPLATES = false; 
    bool DOINSTANCES = false;
    
    bool DOGRAPHVIZ = false; 
    bool DEBUG = false;
    
    char c;    
    while ((c = getopt(argc, argv, "ad:c:f:ho:PCXrGq:LRSsgITDO")) != -1)
        switch (c)
        {
        case 'a':
        	DOTEMPLATES = true;
        	DOINSTANCES = true;
        	DOPARENTS = true;
        	DOCHILDREN = true;
        	DOXCOUNT = true;
        	DOSUCCESSORS = true;
        	DOPREDECESSORS = true;
        	DOSKIPCONTENT = true;
        	DOSKIPUSAGE = true;
        	break;
        case 'c':
            classfile = optarg;
            break;
        case 'f':
            patternmodelfile = optarg;
            break;
        case 'd':
            modelfile = optarg;
            break;
        case 'D':
            DEBUG = true;
            break;               
        case 'o': 
            outputprefix = optarg;
            break;
        case 'P': 
            DOPARENTS = true;
            break;
        case 'C': 
            DOCHILDREN = true;
            break;
        case 'R':
        	DOSUCCESSORS = true;
        	break;
		case 'L':
        	DOPREDECESSORS = true;
        	break;
		case 'S':
        	DOSKIPCONTENT = true;
        	break;       
		case 's':
        	DOSKIPUSAGE = true;
        	break;     
        case 'I':
        	DOINSTANCES = true;
        	break;
        case 'T':
        	DOTEMPLATES = true;
        	break;  
        case 'X': 
            DOXCOUNT = true;
            break;
        case 'G':
        	DOGRAPHVIZ = true;
        	break;
        case 'g':
        	DOOUTPUTRELATIONS = true;
        	break;
        case 'r': 
        	TRANSITIVEREDUCTION = true;
        	break;
        case 'h':
            usage();
            break;
        case 'q':
        	querystring = optarg;
        	break; 
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            usage();
            exit(2);
        }
    

    
    if (patternmodelfile.empty() && modelfile.empty()) {
            cerr << "ERROR: Need to specify -f corpusfile to compute pattern, or -d graphmodelfile -c classfile to decode an existing model" << endl;
            usage();
            exit(2);
    }

    /*if ((!DOPARENTS) && (!DOCHILDREN) && (!DOXCOUNT)) {
        cerr << "No options selected, no relations to load or construct" << endl;
        usage();
        exit(2);
    }*/

    
    

    
    
    
    if (modelfile.empty()) {
        if (outputprefix.empty()) {
            outputprefix = patternmodelfile; //TODO: strip .clsenc. .bin?
            strip_extension(outputprefix, string("bin"));
            strip_extension(outputprefix, string("colibri"));
            strip_extension(outputprefix, string("indexedpatternmodel"));
            strip_extension(outputprefix, string("clsenc"));
            strip_extension(outputprefix, string("txt"));   
        }
        
        cerr << "Loading pattern model " << patternmodelfile << endl;
        IndexedPatternModel patternmodel = IndexedPatternModel(patternmodelfile, DEBUG);
    
        cerr << "Loaded " << patternmodel.types() << " types, " << patternmodel.tokens() << " tokens" << endl;
            
        cerr << "Constructing graph " << endl;
        GraphPatternModel graphmodel = GraphPatternModel(&patternmodel, DOPARENTS, DOCHILDREN, DOXCOUNT, DOTEMPLATES, DOINSTANCES,DOSKIPUSAGE,DOSKIPCONTENT,DOSUCCESSORS,DOPREDECESSORS);
        
        
        
        if (TRANSITIVEREDUCTION) {
        	cerr << "Pruning and keeping only transitive reduction " << endl;
        	const int pruned = graphmodel.transitivereduction();
        	cerr << "\tPruned " << pruned << " edges" << endl;
        }
        
        graphmodel.stats( (ostream*) &cerr );
        
        cerr << "Saving graph " << outputprefix << ".graphpatternmodel.colibri" << endl;
        graphmodel.save(outputprefix + ".graphpatternmodel.colibri");
        
        if (!classfile.empty()) {
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
           
            cerr << "Decoding graph" << endl;
            if (DOGRAPHVIZ) {
            	graphmodel.outputgraph(classdecoder, (ostream*) &cout );
            } else {
            	graphmodel.decode(classdecoder, (ostream*) &cout, DOOUTPUTRELATIONS);
            }
            
        }        
    } else {
        if (!classfile.empty()) {           
            cerr << "Loading graph model " << modelfile << endl;
            GraphPatternModel graphmodel = GraphPatternModel(modelfile, DOPARENTS, DOCHILDREN, DOXCOUNT, DOTEMPLATES, DOINSTANCES,DOSKIPUSAGE,DOSKIPCONTENT,DOSUCCESSORS,DOPREDECESSORS, DEBUG);  
            
            graphmodel.stats( (ostream*) &cerr );
                                         
        
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
            cerr << "Decoding graph" << endl;
            if (DOGRAPHVIZ) {
            	if (!querystring.empty()) {
            		cerr << "Loading class encoder " << classfile << endl;
            		ClassEncoder encoder = ClassEncoder(classfile);
            		
            		
            		const EncAnyGram * anygram = encoder.input2anygram(querystring, true);
            		cerr << "Outputting graph for \"" << anygram->decode(classdecoder) << "\"" << endl;
            		//unsigned char buffer[65536];
            		//char buffersize = encoder.encodestring(querystring, buffer);
            		//EncNGram ngram = EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte
            		graphmodel.outputgraph(classdecoder,(ostream*) &cout, anygram); // (EncAnyGram*) &ngram);
            	} else {            	
            		graphmodel.outputgraph(classdecoder, (ostream*) &cout );
            	}
            } else if (DOOUTPUTRELATIONS) {
            	if (!querystring.empty()) {
            		cerr << "Loading class encoder " << classfile << endl;
            		ClassEncoder encoder = ClassEncoder(classfile);
            		
            		
            		const EncAnyGram * anygram = encoder.input2anygram(querystring, true);
            		cerr << "Outputting graph for \"" << anygram->decode(classdecoder) << "\"" << endl;
            		graphmodel.outputrelations(classdecoder,(ostream*) &cout, anygram); //(EncAnyGram*) &ngram);
            	} else {            	
            		cerr << "Specify a query string with -q" << endl;
            	}             
            } else {
            	graphmodel.decode(classdecoder, (ostream*) &cout);
            }
        } else {
            cerr << "No classer specified" << endl;
            exit(4);
        }
    }
    return 0;
}
