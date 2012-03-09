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
    cerr << "\t-X               Compute/load exclusive count" << endl;
    cerr << "\t------------------------------------------------------------------------------" << endl;
    cerr << "\t-r               Keep only transitive reduction (sizes down the model)" << endl;
    cerr << "\t-d filename.graphpatternmodel.colibri		Graph pattern model to load (for decoding an existing model, use with -c)" << endl;
    cerr << "\t-c classfile     The classfile to use for decoding. If specified, decoded output will be produced (use with -d)" << endl;
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
    bool TRANSITIVEREDUCTION = false;
    
    bool DOGRAPHVIZ = false; 
    
    char c;    
    while ((c = getopt(argc, argv, "d:c:f:ho:PCXrGq:")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
            break;
        case 'f':
            patternmodelfile = optarg;
            break;
        case 'd':
            modelfile = optarg;
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
        case 'X': 
            DOXCOUNT = true;
            break;
        case 'G':
        	DOGRAPHVIZ = true;
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

    if ((!DOPARENTS) && (!DOCHILDREN) && (!DOXCOUNT)) {
        cerr << "No options selected, no relations to load or construct" << endl;
        usage();
        exit(2);
    }

    
    

    
    
    
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
        IndexedPatternModel patternmodel = IndexedPatternModel(patternmodelfile);
    
        cerr << "Loaded " << patternmodel.types() << " types, " << patternmodel.tokens() << " tokens" << endl;
            
        cerr << "Constructing graph " << endl;
        GraphPatternModel graphmodel = GraphPatternModel(&patternmodel, DOPARENTS, DOCHILDREN, DOXCOUNT);
        
        if (TRANSITIVEREDUCTION) {
        	cerr << "Pruning and keeping only transitive reduction " << endl;
        	const int pruned = graphmodel.transitivereduction();
        	cerr << "\tPruned " << pruned << " edges" << endl;
        }
        
        cerr << "Saving graph " << outputprefix << ".graphpatternmodel.colibri" << endl;
        graphmodel.save(outputprefix + ".graphpatternmodel.colibri");
        
        if (!classfile.empty()) {
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
           
            cerr << "Decoding graph" << endl;
            if (DOGRAPHVIZ) {
            	graphmodel.outputgraph(classdecoder, (ostream*) &cout );
            } else {
            	graphmodel.decode(classdecoder, (ostream*) &cout, (ostream*) &cout);
            }
            
        }        
    } else {
        if (!classfile.empty()) {           
            cerr << "Loading graph model " << modelfile << endl;
            GraphPatternModel graphmodel = GraphPatternModel(modelfile);  
            
            cerr << "  Loaded " << graphmodel.model->types() << " types, " << graphmodel.model->tokens() << " tokens" << endl;
            cerr << "  Parent relations available for " << graphmodel.rel_subsumption_parents.size() << " patterns" << endl;
            cerr << "  Child relation available for " << graphmodel.rel_subsumption_children.size() << " patterns" << endl;      
        
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
            cerr << "Decoding graph" << endl;
            if (DOGRAPHVIZ) {
            	if (!querystring.empty()) {
            		cerr << "Loading class encoder " << classfile << endl;
            		ClassEncoder encoder = ClassEncoder(classfile);
            		
            		cerr << "Outputting graph for \"" << querystring << "\"" << endl;
            		unsigned char buffer[65536];
            		char buffersize = encoder.encodestring(querystring, buffer);
            		EncNGram ngram = EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte
            		graphmodel.outputgraph(classdecoder,(ostream*) &cout, (EncAnyGram*) &ngram);
            	} else {            	
            		graphmodel.outputgraph(classdecoder, (ostream*) &cout );
            	}
            } else {
            	graphmodel.decode(classdecoder, (ostream*) &cout, (ostream*) &cout);
            }
        } else {
            cerr << "No classer specified" << endl;
            exit(4);
        }
    }
    return 0;
}
