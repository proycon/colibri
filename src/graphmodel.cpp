#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>
#include <patternmodel.h>
#include <common.h>


using namespace std;


void usage() {
    cerr << "Syntax: graphmodel -f filename.indexedpatternmodel.colibri" << endl;        
    cerr << "Constructs a graph model" << endl;    
    cerr << "\t-P               Compute/load subsumption relations from children to parents (reverse of -C)" << endl;
    cerr << "\t-C               Compute/load subsumption relations from parents to children (reverse of -P)" << endl;
    cerr << "\t-X               Compute/load exclusive count" << endl;
    cerr << "\t-c classfile     The classfile to use for decoding. If specified, decoded output will be produced" << endl;
    
}

int main( int argc, char *argv[] ) {
    string classfile = "";
    string patternmodelfile = "";
    string modelfile = "";
    string outputprefix = "";
    bool DOPARENTS = false;
    bool DOCHILDREN = false;
    bool DOXCOUNT = false;
    
    char c;    
    while ((c = getopt(argc, argv, "d:c:f:ho:PCX")) != -1)
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
        case 'h':
            usage();
            break;
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            usage();
            exit(2);
        }
    

    
    if (patternmodelfile.empty()) {
            cerr << "ERROR: Need to specify -f corpusfile to compute pattern, or -f patternmodelfile -d graphmodelfile -c classfile to decode an existing model" << endl;
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
        
        cerr << "Saving graph " << outputprefix << ".graphpatternmodel.colibri" << endl;
        graphmodel.save(outputprefix + ".graphpatternmodel.colibri");
        
        if (!classfile.empty()) {
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
            cerr << "Decoding graph" << endl;
            graphmodel.decode(classdecoder, (ostream*) &stdout, (ostream*) &stdout);
        }        
    } else {
        cerr << "Loading graph model " << modelfile << endl;
        GraphPatternModel graphmodel = GraphPatternModel(modelfile, patternmodelfile);
        
        if (!classfile.empty()) {
            cerr << "Loading class decoder " << classfile << endl;
            ClassDecoder classdecoder = ClassDecoder(classfile);
            
            cerr << "Decoding graph" << endl;
            graphmodel.decode(classdecoder, (ostream*) &stdout, (ostream*) &stdout);
        }
    }
    return 0;
}
