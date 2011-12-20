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
    cerr << "\t-P       Compute/load subsumption relations from children to parents (reverse of -C)" << endl;
    cerr << "\t-C       Compute/load subsumption relations from parents to children (reverse of -P)" << endl;
    cerr << "\t-X       Compute/load exclusive count" << endl;
    
}

int main( int argc, char *argv[] ) {
    string classfile = "";
    string modelfile = "";
    string outputprefix = "";
    bool DOPARENTS = false;
    bool DOCHILDREN = false;
    bool DOXCOUNT = false;
    
    char c;    
    while ((c = getopt(argc, argv, "c:f:ho:PCX")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
            break;
        case 'f':
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
    

    
    if ((modelfile.empty()) || (classfile.empty())) {
        usage();
        exit(2);
    }

    if ((!DOPARENTS) && (!DOCHILDREN) && (!DOXCOUNT)) {
        usage();
        exit(2);
    }

    if (outputprefix.empty()) {
        outputprefix = modelfile; //TODO: strip .clsenc. .bin?
        strip_extension(outputprefix, string("bin"));
        strip_extension(outputprefix, string("colibri"));
        strip_extension(outputprefix, string("indexedpatternmodel"));
        strip_extension(outputprefix, string("clsenc"));
        strip_extension(outputprefix, string("txt"));   
    }
    cerr << "Loading class decoder " << classfile << endl;
    ClassDecoder classdecoder = ClassDecoder(classfile);
    
    cerr << "Loading model " << modelfile << endl;
    IndexedPatternModel model = IndexedPatternModel(modelfile);
    
    cerr << "Loaded " << model.types() << " types, " << model.tokens() << " tokens" << endl;
    
    cerr << "Constructing graph " << endl;
    GraphPatternModel graph = GraphPatternModel(&model, DOPARENTS, DOCHILDREN, DOXCOUNT);
    
    cerr << "Saving graph " << endl;
    graph.save(outputprefix + ".graphpatternmodel.colibri");
    
    return 0;
}
