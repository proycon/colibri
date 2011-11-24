#include <ngram.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>


using namespace std;


void usage() {
    cerr << "Syntax: graphmodel -f encoded-model -c classfile" << endl;        
}

int main( int argc, char *argv[] ) {
    string classfile = "";
    string modelfile = "";
    string outputprefix = "";
    char c;    
    while ((c = getopt(argc, argv, "c:f:ho:")) != -1)
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

    if (outputprefix.empty()) {
        outputprefix = modelfile; //TODO: strip .clsenc. .bin?
    }
    cerr << "Loading class decoder " << classfile << endl;
    ClassDecoder classdecoder = ClassDecoder(classfile);
    
    cerr << "Loading model " << modelfile << endl;
    EncGramModel model = EncGramModel(modelfile);
    
    cerr << "Loaded " << model.types() << " types, " << model.tokens() << " tokens" << endl;
    
    cerr << "Constructing graph " << endl;
    EncGramGraphModel graph = EncGramGraphModel(model);    
    
    cerr << "Saving graph " << endl;
    graph.save(outputprefix + ".graph.bin");
    
    return 0;

}
