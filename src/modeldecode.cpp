#include <ngram.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;


void usage() {
    cerr << "Syntax: modeldecode -f encoded-model -c classfile" << endl;        
    cerr << "Options: -o <string>      Output prefix" << endl;
}

int main( int argc, char *argv[] ) {
    string classfile = "";
    string modelfile = "";
    string outputprefix = "";
    char c;    
    while ((c = getopt(argc, argv, "c:f:h")) != -1)
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
    
    cerr << "Decoding " << modelfile << endl;
    const string ngramoutputfile = outputprefix + ".ngrams";
    ofstream *NGRAMSOUT =  new ofstream( ngramoutputfile.c_str() );      
    const string skipgramoutputfile = outputprefix + ".skipgrams";
    ofstream *SKIPGRAMSOUT = NULL;
    SKIPGRAMSOUT = new ofstream( skipgramoutputfile.c_str() );      
    model.decode(classdecoder, (ostream*) NGRAMSOUT, (ostream*) SKIPGRAMSOUT);    
    
    
    return 0;
}
