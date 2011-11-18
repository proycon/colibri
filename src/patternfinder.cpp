#include <ngram.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <limits>
#include <iomanip> // contains setprecision()
#include <map>
#include <unistd.h>
#include <cmath>
#include <set>
#include <algorithms.h>

using namespace std;



void usage() {
    cerr << "Syntax: patternfinder -f encoded-corpus" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-c classfile     The classfile to use for decoding. If specified, decoded output will be produced" << endl;
    cerr << "\t-t <number>      Token threshold: n-grams and skipgrams occuring less than this will be pruned (default: 2)" << endl;
    cerr << "\t-l <number>      Maximum n-gram/skipgram length (in words, default: 9)" << endl;
    cerr << "\t-s               Compute skip-grams" << endl;    
    cerr << "\t-T <number>      Skip threshold: only skip content that occurs at least x times will be considered (default: 2) " << endl;
    cerr << "\t-i               Compute and output index" << endl;
    cerr << "\t-L               Compute and output content of skipgrams" << endl;
    cerr << "\t-S <number>      Skip type threshold: only skipgrams with x possible types for the skip will be considered, otherwise the skipgram will be pruned  (default: 2, works only with -L enabled)" << endl;
    cerr << "\t-B               Do NOT consider skipgrams that begin with a skip and have no further skips" << endl;
    cerr << "\t-E               Do NOT consider skipgrams that end in a skip and have no further skips" << endl;
    cerr << "\t-o <string>      Output prefix" << endl;
        
}


int main( int argc, char *argv[] ) {
    
    string classfile = "";
    string corpusfile = "";
    string outputprefix = "";
    
    int MINTOKENS = 2;
    int MINSKIPTOKENS = 2;
    unsigned int MINSKIPTYPES = 2;
    int MAXLENGTH = 8;
    bool DOSKIPGRAMS = false;
    bool DOINDEX = false;
    bool DOSKIPCONTENT = false;
    bool DOINITIALONLYSKIP = true;
    bool DOFINALONLYSKIP = true;
    //bool DOCOMPOSITIONALITY = false;
    
    char c;    
    while ((c = getopt(argc, argv, "c:f:t:T:S:l:o:siLhnBE")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
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
        case 'i':
            DOINDEX = true;
            break;
        case 'L':
            DOSKIPCONTENT = true;
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
    
    if (corpusfile.empty()) {
        cerr << "ERROR: Need to specify -f corpusfile" << endl;
        usage();
        exit(2);
    }
    
    if (outputprefix.empty()) {
        outputprefix = corpusfile;
    }

    
    
    cerr << "Computing model on " << corpusfile << endl;
    EncGramModel model = EncGramModel(corpusfile, MAXLENGTH, MINTOKENS, DOSKIPGRAMS, MINSKIPTOKENS, MINSKIPTYPES, DOINDEX, DOSKIPCONTENT,DOINITIALONLYSKIP,DOFINALONLYSKIP);
        
    cerr << "Saving "  << endl;
    const string outputfile = outputprefix + ".bin";    
    model.save(outputfile);
    
    
    if (!classfile.empty()) {
        cerr << "Loading class decoder " << classfile << endl;
        ClassDecoder classdecoder = ClassDecoder(classfile);
        
        const string ngramoutputfile = outputprefix + ".ngrams";
        ofstream *NGRAMSOUT =  new ofstream( ngramoutputfile.c_str() );      
        const string skipgramoutputfile = outputprefix + ".skipgrams";
        ofstream *SKIPGRAMSOUT = NULL;
        if (DOSKIPGRAMS) SKIPGRAMSOUT = new ofstream( skipgramoutputfile.c_str() );      
        cerr << "Decoding" << endl;
        model.decode(classdecoder, (ostream*) NGRAMSOUT, (ostream*) SKIPGRAMSOUT);    
    }    

}
