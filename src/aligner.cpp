#include "ngram.h"
#include <unordered_map>

using namespace std;

void usage() {
    cerr << "Aligner: aligner -s source-model -S source-class-file -t target-model -T target-class-file" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-s sourcemodelfile       Source model file" << endl;
    cerr << "\t-S sourceclassfile       Source class file (for decoding)" << endl;
    cerr << "\t-t targetmodelfile       Target model file" << endl;
    cerr << "\t-T targetclassfile       Target class file (for decoding)" << endl;
}

int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    char c;    
    while ((c = getopt(argc, argv, "s:S:t:T:")) != -1)
        switch (c)
        {
        case 's':
            sourcemodelfile = optarg;
            break;
        case 't':
            targetmodelfile = optarg;
            break;
        case 'S':
            sourceclassfile = optarg;
            break;
        case 'T':
            targetclassfile = optarg;
            break;        
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
    if (sourcemodelfile.empty() || sourceclassfile.empty() || targetmodelfile.empty() || targetclassfile.empty()) {
        usage();
        exit(2);
    }
    
    
    cerr << "Loading source class decoder " << sourceclassfile << endl;
    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
    cerr << "Loading target class decoder " << targetclassfile << endl;
    ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    
    
    cerr << "Loading source model " << sourcemodelfile << endl;
    EncGramModel sourcemodel = EncGramModel(sourcemodelfile,true,false,true);
    cerr << "  Loaded " << sourcemodel.types() << " types, " << sourcemodel.tokens() << " tokens" << endl;
    //cerr << "  Reverse index has " << sourcemodel.reverse_index_size() << " sentences" << endl;    
    
    cerr << "Loading target model " << targetmodelfile << endl;
    EncGramModel targetmodel = EncGramModel(targetmodelfile,true,false,true);
    cerr << "  Loaded " << targetmodel.types() << " types, " << targetmodel.tokens() << " tokens" << endl;
    //cerr << "  Reverse index has " << targetmodel.reverse_index_size() << " sentences" << endl;
    
    cerr << "Computing alignment model..." << endl;
    CoocAlignmentModel(sourcemodel,targetmodel);    
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
