#include "alignmodel.h"

using namespace std;

void usage() {
    cerr << "Aligner: aligner -C -s source-model -S source-class-file -t target-model -T target-class-file" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-s sourcemodelfile       Source graph model file (*.graphmodel.colibri)" << endl;
    cerr << "\t-S sourceclassfile       Source class file (for decoding)" << endl;
    cerr << "\t-t targetmodelfile       Target model file (*.graphmodel.colibri)"  << endl;
    cerr << "\t-T targetclassfile       Target class file (for decoding)" << endl;
    cerr << "\t-p pruning-threshold     Prune all alignments with a lower score than specified (0 <= x <= 1)" << endl;
    cerr << "\t-C                       Use Jaccard co-occurrence method" << endl;
    cerr << "\t-d                       Decode results" << endl;
}

int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    double prunevalue = 0.8;
    bool DOCOOC = true;
    bool DODECODE = false;
    char c;    
    while ((c = getopt(argc, argv, "s:S:t:T:p:dC")) != -1)
        switch (c)
        {
        case 'p':
            prunevalue = atof(optarg);
            cerr << "Prune value set to: " << prunevalue << endl; 
            break;
        case 'd':
            DODECODE = true;
            break;
        case 'C':
            DOCOOC = true;
            break;
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
    
    if (!DOCOOC) {
    	cerr << "Error: No alignment method selected (add -C for Jaccard-based co-occurence)" << endl;
    	usage();
    	exit(3);
    }
    
    
    cerr << "Loading source class decoder " << sourceclassfile << endl;
    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
    cerr << "Loading target class decoder " << targetclassfile << endl;
    ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    
    
    cerr << "Loading source model " << sourcemodelfile << endl;
    DoubleIndexedGraphPatternModel sourcemodel = DoubleIndexedGraphPatternModel(sourcemodelfile);
    cerr << "  Loaded " << sourcemodel.types() << " types, " << sourcemodel.tokens() << " tokens" << endl;
    cerr << "  Reverse index has " << sourcemodel.reverseindex.size() << " sentences" << endl;    
    
    cerr << "Loading target model " << targetmodelfile << endl;
    DoubleIndexedGraphPatternModel targetmodel = DoubleIndexedGraphPatternModel(targetmodelfile);
    cerr << "  Loaded " << targetmodel.types() << " types, " << targetmodel.tokens() << " tokens" << endl;
    cerr << "  Reverse index has " << targetmodel.reverseindex.size() << " sentences" << endl;
    
    if (DOCOOC) {
		cerr << "Computing alignment model..." << endl;
		CoocAlignmentModel alignmodel = CoocAlignmentModel(sourcemodel,targetmodel, prunevalue);    
		if (DODECODE) {
		    cerr << "Decoding..." << endl;
		    alignmodel.decode(sourceclassdecoder, targetclassdecoder, &cout);    
		}
	}	
	
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
