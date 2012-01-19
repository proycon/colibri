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
    cerr << "\t-o occurence-threshold   Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold        Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold      Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio          Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
}

int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    double prunevalue = 0.8;
    bool DOCOOC = true;
    bool DODECODE = false;
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    
    char c;    
    while ((c = getopt(argc, argv, "s:S:t:T:p:dCo:F:x:X:")) != -1)
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
        case 'o':
            COUNTTHRESHOLD = atoi(optarg);
            break;
		case 'x':
            XCOUNTTHRESHOLD = atoi(optarg);
            break;
		case 'F':
            FREQTHRESHOLD = atof(optarg);
            break;
        case 'X':
            XCOUNTRATIOTHRESHOLD = atof(optarg);
            break;
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
    if (sourcemodelfile.empty()  || targetmodelfile.empty()) {
  	    cerr << "Error: Specify at least a source model, target model, and alignment method!" << endl;
        usage();
        exit(2);
    }
    
    if (!DOCOOC) {
    	cerr << "Error: No alignment method selected (add -C for Jaccard-based co-occurrence)" << endl;
    	usage();
    	exit(3);
    }
    

    
    cerr << "Loading source model " << sourcemodelfile << endl;
    SelectivePatternModel sourcemodel = SelectivePatternModel(sourcemodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD);
    cerr << "  Loaded " << sourcemodel.types() << " types, " << sourcemodel.tokens() << " tokens" << endl;
    cerr << "  Reverse index has " << sourcemodel.reverseindex.size() << " sentences" << endl;    
    
    cerr << "Loading target model " << targetmodelfile << endl;
    SelectivePatternModel targetmodel = SelectivePatternModel(targetmodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD);
    cerr << "  Loaded " << targetmodel.types() << " types, " << targetmodel.tokens() << " tokens" << endl;
    cerr << "  Reverse index has " << targetmodel.reverseindex.size() << " sentences" << endl;
    
    
    AlignmentModel * alignmodel = NULL;
    
    if (DOCOOC) {
		cerr << "Computing alignment model..." << endl;
		alignmodel = new CoocAlignmentModel(sourcemodel,targetmodel, prunevalue);
		cerr << "   Found alignment targets for  " << alignmodel->alignprob.size() << " source constructions" << endl;
		cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;		
	}
	
	
	if ((!sourceclassfile.empty()) && (!targetclassfile.empty())) {
		cerr << "Loading source class decoder " << sourceclassfile << endl;
		ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
	
		cerr << "Loading target class decoder " << targetclassfile << endl;
		ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    	
	
	    cerr << "Decoding..." << endl;
	    alignmodel->decode(sourceclassdecoder, targetclassdecoder, &cout);    
	}	

	if (alignmodel != NULL) {
		delete alignmodel;
	}
	
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
