#include "alignmodel.h"

using namespace std;

void usage() {
    cerr << "Usage: aligner -J -s source-model -t target-model [-S source-class-file -T target-class-file]" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-s sourcemodelfile        Source graph model file (*.graphmodel.colibri)" << endl;    
    cerr << "\t-t targetmodelfile        Target model file (*.graphmodel.colibri)"  << endl;
    cerr << "\t-S sourceclassfile        Source class file (for decoding)" << endl;
    cerr << "\t-T targetclassfile        Target class file (for decoding)" << endl;
    //cerr << "\t-d model                 Load and decode an existing model" << endl; //TODO
    //cerr << "\t-B                       Do a bi-directional alignment and compute intersection of results" << endl; //TODO
	cerr << "\t-l n                      Minimum N length" << endl; //TODO
    cerr << "\t-L n                      Maximum N length" << endl; //TODO
    cerr << "\t-N                        No skip-grams" << endl; //TODO
    cerr << "\t-J                        Use Jaccard co-occurrence method (simplest)" << endl;
    cerr << "\t-D                        Use Dice co-occurrence method" << endl;
    //cerr << "\t-E                       Use EM alignment method" << endl; //TODO
    cerr << "\t-B                        Best alignment only" << endl;    
    cerr << "\t-P probability-threshold  Prune all alignments with an alignment probability lower than specified (0 <= x <= 1)" << endl;    
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;
    cerr << "\t-o occurence-threshold    Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold         Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold       Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio           Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
    cerr << "\t-Z				         No normalisation; return actual co-occurence scores instead of probabilities" << endl;
    cerr << "\t-V				         Verbose debugging output" << endl;
}

int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    double coocprunevalue = 0.0;
    double probprunevalue = 0.0;
    CoocMode COOCMODE = NOCOOC;
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    bool DOBIDIRECTIONAL = false;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
    bool DODEBUG = false;
    bool DONORM = true;
    bool BESTONLY = false;
    
    char c;    
    while ((c = getopt(argc, argv, "s:S:t:T:p:P:JDo:F:x:X:Bl:L:NV")) != -1)
        switch (c)
        {
        case 'B':
        	BESTONLY = true;
        	break;
        case 'p':
            coocprunevalue = atof(optarg);
            break;
        case 'P':
            probprunevalue = atof(optarg);
            break;            
        case 'J':
            COOCMODE = JACCARD;
            break;
        case 'D':
            COOCMODE = DICE;
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
		case 'l':
            MINLENGTH = atoi(optarg);
            break;            
		case 'L':
            MAXLENGTH = atoi(optarg);
            break;    
        case 'N':
            DOSKIPGRAMS = false;
            break;   
        case 'V':
        	DODEBUG = true;
        	break;
        case 'Z':
        	DONORM = false;
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
    
    if (!COOCMODE) {
    	cerr << "Error: No alignment method selected (select -J or -D)" << endl;
    	usage();
    	exit(3);
    }
    
    cerr << "Configuration: " << endl;
    if (COOCMODE == JACCARD) {
  		cerr << "\tCo-occcurrence metric : JACCARD (-J)" << endl;	
    } else if (COOCMODE == DICE) {
    	cerr << "\tCo-occcurrence metric : DICE (-D)" << endl;
    }
    cerr << "\tAlig. prob prune  (-P): " << probprunevalue << endl;
    cerr << "\tCo-oc prune value (-p): " << coocprunevalue << endl;    
	cerr << "\tCount threshold   (-o): " << COUNTTHRESHOLD << endl;
	cerr << "\tFreq threshold    (-F): " << FREQTHRESHOLD << endl;
	cerr << "\tXcount threshold  (-x): " << XCOUNTTHRESHOLD << endl;
	cerr << "\tXcount ratio      (-X): " << XCOUNTRATIOTHRESHOLD << endl;
	cerr << "\tMinimum N length  (-l): " << MINLENGTH << endl;
	cerr << "\tMaximum N length  (-L): " << MAXLENGTH << endl;
	if (!DOSKIPGRAMS) {
		cerr << "\tSKIPGRAMS DISABLED! (-N)";
	}
	if (!DONORM) {
		cerr << "\tSKIPGRAMS DISABLED! (-N)";
	}
	cerr << endl;
	
    
    cerr << "Loading source model " << sourcemodelfile << endl;
    SelectivePatternModel sourcemodel = SelectivePatternModel(sourcemodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS, MINLENGTH, MAXLENGTH);
    cerr << "  Loaded " << sourcemodel.types() << " types, " << sourcemodel.tokens() << " tokens" << endl;
    cerr << "  Ignored " << sourcemodel.ignoredtypes << " types, " << sourcemodel.ignoredtokens << " tokens due to set thresholds" << endl;
    if (sourcemodel.has_xcount()) {
    	cerr << "  Exclusive count available? YES" << endl;
    } else {
    	cerr << "  Exclusive count available? NO" << endl;
    }
    if (sourcemodel.has_index()) {
    	cerr << "  Reverse index has " << sourcemodel.reverseindex.size() << " sentences" << endl;
    } else {
    	cerr << "ERROR: Model " + sourcemodelfile + " contains no indexing information! Unable to align without!" << endl;
    	exit(3);
    }    
    
    cerr << "Loading target model " << targetmodelfile << endl;
    SelectivePatternModel targetmodel = SelectivePatternModel(targetmodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS, MINLENGTH, MAXLENGTH);
    cerr << "  Loaded " << targetmodel.types() << " types, " << targetmodel.tokens() << " tokens" << endl;
    cerr << "  Ignored " << targetmodel.ignoredtypes << " types, " << targetmodel.ignoredtokens << " tokens due to set thresholds" << endl;
    if (targetmodel.has_xcount()) {
    	cerr << "  Exclusive count available? YES" << endl;
    } else {
    	cerr << "  Exclusive count available? NO" << endl;
    }
    if (targetmodel.has_index()) {
    	cerr << "  Reverse index has " << targetmodel.reverseindex.size() << " sentences" << endl;
    } else {
    	cerr << "ERROR: Model " + targetmodelfile + " contains no indexing information! Unable to align without!" << endl;
    	exit(3);
    }       
    
    AlignmentModel * alignmodel = NULL;
    
    if (COOCMODE) {
		cerr << "Computing alignment model..." << endl;
		alignmodel = new CoocAlignmentModel(COOCMODE, sourcemodel,targetmodel, coocprunevalue, probprunevalue, BESTONLY, DONORM, DODEBUG);
		cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
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
