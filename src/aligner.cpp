#include "alignmodel.h"

using namespace std;

void usage() {
    cerr << "Aligner: aligner -J -s source-model -t target-model [-S source-class-file -T target-class-file]" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-s sourcemodelfile       Source graph model file (*.graphmodel.colibri)" << endl;    
    cerr << "\t-t targetmodelfile       Target model file (*.graphmodel.colibri)"  << endl;
    cerr << "\t-S sourceclassfile       Source class file (for decoding)" << endl;
    cerr << "\t-T targetclassfile       Target class file (for decoding)" << endl;
    //cerr << "\t-d model                 Load and decode an existing model" << endl; //TODO
    //cerr << "\t-B                       Do a bi-directional alignment and compute intersection of results" << endl; //TODO
    	cerr << "\t-l n                     Minimum N length" << endl; //TODO
    cerr << "\t-L n                     Maximum N length" << endl; //TODO
    cerr << "\t-N                       No skip-grams" << endl; //TODO
    cerr << "\t-J                       Use Jaccard co-occurrence method (simplest)" << endl;
    cerr << "\t-D                       Use Dice co-occurrence method" << endl;
    //cerr << "\t-E                       Use EM alignment method" << endl; //TODO
    cerr << "\t-p pruning-threshold     Prune all alignments with a lower score than specified (0 <= x <= 1)" << endl;
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
    CoocMode COOCMODE = NOCOOC;
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    bool DOBIDIRECTIONAL = false;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
    
    char c;    
    while ((c = getopt(argc, argv, "s:S:t:T:p:JDo:F:x:X:B:l:L:N")) != -1)
        switch (c)
        {
        case 'B':
        	DOBIDIRECTIONAL = true;
        	break;
        case 'p':
            prunevalue = atof(optarg);
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
    cerr << "\tPrune value       (-P): " << prunevalue << endl;
	cerr << "\tCount threshold   (-o): " << COUNTTHRESHOLD << endl;
	cerr << "\tFreq threshold    (-F): " << FREQTHRESHOLD << endl;
	cerr << "\tXcount threshold  (-x): " << XCOUNTTHRESHOLD << endl;
	cerr << "\tXcount ratio      (-X): " << XCOUNTRATIOTHRESHOLD << endl;
	cerr << "\tMinimum N length  (-l): " << MINLENGTH << endl;
	cerr << "\tMaximum N length  (-L): " << MAXLENGTH << endl;
	if (!DOSKIPGRAMS) {
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
    cerr << "  Ignored " << sourcemodel.ignoredtypes << " types, " << sourcemodel.ignoredtokens << " tokens due to set thresholds" << endl;
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
		alignmodel = new CoocAlignmentModel(COOCMODE, sourcemodel,targetmodel, prunevalue);
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
