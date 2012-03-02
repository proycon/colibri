#include <getopt.h>
#include "alignmodel.h"


using namespace std;

void usage() {
    cerr << "Usage: aligner [-J|-D|-E] -s source-model -t target-model [-S source-class-file -T target-class-file]" << endl;
    cerr << " Input:" << endl;
    cerr << "\t-s sourcemodelfile        Source graph model file (*.graphmodel.colibri)" << endl;    
    cerr << "\t-t targetmodelfile        Target model file (*.graphmodel.colibri)"  << endl;
    cerr << "\t-S sourceclassfile        Source class file (for decoding)" << endl;
    cerr << "\t-T targetclassfile        Target class file (for decoding)" << endl;
    cerr << "\t-d alignmodelfile         Decode an existing alignment model (*.alignmodel.colibri), specify -S and -T as well" << endl;
    cerr << " Alignment method (choose one):" << endl;
    cerr << "\t-J                        Use Jaccard co-occurrence method (simplest)" << endl;
    cerr << "\t-D                        Use Dice co-occurrence method" << endl;
    cerr << "\t-E                        Use EM alignment method" << endl;       
    cerr << " Generic alignment options:" << endl;    
    cerr << "\t-V				         Verbose debugging output" << endl;
    cerr << "\t-b                        Best alignment only" << endl;
    cerr << "\t-G 			             Weigh alignment results based on graph information (subsumption relations)" << endl;
    cerr << "\t-B probability-threshold  Compute bidirectional alignment (intersection), using given probability threshold" << endl;
    cerr << " Co-occurrence alignment options:" << endl;       
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;
    cerr << "\t-Z				         Do normalisation; return probabilities instead of co-occurrence scores" << endl;   
    cerr << " EM Alignment Options:" << endl;
    cerr << "\t-P probability-threshold  Prune all alignments with an alignment probability lower than specified (0 <= x <= 1)" << endl;
    cerr << "\t-I n				         Maximum number of iterations (for EM method, default: 10000)" << endl;
    cerr << "\t-v n				         Convergence delta value (for EM method, default: 0.001)" << endl;
    cerr << " Input filtering:" << endl;
    cerr << "\t-O occurence-threshold    Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold         Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold       Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio           Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
    cerr << "\t-l n                      Minimum N length" << endl; 
    cerr << "\t-L n                      Maximum N length" << endl; 
    cerr << "\t-N                        No skip-grams" << endl;
    cerr << " Output options:" << endl;
    cerr << "\t--simplelex               Output simple word-based lexicon" << endl;
    cerr << "\t--simpletable             Output simple phrase-based translation table" << endl;
    cerr << "\t--targetfirst             Output target before source in simple lexicon and simple translation table output" << endl;
}

int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    string modelfile="";
    double coocprunevalue = 0.0;
    double probprunevalue = 0.0;
    double graphweightfactor = 1.2;
    CoocMode COOCMODE = NOCOOC;
    bool DO_EM = false;
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    bool DOBIDIRECTIONAL = false;
    bool DOGRAPHALIGN = false;
    double bidirprobthreshold = 0.0;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
    bool DODEBUG = false;
    bool DONORM = false;
    bool BESTONLY = false;
    int MAXROUNDS = 10000;
    double CONVERGENCE = 0.001;
    int DOSIMPLELEX = 0;
    int DOSIMPLETABLE = 0;
    int TARGETFIRST = 0;
    string outputprefix = "";
    
    static struct option long_options[] = {      
       {"simplelex", no_argument,       &DOSIMPLELEX, 1},
       {"simpletable", no_argument,       &DOSIMPLETABLE, 1},
       {"targetfirst", no_argument,       &TARGETFIRST, 1},            
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    char c;    
    while ((c = getopt_long(argc, argv, "hd:s:S:t:T:p:P:JDo:O:F:x:X:B:bl:L:NVZEI:v:G",long_options,&option_index)) != -1)
        switch (c)
        {
        case 0:
            if (long_options[option_index].flag != 0)
               break;
        case 'd':
        	modelfile = optarg;
        	break;
        case 'h':
        	usage();
        	exit(0);
        case 'B':
        	bidirprobthreshold = atof(optarg);
        	DOBIDIRECTIONAL = true;
        	break;
        case 'b':
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
        case 'E':
        	DO_EM = true;
        	break;
		case 'G':
			DOGRAPHALIGN = true;
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
        case 'O':
            COUNTTHRESHOLD = atoi(optarg);            
            break;
        case 'o':
        	outputprefix = optarg;
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
        case 'v':       
            CONVERGENCE = atof(optarg);
            break;
		case 'L':
            MAXLENGTH = atoi(optarg);
            break;
        case 'I':
            MAXROUNDS = atoi(optarg);
            break;        
        case 'N':
            DOSKIPGRAMS = false;
            break;   
        case 'V':
        	DODEBUG = true;
        	break;
        case 'Z':
        	DONORM = true;
        	break;        	
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
        
	
    if (modelfile.empty() && (sourcemodelfile.empty()  || targetmodelfile.empty())) {
  	    cerr << "Error: Specify at least a source model, target model, and alignment method to build an alignment model! Or load a pre-existing model" << endl;
        usage();
        exit(2);
    }
    
    if (modelfile.empty() && ((!DO_EM) && (!COOCMODE))) {
    	cerr << "Error: No alignment method selected (select -J or -D)" << endl;
    	usage();
    	exit(3);
    }

	AlignmentModel * alignmodel = NULL;

	if (modelfile.empty()) {
		cerr << "Configuration: " << endl;
		if (DO_EM) {
			cerr << "\tEM-alignment" << endl;
		} else if (COOCMODE == JACCARD) {
			cerr << "\tCo-occurrence metric : JACCARD (-J)" << endl;	
		} else if (COOCMODE == DICE) {
				cerr << "\tCo-occurrence metric : DICE (-D)" << endl;
		}
		cerr << "\tAlig. prob prune  (-P): " << probprunevalue << endl;
		if (DO_EM) cerr << "\tCo-oc prune value (-p): " << coocprunevalue << endl;    
		cerr << "\tCount threshold   (-o): " << COUNTTHRESHOLD << endl;
		cerr << "\tFreq threshold    (-F): " << FREQTHRESHOLD << endl;
		cerr << "\tXcount threshold  (-x): " << XCOUNTTHRESHOLD << endl;
		cerr << "\tXcount ratio      (-X): " << XCOUNTRATIOTHRESHOLD << endl;
		cerr << "\tMinimum N length  (-l): " << MINLENGTH << endl;
		cerr << "\tMaximum N length  (-L): " << MAXLENGTH << endl;
		if (!DOSKIPGRAMS) {
			cerr << "\tSKIPGRAMS DISABLED! (-N)"  << endl;
		}
		if (DONORM) {
			cerr << "\tNormalisation enabled (-Z)"  << endl;
		}
		if (DOBIDIRECTIONAL) {
			cerr << "\tBidirectional alignment enabled (-B)" << endl;
		}
		if (DOGRAPHALIGN) {
			cerr << "\tGraph weighting enabled (-G), weight factor: " << graphweightfactor << endl;
		}
		cerr << endl;
	
		
		cerr << "Loading source model " << sourcemodelfile << endl;
		SelectivePatternModel sourcemodel = SelectivePatternModel(sourcemodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS, MINLENGTH, MAXLENGTH, DOGRAPHALIGN);
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
		if (sourcemodel.has_parents()) {
			cerr << "  Parent relations available for  " << sourcemodel.rel_subsumption_parents.size() << " patterns" << endl;
		}
		
		cerr << "Loading target model " << targetmodelfile << endl;
		SelectivePatternModel targetmodel = SelectivePatternModel(targetmodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS, MINLENGTH, MAXLENGTH, DOGRAPHALIGN);
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
		if (targetmodel.has_parents()) {
			cerr << "  Parent relations available for  " << targetmodel.rel_subsumption_parents.size() << " patterns" << endl;
		}		
		
		if (DO_EM) {
			cerr << "Computing alignment model..." << endl;
			alignmodel = new EMAlignmentModel(&sourcemodel,&targetmodel, MAXROUNDS,  CONVERGENCE, probprunevalue, BESTONLY, DODEBUG);
			cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
			cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;
		
			if (DOBIDIRECTIONAL) {
				cerr << "Computing reverse alignment model (for bidirectional alignment)..." << endl;
				AlignmentModel reversealignmodel = EMAlignmentModel(&targetmodel,&sourcemodel, MAXROUNDS,  CONVERGENCE, probprunevalue, BESTONLY, DODEBUG);
				cerr << "   Found alignment targets for  " << reversealignmodel.alignmatrix.size() << " source constructions" << endl;
				cerr << "   Total of alignment possibilies in matrix: " << reversealignmodel.totalsize() << endl;
				cerr << "Computing intersection of both alignment models..." << endl;
				alignmodel->intersect(&reversealignmodel, bidirprobthreshold);
			}	    
		} else if (COOCMODE) {
			cerr << "Computing alignment model..." << endl;
			alignmodel = new CoocAlignmentModel(COOCMODE, &sourcemodel,&targetmodel, coocprunevalue, probprunevalue, BESTONLY, DONORM, DODEBUG);
			cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
			cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;
		
			if (DOBIDIRECTIONAL) {
				cerr << "Computing reverse alignment model (for bidirectional alignment)..." << endl;
				AlignmentModel reversealignmodel = CoocAlignmentModel(COOCMODE, &targetmodel,&sourcemodel, coocprunevalue, probprunevalue, BESTONLY, DONORM, DODEBUG);
				cerr << "   Found alignment targets for  " << reversealignmodel.alignmatrix.size() << " source constructions" << endl;
				cerr << "   Total of alignment possibilies in matrix: " << reversealignmodel.totalsize() << endl;
				cerr << "Computing intersection of both alignment models..." << endl;
				alignmodel->intersect(&reversealignmodel, bidirprobthreshold);	
			}	    				
		}

		if (DOGRAPHALIGN) {
			cerr << "Weighting based on graph subsumption relations..." << endl;
			const int adjustments = alignmodel->graphalign(sourcemodel, targetmodel, graphweightfactor);
			cerr << "   Made " << adjustments << " adjustments" << endl;			
		}

		if (!outputprefix.empty()) {
		    cerr << "Saving alignment model..." << endl;
			alignmodel->save(outputprefix + ".alignmodel.colibri");
		}


		if ((!sourceclassfile.empty()) && (!targetclassfile.empty())) {
			cerr << "Loading source class decoder " << sourceclassfile << endl;
			ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
	
			cerr << "Loading target class decoder " << targetclassfile << endl;
			ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    	
	
			cerr << "Decoding..." << endl;
			if (DOSIMPLETABLE) {
				alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout);
			} else if (DOSIMPLELEX) {
				alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, true);
			} else { 
				alignmodel->decode(sourceclassdecoder, targetclassdecoder, &cout);
			}       
		}	
		
		
		       
    } else {
   		cerr << "Loading alignment model..." << endl;
    	alignmodel = new AlignmentModel(modelfile);
    	
    	if ((!sourceclassfile.empty()) && (!targetclassfile.empty())) {
			cerr << "Loading source class decoder " << sourceclassfile << endl;
			ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
	
			cerr << "Loading target class decoder " << targetclassfile << endl;
			ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    	
	
			cerr << "Decoding..." << endl;
			if (DOSIMPLETABLE) {
				alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST);
			} else if (DOSIMPLELEX) {
				alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST, true);
			} else { 
				alignmodel->decode(sourceclassdecoder, targetclassdecoder, &cout);
			}        
		}	    	
    } 
    
	


	if (alignmodel != NULL) {
		delete alignmodel;
	}
	
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
