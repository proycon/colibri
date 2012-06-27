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
    cerr << "\t-i inv-alignmodelfile     Use inverse alignment model as well (*.alignmodel.colibri), specify -S and -T as well" << endl;
    cerr << " Alignment method (choose one, though some may be combined):" << endl;
    cerr << "\t-J                        Use Jaccard co-occurrence method (simplest)" << endl;
    //cerr << "\t-D                        Use Dice co-occurrence method" << endl;
    cerr << "\t-E                        Use EM alignment method (sentence-based)" << endl;
    cerr << "\t-2                        Use Alternative EM alignment method (type-based)" << endl;
    cerr << "\t-3                        Use Iterative EM alignment method" << endl;
    cerr << "\t-W giza-s-t.A3:giza-t-s.A3   Extract phrases by matching giza word-alignments with pattern models" << endl;       
    cerr << " Generic alignment options:" << endl;    
    cerr << "\t-V				         Verbose debugging output" << endl;
    cerr << "\t-b n                      Best n alignments only" << endl;
    cerr << "\t-G weight-factor          Weigh alignment results based on graph information (subsumption relations)" << endl;
    cerr << "\t-B probability-threshold  Compute bidirectional alignment (intersection), using given probability threshold (0 <= x < 1). Will automatically enable normalisation (-Z)" << endl;
    cerr << " Co-occurrence alignment options:" << endl;       
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;
    cerr << "\t-Z				         Do normalisation" << endl;
    cerr << "\t-U                        Extract skip-grams from n-grams (requires source and target models to be graph models with template and instance relations)" << endl;   
    cerr << " EM Alignment Options:" << endl;
    cerr << "\t-P probability-threshold  Prune all alignments with an alignment probability lower than specified (0 <= x <= 1)" << endl;
    cerr << "\t-I n				         Maximum number of iterations (for EM method, default: 10000)" << endl;
    cerr << "\t-v n				         Convergence delta value (for EM method, default: 0.001)" << endl;
    cerr << "\t-N                        Do not extract skip-grams in EM-process" << endl;
    cerr << " GIZA Alignment Options:" << endl;
    cerr << "\t-a                        Alignment threshold (0 <= x <= 1). Specifies how strong word alignments have to be if phrases are to be extracted from them (default 0.5)" << endl;
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a jaccard co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;
    cerr << "\t-c pair-count-threshold   Prune phrase pairs that occur less than specified" << endl;
    cerr << " Input filtering:" << endl;
    cerr << "\t-O occurence-threshold    Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold         Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold       Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio           Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
    cerr << "\t-l n                      Minimum N length" << endl; 
    cerr << "\t-L n                      Maximum N length" << endl;     
    cerr << "\t--null                    Take into account zero-fertility words (null alignments) in EM" << endl;
    cerr << " Output options:" << endl;
    cerr << "\t--simplelex               Output simple word-based lexicon" << endl;
    cerr << "\t--simpletable             Output simple phrase-based translation table" << endl;
    cerr << "\t--targetfirst             Output target before source in simple lexicon and simple translation table output (use with --simplelex, --simpletable)" << endl;
    cerr << "\t--moses                   Output phrase-translation table in Moses format (use with --simpletable, --simplelex)" << endl;
}




int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    string modelfile="";
    string invmodelfile="";
    double coocprunevalue = 0.0;
    double probprunevalue = 0.0;
    double graphweightfactor = 0.0; 
    CoocMode COOCMODE = NOCOOC;
    bool DO_EM = false;
    bool DO_EM2 = false;
    bool DO_ITEREM = false;
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    bool DOBIDIRECTIONAL = false;
    double bidirprobthreshold = 0.0;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
    bool EXTRACTSKIPGRAMS = false;
    bool DODEBUG = false;
    bool DONORM = false;
    bool DOGIZA = false;
    string gizast = "";
    string gizats = "";
    int EM_NULL = 0;
    int MAXROUNDS = 10000;
    double CONVERGENCE = 0.001;
    int DOSIMPLELEX = 0;
    int DOSIMPLETABLE = 0;
    int TARGETFIRST = 0;
    int MOSESFORMAT = 0;
    int bestn = 0;
    bool DEBUG = false;
    
    double alignthreshold = 0.5;
    int pairthreshold = 1;
    
    string outputprefix = "";
    
    static struct option long_options[] = {      
       {"simplelex", no_argument,       &DOSIMPLELEX, 1},
       {"simpletable", no_argument,       &DOSIMPLETABLE, 1},
       {"targetfirst", no_argument,       &TARGETFIRST, 1},
       {"moses", no_argument,             &MOSESFORMAT, 1},
       {"null", no_argument,             &EM_NULL, 1}, 
                      
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    
    string raw;
    string::size_type pos;
    
    
    
    char c;    
    while ((c = getopt_long(argc, argv, "hd:s:S:t:T:p:P:JDo:O:F:x:X:B:b:l:L:NVZEI:v:G:i:23W:a:c:U",long_options,&option_index)) != -1)
        switch (c)
        {
        case 0:
            if (long_options[option_index].flag != 0)
               break;
        case 'd':
        	modelfile = optarg;
        	break;
        case 'i':
        	invmodelfile = optarg;
        	break;
        case 'h':
        	usage();
        	exit(0);
        case 'B':
        	bidirprobthreshold = atof(optarg);        	
        	DOBIDIRECTIONAL = true;
        	DONORM = true;
        	break;
        case 'b':
        	bestn = atoi(optarg);
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
        	DEBUG = true;
        	break;
        /*case 'D':
            COOCMODE = DICE;
            break;*/            
        case 'E':
        	DO_EM = true;
        	break;
        case '2':
        	DO_EM2 = true;
        	break;               	
        case '3':
        	DO_ITEREM = true;
        	break;        	
		case 'G':
			graphweightfactor = atof(optarg);
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
        case 'U':
            EXTRACTSKIPGRAMS = true;
            break;   
        case 'V':
        	DODEBUG = true;
        	break;        	
        case 'W':
            DOGIZA = true;
            raw = optarg;
            pos = raw.find(':');
            if (pos == string::npos) {
                cerr << "ERROR: -W expects two giza filenames separated by a colon (:)" << endl;
                usage();
                exit(2);
            }
            gizast = raw.substr(0, pos);
            gizats = raw.substr(pos+1); 
            break;
        case 'a':
            alignthreshold = atoi(optarg);
            break;
        case 'c':
            pairthreshold = atoi(optarg);
            break;
        case 'Z':
        	DONORM = true;
        	break;        	
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
    //if (MOSESFORMAT) TARGETFIRST = 1;    
	
    if (modelfile.empty() && (sourcemodelfile.empty()  || targetmodelfile.empty())) {
  	    cerr << "Error: Specify at least a source model, target model, and alignment method to build an alignment model! Or load a pre-existing model" << endl;
        usage();
        exit(2);
    }
    
    if (modelfile.empty() && ((!DO_EM) && (!DO_EM2) && (!DO_ITEREM) && (!COOCMODE) && (!DOGIZA) )) {
    	cerr << "Error: No alignment method selected (select -J or -D)" << endl;
    	usage();
    	exit(3);
    }

	AlignmentModel * alignmodel = NULL;

	if (modelfile.empty()) {
		cerr << "Configuration: " << endl;
		if (DO_EM) {
			cerr << "\tEM-alignment (-E)" << endl;
		} else if (DO_EM2) {
			cerr << "\tAlternative EM-alignment (-2)" << endl;
		} else if (DO_ITEREM) {
			cerr << "\tIterative EM-alignment (-3)" << endl;
		} else if (COOCMODE == JACCARD) {
			cerr << "\tCo-occurrence metric : JACCARD (-J)" << endl;	
		} else if (COOCMODE == DICE) {
				cerr << "\tCo-occurrence metric : DICE (-D)" << endl;
        } else if (DOGIZA) {
            cerr << "\tGIZA++ source-to-target word-alignment: " << gizast << endl;
            cerr << "\tGIZA++ target-to-source word-alignment: " << gizats << endl;
		}
		
		if (bestn) {
	    	cerr << "\tBest-n            (-b): " << bestn << endl;
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
			cerr << "\tBidirectional alignment enabled (-B), threshold: " <<  bidirprobthreshold << endl;
		}
		if (graphweightfactor > 0) {
			cerr << "\tGraph weighting enabled (-G), weight factor: " << graphweightfactor << endl;
		}
		if ((DO_EM) || (DO_ITEREM) || (DO_EM2)) {
			if (EM_NULL) {
				cerr << "\tNull alignments in EM?  yes" << endl;
			} else {
				cerr << "\tNull alignments in EM?  no" << endl;
			}
		}
		cerr << endl;
	
		
		cerr << "Loading source model " << sourcemodelfile << endl;
		SelectivePatternModel sourcemodel = SelectivePatternModel(sourcemodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS || EXTRACTSKIPGRAMS, MINLENGTH, MAXLENGTH, (graphweightfactor > 0),false, NULL,false, DEBUG);
		cerr << "  Loaded " << sourcemodel.types() << " types, " << sourcemodel.tokens() << " tokens" << endl;
	 	cerr << "  Ignored " << sourcemodel.ignoredtypes << " types, " << sourcemodel.ignoredoccurrences << " occurrences due to set thresholds" << endl;
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
		SelectivePatternModel targetmodel = SelectivePatternModel(targetmodelfile, true, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS || EXTRACTSKIPGRAMS, MINLENGTH, MAXLENGTH, (graphweightfactor > 0),false, NULL,false, DEBUG);
		cerr << "  Loaded " << targetmodel.types() << " types, " << targetmodel.tokens() << " tokens" << endl;
		cerr << "  Ignored " << targetmodel.ignoredtypes << " types, " << targetmodel.ignoredoccurrences << " occurrences due to set thresholds" << endl;
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
		
		
		
		alignmodel = new AlignmentModel(&sourcemodel,&targetmodel, DODEBUG);
		AlignmentModel * reversealignmodel = new AlignmentModel(&targetmodel,&sourcemodel, DODEBUG); 				
		bool EM_INIT = true;
		
		if (COOCMODE) {
			int tmpbestn;
			if ((DO_EM) || (DO_EM2)) {
			 	tmpbestn = 0;
			} else {
				tmpbestn = bestn;
			}
		
			cerr << "Computing Cooc alignment model..." << endl;
			alignmodel->trainCooc(COOCMODE, tmpbestn, coocprunevalue, 0);			
			cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
			cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;			
			if ((DONORM) || (DO_EM) || (DO_EM2)) {
				cerr << "   Normalizing... " << endl;
				alignmodel->normalize();
			}

			if (DOBIDIRECTIONAL) {
				cerr << "Computing reverse alignment model (for bidirectional alignment)..." << endl;
				reversealignmodel->trainCooc(COOCMODE, tmpbestn, coocprunevalue, 0);
				cerr << "   Found alignment targets for  " << reversealignmodel->alignmatrix.size() << " source constructions" << endl;
				cerr << "   Total of alignment possibilies in matrix: " << reversealignmodel->totalsize() << endl;
				if ((DONORM) || (DO_EM) || (DO_EM2)) {
					cerr << "   Normalizing... " << endl;
					reversealignmodel->normalize();
				}			
		
			}
				    					
			EM_INIT = false;
			
		}		
		if ((DO_EM) || (DO_EM2)) {
			cerr << "Computing EM alignment model..." << endl;
			if (DO_EM2) {
			    alignmodel->trainEM2(MAXROUNDS,  CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT); //EM2 (experimental)
			} else {		
			    alignmodel->trainEM(MAXROUNDS,  CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT);
			}
			if (DONORM) alignmodel->normalize();	
			cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
			cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;
						
			if (DOBIDIRECTIONAL) {
				cerr << "Computing reverse alignment model (for bidirectional alignment)..." << endl;
				if (DO_EM2) {				
				    reversealignmodel->trainEM2(MAXROUNDS, CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT);
				} else {
				    reversealignmodel->trainEM(MAXROUNDS, CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT);
				}
				if (DONORM) reversealignmodel->normalize();
				cerr << "   Found alignment targets for  " << reversealignmodel->alignmatrix.size() << " source constructions" << endl;
				cerr << "   Total of alignment possibilies in matrix: " << reversealignmodel->totalsize() << endl;						
			}	    			    		
		}		
		if (DOGIZA) {
			cerr << "Loading source class encoder " << sourceclassfile << endl;
		    ClassEncoder sourceclassencoder = ClassEncoder(sourceclassfile);
    
		    cerr << "Loading target class encoder " << targetclassfile << endl;
    		ClassEncoder targetclassencoder = ClassEncoder(targetclassfile);    
				
		    cerr << "Initialising GIZA++ Word Alignments" << endl;
		    GizaModel gizamodels2t = GizaModel(gizast, &sourceclassencoder, &targetclassencoder);
		    GizaModel gizamodelt2s = GizaModel(gizats, &targetclassencoder, &sourceclassencoder);
		    
		    alignmodel->extractgizapatterns(gizamodels2t, gizamodelt2s, pairthreshold, coocprunevalue, alignthreshold);
		}

	    if (EXTRACTSKIPGRAMS) {
	        alignmodel->extractskipgrams();
	    }
				
		if (DOBIDIRECTIONAL) {
			cerr << "Computing intersection of both alignment models..." << endl;
			alignmodel->intersect(reversealignmodel, bidirprobthreshold, bestn);
			if (DONORM) alignmodel->normalize();
		}
			
			
		if (graphweightfactor > 0) {
			cerr << "Weighting based on graph subsumption relations..." << endl;
			const int adjustments = alignmodel->graphalign(sourcemodel, targetmodel, graphweightfactor);
			cerr << "   Made " << adjustments << " adjustments" << endl;			
		}




		if (!outputprefix.empty()) {
		    cerr << "Saving alignment model..." << endl;
			alignmodel->save(outputprefix);
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
    	if ((sourceclassfile.empty()) || (targetclassfile.empty())) {
    	    if ((EXTRACTSKIPGRAMS) && (!outputprefix.empty())) {
    	    	cerr << "Loading alignment model..." << endl;
        		alignmodel = new AlignmentModel(modelfile, bestn);
	            alignmodel->extractskipgrams();
		        cerr << "Saving alignment model..." << endl;
			    alignmodel->save(outputprefix);	    
    	    } else {
    		    cerr << "Error: Specify -S and -T to decode, or -U and -o to extract skipgrams from an existing model" << endl; 
    		    usage();
    		    exit(2);
    		  }
    	} else {
    	
		    cerr << "Loading source class decoder " << sourceclassfile << endl;
		    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);

		    cerr << "Loading target class decoder " << targetclassfile << endl;
		    ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    	

        	
        	if (invmodelfile.empty()) {
        	   	cerr << "Loading alignment model..." << endl;
        		alignmodel = new AlignmentModel(modelfile, bestn);
        		
        		if  (EXTRACTSKIPGRAMS) alignmodel->extractskipgrams();
        		    	    	
			    cerr << "Decoding..." << endl;
			    if (DOSIMPLETABLE) {
				    alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST, false, MOSESFORMAT);
			    } else if (DOSIMPLELEX) {
				    alignmodel->simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST, true);
			    } else { 
				    alignmodel->decode(sourceclassdecoder, targetclassdecoder, &cout);
			    }            		
        	} else {
        		cerr << "Loading alignment models..." << endl;
        		BiAlignmentModel bialignmodel = BiAlignmentModel(modelfile, invmodelfile);
			    cerr << "Decoding..." << endl;
			    if (DOSIMPLETABLE) {
				    bialignmodel.simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST, false, MOSESFORMAT);
			    } else if (DOSIMPLELEX) {
				    bialignmodel.simpletableoutput(sourceclassdecoder, targetclassdecoder, &cout, TARGETFIRST, true);
			    } else { 
				    bialignmodel.decode(sourceclassdecoder, targetclassdecoder, &cout);
			    }            		    		
		    }
			
		}    	
    } 
    

	if (alignmodel != NULL) {
		delete alignmodel;
	}
	
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
