#include <getopt.h>
#include "alignmodel.h"


using namespace std;

void usage() {
    cerr << "Usage: aligner [-J|-E] -s source-model -t target-model [-S source-class-file -T target-class-file]" << endl;
    cerr << " Input:" << endl;
    cerr << "\t-s sourcemodelfile        Source model file (patternmodel/graphmodel)" << endl;    
    cerr << "\t-t targetmodelfile        Target model file (patternmodel/graphmodel)"  << endl;
    cerr << "\t-S sourceclassfile        Source class file (for model decoding)" << endl;
    cerr << "\t-T targetclassfile        Target class file (for model decoding)" << endl;
    cerr << "\t-d alignmodelfile         Load an existing alignment model (*.alignmodel.colibri), for model decoding specify with -S and -T" << endl;
    cerr << "\t-i alignmodelfile         Load inverse alignment model as well (*.alignmodel.colibri) and compute intersection, for model decoding specify with -S and -T" << endl;
    //cerr << "\t-H translationtable       Load a translation table model (*.transtable.colibri) for decoding, specify with -S and -T" << endl;    
    cerr << " Alignment method (choose one, though some may be combined):" << endl;
    cerr << "\t-J                        Use Jaccard co-occurrence method" << endl;
    //cerr << "\t-D                        Use Dice co-occurrence method" << endl;
    cerr << "\t-E                        Use EM alignment method" << endl;
    //cerr << "\t-2                        Use Alternative EM alignment method (type-based)" << endl;
    //cerr << "\t-3                        Use Iterative EM alignment method" << endl;
    cerr << "\t-W giza-s-t.A3:giza-t-s.A3   Extract phrases by matching giza word-alignments with pattern models (supervised, if -s and -t are specified), or using heuristic methods (unsupervised, if -H is specified). Specify two GIZA alignment models (one for each direction), separated by a colon" << endl;
    cerr << "\t-H method                 Heuristic method for unsupervised phrase alignment. Choose from: growdiag, growdiagfinal (default), s2t, intersection, union" << endl;    
    cerr << " Generic alignment options:" << endl;    
    cerr << "\t-I 1                      Compute intersection (one single joint score)" << endl;
    cerr << "\t-I 2                      Compute intersection (two scores p(s|t),p(t|s) )" << endl;
    cerr << "\t-V				         Verbose debugging output" << endl;
    cerr << "\t-b n                      Best n alignments only" << endl;
    cerr << "\t-G weight-factor          Weigh alignment results based on graph information (subsumption relations)" << endl;
    cerr << "\t-B probability-threshold  Probability threshold used when computing bidirectional alignment (intersection) of alignment model and reverse alignment model (-I, -i)" << endl;
    cerr << "\t-z				         No normalisation" << endl;
    cerr << "\t-U                        Extract skip-grams from n-grams (requires source and target models to be graph models with template and instance relations)" << endl;    
    cerr << " Co-occurrence alignment options:" << endl;       
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;   
    cerr << " EM Alignment Options:" << endl;
    cerr << "\t-P probability-threshold  Prune all alignments with an alignment probability lower than specified (0 <= x <= 1)" << endl;
    cerr << "\t-M n				         Maximum number of iterations (for EM method, default: 10000)" << endl;
    cerr << "\t-v n				         Convergence delta value (for EM method, default: 0.001)" << endl;
    cerr << "\t-N                        Do not extract skip-grams in EM-process" << endl;
    cerr << "\t--null                    Take into account zero-fertility words (null alignments) in EM" << endl;    
    cerr << " GIZA Alignment Options:" << endl;
    cerr << "\t-a                        Alignment threshold (0 <= x <= 1). Specifies how strong word alignments have to be if phrases are to be extracted from them (default 0.5)" << endl;
    cerr << "\t-p cooc-pruning-threshold Prune all alignments with a jaccard co-occurence score lower than specified (0 <= x <= 1). Uses heuristics to prune, final probabilities may turn out lower than they would otherwise be" << endl;
    cerr << "\t-c pair-count-threshold   Prune phrase pairs that occur less than specified" << endl;
    cerr << " Input filtering:" << endl;
    cerr << "\t-O occurence-threshold    Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold         Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold       Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio           Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
    cerr << "\t-l n                      Left context size (in words, default 0)" << endl;
    cerr << "\t-r n                      Right context size (in words, default 0)" << endl; 
    cerr << "\t-L n                      Maximum N length" << endl;         
    cerr << " Output options:" << endl;    
    cerr << "\t-o filename               Write an alignment model to file using this filename (extension *.alignmodel.colibri will be automatically added)",
    cerr << "\t--moses                   Output phrase-translation table in Moses format" << endl;
    cerr << "\t--stats                   Output statistics only (use with -d)" << endl;
    cerr << "\t--removecontext           Generate a model without context from a model with context" << endl;    
}




int main( int argc, char *argv[] ) {
    string sourcemodelfile = "";
    string targetmodelfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    string modelfile="";
    string ttablefile="";
    string ttableoutfile="";
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
    int DOBIDIRECTIONAL = 0;
    double bidirprobthreshold = 0.0;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
    bool EXTRACTSKIPGRAMS = false;
    bool DODEBUG = false;
    bool DONORM = true;
    bool DOGIZA = false;
    string gizast = "";
    string gizats = "";
    int EM_NULL = 0;
    int MAXROUNDS = 10000;
    double CONVERGENCE = 0.001;
    //int DOSIMPLELEX = 0;
    //int DOSIMPLETABLE = 0;
    //int TARGETFIRST = 0;
    int MOSESFORMAT = 0;
    int DOSTATS = 0;
    int REMOVECONTEXT = 0;
    
    int REMOVECONTEXT = 0;
    int bestn = 0;
    
    bool DOPARENTS = false;
    bool DOCHILDREN = false;
    bool DOXCOUNT = false;
    bool DOSUCCESSORS = false;
    bool DOPREDECESSORS = false;
    bool DOSKIPCONTENT = false;
    bool DOSKIPUSAGE = false;    
    bool DOTEMPLATES = false; 
    bool DOINSTANCES = false;    
    PhraseAlignHeuristic phrasealignheuristic = PAH_GROWDIAGFINAL;
    
    unsigned char LEFTCONTEXTSIZE = 0;
    unsigned char RIGHTCONTEXTSIZE = 0;
    
    double alignthreshold = 0.5;
    int pairthreshold = 1;
    
    
    string outputprefix = "";
    
    static struct option long_options[] = {      
       //{"simplelex", no_argument,       &DOSIMPLELEX, 1},
       //{"simpletable", no_argument,       &DOSIMPLETABLE, 1},
       //{"targetfirst", no_argument,       &TARGETFIRST, 1},
       {"moses", no_argument,             &MOSESFORMAT, 1},
       {"stats", no_argument,             &DOSTATS, 1},
       {"null", no_argument,             &EM_NULL, 1}, 
        {"removecontext", no_argument,             &REMOVECONTEXT, 1},                      
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    
    string raw;
    string::size_type pos;
    string a;
    
    
    char c;    
    while ((c = getopt_long(argc, argv, "hd:s:S:t:T:p:P:JDo:O:F:x:X:B:b:l:r:L:NVzEM:v:G:i:23W:a:c:UI:H:",long_options,&option_index)) != -1)
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
        case 'I':
            DOBIDIRECTIONAL = atoi(optarg);
            if ((DOBIDIRECTIONAL != 1) && (DOBIDIRECTIONAL != 2 )) {
                cerr << "Value for -I must be 1 or 2" << endl;
            }
            break;  	
        case 'B':
        	bidirprobthreshold = atof(optarg);        	
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
            LEFTCONTEXTSIZE = atoi(optarg);
            break;     
        case 'r':
            RIGHTCONTEXTSIZE = atoi(optarg);
            break;                 
        case 'v':       
            CONVERGENCE = atof(optarg);
            break;
		case 'L':
            MAXLENGTH = atoi(optarg);
            break;
        case 'M':
            MAXROUNDS = atoi(optarg);
            break;        
        case 'N':
            DOSKIPGRAMS = false;
            break;
        case 'U':
            EXTRACTSKIPGRAMS = true;
            DOTEMPLATES = true;
            DOINSTANCES = true;
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
        case 'z':
        	DONORM = false;
        	break;    
        case 'H':
            a = string(optarg);
            if (a == "growdiagfinal") {
                phrasealignheuristic = PAH_GROWDIAGFINAL;
            } else if (a == "growdiag") {
                phrasealignheuristic = PAH_GROWDIAG;
            } else if (a == "intersection") {
                phrasealignheuristic = PAH_INTERSECTION;
            } else if (a == "union") {
                phrasealignheuristic = PAH_UNION;                
            } else if (a == "s2t") {
                phrasealignheuristic = PAH_S2T;
            } else {
                cerr << "ERROR: Unknown phrase alignment heuristic: '" << a << "'" << endl;
                exit(2);
            }
            break;
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
	
	AlignmentModel * alignmodel = NULL;
    AlignmentModel * reversealignmodel = NULL; 
    SelectivePatternModel * sourcemodel = NULL;
    SelectivePatternModel * targetmodel = NULL;
    GraphFilter filter;
    filter.DOPARENTS = (DOPARENTS || (graphweightfactor > 0));
    filter.DOCHILDREN = DOCHILDREN;
    filter.DOXCOUNT = (DOXCOUNT || (XCOUNTRATIOTHRESHOLD > 0) || (XCOUNTTHRESHOLD > 0));
    filter.DOTEMPLATES = DOTEMPLATES;
    filter.DOINSTANCES = DOINSTANCES;
    filter.DOSKIPUSAGE = DOSKIPUSAGE;
    filter.DOSKIPCONTENT = DOSKIPCONTENT;
    filter.DOSUCCESSORS = DOSUCCESSORS;
    filter.DOPREDECESSORS = DOPREDECESSORS;    
    
    if (!sourcemodelfile.empty()) {
	    cerr << "Loading source model " << sourcemodelfile << endl;   
	    sourcemodel = new SelectivePatternModel(sourcemodelfile, filter, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS || EXTRACTSKIPGRAMS, MINLENGTH, MAXLENGTH, NULL,false, DODEBUG);
	    cerr << "  Loaded " << sourcemodel->types() << " types, " << sourcemodel->tokens() << " tokens" << endl;
     	cerr << "  Ignored " << sourcemodel->ignoredtypes << " types, " << sourcemodel->ignoredoccurrences << " occurrences due to set thresholds" << endl;
	    if (sourcemodel->has_xcount()) {
		    cerr << "  Exclusive count available? YES" << endl;
	    } else {
		    cerr << "  Exclusive count available? NO" << endl;
	    }		
	    if (sourcemodel->has_index()) {
		    cerr << "  Reverse index has " << sourcemodel->reverseindex.size() << " sentences" << endl;
	    } else {
		    cerr << "ERROR: Model " + sourcemodelfile + " contains no indexing information! Unable to align without!" << endl;
		    exit(3);
	    }    
	    if (sourcemodel->has_parents()) {
		    cerr << "  Parent relations available for  " << sourcemodel->rel_subsumption_parents.size() << " patterns" << endl;
	    }
	}
	if (!targetmodelfile.empty()) {
	    cerr << "Loading target model " << targetmodelfile << endl; 		
	    targetmodel = new SelectivePatternModel(targetmodelfile, filter, true, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS || EXTRACTSKIPGRAMS, MINLENGTH, MAXLENGTH, NULL,false, DODEBUG);
	    cerr << "  Loaded " << targetmodel->types() << " types, " << targetmodel->tokens() << " tokens" << endl;
	    cerr << "  Ignored " << targetmodel->ignoredtypes << " types, " << targetmodel->ignoredoccurrences << " occurrences due to set thresholds" << endl;
	    if (targetmodel->has_xcount()) {
		    cerr << "  Exclusive count available? YES" << endl;
	    } else {
		    cerr << "  Exclusive count available? NO" << endl;
	    }
	    if (targetmodel->has_index()) {
		    cerr << "  Reverse index has " << targetmodel->reverseindex.size() << " sentences" << endl;
	    } else {
		    cerr << "ERROR: Model " + targetmodelfile + " contains no indexing information! Unable to align without!" << endl;
		    exit(3);
	    }
	    if (targetmodel->has_parents()) {
		    cerr << "  Parent relations available for  " << targetmodel->rel_subsumption_parents.size() << " patterns" << endl;
	    }		
	}
		
    
	
	if ((DO_EM) || (COOCMODE) || (DOGIZA)) {	
	    //************ BUILD ****************
	    if ((!DOGIZA) && (sourcemodelfile.empty()  || targetmodelfile.empty())) {
      	    cerr << "Error: Specify at least a source model (-s) and target model (-t) to build an alignment model" << endl;
            exit(2);	    
	    } else if ((DOGIZA) && ((DO_EM) || (COOCMODE))) {
	        cerr << "Error: Phrase extraction from GIZA word alignments (-W) can not be combined with EM (-E) or Jaccard (-J)" << endl;
            exit(2);
        } else if ((DOGIZA) && (sourceclassfile.empty() || targetclassfile.empty()) ) {
            cerr << "Error: Phrase extraction from GIZA word alignments (-W) requires source and target classers to be specified (-S and -T)" << endl;
            exit(2);
	    } 
	
	    //no alignment model loaded, build one
	
		cerr << "Configuration: " << endl;
		if (DO_EM) {
			cerr << "\tEM-alignment (-E)" << endl;
		//} else if (DO_EM2) {
        // cerr << "\tAlternative EM-alignment (-2)" << endl;
		//} else if (DO_ITEREM) {
		//	cerr << "\tIterative EM-alignment (-3)" << endl;
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
	
		
        if (outputprefix.empty()) {
            outputprefix = "alignmodel.colibri";
        }
		
		bool usepatternmodels = (!sourcemodelfile.empty() && !targetmodelfile.empty());
		
		if (usepatternmodels) {
		    alignmodel = new AlignmentModel(sourcemodel,targetmodel, LEFTCONTEXTSIZE, RIGHTCONTEXTSIZE, DODEBUG);
		    reversealignmodel = new AlignmentModel(targetmodel,sourcemodel, 0,0, DODEBUG);
		} else {
		    alignmodel = new AlignmentModel(LEFTCONTEXTSIZE, RIGHTCONTEXTSIZE, DODEBUG);
		}    
		
		 				
		bool EM_INIT = true;
		
		if (COOCMODE) {
		    //************ BUILD / COOC ****************
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
				    					
			EM_INIT = false; //no init if EM is done after COOC			
		}		
		if (DO_EM) {
		    //************ BUILD / EM ****************
			cerr << "Computing EM alignment model..." << endl;
			alignmodel->trainEM(MAXROUNDS,  CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT);
			
			if (DONORM) alignmodel->normalize();	
			cerr << "   Found alignment targets for  " << alignmodel->alignmatrix.size() << " source constructions" << endl;
			cerr << "   Total of alignment possibilies in matrix: " << alignmodel->totalsize() << endl;
						
			if ((DOBIDIRECTIONAL) || (!ttableoutfile.empty())) {
				cerr << "Computing reverse alignment model (for bidirectional alignment)..." << endl;
			    reversealignmodel->trainEM(MAXROUNDS, CONVERGENCE, probprunevalue, bestn, EM_NULL, EM_INIT);			
				if (DONORM) reversealignmodel->normalize();
				cerr << "   Found alignment targets for  " << reversealignmodel->alignmatrix.size() << " source constructions" << endl;
				cerr << "   Total of alignment possibilies in matrix: " << reversealignmodel->totalsize() << endl;						
			}	    			    		
		}	
			
		if (DOGIZA) {
		    //************ BUILD / GIZA ****************
			cerr << "Loading source class encoder " << sourceclassfile << endl;
		    ClassEncoder sourceclassencoder = ClassEncoder(sourceclassfile);
    
		    cerr << "Loading target class encoder " << targetclassfile << endl;
    		ClassEncoder targetclassencoder = ClassEncoder(targetclassfile);    
				
		    cerr << "Initialising GIZA++ Word Alignments" << endl;
	        GizaModel gizamodels2t = GizaModel(gizast, &sourceclassencoder, &targetclassencoder);
	        GizaModel gizamodelt2s = GizaModel(gizats, &targetclassencoder, &sourceclassencoder);


		    if (usepatternmodels) {   
		        		    
		        cerr << "Extracting phrases based on GIZA++ Word Alignments and pattern models (semi-supervised)" << endl;
		        int found = alignmodel->extractgizapatterns(gizamodels2t, gizamodelt2s, pairthreshold, coocprunevalue, alignthreshold, DOBIDIRECTIONAL);
		        cerr << "\tFound " << found << " pairs, " << alignmodel->size()  << " source patterns." << endl;
		    } else {
		        cerr << "Extracting phrases based on GIZA++ Word Alignments and heuristics (unsupervised)" << endl;
		        int found = alignmodel->extractgizapatterns_heur(gizamodels2t, gizamodelt2s, phrasealignheuristic, DOBIDIRECTIONAL);
		        cerr << "\tFound " << found << " pairs, " << alignmodel->size()  << " source patterns." << endl;
		    }
		    
		    /*if (DOBIDIRECTIONAL) {
		        GizaModel gizamodels2t = GizaModel(gizast, &sourceclassencoder, &targetclassencoder);
		        GizaModel gizamodelt2s = GizaModel(gizats, &targetclassencoder, &sourceclassencoder);
		        cerr << "Extracting phrases based on GIZA++ Word Alignments (reverse)" << endl;
		        int found = reversealignmodel->extractgizapatterns(gizamodelt2s, gizamodels2t, pairthreshold, coocprunevalue, alignthreshold, DOBIDIRECTIONAL);
		        cerr << "\tFound " << found << " pairs, " << reversealignmodel->size()  << " source patterns." << endl;
		    }*/
		    
		    		  
		}
	} else if (!modelfile.empty()) { //modelfile not empty
	    //************ LOAD ****************
	    cerr << "Loading alignment model..." << endl;
	    if ((sourcemodel != NULL) && (targetmodel != NULL)) {
	        alignmodel = new AlignmentModel(sourcemodel,targetmodel, LEFTCONTEXTSIZE, RIGHTCONTEXTSIZE,DODEBUG);
	        alignmodel->load(modelfile, false, DOSKIPGRAMS, bestn);	    
	        
        } else {
            alignmodel = new AlignmentModel(modelfile, false, DOSKIPGRAMS, bestn, DODEBUG);            
        }  
        cerr << "\tLoaded " << alignmodel->size() << " source patterns." << endl;	        
        
        if (!invmodelfile.empty()) {
        	cerr << "Loading inverse alignment model..." << endl;
	        if ((sourcemodel != NULL) && (targetmodel != NULL)) {
	            reversealignmodel = new AlignmentModel(targetmodel,sourcemodel, 0,0, DODEBUG);
	            reversealignmodel->load(modelfile, false, DOSKIPGRAMS, bestn);	    	        
            } else {
                reversealignmodel = new AlignmentModel(invmodelfile, false, DOSKIPGRAMS, bestn, DODEBUG); 
            }     
            cerr << "\tLoaded " << reversealignmodel->size() << " source patterns for inverse model" << endl;
        }
        
        if (DOSTATS) {
            alignmodel->stats();
        }
        
        if (REMOVECONTEXT) {
            if (!alignmodel->leftsourcecontext && !alignmodel->rightsourcecontext) {
                cerr << "ERROR: Model has no context" << endl;
                exit(2);
            }
            if (outputprefix.empty()) {
                cerr << "ERROR: Specify an explicit output prefix using -o" << endl;
                exit(2);
            }             
            cerr << "Removing context" << endl;
            AlignmentModel * newalignmodel = alignmodel->removecontext();
            if (!outputprefix.empty()) {
                newalignmodel->save(outputprefix);
            }
            exit(0);
        }

        
	} else {
	    cerr << "Error: Don't know what to do.. No model to load or build?" << endl;
	    exit(2);
	}
	
	

    if (EXTRACTSKIPGRAMS) {
        //************ EXTRACTSKIPGRAMS ****************
        cerr << "Extracting skipgrams" << endl;
        int found = alignmodel->extractskipgrams();
        cerr << "\tExtracted " << found << " skipgrams, total " << alignmodel->size() << " source patterns." << endl;
        
        if (DOBIDIRECTIONAL) {
            cerr << "Extracting skipgrams (reverse)" << endl;
            found = reversealignmodel->extractskipgrams();
            cerr << "\tExtracted " << found << " skipgrams, total " << reversealignmodel->size() << " source patterns." << endl;
        }
    }
				
				
	if ((DOBIDIRECTIONAL) && (!DOGIZA)) {
	    //************ INTERSECTION OF MODELS for other methods than giza ****************		    
					
		
		if (DOBIDIRECTIONAL == 1) {
		    cerr << "Computing intersection of both alignment models (joint score)..." << endl;
		    alignmodel->intersect(reversealignmodel, bidirprobthreshold, bestn);
		} else if (DOBIDIRECTIONAL == 2) {			    
		    cerr << "Computing intersection of both alignment models (split scores)..." << endl;
		    alignmodel = new AlignmentModel(*alignmodel, *reversealignmodel);		    
		}	
		cerr << "\t" << alignmodel->size() << " source patterns."	 << endl;
		if (DONORM) {
		    cerr << "Normalizing..." << endl;
		    alignmodel->normalize();
		}    
	}			
     
	// post intersection:		
			
	if (graphweightfactor > 0) {
	    //************ GRAPH WEIGHTING ****************
		cerr << "Weighting based on graph subsumption relations..." << endl;
		const int adjustments = alignmodel->graphalign(*sourcemodel, *targetmodel, graphweightfactor);
		cerr << "   Made " << adjustments << " adjustments" << endl;			
	}

	if (!outputprefix.empty()) {
	    //************ SAVING ****************
	    cerr << "Saving alignment model (" << outputprefix << ")" << endl;
		alignmodel->save(outputprefix);
	}

    if ((!sourceclassfile.empty()) && (!targetclassfile.empty())) {
        //************ DECODING ****************
		cerr << "Loading source class decoder " << sourceclassfile << endl;
		ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);

		cerr << "Loading target class decoder " << targetclassfile << endl;
		ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    	

		cerr << "Decoding..." << endl;
		alignmodel->decode(sourceclassdecoder, targetclassdecoder, &cout, (MOSESFORMAT == 1));
	} else {
	    cerr << "Done without decoding (no -S and -T specified)..." << endl;
	}		

	if (alignmodel != NULL) {
		delete alignmodel;
	}
	
    //EMAlignmentModel(sourcemodel,targetmodel,10000,0.001);    
    
}
