#include "patternmodel.h"
#include <getopt.h>

using namespace std;

void usage() {
    cerr << "Usage: querier -d model -c class-file" << endl;
    cerr << "Options:" << endl;    
    cerr << "\t-d model                 Load an existing model" << endl; //TODO
    cerr << "\t-c class-file            The class file (for encoding and decoding)" << endl; //TODO
    //cerr << "\t-B                       Do a bi-directional alignment and compute intersection of results" << endl; //TODO
    cerr << "\t-l n                     Minimum N length" << endl;
    cerr << "\t-L n                     Maximum N length" << endl;
    cerr << "\t-N                       No skip-grams" << endl;  
    cerr << "\t-o occurence-threshold   Consider only patterns occuring more than specified (absolute occurrence). Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-F freq-threshold        Consider only patterns occuring more than specified (relative frequency of all patterns).  Note: The model you load may already be pruned up to a certain value, only higher numbers have effect." << endl;
    cerr << "\t-x xcount-threshold      Consider only patterns with an *exclusive* count over this threshold" << endl;
    cerr << "\t-X xcount-ratio          Consider only patterns with an *exclusivity ratio* over this threshold (between 0.0 [not exclusive] and 1.0 [entirely exclusive])" << endl;
}


int main( int argc, char *argv[] ) {
	string modelfile = "";
	string classfile = "";
    int COUNTTHRESHOLD = 0;
    int FREQTHRESHOLD = 0; 
    double XCOUNTTHRESHOLD = 0;
    double XCOUNTRATIOTHRESHOLD = 0;
    int MINLENGTH = 0;
    int MAXLENGTH = 99;
    bool DOSKIPGRAMS = true;
	bool DEBUG = false;
	
    char c;    
    while ((c = getopt(argc, argv, "d:c:l:L:No:F:x:X:D")) != -1)
        switch (c)
        {
        case 'd':
            modelfile = optarg;
            break;        
        case 'D':
        	DEBUG = true;
        	break;
        case 'c':
            classfile = optarg;
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
        
    if (modelfile.empty()  || classfile.empty()) {
  	    cerr << "Error: Specify a model and a class file" << endl;
        usage();
        exit(2);
    }
    
    cerr << "Configuration: " << endl;
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
	
	cerr << "Loading model " << modelfile << endl;
	GraphFilter filter;
	filter.DOXCOUNT = ((XCOUNTRATIOTHRESHOLD > 0) || (XCOUNTTHRESHOLD > 0));
    SelectivePatternModel model = SelectivePatternModel(modelfile, filter,false, true, COUNTTHRESHOLD, FREQTHRESHOLD, XCOUNTRATIOTHRESHOLD, XCOUNTTHRESHOLD, DOSKIPGRAMS, MINLENGTH, MAXLENGTH,NULL, false ,DEBUG);
    cerr << "  Loaded " << model.types() << " types, " << model.tokens() << " tokens" << endl;
    cerr << "  Ignored " << model.ignoredtypes << " types, " << model.ignoredoccurrences << " occurrences due to set thresholds" << endl;
    if (model.has_xcount()) {
    	cerr << "  Exclusive count available? YES" << endl;
    } else {
    	cerr << "  Exclusive count available? NO" << endl;
    }
    if (model.has_index()) {
    	cerr << "  Index available? YES" <<  endl;    	
    } else {
	    cerr << "  Index available? NO" <<  endl;
    }    

	cerr << "Loading class decoder " << classfile << endl;
	ClassDecoder classdecoder = ClassDecoder(classfile);
	cerr << "Loading class encoder " << classfile << endl;
	ClassEncoder classencoder = ClassEncoder(classfile);
	
	cerr << "Starting query mode:" << endl;
	model.querier(classencoder, classdecoder, false, true, model.getminn(), model.getmaxn());	
}
