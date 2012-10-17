#include <getopt.h>
#include "alignmodel.h"
#include "classifiers.h"

using namespace std;

void usage() {
    cerr << "Usage: trainclassifiers [-N|-X] -d alignmentmodel -S source-class-file -T target-class-file" << endl;
    cerr << "Classifier types: (pick one)" << endl;
    cerr << " -N           N-Classifier Array, one classifier per pattern size group" << endl;
    cerr << " -X           Construction experts, one classifier per construction" << endl;
    cerr << "Options:" << endl;
    cerr << " -C [id]      Classifier output prefix. The decoder takes this same ID to load your classifier." << endl;
    cerr << " -c [int]     Context threshold. Only create a classifier when at least this many different contexts exist. Defaults to 1." << endl;
    cerr << " -x           disable exemplar weighting" << endl;
    cerr << " -O [options] Timbl options" << endl;
    //cerr << "\t-C number                 Classifier mode" << endl;
    //cerr << "\t   1 - Local context with Classifier Array" << endl;
}

int main( int argc, char *argv[] ) {
    string sourceclassfile = "";
    string targetclassfile = "";
    string modelfile="";
    string outputprefix="classifier";
    bool exemplarweights = true;
    
    static struct option long_options[] = {                        
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    string timbloptions;
    if (exemplarweights) {
        timbloptions = "-a 1 -s -F Tabbed";
    } else {
        timbloptions = "-a 1 -F Tabbed";
    }
    
    
    ClassifierType mode = CLASSIFIERTYPE_NONE;
    
    int contextthreshold = 1;
    
    char c;    
    while ((c = getopt_long(argc, argv, "hd:S:T:C:xO:XNc:",long_options,&option_index)) != -1) {
        switch (c) {
        case 0:
            if (long_options[option_index].flag != 0)
               break;
        case 'd':
        	modelfile = optarg;
        	break;
        case 'h':
        	usage();
        	exit(0);
        case 'S':
            sourceclassfile = optarg;
            break;
        case 'T':
            targetclassfile = optarg;
            break;    
        case 'C':
            outputprefix = optarg;
            break;
        case 'X':
            mode = CLASSIFIERTYPE_CONSTRUCTIONEXPERTS;
            break;
        case 'N':
            mode = CLASSIFIERTYPE_NARRAY;
            break;
        case 'O':
            if (exemplarweights) {
                timbloptions = "-s -F Tabbed " + std::string(optarg);
            } else {
                timbloptions = "-F Tabbed " + std::string(optarg);
            }
            break;
        case 'c':
            contextthreshold = atoi(optarg);
            break; 
        case 'x':
            exemplarweights = false;
            break; 	
        }
    }
    
    if ((modelfile.empty()) || (sourceclassfile.empty()) || (targetclassfile.empty()) || (mode == CLASSIFIERTYPE_NONE)) {
        usage();
        exit(2);
    }
    
    cerr << "Loading source class decoder " << sourceclassfile << endl;
	ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);

	cerr << "Loading target class decoder " << targetclassfile << endl;
	ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);   
		
    cerr << "Loading alignment model " << modelfile << endl;
    AlignmentModel alignmodel = AlignmentModel(modelfile,false,true,0, false);
    
    if (mode == CLASSIFIERTYPE_NARRAY) {
    
        cerr << "Building N-Array classifiers" << endl;
        NClassifierArray classifiers = NClassifierArray(outputprefix, alignmodel.leftsourcecontext, alignmodel.rightsourcecontext, contextthreshold);    
        classifiers.build(&alignmodel, &sourceclassdecoder, &targetclassdecoder, exemplarweights);

        cerr << "Training classifiers" << endl;
        cerr << "   Timbl options: " << timbloptions << endl;
        classifiers.train(timbloptions);
    
    } else if (mode == CLASSIFIERTYPE_CONSTRUCTIONEXPERTS) {
    
        
        cerr << "Building construction expert classifiers" << endl;
        ConstructionExperts classifiers = ConstructionExperts(outputprefix, alignmodel.leftsourcecontext, alignmodel.rightsourcecontext, contextthreshold);    
        classifiers.build(&alignmodel, &sourceclassdecoder, &targetclassdecoder, exemplarweights);
        
        cerr << "Training classifiers" << endl;
        cerr << "   Timbl options: " << timbloptions << endl;
        classifiers.train(timbloptions);
    
    }
    
    writeclassifierconf(outputprefix, mode, contextthreshold, exemplarweights);
    

    
    cerr << "All done" << endl;
}
    
