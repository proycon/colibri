#include <getopt.h>
#include "alignmodel.h"
#include "classifiers.h"

using namespace std;

void usage() {
    cerr << "Usage: trainclassifiers -d alignmentmodel -S source-class-file -T target-class-file" << endl;
    cerr << "Options:" << endl;
    cerr << " -C [id]      Classifier output prefix" << endl;
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
        timbloptions = "-a 1 -s";
    } else {
        timbloptions = "-a 1";
    }
    
    
    
    char c;    
    while ((c = getopt_long(argc, argv, "hd:S:T:C:xt:",long_options,&option_index)) != -1) {
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
        case 'O':
            timbloptions = optarg;
            break;
        case 'x':
            exemplarweights = false;
            break; 	
        }
    }
    
    if ((modelfile.empty()) || (sourceclassfile.empty()) || (targetclassfile.empty())) {
        usage();
        exit(2);
    }
    
    cerr << "Loading source class decoder " << sourceclassfile << endl;
	ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);

	cerr << "Loading target class decoder " << targetclassfile << endl;
	ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);   
		
    cerr << "Loading alignment model " << modelfile << endl;
    AlignmentModel alignmodel = AlignmentModel(modelfile,false,true,0, false);
    
    cerr << "Building classifiers" << endl;
    NClassifierArray classifiers = NClassifierArray(outputprefix, alignmodel.leftsourcecontext, alignmodel.rightsourcecontext);    
    classifiers.build(&alignmodel, &sourceclassdecoder, &targetclassdecoder, exemplarweights);
    
    cerr << "Training classifiers" << endl;
    cerr << "   Timbl options: " << timbloptions << endl;
    classifiers.train(timbloptions);
    
    cerr << "All done" << endl;
}
    
