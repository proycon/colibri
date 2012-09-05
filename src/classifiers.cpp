#include <classifiers.h>

using namespace std;
using namespace Timbl;

BuildClassifier::BuildClassifier(const std::string & id, bool appendmode = false, bool exemplarweights = false) {
    ID = id;
    trainfile = string(id + ".train");
    opened = false;        
    featurevectorsize = 0;    
    this->appendmode = appendmode;
    this->exemplarweights = exemplarweights;
}        

BuildClassifier::~BuildClassifier() {
    if (opened) outputfile.close();    
}

void BuildClassifier::addinstance(vector<const string> featurevector, const string & label, double exemplarweight) {
    if (!opened) {
        if (appendmode) {
            outputfile.open(trainfile, ios::app);
        } else {
            outputfile.open(trainfile);
        }
    }

    if (featurevectorsize == 0) {
        featurevectorsize = featurevector.size();
    } else if (featurevector.size() != featurevectorsize) {
        cerr << "INTERNAL ERROR: Passed feature vector of size " << featurevector.size() << ", expected: " << featurevectorsize << endl;
        exit(6); 
    } 
    
    for (vector<const string>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {    
        const string feature = *iter;
        outputfile << feature << ' ';  
    } 
    outputfile << label;
    if (examplarweights) {
        outputfile << ' ' << exemplarweight;
    }
    outputfile << endl;
}

void BuildClassifier::train(const string & timbloptions) {
    ibasefile = string(id + ".ibase");
    TimblAPI * timblexp = new TimblAPI( timbloptions , ID );
    timblexp->Learn(trainfile);   
    timblexp->WriteInstanceBase( ibasefile );
    delete timblexp;    
}
