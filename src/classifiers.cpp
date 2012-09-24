#include <classifiers.h>

using namespace std;
using namespace Timbl;

Classifier::Classifier(const std::string & _id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool appendmode = false, bool exemplarweights = false) {
    ID = _id;
    trainfile = string(_id + ".train");
    opened = false;        
    featurevectorsize = 0;    
    this->appendmode = appendmode;
    this->exemplarweights = exemplarweights;
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassdecoder = targetclassdecoder;
}        

Classifier::~Classifier() {
    if (opened) outputfile.close();    
}

void Classifier::addinstance(vector<const EncAnyGram *> featurevector, const EncAnyGram * label, double exemplarweight) {
    vector<const string> featurevector_s;
    for (vector<const EncAnyGram *>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram;
        featurevector_s.push_back(anygram->decode(sourceclassdecoder));        
    }
    const string label_s = label->decoder(targetclassdecoder);
    addinstance(featurevector_s, label_s, exemplarweight);
}

void Classifier::addinstance(vector<const string> & featurevector, const string & label, double exemplarweight) {
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

void Classifier::train(const string & timbloptions) {
    ibasefile = string(id() + ".ibase");
    TimblAPI * timbltrainexp = new TimblAPI( timbloptions , ID );
    timbltrainexp->Learn(trainfile);   
    timbltrainexp->WriteInstanceBase( ibasefile );
    delete timbltrainexp;    
}




NClassifierArray::NClassifierArray(const string & id, int maxn, int leftcontextsize, int rightcontextsize) {
    this->leftcontextsize = leftcontextsize;
    this->rightcontextsize = rightcontextsize;
    this->maxn = maxn;
}

NClassifierArray::build(const TranslationTable * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights) {
    if (ttable->leftcontextsize != leftcontextsize) {
        cerr << "Translation table has left context size: " << ttable->leftcontextsize << ", not " << leftcontextsize << endl;
        exit(3); 
    } else {
        cerr << "Translation table has right context size: " << ttable->rightcontextsize << ", not " << rightcontextsize << endl;
        exit(3);
    } 
    for (t_contexts::iterator iter = ttable->sourcecontexts.begin(); ttable->sourcecontexts.end(); iter++) {
        const EncAnyGram * focus = iter->first;
        const int n = focus->n();
        stringstream newid;
        newid << this->id() << ".n" << n;
        if (!classifierarray.count(n) {
            classifierarray[n] = new Classifier(newid.str(), sourceclassdecoder, targetclassdecoder);
            vector<const EncAnyGram *> featurevector;
            for (int i = 0; i < n; i++) {
                const EncNGram * unigram = gettoken(i);
                featurevector.push_back(unigram);
            }
            for (t_aligntargets::iterator iter = alignmatrix[focus].begin(); alignmatrix[focus].end(); iter++) {
                const EncAnyGram * label = iter->first;
                //TODO: add exemplar weight iter->second
                classifierarray[n]->addinstance(featurevector, label);
            }            
            //cleanup
            for (int i = 0; i < n; i++) {
                delete featurevector[i];
            }
        }
    }
    
}





