#include <classifiers.h>

using namespace std;
using namespace Timbl;

Classifier::Classifier(const std::string & _id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool appendmode, bool exemplarweights ) {
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
    vector<string> featurevector_s;
    for (vector<const EncAnyGram *>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram;
        featurevector_s.push_back(anygram->decode(*sourceclassdecoder));        
    }
    const string label_s = label->decode(*targetclassdecoder);
    addinstance(featurevector_s, label_s, exemplarweight);
}

void Classifier::addinstance(vector<string> & featurevector, const string & label, double exemplarweight) {
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
    
    for (vector<string>::iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {    
        const string feature = *iter;
        outputfile << feature << ' ';  
    } 
    outputfile << label;
    if (exemplarweights) {
        outputfile << ' ' << exemplarweight;
    }
    outputfile << endl;
}

void Classifier::train(const string & timbloptions) {
    const string ibasefile = string(id() + ".ibase");
    TimblAPI * timbltrainexp = new TimblAPI( timbloptions , ID );
    timbltrainexp->Learn(trainfile);   
    timbltrainexp->WriteInstanceBase( ibasefile );
    delete timbltrainexp;    
}




NClassifierArray::NClassifierArray(const string & id, int leftcontextsize, int rightcontextsize): ClassifierInterface(id) {
    this->leftcontextsize = leftcontextsize;
    this->rightcontextsize = rightcontextsize;
}

void NClassifierArray::build(const std::string & enctraincorpusfile, AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights) {
    //enctraincorpusfile is ignored
    if (ttable->leftsourcecontext != leftcontextsize) {
        cerr << "Translation table has left context size: " << ttable->leftsourcecontext << ", not " << leftcontextsize << endl;
        exit(3); 
    } else if (ttable->rightsourcecontext != rightcontextsize)  {
        cerr << "Translation table has right context size: " << ttable->rightsourcecontext << ", not " << rightcontextsize << endl;
        exit(3);
    } 
    for (t_contexts::const_iterator iter = ttable->sourcecontexts.begin(); iter != ttable->sourcecontexts.end(); iter++) {
        const EncAnyGram * focus = iter->first;
        if (iter->second.size() > 1) {
            const int n = focus->n();
            stringstream newid;
            newid << this->id() << ".n" << n;
            if (!classifierarray.count(n)) {
                classifierarray[n] = new Classifier(newid.str(), sourceclassdecoder, targetclassdecoder);
                vector<const EncAnyGram *> featurevector;
                for (int i = 0; i < n; i++) {
                    const EncAnyGram * unigram = (const EncAnyGram *) focus->slice(i,1);
                    featurevector.push_back(unigram);                    
                }
                for (t_aligntargets::const_iterator iter2 = ttable->alignmatrix[focus].begin(); iter2 != ttable->alignmatrix[focus].end(); iter2++) {
                    const EncAnyGram * label = iter2->first;
                    //TODO: add exemplar weight iter2->second
                    classifierarray[n]->addinstance(featurevector, label);
                }            
                //cleanup
                for (int i = 0; i < n; i++) {
                    delete featurevector[i];
                }
            }
        }
    }
    
}





