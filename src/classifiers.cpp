#include <classifiers.h>

using namespace std;
using namespace Timbl;

Classifier::Classifier(const std::string & _id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool appendmode, bool exemplarweights ) {
    //for training
    ID = _id;
    trainfile = string(_id + ".train");        
    featurevectorsize = 0;    
    this->appendmode = appendmode;
    this->exemplarweights = exemplarweights;
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassdecoder = targetclassdecoder;
}        

Classifier::Classifier(const std::string & _id, const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder) {
    //for testing
    ID = _id;
    ibasefile = string(_id + ".ibase");        
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassencoder = targetclassencoder;
    testexp = new TimblAPI( timbloptions , ID );
}

Classifier::~Classifier() {
    if (outputfile.is_open()) outputfile.close();    
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
    if (!outputfile.is_open()) {
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
        outputfile << feature << '\t';  
    } 
    outputfile << label;
    if (exemplarweights) {
        outputfile << '\t' << exemplarweight;
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

t_aligntargets Classifier::classify(std::vector<const EncAnyGram *> & featurevector) {
    vector<string> featurevector_s;
    for (vector<const EncAnyGram *>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram;
        featurevector_s.push_back(anygram->decode(*sourceclassdecoder));        
    }
    return classify(featurevector_s);
}

t_aligntargets Classifier::classify(std::vector<string> & featurevector) {
    stringstream features_ss;
    for (vector<string>::iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        features_ss << *iter << "\t";
    }
    features_ss << "?"; //class dummy    
    const ValueDistribution * valuedistribution;
    double distance;
    const string features = features_ss.str();
    testexp->Classify(features, valuedistribution, distance);
    
    //convert valuedistribution to t_aligntargets
    t_aligntargets result;
    for (ValueDistribution::dist_iterator iter = valuedistribution->begin(); iter != valuedistribution->end(); iter++) {
        const string data = iter->second->Value()->Name();
        const double weight = iter->second->Weight();
        const EncAnyGram * target = targetclassencoder->input2anygram(data, false);
        result[target].push_back(weight);         
    }
    return result;
}


NClassifierArray::NClassifierArray(const string & id, int leftcontextsize, int rightcontextsize): ClassifierInterface(id) {
    this->leftcontextsize = leftcontextsize;
    this->rightcontextsize = rightcontextsize;
}

void NClassifierArray::build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights) {
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
            }
            for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * withcontext = *iter2;
        
                const int nwithcontext = withcontext->n();
                vector<const EncAnyGram *> featurevector;
                for (int i = 0; i < nwithcontext; i++) {
                    const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                    featurevector.push_back(unigram);                    
                }                
                for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[focus].begin(); iter3 != ttable->alignmatrix[focus].end(); iter3++) {
                    const EncAnyGram * label = iter3->first;                                                            
                    if (exemplarweights) {
                        //add exemplar weight         
                        double exemplarweight = iter3->second[0]; //first from score vector, conventionally corresponds to p(t|s) //TODO: Additional methods of weight computation?                    
                        classifierarray[n]->addinstance(featurevector, label, exemplarweight);
                    } else {
                        classifierarray[n]->addinstance(featurevector, label);
                    }
                }                        
                //cleanup
                for (int i = 0; i < nwithcontext; i++) {
                    delete featurevector[i];
                }
            }        
        }
    }    
}



t_aligntargets NClassifierArray::classify(std::vector<const EncAnyGram *> & featurevector) {
    const int n = featurevector.size() - 1 - leftcontextsize - rightcontextsize; // - 1 for dummy
    if (classifierarray.count(n)) {
        return classifierarray[n]->classify(featurevector);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        exit(6);
    } 
}

t_aligntargets NClassifierArray::classify(std::vector<string> & featurevector) {
    const int n = featurevector.size() - 1 - leftcontextsize - rightcontextsize; // - 1 for dummy
    if (classifierarray.count(n)) {
        return classifierarray[n]->classify(featurevector);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        exit(6);
    } 
}

void NClassifierArray::classifyfragments(const EncData & input, AlignmentModel * translationtable, t_sourcefragments sourcefragments) {    
     //decoder will call this, sourcefragments and newttable will be filled for decoder
     
     const int maxn = 9; //TODO: make dynamic?
     
     vector<pair<const EncAnyGram*, CorpusReference> > tmpsourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn); //will work on context-informed alignment model, always returns patterns without context
     for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = tmpsourcefragments.begin(); iter != tmpsourcefragments.end(); iter++) {
        //returned anygram will be contextless:
        const EncAnyGram * anygram = iter->first;
        const CorpusReference ref = iter->second;
        
        t_aligntargets translationoptions;
        
        const int contextcount = translationtable->sourcecontexts[anygram].size(); //in how many different contexts does this occur?
        const int n = anygram->n();
        
        if (contextcount == 1) {
            //only one? no need for classifier, just copy from translation table
            const EncAnyGram * anygramwithcontext = *(translationtable->sourcecontexts[anygram].begin());
            translationoptions = translationtable->alignmatrix[anygramwithcontext]; 
        } else {
            //more context possible? classify!

            //extract anygram in context for classifier test input
            const EncAnyGram * withcontext = translationtable->addcontext(&input,anygram, (int) ref.token);
            const int nwithcontext = withcontext->n();
             
            vector<const EncAnyGram *> featurevector;
            for (int i = 0; i < nwithcontext; i++) {
                const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                featurevector.push_back(unigram);                    
            }
            translationoptions = classify(featurevector);
            //cleanup
            for (int i = 0; i < nwithcontext; i++) {
                delete featurevector[i];
            }            
        }

        sourcefragments.push_back(SourceFragmentData(anygram, ref, translationoptions));
     }     
}


void NClassifierArray::train(const string & timbloptions) {
    for (map<int,Classifier*>::iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        cerr << "Training classifier n=" << iter->first << "..." << endl;
        iter->second->train(timbloptions);
    }

    
}
