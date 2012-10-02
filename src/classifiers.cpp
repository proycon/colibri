#include <classifiers.h>
#include <glob.h>
#include "timbl/StringOps.h"

using namespace std;
using namespace Timbl;



vector<string> globfiles(const string& pat){ //from http://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

Classifier::Classifier(const std::string & _id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool appendmode, bool exemplarweights, bool debug ) {
    //for training
    ID = _id;
    trainfile = string(_id + ".train");        
    featurevectorsize = 0;    
    this->appendmode = appendmode;
    this->exemplarweights = exemplarweights;
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassdecoder = targetclassdecoder;
    this->DEBUG = debug;
}        

Classifier::Classifier(const std::string & _id, const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, bool debug) {
    //for testing
    ID = _id;
    ibasefile = string(_id + ".ibase");
    wgtfile = string(_id + ".ibase.wgt");           
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassencoder = targetclassencoder;
    this->DEBUG = true;
    
    //const string moretimbloptions = "-F Tabbed -i " + ibasefile + " -w " + wgtfile + " " + timbloptions + " +D +vdb";
    const string moretimbloptions = "-F Tabbed " + timbloptions + " +D +vdb -G 0";
    if (DEBUG) cerr << "    Instantiating Timbl API: "  << moretimbloptions << endl; 
    testexp = new TimblAPI( moretimbloptions , ID );
    if (!testexp->Valid()) {
        cerr << "Error Instantiating Timbl API: "  << moretimbloptions << endl;
        throw InternalError();
    }
    if (!testexp->GetInstanceBase(ibasefile)) {
        cerr << "Error getting instance base " << ibasefile << endl;
        throw InternalError();
    }
    if (!testexp->GetWeights(wgtfile)) {
        cerr << "Error getting weights " << wgtfile << endl;
        throw InternalError();
    }
}

Classifier::~Classifier() {
    if (outputfile.is_open()) outputfile.close();    
}

void Classifier::addinstance(vector<const EncAnyGram *> featurevector, const EncAnyGram * label, double exemplarweight) {
    vector<string> featurevector_s;
    for (vector<const EncAnyGram *>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram = *iter;
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
        throw InternalError();
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
    const string moretimbloptions = " -F Tabbed " + timbloptions + " +D +vdb";
    TimblAPI * timbltrainexp = new TimblAPI( moretimbloptions , ID );
    timbltrainexp->Learn(trainfile);   
    timbltrainexp->WriteInstanceBase( ibasefile );
    delete timbltrainexp;    
}

t_aligntargets Classifier::classify(std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    if (DEBUG) cerr << "\t\t\tConverting classifier input to strings..." << endl;
    vector<string> featurevector_s;
    for (vector<const EncAnyGram *>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram = *iter;
        featurevector_s.push_back(anygram->decode(*sourceclassdecoder));        
    }
    return classify(featurevector_s, scorehandling, originaltranslationoptions);
}

t_aligntargets Classifier::classify(std::vector<string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    if (DEBUG) cerr << "\t\t\tClassifier input: " << endl;
    stringstream features_ss;
    for (vector<string>::iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        features_ss << *iter << "\t";
        if (DEBUG) cerr << *iter << "\t"; 
    }
    if (DEBUG) cerr << endl;
    features_ss << "?"; //class dummy    
    
    double distance;
    const string features = features_ss.str();
    const ValueDistribution * valuedistribution;    
    const TargetValue * targetvalue = testexp->Classify(features, valuedistribution, distance);
    

    if (targetvalue == NULL) {
            cerr << "INTERNAL ERROR: Classifier::classify: No targetvalue returned" << endl; 
            throw InternalError();
    } else if (valuedistribution == NULL) {
            cerr << "INTERNAL ERROR: Classifier::classify: No value distribution returned" << endl; 
            throw InternalError();
    }
     
    
    //get amount of scores:
    const double epsilon = -500;
    t_aligntargets::iterator tmpiter1 = originaltranslationoptions.begin();    
    const int scorecount = tmpiter1->second.size();   
   
    
    //convert valuedistribution to t_aligntargets    
    t_aligntargets result;
    for (ValueDistribution::dist_iterator iter = valuedistribution->begin(); iter != valuedistribution->end(); iter++) {
        const string data = CodeToStr(iter->second->Value()->Name());        
        const double weight = log(iter->second->Weight()); //convert into logprob
        if (DEBUG) cerr << "\t\t\tGot solution \"" << data << "\" with weight " << iter->second->Weight() << endl;
        const EncAnyGram * target = targetclassencoder->input2anygram(data, false);
        if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_APPEND)) {
            if (originaltranslationoptions.count(target)) {
                result[target] = originaltranslationoptions[target];               
            } else {
                //translation option did not exist yet:
                for (int i = 0; i < scorecount; i++) {
                        result[target].push_back(epsilon);
                }
            }          
        } 
        if (scorehandling == SCOREHANDLING_WEIGHED) {
            for (int i = 0; i < originaltranslationoptions[target].size(); i++) {
                result[target][i] = originaltranslationoptions[target][i] + weight;
            }                        
        }
        if ((scorehandling == SCOREHANDLING_APPEND) || (scorehandling == SCOREHANDLING_REPLACE)) {
            result[target].push_back(weight);
        }         
    }

    //note: any targets not present in classifier output will be pruned!
        
    return result;
}


NClassifierArray::NClassifierArray(const string & id, int leftcontextsize, int rightcontextsize): ClassifierInterface(id) {
    this->leftcontextsize = leftcontextsize;
    this->rightcontextsize = rightcontextsize;
}

void NClassifierArray::load( const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG) {    
    vector<string> files = globfiles(ID + ".n*.ibase");
    for (vector<string>::iterator iter = files.begin(); iter != files.end(); iter++) {
    
        const string filename = *iter;
        
        //grab n
        stringstream nss;
        for (int i = filename.size() - 7; i >= 0; i--) {
            if (filename[i] == 'n') break;
            nss << filename[i];            
        }
        const int n = atoi(nss.str().c_str());
        
        const string nclassifierid = filename.substr(0, filename.size() - 6);
        cerr << "   Loading classifier n=" << n << " id=" << nclassifierid << endl;   
        classifierarray[n] = new Classifier(nclassifierid, timbloptions,  sourceclassdecoder, targetclassencoder, (DEBUG >= 3));
    }    
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
                classifierarray[n] = new Classifier(newid.str(), sourceclassdecoder, targetclassdecoder, false, exemplarweights);
            }
            for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * withcontext = *iter2;
        
                const int nwithcontext = withcontext->n();
                vector<const EncAnyGram *> featurevector;
                for (int i = 0; i < nwithcontext; i++) {
                    const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                    featurevector.push_back(unigram);                    
                }                
                for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
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



t_aligntargets NClassifierArray::classify(std::vector<const EncAnyGram *> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const int n = featurevector.size() - leftcontextsize - rightcontextsize;  
    if ((classifierarray.count(n)) && (classifierarray[n] != NULL)) {
        cerr << "Invoking classifier " << n << endl;
        return classifierarray[n]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        throw InternalError();
    } 
}

t_aligntargets NClassifierArray::classify(std::vector<string> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const int n = featurevector.size() - leftcontextsize - rightcontextsize; 
    if ((classifierarray.count(n)) && (classifierarray[n] != NULL)) {
        cerr << "Invoking classifier " << n << endl;
        return classifierarray[n]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        throw InternalError();
    } 
}

void NClassifierArray::classifyfragments(const EncData & input, AlignmentModel * translationtable, t_sourcefragments & sourcefragments, ScoreHandling scorehandling) {    
     
          
     //decoder will call this
     t_sourcefragments newsourcefragments;
     
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
            t_aligntargets originaltranslationoptions = translationtable->alignmatrix[anygramwithcontext];
            
            for (t_aligntargets::iterator iter = originaltranslationoptions.begin(); iter != originaltranslationoptions.end(); iter++) {
                const EncAnyGram * target = iter->first;
                const double weight = 0; //== log(1.0)
                
                if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_APPEND)) {
                    if (originaltranslationoptions.count(target) == 0) {
                        cerr << "INTERNAL ERROR: Classifier::classify: Original translation option not found" << endl; 
                        throw InternalError();
                    }
                    translationoptions[target] = originaltranslationoptions[target];          
                } 
                if (scorehandling == SCOREHANDLING_WEIGHED) {
                    for (int i = 0; i <= originaltranslationoptions[target].size(); i++) {
                        translationoptions[target][i] = originaltranslationoptions[target][i] + weight;
                    }                        
                }
                if ((scorehandling == SCOREHANDLING_APPEND) || (scorehandling == SCOREHANDLING_REPLACE)) {
                    translationoptions[target].push_back(weight);
                }
                               
                                
            } 
             
        } else {
            //more context possible? classify!
            
            t_aligntargets reftranslationoptions;
            
            if (scorehandling != SCOREHANDLING_REPLACE) {                
                //first aggregate original translation options for all training contexts and renormalize.                
                for (unordered_set<const EncAnyGram *>::iterator iter = translationtable->sourcecontexts[anygram].begin(); iter != translationtable->sourcecontexts[anygram].end(); iter++) {
                    const EncAnyGram * sourceincontext = *iter;
                    for (t_aligntargets::iterator targetiter = translationtable->alignmatrix[sourceincontext].begin(); targetiter != translationtable->alignmatrix[sourceincontext].end(); targetiter++) {
                        const EncAnyGram * target = targetiter->first;
                        if (reftranslationoptions.count(target) == 0) {
                            reftranslationoptions[target] = targetiter->second;
                            for (int i = 0; i < targetiter->second.size(); i++) {
                                reftranslationoptions[target][i] = reftranslationoptions[target][i] / contextcount;  
                            } 
                        } else {
                            for (int i = 0; i < targetiter->second.size(); i++) {
                                reftranslationoptions[target][i] += targetiter->second[i] / contextcount;  
                            }
                        }
                    }
                }
            } 
                        
            //extract anygram in context for classifier test input
            const EncAnyGram * withcontext = translationtable->addcontext(&input,anygram, (int) ref.token);
            const int nwithcontext = withcontext->n();

             
            vector<const EncAnyGram *> featurevector;
            for (int i = 0; i < nwithcontext; i++) {
                const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                featurevector.push_back(unigram);                    
            }
            translationoptions = classify(featurevector, scorehandling, reftranslationoptions);
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


