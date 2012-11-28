#include <classifiers.h>
#include <glob.h>
#include "timbl/StringOps.h"
#include <stdio.h>

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

Classifier::Classifier(const std::string & _id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights, bool debug ) {
    //for training
    ID = _id;
    trainfile = string(_id + ".train");
    std::remove(trainfile.c_str()); //remove pre-existing trainfile (no exception if file does not exist)         
    featurevectorsize = 0;    
    this->exemplarweights = exemplarweights;
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassdecoder = targetclassdecoder;
    this->DEBUG = debug;
    testexp = NULL;
    added = false;    
    used = false;
}        

Classifier::Classifier(const std::string & _id, const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, bool loadondemand, bool debug) {
    //for testing
    ID = _id;

    ibasefile = string(_id + ".ibase");
    wgtfile = string(_id + ".ibase.wgt");           
    this->sourceclassdecoder = sourceclassdecoder;
    this->targetclassencoder = targetclassencoder;
    this->DEBUG = debug;
    this->timbloptions = timbloptions;
    added = false;
    
    used = false;
    if (loadondemand) {
        loaded = false;
    } else {
        load();
    }
}

void Classifier::load() {
    if (DEBUG) cerr << "    Loading classifier " << ID << " from " << ibasefile << endl;
    //const string moretimbloptions = "-F Tabbed -i " + ibasefile + " -w " + wgtfile + " " + timbloptions + " +D +vdb";
    const string moretimbloptions = "-F Tabbed " + timbloptions + " +D +vdb -G 0 +vS";
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
    loaded = true;
}

void Classifier::unload() {
    if ((loaded) && (testexp != NULL)) {
        delete testexp;
        testexp = NULL;
        loaded = false; 
    }    
}

void Classifier::remove() {
    unload();
    trainfile = string(ID + ".train");
    ibasefile = string(ID + ".ibase");
    wgtfile = string(ID + ".ibase.wgt");
    const string outfile = string(ID + ".train");      
    std::remove(trainfile.c_str());    
    std::remove(ibasefile.c_str());
    std::remove(wgtfile.c_str());
    std::remove(outfile.c_str());
}

Classifier::~Classifier() {
    unload();
}

void Classifier::addinstance(vector<const EncAnyGram *> & featurevector, const EncAnyGram * label, double exemplarweight) {
    vector<string> featurevector_s;
    for (vector<const EncAnyGram *>::iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {
        const EncAnyGram * anygram = *iter;
        const string feature = anygram->decode(*sourceclassdecoder);
        featurevector_s.push_back(feature);        
    }
    const string label_s = label->decode(*targetclassdecoder);
    addinstance(featurevector_s, label_s, exemplarweight);
}

void Classifier::addinstance(vector<string> & featurevector, const string & label, double exemplarweight) {
    //if (!outputfile.is_open()) {
        //cerr << "Opening (append) " << trainfile << endl;
        outputfile.open(trainfile, ios::app);
        /*} else {
            cerr << "Opening " << trainfile << endl;
            outputfile.open(trainfile);
        }
        */
        if (!outputfile.good()) {
            cerr << "Unable to write to file " << trainfile << endl;
            throw InternalError();
        }
    //}
    

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
    added = true;
    outputfile.close();
}

void Classifier::train(const string & timbloptions) {
    const string ibasefile = string(id() + ".ibase");
    const string moretimbloptions = " -F Tabbed " + timbloptions + " +D +vdb";
    TimblAPI * timbltrainexp = new TimblAPI( moretimbloptions , ID );    
    ifstream * f  = new ifstream(trainfile);
    if (!f->good()) {
        cerr << "Training file not found: " << trainfile << endl;
        throw InternalError();
    }
    f->close();
    delete f;     
    timbltrainexp->Learn(trainfile);   
    timbltrainexp->WriteInstanceBase( ibasefile );
    timbltrainexp->SaveWeights( ibasefile + ".wgt");
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
    /*
concerning p(t|s) only

        C - Set of translation options from classifier
        S - Set of translation options from statistical phrase table

        WEIGHTED, APPEND, IGNORE METHODS:

        |--0--|-------x--------|---(1-x)--|
          C        C ∩ S           S

        probability mass: x:1-x
        |C ∩ S| can be 0 if classifier & statistical model disagree completely)
        WEIGHTED: redistribute within x according to classifier weights

        REPLACE METHOD:

        |----------x--------|---0----|	
          C        C ∩ S           S
          
          
          FILTEREDWEIGHED METHOD:
          [--0--|----x-----|---0--]
             C      C ∩ S     S
    */


    if (!loaded) load();
    used = true;
    
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
    
    //invoke classifier, obtain valuedistribution
    const TargetValue * targetvalue = testexp->Classify(features, valuedistribution, distance);
    

    if (targetvalue == NULL) {
            cerr << "INTERNAL ERROR: Classifier::classify: No targetvalue returned" << endl; 
            throw InternalError();
    } else if (valuedistribution == NULL) {
            cerr << "INTERNAL ERROR: Classifier::classify: No value distribution returned" << endl; 
            throw InternalError();
    }
     
    
    //get amount of scores:
    const double epsilon = -900;
    int scorecount = 0;
    
    
    
    if (scorehandling != SCOREHANDLING_REPLACE) {
        //how long is the score vector?
        t_aligntargets::iterator tmpiter1 = originaltranslationoptions.begin();
        if (tmpiter1 == originaltranslationoptions.end()) {
                cerr << "INTERNAL ERROR: Classifier::classify: No translation options passed!" << endl; 
                throw InternalError();        
        }
        scorecount = tmpiter1->second.size();   
    }
   
    
    //convert valuedistribution to t_aligntargets
    double cis_oldtotal = 0; //|C ∩ S|   
    double cis_total = 0; //|C ∩ S| 
    t_aligntargets result;
    for (ValueDistribution::dist_iterator iter = valuedistribution->begin(); iter != valuedistribution->end(); iter++) {
        const string data = CodeToStr(iter->second->Value()->Name());        
        const double weight = log(iter->second->Weight()); //convert into logprob
        if (DEBUG) cerr << "\t\t\tGot solution \"" << data << "\" with weight " << iter->second->Weight() << " (log=" << weight << ") ";
        const EncAnyGram * target = targetclassencoder->input2anygram(data, false);
        if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_APPEND) || (scorehandling == SCOREHANDLING_IGNORE)) {
            if (originaltranslationoptions.count(target)) {
                //this target occurs in the original statistical model        
                result[target] = originaltranslationoptions[target];               
            } else {
                //this target only occurs in the classifier, ignore
                break;
            }          
        } 
        
        if (scorehandling == SCOREHANDLING_WEIGHED) {
            //ONLY WEIGH FOR p(t|s)            
            //for (int i = 0; i < originaltranslationoptions[target].size(); i++) {
                if (DEBUG) cerr << " [" << result[target][0] << "+" << weight << "=" << result[target][0] + weight << "] ";
                cis_oldtotal += pow(exp(1),result[target][0]); //counting for normalisation later ..prior
                result[target][0] = result[target][0] + weight;
                cis_total += pow(exp(1), result[target][0]);   //counting for normalisation later ..post
            //}                        
        }
        if ((scorehandling == SCOREHANDLING_APPEND) || (scorehandling == SCOREHANDLING_REPLACE)) {
            result[target].push_back(weight);
        }         
        if (DEBUG) cerr << endl; 
    }




    
    
    if ((scorehandling != SCOREHANDLING_REPLACE))  {

        //iterate over original statistical distribution and add translation options not in classifier (S).. renormalise   C | S for WEIGHED methods
            
        for (t_aligntargets::iterator iter = originaltranslationoptions.begin(); iter != originaltranslationoptions.end(); iter++) {
            const EncAnyGram * target = iter->first;
            if ((result.count(target) == 0) && (scorehandling != SCOREHANDLING_FILTEREDWEIGHED)) {
                //translation option is not in classifier: S only, add to results
                if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_IGNORE)) {        
                    result[target] = originaltranslationoptions[target];
                    //s_total += pow(exp(1), iter->second[0]);
                }                 
            } else {
                //translation option is in classifier as well: C | S
                if (scorehandling == SCOREHANDLING_WEIGHED) {
                    //renormalise within C | S only, S remains as is 
                    result[target][0] = log( pow(exp(1), result[target][0]) * (cis_oldtotal/cis_total) ); 
                } else if (scorehandling == SCOREHANDLING_FILTEREDWEIGHED) {
                    result[target][0] = log( pow(exp(1), result[target][0]) * cis_oldtotal) ; 
                }
            }
            
        }  
        
    }
    
    
    if (result.empty()) {
        //should be possible for SCOREHANDLING_FILTEREDWEIGHED
        cerr << "WARNING: NO TRANSLATION OPTIONS AFTER CLASSIFICATION! FALLING BACK TO ORIGINALS FROM STATISTICAL MODEL" << endl;        
        return originaltranslationoptions;        
    }

    
    //no p(t|s) normalisation needed for APPEND method, neither for score 0 nor for last score
    
     
            
    /*        
        for (t_aligntargets::iterator iter = result.begin(); iter != result.end(); iter++) {
            const EncAnyGram * target = iter->first;
            if (originaltranslationoptions.count(target) != 0) {
                //
                result[target] = 
            }
        }    

        //renormalise (only p(t|s) can be normalised) 
        map<int,double> total;
        for (t_aligntargets::iterator iter = result.begin(); iter != result.end(); iter++) {
            for (int i = 0; i < scorecount; i++) {
                if ((i == 0) || ((i == scorecount -1) && (scorehandling == SCOREHANDLING_APPEND))) {
                    if (total.count(i) == 0) total[i] = 0;
                    total[i] += pow(exp(1), iter->second[i]);
                }
            }
        }
        
        for (t_aligntargets::iterator iter = result.begin(); iter != result.end(); iter++) {
            for (int i = 0; i < scorecount; i++) {
                if ((i == 0) || ((i == scorecount -1) && (scorehandling == SCOREHANDLING_APPEND))) {
                    result[iter->first][i] = log(pow(exp(1),iter->second[i]) / total[i]);
                }
            } 
        }
    */   
        
             
            
    
    

    return result;
}


double Classifier::crossvalidate(const std::string & timbloptions) {

    trainfile = string(ID + ".train");
    string outfile = string(ID + ".out");
    unload();
    const string moretimbloptions = "-F Tabbed " + timbloptions + " +D +vdb -G 0 +vS -t leave_one_out";
    if (DEBUG) cerr << "    Instantiating Timbl API: "  << moretimbloptions << endl; 
    testexp = new TimblAPI( moretimbloptions , ID );
    testexp->Learn(trainfile);   
    testexp->Test(trainfile, outfile);
    double a = testexp->GetAccuracy();
    delete testexp;
    testexp = NULL;     
    std::remove(outfile.c_str());
    return a;
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
        classifierarray[n] = new Classifier(nclassifierid, timbloptions,  sourceclassdecoder, targetclassencoder, false, (DEBUG >= 3));
    }    
}

void NClassifierArray::build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) {
    if (ttable->leftsourcecontext != leftcontextsize) {
        cerr << "Translation table has left context size: " << ttable->leftsourcecontext << ", not " << leftcontextsize << endl;
        exit(3); 
    } else if (ttable->rightsourcecontext != rightcontextsize)  {
        cerr << "Translation table has right context size: " << ttable->rightsourcecontext << ", not " << rightcontextsize << endl;
        exit(3);
    } 
    for (t_contexts::const_iterator iter = ttable->sourcecontexts.begin(); iter != ttable->sourcecontexts.end(); iter++) {
        const EncAnyGram * focus = iter->first;
        //cerr << "DEBUG: " << focus->decode(*sourceclassdecoder) << endl;
                        
        if (iter->second.size() >= contextthreshold) { //only use classifier if contextsthreshold is met (by default 1, so it always is)
            const int n = focus->n();
            stringstream newid;
            newid << this->id() << ".n" << n;
            if (!classifierarray.count(n)) {
                classifierarray[n] = new Classifier(newid.str(), sourceclassdecoder, targetclassdecoder, exemplarweights);
            }
            
            unordered_set<const EncAnyGram *> targets; //used to check if enough different targets exists, passing targetthreshold
            for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * withcontext = *iter2;
                for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
                    const EncAnyGram * label = iter3->first;
                    targets.insert(label);
                }
            }            
            if (targets.size() >= targetthreshold) {
                for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                    const EncAnyGram * withcontext = *iter2;
                    const int nwithcontext = withcontext->n();
                    
                    vector<const EncAnyGram *> featurevector;

                    if (singlefocusfeature) {                    
                        //left context
                        for (int i = 0; i < leftcontextsize; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                        
                        featurevector.push_back(focus);
                        
                        //right context
                        for (int i = nwithcontext - rightcontextsize; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }                             
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
                    if (singlefocusfeature) {    
                        for (int i = 0; i < leftcontextsize; i++) {
                            delete featurevector[i];
                        }
                        for (int i = leftcontextsize + 1; i < leftcontextsize + rightcontextsize + 1; i++) {
                            delete featurevector[i];
                        }
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            delete featurevector[i];
                        }
                    }
                }
            }   
            if (classifierarray[n]->empty()) {
                delete classifierarray[n];
                classifierarray.erase(n);
            }             
        }
    }    
    /*for (map<int, Classifier*>::iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        iter->second->flush();
        iter->second->close();
    }*/
}






t_aligntargets NClassifierArray::classify(const EncAnyGram * focus, std::vector<const EncAnyGram *> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const int n = featurevector.size() - leftcontextsize - rightcontextsize;  
    if ((classifierarray.count(n)) && (classifierarray[n] != NULL)) {
        return classifierarray[n]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        throw InternalError();
    } 
}

t_aligntargets NClassifierArray::classify(const EncAnyGram * focus, std::vector<string> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const int n = featurevector.size() - leftcontextsize - rightcontextsize; 
    if ((classifierarray.count(n)) && (classifierarray[n] != NULL)) {
        return classifierarray[n]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        cerr << "INTERNAL ERROR: NClassifierArray::classify invokes classifier " << n << ", but it does not exist" << endl;
        throw InternalError();
    } 
}

void ClassifierInterface::classifyfragments(const EncData & input, AlignmentModel * translationtable, t_sourcefragments & sourcefragments, ScoreHandling scorehandling) {    
     
          
     //decoder will call this
     t_sourcefragments newsourcefragments;
     
     const int maxn = 9; //TODO: make dynamic?
     

     
     vector<pair<const EncAnyGram*, CorpusReference> > tmpsourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn); //will work on context-informed alignment model, always returns patterns without context
     for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = tmpsourcefragments.begin(); iter != tmpsourcefragments.end(); iter++) {
        //returned anygram will be contextless:
        const EncAnyGram * anygram = iter->first;
        const CorpusReference ref = iter->second;
        
        if (DEBUG >= 2) {
            cerr << "\t\t\tClassifying " << anygram->decode(*sourceclassdecoder) << " @" <<  (int) ref.token << ":" << (int) anygram->n() << ": ";     
        }
        
        t_aligntargets translationoptions;
        
        const int contextcount = translationtable->sourcecontexts[anygram].size(); //in how many different contexts does this occur?
        const int n = anygram->n();

        t_aligntargets reftranslationoptions;            
                                    
        //first aggregate original translation options for all training contexts and renormalize.    
        reftranslationoptions = translationtable->sumtranslationoptions(anygram);
        
        bool bypass = false;
        if (contextcount >= contextthreshold) {
            //classify!
            
            //are there enough targets for this source to warrant a classifier?
            if (reftranslationoptions.size() >= targetthreshold) {
                
                //yes
            
                //extract anygram in context for classifier test input
                const EncAnyGram * withcontext = translationtable->addcontext(&input,anygram, (int) ref.token);
                const int nwithcontext = withcontext->n();

                 
                vector<const EncAnyGram *> featurevector;

                if (singlefocusfeature) {                    
                    //left context
                    for (int i = 0; i < translationtable->leftsourcecontext; i++) {
                        const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                        featurevector.push_back(unigram);                    
                    }         
                    
                    featurevector.push_back(anygram);
                    
                    //right context
                    for (int i = nwithcontext - translationtable->rightsourcecontext; i < nwithcontext; i++) {
                        const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                        featurevector.push_back(unigram);                    
                    }         
                } else {
                    for (int i = 0; i < nwithcontext; i++) {
                        const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                        featurevector.push_back(unigram);                    
                    }                             
                }
                
                translationoptions = classify(anygram, featurevector, scorehandling, reftranslationoptions);
                
                
                
                if (DEBUG >= 2) {
                    multimap<double,string> ordered;
                
                    for (t_aligntargets::iterator iter2 = translationoptions.begin(); iter2 != translationoptions.end(); iter2++) {
                        stringstream optionsout;
                        const EncAnyGram * targetkey = iter2->first;                       
                        if (targetkey->isskipgram()) {
                            optionsout << ((const EncSkipGram*) targetkey)->decode(*targetclassdecoder) << " [ ";                            
                        } else {
                            optionsout <<  ((const EncNGram*) targetkey)->decode(*targetclassdecoder) << " [ ";
                        }  
                        
                        if (reftranslationoptions.count(targetkey)) {
                            for (vector<double>::iterator iter3 = reftranslationoptions[targetkey].begin(); iter3 != reftranslationoptions[targetkey].end(); iter3++) {
                                optionsout << *iter3 << " ";
                            }                             
                        }
                        optionsout << "] --> [ ";
                        double totalscore = 0;
                        for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                            optionsout << *iter3 << " ";
                            totalscore += *iter3;
                        }                                                
                        optionsout << "]; ";
                        
                                       
                        ordered.insert( pair<double,string>( -1 * totalscore, optionsout.str() ));                                                                  
                    }
                    
                    for (multimap<double,string>::iterator iter2 = ordered.begin(); iter2 != ordered.end(); iter2++) {
                        cerr << iter2->second;
                    }
                    cerr << endl;                
                }
                
                //cleanup                
                for (int i = 0; i < featurevector.size(); i++) {
                    if (featurevector[i] != anygram) delete featurevector[i];
                }
                                                
                delete withcontext;
            } else {
                if (DEBUG >= 2) cerr << "\t\t\t\tBypassing classifier... number of reference target options is less than set threshold: " << reftranslationoptions.size() << " < " << targetthreshold << endl;
                bypass = true;
            }                    
        } else {
            if (DEBUG >= 2) cerr << "\t\t\t\tBypassing classifier... number of contexts is less than set threshold: " << contextcount << " < " << contextthreshold << endl;
            bypass = true;
        }
        
        
        if (bypass) {
            //bypass classifier; copy from translation table      
        
            for (t_aligntargets::iterator iter = reftranslationoptions.begin(); iter != reftranslationoptions.end(); iter++) {
                const EncAnyGram * target = iter->first;
                
                if (scorehandling == SCOREHANDLING_REPLACE) {
                    translationoptions[target].push_back(reftranslationoptions[target][0]); //fall back to first statistical value    
                } else if (scorehandling == SCOREHANDLING_APPEND) {
                    translationoptions[target] = reftranslationoptions[target];
                    translationoptions[target].push_back(reftranslationoptions[target][0]); //fall back to first statistical value
                } else {
                    translationoptions[target] = reftranslationoptions[target];
                }
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

void writeclassifierconf(const string & ID, ClassifierType type, int contextthreshold, int targetthreshold, bool exemplarweight, bool singlefocusfeature) {
    const string filename = ID + ".classifierconf";
    ofstream * f = new ofstream(ID + ".classifierconf");
    if (f->good()) {
        *f << (int) type << endl;
        *f << contextthreshold << endl;
        *f << targetthreshold << endl;
        *f << exemplarweight << endl;
        *f << singlefocusfeature << endl;
    } else {
        cerr << "ERROR: No classifier configuration found!" << endl;
        throw InternalError();
    }
    f->close();
    delete f;
}


ClassifierType getclassifierconf(const string & ID, int & contextthreshold, int & targetthreshold, bool & exemplarweight, bool & singlefocusfeature) {
    const string filename = ID + ".classifierconf";
    ifstream * f = new ifstream(ID + ".classifierconf");
    string input;
    int type;
    if (f->good()) {        
        *f >> input; 
        istringstream buffer(input);
        buffer >> type;
        *f >> input; 
        istringstream buffer2(input);
        buffer2 >> contextthreshold;
        *f >> input; 
        istringstream buffer3(input);
        buffer3 >> targetthreshold;  
         *f >> input;
        istringstream buffer4(input);
        buffer4 >> exemplarweight;
        *f >> input;
        istringstream buffer5(input);
        buffer5 >> singlefocusfeature;
    }
    f->close();
    delete f;
    return (ClassifierType) type;
}



void MonoClassifier::load( const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG) {    
    vector<string> files = globfiles(ID + ".ibase");
    classifier = new Classifier(ID, timbloptions,  sourceclassdecoder, targetclassencoder, false, (DEBUG >= 3));    
}

void MonoClassifier::build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) {
    if (ttable->leftsourcecontext != leftcontextsize) {
        cerr << "Translation table has left context size: " << ttable->leftsourcecontext << ", not " << leftcontextsize << endl;
        exit(3); 
    } else if (ttable->rightsourcecontext != rightcontextsize)  {
        cerr << "Translation table has right context size: " << ttable->rightsourcecontext << ", not " << rightcontextsize << endl;
        exit(3);
    } 
    classifier = new Classifier(this->id(), sourceclassdecoder, targetclassdecoder, exemplarweights);
    for (t_contexts::const_iterator iter = ttable->sourcecontexts.begin(); iter != ttable->sourcecontexts.end(); iter++) {
        const EncAnyGram * focus = iter->first;
                        
        if (iter->second.size() >= contextthreshold) { //only use classifier if contextsthreshold is met (by default 1, so it always is)
            const int n = focus->n();
            stringstream newid;

            
            unordered_set<const EncAnyGram *> targets; //used to check if enough different targets exists, passing targetthreshold
            for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * withcontext = *iter2;
                for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
                    const EncAnyGram * label = iter3->first;
                    targets.insert(label);
                }
            }            
            if (targets.size() >= targetthreshold) {
                for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                    const EncAnyGram * withcontext = *iter2;
                    const int nwithcontext = withcontext->n();
                    
                    vector<const EncAnyGram *> featurevector;
                    
                    if (singlefocusfeature) {                    
                        //left context
                        for (int i = 0; i < leftcontextsize; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                        
                        featurevector.push_back(focus);
                        
                        //right context
                        for (int i = nwithcontext - rightcontextsize ; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }                             
                    }

                    for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
                        const EncAnyGram * label = iter3->first;
                        
                        if (exemplarweights) {
                            //add exemplar weight         
                            double exemplarweight = iter3->second[0]; //first from score vector, conventionally corresponds to p(t|s) //TODO: Additional methods of weight computation?                    
                            classifier->addinstance(featurevector, label, exemplarweight);
                        } else {
                            classifier->addinstance(featurevector, label);
                        }
                    }                        
                    //cleanup
                    if (singlefocusfeature) {    
                        for (int i = 0; i < leftcontextsize; i++) {
                            delete featurevector[i];
                        }
                        for (int i = leftcontextsize + 1; i < leftcontextsize + rightcontextsize + 1; i++) {
                            delete featurevector[i];
                        }
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            delete featurevector[i];
                        }
                    }
                }
            }   
        }
    }    
}






t_aligntargets MonoClassifier::classify(const EncAnyGram * focus, std::vector<const EncAnyGram *> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {  
        return classifier->classify(featurevector, scorehandling, originaltranslationoptions);
}

t_aligntargets MonoClassifier::classify(const EncAnyGram * focus, std::vector<string> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
        return classifier->classify(featurevector, scorehandling, originaltranslationoptions);
}

void MonoClassifier::train(const string & timbloptions) {
    classifier->train(timbloptions);
}


void ConstructionExperts::train(const string & timbloptions) {
    unsigned int accepted = 0;
    unsigned int discarded = 0;
    unsigned int count = 0;
    unsigned int total = classifierarray.size();
    map<uint64_t,Classifier*>::iterator iter = classifierarray.begin();
    while (iter != classifierarray.end()) {
        double p = ((double) count / (double) total) * 100; 
        count++;        
        double accuracy = 0;
        if (accuracythreshold > 0) {
            cerr << "Cross-validating classifier #" << count << "/" << total << " -- " << p << "% hash=" << iter->first << ": ";
            accuracy =  iter->second->crossvalidate(timbloptions);
            cerr << accuracy << endl;
        }
        if (accuracy >= accuracythreshold) {
            cerr << "Training classifier #" << count << "/" << total << " -- " << p << "% hash=" << iter->first << "... " << endl;                    
            iter->second->train(timbloptions);
            accepted++;
            iter++;
        } else {
            iter->second->remove();
            iter = classifierarray.erase(iter);
            discarded++;            
        }
    }
    cerr << "Training complete: " << accepted << " classifiers built, " << discarded << " discarded due to not making the accuracy threshold." << endl;        
}

void ConstructionExperts::build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) {
    if (ttable->leftsourcecontext != leftcontextsize) {
        cerr << "Translation table has left context size: " << ttable->leftsourcecontext << ", not " << leftcontextsize << endl;
        exit(3); 
    } else if (ttable->rightsourcecontext != rightcontextsize)  {
        cerr << "Translation table has right context size: " << ttable->rightsourcecontext << ", not " << rightcontextsize << endl;
        exit(3);
    } 
    for (t_contexts::const_iterator iter = ttable->sourcecontexts.begin(); iter != ttable->sourcecontexts.end(); iter++) {
        const EncAnyGram * focus = iter->first;                
        if (iter->second.size() >= contextthreshold) { //only use classifier if contextssthreshold is met (by default 1, so it always is)
            stringstream newid;
            const uint64_t hash = focus->hash();
            newid << this->id() << "." << focus->hash();
            if (!classifierarray.count(hash)) {
                classifierarray[hash] = new Classifier(newid.str(), sourceclassdecoder, targetclassdecoder, exemplarweights);
            }
            
            unordered_set<const EncAnyGram *> targets; //used to check if enough different targets exists, passing targetthreshold
            for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * withcontext = *iter2;
                for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
                    const EncAnyGram * label = iter3->first;
                    targets.insert(label);
                }
            }    
            
            if (targets.size() >= targetthreshold) {
                cerr << "Building classifier hash=" << hash << "..." << endl;
                for (unordered_set<const EncAnyGram *>::const_iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                    const EncAnyGram * withcontext = *iter2;
            
                    const int nwithcontext = withcontext->n();
                    vector<const EncAnyGram *> featurevector;
    
                    if (singlefocusfeature) {                    
                        //left context
                        for (int i = 0; i < leftcontextsize; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                        
                        featurevector.push_back(focus);
                        
                        //right context
                        for (int i = nwithcontext - rightcontextsize; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }         
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            const EncAnyGram * unigram = (const EncAnyGram *) withcontext->slice(i,1);
                            featurevector.push_back(unigram);                    
                        }                             
                    }
                    
                                        
                    for (t_aligntargets::const_iterator iter3 = ttable->alignmatrix[withcontext].begin(); iter3 != ttable->alignmatrix[withcontext].end(); iter3++) {
                        const EncAnyGram * label = iter3->first;
                        cerr << "Adding to classifier hash=" << hash << "..." << endl;
                        if (exemplarweights) {
                            //add exemplar weight         
                            double exemplarweight = iter3->second[0]; //first from score vector, conventionally corresponds to p(t|s) //TODO: Additional methods of weight computation?                    
                            classifierarray[hash]->addinstance(featurevector, label, exemplarweight);
                        } else {
                            classifierarray[hash]->addinstance(featurevector, label);
                        }
                    }                        
                    //cleanup
                    if (singlefocusfeature) {    
                        for (int i = 0; i < leftcontextsize; i++) {
                            delete featurevector[i];
                        }
                        for (int i = leftcontextsize + 1; i < leftcontextsize + rightcontextsize + 1; i++) {
                            delete featurevector[i];
                        }
                    } else {
                        for (int i = 0; i < nwithcontext; i++) {
                            delete featurevector[i];
                        }
                    }
                }
            }
            
            if (classifierarray[hash]->empty()) {
                cerr << "Deleting empty classifier hash=" << hash << "..." << endl;
                delete classifierarray[hash];
                classifierarray.erase(hash);
            }  
        }
    }    
    
    /*for (map<uint64_t, Classifier*>::iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        cerr << "Closing classifier hash=" << iter->first << "..." << endl;
        iter->second->flush();
        iter->second->close();
    }*/
}



t_aligntargets ConstructionExperts::classify(const EncAnyGram * focus, std::vector<const EncAnyGram *> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const uint64_t hash = focus->hash();   
    if ((classifierarray.count(hash)) && (classifierarray[hash] != NULL)) {
        return classifierarray[hash]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        if (DEBUG >= 2) cerr << "Classifier " << hash << " does not exist... falling back to statistical model   ";
        if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_IGNORE) || (scorehandling == SCOREHANDLING_FILTEREDWEIGHED)) {
            return originaltranslationoptions;
        } else {
            t_aligntargets translationoptions; 
            for (t_aligntargets::iterator iter = originaltranslationoptions.begin(); iter != originaltranslationoptions.end(); iter++) {
                const EncAnyGram * target = iter->first;
                if (scorehandling == SCOREHANDLING_APPEND) translationoptions[target] = originaltranslationoptions[target];                    
                translationoptions[target].push_back(iter->second[0]); //first score                  
            }
            return translationoptions;
              
        }                 
    } 
}

t_aligntargets ConstructionExperts::classify(const EncAnyGram * focus, std::vector<string> & featurevector,  ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) {
    const uint64_t hash = focus->hash(); 
    if ((classifierarray.count(hash)) && (classifierarray[hash] != NULL)) {
        return classifierarray[hash]->classify(featurevector, scorehandling, originaltranslationoptions);
    } else {
        if (DEBUG >= 2) cerr << "Classifier " << hash << " does not exist... falling back to statistical model    ";
        if ((scorehandling == SCOREHANDLING_WEIGHED) || (scorehandling == SCOREHANDLING_IGNORE) || (scorehandling == SCOREHANDLING_FILTEREDWEIGHED)) {
            return originaltranslationoptions;
        } else {
            t_aligntargets translationoptions; 
            for (t_aligntargets::iterator iter = originaltranslationoptions.begin(); iter != originaltranslationoptions.end(); iter++) {
                const EncAnyGram * target = iter->first;
                if (scorehandling == SCOREHANDLING_APPEND) translationoptions[target] = originaltranslationoptions[target];                    
                translationoptions[target].push_back(iter->second[0]); //first score                  
            }
            return translationoptions;
        }    
    } 
}


void ConstructionExperts::load( const string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG) {    
    vector<string> files = globfiles(ID + ".*.ibase");
    for (vector<string>::iterator iter = files.begin(); iter != files.end(); iter++) {
    
        const string filename = *iter;
    
        stringstream hash_s;
        const string filenamenoext = filename.substr(0,filename.size() - 6);
        
        for (int i = filenamenoext.find_last_of('.')+1; i < filenamenoext.size(); i++) {
            hash_s << filenamenoext[i];            
        }
        
        const uint64_t hash = atoi(hash_s.str().c_str());
        
        cerr << "   Preparing classifier hash=" << hash << " id=" << filenamenoext << endl;   
        classifierarray[hash] = new Classifier(filenamenoext, timbloptions,  sourceclassdecoder, targetclassencoder, true, (DEBUG >= 3));
    }    
}

void ConstructionExperts::reset() {
    //unload unused classifiers
    for (map<uint64_t, Classifier*>::const_iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        if (!iter->second->used) iter->second->unload();
    }
}

NClassifierArray::~NClassifierArray() {
    for (map<int, Classifier*>::const_iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        Classifier * c = iter->second;
        if (c != NULL) delete c;
    }
}

ConstructionExperts::~ConstructionExperts() {
    for (map<uint64_t, Classifier*>::const_iterator iter = classifierarray.begin(); iter != classifierarray.end(); iter++) {
        Classifier * c = iter->second;
        delete c;
    }
}


MonoClassifier::~MonoClassifier() {
    if (classifier != NULL) delete classifier;
}
