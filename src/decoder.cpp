#include <decoder.h>
#include <algorithm> 
#include <getopt.h>
#include <iostream>
#include <sstream>


using namespace std;

unsigned char UNKNOWNCLASS = 2;
unsigned char BOSCLASS = 3;
unsigned char EOSCLASS = 4;

const EncNGram UNKNOWNUNIGRAM = EncNGram(&UNKNOWNCLASS,1);

StackDecoder::StackDecoder(const EncData & input, AlignmentModel * translationtable, LanguageModel * lm, int stacksize, double prunethreshold, vector<double> tweights, double dweight, double lweight, int dlimit, int maxn, int debug, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, ClassifierInterface * classifier, ScoreHandling scorehandling, bool globalstats) {
        this->input = input;
        this->inputlength = input.length();
        this->translationtable = translationtable;
        this->lm = lm;
        this->stacksize = stacksize;
        this->prunethreshold = prunethreshold;
        this->DEBUG = debug;
        this->sourceclassdecoder = sourceclassdecoder;
        this->targetclassdecoder = targetclassdecoder;        
        this->globalstats = globalstats;
        
        //init stacks
        for (unsigned int i = 0; i <= inputlength; i++) {
            stacks.push_back( Stack(this, i, stacksize, prunethreshold) );
            gappystacks.push_back( Stack(this, i, stacksize, prunethreshold) );
        }

        if (DEBUG >= 3) cerr << "Gathering source fragments:" << endl;        
        
        if (classifier != NULL) {
            //Use classifier
            if (DEBUG >= 3) cerr << "  Calling classifier" << endl;
            classifier->classifyfragments(input, translationtable, sourcefragments, scorehandling);
        } else {            
            //Collect source fragments and translation options straight from translation table
            
            vector<pair<const EncAnyGram*, CorpusReference> > tmpsourcefragments;  
            tmpsourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn);
            if (DEBUG >= 3) cerr << "  " << tmpsourcefragments.size() << " source-fragments found in translation table" << endl;
            for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = tmpsourcefragments.begin(); iter != tmpsourcefragments.end(); iter++) {
                const EncAnyGram * sourcekey = iter->first;
                if (translationtable->leftsourcecontext || translationtable->rightsourcecontext) {
                    //translation table context information, aggregate scores  
                    t_aligntargets translationoptions = translationtable->sumtranslationoptions(sourcekey); 
                    sourcefragments.push_back(SourceFragmentData(iter->first, iter->second, translationoptions));
                } else {
                    //no context information, do simply copy:
                    t_aligntargets translationoptions;
                    for (t_aligntargets::iterator iter2 = translationtable->alignmatrix[sourcekey].begin(); iter2 != translationtable->alignmatrix[sourcekey].end(); iter2++) {
                            const EncAnyGram * targetgram = iter2->first;
                            translationoptions[targetgram] = iter2->second;
                    }
                    sourcefragments.push_back(SourceFragmentData(iter->first, iter->second, translationoptions)); 
                }
            }
            if (DEBUG >= 3) cerr << "  " << sourcefragments.size() << " source-fragments registered" << endl;
        }
        
        
        
        
        //Build a coverage mask, this will be used to check if their are words uncoverable by translations, these will be added as unknown words 
        std::vector<bool> inputcoveragemask;
        inputcoveragemask.reserve(inputlength);
        for (unsigned int i = 0; i < inputlength; i++) inputcoveragemask.push_back(false);
        
        
                        
        for (t_sourcefragments::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {           
            const EncAnyGram * sourcekey = iter->sourcefragment;
            int n = sourcekey->n();
            for (int j = iter->ref.token; j < iter->ref.token + n; j++) inputcoveragemask[j] = true;
                        
                        
            //Output translation options
            if ((DEBUG >= 3) && (sourceclassdecoder != NULL) && (targetclassdecoder != NULL)) {
                if (sourcekey->isskipgram()) {
                    cerr << "\t" << (int) iter->ref.token << ':' << n << " -- " << ((const EncSkipGram*) sourcekey)->decode(*sourceclassdecoder) << " ==> ";                    
                } else {
                    cerr << "\t" << (int) iter->ref.token << ':' << n << " -- " << ((const EncNGram*) sourcekey)->decode(*sourceclassdecoder) << " ==> ";
                }
                if (iter->translationoptions.empty()) {
                    cerr << endl;
                    cerr << "ERROR: No translation options found!!!!" << endl;
                    exit(6);
                } else {
                    for (t_aligntargets::iterator iter2 = iter->translationoptions.begin(); iter2 != iter->translationoptions.end(); iter2++) {
                        const EncAnyGram * targetkey = iter2->first;
                        if (targetkey->isskipgram()) {
                            cerr << " SKIPGRAM!!! "; //debug
                            cerr << ((const EncSkipGram*) targetkey)->decode(*targetclassdecoder) << " [ ";                            
                        } else {
                            cerr <<  ((const EncNGram*) targetkey)->decode(*targetclassdecoder) << " [ ";
                        }                        
                        for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                            cerr << *iter3 << " ";
                        }
                        if (!targetkey->isskipgram()) {
                           cerr << "LM=" << lm->score((const EncNGram*) targetkey) << " ";
                        }
                        cerr << "]; ";                        
                    }
                }
                cerr << endl;            
            }
        }
        
        //Check for uncoverable words               
        for (unsigned int i = 0; i < inputlength; i++) {
            if (!inputcoveragemask[i]) {
                //found one
                
                bool keepunigram = false;
                EncNGram * unigram = input.slice(i, 1);
                const string word = unigram->decode(*sourceclassdecoder);
                cerr << "NOTICE: UNTRANSLATABLE WORD: '" << word << "' (adding)" << endl;  
                
                
                const unsigned int targetcls = targetclassdecoder->gethighestclass() + 1;             
                targetclassdecoder->add(targetcls, word);
                
                //try to make do without an explicit encoder
                int size = 0;
                unsigned char * buffer = inttobytes(targetcls, size);                
                EncNGram targetunigram = EncNGram(buffer, size);
                delete buffer; 
                lm->ngrams[targetunigram] = lm->ngrams[UNKNOWNUNIGRAM];
                
                const EncAnyGram * sourcekey = translationtable->getsourcekey((const EncAnyGram *) unigram);                
                if (sourcekey == NULL) {
                    translationtable->alignmatrix[(const EncAnyGram *) unigram];
                    keepunigram = true;
                    sourcekey = translationtable->getsourcekey((const EncAnyGram *) unigram);
                }
                const EncAnyGram * targetkey = translationtable->gettargetkey((const EncAnyGram *) &targetunigram);
                if (targetkey == NULL) {
                    translationtable->targetngrams.insert(targetunigram);
                    targetkey = translationtable->gettargetkey((const EncAnyGram *) &targetunigram);
                }
                vector<double> scores;
                for (unsigned int j = 0; j < tweights.size(); j++) scores.push_back(1);                 
                
                t_aligntargets translationoptions;                
                translationoptions[targetkey] = scores;
                                                
                sourcefragments.push_back(SourceFragmentData( sourcekey, CorpusReference(0,i), translationoptions) );
                if (!keepunigram) delete unigram; 
                
               /* if (DEBUG >= 3) {
                    cerr << "\t" << i << ":1" << " -- " << sourcekey->decode(*sourceclassdecoder) << " ==> " << targetkey->decode(*targetclassdecoder) << " [ ";
                    for (unsigned int j = 0; j < tweights.size(); j++) cerr << translationtable->alignmatrix[sourcekey][targetkey][j] << " ";
                    cerr << "];" << endl;                     
                }*/
            }
        }
        
        
        this->tweights = vector<double>(tweights.begin(), tweights.end());
        this->dweight = dweight;
        this->lweight = lweight;
        this->dlimit = dlimit;
        
        if (DEBUG >= 3) cerr << "\tComputing future cost:" << endl;
        computefuturecost();
        
        

}


void StackDecoder::computefuturecost() {
        map<pair<int,int>, double> sourcefragments_costbyspan;
        //reorder source fragments by span for more efficiency
        for (t_sourcefragments::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {
            
            const EncAnyGram * candidate = iter->sourcefragment;
            const CorpusReference ref = iter->ref;
            const int n = candidate->n();    
            const pair<int,int> span = make_pair<int,int>((int) ref.token, (int) n);
             
            //cerr << "DEBUG: " << span.first << ':' << span.second << endl;             
            if (iter->translationoptions.size() == 0) {
                    cerr << "INTERNAL ERROR: No translation options" << endl;
                    throw InternalError();
            }
            
            //find cheapest translation option
            double bestscore = -INFINITY;            
            for (t_aligntargets::iterator iter2 = iter->translationoptions.begin(); iter2 != iter->translationoptions.end(); iter2++) {
                if (tweights.size() > iter2->second.size()) {
                    cerr << "ERROR: Too few translation scores specified for an entry in the translation table. Expected at least "  << tweights.size() << ", but got " << iter2->second.size() << " instead. Did you set -W correctly for the specified translation table?" << endl;
                    throw InternalError();
                }
                double score = 0; 
                for (unsigned int i = 0; i < tweights.size(); i++) {
                    double p = iter2->second[i];
                    if (p > 0) p = log(p); //turn into logprob, base e 
                    score += tweights[i] * p;
                }
                const EncAnyGram * translationoption = iter2->first;
                if (translationoption->isskipgram()) {
                    if (DEBUG) cerr << "debug mark" << endl;
                    vector<EncNGram*> parts;
                    (*((const EncSkipGram *) translationoption)).parts(parts);
                    for (vector<EncNGram*>::iterator iter3 = parts.begin(); iter3 != parts.end(); iter3++) {                        
                        EncNGram * part = *iter3;
                        score += lweight * lm->score(part);
                        delete part;
                    }  
                } else {
                    const EncNGram * ngram = (const EncNGram *) translationoption;
                    score += lweight * lm->score(ngram);
                }
                if (score > bestscore) {
                    bestscore = score;
                }                            
            } 
            sourcefragments_costbyspan[span] = bestscore; 
        }   
        //compute future cost
        for (unsigned int length = 1; length <= inputlength; length++) {
            for (unsigned int start = 0; start < inputlength - length + 1; start++) {
                const pair<int,int> span = make_pair((int) start,(int) length);
                map<pair<int,int>, double>::iterator iter = sourcefragments_costbyspan.find(span);
                if (iter != sourcefragments_costbyspan.end()) {
                    futurecost[span] = sourcefragments_costbyspan[span];
                } else {
                    /*if (length == 1){                       
                        cerr << "INTERNAL ERROR: No sourcefragment covers " << span.first << ":" << span.second << " ! Unable to compute future cost!" << endl;
                        exit(6);
                    } else {*/
                        futurecost[span] = -INFINITY;
                    /*}*/
                }
                for (unsigned int l = 1; l < length; l++) {
                    double spanscore = futurecost[make_pair((int) start,(int) l)] + futurecost[make_pair((int) start+l,(int) length - l)];
                    if (spanscore > futurecost[span]) { //(higher score -> lower cost)
                        if (DEBUG >= 3) {
                            cerr << "[" << span.first << ":" << span.second << "]" << " = [" << start << ":" << l << "] + [" << start + l << ":" << length -l << "] = " << futurecost[make_pair((int) start,(int) l)] << " + " << futurecost[make_pair((int) start+l,(int) length - l)] << " = " << spanscore << endl;
                        }
                        futurecost[span] = spanscore;
                    }
                }
            }
        }
        
        if (DEBUG >= 3) {
            cerr << "\tFuture cost precomputation:" << endl;
            for (map<std::pair<int, int>, double>::iterator iter = futurecost.begin(); iter != futurecost.end(); iter++) { 
                cerr << "\t  " << iter->first.first << ":" << iter->first.second << " = " << iter->second; // << "\tbase_e=" << log_e;
                if (sourcefragments_costbyspan.find(iter->first) != sourcefragments_costbyspan.end()) {
                     cerr << " *";
                }
                cerr << endl;                
            }            
        }
}

TranslationHypothesis * StackDecoder::decodestack(Stack & stack) {
        bool dead = false;
        unsigned int totalexpanded = 0;
        bool first = true;        
        TranslationHypothesis * fallbackhyp = NULL;
        
        while (!stack.empty()) {
            if (DEBUG >= 1) cerr << "\t Expanding hypothesis off stack " << stack.index << " -- " << stack.size() -1 << " left:" << endl;
            
            //pop from stack
            TranslationHypothesis * hyp = stack.pop();
            if (first) {                
                fallbackhyp = hyp;
                fallbackhyp->keep = true; //prevent fallback hypothesis from getting pruned
                first = false;
                if (DEBUG == 99) cerr << "DEBUG: FALLBACKHYP=" << (size_t) fallbackhyp << endl;
            }
            if (DEBUG >= 2) {
                cerr << "\t  Popped from stack:" << endl;
                hyp->report();
                cerr << "\t  Expands to:" << endl;
            }
                         
            
            unsigned int expanded = hyp->expand(); //will automatically add to appropriate stacks
            if (DEBUG >= 1) cerr << "\t  Expanded " << expanded << " new hypotheses" << endl;
            totalexpanded += expanded;            
            if ((hyp->deletable()) && (hyp != fallbackhyp)) {                
                if (DEBUG == 99) {
                    cerr << "DEBUG: DONE EXPANDING, DELETING HYPOTHESIS " << (size_t) hyp << endl;
                    hyp->cleanup();
                } else {
                    delete hyp; //if the hypothesis failed to expand it will be pruned
                }                
            }
        }
        if ((totalexpanded == 0) && ((unsigned int) stack.index != inputlength)) {
            dead = true;
            for (unsigned int j = stack.index + 1; j <= inputlength; j++) {
                if (!stacks[j].empty() || (!gappystacks[j].empty())) {
                    dead = false;
                    break;
                }                
            }
            if (dead) {
                cerr << "WARNING: DECODER ENDED PREMATURELY AFTER STACK " << stack.index << " of " << inputlength <<" , NO FURTHER EXPANSIONS POSSIBLE !" << endl;
                //return fallbackhyp;
            }
        }                
        if (fallbackhyp != NULL) fallbackhyp->keep = false;
        if ((!dead) && (fallbackhyp != NULL) && (fallbackhyp->deletable())) {
            if (DEBUG == 99) {
                cerr << "DEBUG: DELETING FALLBACK HYPOTHESIS " << (size_t) fallbackhyp << endl;
                fallbackhyp->cleanup();                
            } else {
                delete fallbackhyp;                
            }            
        }        
        if (DEBUG >= 1) cerr << "\t Expanded " << totalexpanded << " new hypotheses in total after processing stack " << stack.index << endl;
        unsigned int totalpruned = 0;
        for (int j = inputlength; j >= stack.index+1; j--) { //prune further stacks (hypotheses may have been added to any of them).. always prune superior stack first
            unsigned int pruned = stacks[j].prune();
            if ((DEBUG >= 1) && (pruned > 0)) cerr << "\t  Pruned " << pruned << " hypotheses from gapless stack " << j << endl;
            totalpruned += pruned; 
            
            pruned = gappystacks[j].prune();
            if ((DEBUG >= 1) && (pruned > 0)) cerr << "\t  Pruned " << pruned << " hypotheses from gappy stack " << j << endl;
            totalpruned += pruned; 
            
        }
        if (DEBUG >= 1) cerr << "\t Pruned " << totalpruned << " hypotheses from all superior stacks" << endl;        
        if (!dead) {
            fallbackhyp = NULL;
            stack.clear(); //stack is in itself no longer necessary, included pointer elements may live on though! Unnecessary hypotheses will be cleaned up automatically when higher-order hypotheses are deleted
        }
        stats.expanded += totalexpanded;
        stats.pruned += totalpruned;
        return fallbackhyp;
}
    
TranslationHypothesis * StackDecoder::decode() {
   /*
    initialize hypothesisStack[0 .. nf];
    create initial hypothesis hyp_init;
    add to stack hypothesisStack[0];
    for i=0 to nf-1:
    for each hyp in hypothesisStack[i]:
     for each new_hyp that can be derived from hyp:
       nf[new_hyp] = number of source words covered by new_hyp;
       add new_hyp to hypothesisStack[nf[new_hyp]];
       prune hypothesisStack[nf[new_hyp]];
    find best hypothesis best_hyp in hypothesisStack[nf];
    output best path that leads to best_hyp;
  */

    
    
    //create initial hypothesis and add to 0th stack
    TranslationHypothesis * initialhypothesis = new TranslationHypothesis(NULL, this, NULL, 0,NULL,0, vector<double>() );   
    stacks[0].add(initialhypothesis);
    
    TranslationHypothesis * fallbackhyp = NULL;
    
    //for each stack
    for (unsigned int i = 0; i <= inputlength - 1; i++) {
        if (DEBUG >= 1) cerr << "\tDecoding gapless stack " << i << " -- " << stacks[i].size() << " hypotheses" << endl;
        stats.stacksizes[i] = stacks[i].size();
        fallbackhyp = decodestack(stacks[i]);
        if (fallbackhyp != NULL) break;
        
        if (DEBUG >= 1) cerr << "\tDecoding gappy stack " << i << " -- " << gappystacks[i].size() << " hypotheses" << endl;
        stats.gappystacksizes[i] = gappystacks[i].size();
        fallbackhyp = decodestack(gappystacks[i]);
        if (fallbackhyp != NULL) break;
    }   
    
    //solutions are now in last stack: stacks[inputlength]       
     
    
    if (!stacks[inputlength].empty()) {
        TranslationHypothesis * solution = stacks[inputlength].pop();
        return solution;
    } else if (fallbackhyp != NULL) {
        return fallbackhyp;
    } else {
        return NULL;
    }
}




StackDecoder::~StackDecoder() {
    for (int i = inputlength; i >= 0; i--) {
        stacks[i].clear();
    }
}

/*
string StackDecoder::solution(ClassDecoder & targetclassdecoder) {
    for (int stackindex = inputlength; stackindex >= 1; stackindex--) {
        if (!stacks[stackindex].empty()) {
            if (stackindex < inputlength) {
                cerr << "WARNING: INCOMPLETE SOLUTION (LOCAL MAXIMUM IN SEARCH), FROM STACK " << stackindex << " instead of " << inputlength << endl;
            }
            TranslationHypothesis * sol = stacks[stackindex].pop();
            EncData s = sol->getoutput();
            cerr << "SCORE=" << sol->score() << endl;
            return s.decode(targetclassdecoder);
        }
    }
    cerr << "NO SOLUTION FOUND!";
    exit(24);
}
*/

Stack::Stack(StackDecoder * decoder, int index, int stacksize, double prunethreshold) {
    this->decoder = decoder;
    this->index = index;
    this->stacksize = stacksize;
    this->prunethreshold = prunethreshold;
}

Stack::~Stack() {
    //note: superior stacks always have to be deleted before lower stacks!
    for (list<TranslationHypothesis*>::iterator iter = contents.begin(); iter != contents.end(); iter++) {
        TranslationHypothesis * h = *iter;
        if (h->deletable()) {
            if (decoder->DEBUG == 99) {
                cerr << "DEBUG: STACK DESTRUCTION. DELETING " << (size_t) h << endl;
                h->cleanup();
            } else {
                delete h;
            }
        }
    }
}

Stack::Stack(const Stack& ref) { //limited copy constructor
    decoder = ref.decoder;
    index = ref.index;
    stacksize = ref.stacksize;
    prunethreshold = ref.prunethreshold;
}

void Stack::clear() {
    if (decoder->DEBUG >= 3) cerr << "\tClearing stack " << index << endl; 
    for (list<TranslationHypothesis*>::iterator iter = contents.begin(); iter != contents.end(); iter++) {
        TranslationHypothesis * h = *iter;
        if (h->deletable()) {
            if (decoder->DEBUG == 99) {
                cerr << "DEBUG: CLEARING STACK " << index << ", DELETING HYPOTHESIS" << (size_t) this << endl;
                h->cleanup();
            } else {
                delete h;
            }
        }
    }
    contents.clear();
}

TranslationHypothesis * Stack::pop() {
    if (contents.empty()) return NULL;
    TranslationHypothesis * h = contents.front();
    contents.pop_front();
    return h;
}

double Stack::bestscore() {
    if (contents.empty()) return -999999999;
    TranslationHypothesis * h = contents.front();
    return h->score();
}

double Stack::worstscore() {
    if (contents.empty()) return -999999999;
    TranslationHypothesis * h = contents.back();
    return h->score();
}

bool Stack::add(TranslationHypothesis * candidate) {
    //IMPORTANT NOTE: the task of deleting the hypothesis when discarded is left to the caller!

    if (contents.empty()) {
        //empty, add:
        contents.push_back(candidate);
        return true;
    }
    
    double score = candidate->score();
    if (contents.size() >= (size_t) stacksize) {
        if (score < worstscore()) {
            return false;
        } 
    }
    
    //insert at right position and do histogram pruning
    bool added = false;
    int count = 0;
    for (list<TranslationHypothesis*>::iterator iter = contents.begin(); iter != contents.end(); iter++) {
        count++;
        TranslationHypothesis * h = *iter;
        if ( (!added) && (count < stacksize) && (score >= h->score()) ) {
            //insert here
            count++;
            contents.insert(iter, candidate);
            added = true;
        }
        if (count > stacksize) {
            if (h->deletable()) {
                if (decoder->DEBUG == 99) {
                    cerr << endl << "DEBUG: STACK OVERFLOW, PRUNING BY DELETING " << (size_t) h << endl;
                    h->cleanup();
                } else {                
                    delete h;
                }
            }
            contents.erase(iter);
            break;
        }
    }
    if ((!added) && (contents.size() < (size_t) stacksize)) {
        added = true;
        contents.push_back(candidate);
    }
    return added;    
}


int Stack::prune() {
    int pruned = 0;
    if ((prunethreshold != 1) && (prunethreshold != 0)) {
        //pruning based on prunethreshold
        double cutoff = bestscore() + log(prunethreshold);
        for (list<TranslationHypothesis*>::iterator iter = contents.begin(); iter != contents.end(); iter++) {
            TranslationHypothesis*  h = *iter;
            if (h->score() < cutoff) {
                pruned++;
                iter = contents.erase(iter);
                if (h->deletable()) {
                    if (decoder->DEBUG == 9) {
                        cerr << "DEBUG: SCORE EXCEEDS CUTOFF, PRUNING BY DELETING " << (size_t) h << endl;                    
                        h->cleanup();
                    } else {
                        delete h; //TODO: REENABLE AND FIX -- MEMORY LEAK?
                    }
                }
            }
        }
    }
    return pruned;
}


TranslationHypothesis::TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const vector<double> & tscores) {
    this->keep = false;
    this->parent = parent;
    this->decoder = decoder;            
    if (parent != NULL) parent->children.push_back(this);        
    this->sourcegram = sourcegram;
    this->targetgram = targetgram;
    this->sourceoffset = sourceoffset;
    this->targetoffset = targetoffset; 
    this->deleted = false; //debug only
    
    if ((sourcegram != NULL) && (sourcegram->isskipgram())) {
        vector<pair<int, int> > gaps;
        (*((const EncSkipGram *) sourcegram)).getgaps(gaps);
        for (vector<pair<int, int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) sourcegaps.push_back( make_pair<unsigned char, unsigned char>(iter->first, iter->second) );        
    }
    

    vector<bool> targetcoverage; //intermediary data structure used for computing target gaps
    const TranslationHypothesis * h = this;
    while (h->parent != NULL) {
        while (targetcoverage.size() < (size_t) (h->targetoffset + h->targetgram->n()) ) {
           targetcoverage.push_back(false);   
        }        
        if (!h->targetgram->isskipgram()) {
            for (int i = h->targetoffset; i < h->targetoffset + h->targetgram->n(); i++) {
                 targetcoverage[i] = true;
            } 
        } else {
            vector<pair<int, int> > parts;
            (*((const EncSkipGram *) h->targetgram)).getparts(parts);
            /*if (decoder->DEBUG >= 4) {
                    (*((const EncSkipGram *) h->targetgram)).out();
                    cerr << endl;
            }*/
            for (vector<pair<int, int> >::iterator iter = parts.begin(); iter != parts.end(); iter++) {
                if (decoder->DEBUG >= 4) cerr << "DEBUG: part: " << iter->first << ":" << iter->second << endl;
                for (int i = h->targetoffset + iter->first; i < h->targetoffset + iter->first + iter->second; i++) {
                    targetcoverage[i] = true;
                } 
            }               
        }        
        h = h->parent;
    }
        
    int gapbegin = 0;
    int gaplength = 0;
    if (decoder->DEBUG >= 4) cerr << "DEBUG: Targetcoverage: ";
    for (unsigned int i = 0; i <= targetcoverage.size(); i++) {
        if ((decoder->DEBUG >= 4) && (i < targetcoverage.size())) {
            cerr << (int) targetcoverage[i];                
        }
        if (i == targetcoverage.size() || targetcoverage[i]) {
            if (gaplength > 0) {
                targetgaps.push_back( make_pair<unsigned char, unsigned char>( (unsigned char) gapbegin, (unsigned char) gaplength) );        
            }
            gapbegin = i + 1;
            gaplength= 0;
        } else {
            gaplength++;
        }
    }        
    if (decoder->DEBUG >= 4) cerr << endl;
 
    
           
    //compute input coverage
    if (parent == NULL) {
        //initial hypothesis: no coverage at all
        for (unsigned int i = 0; i < decoder->inputlength; i++) inputcoveragemask.push_back(false);
    } else {
        //inherit from parent
        for (unsigned int i = 0; i < decoder->inputlength; i++) {
             inputcoveragemask.push_back(parent->inputcoveragemask[i]);
        }
        //add sourcegram coverage
        bool isskipgram = sourcegram->isskipgram();
        for (int i = sourceoffset; i < sourceoffset + sourcegram->n(); i++) {
            if (isskipgram) {
                //cerr << "DEBUG: processing " << i << " length=" << (int) sourcegram->n() << ": ";
                bool ingap = false;
                for (vector<pair<unsigned char, unsigned char> >::iterator iter = sourcegaps.begin(); iter != sourcegaps.end(); iter++) {
                    if ( (i >= sourceoffset + iter->first) &&  ( i < sourceoffset + iter->first + iter->second ) ) {
                        ingap = true;
                        break;       
                    }             
                }
                if (!ingap) {
                    //cerr << "no gap" << endl;
                    inputcoveragemask[i] = true;
                } else {
                    //cerr << "in a gap" << endl;
                }
            } else {
                inputcoveragemask[i] = true;
            }        
        }    
    }
    
    //find history (last order-1 words) for computation for LM score    
    history = NULL;
    
    const int order = decoder->lm->getorder();
    int begin = targetoffset - (order - 1);
    if (begin < 0) {        
        history = new EncNGram((const unsigned char*) &BOSCLASS, 1);
        begin = 0;
    }
     
    for (int i = begin; i < targetoffset; i++) {
        EncNGram * unigram; 
        //cerr << "DEBUG: history->getoutputtoken " << i << endl;
        unigram = getoutputtoken(i);
            
        if (!unigram->unknown()) {
            if (history != NULL) {
                const EncNGram * old = history;
                
                history = new EncNGram(*history + *unigram);
                delete old;
            } else {
                history = new EncNGram(*unigram);
            }
        } else if (history != NULL) {
            //we have an unknown unigram, erase history and start over from this point on
            delete history;
            history = NULL;
        }
        delete unigram;
    }
    
    
    //Precompute score
    if ((parent != NULL) && (decoder->tweights.size() > tscores.size())) {
        cerr << "Too few translation scores specified for an entry in the translation table. Expected at least "  << decoder->tweights.size() << ", but got " << tscores.size() << " instead. Did you set -W correctly for the specified translation table?" << endl;
        throw InternalError();
    }
  
    this->tscores = tscores;
  
    tscore = 0; 
    for (unsigned int i = 0; i < tscores.size(); i++) {
        double p = tscores[i];
        if (p > 0) p = log(p); //turn into logprob, base e 
        tscore += decoder->tweights[i] * p;
    }
    
    
    lmscore = 0;
    
    if (parent != NULL) {
        if (targetgram->isskipgram()) {
            vector<EncNGram*> parts;
            ((const EncSkipGram *) targetgram)->parts(parts);
            int partcount = 0;
            for (vector<EncNGram*>::iterator iter = parts.begin(); iter != parts.end(); iter++) {
                EncNGram * part = *iter;
                if ((partcount == 0) && (sourcegaps[0].first != 0)) {
                    //first part, no initial gap, call with history
                    //EncNGram ngram = *history + *part;
                    lmscore += decoder->lweight * decoder->lm->score(part, history);
                } else {     
                    lmscore += decoder->lweight * decoder->lm->score(part);
                }
                partcount++;
                delete part;
            } 
        } else {                
            //EncNGram ngram = EncNGram(*history + *((const EncNGram* ) targetgram) );
            //if (decoder->DEBUG >= 4) cerr << "DEBUG4: Calling LM with history for: '" << ngram.decode( *decoder->targetclassdecoder) << "'";
            lmscore += decoder->lweight * decoder->lm->score( (const EncNGram *) targetgram, history);
            //if (decoder->DEBUG >= 4) cerr << " ... DONE" << endl; 
        }
    }
    
       
    if (final()) {
       EncNGram * terminator = NULL;
       
       //find length
       int targetlength = 0;
       TranslationHypothesis * h = this;
       while ((h != NULL) && (h->targetgram != NULL)) {
            int n = h->targetgram->n();
            if (h->targetoffset + n > targetlength) targetlength = h->targetoffset + n;
            h = h->parent; 
       } 
         
       for (int i = targetlength - (order - 1); (i >= 0) && (i < targetlength); i++) {
            if (terminator == NULL) {
                EncNGram * unigram; 
                //cerr << "DEBUG: final1->getoutputtoken " << i << " (" << targetlength << " )" << endl;
                unigram = getoutputtoken(i);
                terminator = new EncNGram(*unigram);
                delete unigram;               
            } else {
                EncNGram * unigram; 
                //cerr << "DEBUG: final2->getoutputtoken " << i << endl;
                unigram = getoutputtoken(i);
                terminator = new EncNGram(*terminator + *unigram);
                delete unigram;
            }
       } 
       EncNGram eos = EncNGram((const unsigned char*) &EOSCLASS, 1);
       //if (decoder->DEBUG >= 4) cerr << "DEBUG4: Calling LM for terminator: " << terminator->decode( *decoder->targetclassdecoder);
       lmscore += decoder->lweight * decoder->lm->scoreword(&eos , terminator );
       //if (decoder->DEBUG >= 4) cerr << "...DONE" << endl;
       delete terminator;
    }
    
    
    //TODO: deal differently with filling in gaps in skips?
    //Total reordering cost is computed by D(e,f) = - Σi (d_i) where d for each phrase i is defined as d = abs( last word position of previously translated phrase + 1 - first word position of newly translated phrase ).
    int prevpos = 0;
    if ((parent != NULL) && (!parent->initial())) {
        prevpos = parent->sourceoffset + parent->sourcegram->n();        
    } 
    double distance = abs( prevpos - sourceoffset);
    
    dscore = decoder->dweight * -distance; 
        
    //Compute estimate for future score
    
    futurecost = 0;
    begin = -1;    
    for (unsigned int i = 0; i <= inputcoveragemask.size() ; i++) {
        if ((!inputcoveragemask[i]) && (begin == -1) && (i < inputcoveragemask.size())) {
            begin = i;
        } else if (((i == inputcoveragemask.size()) || (inputcoveragemask[i])) && (begin != -1)) {
            const double c = decoder->futurecost[make_pair((int) begin,(int) i-begin)];
            if (c == 0) {
                cerr << "INTERNAL ERROR: Future cost for " << begin << ":" << i - begin << " is 0! Not possible!" << endl;
                report();
                throw InternalError();
            }           
            //cerr << "DEBUG: Adding futurecost for " << begin << ":" << i - begin << " = " <<c << endl;
            futurecost += c;
            begin = -1;
        }
    }        
    if (decoder->DEBUG >= 2) {
        report();
    }
  
}

void TranslationHypothesis::report() {
        const double _score = tscore + lmscore + dscore;
        if (decoder->DEBUG == 99) {
            cerr << "\t   Translation Hypothesis " << (size_t) this << endl;
        } else {
            cerr << "\t   Translation Hypothesis "  << endl;  //<< (size_t) this << endl;
        }
        if (decoder->sourceclassdecoder != NULL) { 
            cerr << "\t    Source Fragment: ";
            if (sourcegram == NULL) {
                cerr << "NONE (Initial Hypothesis!)" << endl;
            } else if (sourcegram->isskipgram()) {
                cerr << ((const EncSkipGram *) sourcegram)->decode(*(decoder->sourceclassdecoder)) << endl;
            } else {
                cerr << ((const EncNGram *) sourcegram)->decode(*(decoder->sourceclassdecoder)) << endl;
            }
        }
        if (decoder->targetclassdecoder != NULL) { 
            cerr << "\t    Target Fragment: ";
            if (targetgram == NULL) {
                cerr << "NONE (Initial Hypothesis!)" << endl;
            } else if (targetgram->isskipgram()) {
                cerr << ((const EncSkipGram *) targetgram)->decode(*(decoder->targetclassdecoder)) << endl;
            } else {
                cerr << ((const EncNGram *) targetgram)->decode(*(decoder->targetclassdecoder)) << endl;
            }
            if ((targetgram != NULL)) {
                EncData s = getoutput();
                cerr << "\t    Translation: " << s.decode(*(decoder->targetclassdecoder)) << endl;
            }
        }
        cerr << "\t    targetoffset: " << (int) targetoffset << endl;
        cerr << "\t    History: ";
        if (history == NULL) {
            cerr << "NONE" << endl;
        } else {
            cerr << history->decode(*(decoder->targetclassdecoder)) << endl;
        }
        cerr << "\t    tscore = ";
        for (unsigned int i = 0; i < tscores.size(); i++) {
            if (i > 0) cerr << " + ";
            cerr << decoder->tweights[i] << " * " << tscores[i];
        }
        cerr << " = " << tscore << endl;
        cerr << "\t    fragmentscore = tscore + lmscore + dscore = " << tscore << " + " << lmscore << " + " << dscore << " = " << _score << endl;
        cerr << "\t    totalscore = basescore + fragmentscore + futurecost = " << basescore() << " + " << _score << " + " << futurecost << " = " << score() << endl;
        cerr << "\t    coverage: ";
        for (unsigned int i = 0; i < inputcoveragemask.size(); i++) {
            if (inputcoveragemask[i]) {
                cerr << "1";
            } else {
                cerr << "0";
            }            
        }
        cerr << endl;
        if (!sourcegaps.empty()) {
            cerr << "\t    sourcegaps: ";
            for (vector<pair<unsigned char,unsigned char>>::iterator iter = sourcegaps.begin(); iter != sourcegaps.end(); iter++) {
                cerr << (int) iter->first << ':' << (int) iter->second << ' ';
            }
            cerr << endl;
        }               
        if (!targetgaps.empty()) {
            cerr << "\t    targetgaps: ";
            for (vector<pair<unsigned char,unsigned char>>::iterator iter = targetgaps.begin(); iter != targetgaps.end(); iter++) {
                cerr << (int) iter->first << ':' << (int) iter->second << ' ';
            }
            cerr << endl;
        }        
}

double TranslationHypothesis::basescore() const {
   double s = 0;
    const TranslationHypothesis * h = this->parent;
    while (h != NULL) {
        s += h->tscore + h->lmscore + h->dscore;
        h = h->parent;
    } 
    return s;    
}

double TranslationHypothesis::score() const {
    return basescore() + tscore + lmscore + dscore + futurecost;
} 


bool TranslationHypothesis::deletable() {
    return ((!keep) && (children.empty()));
}

void TranslationHypothesis::cleanup() {  
    if (decoder->DEBUG == 99) {
        cerr << "DEBUG: DELETING HYPOTHESIS " << (size_t) this << endl;
        if (deleted) {
             cerr << "INTERNAL ERROR: DELETING AN ALREADY DELETED HYPOTHESIS!!! THIS SHOULD NOT HAPPEN!!!!!" << (size_t) this << endl;
             throw InternalError();
        }
        deleted = true;
    }
     
    if (history != NULL) {
        delete history;
    }
    if (parent != NULL) {                 
        vector<TranslationHypothesis *>::iterator iter = find(parent->children.begin(), parent->children.end(), this); 
        parent->children.erase(iter);
        if (parent->deletable()) {
            if (decoder->DEBUG == 99) {
                cerr << "DEBUG: DELETING PARENT HYPOTHESIS " << (size_t) parent << endl;
                parent->cleanup();
            } else {
                delete parent;
            }
        }
    } 
}


TranslationHypothesis::~TranslationHypothesis() {  
    cleanup();
}

bool TranslationHypothesis::hasgaps() const {
    return (!targetgaps.empty());
}

bool TranslationHypothesis::final(){
        if (!targetgaps.empty()) return false;
        return ((unsigned int) inputcoverage() == decoder->inputlength); 
}

int TranslationHypothesis::fitsgap(const EncAnyGram * candidate, const int offset) {
    //returns -1 if no fit, begin index of gap otherwise
    int i = 0;
    for (vector<pair<unsigned char, unsigned char> >::iterator iter = targetgaps.begin(); iter != targetgaps.end(); iter++) {
        if (i >= offset) { //offset is gap index (for skipping gaps), NOT begin offset        
            if ((candidate->variablewidth()) || (candidate->n() <= iter->second)) {
                return iter->first;
            }
        }
        i++; 
    }     
    return -1;
}

unsigned int TranslationHypothesis::expand() {
    bool oldkeep = this->keep;
    this->keep = true; //lock this hypothesis, preventing it from being deleted when expanding and rejecting its last dying child
    if (deleted) { 
        //only in debug==99
        cerr << "ERROR: Expanding an already deleted hypothesis! This should never happen! " << (size_t) this << endl;
        throw InternalError();
    }
    
    
    unsigned int expanded = 0;
    int thiscov = inputcoverage();
    //expand directly in decoder.stack()

    //find new source fragment
    for (t_sourcefragments::iterator iter = decoder->sourcefragments.begin(); iter != decoder->sourcefragments.end(); iter++) {
        const EncAnyGram * sourcecandidate = iter->sourcefragment;
        const CorpusReference ref = iter->ref; 
        if (!conflicts(sourcecandidate, ref)) {
            //find target fragments for this source fragment
            if (iter->translationoptions.empty()) {
                cerr << "ERROR: Translation options are empty! This should never happen!" << endl;
                throw InternalError();
            }
            int c = 0;                
            for (t_aligntargets::const_iterator iter2 =  iter->translationoptions.begin(); iter2 != iter->translationoptions.end(); iter2++) {
                c++;

                if ((decoder->dlimit >= 0) && (decoder->dlimit < 999)) {                             
                    int prevpos = 0;
                    if (sourcegram != NULL) {                                                                                    
                        prevpos = sourceoffset + sourcegram->n();        
                    } 
                    double distance = abs( prevpos - sourceoffset);
                    if (distance > decoder->dlimit) continue; //skip
                }
                
                //cerr << "DEBUG: " << c << " of " << decoder->translationtable->alignmatrix[sourcecandidate].size() << endl;  
                //create hypothesis for each target fragment
                const EncAnyGram * targetcandidate = iter2->first;
                int length;
                
                if (targetgram != NULL) { 
                    length = targetgram->n();
                } else {    
                    length = 0;
                }
                int newtargetoffset = targetoffset + length;
                
                
                //If there are target-side gaps, this expansion must fill (one of) those first
                int gapoffset = 0;
                int fitsgapindex;
                do {
                    fitsgapindex = fitsgap(targetcandidate, gapoffset++);
                    if (hasgaps()) {
                        if (fitsgapindex == -1) {
                            //can't fill a gap, don't create hypothesis, break
                            break;
                        } else{
                            newtargetoffset = fitsgapindex;
                        }             
                    }                    
                    
                    TranslationHypothesis * newhypo = new TranslationHypothesis(this, decoder, sourcecandidate, ref.token, targetcandidate, newtargetoffset , iter2->second);                    
                    if ((!newhypo->fertile()) && (!newhypo->final())) {
                        if (decoder->DEBUG >= 3) cerr << "\t    Hypothesis not fertile, discarding..."<< endl;
                        if (decoder->DEBUG == 99) {
                            newhypo->cleanup();
                        } else {
                            delete newhypo;
                        }                         
                        continue;
                    }                    
                    //add to proper stack
                    int cov = newhypo->inputcoverage();
                    if (thiscov >= cov) {
                        cerr << "INTERNAL ERROR: Hypothesis expansion did not lead to coverage expansion! This should not happen. New hypo has coverage " << cov << ", parent: " << thiscov << endl;
                        continue;        
                    }
                    
                    bool accepted;
                    if (newhypo->hasgaps()) {
                        if (decoder->DEBUG >= 2) cerr << "\t    Adding to gappy stack " << cov;
                        accepted = decoder->gappystacks[cov].add(newhypo);
                    } else {
                        if (decoder->DEBUG >= 2) cerr << "\t    Adding to gapless stack " << cov;
                        accepted = decoder->stacks[cov].add(newhypo);
                        if (hasgaps()) decoder->stats.gapresolutions++;
                    }
                    if (decoder->globalstats && accepted) {
                        newhypo->stats();
                    }                    
                    if (decoder->DEBUG >= 2) {
                        if (accepted) {
                            cerr << " ... ACCEPTED" << endl;
                        } else {
                            cerr << " ... REJECTED" << endl;
                        }
                    }
                    expanded++;                    
                    
                    if (!accepted) {
                        decoder->stats.discarded++;
                        if (decoder->DEBUG == 99) {
                            cerr << "DEBUG: IMMEDIATELY DELETING NEWLY CREATED HYPOTHESIS (REJECTED BY STACK) " << (size_t) newhypo << endl;
                            newhypo->cleanup();
                        } else {
                            delete newhypo;
                        } 
                    }
        
                } while (fitsgapindex != -1); //loop until no gaps can be filled (implies only a single iteration if there are no gaps at all)  
            
            }                 
        }        
    }

    this->keep = oldkeep; //release lock
    return expanded;
}

bool TranslationHypothesis::conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref, bool skipduplicatecheck) {
    if ((sourcegram == NULL) && (parent == NULL)) return false; //no sourcegram, no conflict (this is an empty initial hypothesis)
    
    if (skipduplicatecheck) {
        const TranslationHypothesis * h = this;
        size_t candidatehash = sourcecandidate->hash(); 
        while  (h->parent != NULL) {
            if (candidatehash == h->sourcegram->hash()) return true; //source was already added, can not add twice
            h = h->parent;
        }
    }        
    
    if ( (ref.token + sourcecandidate->n() > sourceoffset) && (ref.token < sourceoffset + sourcegram->n()  ) ) { 
        //conflict    
        if (sourcegram->isskipgram()) {
            //if this falls nicely into a gap than it may not be a conflict after all
            bool ingap = false;
            for (vector<pair<unsigned char, unsigned char> >::iterator iter = sourcegaps.begin(); iter < sourcegaps.end(); iter++) {
                if ( (ref.token + sourcecandidate->n() > sourceoffset + iter->first) &&  ( ref.token < sourceoffset +iter->first + iter->second ) ) {
                    ingap = true;
                    break;       
                }
            }            
            if (!ingap) {
                return true;
            }
        } else {
            //MAYBE TODO: deal with partial overlap?
            return true;
        }
    }
    
    //no confict, check parents    
    if (parent != NULL) {
        return parent->conflicts(sourcecandidate, ref, true);
    } else {
        return false;
    }
}

bool TranslationHypothesis::fertile() {
    vector<int> fertilitymask;
    for (unsigned int i = 0; i < decoder->inputlength; i++) {
        if (inputcoveragemask[i]) {
            fertilitymask.push_back(-1); //for tokens already covered
        } else {
            fertilitymask.push_back(0);
        }
    }
    for (t_sourcefragments::iterator iter = decoder->sourcefragments.begin(); iter != decoder->sourcefragments.end(); iter++) {
        const EncAnyGram * sourcecandidate = iter->sourcefragment;
        const CorpusReference ref = iter->ref;
        
        size_t candidatehash = sourcecandidate->hash();
        bool alreadyused = false;
        const TranslationHypothesis * h = this;
        while (h != NULL) {
            if ((h->sourcegram != NULL) && (candidatehash == h->sourcegram->hash())) {
                alreadyused = true;
                break;
            } 
            h = h->parent;            
        } 
        if (alreadyused) continue;
                
        int length = sourcecandidate->n();
    
        bool applicable = true;
        for (int i = ref.token; i < ref.token + length; i++) {
            if (fertilitymask[i] < 0) {
                applicable = false;
                break;
            } 
        }
        if (applicable) {
            for (int i = ref.token; i < ref.token + length; i++) {
                fertilitymask[i]++;
            }
        }    
    }
    if (decoder->DEBUG >= 3) {
        cerr << "\t     fertility: ";
        for (unsigned int i = 0; i < decoder->inputlength; i++) {
            if (fertilitymask[i] >= 0)
                cerr << fertilitymask[i];
            else 
                cerr << 'X';
        }
        cerr << endl;
    }
    //check fertility mask for 0s
    for (unsigned int i = 0; i < decoder->inputlength; i++) {
        if (fertilitymask[i] == 0) return false;
    }     
    return true;
}

int TranslationHypothesis::inputcoverage() {
    int c = 0;
    for (unsigned int i = 0; i < inputcoveragemask.size(); i++) {
        if (inputcoveragemask[i]) c++; 
    }
    return c;
}


EncNGram * TranslationHypothesis::getoutputtoken(int index) {
    if (parent == NULL) {
        cerr << "ERROR: getoutputtoken left unresolved!" << endl;
        return NULL;
    } else if ((index >= targetoffset) && (index < targetoffset + targetgram->n())) {
        index = index - targetoffset;
        if ((index < 0) || (index >= targetgram->n())) {
            cerr << "ERROR: TranslationHypothesis::getoutputtoken() with index " << index << " is out of bounds" << endl;
            throw InternalError();
        }
        if (targetgram->isskipgram()) {
            const EncSkipGram * targetskipgram = (const EncSkipGram *) targetgram;
            return targetskipgram->gettoken(index);
        } else {
            const EncNGram * targetngram = (const EncNGram *) targetgram;        
            return targetngram->slice(index,1);
        }
    } else {
        return parent->getoutputtoken(index);
    }
}

EncData TranslationHypothesis::getoutput(deque<TranslationHypothesis*> * path) { //get output        
    //backtrack
    if (path == NULL) path = new deque<TranslationHypothesis*>;
    if (parent != NULL) {
            path->push_front(this);
            return parent->getoutput(path);
    } else {
        //we're back at the initial hypothesis, now we construct the output forward again
        map<int, EncNGram *> outputtokens; //unigrams, one per index                                  
        while (!path->empty()) {
            TranslationHypothesis * hyp = path->front();
            path->pop_front();
            for (int i = hyp->targetoffset; i < hyp->targetoffset + hyp->targetgram->n(); i++) {
                outputtokens[i] = hyp->getoutputtoken(i);
            }                          
        }
        unsigned char buffer[8192];
        int cursor = 0;
        for (map<int, EncNGram *>::iterator iter = outputtokens.begin(); iter != outputtokens.end(); iter++) {
            EncNGram * unigram = iter->second;
            for (int i = 0; i < unigram->size(); i++) {
                buffer[cursor++] = unigram->data[i];            
            }
            buffer[cursor++] = 0;
            delete unigram; //important cleanup            
        }     
        delete path; //important cleanup
        return EncData(buffer, cursor); 
    }    
}

void TranslationHypothesis::stats() {
    TranslationHypothesis * h = this;
    int stepcount = 0;    
    while ((h != NULL) && (h->parent != NULL)) {
        stepcount++;
        const int source_n = h->sourcegram->n();
        if (h->sourcegram->isskipgram()) {            
            decoder->stats.sourceskipgramusage[source_n]++;
        } else {
            decoder->stats.sourcengramusage[source_n]++;
        }
        const int target_n = h->targetgram->n();
        if (h->targetgram->isskipgram()) {            
            decoder->stats.targetskipgramusage[target_n]++;
        } else {
            decoder->stats.targetngramusage[target_n]++;
        }
        h = h->parent;
    }
    decoder->stats.steps.push_back(stepcount);        
}
    
 

void usage() {
    cerr << "Colibri MT Decoder" << endl;
    cerr << "   Maarten van Gompel, Radboud University Nijmegen" << endl;
    cerr << "----------------------------------------------------" << endl;
    cerr << "Usage: decoder -t translationtable -l lm-file [-S source-class-file -T target-class-file]" << endl;
    cerr << " Input:" << endl;
    cerr << "   (input will be read from standard input)" << endl;
    cerr << "\t-t translationtable       Translation table (*.translationtable.colibri)" << endl;
    cerr << "\t-l lm-file                Language model file (in ARPA format)" << endl;    
    cerr << "\t-S sourceclassfile        Source class file" << endl;
    cerr << "\t-T targetclassfile        Target class file" << endl;
    cerr << "\t-s stacksize              Stacksize" << endl;
    cerr << "\t-p prune threshold        Prune threshold" << endl;
    cerr << "\t-W weight                 Translation model weight (issue multiple times for multiple weights: -W 1 -W 1" << endl;    
    cerr << "\t-L weight                 Language model weight" << endl;
    cerr << "\t-D weight                 Distortion model weight" << endl;
    cerr << "\t-M distortion-limit       Distortion limit (max number of source words skipped, default: unlimited)" << endl;
    cerr << "\t-N                        No skipgrams" << endl;
    cerr << "\t-C id                     Use classifier" << endl;            
    cerr << "\t-O options                Timbl options for classification" << endl;
    cerr << "\t--moses                   Translation table is in Moses format" << endl;
    cerr << "\t-v verbosity              Verbosity/debug level" << endl;
    cerr << "\t--stats                   Compute and output decoding statistics for each solution" << endl;
    cerr << "\t--globalstats             Compute and output decoding statistics for all hypothesis accepted on a stack" << endl;
    cerr << "\t-x mode                   How to handle classifier scores? (Only with -C). Choose from:" << endl;
    cerr << "\t       weighed            Apply classifier score as weight to original scores (default)" << endl;
    cerr << "\t       append             Append classifier score to translation score vector (make sure to specify an extra weight using -W)" << endl;
    cerr << "\t       replace            Use only classifier score, replacing translation table scores (make to specify only one weight using -W)" << endl;        
}

void addsentencemarkers(ClassDecoder & targetclassdecoder, ClassEncoder & targetclassencoder) {
    targetclassdecoder.add(BOSCLASS,"<s>");
    targetclassdecoder.add(EOSCLASS,"</s>");
    targetclassencoder.add("<s>", BOSCLASS);
    targetclassencoder.add("</s>", EOSCLASS);
}

int addunknownwords( AlignmentModel & ttable, LanguageModel & lm, ClassEncoder & sourceclassencoder, ClassDecoder & sourceclassdecoder,  ClassEncoder & targetclassencoder, ClassDecoder & targetclassdecoder, int tweights_size) {
    int added = 0;
    if (sourceclassencoder.gethighestclass() > sourceclassdecoder.gethighestclass()) {
        for (unsigned int i = sourceclassdecoder.gethighestclass() + 1; i <= sourceclassencoder.gethighestclass(); i++) {
            added++;
            unsigned int cls = i;
            const string word = sourceclassencoder.added[cls];
            sourceclassdecoder.add(cls, word);            
            
            unsigned int targetcls = targetclassencoder.gethighestclass() + 1;
            targetclassencoder.add(word, targetcls);
            targetclassdecoder.add(targetcls, word);
            
            cerr << "NOTICE: Unknown word in input: " << word << " (" << cls << ", " << targetcls << ")" << endl;
            
            EncNGram sourcegram = sourceclassencoder.input2ngram( word,false,false);
            EncNGram targetgram = targetclassencoder.input2ngram( word,false,false);
            
            
            //ttable.sourcengrams.insert(sourcegram);
            ttable.targetngrams.insert(targetgram);            
                
            lm.ngrams[targetgram] = lm.ngrams[UNKNOWNUNIGRAM];
            
            
            
            const EncAnyGram * sourcekey = ttable.getsourcekey((const EncAnyGram*) &sourcegram );
            if (sourcekey == NULL) {
                sourcekey = new EncNGram(sourcegram); //new pointer
            }
            const EncAnyGram * targetkey = ttable.gettargetkey((const EncAnyGram*) &targetgram );
            
            vector<double> scores;
            for (int j = 0; j < tweights_size; j++) scores.push_back(1);
            ttable.alignmatrix[sourcekey][targetkey] = scores;
            
            
        }
    }
    sourceclassencoder.added.clear();
    return added;    
}

void DecodeStats::output() {
    unsigned int totalsourcengrams = 0;
    unsigned int totaltargetngrams = 0;
    unsigned int totalsourceskipgrams = 0;
    unsigned int totaltargetskipgrams = 0;
    int maxn = 0;
    cerr << "N\tsource-ngrams\tsource-skipgrams\ttarget-ngrams\ttarget-skipgrams" << endl;
    cerr << "-------------------------------------------------------------------------------------------------------" << endl;
    for (int i = 0; i <= 9; i++) {
        if (sourcengramusage[i] ||sourceskipgramusage[i] || targetngramusage[i]  || targetskipgramusage[i]) {
            maxn = i;
        }
        totalsourcengrams += sourcengramusage[i];
        totalsourceskipgrams += sourceskipgramusage[i];
        totaltargetngrams += targetngramusage[i];
        totaltargetskipgrams += targetskipgramusage[i];
    }
    for (int i = 0; i <= maxn; i++) {
        double p;                    
        if (sourcengramusage[i] ||  sourceskipgramusage[i] || targetngramusage[i]  || targetskipgramusage[i]) {                    
            cerr << i << "\t";
            if (totalsourcengrams > 0) {
                p = (double) sourcengramusage[i] / totalsourcengrams;
            } else {
                p = 0;
            }
            cerr << sourcengramusage[i] << " (" << p << ")\t";
            
            if (totalsourceskipgrams > 0) {            
                p = (double)  sourceskipgramusage[i] / totalsourceskipgrams;
            } else {
                p = 0;
            }                                                
            cerr << sourceskipgramusage[i] << " (" << p << ")\t";
            
            if (totaltargetngrams > 0) {
                p = (double)  targetngramusage[i] / totaltargetngrams;
            } else {
                p = 0;
            }
            cerr << targetngramusage[i] << " (" << p << ")\t";
            
            if (totaltargetskipgrams > 0) {
                p = (double)  targetskipgramusage[i] / totaltargetskipgrams;
            } else {
                p = 0;
            }                                                
            cerr << targetskipgramusage[i] << " (" << p << ")" << endl;                        
        }                     
    }
        cerr << "-------------------------------------------------------------------------------------------------------" << endl;
        cerr << " \t" << totalsourcengrams << "\t" << totalsourceskipgrams << "\t" << totaltargetngrams << "\t" << totaltargetskipgrams << endl;    
}


int main( int argc, char *argv[] ) {
    int MOSESFORMAT = 0;
    int STATS = 0;
    int GLOBALSTATS = 0;
    vector<double> tweights;
    double lweight = 1.0;
    double dweight = 1.0;
    string transtablefile = "";
    string lmfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    int stacksize = 10;
    int dlimit = 999;
    double prunethreshold = 0.5;
    int maxn = 9;
    static struct option long_options[] = {      
       {"moses", no_argument,             &MOSESFORMAT, 1},
       {"stats", no_argument,             &STATS, 1},
       {"globalstats", no_argument,             &GLOBALSTATS, 1},                       
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    bool DOSKIPGRAMS = true;
    string timbloptions = "-a 1";
    string classifierid = "";
    ScoreHandling scorehandling = SCOREHANDLING_WEIGHED;
    
    //temp:
    string raw;
    
    string ws;
    stringstream linestream;
    string s;
    
    int debug = 0;
    char c;    
    while ((c = getopt_long(argc, argv, "ht:S:T:s:p:l:W:L:D:v:M:NC:x:O:",long_options,&option_index)) != -1) {
        switch (c) {
        case 0:
            if (long_options[option_index].flag != 0)
               break;
        case 'h':
        	usage();
        	exit(0);
        case 'l':
            lmfile = optarg;
            break;       
        case 't':
            transtablefile = optarg;
            break;
        case 'S':
            sourceclassfile = optarg;
            break;
        case 'T':
            targetclassfile = optarg;
            break;
        case 'L':
            lweight = atof(optarg);
            break;
        case 'D':
            dweight = atof(optarg);
            break;
        case 'M':
            dlimit = atoi(optarg);
            break;
        case 's':
            stacksize = atoi(optarg);
            break;
        case 'p':
            prunethreshold = atof(optarg);
            break;            
        case 'W':
            tweights.push_back(atof(optarg));
            break;       
        case 'v':
            debug = atoi(optarg);
            break;     
        case 'N':
            DOSKIPGRAMS = false;
            break;
        case 'C':
            classifierid = optarg;
            break;
        case 'O':
            timbloptions = optarg;
            break;            
        case 'x':
            s = optarg;
            if (s == "weighed") {
                scorehandling = SCOREHANDLING_WEIGHED;
            } else if (s == "append") {
                scorehandling = SCOREHANDLING_APPEND;
            } else if (s == "replace") {
                scorehandling = SCOREHANDLING_REPLACE;
            } else if (s == "ignore") {
                scorehandling = SCOREHANDLING_IGNORE;                
            } else {
                cerr << "Invalid value for -x: '" << s << "'" << endl;
                exit(2);    
            }
            break;        
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }        
    }
    
    if (tweights.empty()) {
        cerr << "WARNING: No translation weights specified, assuming default of 1.0,1.0. May clash with translation table!" << endl;
        tweights.push_back(1.0);
        tweights.push_back(1.0);        
    }
    
    if (transtablefile.empty()) {
        cerr << "ERROR: No translation table specified (-t)" << endl;
        exit(2);
    }
    if (lmfile.empty()) {
        cerr << "ERROR: No language model specified (-l)" << endl;
        exit(2);
    }    
    if (sourceclassfile.empty()) {
        cerr << "ERROR: No source class file specified (-S)" << endl;
        exit(2);
    }
    if (targetclassfile.empty()) {
        cerr << "ERROR: No target class file specified (-T)" << endl;
        exit(2);
    }        
    
    cerr << "Colibri MT Decoder" << endl;
    cerr << "   Maarten van Gompel, Radboud University Nijmegen" << endl;
    cerr << "----------------------------------------------------" << endl;
    cerr << "Translation weights:  ";
    for (unsigned int i = 0; i < tweights.size(); i++) cerr << tweights[i] << " "; 
    cerr << endl;
    cerr << "Distortion weight:    " << dweight << endl;
    cerr << "LM weight:            " << lweight << endl;    
    cerr << "Stacksize:            " << stacksize << endl;
    cerr << "Prune threshold:      " << prunethreshold << endl;
    cerr << "Distortion limit:     " << dlimit << endl;
    cerr << "----------------------------------------------------" << endl;
    cerr << "Source classes:       " << sourceclassfile << endl;
    ClassEncoder sourceclassencoder = ClassEncoder(sourceclassfile);
    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);       
    cerr << "Target classes:       " << targetclassfile << endl;
    ClassEncoder targetclassencoder = ClassEncoder(targetclassfile);    
    ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);
    addsentencemarkers(targetclassdecoder, targetclassencoder);        
    cerr << "Language model:       " << lmfile << endl;
    LanguageModel lm = LanguageModel(lmfile, targetclassencoder, &targetclassdecoder, (debug >= 4) );
    cerr << "   loaded " << lm.size() << " n-grams, order=" << lm.getorder() << endl;
    
    
    cerr << "Translation table:    " << transtablefile << endl;
    
    AlignmentModel * transtable;
    if (MOSESFORMAT) {
        transtable = new AlignmentModel(transtablefile, &sourceclassencoder, &targetclassencoder);
    } else {
        transtable = new AlignmentModel(transtablefile, true, DOSKIPGRAMS);        
    }
    
    cerr << "   loaded translations for " << transtable->size() << " patterns" << endl;
    
    ClassifierInterface * classifier = NULL;
    if (!classifierid.empty()) {
        //cerr << "Computing reverse index for translation table" << endl;
        //transtable->computereverse(); //not necessary 
        cerr << "Loading classifiers" << endl;
        cerr << "   ID: " << classifierid << endl;
        cerr << "   Timbl options: " << timbloptions << endl;
        cerr << "   Score handling: ";
        if (scorehandling == SCOREHANDLING_WEIGHED) {
            cerr << "weighed" << endl;
        } else if (scorehandling == SCOREHANDLING_APPEND) {
            cerr << "append" << endl;
        } else if (scorehandling == SCOREHANDLING_REPLACE) {
            cerr << "replace" << endl;
        }
        bool exemplarweights = true; //read from config
        classifier = (ClassifierInterface*) new NClassifierArray(classifierid, (int) transtable->leftsourcecontext, (int) transtable->rightsourcecontext);
        classifier->load(timbloptions, &sourceclassdecoder, &targetclassencoder, debug);        
    }   
        
    //const int firstunknownclass_source = sourceclassencoder.gethighestclass()+1;    
    //const int firstunknownclass_target = targetclassencoder.gethighestclass()+1;

    DecodeStats overallstats;
    
    string input;
    unsigned char buffer[8192]; 
    int size;
    while (getline(cin, input)) {        
        if (input.length() > 0) {                    
            cerr << "INPUT: " << input << endl;
            if (debug >= 1) cerr << "Processing input" << endl;        
            size = sourceclassencoder.encodestring(input, buffer, true, true) - 1; //weird patch: - 1  to get n() right later               
            const EncData * const inputdata = new EncData(buffer,size);
            if (debug >= 1) cerr << "Processing unknown words" << endl; 
            addunknownwords(*transtable, lm, sourceclassencoder, sourceclassdecoder, targetclassencoder, targetclassdecoder, tweights.size());
            if (debug >= 1) cerr << "Setting up decoder" << endl;
            //TODO: add classifier support
            StackDecoder * decoder = new StackDecoder(*inputdata, transtable, &lm, stacksize, prunethreshold, tweights, dweight, lweight, dlimit, maxn, debug, &sourceclassdecoder, &targetclassdecoder, classifier, scorehandling, (bool) GLOBALSTATS);
            if (debug >= 1) cerr << "Decoding..." << endl;
            TranslationHypothesis * solution = decoder->decode();                    
            if (solution != NULL) {
                EncData s = solution->getoutput();
                if (decoder->DEBUG >= 1) {
                    TranslationHypothesis * h = solution;            
                    cerr << "HISTORY=" << endl;
                    while (h != NULL) {
                        cerr << "\t";
                        for (unsigned int i = 0; i < decoder->inputlength; i++) {
                            cerr << h->inputcoveragemask[i];
                        }
                        cerr << "  [" << h->score() << "]" << endl;
                        h = h->parent;
                    }                
                }
                cerr << "GAPLESS STACKSIZES: ";
                for (map<int,int>::iterator iter = decoder->stats.stacksizes.begin(); iter != decoder->stats.stacksizes.end(); iter++) {
                    if (iter->first > 0) cerr << iter->second << ' ';
                }
                cerr << endl;
                if (DOSKIPGRAMS) {
                    cerr << "GAPPY STACKSIZES: ";
                    for (map<int,int>::iterator iter = decoder->stats.gappystacksizes.begin(); iter != decoder->stats.gappystacksizes.end(); iter++) {
                        if (iter->first > 0) cerr << iter->second << ' ';
                    }
                }
                cerr << endl;                
                cerr << "STATS:\tTotal expansions: " << decoder->stats.expanded << endl;
                cerr << "       \tof which rejected: " << decoder->stats.discarded << endl;
                cerr << "       \tof which pruned: " << decoder->stats.pruned << endl;
                cerr << "       \tof which gapresolutions: " << decoder->stats.gapresolutions << endl;
                cerr << "SCORE=" << solution->score() << endl;                
                cerr << "DONE. OUTPUT:" << endl;
                cout << s.decode(targetclassdecoder) << endl;
                if ((STATS) && (!GLOBALSTATS)) solution->stats(); 
                if (decoder->DEBUG == 99) {
                    cerr << "DEBUG: SOLUTION DESTRUCTION. DELETING " << (size_t) solution << endl;
                    solution->cleanup();
                } else {
                    delete solution;
                }
            } else {
                cerr << "ERROR: NO SOLUTION FOUND!!!" << endl;
                exit(12);
            }                
            

            
            
            //output statistics
            if (STATS) {
                decoder->stats.output();
                overallstats.add(decoder->stats);
            }
            
            
            //delete inputdata; //TODO: REENABLE, MEMORY LEAK
            if (decoder->DEBUG == 99) cerr << "DEALLOCATING DECODER" << endl;
            delete decoder;
       }
    }
    
    if (STATS) {
        cerr << "==========================================================================" << endl;
        cerr << "OVERALL STATISTICS FOR ALL SENTENCES:" << endl;
        overallstats.output();
        cerr << endl; 
    }
    
    delete transtable;
}
