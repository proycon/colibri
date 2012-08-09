#include <decoder.h>
#include <algorithm> 
#include <getopt.h>
#include <iostream>
#include <sstream>


using namespace std;

unsigned char UNKNOWNCLASS = 2;
unsigned char BOSCLASS = 3;
unsigned char EOSCLASS = 4;

StackDecoder::StackDecoder(const EncData & input, TranslationTable * translationtable, LanguageModel * lm, int stacksize, double prunethreshold, vector<double> tweights, double dweight, double lweight, int maxn, int debug, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) {
        this->input = input;
        this->inputlength = input.length();
        this->translationtable = translationtable;
        this->lm = lm;
        this->stacksize = stacksize;
        this->prunethreshold = prunethreshold;
        this->DEBUG = debug;
        this->sourceclassdecoder = sourceclassdecoder;
        this->targetclassdecoder = targetclassdecoder;        
        
        
        //init stacks
        for (int i = 0; i <= inputlength; i++) {
            stacks.push_back( Stack(this, i, stacksize, prunethreshold) );
        }

        if (DEBUG >= 3) cerr << "Gathering source fragments:" << endl;        
        sourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn);
        
        //Build a coverage mask, this will be used to check if their are words uncoverable by translations, these will be added as unknown words 
        std::vector<bool> inputcoveragemask;
        inputcoveragemask.reserve(inputlength);
        for (int i = 0; i < inputlength; i++) inputcoveragemask.push_back(false);                 
        for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {           
            const EncAnyGram * anygram = iter->first;
            int n;
            if (anygram->isskipgram()) {
                n = ((const EncSkipGram*) anygram)->n();
            } else {
                n = (int) ((const EncNGram*) anygram)->n(); 
            }
            for (int j = iter->second.token; j < iter->second.token + n; j++) inputcoveragemask[j] = true;
            
            //Output translation options
            if ((DEBUG >= 3) && (sourceclassdecoder != NULL) && (targetclassdecoder != NULL)) {
                if (anygram->isskipgram()) {
                    cerr << "\t" << (int) iter->second.token << ':' << n << " -- " << ((const EncSkipGram*) anygram)->decode(*sourceclassdecoder) << " ==> ";                    
                } else {
                    cerr << "\t" << (int) iter->second.token << ':' << n << " -- " << ((const EncNGram*) anygram)->decode(*sourceclassdecoder) << " ==> ";
                }
                const EncAnyGram * sourcekey = translationtable->getsourcekey(anygram);
                if (sourcekey == NULL) {
                    cerr << endl;
                    cerr << "ERROR: No translation options found!!!!" << endl;
                    exit(6);
                } else {
                    for (std::unordered_map<const EncAnyGram*, std::vector<double> >::iterator iter2 = translationtable->alignmatrix[sourcekey].begin(); iter2 != translationtable->alignmatrix[sourcekey].end(); iter2++) {
                        const EncAnyGram * anygram2 = iter2->first;
                        if (anygram2->isskipgram()) {
                            cerr << ((const EncSkipGram*) anygram2)->decode(*targetclassdecoder) << " [ ";                            
                        } else {
                            cerr <<  ((const EncNGram*) anygram2)->decode(*targetclassdecoder) << " [ ";
                        }                        
                        for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                            cerr << *iter3 << " ";
                        }
                        cerr << "]; ";                        
                    }
                }
                cerr << endl;            
            }
        }
        
        //Check for uncoverable words
        for (int i = 0; i < inputlength; i++) {
            if (!inputcoveragemask[i]) {
                //found one
                
                EncNGram * unigram = input.slice(i, 1);
                const string word = unigram->decode(*sourceclassdecoder);
                cerr << "NOTICE: UNTRANSLATABLE WORD: '" << word << "' (adding)" << endl;  
                                
                targetclassdecoder->add(targetclassdecoder->gethighestclass() + 1, word);
                lm->ngrams[*unigram] = -99; //TODO: more sane LM value
                
                const EncAnyGram * sourcekey = translationtable->getsourcekey((const EncAnyGram *) unigram);                
                if (sourcekey == NULL) {
                    translationtable->sourcengrams.insert(*unigram);
                    sourcekey = translationtable->getsourcekey((const EncAnyGram *) unigram);
                }
                const EncAnyGram * targetkey = translationtable->gettargetkey((const EncAnyGram *) unigram);
                if (targetkey == NULL) {
                    translationtable->targetngrams.insert(*unigram);
                    targetkey = translationtable->gettargetkey((const EncAnyGram *) unigram);
                }
                vector<double> scores;
                for (int j = 0; j < tweights.size(); j++) scores.push_back(1);                 
                translationtable->alignmatrix[sourcekey][targetkey] = scores;
                                                
                sourcefragments.push_back(make_pair( sourcekey, CorpusReference(0,i) ));
                delete unigram; 
                
                if (DEBUG >= 3) {
                    cerr << "\t" << i << ":1" << " -- " << word << " ==> " << word << " [ 1 1 ];" << endl;                    
                }
            }
        }
        
        
        this->tweights = vector<double>(tweights.begin(), tweights.end());
        this->dweight = dweight;
        this->lweight = lweight;
        
        computefuturecost();
        
        

}


void StackDecoder::computefuturecost() {
        map<pair<int,int>, double> sourcefragments_costbyspan;
        //reorder source fragments by span for more efficiency
        for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {
            
            const EncAnyGram * anygram = iter->first;
            const CorpusReference ref = iter->second;
            const int n = anygram->n();    
            const pair<int,int> span = make_pair<int,int>((int) ref.token, (int) n);
             
            //cerr << "DEBUG: " << span.first << ':' << span.second << endl;
             
            const EncAnyGram * candidate = translationtable->getsourcekey(anygram);
            if (translationtable->alignmatrix[candidate].size() == 0) {
                    cerr << "INTERNAL ERROR: No translation options" << endl;
                    exit(6);  
            }
            
            //find cheapest translation option
            double bestscore = -INFINITY;            
            for (std::unordered_map<const EncAnyGram*, std::vector<double> >::iterator iter2 = translationtable->alignmatrix[candidate].begin(); iter2 != translationtable->alignmatrix[candidate].end(); iter2++) {
                if (tweights.size() > iter2->second.size()) {
                    cerr << "Too few translation scores specified for an entry in the translation table. Expected at least "  << tweights.size() << ", got " << iter2->second.size() << endl;
                    exit(6);  
                }
                double score = 0; 
                for (int i = 0; i < tweights.size(); i++) {
                    double p = iter2->second[i];
                    if (p > 0) p = log10(p); //turn into logprob, base10 
                    score += tweights[i] * p;
                }
                const EncAnyGram * translationoption = iter2->first;
                if (translationoption->isskipgram()) {
                    vector<EncNGram*> parts;
                    (*((const EncSkipGram *) translationoption)).parts(parts);
                    for (vector<EncNGram*>::iterator iter = parts.begin(); iter != parts.end(); iter++) {
                        EncNGram * part = *iter; 
                        score += lweight * lm->score(*part);
                        delete part;
                    }  
                } else {
                    const EncNGram * ngram = (const EncNGram *) translationoption;
                    score += lweight * lm->score(*ngram);
                }
                if (score > bestscore) {
                    bestscore = score;
                }                            
            } 
            sourcefragments_costbyspan[span] = bestscore; 
        }   
        
        //compute future cost
        for (int length = 1; length <= inputlength; length++) {
            for (int start = 0; start < inputlength - length + 1; start++) {
                const pair<int,int> span = make_pair((int) start,(int) length);
                bool found = false;
                for (map<pair<int,int>, double>::iterator iter = sourcefragments_costbyspan.find(span); iter != sourcefragments_costbyspan.end(); iter++) {
                    found = true;
                    futurecost[span] = sourcefragments_costbyspan[span];
                }
                if (!found) {
                    if (length == 1){                       
                        cerr << "INTERNAL ERROR: No sourcefragment covers " << span.first << ":" << span.second << " ! Unable to compute future cost!" << endl;
                        exit(6);
                    } else {
                        futurecost[span] = -INFINITY;
                    }
                }
                for (int i = 1; i < length; i++) {
                    double spanscore = futurecost[make_pair((int) start,(int) i)] + futurecost[make_pair((int) start+i,(int) length - i)];
                    if (spanscore > futurecost[span]) { //(higher score -> lower cost)
                        futurecost[span] = spanscore;
                    }
                }
            }
        }
        
        if (DEBUG >= 3) {
            cerr << "\tFuture cost precomputation:" << endl;
            for (map<std::pair<int, int>, double>::iterator iter = futurecost.begin(); iter != futurecost.end(); iter++) {
                cerr << "\t  " << iter->first.first << ":" << iter->first.second << " = " << iter->second << endl;                
            }            
        }
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
    bool dead = false;
    
    //for each stack
    for (int i = 0; i <= inputlength - 1; i++) {
        if (DEBUG >= 1) {
            cerr << "\tDecoding Stack " << i << " -- " << stacks[i].size() << " hypotheses" << endl;
        }        
        unsigned int totalexpanded = 0;
        bool first = true;        
        while (!stacks[i].empty()) {
            if (DEBUG >= 1) cerr << "\t Expanding hypothesis off stack " << i << " -- " << stacks[i].size() -1 << " left:" << endl;
            
            //pop from stacks[i]
            TranslationHypothesis * hyp = stacks[i].pop();
            if (first) {                
                fallbackhyp = hyp;
                fallbackhyp->keep = true; //prevent fallback hypothesis from getting pruned
                first = false;
                cerr << "DEBUG: FALLBACKHYP=" << (size_t) fallbackhyp << endl;
            }
            if (DEBUG >= 2) {
                cerr << "\t  Popped from stack:" << endl;
                hyp->report();
                cerr << "\t  Expands to:" << endl;
            }
                         
            bool finalonly = (i == inputlength - 1);
            unsigned int expanded = hyp->expand(finalonly); //will automatically add to appropriate stacks
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
        if ((totalexpanded == 0) && (i != inputlength)) {
            dead = true;
            for (int j = i + 1; j <= inputlength; j++) {
                if (!stacks[j].empty()) dead = false;
            }
            if (dead) {
                cerr << "DECODER ENDED PREMATURELY AFTER STACK " << i << ", NO FURTHER EXPANSIONS POSSIBLE." << endl;
                break;
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
        if (DEBUG >= 1) cerr << "\t Expanded " << totalexpanded << " new hypotheses in total after processing stack " << i << endl;
        unsigned int totalpruned = 0;
        for (int j = inputlength; j >= i+1; j--) { //prune further stacks (hypotheses may have been added to any of them).. always prune superior stack first
            unsigned int pruned = stacks[j].prune();
            if ((DEBUG >= 1) && (pruned > 0)) cerr << "\t  Pruned " << pruned << " hypotheses from stack " << j << endl;
            totalpruned += pruned; 
        }
        if (DEBUG >= 1) cerr << "\t Pruned " << totalpruned << " hypotheses from all superior stacks" << endl;
        if (!dead) {
            fallbackhyp = NULL;
            stacks[i].clear(); //stack is in itself no longer necessary, included pointer elements may live on though! Unnecessary hypotheses will be cleaned up automatically when higher-order hypotheses are deleted
        }
    }   
    
    //solutions are now in last stack: stacks[inputlength]       
     
    TranslationHypothesis * solution;
    if (!stacks[inputlength].empty()) {
        TranslationHypothesis * solution = stacks[inputlength].pop();
        return solution;
    } else if (dead) {
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
    if (contents.size() >= stacksize) {
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
    if ((!added) && (contents.size() < stacksize)) {
        added = true;
        contents.push_back(candidate);
    }
    return added;    
}


int Stack::prune() {
    int pruned = 0;
    if ((prunethreshold != 1) && (prunethreshold != 0)) {
        //pruning based on prunethreshold
        double cutoff = bestscore() / prunethreshold;
        for (list<TranslationHypothesis*>::iterator iter = contents.begin(); iter != contents.end(); iter++) {
            TranslationHypothesis*  h = *iter;
            if (h->score() < cutoff) {
                pruned++;
                iter = contents.erase(iter);
                if (h->deletable()) {
                    if (decoder->DEBUG) {
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
    
    if ((targetgram != NULL) && (targetgram->isskipgram())) {
        vector<pair<int, int> > gaps;
        (*((const EncSkipGram *) targetgram)).getgaps(gaps);
        for (vector<pair<int, int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) targetgaps.push_back( make_pair<unsigned char, unsigned char>(iter->first, iter->second) );
    }    
    
       
    
    //compute input coverage
    if (parent == NULL) {
        //initial hypothesis: no coverage at all
        for (int i = 0; i < decoder->inputlength; i++) inputcoveragemask.push_back(false);
    } else {
        //inherit from parent
        for (int i = 0; i < decoder->inputlength; i++) {
             inputcoveragemask.push_back(parent->inputcoveragemask[i]);
        }
        //add sourcegram coverage
        bool isskipgram = sourcegram->isskipgram();
        for (int i = sourceoffset; i < sourceoffset + sourcegram->n(); i++) {
            if (isskipgram) {
                bool ingap = false;
                for (vector<pair<unsigned char, unsigned char> >::iterator iter = sourcegaps.begin(); iter < sourcegaps.end(); iter++) {
                    if ( (i + sourcegram->n() > sourceoffset + iter->first) &&  ( i < sourceoffset +iter->first + iter->second ) ) {
                        ingap = true;
                        break;       
                    }             
                }
                if (!ingap) {
                    inputcoveragemask[i] = true;
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
    if (begin < 0) begin = -1; //TODO : CHECK
    
    for (int i = begin; i < begin + (order-1); i++) {
        EncNGram * unigram;
        if (begin == -1) {
            unigram = new EncNGram((const unsigned char*) &BOSCLASS, 1); //TODO: CHECK
        } else { 
            unigram = getoutputtoken(i);
        }            
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
        cerr << "Too few translation scores specified for an entry in the translation table. Expected at least "  << decoder->tweights.size() << ", got " << tscores.size() << endl;
        exit(6);  
    }
  
    tscore = 0; 
    for (int i = 0; i < tscores.size(); i++) {
        double p = tscores[i];
        if (p > 0) p = log10(p); //turn into logprob, base10 
        tscore += decoder->tweights[i] * p;
    }
    
    
    lmscore = 0;
    
    if (parent != NULL) {
        if (history != NULL) {
            if (targetgram->isskipgram()) {
                vector<EncNGram*> parts;
                ((const EncSkipGram *) targetgram)->parts(parts);
                int partcount = 0;
                for (vector<EncNGram*>::iterator iter = parts.begin(); iter != parts.end(); iter++) {
                    EncNGram * part = *iter;
                    if ((partcount == 0) && (sourcegaps[0].first != 0)) {
                        //first part, no initial gap, join with history
                        EncNGram ngram = *history + *part;
                        lmscore += decoder->lweight * decoder->lm->score(ngram);
                    } else {     
                        lmscore += decoder->lweight * decoder->lm->score(*part);
                    }
                    partcount++;
                    delete part;
                } 
            } else {
                EncNGram ngram = EncNGram(*history + *((const EncNGram* ) targetgram) );
                lmscore += decoder->lweight * decoder->lm->score(ngram);
            }
        } else {
            if (targetgram->isskipgram()) {
                vector<EncNGram*> parts;
                ( (const EncSkipGram *) targetgram)->parts(parts);
                for (vector<EncNGram*>::iterator iter = parts.begin(); iter != parts.end(); iter++) {
                    EncNGram * part = *iter; 
                    lmscore += decoder->lweight * decoder->lm->score(*part);
                    delete part;
                } 
            } else {
               lmscore += decoder->lweight * decoder->lm->score( *((const EncNGram *) targetgram));
            }
        }
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
    for (int i = 0; i <= inputcoveragemask.size() ; i++) {
        if ((!inputcoveragemask[i]) && (begin == -1) && (i < inputcoveragemask.size())) {
            begin = i;
        } else if (((i == inputcoveragemask.size()) || (inputcoveragemask[i])) && (begin != -1)) {
            const double c = decoder->futurecost[make_pair((int) begin,(int) i-begin)];
            if (c == 0) {
                cerr << "INTERNAL ERROR: Future cost for " << begin << ":" << i - begin << " is 0! Not possible!" << endl;
                report();
                exit(6);
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
        cerr << "\t    fragmentscore = tscore + lmscore + dscore = " << tscore << " + " << lmscore << " + " << dscore << " = " << _score << endl;
        cerr << "\t    futurecost = " << futurecost << endl;        
        cerr << "\t    totalscore = allfragmentscores + futurecost = " << score() << endl;
        cerr << "\t    coverage: ";
        for (int i = 0; i < inputcoveragemask.size(); i++) {
            if (inputcoveragemask[i]) {
                cerr << "1";
            } else {
                cerr << "0";
            }            
        }
        cerr << endl;
}

double TranslationHypothesis::score() const {
    double s = 0;
    const TranslationHypothesis * h = this;
    while (h != NULL) {
        s += h->tscore + h->lmscore + h->dscore;
        h = h->parent;
    } 
    s += futurecost;
    return s;    
} 


bool TranslationHypothesis::deletable() {
    return ((!keep) && (children.empty()));
}

void TranslationHypothesis::cleanup() {  
    if (decoder->DEBUG == 99) {
        cerr << "DEBUG: DELETING HYPOTHESIS " << (size_t) this << endl;
        if (deleted) {
             cerr << "INTERNAL ERROR: DELETING AN ALREADY DELETED HYPOTHESIS!!! THIS SHOULD NOT HAPPEN!!!!!" << (size_t) this << endl;
             exit(6);
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

bool TranslationHypothesis::final(){
        if (!targetgaps.empty()) return false;
        return (inputcoverage() == decoder->inputlength); 
}


unsigned int TranslationHypothesis::expand(bool finalonly) {
    bool oldkeep = this->keep;
    this->keep = true; //lock this hypothesis, preventing it from being deleted when expanding and rejecting its last dying child
    if (deleted) { 
        //only in debug==99
        cerr << "ERROR: Expanding an already deleted hypothesis! This should never happen! " << (size_t) this << endl;
        exit(6);        
    }
    
    
    unsigned int expanded = 0;
    int thiscov = inputcoverage();
    //expand directly in decoder.stack()
    if (targetgaps.empty()) {
        //find new source fragment
        for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = decoder->sourcefragments.begin(); iter != decoder->sourcefragments.end(); iter++) {
            const EncAnyGram * sourcecandidate = iter->first;
            const CorpusReference ref = iter->second; 
            if (!conflicts(sourcecandidate, ref)) {
                //find target fragments for this source fragment
                if (decoder->translationtable->alignmatrix.find(sourcecandidate) == decoder->translationtable->alignmatrix.end()) {
                    cerr << "ERROR: Source candidate not found in translation table. This should never happen!" << endl;
                    exit(6);
                }
                int c = 0;                
                for (std::unordered_map<const EncAnyGram*, vector<double> >::const_iterator iter2 =  decoder->translationtable->alignmatrix[sourcecandidate].begin(); iter2 != decoder->translationtable->alignmatrix[sourcecandidate].end(); iter2++) {
                    c++;
                    //cerr << "DEBUG: " << c << " of " << decoder->translationtable->alignmatrix[sourcecandidate].size() << endl;  
                    //create hypothesis for each target fragment
                    const EncAnyGram * targetcandidate = iter2->first;
                    int length;
                    if (targetgram != NULL) { 
                        length = targetgram->n();
                    } else {    
                        length = 0;
                    }
                    TranslationHypothesis * newhypo = new TranslationHypothesis(this, decoder, sourcecandidate, ref.token, targetcandidate, targetoffset + length , iter2->second);
                    if ((finalonly) && (!newhypo->final())) {
                        if (!newhypo->deletable()) {
                            cerr << "INTERNAL ERROR: Newly created hypothesis not deletable? shouldn't happen" << endl;
                            exit(6);
                        }
                        if (decoder->DEBUG == 99) {
                            cerr << "DEBUG: IMMEDIATELY DELETING NEWLY CREATED HYPOTHESIS (NOT FINAL) " << (size_t) newhypo << endl;
                            newhypo->cleanup();
                        } else {
                            delete newhypo;
                        }
                    } else {
                        //add to proper stack
                        int cov = newhypo->inputcoverage();
                        if (thiscov >= cov) {
                            cerr << "INTERNAL ERROR: Hypothesis expansion did not lead to coverage expansion! This should not happen. New hypo has coverage " << cov << ", parent: " << thiscov << endl;
                            exit(6);        
                        }
                        if (decoder->DEBUG >= 2) cerr << "\t    Adding to stack " << cov;
                        bool accepted = decoder->stacks[cov].add(newhypo);
                        if (decoder->DEBUG >= 2) {
                            if (accepted) {
                                cerr << " ... ACCEPTED" << endl;
                            } else {
                                cerr << " ... REJECTED" << endl;
                            }
                        }
                        expanded++;
                        
                        if (!accepted) {
                            if (decoder->DEBUG == 99) {
                                cerr << "DEBUG: IMMEDIATELY DELETING NEWLY CREATED HYPOTHESIS (REJECTED BY STACK) " << (size_t) newhypo << endl;
                                newhypo->cleanup();
                            } else {
                                delete newhypo;
                            } 
                        }
                    }
                }                 
            }        
        }
    } else {
        //there are target-side gaps, attempt to fill
        //TODO
    }
    this->keep = oldkeep; //release lock
    return expanded;
}

bool TranslationHypothesis::conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref) {
    if ((sourcegram == NULL) && (parent == NULL)) return false; //no sourcegram, no conflict (this is an empty initial hypothesis)
    
    if (sourcecandidate->hash() == sourcegram->hash()) return true; //source was already added, can not add twice
     
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
        return parent->conflicts(sourcecandidate, ref);
    } else {
        return false;
    }
}



int TranslationHypothesis::inputcoverage() {
    int c = 0;
    for (int i = 0; i < inputcoveragemask.size(); i++) {
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
            exit(6);
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
            int index = iter->first;
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
    cerr << "\t-W w1,w2,w3               Translation model weights (comma seperated)" << endl;    
    cerr << "\t-L weight                 Language model weight" << endl;
    cerr << "\t-D weight                 Distortion model weight" << endl;      
    cerr << "\t--moses                   Translation table is in Moses format" << endl;
    cerr << "\t-d debug                  Debug level" << endl;
    
}

void addsentencemarkers(ClassDecoder & targetclassdecoder, ClassEncoder & targetclassencoder) {
    targetclassdecoder.add(BOSCLASS,"<s>");
    targetclassdecoder.add(EOSCLASS,"</s>");
    targetclassencoder.add("<s>", BOSCLASS);
    targetclassencoder.add("</s>", EOSCLASS);
}

int addunknownwords( TranslationTable & ttable, LanguageModel & lm, ClassEncoder & sourceclassencoder, ClassDecoder & sourceclassdecoder,  ClassEncoder & targetclassencoder, ClassDecoder & targetclassdecoder, int tweights_size) {
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
            
            ttable.sourcengrams.insert(sourcegram);
            ttable.targetngrams.insert(targetgram);
            
            
            
            lm.ngrams[targetgram] = -99; //TODO: Use better unknown value from LM?
            
            
            
            const EncAnyGram * sourcekey = ttable.getsourcekey((const EncAnyGram*) &sourcegram );
            const EncAnyGram * targetkey = ttable.gettargetkey((const EncAnyGram*) &targetgram );
            
            vector<double> scores;
            for (int j = 0; j < tweights_size; j++) scores.push_back(1);
            ttable.alignmatrix[sourcekey][targetkey] = scores;
            
            
        }
    }
    sourceclassencoder.added.clear();
    return added;    
}


int main( int argc, char *argv[] ) {
    int MOSESFORMAT = 0;
    vector<double> tweights;
    double lweight = 1.0;
    double dweight = 1.0;
    string transtablefile = "";
    string lmfile = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    int stacksize = 10;
    double prunethreshold = 0.5;
    int maxn = 9;
    static struct option long_options[] = {      
       {"moses", no_argument,             &MOSESFORMAT, 1},                       
       {0, 0, 0, 0}
     };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    
    //temp:
    string raw;
    string::size_type pos;
    
    string ws;
    stringstream linestream;
    double w;
    
    int debug = 0;
    char c;    
    while ((c = getopt_long(argc, argv, "ht:S:T:s:p:l:W:L:D:d:",long_options,&option_index)) != -1) {
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
        case 's':
            stacksize = atoi(optarg);
            break;
        case 'p':
            prunethreshold = atof(optarg);
            break;            
        case 'W':
            ws = optarg;
            while (linestream >> w) tweights.push_back(w);
            break;       
        case 'd':
            debug = atoi(optarg);
            break;     
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }        
    }
    
    if (tweights.empty()) {
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
    for (int i = 0; i < tweights.size(); i++) cerr << tweights[i] << " "; 
    cerr << endl;
    cerr << "Distortion weight:    " << dweight << endl;
    cerr << "LM weight:            " << lweight << endl;    
    cerr << "Stacksize:            " << stacksize << endl;
    cerr << "Prune threshold:      " << prunethreshold << endl;
    cerr << "----------------------------------------------------" << endl;
    cerr << "Source classes:       " << sourceclassfile << endl;
    ClassEncoder sourceclassencoder = ClassEncoder(sourceclassfile);
    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);       
    cerr << "Target classes:       " << targetclassfile << endl;
    ClassEncoder targetclassencoder = ClassEncoder(targetclassfile);    
    ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);
    addsentencemarkers(targetclassdecoder, targetclassencoder);        
    cerr << "Language model:       " << lmfile << endl;
    LanguageModel lm = LanguageModel(lmfile, targetclassencoder);
    cerr << "   loaded " << lm.size() << " n-grams, order=" << lm.getorder() << endl;
    cerr << "Translation table:    " << transtablefile << endl;
    //TODO: Moses format
    TranslationTable transtable = TranslationTable(transtablefile);
    cerr << "   loaded translations for " << transtable.size() << " patterns" << endl;
        
    const int firstunknownclass_source = sourceclassencoder.gethighestclass()+1;    
    const int firstunknownclass_target = targetclassencoder.gethighestclass()+1;

    
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
            addunknownwords(transtable, lm, sourceclassencoder, sourceclassdecoder, targetclassencoder, targetclassdecoder, tweights.size());
            if (debug >= 1) cerr << "Setting up decoder" << endl;
            StackDecoder * decoder = new StackDecoder(*inputdata, &transtable, &lm, stacksize, prunethreshold, tweights, dweight, lweight, maxn, debug, &sourceclassdecoder, &targetclassdecoder);
            if (debug >= 1) cerr << "Decoding..." << endl;
            TranslationHypothesis * solution = decoder->decode();
            cerr << "DONE. OUTPUT:" << endl;        
            if (solution != NULL) {
                EncData s = solution->getoutput();
                cerr << "SCORE=" << solution->score() << endl;
                cout << s.decode(targetclassdecoder) << endl;
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
            //delete inputdata; //TODO: REENABLE, MEMORY LEAK
            if (decoder->DEBUG == 99) cerr << "DEALLOCATING DECODER" << endl;
            delete decoder;
       }
    }
}
