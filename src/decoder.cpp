unsigned char UNKNOWNCLASS = 2;

StackDecoder::StackDecoder(const EncData & input, TranslationTable * translationtable, int stacksize, double prunethreshold, int maxn) {
        this->input = input;
        this->inputlength = input.length();
        this->translationtable = translationtable;
        this->stacksize = stacksize;
        this->prunethreshold = prunethreshold;
        sourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn);
        
        //TODO: Deal with unknown tokens?
}
StackDecoder::setdebug(int debug) {
    this->DEBUG = debug;
}

    
void StackDecoder::decode() {p
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
    const TranslationHypothesis * initialhypothesis = new TranslationHypothesis(NULL, this, NULL, NULL);
        
    stacks[0].insert(initialhypothesis);
    
    
    //for each stack
    for (int i = 0; i <= inputlength - 1; i++) {        
        if (!stacks[i].empty()) {
            if (DEBUG >= 1) {
                cerr << "\tDecoding Stack " << i << endl;
            }
            //pop from stacks[i]
            const TranslationHypothesis * hyp = stacks[i][0]
            stacks[i][0].erase(hyp);
            
            
            bool finalonly = (i == inputlength - 1) 
            unsigned int expanded = hyp->expand(finalonly); //will automatically add to appropriate stacks
            if (DEBUG >= 1) cerr << "\t\tExpanded " << expanded << " new hypotheses" << endl;
            unsigned int pruned = 0;
            for (int j = i+1; j <= inputlength; j++){ //prune further stacks (hypotheses may have been added to any of them)
                pruned += prune(j); 
            }
            if (DEBUG >= 1) cerr << "\t\tPruned " << pruned << " hypotheses" << endl;            
            if (hyp->children.empty()) delete hyp; //if the hypothesis failed to expand it will be pruned
            stacks[i].clear(); //stack is in itself no longer necessary, included pointer elements may live on though! Unnecessary hypotheses will be cleaned up automatically when higher-order hypotheses are deleted             
        }
    }   
    
    //solutions are now in last stack: stacks[inputlength]         
}


unsigned int StackDecoder::prune(int stackindex) {    
    //TODO: Add consolidation stage prior to pruning stage? (merge hypotheses with same output)    
    unsigned int pruned = 0;
    unsigned int count = 0;
    unsigned double best = 0;
    for (multiset<const TranslationHypothesis*>::const_iterator iter = stacks[stackindex].begin(); iter != stacks[stackindex].end(); iter++) {
        const TranslationHypothesis*  h = iter;
        count++;
        if (best == 0) best = h->score(); //will only be set once, first is always best
        if ((count > stacksize) || (h->score() < best * prunethreshold)) {
            pruned++;
            stacks[stackindex].erase(iter); //delete form here onwards
            delete h;            
        }
    }
    return pruned;
}


StackDecoder::~StackDecoder() {
    for (multiset<const TranslationHypothesis*>::const_iterator iter = stacks[inputlength].begin(); iter != stacks[inputlength].end(); iter++) {
        const TranslationHypothesis*  h = iter;
        delete h;
    }
}

string StackDecoder::solution(ClassDecoder & targetclassdecoder) {
    if (stacks[inputlength].empty()) {
        cerr << "ERROR: No solution found!" << endl;
        exit(6);
    }
    TranslationHypothesis * sol = stacks[inputlength].begin();
    EncData s = sol->getoutput();
    return s.decode(targetclassdecoder);
}


TranslationHypothesis::TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const vector<double> & tscores) {
    this->parent = parent;
    this->decoder = decoder;            
    if (parent != NULL) parent.children.push_back(this);        
    this->sourcegram = sourcegram;
    this->targetgram = targetgram;
    this->sourceoffset = sourceoffset;
    this->targetoffset = targetoffset; 
    
    
    if ((sourcegram != NULL) && (sourcegram.isskipgram())) {
        vector<pair<int, int> > gaps;
        *((const EncSkipGram *) sourcegram).gaps(gaps);
        for (vector<pair<int, int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) sourcegaps.push_back( make_pair<unsigned char, unsigned char>(iter->first, iter->second) );        
    }
    
    if ((targetgram != NULL) && (targetgram.isskipgram())) {
        vector<pair<int, int> > gaps;
        *((const EncSkipGram *) targetgram).gaps(gaps);
        for (vector<pair<int, int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) targetgaps.push_back( make_pair<unsigned char, unsigned char>(iter->first, iter->second) );
    }    
    
    this->tscores = vector(tscores.begin(), tscores.end());    
}

TranslationHypothesis::~TranslationHypothesis() {
    if (parent != NULL) {
        parent.children.erase(this);
        if (parent.children.empty()) delete parent;
    } 
}

bool TranslationHypothesis::final(){
        if (!targetgaps.empty()) return false;
        return (inputcoverage() = decoder.inputlength); 
}


unsigned int TranslationHypothesis::expand(bool finalonly) {
    unsigned int expanded = 0;
    //expand directly in decoder.stack()
    if (targetgaps.empty()) {
        //find new source fragment
        for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = decoder->sourcefragments.begin(); iter != decoder->sourcefragments.end(); iter++) {
            const EncAnyGram * sourcecandidate = iter->first;
            const CorpusReference ref = iter->second; 
            if (!conflicts(sourcecandidate, ref)) {
                //find target fragments for this source fragment
                for (std::unordered_map<const EncAnyGram*, vector<double> >::iterator iter2 =  phrasetable->alignmatrix[sourcecandidate].begin(); iter2 != phrasetable->alignmatrix[sourcecandidate].end(); iter2++) {
                    //create hypothesis for each target fragment
                    const EncAnyGram * targetcandidate = iter2->first;
                    newhypo = new Hypothesis(self, decoder, sourcecandidate, ref.token, targetcandidate, targetoffset + targetgram->n() , iter2->second);
                    if ((finalonly) && (!newhypo->final()) {
                        delete newhypo;
                        break;
                    } 
                    //add to proper stack
                    int cov = newhypo->inputcoverage();                    
                    decoder->stacks[cov].insert(newhypo);
                    expanded++;
                }                 
            }        
        }
        return expanded;
    } else {
        //there are target-side gaps, attempt to fill
        //TODO
        return expanded;
    }
}

bool TranslationHypothesis::conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref) {
    if ((sourcegram == NULL) && (parent == NULL)) return false; //no sourcegram, no conflict (this is an empty initial hypothesis)
    
    if (sourcecandidate->hash() == sourcegram->hash()) return true; //source was already added, can not add twice
     
    if ( (sourcecandidate->ref.token + sourcecandidate->n() > sourceoffset) && ( sourcecandidate->ref.token < sourceoffset + sourcegram->n()  ) ) { 
        //conflict    
        if (sourcegram->isskipgram()) {
            //if this falls nicely into a gap than it may not be a conflict after all
            ingap = false;
            for (vector<pair<unsigned char, unsigned char> >::iterator iter = sourcegaps.begin(); iter < sourcegaps.end(); iter++) {
                if (sourcecandidate->ref.token + sourcecandidate->n() > sourceoffset + iter->first) &&  ( sourcecandidate->ref.token < sourceoffset +iter->first + iter->second ) ) {
                    ingap = true;
                    break;       
                }
            }            
            if (!ingap) {
                return true
            }
        } else {
            //MAYBE TODO: deal with partial overlap?
            return true
        }
    }
    
    //no confict, check parents    
    if (parent != NULL) {
        return parent->uncovered(sourcecandidate);
    } else {
        return false;
    }
}


void TranslationHypothesis::computeinputcoverage(vector<bool> & container) {
    bool isskipgram = sourcegram->isskipgram();
    for (int i = sourceoffset; i < sourceoffset + sourcegram->n(); i++) {
        if (isskipgram) {
            ingap = false;
            for (vector<pair<unsigned char, unsigned char> >::iterator iter = sourcegaps.begin(); iter < sourcegaps.end(); iter++) {
                if (sourcecandidate->ref.token + sourcecandidate->n() > sourceoffset + iter->first) &&  ( sourcecandidate->ref.token < sourceoffset +iter->first + iter->second ) ) {
                    ingap = true;
                    break;       
                }             
            }
            if (!ingap) {
                container[i] = true;
            }
        } else {
            container[i] = true;
        }        
    }    
}

int TranslationHypothesis::inputcoverage() {
    const int L = decoder.input.length()
    vector<bool> container;
    for (int i = 0; i < L; i++) {
        container.push_back(false);
    }
    computesourcecoverage(container);
    int result = 0;
    for (int i = 0; i < L; i++)
        if (container[i]) result++; 
    }
    return result;
}

EncNGram * TranslationHypothesis::getoutputtoken(int index) {
    index = index - targetoffset;
    if ((index < 0) || (index >= targetgram->n()) {
        cerr << "ERROR: TranslationHypothesis::getoutputtoken() with index " << index << " is out of bounds" << endl;
        exit(6);
    }
    if (targetgram->isskipgram()) {
        const EncSkipGram * targetskipgram = (const EncSkipGram *) targetgram;
        return new targetskipgram->gettoken(index);
    } else {
        const EncNGram * targetngram = (const EncNGram *) targetgram;
        const EncNGram * unigram = targetngram->slice(index,1);
        return unigram;
    }
}

EncData TranslationHypothesis::getoutput(deque<const TranslationHypothesis*> * path) { //get output        
    //backtrack
    if (path == NULL) path = new deque<const TranslationHypothesis*>;
    if (parent != NULL) {
            path->push_front(self);
            return parent->getoutput(path);
    } else {
        //we're back at the initial hypothesis, now we construct the output forward again
        map<int, const EncNGram *> outputtokens; //unigrams, one per index                                  
        while (!path->empty()) {
            const TranslationHypothesis * hyp = path->pop_front();
            for (int i = hyp->targetoffset; i < hyp->targetoffset + hyp->targetgram->n(); i++) {
                outputtokens[i] = hyp->getoutputtoken(i);
            }                          
        }
        unsigned char buffer[8192];
        int cursor = 0;
        for (map<int, const EncNGram *>::iterator iter = outputtokens.begin(); iter != outputtokens.end(); iter++) {
            int index = iter->first;
            EncNGram * unigram = iter->second;
            for (int i = 0; i < unigram->size(); i++) {
                buffer[cursor++] = unigram->data[i];
            }
            delete unigram; //important cleanup            
        }         
        delete path; //important cleanup
        return EncData(buffer, cursor); 
    }    
}