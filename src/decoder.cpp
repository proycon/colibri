#include <decoder.h>
#include <algorithm> 
#include <getopt.h>
#include <iostream>
#include <sstream>


using namespace std;

unsigned char UNKNOWNCLASS = 2;
unsigned char BOSCLASS = 3;
unsigned char EOSCLASS = 4;

StackDecoder::StackDecoder(const EncData & input, TranslationTable * translationtable, LanguageModel * lm, int stacksize, double prunethreshold, vector<double> tweights, double dweight, double lweight, int maxn) {
        this->input = input;
        this->inputlength = input.length();
        this->translationtable = translationtable;
        this->lm = lm;
        this->stacksize = stacksize;
        this->prunethreshold = prunethreshold;
        sourcefragments = translationtable->getpatterns(input.data,input.size(), true, 0,1,maxn);
        //sanity check:
        /*for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {
            const EncAnyGram * anygram = iter->first;
            cerr << anygram << endl;            
        }*/
        
        this->tweights = vector<double>(tweights.begin(), tweights.end());
        this->dweight = dweight;
        this->lweight = lweight;
        
        computefuturecost();
        //TODO: Deal with unknown tokens? and <s> </s>
        
        
}

void StackDecoder::setdebug(int debug) {
    this->DEBUG = debug;
}

void StackDecoder::computefuturecost() {
        map<pair<int,int>, double> sourcefragments_costbyspan;
        //reorder source fragments by span for more efficiency
        for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = sourcefragments.begin(); iter != sourcefragments.end(); iter++) {
            const EncAnyGram * anygram = iter->first;
            const CorpusReference ref = iter->second;     
            const int n = anygram->n();    
            const pair<int,int> span = make_pair<int,int>((int) ref.token, (int) n);
             
            const EncAnyGram * candidate = translationtable->getsourcekey(anygram);
            
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
            for (int start = 0; start < inputlength - length; start++) {
                const pair<int,int> span = make_pair((int) start,(int) length);
                bool found = false;
                for (map<pair<int,int>, double>::iterator iter = sourcefragments_costbyspan.find(span); iter != sourcefragments_costbyspan.end(); iter++) {
                    found = true;
                    futurecost[span] = sourcefragments_costbyspan[span];
                }
                if (!found) futurecost[span] = -INFINITY;
                for (int i = 1; i < length; i++) {
                    double spanscore = futurecost[make_pair((int) start,(int) i)] + futurecost[make_pair((int) start+i,(int) length - i)];
                    if (spanscore > futurecost[span]) { //(higher score -> lower cost)
                        futurecost[span] = spanscore;
                    }
                }
            }
        }
}
    
void StackDecoder::decode() {
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
        
    stacks[0].insert(initialhypothesis);
    
    
    //for each stack
    for (int i = 0; i <= inputlength - 1; i++) {        
        if (!stacks[i].empty()) {
            if (DEBUG >= 1) {
                cerr << "\tDecoding Stack " << i << endl;
            }
            //pop from stacks[i]
            TranslationHypothesis * hyp = *(stacks[i].begin());
            stacks[i].erase(hyp);
            
            
            bool finalonly = (i == inputlength - 1); 
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
    double best = 0;
    for (multiset<TranslationHypothesis*>::const_iterator iter = stacks[stackindex].begin(); iter != stacks[stackindex].end(); iter++) {
        const TranslationHypothesis*  h = *iter;
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
    for (multiset<TranslationHypothesis*>::const_iterator iter = stacks[inputlength].begin(); iter != stacks[inputlength].end(); iter++) {
        const TranslationHypothesis*  h = *iter;
        delete h;
    }
}

string StackDecoder::solution(ClassDecoder & targetclassdecoder) {
    if (stacks[inputlength].empty()) {
        cerr << "ERROR: No solution found!" << endl;
        exit(6);
    }
    TranslationHypothesis * sol = *(stacks[inputlength].begin());
    EncData s = sol->getoutput();
    return s.decode(targetclassdecoder);
}


TranslationHypothesis::TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const vector<double> & tscores) {
    this->parent = parent;
    this->decoder = decoder;            
    if (parent != NULL) parent->children.push_back(this);        
    this->sourcegram = sourcegram;
    this->targetgram = targetgram;
    this->sourceoffset = sourceoffset;
    this->targetoffset = targetoffset; 
    
    
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
    double tscore = 0; 
    for (int i = 0; i < decoder->tweights.size(); i++) {
        double p = tscores[i];
        if (p > 0) p = log10(p); //turn into logprob, base10 
        tscore += decoder->tweights[i] * p;
    }
    
    
    double lmscore = 0;
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
    
    //TODO: deal differently with filling in gaps in skips?
    int prevpos = 0;
    if ((parent != NULL) && (!parent->initial())) {
        prevpos = parent->sourceoffset + parent->sourcegram->n();        
    } 
    double dscore = pow(decoder->dweight, abs(sourceoffset - prevpos - 1));  
        
    //Compute estimate for future score
    double futurescore = 0;    
    begin = -1;    
    for (int i = 0; i < inputcoveragemask.size(); i++) {
        if (!inputcoveragemask[i]) {
            begin = i;
        } else if (begin != -1) {
            futurescore += decoder->futurecost[make_pair((int) begin,(int) i-begin)];
            begin = -1;
        }
    }
    if (begin != -1) futurescore += decoder->futurecost[make_pair((int) begin,(int) inputcoveragemask.size()-begin)];            
    
    
    //total score            
    _score = tscore + lmscore + dscore + futurescore;
  
}

TranslationHypothesis::~TranslationHypothesis() {
    if (history != NULL) {
        delete history;
    }
    if (parent != NULL) {         
        vector<TranslationHypothesis *>::iterator iter = find(parent->children.begin(), parent->children.end(), this); 
        parent->children.erase(iter);
        if (parent->children.empty()) delete parent;
    } 
}

bool TranslationHypothesis::final(){
        if (!targetgaps.empty()) return false;
        return (inputcoverage() == decoder->inputlength); 
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
                for (std::unordered_map<const EncAnyGram*, vector<double> >::iterator iter2 =  decoder->translationtable->alignmatrix[sourcecandidate].begin(); iter2 != decoder->translationtable->alignmatrix[sourcecandidate].end(); iter2++) {
                    //create hypothesis for each target fragment
                    const EncAnyGram * targetcandidate = iter2->first;
                    TranslationHypothesis * newhypo = new TranslationHypothesis(this, decoder, sourcecandidate, ref.token, targetcandidate, targetoffset + targetgram->n() , iter2->second);
                    if ((finalonly) && (!newhypo->final())) {
                        delete newhypo;
                        break;
                    } 
                    //add to proper stack
                    int cov = newhypo->inputcoverage();
                    decoder->stacks[cov].insert(newhypo);
                    expanded++;
                    //no delete newhypo
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
    /*
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
    */
}

EncNGram * TranslationHypothesis::getoutputtoken(int index) {
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
            path->pop_front(); //WARNING: Calls TranslationHypothesis * destructor??? may have unwanted side-effect
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
            delete unigram; //important cleanup            
        }         
        delete path; //important cleanup
        return EncData(buffer, cursor); 
    }    
}


double TranslationHypothesis::score() const {
    return _score;
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
    
}

void addsentencemarkers(ClassDecoder & targetclassdecoder, ClassEncoder & targetclassencoder) {
    targetclassdecoder.add(BOSCLASS,"<s>");
    targetclassdecoder.add(EOSCLASS,"</s>");
    targetclassencoder.add("<s>", BOSCLASS);
    targetclassencoder.add("</s>", EOSCLASS);
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
    
    char c;    
    while ((c = getopt_long(argc, argv, "ht:S:T:s:p:l:W:L:D:",long_options,&option_index)) != -1) {
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
        
    string input;
    unsigned char buffer[8192]; 
    int size;
    while (getline(cin, input)) {        
        cerr << "INPUT: " << input << endl;        
        size = sourceclassencoder.encodestring(input, buffer, true);
        
        
                
        StackDecoder decoder = StackDecoder(EncData(buffer, size), &transtable, &lm, stacksize, prunethreshold, tweights, dweight, lweight, maxn);
        decoder.decode();
        cout << decoder.solution(targetclassdecoder) << endl;        
    }
}
