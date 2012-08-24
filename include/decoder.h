#ifndef DECODER_H
#define DECODER_H


#include <alignmodel.h>
#include <lm.h>
#include <math.h>
#include <deque>

class StackDecoder;
class Stack;


class TranslationHypothesis {
    private:        
        
        bool deleted; //for debug only
    public:
        bool keep; //if true, prevents this hypothesis from being automatically deleted by its last dying child
        std::vector<double> tscores;
        double tscore;
        double lmscore;
        double dscore;
        double futurecost;
        StackDecoder * decoder;
        
        TranslationHypothesis * parent;
        std::vector<TranslationHypothesis *> children; //reference counting
         
        std::vector<bool> inputcoveragemask; 
    
        const EncAnyGram * sourcegram;
        const EncAnyGram * targetgram;
        unsigned char sourceoffset;
        unsigned char targetoffset;
        
        const EncNGram * history; //order-1 n-gram in targetlanguage, used for computation of LM score 
                
        std::vector<std::pair<unsigned char, unsigned char> > sourcegaps; //offsets are relative to sourcegram
        std::vector<std::pair<unsigned char, unsigned char> > targetgaps; //filling gaps will take priority upon expansion, offsets are relative to the entire hypothesis  
        
        TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const std::vector<double> & tscores);
        ~TranslationHypothesis();
        void cleanup();
        
        unsigned int expand(); //expands directly in the appropriate stack of the decoder. If finalonly is set, new hypotheses are expected to be final/complete, without gaps
        
        double basescore() const;    
        double score() const;        
        void report();  
        bool initial() const { return (parent == NULL); }
        
        bool conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref, bool skipduplicatecheck = false); //checks if a source pattern's coverage conflicts with this hypothesis, if not, it is a candidate to be added upon hypothesis expansion. Latter parameter is for internal usage
        bool final();     
        bool hasgaps() const;   
        
        void stats();
        
        bool fertile();
                
        int fitsgap(const EncAnyGram *, const int offset = 0); //returns -1 if no fit, index of gap begin otherwise
                
        //void computeinputcoverage(vector<bool> & container); //compute source coverage        
        int inputcoverage(); //return input coverage (absolute number of tokens covered)
        
        bool deletable(); 
        
        EncNGram * getoutputtoken(int index); //get output token (unigram, will return unknown class if in gap)                    
        EncData getoutput(std::deque<TranslationHypothesis*> * path = NULL); //get output
                        
};


struct HypCompare : public std::binary_function<const TranslationHypothesis*, const TranslationHypothesis*, bool>
{
    bool operator()(const TranslationHypothesis* x, TranslationHypothesis* y) const
    {   
        return x->score() >= y->score();
    }
};

class DecodeStats {
  public:
    //statistics
    unsigned int discarded;
    unsigned int expanded; 
    unsigned int pruned;
    unsigned int gapresolutions;
    
    std::map<int,int> stacksizes;
    std::map<int,int> gappystacksizes;


    std::map<int,int> sourcengramusage;
    std::map<int,int> sourceskipgramusage;
    std::map<int,int> targetngramusage;
    std::map<int,int> targetskipgramusage;
    std::vector<int> steps;
     
    DecodeStats() { reset(); } 
    void reset() {

        discarded = 0;
        expanded = 0; 
        pruned = 0;
        gapresolutions = 0;
    
        sourcengramusage.clear();
        sourceskipgramusage.clear();
        targetngramusage.clear();
        targetskipgramusage.clear();
        stacksizes.clear();
        gappystacksizes.clear();
        steps.clear();
    }

    void output();
    
    void add(DecodeStats& other) {
        for (std::map<int,int>::iterator iter = other.sourcengramusage.begin(); iter != other.sourcengramusage.end(); iter++) sourcengramusage[iter->first] += iter->second; 
        for (std::map<int,int>::iterator iter = other.sourceskipgramusage.begin(); iter != other.sourceskipgramusage.end(); iter++) sourceskipgramusage[iter->first] += iter->second;
        for (std::map<int,int>::iterator iter = other.targetngramusage.begin(); iter != other.targetngramusage.end(); iter++) targetngramusage[iter->first] += iter->second; 
        for (std::map<int,int>::iterator iter = other.targetskipgramusage.begin(); iter != other.targetskipgramusage.end(); iter++) targetskipgramusage[iter->first] += iter->second; 
         
        discarded += other.discarded;
        expanded += other.expanded;
        pruned += other.pruned;
        gapresolutions += other.gapresolutions;
    } 

};

class StackDecoder;

class Stack {
   private:
    int stacksize;
    double prunethreshold;
    StackDecoder * decoder;
   public:
   int index;
    Stack(StackDecoder * decoder, int index, int stacksize, double prunethreshold);
    ~Stack();
    Stack(const Stack& ref); //limited copy constructor
    std::list<TranslationHypothesis *> contents;
    bool add(TranslationHypothesis *); 
    double bestscore();
    double worstscore();
    size_t size() { return contents.size(); }
    bool empty() { return contents.empty(); }
    void clear();   
    int prune();
    TranslationHypothesis * pop();
};

class StackDecoder {
    private:
        int stacksize;
        double prunethreshold;
    public:
        int DEBUG;
        std::map<std::pair<int, int>, double> futurecost; //(start, end) => cost    
        EncData input;
        unsigned int inputlength;
        TranslationTable * translationtable;
        ClassDecoder * sourceclassdecoder;
        ClassDecoder * targetclassdecoder;
        LanguageModel * lm;    
        bool globalstats;
        
        
        
        std::vector<double> tweights; //translation model weights
        double dweight; //distortion model weight
        double lweight; //language model weight
        int dlimit; //distortion limit
        
                
        std::vector<std::pair<const EncAnyGram*, CorpusReference> >  sourcefragments;        
        
        std::vector<Stack> stacks;
        std::vector<Stack> gappystacks;
                
        StackDecoder(const EncData & input, TranslationTable * translationtable, LanguageModel * lm, int stacksize, double prunethreshold, std::vector<double> tweights, double dweight, double lweight, int dlimit, int maxn, int debug, ClassDecoder *, ClassDecoder *, bool globalstats = false);
        ~StackDecoder();
        
        TranslationHypothesis * decodestack(Stack & stack); //returns fallback hypothesis if dead, NULL otherwise
        TranslationHypothesis * decode(); //returns solution
                
        unsigned int getsourcefragments(int maxn);
        
        unsigned int prune(int stackindex);    
        
        void computefuturecost();      

        DecodeStats stats;


};


#endif
