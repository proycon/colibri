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
        bool keep; //if true, prevents this hypothesis from being automatically deleted by its last dying child
    public:
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
                
        std::vector<std::pair<unsigned char, unsigned char> > sourcegaps;
        std::vector<std::pair<unsigned char, unsigned char> > targetgaps; //filling gaps will take priority upon expansion 
        
        TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const std::vector<double> & tscores);
        ~TranslationHypothesis();
        void cleanup();
        
        unsigned int expand(bool finalonly=false); //expands directly in the appropriate stack of the decoder. If finalonly is set, new hypotheses are expected to be final/complete, without gaps
        
        double score() const;        
        void report();  
        bool initial() const { return (parent == NULL); }
        
        bool conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref); //checks if a source pattern's coverage conflicts with this hypothesis, if not, it is a candidate to be added upon hypothesis expansion
        bool final();        
        
                
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

class StackDecoder;

class Stack {
   private:
    int index;
    int stacksize;
    double prunethreshold;
   public:
    Stack(int index, int stacksize, double prunethreshold);
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
        
        std::vector<double> tweights; //translation model weights
        double dweight; //distortion model weight
        double lweight; //language model weight
        
                
        std::vector<std::pair<const EncAnyGram*, CorpusReference> >  sourcefragments;        
        
        std::vector<Stack> stacks;
                
        StackDecoder(const EncData & input, TranslationTable * translationtable, LanguageModel * lm, int stacksize, double prunethreshold, std::vector<double> tweights, double dweight, double lweight, int maxn);
        ~StackDecoder();
        
        void decode();
                
        unsigned int getsourcefragments(int maxn);
        
        unsigned int prune(int stackindex);
        
        std::string solution(ClassDecoder & targetclassdecoder); //return solution as string
        
        void computefuturecost();
        void setdebug(int debug, ClassDecoder *, ClassDecoder *);      

};


#endif
