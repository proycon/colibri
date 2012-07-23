#include <alignmodel.h>
#include <bitset>

class TranslationHypothesis {
    public:
        StackDecoder * decoder;
        
        TranslationHypothesis * parent;
        vector<TranslationHypothesis *> children; //reference counting
         
        //std::bitset * sourcecoverage; 
        
        const EncAnyGram * sourcegram;
        const EncAnyGram * targetgram;
        unsigned char sourceoffset;
        unsigned char targetoffset;
                
        vector<pair<unsigned char, unsigned char> > sourcegaps;
        vector<pair<unsigned char, unsigned char> > targetgaps; //filling gaps will take priority upon expansion
        vector<double> tscores; //scores from the translation table 
        
        TranslationHypothesis(TranslationHypothesis * parent, StackDecoder * decoder,  const EncAnyGram * sourcegram , unsigned char sourceoffset,  const EncAnyGram * targetgram, unsigned char targetoffset, const std::vector<double> & tscores);
        ~TranslationHypothesis();
        
        unsigned int expand(bool finalonly=false); //expands directly in the appropriate stack of the decoder. If finalonly is set, new hypotheses are expected to be final/complete, without gaps
        
        double score();
                  
        bool initial() const { return (parent == NULL) }
        bool final();        
        
        bool conflicts(const EncAnyGram * sourcecandidate, const CorpusReference & ref); //checks if a source pattern's coverage conflicts with this hypothesis, if not, it is a candidate to be added upon hypothesis expansion
        
                
        void computeinputcoverage(vector<bool> & container); //compute source coverage        
        int inputcoverage(); //return input coverage (absolute number of tokens covered)
        
        bool operator< (const TranslationHypothesis& other) const {
            return (score() < other.score());
        }
        
        
        EncNGram getoutputtoken(int index); //get output token (unigram, will return unknown class if in gap)                    
        EncData getoutput(deque<const TranslationHypothesis*> * path = NULL); //get output                
}



class StackDecoder {
    private:
        int stacksize;
        int DEBUG;
    public:
        EncData input;
        unsigned int inputlength;
        TranslationTable * translationtable;
                
        vector<pair<const EncAnyGram*, CorpusReference> >  sourcefragments;        
        
        map<unsigned char, multiset<const TranslationHypothesis *> > stacks;
                
        StackDecoder(const EncData & input, TranslationTable * translationtable, int stacksize, double prunethreshold, int maxn);
        ~StackDecoder();
        
        void decode();
                
        unsigned int getsourcefragments(int maxn);
        
        unsigned int prune(int stackindex);
        
        string StackDecoder::solution(ClassDecoder & targetclassdecoder); //return solution as string        
}