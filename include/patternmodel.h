#include <ngram.h>
#include <cmath>
#include <cstdint>

class CorpusReference {
    /* Reference to a position in the corpus */
   public:
    uint32_t sentence;
    unsigned char token;  
    CorpusReference(uint32_t, unsigned char);  
    CorpusReference(std::istream * in);  
    void writeasbinary(std::ostream * out) const; 
    bool operator< (CorpusReference& other) {
        if (sentence < other.sentence) {
            return true;
        } else if (sentence == other.sentence) {
            return (token < other.token);
        } else {
            return false;
        }
    }
};



class NGramData {
   public:
    std::set<CorpusReference> refs;
    uint32_t count() const {
        return (uint32_t) refs.size();
    }
    void writeasbinary(std::ostream * out) const; 
    NGramData() {};
    NGramData(std::istream * in);
};

class SkipGramData {
   public:
    uint32_t _count;
    std::unordered_map<EncSkipGram,NGramData> skipcontent;
    int count() const {
       return _count;
    }
    SkipGramData() { _count = 0; }
    double entropy();
};


class EncGramIndexedModel {    
   private:
    int MINTOKENS; // = 2;
    int MINSKIPTOKENS; // = 2;
    int MINSKIPTYPES; //= 2;
    int MAXLENGTH; //= 8;
    bool DOSKIPGRAMS; //= false;
    bool DOREVERSEINDEX; //= false;
    bool DOSKIPCONTENT; //= false;
    bool DOINITIALONLYSKIP; //= true;
    bool DOFINALONLYSKIP; //= true;

    unsigned long ngramtokencount;
    unsigned long skipgramtokencount; 
    int ngramtypecount;
    int skipgramtypecount;
    
    
    int tokencount[10]; //relative token count
    int skiptokencount[10];
    int typecount[10]; //relative token count
    int skiptypecount[10];
   public:
    std::unordered_map<EncNGram,NGramData > ngrams;
    std::unordered_map<EncSkipGram,SkipGramData > skipgrams;    
    
    std::unordered_map< int,std::vector<EncNGram> > ngram_reverse_index;
    std::unordered_map< int,std::vector<EncSkipGram> > skipgram_reverse_index;
        
    EncGramIndexedModel(const std::string & filename, bool DOREVERSEINDEX = false);
    EncGramIndexedModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOREVERSEINDEX = false, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
    int maxlength() const { return MAXLENGTH; }
    
    int types() const { return ngramtypecount + skipgramtypecount; }
    int tokens() const { return ngramtokencount + skipgramtokencount; }
    
    bool exists(const EncAnyGram* key) const;
    int count(const EncAnyGram* key);
    double freq(const EncAnyGram* key);    
    double relfreq(const EncAnyGram* key);        
    //std::set<int> * index(const EncAnyGram* key);    
    //int index_size() const;
    
    /* Reverse SENTENCE index */
    std::set<int> reverse_index_keys(); 
    bool reverse_index_haskey(const int i) const;    
    int reverse_index_size(const int i);
    int reverse_index_size();
    std::vector<EncAnyGram*> reverse_index(const int i);
    EncAnyGram* get_reverse_index_item(const int, const int);
    
    
    void save(const std::string & filename);
    
    size_t hash();
    
    
    void decode(ClassDecoder & classdecoder, std::ostream *NGRAMSOUT, std::ostream *SKIPGRAMSOUT);    
};


//typedef std::unordered_map<EncNGram,set<CorpusReference> > freqlist;
//typedef std::unordered_map<EncSkipGram,set<CorpusReference> > skipgram_freqlist;

/*

class EncGramGraphModel {    
   private:    
    std::unordered_map<EncAnyGram,std::unordered_set<EncAnyGram> > rel_subsumption_children;
    //std::unordered_map<EncAnyGram,std::unordered_set<EncAnyGram> > rel_subsumption_parents; //reverse
    
    //std::unordered_map<EncNGram,std::unordered_set<EncSkipGram> > rel_inskips; //ngrams pluggable in gaps in skipgrams (reverse index of skipgram.data)
        
   public:
    EncGramGraphModel(EncGramIndexedModel& model); //compute entire model    
    EncGramGraphModel(const std::string & filename);
        

    void save(const std::string & filename);
};

*/


/*

class AlignmentModel {
   public:
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);

};

class EMAlignmentModel: public AlignmentModel {
   public:    
    EMAlignmentModel(EncGramindexedModel & sourcemodel, EncGramindexedModel & targetmodel, const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001);        
    //save(const std::string filename);
};

class CoocAlignmentModel: public AlignmentModel {
   private:
    double absthreshold;
    double relthreshold;
   public:
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
   
    CoocAlignmentModel(EncGramindexedModel & sourcemodel, EncGramindexedModel & targetmodelconst, double absthreshold = 0,  const double relthreshold = 0);         
   
    double cooc( std::set<uint32_t> & sourceindex, std::set<uint32_t> & targetindex); 
    int compute(const EncAnyGram * sourcegram, std::set<int> & sourceindex, EncGramindexedModel & targetmodel);
    
    
    
    //void save(const std::string filename);
};

*/
