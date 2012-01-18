#include <ngram.h>
#include <cmath>
#include <cstdint>

const char MAXN = 20;

class CorpusReference {
    /* Reference to a position in the corpus */
   public:
    uint32_t sentence;
    unsigned char token;  
    CorpusReference(uint32_t, unsigned char);  
    CorpusReference(std::istream * in);  
    CorpusReference(const CorpusReference& other) { //copy constructor
        sentence = other.sentence;
        token = other.token;
    };     
    void writeasbinary(std::ostream * out) const; 
    bool operator< (const CorpusReference& other) const {
        if (sentence < other.sentence) {
            return true;
        } else if (sentence == other.sentence) {
            return (token < other.token);
        } else {
            return false;
        }
    }
    bool operator==(const CorpusReference &other) const { return ( (sentence == other.sentence) && (token == other.token)); };
    bool operator!=(const CorpusReference &other) const { return ( (sentence != other.sentence) || (token != other.token)); };
};

class AnyGramData {
    public: 
     virtual std::set<CorpusReference> get_refs() const =0;
};

class NGramData: public AnyGramData {
   public:
    std::set<CorpusReference> refs;
    uint32_t count() const {
        return (uint32_t) refs.size();
    }
    void writeasbinary(std::ostream * out) const; 
    std::set<CorpusReference> get_refs() const {
        return refs;
    }
    NGramData() {};
    NGramData(std::istream * in);
};

class SkipGramData: public AnyGramData {
   public:
    uint32_t _count;
    std::unordered_map<EncSkipGram,NGramData> skipcontent;
    int count() const {
       return _count;
    }
    SkipGramData() { _count = 0; }
    double entropy() const;
    std::set<CorpusReference> get_refs() const;
};

class ModelReader {
   public:
    virtual uint64_t id() =0;
    virtual void readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes, bool ignore=false) =0;
    virtual void readngram(std::istream * in, const EncNGram & ngram, bool ignore = false) =0;
    virtual void readskipgram(std::istream * in, const EncSkipGram & skipgram, bool ignore = false) =0;
    virtual void readfooter(std::istream * in, bool ignore = false) =0;    
    
    virtual void readfile(const std::string & filename);
};

class ModelWriter {
   public:
    virtual uint64_t id() =0;    
    virtual void writeheader(std::ostream * out) =0;
    virtual void writengrams(std::ostream * out) =0;
    virtual void writengram(std::ostream * out, const EncNGram & ngram) =0;
    virtual void writeskipgrams(std::ostream * out) =0;
    virtual void writeskipgram(std::ostream * out, const EncSkipGram & skipgram) =0;
    virtual void writefooter(std::ostream * out) =0;
    
    virtual uint64_t tokens() const =0;
    virtual uint64_t types() const =0;
    virtual void writefile(const std::string & filename);
};




class IndexedPatternModel: public ModelReader, public ModelWriter {    
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

    
    
    int ngramtypecount;
    int skipgramtypecount;    

   public:

    unsigned long ngramtokencount;
    unsigned long skipgramtokencount; 
    
    int tokencount[MAXN]; //relative token count
    int skiptokencount[MAXN];
    int typecount[MAXN]; //relative token count
    int skiptypecount[MAXN];   
   
    std::unordered_map<const EncNGram,NGramData > ngrams;
    std::unordered_map<const EncSkipGram,SkipGramData > skipgrams;    
    
    std::unordered_map< int,std::vector<EncNGram> > ngram_reverse_index;
    std::unordered_map< int,std::vector<EncSkipGram> > skipgram_reverse_index;
           
    IndexedPatternModel(const std::string & filename = "", bool DOREVERSEINDEX = false);
    IndexedPatternModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOREVERSEINDEX = false, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
    int maxlength() const { return MAXLENGTH; }
    
    uint64_t types() const { return ngrams.size() + skipgrams.size(); }
    uint64_t tokens() const { return ngramtokencount + skipgramtokencount; }
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);
    const AnyGramData* getdata(const EncAnyGram* key);
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
    
    virtual uint64_t id() { return 1; }
    virtual void readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes, bool ignore = false) {};
    virtual void readngram(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgram(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out) {};
    virtual void writengrams(std::ostream * out);
    virtual void writengram(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writeskipgram(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writefooter(std::ostream * out) {}; 
        
    void save(const std::string & filename) { ModelWriter::writefile(filename); }
    
    size_t hash();
    
    
    void decode(ClassDecoder & classdecoder, std::ostream *NGRAMSOUT, std::ostream *SKIPGRAMSOUT);    
};


class UnindexedPatternModel: public ModelReader, public ModelWriter {
   /* unindexed model */    
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

    
    int ngramtypecount;
    int skipgramtypecount;
    

   public:
    unsigned long ngramtokencount;
    unsigned long skipgramtokencount; 
    
    int tokencount[MAXN]; //relative token count
    int skiptokencount[MAXN];
    int typecount[MAXN]; //relative token count
    int skiptypecount[MAXN];   
   
    std::unordered_map<const EncNGram,uint32_t > ngrams;
    std::unordered_map<const EncSkipGram,uint32_t > skipgrams;    
            
    UnindexedPatternModel(const std::string & filename, bool DOREVERSEINDEX = false);
    UnindexedPatternModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
    int maxlength() const { return MAXLENGTH; }
    
    uint64_t types() const { return ngrams.size() + skipgrams.size(); }
    uint64_t tokens() const { return ngramtokencount + skipgramtokencount; }
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);
    int count(const EncAnyGram* key);
    double freq(const EncAnyGram* key);    
    double relfreq(const EncAnyGram* key);        
    //std::set<int> * index(const EncAnyGram* key);    
    //int index_size() const;

    
    virtual uint64_t id() { return 0; }
    virtual void readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes, bool ignore = false) {};
    virtual void readngram(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgram(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out) {};
    virtual void writengrams(std::ostream * out);
    virtual void writengram(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writeskipgram(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writefooter(std::ostream * out) {}; 
        
    void save(const std::string & filename) { ModelWriter::writefile(filename); }
    
    size_t hash();
    
    void decode(ClassDecoder & classdecoder, std::ostream *NGRAMSOUT, std::ostream *SKIPGRAMSOUT);    
};


//typedef std::unordered_map<EncNGram,set<CorpusReference> > freqlist;
//typedef std::unordered_map<EncSkipGram,set<CorpusReference> > skipgram_freqlist;



class GraphPatternModel: public ModelReader, public ModelWriter {               
   protected:
    bool DOPARENTS;
    bool DOCHILDREN;
    bool DOXCOUNT;
    bool DOSKIPCONTENT;
    bool DOINSKIPCONTENT;
    
    bool HASPARENTS;
    bool HASCHILDREN;
    bool HASXCOUNT;
    bool HASSKIPCONTENT;
    bool HASINSKIPCONTENT;
    
    bool DELETEMODEL;
    
    void readrelations(std::istream * in,const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > &, bool ignore = false);
    void writerelations(std::ostream * out, const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & );
    
    bool secondpass;       
   public:
   
    IndexedPatternModel * model;
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_parents;
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_children;        
    std::unordered_map<const EncAnyGram*,int> data_xcount;        
   
    //std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipcontent; //skipgram -> skipcontent       
    //std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_inskipcontent; //skipcontent -> skipgram
    
    uint64_t types() const { return model->types(); }
    uint64_t tokens() const { return model->tokens(); }
    

   
    GraphPatternModel(IndexedPatternModel * model, bool DOPARENTS=true,bool DOCHILDREN=false,bool DOXCOUNT=false); //compute entire model
    GraphPatternModel(const std::string & graphmodelfilename, IndexedPatternModel * model, bool DOPARENTS=true,bool DOCHILDREN=true,bool DOXCOUNT=true ) {
        //do everything (provided that it exists in file)
        this->DOPARENTS = DOPARENTS;
        this->DOCHILDREN = DOCHILDREN;
        this->DOXCOUNT = DOXCOUNT;
        
    	DELETEMODEL = false;        
        this->model = model;
        secondpass = true;
    	readfile(graphmodelfilename);        
    }
    GraphPatternModel(const std::string & graphmodelfilename, bool DOPARENTS=true,bool DOCHILDREN=true,bool DOXCOUNT=true) {
        //do everything (provided that it exists in file)
        this->DOPARENTS = DOPARENTS;
        this->DOCHILDREN = DOCHILDREN;
        this->DOXCOUNT = DOXCOUNT;    
    
    	DELETEMODEL = true;
    	model = new IndexedPatternModel();
    	std::cerr << "Pass one, reading implied indexedpatternmodel..." << std::endl;
    	//reading is done in two passes
    	secondpass = false;
    	readfile(graphmodelfilename);    
        std::cerr << "Pass two, reading graph data..." << std::endl;
        secondpass = true;
        readfile(graphmodelfilename);
    }    
    
    ~GraphPatternModel();
    
    int xcount(const EncAnyGram* anygram); //exclusive count    
        
    void save(const std::string & filename) { writefile(filename); }
    
    virtual uint64_t id() { return 20; }
    
    virtual void readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes, bool ignore = false);
    virtual void readngram(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgram(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out);
    virtual void writengrams(std::ostream * out);
    virtual void writengram(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgram(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writefooter(std::ostream * out) {};    
    
    void decode(ClassDecoder & classdecoder, std::ostream *NGRAMSOUT, std::ostream *SKIPGRAMSOUT);    
};


class IndexCountData {
	public:
	 uint32_t xcount;
	 std::multiset<uint32_t> sentences; //may occur multiple times in same sentence
};

class DoubleIndexedGraphPatternModel: public ModelReader {
    // Read only model, reads graphpatternmodel in a simplified, less memory intensive representation optimised for alignment tasks 
    private:
     bool HASPARENTS;
     bool HASCHILDREN;
     bool HASXCOUNT;
     
    
     int ngramtypecount;
     int skipgramtypecount;    
     void readrelations(std::istream * in,const EncAnyGram*);
    public:

     unsigned long ngramtokencount;
     unsigned long skipgramtokencount;  
     std::unordered_map<const EncNGram, IndexCountData> ngrams;
     std::unordered_map<const EncSkipGram,IndexCountData> skipgrams;
    
     std::unordered_map<uint32_t,std::vector<const EncAnyGram*> > reverseindex;    
     DoubleIndexedGraphPatternModel(const std::string & filename); //read a graph pattern model
  
     uint64_t types() const { return ngrams.size() + skipgrams.size(); }
     uint64_t tokens() const { return ngramtokencount + skipgramtokencount; }
     
    virtual uint64_t id() { return 20; }
    
    virtual void readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes, bool ignore = false);
    virtual void readngram(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgram(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};
};
