#ifndef PATTERNMODEL_H
#define PATTERNMODEL_H

#include <ngram.h>
#include "classencoder.h"
#include <cmath>
#include <cstdint>
#include <map>

const char MAXN = 20;

const short GRAPHPATTERNMODELVERSION = 1;
const short INDEXEDPATTERNMODELVERSION = 2;

enum ModelType {
	UNINDEXEDPATTERNMODEL = 10, 
    INDEXEDPATTERNMODEL = 20,
    GRAPHPATTERNMODEL = 30,        
};


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
   protected:
    uint64_t totaltokens;
    uint64_t totaltypes;
    bool DEBUG; 
    
    int FOUNDMINN; 
    int FOUNDMAXN; 
   public:
    uint64_t model_id;
    virtual uint64_t id() =0;
    virtual void readheader(std::istream * in, bool ignore=false) =0;
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, bool ignore = false) =0;
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore = false) =0;
    virtual void readfooter(std::istream * in, bool ignore = false) =0;    
    
    virtual void readfile(const std::string & filename, const bool DEBUG=false);
    
    int getminn() { return FOUNDMINN; };
    int getmaxn() { return FOUNDMAXN; };
};

class ModelWriter {

   public:
    virtual uint64_t id() =0;    
    virtual void writeheader(std::ostream * out) =0;
    virtual void writengrams(std::ostream * out) =0;
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram) =0;
    virtual void writeskipgrams(std::ostream * out) =0;
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram) =0;
    virtual void writefooter(std::ostream * out) =0;
    

    virtual uint64_t types() const =0;
    virtual uint64_t tokens() const =0;
    virtual void writefile(const std::string & filename);
};

class ModelQuerier {
     std::map<int, std::vector< std::vector< std::pair<int,int> > > > gapconf;
	public:
	 ModelQuerier();
	 virtual const EncAnyGram * getkey(const EncAnyGram * anygram) =0;
	 virtual int maxlength() const =0;
	 virtual int count(const EncAnyGram * anygram) =0;
	 virtual double freq(const EncAnyGram * anygram) =0;
	 virtual void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &) =0;	 	 
	 std::vector<std::pair<const EncAnyGram*, CorpusReference> > getpatterns(const unsigned char * data, const unsigned char datasize, bool doskipgrams=true, uint32_t linenum=0, const int minn = 1, const int maxn = MAXN); 
	 void querier(ClassEncoder & encoder, ClassDecoder & decoder, bool exact = false,bool repeat = true, const int minn = 1, const int maxn = MAXN);	 	  
};


class IndexedPatternModel: public ModelReader, public ModelWriter, public ModelQuerier {    
   private:
    int MINTOKENS; // = 2;
    int MINSKIPTOKENS; // = 2;
    int MINSKIPTYPES; //= 2;
    int MAXLENGTH; //= 8;
    bool DOSKIPGRAMS; //= false;
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
           
    IndexedPatternModel(const std::string & filename = "", const bool DEBUG=false);
    IndexedPatternModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    IndexedPatternModel(const std::string & corpusfile, IndexedPatternModel & refmodel, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
    int maxlength() const { return MAXLENGTH; }
    
    std::vector<unsigned char> sentencesize;
    
    uint64_t types() const { return ngrams.size() + skipgrams.size(); }    
    uint64_t tokens() const { return totaltokens; }    
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);
    const AnyGramData* getdata(const EncAnyGram* key);
    int count(const EncAnyGram* key);
    double freq(const EncAnyGram* key);    
    //double relfreq(const EncAnyGram* key);        
    //std::set<int> * index(const EncAnyGram* key);    
    //int index_size() const;
    
    void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &);	 	 

    
    /* Reverse SENTENCE index */
    std::set<int> reverse_index_keys(); 
    bool reverse_index_haskey(const int i) const;    
    int reverse_index_size(const int i);
    int reverse_index_size();
    std::vector<EncAnyGram*> reverse_index(const int i);
    EncAnyGram* get_reverse_index_item(const int, const int);
    
    virtual uint64_t id() { return INDEXEDPATTERNMODEL + INDEXEDPATTERNMODELVERSION; }
    virtual void readheader(std::istream * in, bool ignore = false) {};
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false);    
    
    virtual void writeheader(std::ostream * out) {};
    virtual void writengrams(std::ostream * out);
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writefooter(std::ostream * out); 
        
    void writeanygram(const EncAnyGram * anygram, std::ostream * out); //write the anygram itself (not its data!)
        
    void save(const std::string & filename) { ModelWriter::writefile(filename); }
    
    size_t hash();
    
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT);
    void decode(IndexedPatternModel & testmodel, ClassDecoder & classdecoder, std::ostream *OUT);
    
    void coveragereport(std::ostream *OUT, int segmentsize = 100000);    
};


class UnindexedPatternModel: public ModelReader, public ModelWriter, public ModelQuerier {
   /* unindexed model */    
   private:
    int MINTOKENS; // = 2;
    int MINSKIPTOKENS; // = 2;
    int MINSKIPTYPES; //= 2;
    int MAXLENGTH; //= 8;
    bool DOSKIPGRAMS; //= false;
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
            
    UnindexedPatternModel(const std::string & filename, const bool DEBUG=false);
    UnindexedPatternModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    UnindexedPatternModel(const std::string & corpusfile, UnindexedPatternModel & refmodel, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
        
    int maxlength() const { return MAXLENGTH; }
    
    uint64_t types() const { return ngrams.size() + skipgrams.size(); }
    uint64_t tokens() const { return totaltokens; }
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);
    int count(const EncAnyGram* key);
    double freq(const EncAnyGram* key);    
    //double relfreq(const EncAnyGram* key);        
    //std::set<int> * index(const EncAnyGram* key);    
    //int index_size() const;

	void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &);
    
    virtual uint64_t id() { return UNINDEXEDPATTERNMODEL; }
    virtual void readheader(std::istream * in,  bool ignore = false) {};
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out) {};
    virtual void writengrams(std::ostream * out);
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writefooter(std::ostream * out) {}; 
        
    void save(const std::string & filename) { ModelWriter::writefile(filename); }
    
    size_t hash();
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT);
    void decode(UnindexedPatternModel & testmodel, ClassDecoder & classdecoder, std::ostream *OUT);    
};


//typedef std::unordered_map<EncNGram,set<CorpusReference> > freqlist;
//typedef std::unordered_map<EncSkipGram,set<CorpusReference> > skipgram_freqlist;


int transitivereduction(std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relations );


class GraphPatternModel: public ModelReader, public ModelWriter {               
   protected:
    bool DOPARENTS;
    bool DOCHILDREN;
    bool DOXCOUNT;
    bool DOTEMPLATES;
    bool DOINSTANCES;
    bool DOSKIPUSAGE;
    bool DOSKIPCONTENT;
    bool DOSUCCESSORS;
    bool DOPREDECESSORS;
    
    
    bool HASPARENTS;
    bool HASCHILDREN;
    bool HASXCOUNT;
    bool HASTEMPLATES;
    bool HASINSTANCES;
    bool HASSKIPUSAGE;
    bool HASSKIPCONTENT;
    bool HASSUCCESSORS;
    bool HASPREDECESSORS;
    
    bool DELETEMODEL;
    bool TRANSITIVE;
    
    void readrelations(std::istream * in,const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > &, bool ignore = false);
    void writerelations(std::ostream * out, const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & );
    
    bool secondpass;       
   public:
   
    IndexedPatternModel * model;
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_parents;
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_children;        
    std::unordered_map<const EncAnyGram*,int> data_xcount;        
   
   
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_templates; //instance -> skipgram
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_instances; //skipgram -> instance
    
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipusage; //skipcontent -> skipgram
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipcontent; //skipgram -> skipcontent       
    
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_successors;  
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_predecessors;
    //std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipcontent; //skipgram -> skipcontent       
    //std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_inskipcontent; //skipcontent -> skipgram
    
    uint64_t types() const { return model->types(); }
    uint64_t tokens() const { return model->tokens(); }
    

   
    GraphPatternModel(IndexedPatternModel * model, bool DOPARENTS=true,bool DOCHILDREN=false,bool DOXCOUNT=false, bool DOTEMPLATES =false, bool DOINSTANCES = false,bool DOSKIPUSAGE=false,bool DOSKIPCONTENT=false,bool DOSUCCESSORS=false,bool DOPREDECESSORS=false); //compute entire model
    GraphPatternModel(const std::string & graphmodelfilename, IndexedPatternModel * model, bool DOPARENTS=true,bool DOCHILDREN=true,bool DOXCOUNT=true, bool DOTEMPLATES =false, bool DOINSTANCES = false,bool DOSKIPUSAGE=false,bool DOSKIPCONTENT=false,bool DOSUCCESSORS=false,bool DOPREDECESSORS=false, const bool DEBUG=false) {
        //do everything (provided that it exists in file)
        this->DOPARENTS = DOPARENTS;
        this->DOCHILDREN = DOCHILDREN;
        this->DOXCOUNT = DOXCOUNT;
        this->DOTEMPLATES = DOTEMPLATES;
		this->DOINSTANCES = DOINSTANCES;
		this->DOSKIPUSAGE = DOSKIPUSAGE;
		this->DOSKIPCONTENT = DOSKIPCONTENT;
		this->DOSUCCESSORS = DOSUCCESSORS;
		this->DOPREDECESSORS = DOPREDECESSORS;
        model->model_id = GRAPHPATTERNMODEL+GRAPHPATTERNMODELVERSION;
    	DELETEMODEL = false;        
        this->model = model;
        secondpass = true;
    	readfile(graphmodelfilename, DEBUG);        
    }
    GraphPatternModel(const std::string & graphmodelfilename, bool DOPARENTS=true,bool DOCHILDREN=true,bool DOXCOUNT=true, bool DOTEMPLATES =false, bool DOINSTANCES = false,bool DOSKIPUSAGE=false,bool DOSKIPCONTENT=false,bool DOSUCCESSORS=false,bool DOPREDECESSORS=false, const bool DEBUG=false) {
        //do everything (provided that it exists in file)
        this->DOPARENTS = DOPARENTS;
        this->DOCHILDREN = DOCHILDREN;
        this->DOXCOUNT = DOXCOUNT;    
        this->DOTEMPLATES = DOTEMPLATES;
		this->DOINSTANCES = DOINSTANCES;
		this->DOSKIPUSAGE = DOSKIPUSAGE;
		this->DOSKIPCONTENT = DOSKIPCONTENT;
		this->DOSUCCESSORS = DOSUCCESSORS;
		this->DOPREDECESSORS = DOPREDECESSORS;
    	
    	DELETEMODEL = true;
    	model = new IndexedPatternModel();
    	model->model_id = GRAPHPATTERNMODEL+GRAPHPATTERNMODELVERSION;
    	std::cerr << "Pass one, reading implied indexedpatternmodel..." << std::endl;
    	//reading is done in two passes
    	secondpass = false;
    	readfile(graphmodelfilename, DEBUG);    
        std::cerr << "Pass two, reading graph data..." << std::endl;
        secondpass = true;
        readfile(graphmodelfilename, DEBUG);
    }    
        
    
    ~GraphPatternModel();
    
    int transitivereduction();
    
    int xcount(const EncAnyGram* anygram); //exclusive count    
        
    void save(const std::string & filename) { writefile(filename); }
    
    virtual uint64_t id() { return GRAPHPATTERNMODEL + GRAPHPATTERNMODELVERSION; }
    
    
    
    virtual void readheader(std::istream * in, bool ignore = false);
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out);
    virtual void writengrams(std::ostream * out);
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writefooter(std::ostream * out) {};    
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT);
    void stats(std::ostream *OUT);
    void outputgraph(ClassDecoder & classdecoder, std::ostream *OUT);
    void outputgraph(ClassDecoder & classdecoder, std::ostream *OUT, const EncAnyGram *);
    void outputgraphvizrelations( const EncAnyGram * anygram, std::ostream *OUT, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> > & relationhash, const std::string & colour);
    void outputgraphvizrelations( const std::unordered_set<const EncAnyGram *> &, std::ostream *OUT, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> > & relationhash, const std::string & colour);    
  
    void outputrelations(ClassDecoder & classdecoder, std::ostream *OUT, const EncAnyGram * focusinput);
    void outputrelations(ClassDecoder & classdecoder, std::ostream *OUT, std::unordered_set<const EncAnyGram*>   & relations );

    void outputcoverage(ClassDecoder & classdecoder, std::ostream *OUT);
    
    void findincomingnodes(const EncAnyGram * focus, std::unordered_set<const EncAnyGram *> & relatednodes);
    void findincomingnodes(const EncAnyGram * focus, const EncAnyGram * anygram, std::unordered_set<const EncAnyGram *> & relatednodes, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> >  & relationhash );
};


class AlignConstraintInterface {
	public: 
	  virtual const EncAnyGram * getsourcekey(const EncAnyGram* key) =0;
      virtual const EncAnyGram * gettargetkey(const EncAnyGram* key) =0;
};

class IndexCountData {
	public:
	 uint32_t count;
	 uint32_t xcount;
	 std::multiset<uint32_t> sentences; //may occur multiple times in same sentence 
};


class SelectivePatternModel: public ModelReader, public ModelQuerier {
    // Read only model, reads graphpatternmodel/indexedmodel/unindexedmodel in a simplified, selective, less memory intensive representation. For for example alignment tasks, supports double indexes if fed an indexed model, and exclusive counts if fed a graph model. Whilst offering more functionality, it is also limited in the sense that it does not offer the full representation the complete model does.
    private:
     bool HASPARENTS;
     bool HASCHILDREN;
     bool HASXCOUNT;
     bool HASTEMPLATES;
     bool HASINSTANCES;
     bool HASSKIPUSAGE;
     bool HASSKIPCONTENT;
     bool HASSUCCESSORS;
     bool HASPREDECESSORS;
 
     
     bool DOXCOUNT;
     bool DOFORWARDINDEX;
     bool DOREVERSEINDEX;
     bool DOPARENTS;
     bool DOCHILDREN;
     bool DOTEMPLATES;
     bool DOINSTANCES;
     bool DOSKIPUSAGE;
     bool DOSKIPCONTENT;
     bool DOSUCCESSORS;
     bool DOPREDECESSORS;
     
     
     
     int COUNTTHRESHOLD;
     double FREQTHRESHOLD;
     double XCOUNTTHRESHOLD;
     double XCOUNTRATIOTHRESHOLD;
     
     int MINLENGTH;
     int MAXLENGTH;
     bool DOSKIPGRAMS;
    
     int ngramtypecount;
     int skipgramtypecount;   
     
     bool secondpass; 
     
     AlignConstraintInterface * alignconstrain;
     bool alignconstrainsource;
     
     void readrelations(std::istream * in, const EncAnyGram * anygram = NULL, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > * relationhash = NULL, bool ignore = false);
    public:
     unsigned long ngramtokencount;
     unsigned long skipgramtokencount;  
     unsigned long ignoredtypes;
     unsigned long ignoredtokens;
     std::unordered_map<const EncNGram, IndexCountData> ngrams;
     std::unordered_map<const EncSkipGram,IndexCountData> skipgrams;
     
     std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_parents;
     std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_children;     
    

              
    
     std::unordered_map<uint32_t,std::vector<const EncAnyGram*> > reverseindex;    
     SelectivePatternModel(const std::string & filename, bool DOFORWARDINDEX = true, bool DOREVERSEINDEX = true, bool DOXCOUNT = true, int COUNTTHRESHOLD = 0, double FREQTHRESHOLD = 0, double XCOUNTRATIOTHRESHOLD = 0, int XCOUNTTHRESHOLD = 0, bool DOSKIPGRAMS = true,  int MINLENGTH = 0, int MAXLENGTH=99, bool DOPARENTS = false, bool DOCHILDREN = false, AlignConstraintInterface * alignconstrain = NULL, bool alignconstrainsource = true , const bool DEBUG=false); //read a graph pattern model
  
     uint64_t types() const { return ngrams.size() + skipgrams.size(); }
     uint64_t tokens() const { return totaltokens; }
     
    virtual uint64_t id() { return 0; }
    
    int maxlength() const { return MAXLENGTH; }
    const EncAnyGram* getkey(const EncAnyGram* key);
    int count(const EncAnyGram* key);
    double freq(const EncAnyGram* key);    
    int xcount(const EncAnyGram* key);
    double xcountratio(const EncAnyGram* key);
    
    int countforsentence(const EncAnyGram* key, const uint64_t sentence);
    
	void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &);    
    
    
    
    int transitivereduction();
    
    bool has_xcount() { return (HASXCOUNT); }
    bool has_index() { return ((model_id != UNINDEXEDPATTERNMODEL)) ; }
    bool has_parents() { return (HASPARENTS) ; }
    bool has_children() { return (HASCHILDREN) ; }
    
    virtual void readheader(std::istream * in, bool ignore = false);
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};
};

namespace std {
    template <>
    struct hash<CorpusReference> {
     public: 
            size_t operator()(CorpusReference ref) const throw() {            
                return (ref.sentence * 1000) + ref.token;
            }
    };    
}

#endif
