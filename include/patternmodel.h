#ifndef PATTERNMODEL_H
#define PATTERNMODEL_H

#include <ngram.h>
#include "classencoder.h"
#include <cmath>
#include <cstdint>
#include <map>

const unsigned char MAXN = 0xff;

const short GRAPHPATTERNMODELVERSION = 2; //unsigned:1
const short INDEXEDPATTERNMODELVERSION = 3; //unsigned:2
const short UNINDEXEDPATTERNMODELVERSION = 2; //unsigned:1

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
    CorpusReference() { sentence=0; token = 0; } //needed this to fix a compiler error, but not sure why.. don't really want it   
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

class ExtCorpusReference {
    public:
     uint32_t sentence;
     unsigned char token;
     char n;
     ExtCorpusReference(const CorpusReference & ref, char n) { sentence = ref.sentence; token = ref.token; this->n = n; };
     
     bool contains(const ExtCorpusReference& other) const {
        if (sentence != other.sentence) return false;        
        return ((other.token >= token) && (other.token + other.n <= token + n));
     }
     
     bool operator< (const ExtCorpusReference& other) const {
        if (sentence < other.sentence) {
            return true;
        } else if (sentence == other.sentence) {
            return (token < other.token);
        } else {
            return false;
        }
    } 
    
    bool operator> (const ExtCorpusReference& other) const {
        if (sentence > other.sentence) {
            return true;
        } else if (sentence == other.sentence) {
            return (token + n > other.token + other.n);
        } else {
            return false;
        }
    } 
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
    
    uint64_t totaltypes;
    bool DEBUG; 
    
    int FOUNDMINN; 
    int FOUNDMAXN; 
   public:
    uint64_t totaltokens; //INCLUDES TOKENS NOT COVERED BY THE MODEL!
   
    uint64_t model_id;
    virtual uint64_t id() =0;
    virtual void readheader(std::istream * in, bool ignore=false) =0;
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, int ngramversion=1, bool ignore = false) =0;
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, int ngramversion=1, bool ignore = false) =0;
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

class ModelQuerierBase {
    std::map<int, std::vector< std::vector< std::pair<int,int> > > > gapconf;
   public:
    ModelQuerierBase();
    virtual const EncAnyGram * getfocuskey(const EncAnyGram * anygram) { return getkey(anygram); }; //without context, defaults to getkey (for models not supporting context). getpatterns uses only this! The extracted patterns will always be without context
    virtual const EncAnyGram * getkey(const EncAnyGram * anygram) =0; //default, includes context if available
    std::vector<std::pair<const EncAnyGram*, CorpusReference> > getpatterns(const unsigned char * data, const unsigned char datasize, bool doskipgrams=true, uint32_t linenum=0, const int minn = 1, const int maxn = MAXN);
};

class ModelQuerier: public ModelQuerierBase {
	public:
	 //ModelQuerier() {} ;
	 virtual int maxlength() const =0;
     virtual int occurrencecount(const EncAnyGram* key) =0;
     virtual int coveragecount(const EncAnyGram* key) =0;    
     virtual double coverage(const EncAnyGram* key) =0;	 
	 virtual void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &) =0;	 	 	  
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


    void computestats(); //compute occurrence count sums
   public:

    //occurence counts
    unsigned long totalngramcount;
    unsigned long totalskipgramcount;     
    unsigned int ngramcount[MAXN]; 
    unsigned int skipgramcount[MAXN];    
    unsigned int ngramtypes[MAXN]; 
    unsigned int skipgramtypes[MAXN];
    
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
    uint64_t occurrences() const { return totalngramcount + totalskipgramcount; }    
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);
    const AnyGramData* getdata(const EncAnyGram* key);
    
    
    int occurrencecount(const EncAnyGram* key);
    int coveragecount(const EncAnyGram* key);    
    double coverage(const EncAnyGram* key);    

    
    void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &);	 	 

    
    /* Reverse SENTENCE index */
    std::set<int> reverse_index_keys(); 
    bool reverse_index_haskey(const int i) const;    
    int reverse_index_size(const int i);
    int reverse_index_size();
    std::vector<EncAnyGram*> reverse_index(const int i);
    EncAnyGram* get_reverse_index_item(const int, const int);
    
    virtual uint64_t id() { return INDEXEDPATTERNMODEL + INDEXEDPATTERNMODELVERSION; }
    virtual void readheader(std::istream * in, bool ignore = false);
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, int ngramversion=1,bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, int ngramversion=1,bool ignore = false);
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
    
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT, bool outputhash=false);
    void decode(IndexedPatternModel & testmodel, ClassDecoder & classdecoder, std::ostream *OUT, bool outputhash=false);
    
    bool skipgramvarietycheck(SkipGramData & skipgramdata, int mintypecount=2);
    void coveragereport(std::ostream *OUT, const std::string & corpusfile = "", std::ostream *HTMLOUT = NULL, ClassDecoder * decoder = NULL, int segmentsize = 100000);
    void histogram(std::ostream *OUT);
    void report(std::ostream *OUT);
    //unsigned int prunebyalignment(std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > & alignmatrix, double threshold = 0.0);
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


    void computestats(); //compute occurrence count sums
   public:

    //occurence counts
    unsigned long totalngramcount;
    unsigned long totalskipgramcount;     
    unsigned int ngramcount[MAXN]; 
    unsigned int skipgramcount[MAXN];    
    unsigned int ngramtypes[MAXN]; 
    unsigned int skipgramtypes[MAXN];
    
    
    std::unordered_map<const EncNGram,uint32_t > ngrams;
    std::unordered_map<const EncSkipGram,uint32_t > skipgrams;    
            
    UnindexedPatternModel(const std::string & filename, const bool DEBUG=false);
    UnindexedPatternModel(const std::string & corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    UnindexedPatternModel(const std::string & corpusfile, UnindexedPatternModel & refmodel, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2,  bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
        
    int maxlength() const { return MAXLENGTH; }
    
    uint64_t types() const { return ngrams.size() + skipgrams.size(); }
    uint64_t tokens() const { return totaltokens; }
    uint64_t occurrences() const { return totalngramcount + totalskipgramcount; }
    
    bool exists(const EncAnyGram* key) const;
    const EncAnyGram* getkey(const EncAnyGram* key);


    int occurrencecount(const EncAnyGram* key);
    int coveragecount(const EncAnyGram* key);    
    double coverage(const EncAnyGram* key);           
    //std::set<int> * index(const EncAnyGram* key);    
    //int index_size() const;

	void outputinstance(const EncAnyGram *, CorpusReference, ClassDecoder &);
    
    virtual uint64_t id() { return UNINDEXEDPATTERNMODEL; }
    virtual void readheader(std::istream * in,  bool ignore = false) {};
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, int ngramversion=1,bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, int ngramversion=1,bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out) {};
    virtual void writengrams(std::ostream * out);
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writefooter(std::ostream * out) {}; 
        
    void save(const std::string & filename) { ModelWriter::writefile(filename); }
    
    size_t hash();
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT, bool outputhash=false);
    void decode(UnindexedPatternModel & testmodel, ClassDecoder & classdecoder, std::ostream *OUT, bool outputhash=false);
    
    void report(std::ostream *OUT);
    void histogram(std::ostream *OUT);
        
};


//typedef std::unordered_map<EncNGram,set<CorpusReference> > freqlist;
//typedef std::unordered_map<EncSkipGram,set<CorpusReference> > skipgram_freqlist;


int transitivereduction(std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relations );


class GraphFilter {
   public:
    bool DOPARENTS;
    bool DOCHILDREN;
    bool DOXCOUNT;
    bool DOTEMPLATES;
    bool DOINSTANCES;
    bool DOSKIPUSAGE;
    bool DOSKIPCONTENT;
    bool DOSUCCESSORS;
    bool DOPREDECESSORS;
  
  GraphFilter() {
    DOPARENTS = false;
    DOCHILDREN = false;
    DOXCOUNT = false;
    DOTEMPLATES = false;
    DOINSTANCES = false; 
    DOSKIPUSAGE = false;
    DOSKIPCONTENT = false;
    DOSUCCESSORS = false;
    DOPREDECESSORS = false;
  }
};    


class GraphRelations {
   public:
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
    
    bool TRANSITIVE;
    bool secondpass;          
  
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_parents;
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_subsumption_children;        
    std::unordered_map<const EncAnyGram*,int> data_xcount;        
   
   
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_templates; //instance -> skipgram
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_instances; //skipgram -> instance
    
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipusage; //skipcontent -> skipgram
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_skipcontent; //skipgram -> skipcontent       
    
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_successors;  
    std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > rel_predecessors;  
    
    void readrelations(std::istream * in,const EncAnyGram* = NULL, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > * = NULL, int ngramversion=1,bool ignore = false);
    
    void getrelations(std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relations, const EncAnyGram * anygram, std::unordered_set<const EncAnyGram*> & container);
    
    int transitivereduction();
    
    bool has_xcount() { return (HASXCOUNT); }
    bool has_parents() { return (HASPARENTS) ; }
    bool has_children() { return (HASCHILDREN) ; }
    bool has_templates() { return (HASTEMPLATES) ; }
    bool has_instances() { return (HASINSTANCES) ; }
    bool has_skipusage() { return (HASSKIPUSAGE) ; }
    bool has_skipcontent() { return (HASSKIPCONTENT) ; }
    bool has_successors() { return (HASSUCCESSORS) ; }
    bool has_predecessors() { return (HASPREDECESSORS) ; }

  void applyfilter(const GraphFilter & model) {
    DOPARENTS = model.DOPARENTS;
    DOCHILDREN = model.DOCHILDREN;
    DOXCOUNT = model.DOXCOUNT;
    DOTEMPLATES = model.DOTEMPLATES;
    DOINSTANCES = model.DOINSTANCES;
    DOSKIPUSAGE = model.DOSKIPUSAGE;
    DOSKIPCONTENT = model.DOSKIPCONTENT;
    DOPREDECESSORS = model.DOPREDECESSORS;
    DOSUCCESSORS = model.DOSUCCESSORS;  
  }
    
   virtual const EncAnyGram* getkey(const EncAnyGram* key) =0;
};

class GraphPatternModel: public ModelReader, public ModelWriter, public GraphRelations {               
   protected:

    
    bool DELETEMODEL;
    
    //void readrelations(std::istream * in,const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > &, bool ignore = false);
    void writerelations(std::ostream * out, const EncAnyGram*, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & );
   public:
   
    IndexedPatternModel * model;

    
    uint64_t types() const { return model->types(); }
    uint64_t tokens() const { return model->tokens(); }
    uint64_t occurrences() const { return model->occurrences(); }
    

   
    GraphPatternModel(IndexedPatternModel * model, const GraphFilter & filter); //compute entire model
    GraphPatternModel(const std::string & graphmodelfilename, IndexedPatternModel * model, const GraphFilter & filter, const bool DEBUG=false) {
        //do everything (provided that it exists in file)
        applyfilter(filter);
        model->model_id = GRAPHPATTERNMODEL+GRAPHPATTERNMODELVERSION;
    	DELETEMODEL = false;        
        this->model = model;
        secondpass = true;
    	readfile(graphmodelfilename, DEBUG);        
    }
    GraphPatternModel(const std::string & graphmodelfilename, const GraphFilter & filter, const bool DEBUG=false) {
        //do everything (provided that it exists in file)
        applyfilter(filter);    	
    	DELETEMODEL = true;
    	model = new IndexedPatternModel();
    	model->model_id = GRAPHPATTERNMODEL+GRAPHPATTERNMODELVERSION;
    	std::cerr << "Pass one, reading implied indexedpatternmodel..." << std::endl;
    	//reading is done in two passes
    	secondpass = false;
    	readfile(graphmodelfilename, DEBUG);
    	model->totaltokens = totaltokens;    
        std::cerr << "Pass two, reading graph data..." << std::endl;
        secondpass = true;
        readfile(graphmodelfilename, DEBUG);
    }    
        
    
    
    ~GraphPatternModel();
    
    int transitivereduction();
    
    int computexcount(const EncAnyGram* anygram); //exclusive count    
        
    void save(const std::string & filename) { writefile(filename); }
    
    virtual uint64_t id() { return GRAPHPATTERNMODEL + GRAPHPATTERNMODELVERSION; }
        
    
    const EncAnyGram* getkey(const EncAnyGram* key) { return model->getkey(key); }
    
    virtual void readheader(std::istream * in, bool ignore = false);
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, int ngramversion=1,bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram, int ngramversion=1,bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};    
    
    virtual void writeheader(std::ostream * out);
    virtual void writengrams(std::ostream * out);
    virtual void writengramdata(std::ostream * out, const EncNGram & ngram);
    virtual void writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram);
    virtual void writeskipgrams(std::ostream * out);
    virtual void writefooter(std::ostream * out) {};    
    
    void decode(ClassDecoder & classdecoder, std::ostream *OUT, bool outputrelations = false);
    void stats(std::ostream *OUT);
    void coveragereport(std::ostream *OUT, int segmentsize = 100000, double xratiothreshold = 0.8);
    
    void outputgraph(ClassDecoder & classdecoder, std::ostream *OUT);
    void outputgraph(ClassDecoder & classdecoder, std::ostream *OUT, const EncAnyGram *);
    void outputgraphvizrelations( const EncAnyGram * anygram, std::ostream *OUT, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> > & relationhash, const std::string & colour);
    void outputgraphvizrelations( const std::unordered_set<const EncAnyGram *> &, std::ostream *OUT, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> > & relationhash, const std::string & colour);    
  
    void outputrelations(ClassDecoder & classdecoder, std::ostream *OUT, const EncAnyGram * focusinput, bool outputquery=false);
    void outputrelations(ClassDecoder & classdecoder, std::ostream *OUT, std::unordered_set<const EncAnyGram*>   & relations );

    void outputcoverage(ClassDecoder & classdecoder, std::ostream *OUT);
    
    void findincomingnodes(const EncAnyGram * focus, std::unordered_set<const EncAnyGram *> & relatednodes);
    void findincomingnodes(const EncAnyGram * focus, const EncAnyGram * anygram, std::unordered_set<const EncAnyGram *> & relatednodes, std::unordered_map<const EncAnyGram *, std::unordered_set<const EncAnyGram*> >  & relationhash );
};


class AlignConstraintInterface {
	public: 
	  virtual const EncAnyGram * getsourcekey(const EncAnyGram* key, bool allowfallback=true) =0;
      virtual const EncAnyGram * gettargetkey(const EncAnyGram* key) =0;
};

class IndexCountData {
	public:
	 uint32_t count;
	 uint32_t xcount;
	 std::multiset<uint32_t> sentences; //may occur multiple times in same sentence 
};


class SelectivePatternModel: public ModelReader, public ModelQuerier, public GraphRelations {
    // Read only model, reads graphpatternmodel/indexedmodel/unindexedmodel in a simplified, selective, less memory intensive representation. For for example alignment tasks, supports double indexes if fed an indexed model, and exclusive counts if fed a graph model. Whilst offering more functionality, it is also limited in the sense that it does not offer the full representation the complete model does.
    private:
     bool DOFORWARDINDEX;
     bool DOREVERSEINDEX;

     
     
     
     int COUNTTHRESHOLD;
     double FREQTHRESHOLD;
     double XCOUNTTHRESHOLD;
     double XCOUNTRATIOTHRESHOLD;
     
     int MINLENGTH;
     int MAXLENGTH;
     bool DOSKIPGRAMS;
    
     int ngramtypecount;
     int skipgramtypecount;   
     
    
     AlignConstraintInterface * alignconstrain;
     bool alignconstrainsource;
     
     //void readrelations(std::istream * in, const EncAnyGram * anygram = NULL, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > * relationhash = NULL, bool ignore = false);
    
     void computestats(); //compute occurrence count sums
   public:

    //occurence counts
    unsigned long totalngramcount;
    unsigned long totalskipgramcount;     
    unsigned int ngramcount[MAXN]; 
    unsigned int skipgramcount[MAXN];    
    unsigned int ngramtypes[MAXN]; 
    unsigned int skipgramtypes[MAXN];
    
      
     unsigned long ignoredtypes;
     unsigned long ignoredoccurrences;
     std::unordered_map<const EncNGram, IndexCountData> ngrams;
     std::unordered_map<const EncSkipGram,IndexCountData> skipgrams;
     
              
    
     std::unordered_map<uint32_t,std::vector<const EncAnyGram*> > reverseindex;    
     SelectivePatternModel(const std::string & filename, const GraphFilter & filter, bool DOFORWARDINDEX = true, bool DOREVERSEINDEX = true, int COUNTTHRESHOLD = 0, double FREQTHRESHOLD = 0, double XCOUNTRATIOTHRESHOLD = 0, int XCOUNTTHRESHOLD = 0, bool DOSKIPGRAMS = true,  int MINLENGTH = 0, int MAXLENGTH=99, AlignConstraintInterface * alignconstrain = NULL, bool alignconstrainsource = true , const bool DEBUG=false); //read a graph pattern model
  
     uint64_t types() const { return ngrams.size() + skipgrams.size(); }
     uint64_t tokens() const { return totaltokens; }
     
    virtual uint64_t id() { return 0; }
    
    int maxlength() const { return MAXLENGTH; }
    const EncAnyGram* getkey(const EncAnyGram* key);
    
    
    int occurrencecount(const EncAnyGram* key);
    int coveragecount(const EncAnyGram* key);    
    double coverage(const EncAnyGram* key);
        
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
    virtual void readngramdata(std::istream * in, const EncNGram & ngram, int ngramversion=1,bool ignore = false);
    virtual void readskipgramdata(std::istream * in, const EncSkipGram & skipgram,int ngramversion=1,bool ignore = false);
    virtual void readfooter(std::istream * in, bool ignore = false) {};
};




void readcorpus( const std::string & corpusfile, std::unordered_map<CorpusReference, const EncNGram *> & tokens);
void replaceAll(std::string& str, const std::string& from, const std::string& to);

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
