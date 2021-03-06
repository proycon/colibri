#ifndef ALIGNMODEL_H
#define ALIGNMODEL_H

#include <patternmodel.h>
#include <common.h>
#include <gizamodel.h>

enum CoocMode {
	NOCOOC = 0,
	JACCARD = 1,
	DICE = 2,
	QUICK = 3,
};

enum PhraseAlignHeuristic {
    PAH_S2T = 0,
    PAH_T2S = 1,
    PAH_INTERSECTION = 2,
    PAH_UNION = 3,
    PAH_GROWDIAG = 4,
    PAH_GROWDIAGFINAL = 5 
};

void orderedinsert(std::list<double> &, double value);
void recompute_token_index(std::unordered_map<const EncAnyGram *, std::vector<int> > & tokenfwindex, std::unordered_map<int, std::vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, const std::vector<const EncAnyGram*> * patterns, bool includeskipgrams = false );
void recompute_token_index(std::unordered_map<const EncAnyGram *, std::vector<int> > & tokenfwindex, std::unordered_map<int, std::vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, std::unordered_set<const EncAnyGram*> * patterns, bool includeskipgrams = false );
size_t get_templates(const EncAnyGram * anygram, SelectivePatternModel * model, std::unordered_set<const EncSkipGram *> & container);
void find_clusters(std::unordered_map<const EncSkipGram*,uint16_t> skipgrams, std::vector<std::unordered_set<const EncSkipGram*> > & clusters , SelectivePatternModel * model );

const short ALIGNMENTMODELVERSION = 7; //unsigned: 4, no-keywords: 5, ordered-keywords,always-without-context: 7
const int ALIGNMENTMODEL = 100;


typedef std::unordered_map<const EncAnyGram*, std::vector<double> > t_aligntargets;
typedef std::unordered_map<const EncAnyGram*, t_aligntargets > t_alignmatrix;


typedef std::unordered_map<const EncAnyGram*, std::unordered_set<const EncAnyGram *> > t_contexts;
typedef std::unordered_map<const EncAnyGram*, std::unordered_map<const EncAnyGram *, double> > t_keywords_source; //target -> key -> double 
typedef std::unordered_map<const EncAnyGram*, t_keywords_source > t_keywords; //source -> target -> key -> double 
typedef std::unordered_map<const EncAnyGram*, std::vector< std::unordered_set<const EncAnyGram *> > > t_keywordflags_source; // target -> key
typedef std::unordered_map<const EncAnyGram*, t_keywordflags_source > t_keywordflags; //source -> target -> [key -> double ]

class AlignmentModel: public AlignConstraintInterface, public ModelQuerierBase {
   protected:
    bool DEBUG;
    
   public:
    int ptsfield;
    //std::unordered_set<EncNGram> sourcengrams;
    //std::unordered_set<EncSkipGram> sourceskipgrams;
    std::unordered_set<EncNGram> targetngrams;
    std::unordered_set<EncSkipGram> targetskipgrams;   
   
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel; 
   
    ClassDecoder * debug_sourceclassdecoder;
    ClassDecoder * debug_targetclassdecoder;
   
    unsigned char leftsourcecontext;
    unsigned char rightsourcecontext;
    
    
    //Initialise new alignment model from pattern models (actual computation needs to be invoked explicitly, different types possible)    
    AlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, unsigned char leftsourcecontext=0, unsigned char rightsourcecontext=0, const bool DEBUG = false); //prepare to compute on the basis of two pattern models

    AlignmentModel(unsigned char leftsourcecontext=0, unsigned char rightsourcecontext=0, int ptsfield = 1, const bool DEBUG = false); //prepare to compute without pattern models (heuristic giza approach), or empty model

    //Load alignment model from (binary) file
    AlignmentModel(const std::string & filename, bool logprobs = true, int ptsfield = 1, bool allowskipgrams = true, const int bestn = 0, bool DEBUG = false, const int bestnkeywords = 100, const double keywordprobthreshold = 0.0); //load from binary file
    void load(const std::string & filename, bool logprobs = true, bool allowskipgrams = true, const int bestn = 0, const int bestnkeywords = 100, const double  keywordprobthreshold = 0.0);
        
    //Load alignment model from Moses phrasetable (text)
    AlignmentModel(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs= true, int ptsfield = 3, bool DEBUG = false); //load from Moses text file
    void load(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs= true, int ptsfield = 3); 
    
    //create on the basis of two alignment models (intersection), scores from both models will be represented in the score array (e.g p(t|s) and p(s|t))
    AlignmentModel(const std::string & s2tfilename, const std::string & t2sfilename, const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0, bool DEBUG = false); //create on the basis of two alignment models, will generate two scores: p(t|s) and p(s|t) 
    AlignmentModel(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0, bool DEBUG = false); 
    void load(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0); //take the intersection of two existing models

        
    ~AlignmentModel();
        
    size_t size() { return alignmatrix.size(); }
    
    t_alignmatrix alignmatrix;
    t_contexts sourcecontexts; //focus pattern -> [patterns_in_context]
    t_keywords keywords;     
    
    void computereverse();
    t_alignmatrix reversealignmatrix; //for computation of 2nd score    
        
    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool mosesformat = false);
    
    
    void enabledebug(ClassDecoder * s = NULL, ClassDecoder * t = NULL) { DEBUG = true; debug_sourceclassdecoder = s; debug_targetclassdecoder = t; }
    
    const EncAnyGram * getsourcekey(const EncAnyGram* key, bool allowfallback=true);
    const EncAnyGram * gettargetkey(const EncAnyGram* key, bool returnselfifnotfound=false, bool forcemodel=false);
    const EncAnyGram * getkey(const EncAnyGram* key) { return getsourcekey(key); } //alias for getsourcekey, needed by ModelQuerier
    const EncAnyGram * getfocuskey(const EncAnyGram* key); //alias for getsourcekey, wITHOUT context, needed by ModelQuerier
    const EncAnyGram * getkeywordkey(const EncAnyGram * key); //for global context keywords 
    
    std::pair<EncNGram, EncNGram> getsourcecontext(const EncAnyGram* key);
    
    unsigned int totalsize() {
    	unsigned int c = 0;
    	for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    		c += iter->second.size();
    	}
    	return c;
	}
	
	void intersect(AlignmentModel * reversemodel, double probthreshold = 0, int bestn = 0); //Compute intersection with reverse model	
	int graphalign(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double impactfactor = 1.2);
	
	
	double cooc( CoocMode mode, const std::multiset<uint32_t> & sourceindex, const std::multiset<uint32_t> & targetindex,  const double threshold = 0); //multiset instead of vector cause we want the ordering to easily compute co-occurence
	
	
	int prune(const double prunethreshold);
	void normalize(t_alignmatrix * matrix = NULL);
	
	unsigned int trainCooc(CoocMode mode, const int bestn = 0, const double absthreshold = 0,  const double relthreshold = 0);
	unsigned int trainCooc(CoocMode mode, const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel, const int bestn = 0, const double absthreshold = 0,  const double relthreshold = 0, const bool normalize = false);
		
	void trainEM(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, const bool DONULL = true, const bool INIT=true);
	void trainEM2(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, const bool DONULL = true, const bool INIT=true);	
	

	int extractgizapatterns(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold = 0.5, int computereverse = 2, ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL); //classdecoders only for verbose output
	int extractgizapatterns(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold=0.5,  int computereverse = 2, ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL);
	
	int extractgizapatterns2(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold = 0.5, int computereverse = 2, bool bestonly = true, bool weighbyalignmentscore = true, bool expandunaligned = false, bool combine = false, double unionweight = 1.05, ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL);
	int extractgizapatterns2(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold=0.5,  int computereverse = 2, bool bestonly = true, bool weighbyalignmentscore = true, bool expandunaligned = false, bool combine = false, double unionweight = 1.05, ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL);
	
	
    void addextractedpattern(const EncAnyGram * sourcepattern, const EncAnyGram * targetpattern, double score = 1.0, int computereverse =2, const EncAnyGram * sourcepatternwithcontext = NULL);	
	
	
	int extractgizapatterns_heur(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, PhraseAlignHeuristic phrasealignheuristic, int computereverse);
	int extractgizapatterns_heur(GizaSentenceAlignment & sentence_a, int sentenceindex, int computereverse);
	
	GizaSentenceAlignment extractgiza_growdiag(GizaSentenceAlignment & sentence_s2t ,GizaSentenceAlignment & sentence_t2s);
	void extractgiza_final(GizaSentenceAlignment & sentence_a ,GizaSentenceAlignment & sentence_s2t , GizaSentenceAlignment & sentence_t2s );
	
	void integratereverse(int computereverse);
	int extractskipgrams(const int absolutecoocthreshold = 2);
	
	EncAnyGram * addcontext(const EncData * sentence, const EncAnyGram * focus, int sourceindex);
	EncAnyGram * addcontext(const EncData * sentence, const EncAnyGram * focus, int sourceindex, int leftsourcecontext, int rightsourcecontext);
	
	
	unsigned int prunepatternmodel(IndexedPatternModel & patternmodel, double threshold);
	
	void save(const std::string & filename, const int bestnkeywords = 99999999);
	
	t_aligntargets sumtranslationoptions(const EncAnyGram * sourcefocus, bool debug = false);	
	AlignmentModel * removecontext();
		

    int computekeywords(ModelQuerier *, const EncAnyGram * sourcekey, const EncAnyGram * sourcegram, std::unordered_map<const EncAnyGram *, std::unordered_map<const EncAnyGram *, int> > & countmap, int absolute_threshold, double probability_threshold , int filter_threshold, int bestnkeywords=1000);

    int computekeywords(IndexedPatternModel & sourcepatternmodel, IndexedPatternModel & targetpatternmodel, int include_threshold = 1, int absolute_threshold = 3, double probability_threshold = 0.0000000001, int filter_threshold = 20, int bestnkeywords=1000);
	int computekeywords(IndexedPatternModel & sourcepatternmodel, IndexedPatternModel & targetpatternmodel, const EncAnyGram * sourcegram, int absolute_threshold = 3, double probability_threshold = 0.0000000001, int filter_threshold = 20,  int bestnkeywords=1000);   
	
    int computekeywords(SelectivePatternModel & sourcepatternmodel, SelectivePatternModel & targetpatternmodel, int include_threshold = 1, int absolute_threshold = 3, double probability_threshold = 0.0000000001, int filter_threshold = 20,  int bestnkeywords=1000);
	int computekeywords(SelectivePatternModel & sourcepatternmodel, SelectivePatternModel & targetpatternmodel, const EncAnyGram * sourcegram, int absolute_threshold = 3, double probability_threshold = 0.0000000001, int filter_threshold = 20,  int bestnkeywords=1000);   
	
	
	void stats();
};

class SourceFragmentData {
  public:
    const EncAnyGram * sourcefragment;
    CorpusReference ref;
    t_aligntargets translationoptions;        
    SourceFragmentData(const EncAnyGram * sourcefragment, const CorpusReference ref, t_aligntargets translationoptions) {
        this->sourcefragment = sourcefragment;
        this->ref = ref;
        this->translationoptions = translationoptions;
    }        
};

typedef std::vector<SourceFragmentData> t_sourcefragments;



/*
class TranslationTable: public ModelQuerierBase {
   protected:
    bool DEBUG;
    void load(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0);
   public:  
    // A translation table is an alignment model with multiple scores associated 
    std::unordered_set<EncNGram> sourcengrams;
    std::unordered_set<EncNGram> targetngrams;
    std::unordered_set<EncSkipGram> sourceskipgrams;
    std::unordered_set<EncSkipGram> targetskipgrams;
    
       
    
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, std::vector<double> > > alignmatrix; 
            
    TranslationTable(const std::string & filename, bool logprobs = true, bool allowskipgrams = true, bool DEBUG = false); //load from binary file
        
    TranslationTable(const std::string & s2tfilename, const std::string & t2sfilename, const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0, bool DEBUG = false); //create on the basis of two alignment models, will generate two scores: p(t|s) and p(s|t)
    
    TranslationTable(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0, bool DEBUG = false); //create on the basis of two already loaded alignment models, will generate two scores: p(t|s) and p(s|t)
    TranslationTable(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs= true, bool DEBUG = false); //load from Moses text file

    
    const EncAnyGram * getsourcekey(const EncAnyGram* key);
    const EncAnyGram * gettargetkey(const EncAnyGram* key);    
    const EncAnyGram * getkey(const EncAnyGram* key) { return getsourcekey(key); } //alias for getsourcekey, needed by ModelQuerier
       
    void save(const std::string & filename); //save as binary 
    //TODO void save(const std::string & filename, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder); //save as moses    
    void decode(ClassDecoder & sourcedecoder, ClassDecoder & targetdecoder, std::ostream * OUT, bool mosesformat=false); //decode
    
    
    size_t size() { return alignmatrix.size(); }
};
*/

/*
class BiAlignmentModel: public AlignmentModel {
   public:
    BiAlignmentModel(const std::string & sourcefilename, const std::string & targetfilename);
    
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignmatrixrev;
        
    virtual void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    virtual void simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool targetfirst = false, bool wordbased = false, bool mosesformat = false);    
	
	
};
*/

/*
class CoocAlignmentModel: public AlignmentModel {
   private:
    double absthreshold; //cooc threshold
    double probthreshold;
    bool normalize;
    int bestn;
   public:   
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel;
    CoocMode mode;
    CoocAlignmentModel(CoocMode mode, SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int bestn = 0, const double absthreshold = 0,  const double relthreshold = 0, bool DONORM = true, bool DEBUG = false);         
   
    void save(const std::string & filename);
   
     
    unsigned int compute(const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel);
    

    //void save(const std::string filename);
};


class EMAlignmentModel: public AlignmentModel {
   protected:
    bool INIT;
    bool DONULL;
   public:
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel;    
    EMAlignmentModel() {};
    EMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, bool INIT=true, bool DONULL=true, bool DEBUG = false);   
    void trainEM(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0);     
    void save(const std::string & filename);
};



class EMAlignmentModel2: public AlignmentModel {
   protected:
    bool INIT;
    bool DONULL;
   public:    
    EMAlignmentModel2() {};
    EMAlignmentModel2(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, bool INIT=true, bool DONULL=true, bool DEBUG = false);   
    void train(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0);     
    void save(const std::string & filename);
};

*/


/*class ItEMAlignmentModel: public EMAlignmentModel {
   public:
    ItEMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, bool DONULL=true, bool DEBUG = false);        
};*/

/*
class EMAlignmentModel3: public EMAlignmentModel {
   // barely-functional EM trial based on a weird idea, will probably be removed 
   public:
    EMAlignmentModel3(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, bool DONULL=true, bool DEBUG = false);        
    unsigned int expectation(const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel, std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > & count, std::unordered_map<const EncAnyGram*, double> & total);
};
*/

#endif
