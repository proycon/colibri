#include <patternmodel.h>
#include <list>
#include <gizamodel.h>

enum CoocMode {
	NOCOOC = 0,
	JACCARD = 1,
	DICE = 2,
	QUICK = 3,
};

void orderedinsert(std::list<double> &, double value);
void recompute_token_index(std::unordered_map<const EncAnyGram *, std::vector<int> > & tokenfwindex, std::unordered_map<int, std::vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, const std::vector<const EncAnyGram*> * patterns );

class AlignmentModel: public AlignConstraintInterface {
   protected:
    bool DEBUG;
    std::unordered_map<const EncNGram,bool> sourcengrams;
    std::unordered_map<const EncNGram,bool> targetngrams;
    std::unordered_map<const EncSkipGram,bool> sourceskipgrams;
    std::unordered_map<const EncSkipGram,bool> targetskipgrams;
   public:
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel; 
   
    AlignmentModel() { DEBUG = false; }
    AlignmentModel(const std::string & filename, const int bestn = 0);
    AlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const bool DEBUG = false);
    
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignmatrix;    
    virtual void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    virtual void simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool targetfirst = false, bool wordbased = false, bool mosesformat = false);
    void enabledebug() { DEBUG = true; }
    
    const EncAnyGram * getsourcekey(const EncAnyGram* key);
    const EncAnyGram * gettargetkey(const EncAnyGram* key);
    
    unsigned int totalsize() {
    	unsigned int c = 0;
    	for (std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    		c += iter->second.size();
    	}
    	return c;
	}
	
	void intersect(AlignmentModel * reversemodel, double probthreshold = 0, int bestn = 0); //Compute intersection with reverse model	
	int graphalign(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double impactfactor = 1.2);
	
	
	double cooc( CoocMode mode, const std::multiset<uint32_t> & sourceindex, const std::multiset<uint32_t> & targetindex,  const double threshold = 0); //multiset instead of vector cause we want the ordering to easily compute co-occurence
	
	
	int prune(const double prunethreshold);
	void normalize();
	
	unsigned int trainCooc(CoocMode mode, const int bestn = 0, const double absthreshold = 0,  const double relthreshold = 0);
	unsigned int trainCooc(CoocMode mode, const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel, const int bestn = 0, const double absthreshold = 0,  const double relthreshold = 0, const bool normalize = false);
		
	void trainEM(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, const bool DONULL = true, const bool INIT=true);
	void trainEM2(const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, const int bestn = 0, const bool DONULL = true, const bool INIT=true);	
	

	int extractgizapatterns(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold = 0.5, ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL); //classdecoders only for verbose output
	int extractgizapatterns(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold=0, const double coocthreshold=0, const double alignscorethreshold=0.5,  ClassDecoder * sourcedecoder = NULL, ClassDecoder * targetdecoder = NULL);
	
	
	void save(const std::string & filename);	
};

class BiAlignmentModel: public AlignmentModel {
   public:
    BiAlignmentModel(const std::string & sourcefilename, const std::string & targetfilename);
    
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignmatrixrev;
        
    virtual void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    virtual void simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool targetfirst = false, bool wordbased = false, bool mosesformat = false);    
	
	
};


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
*/


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
