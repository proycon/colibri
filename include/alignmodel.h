#ifndef ALIGNMODEL_H
#define ALIGNMODEL_H

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
void recompute_token_index(std::unordered_map<const EncAnyGram *, std::vector<int> > & tokenfwindex, std::unordered_map<int, std::vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, const std::vector<const EncAnyGram*> * patterns, bool includeskipgrams = false );
size_t get_templates(const EncAnyGram * anygram, SelectivePatternModel * model, std::unordered_set<const EncSkipGram *> & container);
void find_clusters(std::unordered_map<const EncSkipGram*,uint16_t> skipgrams, std::vector<std::unordered_set<const EncSkipGram*> > & clusters , SelectivePatternModel * model );


class AlignmentModel: public AlignConstraintInterface {
   protected:
    bool DEBUG;
   public:
    std::unordered_set<EncNGram> sourcengrams;
    std::unordered_set<EncNGram> targetngrams;
    std::unordered_set<EncSkipGram> sourceskipgrams;
    std::unordered_set<EncSkipGram> targetskipgrams;   
   
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel; 
   
    AlignmentModel() { DEBUG = false; }
    AlignmentModel(const std::string & filename, const int bestn = 0);
    AlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const bool DEBUG = false);
    
    void load(const std::string & filename, const int bestn = 0);
    
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignmatrix;    
    virtual void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool mosesformat = false);
    //virtual void simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool targetfirst = false, bool wordbased = false, bool mosesformat = false);
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
	int extractskipgrams(const int absolutecoocthreshold = 2);
	
	unsigned int prunepatternmodel(IndexedPatternModel & patternmodel, double threshold);
	
	void save(const std::string & filename);	
};


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
    
    TranslationTable(const std::string & filename, bool multiscore=true); //load from binary file
    
    TranslationTable(const std::string & s2tfilename, const std::string & t2sfilename, const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0); //create on the basis of two alignment models, will generate two scores: p(t|s) and p(s|t)
    
    TranslationTable(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold = 0, const double t2sthreshold = 0, const double productthreshold = 0); //create on the basis of two already loaded alignment models, will generate two scores: p(t|s) and p(s|t)
    //TODO TranslationTable(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder); //load from Moses text file

    
    const EncAnyGram * getsourcekey(const EncAnyGram* key);
    const EncAnyGram * gettargetkey(const EncAnyGram* key);    
    const EncAnyGram * getkey(const EncAnyGram* key) { return getsourcekey(key); } //alias for getsourcekey, needed by ModelQuerier
       
    void save(const std::string & filename); //save as binary 
    //TODO void save(const std::string & filename, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder); //save as moses    
    void decode(ClassDecoder & sourcedecoder, ClassDecoder & targetdecoder, std::ostream * OUT, bool mosesformat=false); //decode
    
    
    size_t size() { return alignmatrix.size(); }
};

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

#endif
