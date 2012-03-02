#include <patternmodel.h>

enum CoocMode {
	NOCOOC = 0,
	JACCARD = 1,
	DICE = 2,
	QUICK = 3,
};


class AlignmentModel: public AlignConstraintInterface {
   protected:
    bool DEBUG;
    std::unordered_map<const EncNGram,bool> sourcengrams;
    std::unordered_map<const EncNGram,bool> targetngrams;
    std::unordered_map<const EncSkipGram,bool> sourceskipgrams;
    std::unordered_map<const EncSkipGram,bool> targetskipgrams;
   public:
    AlignmentModel() { DEBUG = false; }
    AlignmentModel(const std::string & filename);
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignmatrix;    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    void simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT, bool targetfirst = false, bool wordbased = false, bool mosesformat = false);
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
	
	void intersect(AlignmentModel * reversemodel, double probthreshold = 0); //Compute intersection with reverse model
	
	
	int graphalign(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double impactfactor = 1.2);
	
	virtual void save(const std::string & filename) {};	
};


class CoocAlignmentModel: public AlignmentModel {
   private:
    double absthreshold; //cooc threshold
    double probthreshold;
    bool normalize;
    bool bestonly; 
   public:   
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel;
    CoocMode mode;
    CoocAlignmentModel(CoocMode mode, SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, double absthreshold = 0,  const double relthreshold = 0, bool BESTONLY = false, bool DONORM = true, bool DEBUG = false);         
   
    void save(const std::string & filename);
   
    double cooc( const std::multiset<uint32_t> & sourceindex, const std::multiset<uint32_t> & targetindex,  const double threshold = 0); //multiset instead of vector cause we want the ordering to easily compute co-occurence 
    unsigned int compute(const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel);
    

    //void save(const std::string filename);
};


class EMAlignmentModel: public AlignmentModel {
   public:
    SelectivePatternModel * sourcemodel;
    SelectivePatternModel * targetmodel;    
    EMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001, double threshold = 0.0, bool BESTONLY = false, bool DEBUG = false);        
    //save(const std::string filename);
};
