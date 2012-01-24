#include <patternmodel.h>

enum CoocMode {
	NOCOOC = 0,
	JACCARD = 1,
	DICE = 2
};

class AlignmentModel {
   protected:
    bool DEBUG;
   public:
    AlignmentModel() { DEBUG = false; }
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    void enabledebug() { DEBUG = true; }
    
    
    unsigned int totalsize() {
    	unsigned int c = 0;
    	for (std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> >::iterator iter = alignprob.begin(); iter != alignprob.end(); iter++) {
    		c += iter->second.size();
    	}
    	return c;
	}
};

/*class EMAlignmentModel: public AlignmentModel {
   public:    
    EMAlignmentModel(IndexedPatternModel & sourcemodel, IndexedPatternModel & targetmodel, const int MAXROUNDS=10000, const double CONVERGEDTHRESHOLD=0.001);        
    //save(const std::string filename);
};*/

class CoocAlignmentModel: public AlignmentModel {
   private:
    double absthreshold;
    double relthreshold;
   public:   
    CoocMode mode;
    CoocAlignmentModel(CoocMode mode, SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double absthreshold = 0,  const double relthreshold = 0, bool DEBUG = false);         
   
    double cooc( const std::multiset<uint32_t> & sourceindex, const std::multiset<uint32_t> & targetindex); //multiset instead of vector cause we want the ordering to easily compute co-occurence 
    int compute(const EncAnyGram * sourcegram, const std::multiset<uint32_t> & sourceindex, SelectivePatternModel & targetmodel);
    
    
    
    //void save(const std::string filename);
};
