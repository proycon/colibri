#include <patternmodel.h>

class AlignmentModel {
   public:
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);
    
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
    CoocAlignmentModel(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double absthreshold = 0,  const double relthreshold = 0);         
   
    double cooc( std::multiset<uint32_t> & sourceindex, std::multiset<uint32_t> & targetindex); //multiset instead of vector cause we want the ordering to easily compute co-occurence 
    int compute(const EncAnyGram * sourcegram, std::multiset<uint32_t> & sourceindex, SelectivePatternModel & targetmodel);
    
    
    
    //void save(const std::string filename);
};
