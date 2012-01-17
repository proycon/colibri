#include <patternmodel.h>

class AlignmentModel {
   public:
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
    void decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, std::ostream * OUT);

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
    std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > alignprob;    
   
    CoocAlignmentModel(DoubleIndexedGraphPatternModel & sourcemodel, DoubleIndexedGraphPatternModel & targetmodel, double absthreshold = 0,  const double relthreshold = 0);         
   
    double cooc( std::multiset<uint32_t> & sourceindex, std::multiset<uint32_t> & targetindex); 
    int compute(const EncAnyGram * sourcegram, std::multiset<uint32_t> & sourceindex, DoubleIndexedGraphPatternModel & targetmodel);
    
    
    
    //void save(const std::string filename);
};
