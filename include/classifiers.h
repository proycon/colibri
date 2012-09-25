#include <alignmodel.h>
#include "timbl/TimblAPI.h"


/*
class Classifier {
  protected:
    ClassEncoder * targetclassencoder;        
    TimblAPI * timblexp;
  public:
    Classifier(const std::string & id, ClassEncoder * targetclassencoder);
    ~Classifier();    
    TranslationTable * classify(const EncAnyGram * sourcegram, vector<const std::string> featurevector);         
};
*/

class Classifier {
  private:
    int featurevectorsize; //nr of unigram features (each feature is a single word) in the feature vector
    bool exemplarweights;
    std::ofstream outputfile;
    std::string trainfile;
    std::string ibasefile;       
    std::string ID;
    bool appendmode;
    ClassDecoder * sourceclassdecoder;
    ClassDecoder * targetclassdecoder;
    ClassEncoder * targetclassencoder;
    Timbl::TimblAPI * testexp;
  public:
    Classifier(const std::string & id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool append = false, bool examplarweights = false); //for building
    Classifier(const std::string & id, const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder); //for testing            
    ~Classifier();
    void addinstance(std::vector<const EncAnyGram *> featurevector, const EncAnyGram * label, double exemplarweight = 1);
    void addinstance(std::vector<std::string> & featurevector, const std::string & label, double exemplarweight = 1);
    void train(const std::string & timbloptions);
    const std::string id() { return ID; };
    t_aligntargets classify(std::vector<const EncAnyGram *> featurevector);
    t_aligntargets classify(std::vector<std::string> featurevector);     
};



class ClassifierInterface {
    protected:
        std::string ID;
    public:
        ClassifierInterface(const std::string & _id) {
            ID = _id;
        }
        const std::string id() { return ID; };
        virtual void build(const std::string & traincorpusfile, const AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) =0;
        virtual void load() =0;
        virtual void train() =0;    
};

class NClassifierArray: public ClassifierInterface {
    protected:
        int leftcontextsize;
        int rightcontextsize;
    public:
        std::map<int, Classifier*> classifierarray;    
        NClassifierArray(const std::string & id, int leftcontextsize, int rightcontextsize);
        void build(const std::string & enctraincorpusfile, AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights = true);        
        void load();        
        void train();     
           
};

/*class ConstructionExperts: public ClassifierInterface {
        
};*/
