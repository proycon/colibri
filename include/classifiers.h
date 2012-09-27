#include <alignmodel.h>
#include "timbl/TimblAPI.h"

/*enum ClassifierMode {
    LOCALCONTEXT_ARRAY = 1
};*/

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
    Classifier(const std::string & id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool append = false, bool exemplarweights = true); //for building
    Classifier(const std::string & id, const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder); //for testing            
    ~Classifier();
    void addinstance(std::vector<const EncAnyGram *> featurevector, const EncAnyGram * label, double exemplarweight = 1);
    void addinstance(std::vector<std::string> & featurevector, const std::string & label, double exemplarweight = 1);
    void train(const std::string & timbloptions);
    const std::string id() { return ID; };
    t_aligntargets classify(std::vector<const EncAnyGram *> & featurevector);
    t_aligntargets classify(std::vector<std::string> & featurevector);     
};



class ClassifierInterface {
    protected:
        std::string ID;
    public:
        ClassifierInterface(const std::string & _id) {
            ID = _id;
        }
        const std::string id() { return ID; };
        virtual void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,  bool exemplarweights = true) =0;
        virtual void train(const std::string & timbloptions) =0;
        
        virtual void classifyfragments(const EncData & input, AlignmentModel * original, t_sourcefragments sourcefragments) =0; //decoder will call this, sourcefragments and newtable will be filled for decoder  
};

class NClassifierArray: public ClassifierInterface {
    protected:
        int leftcontextsize;
        int rightcontextsize;
    public:
        std::map<int, Classifier*> classifierarray;    
        NClassifierArray(const std::string & id, int leftcontextsize, int rightcontextsize);
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights = true);       
        void train(const std::string & timbloptions);     

        void classifyfragments(const EncData & input, AlignmentModel * original, t_sourcefragments sourcefragments); //decoder will call this, sourcefragments will be filled for decoder
        t_aligntargets classify(std::vector<const EncAnyGram *> & featurevector);
        t_aligntargets classify(std::vector<std::string> & featurevector);          
        AlignmentModel * classify(AlignmentModel * original);
};

/*class ConstructionExperts: public ClassifierInterface {
        
};*/
