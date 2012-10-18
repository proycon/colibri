#include <alignmodel.h>
#include "timbl/TimblAPI.h"


/*enum ClassifierMode {
    LOCALCONTEXT_ARRAY = 1
};*/

enum ScoreHandling {
    SCOREHANDLING_IGNORE = 0, //ignore classifier score
    SCOREHANDLING_WEIGHED = 1, //add classifier output as weight to original translation model scores
    SCOREHANDLING_APPEND = 2, //append additional score to score vector
    SCOREHANDLING_REPLACE = 3 //completely replace translation score with classifier score
};


enum ClassifierType { 
    CLASSIFIERTYPE_NONE = 0, 
    CLASSIFIERTYPE_NARRAY = 1, 
    CLASSIFIERTYPE_CONSTRUCTIONEXPERTS = 2
};

ClassifierType getclassifierconf(const std::string & ID, int & contextthreshold, int & targetthreshold, bool & exemplarweight);
void writeclassifierconf(const std::string & ID, ClassifierType, int contextthreshold, int targetthreshold, bool exemplarweight);

class Classifier {
  private:
    int featurevectorsize; //nr of unigram features (each feature is a single word) in the feature vector
    bool exemplarweights;
    std::ofstream outputfile;
    std::string trainfile;
    std::string ibasefile;
    std::string wgtfile;           
    std::string ID;
    bool appendmode;
    ClassDecoder * sourceclassdecoder;
    ClassDecoder * targetclassdecoder;
    ClassEncoder * targetclassencoder;
    Timbl::TimblAPI * testexp;
    bool added;
    bool DEBUG;
  public:
    Classifier(const std::string & id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder,bool append = false, bool exemplarweights = true, bool debug=false); //for building
    Classifier(const std::string & id, const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, bool debug); //for testing            
    ~Classifier();
    void addinstance(std::vector<const EncAnyGram *> featurevector, const EncAnyGram * label, double exemplarweight = 1);
    void addinstance(std::vector<std::string> & featurevector, const std::string & label, double exemplarweight = 1);
    void train(const std::string & timbloptions);
    const std::string id() { return ID; };
    bool empty() { return !added; }
    t_aligntargets classify(std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions );
    t_aligntargets classify(std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);     
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
        virtual void load( const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG =0) =0;
        virtual void classifyfragments(const EncData & input, AlignmentModel * original, t_sourcefragments & sourcefragments, ScoreHandling scorehandling, int contextthreshold = 2, int targetthreshold = 2); //decoder will call this  
        virtual t_aligntargets classify(const EncAnyGram * focus,  std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) =0;
        virtual t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) =0;          
};

class NClassifierArray: public ClassifierInterface {
    protected:
        int leftcontextsize;
        int rightcontextsize;
        int contextthreshold;
        int targetthreshold;
    public:
        std::map<int, Classifier*> classifierarray;    
        NClassifierArray(const std::string & id, int leftcontextsize, int rightcontextsize, int contextthreshold = 2, int targetthreshold = 2);
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights = true);       
        void train(const std::string & timbloptions);     
        void load(const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG=0);
        t_aligntargets classify(const EncAnyGram * focus,  std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);          
        AlignmentModel * classify(AlignmentModel * original);
};

class ConstructionExperts: public ClassifierInterface {
    protected:
        int leftcontextsize;
        int rightcontextsize;
        int contextthreshold;
        int targetthreshold;
    public:
        std::map<uint64_t, Classifier*> classifierarray;    
        ConstructionExperts(const std::string & id, int leftcontextsize, int rightcontextsize, int contextthreshold = 2, int targetthreshold = 2);
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights = true);       
        void train(const std::string & timbloptions);     
        void load(const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG=0);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);          
        AlignmentModel * classify(AlignmentModel * original);
};

/*class ConstructionExperts: public ClassifierInterface {
        
};*/
