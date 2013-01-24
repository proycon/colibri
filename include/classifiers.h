#include <alignmodel.h>
#include "timbl/TimblAPI.h"


/*enum ClassifierMode {
    LOCALCONTEXT_ARRAY = 1
};*/

enum ScoreHandling { //C = set of distributions from classifier, S = set of distributions from statistical model
    SCOREHANDLING_IGNORE = 0, //ignore classifier score: (C ∩ S) ∪ S
    SCOREHANDLING_WEIGHED = 1, //add classifier output as weight to original translation model scores: (C ∩ S) ∪ S
    SCOREHANDLING_APPEND = 2, //append additional score to score vector.  (C ∩ S) ∪ S
    SCOREHANDLING_REPLACE = 3, //completely replace translation score with classifier score: C
    SCOREHANDLING_FILTEREDWEIGHED = 4 //completely replace translation score with classifier score: (C ∩ S)
};


enum ClassifierType { 
    CLASSIFIERTYPE_NONE = 0, 
    CLASSIFIERTYPE_NARRAY = 1, 
    CLASSIFIERTYPE_CONSTRUCTIONEXPERTS = 2,
    CLASSIFIERTYPE_MONO = 3
};

ClassifierType getclassifierconf(const std::string & ID, int & contextthreshold, int & targetthreshold, bool & exemplarweight, bool & singlefocusfeature);
void writeclassifierconf(const std::string & ID, ClassifierType, int contextthreshold, int targetthreshold, bool exemplarweight, bool singlefocusfeature);

class Classifier {
  private:
    int featurevectorsize; //nr of unigram features (each feature is a single word) in the feature vector
    bool exemplarweights;
    std::ofstream outputfile;
    std::string trainfile;
    std::string ibasefile;
    std::string wgtfile;           
    std::string ID;
    std::string timbloptions;
    bool appendmode;
    ClassDecoder * sourceclassdecoder;
    ClassDecoder * targetclassdecoder;
    ClassEncoder * targetclassencoder;
    Timbl::TimblAPI * testexp;
    bool added;
    bool loaded;
    bool DEBUG;
  public:
    bool used;
    Classifier(const std::string & id, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights = true, bool debug=false); //for building
    Classifier(const std::string & id, const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, bool loadondemand = false, bool debug = false); //for testing
    ~Classifier();            
    void addinstance(std::vector<const EncAnyGram *> & featurevector, const EncAnyGram * label, double exemplarweight = 1);
    void addinstance(std::vector<std::string> & featurevector, const std::string & label, double exemplarweight = 1);
    void train(const std::string & timbloptions);
    double crossvalidate(const std::string & timbloptions); //returns accuracy (measured using leave one out cross validation)
    const std::string id() { return ID; };
    bool empty() { return !added; }
    bool enabledebug() { DEBUG = true; }
    void close() { outputfile.close(); }
    void flush() { outputfile.flush(); }
    void load();
    void unload();
    void remove();
    t_aligntargets classify(std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions );
    t_aligntargets classify(std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);     
};



class ClassifierInterface {
    protected:
        std::string ID;
        int leftcontextsize;
        int rightcontextsize;
        int contextthreshold;
        int targetthreshold;  
        bool singlefocusfeature;
        bool exemplarweights;
        int DEBUG;
        ClassDecoder * sourceclassdecoder;
        ClassDecoder * targetclassdecoder;
    public:
        ClassifierInterface(const std::string & _id, int leftcontextsize, int rightcontextsize, int contextthreshold, int targetthreshold, bool exemplarweights, bool singlefocusfeature) {
            DEBUG = 0;
            sourceclassdecoder = NULL;
            targetclassdecoder = NULL;
            ID = _id;
            this->leftcontextsize = leftcontextsize;
            this->rightcontextsize = rightcontextsize;
            this->contextthreshold = contextthreshold;
            this->targetthreshold = targetthreshold;
            this->exemplarweights = exemplarweights;
            this->singlefocusfeature = singlefocusfeature;            
        }
        void enabledebug(int debug, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) {
            this->DEBUG = debug;
            this->sourceclassdecoder = sourceclassdecoder;
            this->targetclassdecoder = targetclassdecoder;
        
        }
        const std::string id() { return ID; };
        virtual void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) =0;
        virtual void add(const EncAnyGram * focus, const EncAnyGram * withcontext, t_aligntargets & targets, int leftcontextsize, int rightcontextsize, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) =0;
        virtual void train(const std::string & timbloptions) =0;
        virtual void load( const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG =0) =0;
        virtual void classifyfragments(const EncData & input, AlignmentModel * original, t_sourcefragments & sourcefragments, ScoreHandling scorehandling); //decoder will call this
        virtual t_aligntargets classifyfragment(const EncAnyGram * focus, const EncAnyGram * withcontext, t_aligntargets & reftranslationfragments, ScoreHandling scorehandling, int leftcontextsize, int rightcontextsize);
        virtual t_aligntargets classify(const EncAnyGram * focus,  std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) =0;
        virtual t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions) =0; 
        virtual void reset() {};         
};


/* Monolithic single classifier, in which focus part is joined as a single feature */
class MonoClassifier: public ClassifierInterface {      
    public:
        Classifier * classifier;  
        MonoClassifier(const std::string & id, int leftcontextsize, int rightcontextsize, int contextthreshold, int targetthreshold, bool exemplarweights, bool singlefocusfeature): ClassifierInterface(id, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature) {};
        ~MonoClassifier(); 
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder); 
        void add(const EncAnyGram * focus, const EncAnyGram * withcontext, t_aligntargets & targets, int leftcontextsize, int rightcontextsize, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);
        void train(const std::string & timbloptions);     
        void load(const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG=0);
        t_aligntargets classify(const EncAnyGram * focus,  std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);          
        AlignmentModel * classify(AlignmentModel * original);
};

class NClassifierArray: public ClassifierInterface {
    public:
        std::map<int, Classifier*> classifierarray;    
        NClassifierArray(const std::string & id, int leftcontextsize, int rightcontextsize, int contextthreshold, int targetthreshold, bool exemplarweights, bool singlefocusfeature): ClassifierInterface(id, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature) {};
        ~NClassifierArray(); 
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);       
        void add(const EncAnyGram * focus, const EncAnyGram * withcontext, t_aligntargets & targets, int leftcontextsize, int rightcontextsize, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);
        void train(const std::string & timbloptions);     
        void load(const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG=0);
        t_aligntargets classify(const EncAnyGram * focus,  std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);          
        AlignmentModel * classify(AlignmentModel * original);
};

class ConstructionExperts: public ClassifierInterface {
    public:
        double accuracythreshold;
        std::map<uint64_t, Classifier*> classifierarray;    
        ConstructionExperts(const std::string & id, int leftcontextsize, int rightcontextsize, int contextthreshold, int targetthreshold, bool exemplarweights, bool singlefocusfeature): ClassifierInterface(id, leftcontextsize, rightcontextsize, contextthreshold, targetthreshold, exemplarweights, singlefocusfeature) { accuracythreshold = 0; }; 
        ~ConstructionExperts();
        void build(AlignmentModel * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);
        void add(const EncAnyGram * focus, const EncAnyGram * withcontext, t_aligntargets & targets, int leftcontextsize, int rightcontextsize, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);
        void train(const std::string & timbloptions);     
        void load(const std::string & timbloptions, ClassDecoder * sourceclassdecoder, ClassEncoder * targetclassencoder, int DEBUG=0);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<const EncAnyGram *> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);
        t_aligntargets classify(const EncAnyGram * focus, std::vector<std::string> & featurevector, ScoreHandling scorehandling, t_aligntargets & originaltranslationoptions);          
        AlignmentModel * classify(AlignmentModel * original);
        void reset();
};

/*class ConstructionExperts: public ClassifierInterface {
        
};*/
