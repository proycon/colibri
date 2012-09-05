#include <alignmodel.h>
#include <TimblAPI.h>

class Classifier {
  protected:
    ClassEncoder * targetclassencoder;        
    TimblAPI * timblexp;
  public:
    Classifier(const std::string & id, ClassEncoder * targetclassencoder);
    ~Classifier();    
    TranslationTable * classify(const EncAnyGram * sourcegram, vector<const std::string> featurevector);         
};

class BuildClassifier {
  private:
    int featurevectorsize; //nr of unigram features (each feature is a single word) in the feature vector
    bool exemplarweights;
    std::ofstream outputfile;
    string trainfile;    
    std::string ID;
    bool opened;
    bool append;
  public:
    BuildClassifier(const std::string & trainfile, bool append = false, bool examplarweights = false);        
    ~BuildClassifier();
    void addinstance(vector<const std::string> featurevector, const std::string & label, double exemplarweight = 1);
    void train(const std::string & timbloptions);
    const std::string id() { return ID; };
};


class ClassifierInterface {
    protected:
        std::string ID;
    public:
        ClassifierInterface(const std::string & id) {
            ID = id;
        }
        const std::string id() { return ID; };
        virtual void build(const string & traincorpusfile, const TranslationTable * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder) =0;
        virtual void load() =0;
        virtual void train() =0;    
};

class NClassifierArray: public ClassifierInterface {
    public:
        NClassifierArray(const std::string & id): ClassifierInterface(id) {};
        void build(const string & traincorpusfile, const TranslationTable * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder);        
        void load();        
        void train();        
};

class ConstructionExperts: public ClassifierInterface {
        
};
