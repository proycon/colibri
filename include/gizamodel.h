#ifndef GIZAMODEL_H
#define GIZAMODEL_H

#include "ngram.h"
#include "classencoder.h"
//#include "classdecoder.h"
#include <map>

class GizaSentenceAlignment {
    public:
     int index;
     EncData * source; //holds entire sentence 
     EncData * target; //holds entire sentence
     //std::unordered_map<unsigned char , std::vector<unsigned char> > alignment;
     std::multimap<const unsigned char,const unsigned char> alignment;

     GizaSentenceAlignment(const std::string & sourceline,const std::string & targetline, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, const int index);     
     //GizaSentenceAlignment(const EncNGram * source, const EncNGram * target, const int index = 0);
     GizaSentenceAlignment(const GizaSentenceAlignment& ref);
     ~GizaSentenceAlignment();
     
     void parsesource(const std::string & line, ClassEncoder * sourceencoder);
     void parsetarget(const std::string & line, ClassEncoder * targetencoder);
     
     GizaSentenceAlignment intersect(const GizaSentenceAlignment & other);
     GizaSentenceAlignment unify(const GizaSentenceAlignment & other);  
     
     void htmlout(std::ostream*, ClassDecoder & sourcedecoder, ClassDecoder & targetdecoder);
};

class GizaModel {
    private:
        std::ifstream * IN;
        ClassEncoder * sourceencoder;
        ClassEncoder * targetencoder;
        int sentenceindex; 
        std::string nextlinebuffer;
    public:
        
        
        GizaModel(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder); 
        ~GizaModel();
                
        bool eof() const { return IN->eof(); };  
        int index() const { return sentenceindex; }
              
        GizaSentenceAlignment readsentence();
                        
};

#endif
