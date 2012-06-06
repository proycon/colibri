#ifndef GIZAMODEL_H
#define GIZAMODEL_H

#include "ngram.h"
#include "classencoder.h"
#include "classdecoder.h"
#include <map>

class GizaSentenceAlignment {
    public:
     int index;
     EncNGram * source; //holds entire sentence 
     EncNGram * target; //holds entire sentence
     //std::unordered_map<unsigned char , std::vector<unsigned char> > alignment;
     std::multimap<const unsigned char,const unsigned char> alignment;

     GizaSentenceAlignment(const std::string & sourceline,const std::string & targetline, ClassEncoder * classencoder, const int index);     
     //GizaSentenceAlignment(const EncNGram * source, const EncNGram * target, const int index = 0);
     GizaSentenceAlignment(const GizaSentenceAlignment& ref);
     ~GizaSentenceAlignment();
     
     void parsesource(const std::string & line, ClassEncoder * classencoder);
     void parsetarget(const std::string & line, ClassEncoder * classencoder);
     
     GizaSentenceAlignment intersect(const GizaSentenceAlignment & other);  
};

class GizaModel {
    private:
        std::ifstream * IN;
        ClassEncoder * classencoder;
        ClassDecoder * classdecoder;
        int sentenceindex; 
    public:
                
        GizaModel(const std::string & filename, ClassEncoder * classencoder, ClassDecoder * classdecoder = NULL);
        ~GizaModel();
                
        bool eof() const { return IN->eof(); };
        
        GizaSentenceAlignment readsentence();
        
        
};

#endif
