#ifndef LM_H
#define LM_H

#include <ngram.h>
#include <classencoder.h>
#include <map>

class LanguageModel {
    private:
        bool DEBUG;
        int order;
    public:
        std::unordered_map<EncNGram, double> ngrams;
        std::unordered_map<EncNGram, double> backoff; //MAYBE TODO: merge with ngrams? <EncNGram, pair<double,double> > ?
        std::map<int,unsigned int> total;
        
        LanguageModel(const std::string & filename,  ClassEncoder & encoder, bool debug = false);
        
        
        double score(const EncNGram * ngram, const EncNGram * history = NULL); //returns logprob (base 10)        
        double scoreword(const EncNGram * word, const EncNGram * history = NULL); //returns logprob (base 10)
         
        
        //double score(EncNGram ngram); //returns logprob (base 10)
        //double score(EncData & data, bool fullsentence = false); //returns logprob (base 10)
        
        int getorder() { return order; }
        size_t size() { return ngrams.size(); }
};

#endif
