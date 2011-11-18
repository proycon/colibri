#include <string>
#include <iostream>
#include <ostream>
#include "classdecoder.h"
#include <unordered_map>
#include <vector>
#include <set>
#include <iomanip> // contains setprecision()

const int MAXSKIPS = 4;


class EncAnyGram {
    protected:
     char _size;    
    public:    
     unsigned char* data;
    
     EncAnyGram();
     EncAnyGram(const unsigned char* dataref, const char size);
     EncAnyGram(const EncAnyGram& ref);     
     virtual ~EncAnyGram();
     
     virtual const char n() const;
     virtual const char size() const;
     
     virtual std::string decode(ClassDecoder& classdecoder) const;
     virtual bool out() const;
     
    
     virtual bool operator==(const EncAnyGram &other) const;
     virtual bool operator!=(const EncAnyGram &other) const;
     virtual EncAnyGram & operator =(EncAnyGram other);    
    
     const size_t hash() const;
    
     //EncNGram * slice(const int begin,const int length) const;
    
     virtual bool isskipgram() const { 
         return false; 
     }
     virtual const char gapcount() const {
         return 0;
     }
     virtual const char gapsize(char i) const {
         return 0;
     }
     
};


class EncNGram: public EncAnyGram {
   public:
    EncNGram(): EncAnyGram() {}; 
    EncNGram(const unsigned char* dataref, const char size): EncAnyGram(dataref, size) {};
    EncNGram(const EncNGram& ref): EncAnyGram(ref) {};     
    
    EncNGram * slice(const int begin,const int length) const;    
};

/*

class EncNGram {
    protected:
     char _size;    
    public:    
     unsigned char* data;
    
     EncNGram();
     EncNGram(const unsigned char* dataref, const char size);
     
     EncNGram(const EncNGram& ref);     
     ~EncNGram();
     
     const char n() const;
     const char size() const;
     
    std::string decode(ClassDecoder& classdecoder) const;
    bool out() const;
    
    bool operator==(const EncNGram &other) const;
    bool operator!=(const EncNGram &other) const;
    EncNGram & operator =(EncNGram other);    
    
    EncNGram * slice(const int begin,const int length) const;
    
    virtual bool isskipgram() const { 
        return false; 
    }
};
*/

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size);


 
 
class EncSkipGram: public EncAnyGram {
    public:
      char skipsize[MAXSKIPS]; //4 bytes reserved for skip size
      char skipcount; //number of skips
      const char n() const;
      
      EncSkipGram(const std::vector<EncNGram*> & dataref, const std::vector<int> & skipref, bool initialskip = false, bool finalskip = false);
      EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount);
      EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn);
    
      std::string decode(ClassDecoder& classdecoder) const;
      bool out() const;
      bool isskipgram() const { return (skipcount > 0); }
      
      const char gapcount() const {
         return skipcount;
      }
      const char gapsize(char i) const {
            return skipsize[i]; //TODO: Add proper exceptions
      }
      
    //EncSkipGram(const unsigned char* dataref, const char size): EncNGram(dataref, size) {      
    //}     
  
};


//size_t jenkinshash(unsigned char * data, char size);

namespace std {

    template <>
    struct hash<EncAnyGram> {
     public: 
            size_t operator()(EncAnyGram anygram) const throw() {            
                return anygram.hash();
            }
    };

    
    template <>
    struct hash<EncNGram> {
     public: 
            size_t operator()(EncNGram ngram) const throw() {            
                return ngram.hash();
            }
    };
    

    template <>
    struct hash<EncSkipGram> {
     public: 
          size_t operator()(EncSkipGram skipgram) const throw() {                            
              return skipgram.hash();              
          }
           /* size_t operator()(EncSkipGram ngram) const throw() {                            
                //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
                unsigned long h;
                int i;
                bool prevnull = false;
                int skipnum = 0;
                for(h = i = 0; i < ngram.size(); ++i)
                {                
                    h += ngram.data[i];
                    h += (h << 10);
                    h ^= (h >> 6);
                    if (ngram.data[i] == 0) {
                        if (prevnull) {
                            h += 46021 + ngram.skipsize[skipnum];                        
                            h += (h << 10);
                            h ^= (h >> 6);
                            skipnum++;
                        } else {
                            prevnull = true;
                        }                    
                    } else {
                        prevnull = false;
                    }                
                }
                h += (h << 3);
                h ^= (h >> 11);
                h += (h << 15);
                return h;
            }*/
    };

}



typedef std::unordered_map<EncNGram,int> freqlist;
typedef std::unordered_map<EncSkipGram,int> skipgram_freqlist;


class skipgramdata {
   public:
    int count;
    skipgram_freqlist skips;
    skipgramdata() {
        count = 0;
    }
};

typedef std::unordered_map<EncSkipGram,skipgramdata> skipgrammap;


class EncGramModel {    
   private:
    int MINTOKENS; // = 2;
    int MINSKIPTOKENS; // = 2;
    int MINSKIPTYPES; //= 2;
    int MAXLENGTH; //= 8;
    bool DOSKIPGRAMS; //= false;
    bool DOINDEX; //= false;
    bool DOSKIPCONTENT; //= false;
    bool DOINITIALONLYSKIP; //= true;
    bool DOFINALONLYSKIP; //= true;

    unsigned long ngramtokencount;
    unsigned long skipgramtokencount; 
    int ngramtypecount;
    int skipgramtypecount;
    
    
    int tokencount[10]; //relative token count
    int skiptokencount[10];
   public:   
    std::vector<freqlist> ngrams;
    std::vector<skipgrammap> skipgrams;
    
    std::unordered_map< EncNGram,std::set<int>  > ngram_index;
    std::unordered_map< EncSkipGram,std::set<int> > skipgram_index;
        
    EncGramModel(std::string filename);
    EncGramModel(const std::string corpusfile, int MAXLENGTH, int MINTOKENS = 2, bool DOSKIPGRAMS = true, int MINSKIPTOKENS = 2, int MINSKIPTYPES = 2, bool DOINDEX = false, bool DOSKIPCONTENT = false, bool DOINITIALONLYSKIP= true, bool DOFINALONLYSKIP = true);
    
    int maxlength() const { return MAXLENGTH; }
    
    int types() const { return ngramtypecount + skipgramtypecount; }
    int tokens() const { return ngramtokencount + skipgramtokencount; }
    
    bool exists(EncNGram* key) const;
    bool exists(EncSkipGram* key) const;
    int count(EncNGram* key);
    int count(EncSkipGram* key);
    double freq(EncNGram* key);    
    double freq(EncSkipGram* key);    
    double relfreq(EncNGram* key);    
    double relfreq(EncSkipGram* key);    
    
    void save(std::string filename);
    
    size_t hash();
    
    
    void decode(ClassDecoder & classdecoder, std::ostream *NGRAMSOUT, std::ostream *SKIPGRAMSOUT);    
};



/*
class EncGramGraphModel {
    
    
   private:
    
    unordered_map<size_t,std::vector<EncNGram*> > rel_subsumption_children;
    
   public:
    EncGramGraphModel(EncGramModel& model); //compute entire model    
    
    
    
    EncGramGraphModel(std:string filename);
    save(std:string filename);
};
*/

double compute_entropy(freqlist & data, const int total);
double compute_entropy(skipgram_freqlist & data, const int total);
