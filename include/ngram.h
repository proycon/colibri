#include <string>
#include <iostream>
#include "classdecoder.h"
#include <unordered_map>
#include <vector>

const int MAXSKIPS = 4;

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
            
    
    bool operator==(const EncNGram &other) const;
    bool operator!=(const EncNGram &other) const;
    EncNGram & operator =(EncNGram other);    
    
    EncNGram * slice(const int begin,const int length) const;
};

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size);


class EncSingleSkipGram: public EncNGram {
    private:
      char _n;
    public:
      const char n() const { return _n; };
      
    EncSingleSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn);
      
    EncSingleSkipGram(const unsigned char* dataref, const char size, const char n): EncNGram(dataref, size) {
        _n = n;        
    }     
  
};
 
 
class EncSkipGram: public EncNGram {
    public:
      char skipsize[MAXSKIPS]; //4 bytes reserved for skip size
      char skipcount; //number of skips
      const char n() const;
      
      EncSkipGram(const std::vector<EncNGram*> & dataref, const std::vector<int> & skipref, bool initialskip = false, bool finalskip = false);
      EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn);
    
      
    //EncSkipGram(const unsigned char* dataref, const char size): EncNGram(dataref, size) {      
    //}     
  
};


size_t jenkinshash(unsigned char * data, char size);

namespace std {

template <>
struct hash<EncNGram> {
 public: 
        size_t operator()(EncNGram ngram) const throw() {            
            //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
            return jenkinshash(ngram.data, ngram.size());
        }
};



template <>
struct hash<EncSingleSkipGram> {
 public: 
        size_t operator()(EncSingleSkipGram ngram) const throw() {            
            //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
            return jenkinshash(ngram.data, ngram.size());
        }
};


template <>
struct hash<EncSkipGram> {
 public: 
        size_t operator()(EncSkipGram ngram) const throw() {            
            //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
            return jenkinshash(ngram.data, ngram.size());
        }
};

}
