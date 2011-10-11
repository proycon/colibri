#include <string>
#include <iostream>
#include "classdecoder.h"
#include <unordered_map>

class EncNGram {
    protected:
     char _size;
    public:
     unsigned char* data;
    
     EncNGram();
     EncNGram(const unsigned char* dataref, const char size);
     
     EncNGram(const EncNGram& ref);
     
     ~EncNGram();
     
     const int n() const;
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
      const int n() const { return (int) _n; };
      
    EncSingleSkipGram(const EncNGram & pregap, const EncNGram & postgap);
      
    EncSingleSkipGram(const unsigned char* dataref, const char size, const char n): EncNGram(dataref, size) {
        _n = n;        
    }     
    
    //const EncNGram preskippart() const;
    
    //const EncNGram postskippart() const;
};


namespace std {

template <>
struct hash<EncNGram> {
 public: 
        size_t operator()(EncNGram ngram) const throw() {            
            //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
            unsigned long h;
            int i;
            for(h = i = 0; i < ngram.size(); ++i)
            {
                h += ngram.data[i];
                h += (h << 10);
                h ^= (h >> 6);
            }
            h += (h << 3);
            h ^= (h >> 11);
            h += (h << 15);
            return h;
        }
};



template <>
struct hash<EncSingleSkipGram> {
 public: 
        size_t operator()(EncSingleSkipGram ngram) const throw() {            
            //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
            unsigned long h;
            int i;
            for(h = i = 0; i < ngram.size(); ++i)
            {
                h += ngram.data[i];
                h += (h << 10);
                h ^= (h >> 6);
            }
            h += (h << 3);
            h ^= (h >> 11);
            h += (h << 15);
            return h;
        }
};

}
