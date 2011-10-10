#include <string>
#include <iostream>
#include "classdecoder.h"

class EncNGram {
    private:
     char _size;
    public:
     unsigned char* data;
    
     EncNGram(const unsigned char* dataref, char size);
     
     EncNGram(const EncNGram& ref);
     
     ~EncNGram();
     
     const int n() const;
     const char size() const;
     
    std::string decode(ClassDecoder& classdecoder) const;
            
    
    bool operator==(const EncNGram &other) const;
    bool operator!=(const EncNGram &other) const;
    EncNGram & operator =(EncNGram other);    
    
    EncNGram slice(const int begin,const int length) const;
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

}

EncNGram getencngram(const int index, const int n, const unsigned char *line, const int size);
