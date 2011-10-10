#include <string>
#include "classdecoder.h"

class EncNGram {
    private:
     char _size;
    public:
     unsigned char* data;
    
     EncNGram(unsigned char* dataref, char size);
     
     ~EncNGram() {
         free(data);
     }
     
     const int n() const;
     const char size() const;
     
    std::string decode(ClassDecoder& classdecoder);
            
    EncNGram slice(const int begin,const int length);
    bool operator==(const EncNGram &other) const;
    bool operator!=(const EncNGram &other) const;
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

EncNGram getencngram(const int index, const int n, unsigned char *line, const int size);
