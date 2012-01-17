#include <string>
#include <iostream>
#include <ostream>
#include <istream>
#include "classdecoder.h"
#include <unordered_map>
#include <vector>
#include <set>
#include <unordered_set>
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
     
     virtual void writeasbinary(std::ostream * out) const; //write binary output
             
     //virtual bool is_subgram() const; //TODO?
     
     //virtual bool is_supergram() const; //TODO?
     
    
};


class EncNGram: public EncAnyGram {
   public:
    EncNGram(): EncAnyGram() {}; 
    EncNGram(const unsigned char* dataref, const char size): EncAnyGram(dataref, size) {};
    EncNGram(const EncNGram& ref): EncAnyGram(ref) {};     
    EncNGram(std::istream * in);
    EncNGram * slice(const int begin,const int length) const;    
    
    
    int subngrams(std::vector<EncNGram*> & container) const;     
};

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size, const unsigned int linenum = 0);


 
 
class EncSkipGram: public EncAnyGram {
    public:
      char skipsize[MAXSKIPS]; //4 bytes reserved for skip size
      char skipcount; //number of skips
      const char n() const;
      
      EncSkipGram(const std::vector<EncNGram*> & dataref, const std::vector<int> & skipref, bool initialskip = false, bool finalskip = false);
      EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount);
      EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn);
      EncSkipGram(std::istream * in, const char gapcount = -1);
      std::string decode(ClassDecoder& classdecoder) const;
      bool out() const;
      bool isskipgram() const { return (skipcount > 0); }
      
      const char gapcount() const {
         return skipcount;
      }
      const char gapsize(char i) const {
            return skipsize[i]; //TODO: Add proper exceptions
      }
      
      int parts(std::vector<EncNGram*> & container) const; //returns all consecutive parts
      
      void writeasbinary(std::ostream * out) const; //write binary output
      
};


namespace std {

    template <>
    struct hash<const EncAnyGram> {
     public: 
            size_t operator()(const EncAnyGram anygram) const throw() {            
                return anygram.hash();
            }
    };
    
    template <>
    struct hash<const EncAnyGram*> {
     public: 
            size_t operator()(const EncAnyGram * anygram) const throw() {            
                return anygram->hash();
            }
    };
    

    
    
    template <>
    struct hash<const EncNGram> {
     public: 
            size_t operator()(const EncNGram ngram) const throw() {            
                return ngram.hash();
            }
    };
    

    template <>
    struct hash<const EncSkipGram> {
     public: 
          size_t operator()(const EncSkipGram skipgram) const throw() {                            
              return skipgram.hash();              
          }
    };

    template <>
    struct hash<EncSkipGram> {
     public: 
          size_t operator()(EncSkipGram skipgram) const throw() {                            
              return skipgram.hash();              
          }
    };    

}
