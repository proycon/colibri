#ifndef NGRAM_H
#define NGRAM_H

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


const char MAXSKIPS = 4; //Maximum number of skips - THESE NUMBERS CAN NOT BE CHANGED WITHOUT ALSO CHANGING THE SkipConf IMPLEMENTATION!
const char MAXSKIPSIZE = 16; //Maximum length of each skip   - THESE NUMBERS CAN NOT BE CHANGED WITHOUT ALSO CHANGING THE SkipConf IMPLEMENTATION!


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
         for (int i = 0; i < _size; i++) {
         	if ((i == 0) && (data[i] == 0)) return true;
         	if ((i == _size - 1) && (data[i] == 0)) return true;
         	if ((i > 0) && (data[i] == 0) && (data[i-1] == 0)) return true;
         }
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
    int splits(std::vector<std::pair<EncNGram*, EncNGram*> > & container) const;     
};

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size, const unsigned int linenum = 0);

/*
class SkipConf { //Skip configuration
	public:
	 uint16_t value;
	 SkipConf(const uint16_t value);
	 SkipConf(const unsigned char * skipref, const char skipcount);	 
	 char count(); //returns number of skips
	 char skipsize(const char index); //returns size of the specified skip
};
*/
 
 
class EncSkipGram: public EncAnyGram {
    public:
      char skipsize[MAXSKIPS]; //4 bytes reserved for skip size
      char skipcount; //number of skips
      
      //SkipConf skipconf; //skip configuration (encoded in 16 bits, 4 bit per gap: max 4 gaps, max 16 spaces in a gap)
      
      const char n() const;
      
      EncSkipGram(const std::vector<EncNGram*> & dataref, const std::vector<int> & skipref, bool initialskip = false, bool finalskip = false);
      EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount);
      //EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn);
      EncSkipGram(std::istream * in, const char gapcount = -1);
      std::string decode(ClassDecoder& classdecoder) const;
      bool out() const;
      bool isskipgram() const { return (skipcount > 0); }
      
      const char gapcount() const {
      		return skipcount;
         	//return skipconf.count();
      }
      const char gapsize(char i) const {
      		if (i >= MAXSKIPS) {
      			std::cerr << "ERROR: Request for size of gap " << i << " exceed maximum amount of gaps" << std::endl;
      			exit(14);
      		}
      		return skipsize[i]; 
            //return skipconf.skipsize(); 
      }
      
      int parts(std::vector<EncNGram*> & container) const; //returns all consecutive parts
      //int instantiate(const EncSkipGram * skipcontent, std::vector<EncSkipGram*> & container) const;
      //int instantiate(const EncSkipGram * skipcontent, std::vector<EncSkipGram*> & container, const std::vector<EncNGram*> & parts, const std::vector<EncNGram*> & contentparts) const; //returns all instances (as skipgrams cause instances may be skipgrams) 
	  
	  EncNGram instantiate(const EncSkipGram * skipcontent, const std::vector<EncNGram*> & contentparts) const;      
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
#endif
