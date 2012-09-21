#ifndef CLASSDECODER_H
#define CLASSDECODER_H

#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>

class ClassDecoder {
    private:
     std::unordered_map<unsigned int,std::string> classes;
     unsigned int unknownclass;
     unsigned int bosclass;
     unsigned int eosclass;
     unsigned int highestclass;
    public:
    
    ClassDecoder(const std::string & filename);
    
    std::vector<std::string> decodeseq(const std::vector<int> & seq);
    
    void decodefile(const std::string & filename, unsigned int start = 0, unsigned int end = 0);
    //std::string decodestring(const unsigned char * data, unsigned char datasize); 
    
    int size() const {
        return classes.size();
    }
    
    std::string operator[](unsigned int key) {
         return classes[key];
    }
    
    void add( unsigned int, std::string); 
    unsigned int gethighestclass() { return highestclass; }
    bool hasclass(unsigned int key) const { return (classes.count(key) > 0); } 
};

unsigned int bytestoint(const unsigned char* a, const int l);
int readline(std::istream* IN, unsigned char* buffer, const int);

const int countwords(const unsigned char* data, const int l);
std::pair<int,int> getwords(const unsigned char* data, const int datasize, const int n, const int begin = 0);
#endif
