#ifndef CLASSDECODER_H
#define CLASSDECODER_H

#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>

class ClassDecoder {
    private:
     std::unordered_map<unsigned int,std::string> classes;
     int unknownclass;
     int highestclass;
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
};

unsigned int bytestoint(const unsigned char* a, const int l);
int readline(std::istream* IN, unsigned char* buffer, const int);

const int countwords(const unsigned char* data, const int l);
#endif
