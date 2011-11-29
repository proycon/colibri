#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>

class ClassDecoder {
    private:
     std::unordered_map<unsigned int,std::string> classes;
    public:
    
    ClassDecoder(const std::string filename);
    
    std::vector<std::string> decodeseq(std::vector<int> seq);
    
    void decodefile(const std::string filename);
    
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
