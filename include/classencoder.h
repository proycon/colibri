#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>

class ClassEncoder {
    private:
     std::unordered_map<std::string,unsigned int> classes;
    public:
    ClassEncoder();
    ClassEncoder(const std::string &); //load an existing classer
    void build(const std::string & filename); //build a class from this dataset
    
    
    std::vector<unsigned int> encodeseq(const std::vector<std::string> & seq);
    void encodefile(const std::string &, const std::string &);
    void save(const std::string & filename);
    
    int size() const {
        return classes.size();
    }
    
    unsigned int operator[](const std::string & key) {
         return classes[key];
    }
};    

bool validclass(unsigned int cls);
unsigned char * inttobytes(unsigned int, int & length);
int readline(std::istream* IN, unsigned char* buffer, const int);

const int countwords(const unsigned char* data, const int l);

