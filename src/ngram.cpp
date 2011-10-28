#include "ngram.h"
#include <iostream>

using namespace std;

EncNGram::EncNGram() {
    _size = 0;
    data = NULL;
}

EncNGram::EncNGram(const unsigned char* dataref, const char size) {
   //create a copy of the character data (will take less space than storing pointers anyhow!)
   _size = size;
   data = new unsigned char[size];
   for (int i = 0; i < size; i++) {
        data[i] = dataref[i];
   }
}

EncNGram::EncNGram(const EncNGram& ref) {
    _size = ref.size();
    data = new unsigned char[_size];   
    for (int i = 0; i < _size; i++) {
        data[i] = ref.data[i];
    }    
}

EncNGram::~EncNGram() {     
    if (data != NULL) delete [] data;        
    data = NULL;
}

const char EncNGram::size() const {
    return _size;
}

const char EncNGram::n() const {
    char count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    return count;
}

EncNGram * EncNGram::slice(const int begin,const int length) const {    
    return getencngram(begin, length, data, _size);
}

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size) {
    int currentindex = 0;
    int beginpos = 0;
    int endpos = -1;
    for (int i = 0; i < size; i++) {
        if (line[i] == 0) {
            currentindex++;
            if (currentindex == index) {
                beginpos = i+1;
            } else if (currentindex == index + n) {
                endpos = i - 1;
            }
        }        
    }
    if (endpos == -1) {
        endpos = size - 1;
    }
    const char bytesize = (char) (endpos - beginpos + 1);    
    return new EncNGram(line + beginpos, bytesize);
}

std::string EncNGram::decode(ClassDecoder& classdecoder) const {
    //cout << "DECODING NGRAM size=" << (int) _size << " n=" << n() << " data=" << data << endl;
    std::string result = ""; 
    int begin = 0;
    int l = 0;;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {            
            //cout << "N: " << n << endl;        
            if ((i > 0) && (data[i-1] == 0)) {
                //two 0 bytes in a row, indicates a gap:
                result += "* ";
            } else {            
                const unsigned int cls = bytestoint(data + begin, l);              
                if (cls == 1) {
                    //cout << "EOL FOUND" << endl;
                    return result;
                } else {  
                    //cout << " CLASS " << cls << " (length " << l << ") DECODES TO " << classdecoder[cls] << endl;
                    result += classdecoder[cls] + ' ';
                }
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        result += classdecoder[cls];
        //cout << "FINAL CLASS " << cls << " DECODES TO " << classdecoder[cls] << endl;
    }    
    return result;
}

bool EncNGram::operator==(const EncNGram &other) const {
        const char othersize = other.size();
        if (_size == othersize) {
            for (int i = 0; i < _size; i++) {
                if (data[i] != other.data[i]) return false;
            }
            return true;
        } else {
            return false;
        }        
}
bool EncNGram::operator!=(const EncNGram &other) const {
    return !(*this == other);
}

EncNGram & EncNGram::operator =(EncNGram other) { //(note: argument passed by value!
        //delete old data
        if (data != NULL) delete [] data;
        
        //set new data
        _size = other.size();        
        data = new unsigned char[_size];   
        for (int i = 0; i < _size; i++) {
            data[i] = other.data[i];
        }  
 
        // by convention, always return *this (for chaining)
        return *this;
}

EncSingleSkipGram::EncSingleSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn): EncNGram() {
    const char pregapsize = pregap.size();
    const char postgapsize = postgap.size();
    _size = pregapsize + postgapsize + 2;
    data = new unsigned char[_size];    
    int cursor = 0;
    for (int i = 0; i < pregapsize; i++) {
        data[cursor++] = pregap.data[i];
    }
    data[cursor++] = '\0'; //double \0 byte indicates gap
    data[cursor++] = '\0'; //double \0 byte indicates gap
    for (int i = 0; i < postgapsize; i++) {
        data[cursor++] = postgap.data[i];
    }        
    _n = refn;
}



EncSkipGram::EncSkipGram(const vector<EncNGram*> & dataref, const vector<int> & skipref, bool initialskip, bool finalskip): EncNGram() {
    //compute size
    _size = 0;
    skipcount = 0;
    for (char i = 0; i < (char) dataref.size(); i++) {
        _size += (dataref[i])->size();
        if (i < (char) dataref.size() - 1) skipcount++;
    }
    if (initialskip) {
        _size += 2; //two null bytes at start
        skipcount += 1;
    }
    if (finalskip) {
        _size += 2; //extra null byte at end
        skipcount += 1;
    }
    
    
    
    if ((char) skipref.size() != skipcount) {
        cerr << "ENCSKIPGRAM ERROR: Expected " << skipcount << " skips, but configuration specifies " << skipref.size() << endl;
        exit(1);
    } else if (skipcount > MAXSKIPS) {
        cerr << "ENCSKIPGRAM ERROR: Too many skips for memory! " << skipcount << endl;
        exit(1);
    }
                    
    data = new unsigned char[_size]; 
    
    for (unsigned int i = 0; i < skipref.size(); i++) {
        skipsize[i] = (char) skipref[i];
    }
    
    //good, now fill the data buffer
    int cursor = 0;
    if (initialskip) {
        data[cursor++] = '\0';
        data[cursor++] = '\0';
    }
    for (unsigned int i = 0; i < dataref.size(); i++) {
        for (int j = 0; j < dataref[i]->size(); j++) {
            data[cursor++] = dataref[i]->data[j];
            data[cursor++] = '\0';
        }
        if (i < dataref.size() - 1) {            
            data[cursor++] = '\0';
        }                
    }    
    if (finalskip) {
        data[cursor++] = '\0';
        data[cursor++] = '\0';
    }        
}

const char EncSkipGram::n() const {
    char count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    for (int i = 0; i < skipcount; i++) {
        count += skipsize[i] - 1; //minus one because each skip already counted for one in the previous loop
    }
    return count;    
}


EncSkipGram::EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn): EncNGram() {
    const char pregapsize = pregap.size();
    const char postgapsize = postgap.size();
    skipcount = 1;
    skipsize[0] = refn - pregap.n() - postgap.n();
            
    _size = pregapsize + postgapsize + 2;
    data = new unsigned char[_size];    
    int cursor = 0;
    for (int i = 0; i < pregapsize; i++) {
        data[cursor++] = pregap.data[i];
    }
    data[cursor++] = '\0'; //double \0 byte indicates gap
    data[cursor++] = '\0'; //double \0 byte indicates gap
    for (int i = 0; i < postgapsize; i++) {
        data[cursor++] = postgap.data[i];
    }        
}

size_t jenkinshash(unsigned char * data, char size) {
    //jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
    unsigned long h;
    int i;
    for(h = i = 0; i < size; ++i)
    {
        h += data[i];
        h += (h << 10);
        h ^= (h >> 6);
    }
    h += (h << 3);
    h ^= (h >> 11);
    h += (h << 15);
    return h;
}

