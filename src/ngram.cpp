#include "ngram.h"

EncNGram::EncNGram(unsigned char* dataref, char size) {
       //create a copy of the character data (will take less space than storing pointers anyhow!)
       _size = size;
       data = new unsigned char[size];
       for (int i = 0; i < size; i++) {
            data[i] = dataref[i];
       }
}

const char EncNGram::size() const {
    return _size;
}

const int EncNGram::n() const {
    int count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    return count;
}

EncNGram EncNGram::slice(const int begin,const int length) {
    return getencngram(begin, length, data, _size);
}


EncNGram getencngram(const int index, const int n, unsigned char *line, const int size) {
    int currentindex = 0;
    int beginpos = -1;
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
    return EncNGram(line + beginpos, bytesize); 
}

std::string EncNGram::decode(ClassDecoder& classdecoder) {
    std::string result = ""; 
    unsigned char* begin = data;
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) {
            //cout << "N: " << n << endl;
            const unsigned int cls = bytestoint(begin, i - 1);  
            if (cls == 1) {
                return result;
            } else {
                result += classdecoder[cls];
            }
            begin = data + i;
        }
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

