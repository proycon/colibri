#include "ngram.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <algorithms.h>
#include <limits>

using namespace std;


//TODO: Computation of skipgramtokens[] is still bugged upon processing a corpus

EncAnyGram::EncAnyGram() {
    _size = 0;
    data = NULL;
}

EncAnyGram::EncAnyGram(const unsigned char* dataref, const char size) {
   //create a copy of the character data (will take less space than storing pointers anyhow!)
   if (size <= 0) {
       cerr << "EncAnyGram: Data size must be > 0! Parameter says " << (int) size << "!" << endl;
       exit(3);
   }
   _size = size;
   
   data = new unsigned char[size];
   for (int i = 0; i < size; i++) {
        data[i] = dataref[i];
   }
}

EncAnyGram::EncAnyGram(const EncAnyGram& ref) {
    _size = ref.size();
    if (_size <= 0) {
        cerr << "EncAnyGram: Data size must be > 0, reference n-gram has " << (int) _size << "!" << endl;
        exit(3);    
    }
    data = new unsigned char[_size];   
    for (int i = 0; i < _size; i++) {
        data[i] = ref.data[i];
    }    
}

EncAnyGram::~EncAnyGram() {     
    if (data != NULL) delete [] data;        
    data = NULL;
}




const char EncAnyGram::size() const {
    return _size;
}

const char EncAnyGram::n() const {
    char count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    return count;
}

const size_t EncAnyGram::hash() const {
    //adapted from jenkins hash: http://en.wikipedia.org/wiki/Jenkins_hash_function
    unsigned long h;
    int i;
    bool prevnull = false;
    int skipnum = 0;
    for(h = i = 0; i < size(); ++i)
    {                
        h += data[i];
        h += (h << 10);
        h ^= (h >> 6);
        if (data[i] == 0) {
            if (prevnull) {
                h += 46021 + gapsize(skipnum);                        
                h += (h << 10);
                h ^= (h >> 6);
                skipnum++;
            } else {
                prevnull = true;
            }                    
        } else {
            prevnull = false;
        }                
    }
    h += (h << 3);
    h ^= (h >> 11);
    h += (h << 15);
    return h;
}

EncNGram * EncNGram::slice(const int begin,const int length) const {    
    if (length <= 0) {
        cerr << "slice got length argument <= 0! Not possible!" << endl;
        exit(3);
    }
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
    if (bytesize <= 0) {
        cerr << "getencgram yields ngram with size <= 0! Not possible!" << " index="<<index << " n="<<n <<" size="<< size << endl;
        exit(3);
    }
    return new EncNGram(line + beginpos, bytesize);
}

int EncNGram::subngrams(vector<EncNGram*> & container) const {
    int count = 0;
    const int N = n();
    for (int begin = 0; begin < N; begin++) {
        for (int length = 1; length < N-begin; length++)
            if (length < N) {
                count++; 
                container.push_back( slice(begin,length) );
            }
    }        
    return count;
}



std::string EncAnyGram::decode(ClassDecoder& classdecoder) const {
    //cout << "DECODING NGRAM size=" << (int) _size << " n=" << n() << " data=" << data << endl;
    std::string result = ""; 
    int begin = 0;
    int l = 0;;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {             
            const unsigned int cls = bytestoint(data + begin, l);              
            if (cls == 1) {
                //cout << "EOL FOUND" << endl;
                return result;
            } else {  
                //cout << " CLASS " << cls << " (length " << l << ") DECODES TO " << classdecoder[cls] << endl;
                result += classdecoder[cls] + ' ';
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

bool EncAnyGram::out() const {
    //cout << "DECODING NGRAM size=" << (int) _size << " n=" << n() << " data=" << data << endl;
    int begin = 0;
    int l = 0;;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {             
            const unsigned int cls = bytestoint(data + begin, l);              
            if (cls == 1) {
                //cout << "EOL FOUND" << endl;
                return true;
            } else {  
                //cout << " CLASS " << cls << " (length " << l << ") DECODES TO " << classdecoder[cls] << endl;
                cout << cls << ' ';
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        cout << cls << ' ';
        //cout << "FINAL CLASS " << cls << " DECODES TO " << classdecoder[cls] << endl;
    }    
    return true;
}


bool EncAnyGram::operator==(const EncAnyGram &other) const {
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
bool EncAnyGram::operator!=(const EncAnyGram &other) const {
    return !(*this == other);
}

EncAnyGram & EncAnyGram::operator =(EncAnyGram other) { //(note: argument passed by value!
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


void EncAnyGram::writeasbinary(ostream * out) const {
    if (_size <= 0) {
        cerr << "Writing skipgram with size <= 0! Not possible!" << endl;
        exit(3);
    }
    out->write( &_size, sizeof(char) ); //data length
    out->write( (char*) data , (int) _size ); //data
}

void EncSkipGram::writeasbinary(ostream * out) const {
    out->write( &skipcount, sizeof(char) ); //nr of gaps
    for (int j = 0; j < skipcount; j++) { //skip configuration
            out->write( skipsize + j , sizeof(char) );
    }
    if (_size <= 0) {
        cerr << "Writing skipgram with size <= 0! Not possible!" << endl;
        exit(3);
    }
    out->write( &_size, sizeof(char) ); //size
    out->write( (char*) data , (int) _size ); //data
}



EncSkipGram::EncSkipGram(const vector<EncNGram*> & dataref, const vector<int> & skipref, bool initialskip, bool finalskip): EncAnyGram() {
    //compute size
    _size = 0;
    //cerr << "-- INITIAL: " << initialskip << endl;
    //cerr << "-- FINAL: " << finalskip << endl;
    for (char i = 0; i < (char) dataref.size(); i++) {
        _size += (dataref[i])->size();
        if (i < (char) dataref.size() - 1) _size += 2; //double 0 byte delimiting subngrams
        //cerr << "--SUBNGRAM-- n=" << (int) (dataref[i])->n() << ",size=" << (int) (dataref[i])->size()  << endl;        
        //for (int j = 0; j < (dataref[i])->size(); j++) cerr << (int) dataref[i]->data[j] << endl;
    }
    if (initialskip) {
        _size += 2; //two null bytes at start
    }
    if (finalskip) {
        _size += 2; //extra null byte at end
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
        }
        if (i < dataref.size() - 1) {            
            data[cursor++] = '\0';
            data[cursor++] = '\0';
        }                
    }    
    if (finalskip) {
        data[cursor++] = '\0';
        data[cursor++] = '\0';
    }       
    
    //sanity check
    skipcount = 0;
    bool prevnull = false;
    //cerr << "--SKIPGRAM-- (size=" << (int) _size << ")" << endl;
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == '\0') {
            if (prevnull) {
                prevnull = false;
                skipcount++;
            } else {
                prevnull = true;
            }
        } else {
            prevnull = false;
        }        
    }
    if ((char) skipref.size() != skipcount) {
        cerr << "ENCSKIPGRAM ERROR: Skipgram contains " << (int) skipcount << " skips, but configuration specifies " << skipref.size() << endl;      
        cerr << data <<endl;
        exit(1);
    }    
}

EncSkipGram::EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount): EncAnyGram(dataref,size) {
    this->skipcount = skipcount;
    for (int i = 0; i < skipcount; i++) {
        skipsize[i] = skipref[i];
    }
}


const char EncSkipGram::n() const {    
    char count = 0;
    bool item = (data[0] != 0);
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) { 
            if (item) count++;
            item = false;
        } else {
            item = true;
        }
    }    
    if (item) count++;
    for (int i = 0; i < skipcount; i++) {
        count += skipsize[i]; //minus two because each skip already counted for two in the previous loop
    }
    return count;    
}


EncSkipGram::EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn): EncAnyGram() {
    const char pregapsize = pregap.size();
    const char postgapsize = postgap.size();
    skipcount = 1;
    skipsize[0] = refn - pregap.n() - postgap.n();
            
    _size = pregapsize + postgapsize + 2;
    if (_size <= 0) {
        cerr << "EncSkipGram: Data size must be > 0!" << endl;
        exit(3);    
    }
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



std::string EncSkipGram::decode(ClassDecoder& classdecoder) const {
    std::string result = ""; 
    int begin = 0;
    int l = 0;
    int skipnum = 0;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {            
            if ((i > 0) && (data[i-1] == 0)) {
                char chr[2];
                chr[0] = 48 + skipsize[skipnum];
                chr[1] = '\0';                
                skipnum++;
                result += std::string("{*") + std::string(chr) + std::string("*} ");                
            } else {            
                const unsigned int cls = bytestoint(data + begin, l);              
                if (cls == 1) {
                    return result;
                } else if (cls > 0) {  
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
    }    
    return result;
}



bool EncSkipGram::out() const {
    int begin = 0;
    int l = 0;
    int skipnum = 0;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {            
            if ((i > 0) && (data[i-1] == 0)) {
                cout << "{*" << (int) skipsize[skipnum++] << "*} ";                
            } else {            
                const unsigned int cls = bytestoint(data + begin, l);              
                if (cls == 1) {
                    return true;
                } else if (cls > 0) {  
                    cout << cls << ' ';
                }
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        cout << cls;
    }    
    return true;
}


int EncSkipGram::parts(std::vector<EncNGram*> & container) const {
    int begin = 0;
    int count = 0;
    bool prevnull = false;
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == '\0') {
            if (prevnull) {                            
                prevnull = false;                
            } else {                
                prevnull = true;                
                container.push_back( new EncNGram(data + begin,i-begin) );
                begin = i+1;
                count++;
            }
        } else {
            prevnull = false;
        }        
    }
    if (begin < _size - 1) {
        container.push_back( new EncNGram(data + begin,_size-begin) );
    }
    return count;
}




/*size_t jenkinshash(unsigned char * data, char size) {
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
}*/




EncNGram::EncNGram(istream * in) {
    in->read(&_size, sizeof(char));    
    if (_size <= 0) {
        cerr << "EncNGram: data has to have size >0, " << (int) _size << " is not possible!" << endl;;
        exit(5);
    }
    data = new unsigned char[_size];
    in->read((char*) data, (int) _size); //read data                                                
}

EncSkipGram::EncSkipGram(istream * in, const char gapcount) {
    if (gapcount < 0) {
        in->read(&skipcount, sizeof(char));
    } else {
        skipcount = gapcount;
    }    
    for (int j = 0; j < skipcount; j++) {
        in->read(skipsize + j, sizeof(char)); //reads in skipref[j]                
    }
    in->read(&_size, sizeof(char));
    if (_size <= 0) {
        cerr << "EncSkipGram: data has to have size >0, read " << (int) _size << ", not possible!" << endl;;
        exit(5);
    }
    data = new unsigned char[_size];
    in->read((char*) data, (int) _size); //read data                                                
}


