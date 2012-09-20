#include "ngram.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <algorithms.h>
#include <limits>

using namespace std;


EncAnyGram::EncAnyGram() {
    _size = 0;
    data = NULL;
}

EncAnyGram::EncAnyGram(const unsigned char* dataref, const char size) {
   //create a copy of the character data (will take less space than storing pointers anyhow!)
   if (size <= 0) {
       cerr << "INTERNAL ERROR EncAnyGram(): Data size must be > 0! Parameter says " << (int) size << "!" << endl;
       exit(13);
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
        cerr << "INTERNAL ERROR EncAnyGram(): Data size must be > 0, reference n-gram has " << (int) _size << " (n=" << (int) ref.n() << ") !" << endl;
        exit(13);    
    }
    data = new unsigned char[_size];   
    for (int i = 0; i < _size; i++) {
        data[i] = ref.data[i];
    }    
}

EncAnyGram::EncAnyGram(const EncData& ref) {
    _size = ref.size();
    if (_size <= 0) {
        cerr << "INTERNAL ERROR EncAnyGram(): Data size must be > 0, reference encdata has " << (int) _size << " (n=" << (int) ref.length() << ") !" << endl;
        exit(13);    
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
    if (isskipgram()) {
        return ((const EncSkipGram *) this)->n();
    } else {
        char count = 1; 
        for (int i = 0; i < _size; i++) {
            if (data[i] == 0) count++;
        }    
        return count;
    }
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
        cerr << "INTERNAL ERROR EncNGram::slice(): slice got length argument <= 0! Not possible!" << endl;
        exit(13);
    }
    return getencngram(begin, length, data, _size);
}

EncNGram * getencngram(const int index, const int n, const unsigned char *line, const int size, const unsigned int linenum) {
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
        cerr << "INTERNAL ERROR getencgram(): yielding ngram with size <= 0! Not possible!" << " index="<<index << " n="<<n <<" size="<< (int) bytesize << " beginpos=" << beginpos << " endpos=" << endpos << " sentencesize=" << size << endl;
        if (linenum > 0) cerr << "OCCURRED ON LINE " << linenum << endl;
        cerr << "ENCODED DATA DUMP FOR DEBUG (index:byte): ";
        for (int i = 0; i < size; i++) {
        	cerr << i << ':' << (int) line[i] << ' ';        	
        }
        cerr << endl;
        exit(13);
    }
    return new EncNGram(line + beginpos, bytesize);
}
/*
EncAnyGram * EncSkipGram::slice(const int start,const int length) const {
    //low level slice algorithm for skipgrams

    const unsigned char * newskipconf;
    unsigned char * buffer = new unsigned char[4048];
    unsigned char newsize = 0;   

    const unsigned char unknownclass = 2;
    bool prevnull = false;
    int skipnum = 0;
    int cursor = 0;
    bool capture = false;
    int begin = 0;    
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) {
            //if (capture) {
                 //return new EncNGram(data + begin,i-begin);                
            //}        
            if (prevnull) {
                if (skipsize[skipnum] == 0) throw Variablewidthexception();
                for (int j = 0; j < skipsize[skipnum]; j++) {                    
                    if (cursor == start) {
                        begin = i + 1;
                        buffer[newsize++] = 0;
                                             
                        //return new EncNGram(&unknownclass,1);
                    }
                    if (cursor==start+length) (
                        
                    }  
                    cursor++;
                }           
                skipnum++;              
            } else {
                if ((cursor >= start) && (cursor <= start+length)) {  
                    buffer[newsize++] = data[i];
                }            
                cursor++;
            }            
            prevnull = true;
        } else {
            prevnull = false;
            if ((cursor >= start) && (cursor <= start+length)) {            
                buffer[newsize++] = data[i];
            }         
        }        
    }
    if (capture) {
        return new EncNGram(data + begin,_size-begin);
    } else {
        cerr << "ERROR: EncSkipGram::gettoken(): index not found " << index << endl;
        exit(6);
        return NULL;        
    }  
}
*/

int EncNGram::subngrams(vector<EncNGram*> & container) const {
    int count = 0;
    const int N = n();
    //maybe TODO: container.reserve() ? 
    for (int begin = 0; begin < N; begin++) {
        for (int length = 1; length < N-begin+1; length++)
            if (length < N) {
                count++; 
                container.push_back( slice(begin,length) );
            }
    }        
    return count;
}

int EncNGram::splits(vector<pair<EncNGram*, EncNGram*> > & container) const {
	int count = 0;
	const int N = n();
	container.reserve(N-2);
    for (int begin = 1; begin < N - 1; begin++) {
    	int length = N-begin;
    	count++;
    	container.push_back( make_pair<EncNGram *,EncNGram *>(slice(0,begin),slice(begin,length) ) ); 
    }        	
    return count;
}

int EncNGram::getclass(const int index) const {
    int begin = 0;
    int l = 0;
    int curindex = 0;
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {                        
            if (curindex == index) return bytestoint(data + begin, l);
            curindex++;                          
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        if (curindex == index) return bytestoint(data + begin, l);
    }    
    return 0;    
} 

std::string EncAnyGram::decode(ClassDecoder& classdecoder) const {
    //cout << "DECODING NGRAM size=" << (int) _size << " n=" << n() << " data=" << data << endl;    
    if (isskipgram()) {
        return ((const EncSkipGram*) this)->decode(classdecoder);
    } else {
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
                    if (classdecoder.hasclass(cls)) {
                        result += classdecoder[cls] + ' ';
                    } else {
                        result += "{NOTFOUND!} ";       //should never happen             
                    }
                }
                begin = i + 1;            
                l = 0;
            }
        }
        if (l > 0) {
            const unsigned int cls = bytestoint(data + begin, l);
            if (classdecoder.hasclass(cls)) {  
                result += classdecoder[cls];            
            } else {
                result += "{NOTFOUND!} ";
            }
            //cout << "FINAL CLASS " << cls << " DECODES TO " << classdecoder[cls] << endl;
        }    
        return result;
    }
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
                cerr << cls << ' ';
            }
            begin = i + 1;           
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        cerr << cls << ' ';
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
        cerr << "INTERNAL ERROR EncNGram::writeasbinary(): Writing skipgram with size <= 0! Not possible!" << endl;
        exit(13);
    }
    out->write( &_size, sizeof(char) ); //data length
    out->write( (char*) data , (int) _size ); //data
}

void EncSkipGram::writeasbinary(ostream * out) const {
    out->write( &skipcount, sizeof(char) ); //nr of gaps
    if (skipcount > MAXSKIPS) {
    	cerr << "INTERNAL ERROR EncSkipGram::writeasbinary(): Writing skipgram with skipcount > MAXSKIPS! (" << (int) skipcount << " > " << (int) MAXSKIPS << "). Not possible!" << endl;
    	exit(13);
    }
    for (int j = 0; j < skipcount; j++) { //skip configuration
            out->write( skipsize + j , sizeof(char) );
    }
    if (_size <= 0) {
        cerr << "INTERNAL ERROR EncSkipGram::writeasbinary(): Writing skipgram with size <= 0! Not possible!" << endl;
        exit(13);
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
        cerr << "INTERNAL ERROR: EncSkipGram(): Skipgram contains " << (int) skipcount << " skips, but configuration specifies " << (int) skipref.size() << endl;      
        cerr << data <<endl;
        exit(13);
    } 
}

EncSkipGram::EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount): EncAnyGram(dataref,size) {
    this->skipcount = skipcount;
    for (int i = 0; i < skipcount; i++) {
        skipsize[i] = skipref[i];
    }
    
    //sanity check
    char foundskipcount = 0;
    bool prevnull = false;
    //cerr << "--SKIPGRAM-- (size=" << (int) _size << ")" << endl;
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == '\0') {
            if (prevnull) {
                prevnull = false;
                foundskipcount++;
            } else {
                prevnull = true;
            }
        } else {
            prevnull = false;
        }        
    }
    if (skipcount != foundskipcount) {
        cerr << "INTERNAL ERROR: EncSkipGram(): Skipgram contains " << (int) foundskipcount << " skips, but configuration specifies " << (int) skipcount << endl;      
        cerr << data <<endl;
        exit(13);
    }    
}


bool EncSkipGram::operator==(const EncSkipGram &other) const {
        const char othersize = other.size();
        if (_size == othersize) {
            if (skipcount != other.skipcount) return false;
            for (int i = 0; i < _size; i++) {
                if (data[i] != other.data[i]) return false;
            }
            for (int i = 0; i < skipcount; i++) {
               if (skipsize[i] != other.skipsize[i]) return false;
            }
            return true;            
        } else {
            return false;
        }        
}
bool EncSkipGram::operator!=(const EncSkipGram &other) const {
    return !(*this == other);
}

bool EncSkipGram::variablewidth() const {
    for (int i = 0; i < skipcount; i++) {
        if (skipsize[i] == 0) return true;
    }
    return false;
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
        if (skipsize[i] ==0) throw Variablewidthexception();
        count += skipsize[i]; //minus two because each skip already counted for two in the previous loop
    }
    return count;    
}

/*
EncSkipGram::EncSkipGram(const EncNGram & pregap, const EncNGram & postgap, const char refn): EncAnyGram() {
    const char pregapsize = pregap.size();
    const char postgapsize = postgap.size();
    skipcount = 1;
    skipsize[0] = refn - pregap.n() - postgap.n();
            
    _size = pregapsize + postgapsize + 2;
    if (_size <= 0) {
        cerr << "INTERNAL ERROR EncSkipGram(): Data size must be > 0! Got " << (int) _size << endl;
        exit(13);    
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
*/



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
                if (skipsize[skipnum] == 0) {
                    chr[0] = 'V';
                } else {
                    chr[0] = 48 + skipsize[skipnum];
                }
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
                cerr << "{*" << (int) skipsize[skipnum++] << "*} ";                
            } else {            
                const unsigned int cls = bytestoint(data + begin, l);              
                if (cls == 1) {
                    return true;
                } else if (cls > 0) {  
                    cerr << cls << ' ';
                }
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        cerr << cls;
    }    
    return true;
}


int EncSkipGram::parts(std::vector<EncNGram*> & container) const {
    int begin = 0;
    bool prevnull = false;
    container.reserve(skipcount);
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == 0) {
            if (prevnull) {                            
                if (i-begin - 1 > 0) container.push_back( new EncNGram(data + begin,i-begin - 1) );
                begin = i+1;                
            }
            prevnull = true;
        } else {
            prevnull = false;
        }        
    }
    if (_size - begin > 0) {        
        container.push_back( new EncNGram(data + begin,_size-begin) );        
    }
    return container.size();
}


EncNGram * EncSkipGram::gettoken(int index) const {
    const unsigned char unknownclass = 2;
    bool prevnull = false;
    int skipnum = 0;
    int cursor = 0;
    
    int begin = 0;
    bool capture = false;    
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == 0) {
            if ((capture) && (i-begin > 0)) {
                 return new EncNGram(data + begin,i-begin);                
            }        
            if (prevnull) {
                if (skipsize[skipnum] == 0) throw Variablewidthexception();
                for (int j = 0; j < skipsize[skipnum]; j++) {                    
                    if (cursor == index) return new EncNGram(&unknownclass,1); 
                    cursor++;
                }           
                skipnum++;              
            } else {                
                cursor++;              
            }
            if (cursor == index) {
                begin = i+1;
                capture = true;
            }  
            prevnull = true;
        } else {
            prevnull = false;            
        }
    }
    if (capture) {
        return new EncNGram(data + begin,_size-begin);
    } else {
        cerr << "ERROR: EncSkipGram::gettoken(): index not found " << index << endl;
        exit(6);
        return NULL;        
    }  
}

EncNGram * EncNGram::gettoken(int index) const {
    return slice(index,1);
}

bool EncAnyGram::unknown() { //does this anygram have an unknown class in it?
    unsigned char unknownclass= 2;
    for (int i = 0; i < _size; i++) {
        if ((i == 0) && (data[i] == unknownclass) && ((_size == 1) || (data[i+1] == 0))) {
            return true;
        } else if ((data[i - 1] == 0) && (data[i] == unknownclass) && ((_size == i) || (data[i+1] == 0))) {
            return true;
        }
    }  
    return false;
}


void EncSkipGram::mask(std::vector<bool> & container) const { //returns a boolean mask of the skipgram (0 = gap(encapsulation) , 1 = skipgram coverage)
    bool prevnull = false;
    int skipnum = 0;
    container.reserve(n());
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == 0) {
            if (prevnull) {       
                if (skipsize[skipnum] == 0) throw Variablewidthexception();
                for (int j = 0; j < skipsize[skipnum]; j++) {
                    container.push_back(false);
                }           
                skipnum++;              
            } else {
                container.push_back(true);
            }
            prevnull = true;
        } else {
            prevnull = false;            
        }
    }
    if (data[_size-1] != 0) container.push_back(true);
}


void EncSkipGram::getgaps(std::vector<std::pair<int,int> > & gaps) const {
    int pos = 0;
    bool prevnull = false;
    int skipnum = 0;
    gaps.reserve(skipcount);
    for (int i = 0; i < _size; i++) {
        //cerr << (int) data[i] << ':' << prevnull << ':' << skipcount << endl;
        if (data[i] == 0) {
            if (prevnull) {                  
                gaps.push_back(pair<int,int>( pos , skipsize[skipnum] ) );
                pos += skipsize[skipnum];
                skipnum++;              
            } else {
                pos++;
            }
            prevnull = true;        
        } else {
            prevnull = false;            
        }
        
    }  
}

void EncSkipGram::getparts(std::vector<std::pair<int,int> > & p) const {    
    std::vector<std::pair<int,int> > gaps;
    getgaps(gaps);
    
    int beginpart = 0;
    for (vector<pair<int,int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) {
        if (iter->first - beginpart > 0) {
            p.push_back(pair<int,int>(beginpart, iter->first - beginpart  ) );
            beginpart = iter->first + iter->second;
        }        
    } 
    if (beginpart != n()) {
        p.push_back(pair<int,int>(beginpart, n() - beginpart) );
    }    
}




EncSkipGram EncSkipGram::extractskipcontent(EncNGram & instance) const {
    if (instance.n() != n()) {
        cerr << "WARNING: Extractskipcontent(): instance.n() != skipgram.n(), " << (int) instance.n() << " != " << (int) n() << endl;
        cerr << "INSTANCE: " << instance.out() << endl;
        cerr << "SKIPGRAM: " << out() << endl;
    }
    std::vector<std::pair<int,int> > gaps;
    getgaps(gaps);
    vector<EncNGram*> subngrams;
    vector<int> skipcontent_skipref;
    int cursor = 0;
    for (std::vector<std::pair<int,int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) {
        EncNGram * subngram = instance.slice(iter->first,iter->second);
        subngrams.push_back(subngram);
        if (cursor > 0) skipcontent_skipref.push_back(iter->first - cursor);
        cursor = iter->first + iter->second;
    }
    EncSkipGram sc = EncSkipGram(subngrams, skipcontent_skipref, false, false);
    for (std::vector<EncNGram*>::iterator iter = subngrams.begin(); iter != subngrams.end(); iter++) {
        delete *iter;
    } 
    return sc;
}

/*int instantiate(const EncSkipGram * skipcontent,std::vector<EncSkipGram*> & container) const {
	skipcontent->parts(contentparts); //parts
	parts(p); //parts
	int r = instantiate(skipcontent,container,p, contentparts);
	for (vector<EncNGram*>::iterator iter = p.begin(); iter != p.end(); iter++) delete *iter;
	for (vector<EncNGram*>::iterator iter = contentparts.begin(); iter != contentparts.end(); iter++) delete *iter;
	return r;
}

int instantiate(const EncSkipGram * skipcontent,std::vector<EncSkipGram*> & container, const std::vector<EncNGram*> & p, const std::vector<EncNGram*> & contentparts) const {*/

EncNGram EncSkipGram::instantiate(const EncSkipGram * skipcontent) const {
	vector<EncNGram*> contentparts;
	skipcontent->parts(contentparts); 
	return instantiate(skipcontent,contentparts);
}

EncNGram EncSkipGram::instantiate(const EncSkipGram * skipcontent, const std::vector<EncNGram*> & contentparts) const {

	unsigned char buffer[2048]; 
    int l = 0;
    int skipnum = 0;
	int buffercursor = 0;
	if ((size_t) skipcount != contentparts.size()) {
		cerr << "FATAL ERROR: content parts should be equal to skipcount! " <<  contentparts.size() << " content parts, " << (int) skipcount << " skipcount" << endl;
		exit(13);
	}
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {            
            if ((i > 0) && (data[i-1] == 0)) {
		        if ((size_t) skipnum >= contentparts.size()) {
					cerr << "FATAL ERROR: not enough content parts for instantiation! " <<  contentparts.size() << " content parts, i=" << i << endl;
					cerr << "DEBUG: skipgram out:" << endl;
					out();
					cerr << "DEBUG: skipcontent out:" << endl;
					skipcontent->out();
					exit(13);
				} 
				const EncNGram * ngram = contentparts[skipnum];
				for (int j = 0; j < ngram->size(); j++) {
						buffer[buffercursor++] = ngram->data[j];
				}
				skipnum++;
            }            
            //begin = i + 1;            
            l = 0;
       }
       buffer[buffercursor++] = data[i];        
    }

	EncNGram x = EncNGram(buffer, buffercursor); 
	return x;
}


bool EncSkipGram::classvector(vector<int> & container) const {
    int begin = 0;
    int l = 0;
    int skipnum = 0;
    container.reserve(n());
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {            
            if ((i > 0) && (data[i-1] == 0)) {
                for (int j = 0; j < skipsize[skipnum]; j++) container.push_back(0);
                skipnum++;       
            } else {            
                const unsigned int cls = bytestoint(data + begin, l);              
                if (cls == 1) {
                    return true;
                } else if (cls > 0) {
                    container.push_back(cls);  
                }
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        container.push_back(cls);  
    }    
    return true;
}

bool EncNGram::classvector(vector<int> & container) const {
    int begin = 0;
    int l = 0;
    container.reserve(n());
    for (int i = 0; i < _size; i++) {
        l++;
        if ((data[i] == 0) && (l > 0)) {              
            const unsigned int cls = bytestoint(data + begin, l);              
            if (cls == 1) {
                return true;
            } else if (cls > 0) {
                container.push_back(cls);  
            }        
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        container.push_back(cls);  
    }    
    return true;
}


int EncSkipGram::instancetemplaterelation(const EncSkipGram *other) const {
    if (n() != other->n()) return 0;
    vector<int> thisvector; 
    this->classvector(thisvector);
    vector<int> othervector;
    other->classvector(othervector);
    int result = 0;    
    for (int i = 0; i < n(); i++) {
        if (thisvector[i] == othervector[i]) {
            //good.. pass
        } else if (thisvector[i] == 0) {
            if (result == 1) return 0; //mismatch
            result = -1; //this one is the template, other instance
        } else if (othervector[i] == 0) {
            if (result == -1) return 0; //mismatch
            result = 1; //this one is the instance, other template
        } else {
            //mismatch
            return 0;
        }
    }    
    return result;
}	

/*
SkipConf::SkipConf( const  uint16_t value ) {
	this->value = value;
}

SkipConf::SkipConf( const unsigned char * skipref, const char skipcount) {
	value = 0;
	if (skipcount > MAXSKIPS) {
		cerr << "ERROR: Too many skips (" << skipcount << "), the maximum is " << MAXSKIPS << endl;
	}
	for (int i = 0; i < skipcount; i++) {
		 if (skipref[i] >= MAXSKIPSIZE) {
		 		cerr << "ERROR: Skip " <<  i+1 << " is oversized: " << (int) skipref[i] << "), the maximum is " << MAXSKIPSIZE << endl;
		 }
		 value = value | skipref[i] << i*4;
	
}

char SkipConf:skipsize( const char index) {
	char size = (value >> index*4) & 0xf; //shift and mask off
	return (char) size; 
}


char SkipConf::count() {
	if (value == 0) return 0;
	char skipcount = 1;
	for (int i = 1; i < MAXSKIPS; i++) {
		if (value > pow(2,i*4)) {
			skipcount++;
		} else {
			break;
		}		
	}
	return skipcount;
}
*/


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
        cerr << "INTERNAL ERROR: EncNGram(): data has to have size >0, " << (int) _size << " is not possible!" << endl;;
        exit(13);
    }
    data = new unsigned char[_size];
    in->read((char*) data, (int) _size); //read data                                                
}

EncSkipGram::EncSkipGram(istream * in, const char _skipcount) {
    if (_skipcount < 0) {
        in->read(&skipcount, sizeof(char));
    } else {
        skipcount = _skipcount;
    }    
     if (skipcount > MAXSKIPS) {
    	cerr << "INTERNAL ERROR EncSkipGram::EncSkipGram(): Reading skipgram with skipcount > MAXSKIPS! (" << (int) skipcount << " > " << (int) MAXSKIPS << "). Not possible!" << endl;
    	exit(13);
    }
    for (int j = 0; j < skipcount; j++) {
        in->read(skipsize + j, sizeof(char)); //reads in skipref[j]                
    }
    in->read(&_size, sizeof(char));
    if (_size <= 0) {
        cerr << "INTERNAL ERROR: EncSkipGram(): data has to have size >0, read " << (int) _size << ", not possible!" << endl;;
        cerr << "skipcount=" << (int) skipcount << endl;
        exit(13);
    }
    data = new unsigned char[_size];
    in->read((char*) data, (int) _size); //read data                                                
}





EncData::EncData(const unsigned char* dataref, const int size) {
   //create a copy of the character data (will take less space than storing pointers anyhow!)
   if (size <= 0) {
       cerr << "INTERNAL ERROR EncData(): Data size must be > 0! Parameter says " << (int) size << "!" << endl;
       exit(13);
   }
   _size = size;
   
   data = new unsigned char[size];
   for (int i = 0; i < size; i++) {
        data[i] = dataref[i];
   }
}

EncData::EncData(const EncData& ref) {
    _size = ref.size();
    if (_size <= 0) {
        cerr << "INTERNAL ERROR EncAnyGram(): Data size must be > 0, reference data has " << (int) _size << " length=" << (int) ref.length() << ") !" << endl;
        exit(13);    
    }
    data = new unsigned char[_size];   
    for (int i = 0; i < _size; i++) {
        data[i] = ref.data[i];
    }    
}

EncData::~EncData() {     
    if (data != NULL) delete [] data;        
    data = NULL;
}

const int EncData::length() const {
    char count = 1; 
    for (int i = 0; i < _size; i++) {
        if (data[i] == 0) count++;
    }    
    return count;
}

EncNGram * EncData::slice(const int begin,const int length) const {    
    if (length <= 0) {
        cerr << "INTERNAL ERROR EncNGram::slice(): slice got length argument <= 0! Not possible!" << endl;
        exit(13);
    }
    return getencngram(begin, length, data, _size);
}

bool EncData::match(const EncNGram * ngram, const int offset) {
    if (offset + ngram->n() > length()) return false;
    const EncNGram * testpattern = slice(offset, ngram->n() );
    const int s = testpattern->size();
    if (ngram->size() == s) {
        for (int i = 0; i < s; i++) {
            if (ngram->data[i] != testpattern->data[i]) {
                delete testpattern;
                return false;
            }
        }
        delete testpattern;
        return true;
    } else {
        delete testpattern;
        return false; //no match
    }   
}

bool EncData::match(const EncSkipGram * skipgram, const int offset) {
    if (offset + skipgram->n() > length()) return false;
    vector<pair<int,int> > gaps;
    vector<EncNGram *> parts;
    skipgram->getgaps(gaps);
    skipgram->parts(parts);    
    
    int begin = 0;
    int partindex = 0;
    bool result = true;
    for (vector<pair<int,int> >::iterator iter = gaps.begin(); iter != gaps.end(); iter++) {
        const int gapbegin = iter->first;
        const int gapsize = iter->second;        
        if (gapbegin > begin) {
            if ((size_t) partindex >= parts.size()) {
                cerr << "INTERNAL ERROR: EncData::match(), partindex >= parts" << endl;
                exit(6); 
            }
            if (!match(parts[partindex], begin)) {
                result = false;
                break;
            }
            
            //prepare for next round 
            partindex++;
            begin = gapbegin + gapsize; 
        }
    }     
    if ((size_t) partindex < parts.size()) {
        if (!match(parts[partindex], begin)) {
            result = false;
        }
    }
    for (vector<EncNGram *>::iterator iter = parts.begin(); iter != parts.end(); iter++) {
        delete *iter;
    }
    return result;  
}


bool EncData::out() const {
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
                cerr << cls << ' ';
            }
            begin = i + 1;           
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);  
        cerr << cls << ' ';
        //cout << "FINAL CLASS " << cls << " DECODES TO " << classdecoder[cls] << endl;
    }    
    return true;
}


EncNGram EncNGram::operator +(const EncNGram& other) const {
    if (_size + other.size() > 255) {
        cerr << "ERROR: EncNGram::operator+ .. n-grams too big" << endl; 
        exit(6);
    } 
    unsigned char buffer[512];
    unsigned char offset;
    unsigned char newsize = _size + other.size(); 
    for (int i = 0; i < _size; i++) {
        buffer[i] = data[i];
    }
    offset = _size;
    if (buffer[_size -1] != 0) {
        newsize++;
        offset++;
        buffer[(int) _size] = 0;
    }    
    for (int i = 0; i < other.size(); i++) {
        buffer[offset+i] = other.data[i];
    }
    return EncNGram(buffer, newsize);
}


std::string EncData::decode(ClassDecoder& classdecoder) const {
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
                if (classdecoder.hasclass(cls)) {
                    result += classdecoder[cls] + ' ';
                } else {
                    result += "{NOTFOUND!} ";       //should never happen             
                }
            }
            begin = i + 1;            
            l = 0;
        }
    }
    if (l > 0) {
        const unsigned int cls = bytestoint(data + begin, l);
        if (classdecoder.hasclass(cls)) {  
            result += classdecoder[cls];            
        } else {
            result += "{NOTFOUND!} ";
        }
        //cout << "FINAL CLASS " << cls << " DECODES TO " << classdecoder[cls] << endl;
    }    
    return result;
}
