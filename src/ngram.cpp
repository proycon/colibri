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





EncGramModel::EncGramModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  int MINSKIPTYPES, bool DOINDEX, bool DOREVERSEINDEX, bool DOSKIPCONTENT, bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
    this->DOINDEX = DOINDEX;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;
    this->DOINITIALONLYSKIP = DOINITIALONLYSKIP;
    this->DOFINALONLYSKIP = DOFINALONLYSKIP;
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;

    const int BUFFERSIZE = 65536;
    unsigned char line[BUFFERSIZE];

    ngrams.push_back(freqlist());
    skipgrams.push_back(skipgrammap());   

    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;
        ngrams.push_back(freqlist());
        skipgrams.push_back(skipgrammap());            
        
        int linenum = 0;
    
        tokencount[n] = 0;
        skiptokencount[n] = 0;
                
        vector< vector< pair<int,int> > > gaps;
        compute_multi_skips(gaps, vector<pair<int,int> >(), n);    
        
        
        ifstream *IN =  new ifstream( corpusfile.c_str() );    
        vector<unsigned int> words;
        while (IN->good()) {
            const int linesize = readline(IN, line, BUFFERSIZE );            
                    
            linenum++;

            if (linenum % 10000 == 0) {
                cerr << "\t@" << linenum << endl;
            }
                            
            
            const int l = countwords(line, linesize);            
            
            if (linesize > 0) //no { on purpose! applies to next for loop
            for (int i = 0; i < l - n + 1; i++) {
                
                EncNGram * ngram = getencngram(i,n, line, linesize);  
                //cout << "NGRAM("<<ngram->n()<<","<<(int)ngram->size() << ")" << endl;
                              
                if (n > 2) {                    
                    EncNGram * subngram1 = ngram->slice(0, n - 1);
                    if (!(ngrams[n-1].count(*subngram1))) {
                        delete subngram1;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram1;

                                                             
                    EncNGram * subngram2 = ngram->slice(1, n - 1);
                    if (!(ngrams[n-1].count(*subngram2))) {
                        delete subngram2;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram2;                    
                }
                
                
                ngrams[n][*ngram] += 1;            
                tokencount[n]++;
            

                if (DOINDEX) ngram_index[*ngram].insert(linenum);

                if (DOSKIPGRAMS) {
                    for (size_t j = 0; j < gaps.size(); j++) {

                        if (gaps[j].size() == 1) {
                           if (!DOINITIALONLYSKIP) {
                                if (gaps[j][0].first == 0) continue; 
                           }
                           if (!DOFINALONLYSKIP) {
                               if (gaps[j][0].first + gaps[j][0].second == n) continue; 
                           }
                        }

                        
                        vector<EncNGram*> subngrams;
                        vector<int> skipref;
                        bool initialskip = false;
                        bool finalskip = false;
                        int cursor = 0;
                        bool docount = true;    
                        int oc = 0;             
                        
                           
                        //cerr << "INSTANCE SIZE: " << gaps[j].size() << endl;
                        for (size_t k = 0; k < gaps[j].size(); k++) {                                                        
                            const int begin = gaps[j][k].first;  
                            const int length = gaps[j][k].second;                        
                            //cerr << begin << ';' << length << ';' << n << endl;
                            skipref.push_back( length); 
                            if (k == 0) {
                                initialskip = (begin == 0);
                            }                            
                            if (begin > cursor) {
                                EncNGram * subngram = ngram->slice(cursor,begin-cursor);
                                subngrams.push_back(subngram);                                                               
                                oc = ngrams[subngram->n()].count(*subngram);
                                if (oc) oc = ngrams[subngram->n()][*subngram];
                                if ((oc == 0) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) )    {
                                    docount = false;
                                    break;
                                }
                            }
                            cursor = begin + length;
                        }   
                        if (cursor < n) {
                            EncNGram * subngram = ngram->slice(cursor,n-cursor);
                            subngrams.push_back(subngram);
                            oc = ngrams[subngram->n()].count(*subngram);
                            if (oc) oc = ngrams[subngram->n()][*subngram];
                            if ((oc == 0) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) )    {
                                docount = false;
                                break;
                            }
                        } else {
                            finalskip = true;
                        }
                        if (initialskip && finalskip && skipref.size() <= 1) docount = false; //the whole n-gram is a skip, discard
                        if (docount) {
                            EncSkipGram skipgram = EncSkipGram(subngrams, skipref, initialskip, finalskip);                                                    
                            vector<EncNGram*> skipcontent_subngrams;
                            vector<int> skipcontent_skipref;
                            cursor = 0;
                            for (size_t k = 0; k < gaps[j].size(); k++) {
                                const int begin = gaps[j][k].first;  
                                const int length = gaps[j][k].second;
                                EncNGram * subskip = ngram->slice(begin,length);                                
                                skipcontent_subngrams.push_back(subskip);
                                if (cursor > 0) skipcontent_skipref.push_back(begin - cursor);
                                cursor = begin+length;
                            }   
                            EncSkipGram skipcontent = EncSkipGram(skipcontent_subngrams, skipcontent_skipref, false, false);                                                        
                            skipgrams[n][skipgram].count++;                                                            
                            skipgrams[n][skipgram].skips[skipcontent] += 1;
                            skiptokencount[n]++;                            
                            for (size_t k = 0; k < skipcontent_subngrams.size(); k++) {       
                                delete skipcontent_subngrams[k];
                            }                            
                            if (DOINDEX) skipgram_index[skipgram].insert(linenum);
                        }
                        //cleanup
                        for (size_t k = 0; k < subngrams.size(); k++) {
                            delete subngrams[k];
                        }
                    }                    
                }
                
                delete ngram;                 
            }            
            
        };

       cerr << "Found " << ngrams[n].size() << " " << n << "-grams (" << tokencount[n] << " tokens)";
       if (DOSKIPGRAMS) {
        cerr << " and " << skipgrams[n].size() << " skipgrams (" << skiptokencount[n] << " tokens)" << endl;
       } else {
        cerr << endl;
       }
    

       //prune n-grams
       int pruned = 0;
       for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
            //if (DOINDEX) iter2++;
            if (iter->second < MINTOKENS) {
                if (DOINDEX) ngram_index.erase( iter->first);        
                tokencount[n] -= iter->second;
                pruned++;
                ngrams[n].erase(iter->first);                        
            } else {
                ngramtokencount += iter->second;
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << ngrams[n].size() <<  " left (" << tokencount[n] << " tokens)" << endl;
    
       
       if (DOSKIPGRAMS) {       
           //prune skipgrams
           pruned = 0;
           for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {                
                bool pruneskipgram = false;
                if ((iter->second.count < MINTOKENS) || ((DOSKIPCONTENT && iter->second.skips.size() < MINSKIPTYPES)))  {
                    pruneskipgram = true;
                } else if (DOSKIPCONTENT) {                
                    int prunedskiptokens = 0;
                    for(skipgram_freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                        if (iter2->second < MINSKIPTOKENS) {
                            //prune skip
                            prunedskiptokens += iter2->second;
                            iter->second.skips.erase(iter2->first); 
                        }
                    }
                    if ( (prunedskiptokens > 0) && ( (iter->second.skips.size() < MINSKIPTYPES) || (iter->second.count - prunedskiptokens < MINTOKENS) ) ) { //reevaluate
                        pruneskipgram = true;
                    } else {
                        skipgramtokencount += iter->second.count;
                        skiptokencount[n] -= prunedskiptokens;
                    }
                }
                if (pruneskipgram) {
                    if (DOINDEX) skipgram_index.erase(iter->first);
                    skiptokencount[n] -= iter->second.count;
                    pruned++;
                    skipgrams[n].erase(iter->first);
                }
                        
           }
           cerr << "Pruned " << pruned << " skipgrams, " << skipgrams[n].size() <<  " left (" << skiptokencount[n] << " tokens)" << endl;
           
        }
    
        ngramtypecount += ngrams[n].size();
        skipgramtypecount += skipgrams[n].size();     
    }

    if (DOREVERSEINDEX) {
        //TODO: Compute reverse index        
    }
        
}


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
        cerr << "EncSkipGram: data has to have size >0, read " << _size << ", not possible!" << endl;;
        exit(5);
    }
    data = new unsigned char[_size];
    in->read((char*) data, (int) _size); //read data                                                
}

EncGramModel::EncGramModel(const string & filename, bool DOINDEX, bool DOREVERSEINDEX, bool DOSKIPCONTENT) {
    const bool DEBUG = false;
    this->DOINDEX = DOINDEX;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;    
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    ngrams.push_back(freqlist());
    skipgrams.push_back(skipgrammap());    
    
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if (!f) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    unsigned long totaltokens;
    f.read( (char*) &totaltokens, sizeof(unsigned long));        
    unsigned long totaltypes;
    f.read( (char*) &totaltypes, sizeof(unsigned long));            
    
    for (int i = 0; i < totaltypes; i++) {           
        char gapcount;
        /*unsigned char check;
        f.read((char*) &check, sizeof(char));
        if (check != 255) {
            cerr << "Error during read of item " << i + 1 << " , expected check 255, got " << (int) check << endl;            
            exit(2);
        }*/
        if (DEBUG) cerr << "\t@" << i;
        f.read(&gapcount, sizeof(char));
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            ngramtypecount++;
            //NGRAM
            EncNGram ngram = EncNGram(&f); //read from file                        
            if (ngram.n() > MAXLENGTH) {
                tokencount[ngram.n()] = 0;
                skiptokencount[ngram.n()] = 0;
                for (int j = 0; j < ngram.n() - MAXLENGTH; j++) {
                    ngrams.push_back(freqlist());
                    skipgrams.push_back(skipgrammap());
                }
                MAXLENGTH = ngram.n();
            }
            int count;
            f.read((char*) &count, sizeof(int)); //read occurrence count
            ngrams[ngram.n()][ngram] = count; //assign count 
            ngramtokencount += count;
            tokencount[ngram.n()] += count;
            int indexcount;                        
            f.read((char*) &indexcount, sizeof(int));
            for (int j = 0; j < indexcount; j++) {
                int index;
                f.read((char*) &index, sizeof(int));                
                if (DOINDEX) {
                    ngram_index[ngram].insert(index);         
                }
                if (DOREVERSEINDEX) {
                    bool found = false;
                    for (int k = 0; k < ngram_reverse_index[index].size(); k++) if (ngram_reverse_index[index][k] == ngram) { found = true; break; };
                    if (!found) ngram_reverse_index[index].push_back(ngram);
                }
                       
            }
        } else {
            //SKIPGRAM
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            skipgramtypecount++;            
            EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            int count;
            f.read((char*) &count, sizeof(int)); //read occurrence count            
            skipgrams[skipgram.n()][skipgram].count = count; //assign
            skipgramtokencount += count;
            skiptokencount[skipgram.n()] += count;            
            int skipcontentcount;
            f.read((char*) &skipcontentcount, sizeof(int));   
            if (DEBUG) cerr << "\tcontentcount=" << (int) skipcontentcount;
            for (int j = 0; j < skipcontentcount; j++) {                                
                EncSkipGram skipcontent = EncSkipGram(&f);  //also when !DOSKIPCONTENT, bytes have to be read
                f.read((char*) &count, sizeof(int)); //read occurrence count                
                if (DOSKIPCONTENT) {
                    skipgrams[skipgram.n()][skipgram].skips[skipcontent] = count; //skipcontent occurrence 
                }
            }
            int indexcount;
            f.read((char*) &indexcount, sizeof(int));        
            for (int j = 0; j < indexcount; j++) {
                int index;
                f.read((char*) &index, sizeof(int));
                if (DOINDEX) {
                    skipgram_index[skipgram].insert(index);
                }
                if (DOREVERSEINDEX) {
                    bool found = false;
                    for (int k = 0; k < skipgram_reverse_index[index].size(); k++) if (skipgram_reverse_index[index][k] == skipgram) { found = true; break; };
                    if (!found) skipgram_reverse_index[index].push_back(skipgram);                    
                }
            }
        }
        if (DEBUG)  cerr << endl;      //DEBUG  
    }
    f.close();
}



void EncGramModel::save(const std::string & filename) {
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    
    const char czero = 0;
    const int zero = 0;
    //const unsigned char check = 255;
    const unsigned long totaltokens = tokens();
    const unsigned long totaltypes = types();
    
    f.write( (char*) &totaltokens, sizeof(unsigned long) );
    f.write( (char*) &totaltypes, sizeof(unsigned long) );
    for (int n = 1; n <= MAXLENGTH; n++) {
        for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {        
            //f.write((char*) &check, sizeof(char)); 
            f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
            iter->first.writeasbinary(&f);
            f.write( (char*) &iter->second, sizeof(int) ); //occurrence count                                     
            if (DOINDEX) {
                const int indexcount = ngram_index[iter->first].size();
                f.write( (char*) &indexcount, sizeof(int));
                for (set<int>::iterator iter2 = ngram_index[iter->first].begin(); iter2 != ngram_index[iter->first].end(); iter2++) {                    
                    const int c = *iter2;
                    f.write( (char*) &c, sizeof(int) );
                }                
            } else {
                f.write( (char*) &zero,sizeof(int)); //#indices
            }
        }    
    }    
    for (int n = 1; n <= MAXLENGTH; n++) {                
        for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {                            
            //f.write((char*)  &check, sizeof(char)); 
            iter->first.writeasbinary(&f);
            f.write( (char*) &(iter->second.count), sizeof(int) ); //occurrence count                         
            if (DOSKIPCONTENT) {
                const int nrofskips = iter->second.skips.size();                
                f.write( (char*) &nrofskips, sizeof(int) );
                for(skipgram_freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {                    
                    iter2->first.writeasbinary(&f);                    
                    f.write( (char*) &(iter2->second), sizeof(int) ); //occurrence count
                }
            } else {
                f.write( (char*) &zero,sizeof(int) ); //#skips     
            }
            if (DOINDEX) {
                const int indexcount = skipgram_index[iter->first].size();
                f.write( (char*) &indexcount, sizeof(int) ); //#skips
                for (set<int>::iterator iter2 = skipgram_index[iter->first].begin(); iter2 != skipgram_index[iter->first].end(); iter2++) {
                    const int c = *iter2;
                    f.write( (char*) &c, sizeof(int) );
                }                
            } else {
                f.write( (char*) &zero,sizeof(int)); //#indices
            }
        }            
    }
    f.close();    
}


bool EncGramModel::exists(const EncAnyGram* key) const {    
    if (key->gapcount() == 0) {
        return (ngrams[key->n()].count(*( (EncNGram*) key) ) > 0);
    } else {
        return (skipgrams[key->n()].count(*( (EncSkipGram*) key) ) > 0);
    }
    return false;
}



int EncGramModel::count(const EncAnyGram* key) {    
    const int n = key->n();
    if (key->gapcount() == 0) {        
        if (ngrams[n].count(*( (EncNGram*) key) ) > 0) return ngrams[n][*( (EncNGram*) key) ];
    } else {
        const int n = key->n();
        if (skipgrams[n].count(*( (EncSkipGram*) key) ) > 0) return skipgrams[n][*( (EncSkipGram*) key) ].count;   
    }
    return 0;
}


double EncGramModel::freq(const EncAnyGram* key) {    
    const int n = key->n();
    if (key->gapcount() == 0) {        
        if (ngrams[n].count(*( (EncNGram*) key) ) > 0) return ngrams[n][*( (EncNGram*) key) ] / tokens();
    } else {
        if (skipgrams[n].count( *( (EncSkipGram*) key)) > 0) return skipgrams[n][ *( (EncSkipGram*) key)].count / tokens();
    }
    return 0;
}


double EncGramModel::relfreq(const EncAnyGram* key) {    
    const int n = key->n();
    if (key->gapcount() == 0) {        
        if (ngrams[n].count(*( (EncNGram*) key) ) > 0) return ngrams[n][*( (EncNGram*) key) ] / tokencount[n];
    } else {
        if (skipgrams[n].count( *( (EncSkipGram*) key)) > 0) return skipgrams[n][ *( (EncSkipGram*) key)].count / skiptokencount[n];
    }
    return 0;
}


set<int> * EncGramModel::index(const EncAnyGram* key) {
    if (key->gapcount() == 0) {        
        if (ngram_index.count(*( (EncNGram*) key) ) > 0) return &ngram_index[*( (EncNGram*) key) ];
    } else {
        if (skipgram_index.count( *( (EncSkipGram*) key)) > 0) return &skipgram_index[ *( (EncSkipGram*) key)];
    }
}

int EncGramModel::index_size() const {
    return ngram_index.size() + skipgram_index.size();
}




std::set<int> EncGramModel::reverse_index_keys() {
    set<int> keys;
    for (unordered_map<int,vector<EncNGram> >::iterator iter = ngram_reverse_index.begin(); iter != ngram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }   
    for (unordered_map<int,vector<EncSkipGram> >::iterator iter = skipgram_reverse_index.begin(); iter != skipgram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }    
    return keys;
}


int EncGramModel::reverse_index_size(const int i) {
    int s = 0;
    if (ngram_reverse_index.count(i)) s += ngram_reverse_index[i].size();
    //if (skipgram_reverse_index.count(i)) s += skipgram_reverse_index[i].size();
    return s;
    
}

int EncGramModel::reverse_index_size() {
    return ngram_reverse_index.size();//+ skipgram_reverse_index.size();
}

bool EncGramModel::reverse_index_haskey(const int i) const {
    return ((ngram_reverse_index.count(i) > 0)); //|| (skipgram_reverse_index.count(i) > 0));
}


vector<EncAnyGram*> EncGramModel::reverse_index(const int i) {
    vector<EncAnyGram*> revindex;
    if (ngram_reverse_index.count(i) > 0)
        for (vector<EncNGram>::iterator iter = ngram_reverse_index[i].begin(); iter != ngram_reverse_index[i].end(); iter++) {
            revindex.push_back(&(*iter));
        }   
    if (skipgram_reverse_index.count(i) > 0)        
        for (vector<EncSkipGram>::iterator iter = skipgram_reverse_index[i].begin(); iter != skipgram_reverse_index[i].begin(); iter++) {
            revindex.push_back(&(*iter));
        }   
    return revindex;
}


EncAnyGram* EncGramModel::get_reverse_index_item(const int key, const int i) {
    //const int s = ngram_reverse_index[key].size();
    //if (i < s) {                
        vector<EncNGram>::iterator iter = ngram_reverse_index[key].begin() + i;
        return &(*iter);
    //} else {
        //vector<EncSkipGram>::iterator iter = skipgram_reverse_index[key].begin() + (i - s);
        //return &(*iter);
    //}    
}

void EncGramModel::decode(ClassDecoder & classdecoder, ostream *NGRAMSOUT, ostream *SKIPGRAMSOUT) {
    const int grandtotal = ngramtokencount + skipgramtokencount;   
    for (int n = 1; n <= MAXLENGTH; n++) {   
        for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
           const double freq1 = (double) iter->second / tokencount[n];
           const double freq2 = (double) iter->second / ngramtokencount;
           const double freq3 = (double) iter->second / grandtotal;
           const EncNGram ngram = iter->first;
           *NGRAMSOUT << (int) ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq1 << '\t' << freq2 << '\t' << freq3;
            if (DOINDEX) {
                *NGRAMSOUT << '\t';
                for (set<int>::iterator iter2 = ngram_index[iter->first].begin(); iter2 != ngram_index[iter->first].end(); iter2++) {
                    *NGRAMSOUT << *iter2 << ' ';
                }                
            }
            *NGRAMSOUT << endl;
        }
       

       if (SKIPGRAMSOUT != NULL) {
           for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {
               const double freq1 = (double) iter->second.count / skiptokencount[n]; 
               const double freq2 = (double) iter->second.count / skipgramtokencount;           
               const double freq3 = (double) iter->second.count / grandtotal;                          
               const EncSkipGram skipgram = iter->first;                              
               *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';
               const int skiptypes = iter->second.skips.size();               
               const double entropy = compute_entropy(iter->second.skips, iter->second.count);
               *SKIPGRAMSOUT << skiptypes << '\t' << iter->second.count << '\t' << entropy << '\t';
                for(skipgram_freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                    const EncSkipGram skipcontent = iter2->first;
                    *SKIPGRAMSOUT << skipcontent.decode(classdecoder) << '|' << iter2->second << '|';
                }
                if (DOINDEX) {
                    *SKIPGRAMSOUT << '\t';
                    for (set<int>::iterator iter2 = skipgram_index[iter->first].begin(); iter2 != skipgram_index[iter->first].end(); iter2++) {
                        *SKIPGRAMSOUT << *iter2 << ' ';
                    }                
                }
               *SKIPGRAMSOUT << endl;
           }
        }
    }
}


double compute_entropy(freqlist & data, const int total) {
    double entropy = 0;
    for(freqlist::iterator iter = data.begin(); iter != data.end(); iter++ ) {
      double p = iter->second / (double) total;
      //cout << setprecision(numeric_limits<double>::digits10 + 1) << iter->second << " / " << total << " = " << p << endl;
      entropy += p * log2(p);
    }    
    return -1 * entropy;
}

double compute_entropy(skipgram_freqlist & data, const int total) {
    double entropy = 0;
    for(skipgram_freqlist::iterator iter = data.begin(); iter != data.end(); iter++ ) {
      double p = iter->second / (double) total;
      //cout << setprecision(numeric_limits<double>::digits10 + 1) << iter->second << " / " << total << " = " << p << endl;
      entropy += p * log2(p);
    }    
    return -1 * entropy;
}





EncGramGraphModel::EncGramGraphModel(EncGramModel& model) {
    
    for (int n = 2; n <= model.maxlength(); n++) {
        cerr << "Computing subsumption relations on " << n << "-grams" << endl;
        for(freqlist::iterator iter = model.ngrams[n].begin(); iter != model.ngrams[n].end(); iter++ ) {
            const EncNGram * ngram = &(iter->first);
            vector<EncNGram*> subngrams;
            ngram->subngrams(subngrams);
            for (vector<EncNGram*>::iterator iter2 = subngrams.begin(); iter2 != subngrams.end(); iter2++) {
                const EncAnyGram * subngram = *iter2;
                if (model.exists(subngram)) {
                    //subgram exists, add relation:
                    rel_subsumption_children[*ngram].insert(*subngram);
                                        
                    //reverse:
                    //rel_subsumption_parents[*subgram].insert(*ngram);
                }
                //free memory:
                delete subngram;
            }        
        }
    }   
    
    for (int n = 2; n <= model.maxlength(); n++) {
        cerr << "Computing subsumption relations on skip-" << n << "-grams" << endl;
        for(skipgrammap::iterator iter = model.skipgrams[n].begin(); iter != model.skipgrams[n].end(); iter++ ) {        
            const EncSkipGram * skipgram = &(iter->first);
            
            vector<EncNGram*> parts;
            skipgram->parts(parts);
            for (vector<EncNGram*>::iterator iter2 = parts.begin(); iter2 != parts.end(); iter2++) {
                const EncNGram * ngram = *iter2;
                vector<EncNGram*> subngrams;
                ngram->subngrams(subngrams);
                for (vector<EncNGram*>::iterator iter3 = subngrams.begin(); iter3 != subngrams.end(); iter3++) {
                    //subgram exists, add relation:
                    const EncAnyGram * subngram = *iter3;
                    rel_subsumption_children[*skipgram].insert(*subngram);
                    delete subngram;
                }
                delete ngram;
            }    
        }
    }
    
    
     
}

EncGramGraphModel::EncGramGraphModel(const string & filename) {
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if (!f) {
           cerr << "File does not exist: " << filename << endl;
           exit(3);
    }
    
    unsigned long supernodes;
    f.read( (char*) &supernodes, sizeof(unsigned long));        
    
    for (int i = 0; i < supernodes; i++) {
        char gapcount;
        f.read(&gapcount, sizeof(char));
        EncAnyGram * anygram;
        if (gapcount == 0) {
            anygram = new EncNGram(&f);
        } else {
            anygram = new EncSkipGram(&f, gapcount);
        }
        unsigned int relations;
        f.read( (char*) &relations, sizeof(int));        
        for (int j = 0; j < relations; j++) {
            char gapcount2;
            f.read(&gapcount2, sizeof(char));
            EncAnyGram * anygram2;
            if (gapcount2 == 0) {
                anygram2 = new EncNGram(&f);
            } else {
                anygram2 = new EncSkipGram(&f, gapcount);
            }
            //rel_subsumption_children[anygram].insert(anygram2); //TODO: FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            delete anygram2;
        }
        delete anygram;
    }    
    f.close();
}

void EncGramGraphModel::save(const string & filename) {
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    
    const char czero = 0;
    const int zero = 0;
    //const unsigned char check = 255;
    const unsigned long supernodes = rel_subsumption_children.size();
    
    f.write( (char*) &supernodes, sizeof(unsigned long) );    
    for(std::unordered_map<EncAnyGram,std::unordered_set<EncAnyGram> >::iterator iter = rel_subsumption_children.begin(); iter != rel_subsumption_children.end(); iter++ ) {        
        const EncAnyGram * supergram = &(iter->first);
        if (supergram->gapcount() == 0) {
            f.write( &czero, sizeof(char));
        }
        supergram->writeasbinary(&f);
        const int relations = iter->second.size();
        f.write( (char*) &relations, sizeof(int) );        
        for(std::unordered_set<EncAnyGram>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++ ) {            
            const EncAnyGram * subgram = &(*iter2); //ugly, I know
            if (subgram->gapcount() == 0) {
                f.write( &czero, sizeof(char));
            }
            subgram->writeasbinary(&f);
        }
    }
    f.close();
    
}


EMAlignmentModel::EMAlignmentModel(EncGramModel & sourcemodel, EncGramModel & targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD) {
    int round = 0;    
    unsigned long c = 0;
    double totaldivergence = 0;
    double prevavdivergence = 0;
    bool converged = false;
    
    set<int> reverseindexkeys = sourcemodel.reverse_index_keys();
    pair<EncAnyGram*,EncAnyGram*> targetgram_given_sourcegram;
            
    do {       
        round++; 
        cerr << "  EM Round " << round << "... ";
        //use reverse index to iterate over all sentences
        for (set<int>::iterator iter = reverseindexkeys.begin(); iter != reverseindexkeys.end(); iter++) {
            const int key = *iter;        
            const int sourcegrams_size =  sourcemodel.reverse_index_size(key);
            //vector<EncAnyGram*> sourcegrams = sourcemodel.reverse_index(key);
            const int targetgrams_size =  targetmodel.reverse_index_size(key);
            
            if (targetgrams_size > 0 ) { //target model contains sentence?               
                //vector<EncAnyGram*> targetgrams = targetmodel.reverse_index(key);  
                cerr << key << "=" << sourcegrams_size << "x" << targetgrams_size << " ";

                //compute sentencetotal for normalisation later in count step, sum_s(p(t|s))                
                cerr << "A";                            
                unordered_map<EncAnyGram*, double> sentencetotal;                
                for (int i = 0; i < targetgrams_size; i++) {    
                    EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        targetgram_given_sourcegram.first = targetgram;
                        targetgram_given_sourcegram.second = sourcegram;                        
                        sentencetotal[targetgram] += transprob[targetgram_given_sourcegram]; //compute sum over all source conditions for a targetgram under consideration
                    }
                }               
                
                cerr << "B";                            
                //collect counts (for evidence that a targetgram is aligned to a sourcegram)
                alignmentprobmap count;                
                unordered_map<EncAnyGram*, double> total;
                for (int i = 0; i < targetgrams_size; i++) {    
                    EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        targetgram_given_sourcegram.first = targetgram;
                        targetgram_given_sourcegram.second = sourcegram;
                        
                        const double countvalue = transprob[targetgram_given_sourcegram] / sentencetotal[targetgram];
                        count[targetgram_given_sourcegram] += countvalue;
                        total[targetgram] += countvalue;
                    }                    
                }
                
                 
                cerr << "C"; 
                double prevtransprob;
                //estimate new probabilities (normalised count is the new estimated probability)
                for (int i = 0; i < targetgrams_size; i++) {    
                    EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);
                        targetgram_given_sourcegram.first = targetgram;
                        targetgram_given_sourcegram.second = sourcegram;
                        
                        prevtransprob = transprob[targetgram_given_sourcegram];
                        const double newtransprob = count[targetgram_given_sourcegram] / total[targetgram];
                        transprob[targetgram_given_sourcegram] = newtransprob;                        
                        
                        //for computation of convergence
                        const double divergence = abs(newtransprob - prevtransprob);
                        totaldivergence += divergence;
                        c++;
                    }
                }
            }
        }
        const double avdivergence = totaldivergence / c;
        converged = (((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence)) <= CONVERGEDTHRESHOLD);       
        prevavdivergence = avdivergence;
        cerr << " average divergence = " << avdivergence << ", transprob size = " << transprob.size() << endl;
    } while (!converged);    
}


double CoocAlignmentModel::cooc( set<int> & sourceindex, set<int> & targetindex ) {    
    //Jaccard co-occurrence
    
    int intersectioncount = 0;    

    //count union and intersections
    for (set<int>::iterator iter = sourceindex.begin(); iter != sourceindex.end(); iter++) {
        for (set<int>::iterator iter2 = sourceindex.begin(); iter2 != sourceindex.end(); iter2++) {    
            if (*iter == *iter2) intersectioncount++;                        
        }        
    }    
    const int unioncount = (sourceindex.size() + targetindex.size()) - intersectioncount; 
    return intersectioncount / unioncount;
}

CoocAlignmentModel::CoocAlignmentModel(EncGramModel & sourcemodel, EncGramModel & targetmodel) {
    int c = 0;
    //Count co-occurence
    for (unordered_map<EncNGram,set<int> >::iterator iter = sourcemodel.ngram_index.begin();  iter != sourcemodel.ngram_index.end(); iter++) {        
        cerr << "N" << ++c << " ";
        const EncAnyGram * sourcegram = &iter->first;
        for (unordered_map<EncNGram,set<int> >::iterator iter2 = targetmodel.ngram_index.begin();  iter2 != targetmodel.ngram_index.end(); iter2++) {
            const EncAnyGram * targetgram = &iter2->first;
            cerr << iter->second.size() << "n" << iter2->second.size() << " ";
            alignprob[sourcegram][targetgram] = cooc(iter->second, iter2->second);
        }
        for (unordered_map<EncSkipGram,set<int> >::iterator iter2 = targetmodel.skipgram_index.begin();  iter2 != targetmodel.skipgram_index.end(); iter2++) {
            const EncAnyGram * targetgram = &iter2->first;
            cerr << iter->second.size() << "s" << iter2->second.size() << " ";
            alignprob[sourcegram][targetgram] = cooc(iter->second, iter2->second);
        }        
    }
    for (unordered_map<EncSkipGram,set<int> >::iterator iter = sourcemodel.skipgram_index.begin();  iter != sourcemodel.skipgram_index.end(); iter++) {
        cerr << ++c << "S ";
        const EncAnyGram * sourcegram = &iter->first;
        for (unordered_map<EncNGram,set<int> >::iterator iter2 = targetmodel.ngram_index.begin();  iter2 != targetmodel.ngram_index.end(); iter2++) {
            const EncAnyGram * targetgram = &iter2->first;
            cerr << iter->second.size() << "n" << iter2->second.size() << " ";
            alignprob[sourcegram][targetgram] = cooc(iter->second, iter2->second);
        }
        for (unordered_map<EncSkipGram,set<int> >::iterator iter2 = targetmodel.skipgram_index.begin();  iter2 != targetmodel.skipgram_index.end(); iter2++) {
            const EncAnyGram * targetgram = &iter2->first;
            cerr << iter->second.size() << "s" << iter2->second.size() << " ";
            alignprob[sourcegram][targetgram] = cooc(iter->second, iter2->second);
        }        
    }
        
}
    
