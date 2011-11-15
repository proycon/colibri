#include "ngram.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithms.h>
#include <limits>

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

bool EncNGram::out() const {
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

/*
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
}*/



EncSkipGram::EncSkipGram(const vector<EncNGram*> & dataref, const vector<int> & skipref, bool initialskip, bool finalskip): EncNGram() {
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

EncSkipGram::EncSkipGram(const unsigned char *dataref, const char size, const unsigned char* skipref, const char skipcount) {
    _size = size;
    data = new unsigned char[size];
    for (int i = 0; i < size; i++) {
        data[i] = dataref[i];
    }   
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



EncGramModel::EncGramModel(const string corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  int MINSKIPTYPES, bool DOINDEX, bool DOSKIPCONTENT, bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
    this->DOINDEX = DOINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;
    this->DOINITIALONLYSKIP = DOINITIALONLYSKIP;
    this->DOFINALONLYSKIP = DOFINALONLYSKIP;
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;

    unsigned char line[65536];

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
            const int linesize = readline(IN, line);            
                    
            linenum++;

            if (linenum % 10000 == 0) {
                cerr << "\t@" << linenum << endl;
            }
                            
            
            const int l = countwords(line, linesize);            
            
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
                    bool prunedskip = false;
                    for(skipgram_freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                        if (iter2->second < MINSKIPTOKENS) {
                            //prune skip
                            prunedskip = true;
                            iter->second.count -= iter2->second;
                            iter->second.skips.erase(iter2->first); 
                        }
                    }
                    if ( (prunedskip) && ( (iter->second.skips.size() < MINSKIPTYPES) || (iter->second.count < MINTOKENS) ) ) { //reevaluate
                        pruneskipgram = true;
                    } else {
                        skipgramtokencount += iter->second.count;
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
        
        if (DOINDEX) { //TODO: REFACTOR
            /*cerr << "Writing ngram index" << endl;

            for(unordered_map<EncNGram,set<int>>::iterator iter = ngram_index.begin(); iter != ngram_index.end(); iter++ ) {
                const EncNGram ngram = iter->first;
                *NGRAMINDEX << ngram.decode(classdecoder) << '\t';
                for (set<int>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                    *NGRAMINDEX << *iter2 << ' ';
                }
                *NGRAMINDEX << endl;
            }
                        
            if (DOSKIPGRAMS) {
                cerr << "Writing skipgram index" << endl;;

                for(unordered_map<EncSkipGram,set<int>>::iterator iter = skipgram_index.begin(); iter != skipgram_index.end(); iter++ ) {
                    const EncSkipGram skipgram = iter->first;
                    *SKIPGRAMINDEX << (int) skipgram.n() << '\t' << skipgram.decode(classdecoder) << '\t';
                    for (set<int>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                        *SKIPGRAMINDEX << *iter2 << ' ';
                    }
                    *SKIPGRAMINDEX << endl;
                }
                
            }*/
        }
    
        ngramtypecount += ngrams[n].size();
        skipgramtypecount += skipgrams[n].size();     
    }

        
}

EncGramModel::EncGramModel(string filename) {
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    char buffer[1024];
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    
    f.read(buffer, sizeof(unsigned long));
    const unsigned long totaltokens = atol(buffer);
    f.read(buffer, sizeof(unsigned long));
    const unsigned long totaltypes = atol(buffer);
    for (int i = 0; i < totaltypes; i++) {
        f.read(buffer, sizeof(char));
        const char gapcount = buffer[0];
        if (gapcount == 0) {
            //NGRAM
            f.read(buffer, sizeof(char));
            const char size = buffer[0];
            f.read(buffer, (int) size); //read data
            EncNGram ngram = EncNGram( (unsigned char*) buffer, size);
            if (ngram.n() > MAXLENGTH) {
                for (int j = 0; j < ngram.n() - MAXLENGTH; j++) {
                    ngrams.push_back(freqlist());
                    skipgrams.push_back(skipgrammap());
                }
            }
            f.read(buffer, sizeof(int)); //read occurrence count
            ngrams[ngram.n()][ngram] = atoi(buffer); //assign count 
            f.read(buffer, sizeof(int));
            const int indexcount = atoi(buffer);
            for (int j = 0; j < indexcount; j++) {
                f.read(buffer, sizeof(int));
                const int index = atoi(buffer);
                ngram_index[ngram].insert(index);
            }
        } else {
            //SKIPGRAM
            char skipref[4];
            for (int j = 0; j < gapcount; j++) {
                f.read(buffer, sizeof(char));
                skipref[j] = buffer[0];
            }
            f.read(buffer, sizeof(char));
            const char size = buffer[0];
            f.read(buffer, (int) size); //read data
            EncSkipGram skipgram = EncSkipGram( (unsigned char*) buffer, size, (unsigned char*) skipref, gapcount);
            if (skipgram.n() > MAXLENGTH) {
                for (int j = 0; j < skipgram.n() - MAXLENGTH; j++) {
                    ngrams.push_back(freqlist());
                    skipgrams.push_back(skipgrammap());
                }
            }            
            f.read(buffer, sizeof(int)); //read occurrence count
            skipgrams[skipgram.n()][skipgram].count = atoi(buffer); //assign
            f.read(buffer, sizeof(int));   
            const int skipcontentcount = atoi(buffer);
            for (int j = 0; j < skipcontentcount; j++) {
                f.read(buffer, sizeof(char));   
                const char skip_gapcount = buffer[0];
                char skipref2[4];
                for (int j = 0; j < skip_gapcount; j++) {
                    f.read(buffer, sizeof(char));
                    skipref2[j] = buffer[0];
                }
                f.read(buffer, sizeof(char));   
                const char skip_size = buffer[0];
                f.read(buffer, (int) skip_size); //read data
                EncSkipGram skipcontent = EncSkipGram((unsigned char*) buffer, size, (unsigned char*) skipref2, skip_gapcount);
                f.read(buffer, sizeof(int)); //read occurrence count
                skipgrams[skipgram.n()][skipgram].skips[skipcontent] = atoi(buffer); //skipcontent occurrence 
            }
            f.read(buffer, sizeof(int));
            const int indexcount = atoi(buffer);
            for (int j = 0; j < indexcount; j++) {
                f.read(buffer, sizeof(int));
                const int index = atoi(buffer);
                skipgram_index[skipgram].insert(index);
            }
        }        
    }
    f.close();
}



void EncGramModel::save(std::string filename) {
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    
    const char czero = 0;
    const int zero = 0;
    const int totaltokens = tokens();
    const int totaltypes = types();

    
    f.write( (char*) &totaltokens, sizeof(unsigned long) );
    f.write( (char*) &totaltypes, sizeof(unsigned long) );
    for (int n = 1; n <= MAXLENGTH; n++) {
        for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {        
            const int size = iter->first.size();
            f.write( &czero,sizeof(char) ); //#number of gaps (always zero for ngrams)
            f.write( (char*) &size, sizeof(char) ); //data length
            f.write( (char*) iter->first.data , iter->first.size() ); //data
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
            int size = iter->first.size();
            f.write( &iter->first.skipcount, sizeof(char) ); //nr of gaps
            for (int j = 0; j < iter->first.skipcount; j++) { //skip configuration
                    f.write( &iter->first.skipsize[j], sizeof(char) );
            }
            f.write( (char*) &size, sizeof(char) ); //size
            f.write( (char*) iter->first.data , iter->first.size() ); //data
            f.write( (char*) &iter->second.count, sizeof(int) ); //occurrence count                         
            if (DOSKIPCONTENT) {
                const int skipcount = skipgram_index[iter->first].size();                
                f.write( (char*) &skipcount, sizeof(int) );
                for(skipgram_freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                    size = iter2->first.size();
                    f.write( &iter2->first.skipcount, sizeof(char) ); //nr of gaps
                    for (int j = 0; j < iter2->first.skipcount; j++) { //skip configuration
                            f.write( &iter2->first.skipsize[j], sizeof(char) );
                    }
                    f.write( (char*) &size, sizeof(char) ); //size
                    f.write( (char*) iter2->first.data , iter2->first.size() ); //data
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

bool EncGramModel::exists(EncNGram* key) const {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (ngrams[n].count(*key) > 0) return true;
    }
    return false;
}

bool EncGramModel::exists(EncSkipGram* key) const {    
    for (int n = 1; n <= MAXLENGTH; n++) {
            if (skipgrams[n].count(*key) > 0) return true;
    }
    return false;
}

int EncGramModel::count(EncNGram* key) {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (ngrams[n].count(*key) > 0) 
         return ngrams[n][*key];
    }
    return 0;
}

int EncGramModel::count(EncSkipGram* key) {
    for (int n = 1; n <= MAXLENGTH; n++) {
        const EncSkipGram* key2 = (EncSkipGram*) key;
        if (skipgrams[n].count(*key2) > 0) 
         return skipgrams[n][*key2].count;
    }
    return 0;
}

double EncGramModel::freq(EncNGram* key) {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (ngrams[n].count(*key) > 0) return ngrams[n][*key] / tokens();
    }
    return 0;
}

double EncGramModel::freq(EncSkipGram* key) {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (skipgrams[n].count( *( (EncSkipGram*) key)) > 0) return skipgrams[n][ *( (EncSkipGram*) key)].count / tokens();
    }
    return 0;    
}


double EncGramModel::relfreq(EncNGram* key) {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (ngrams[n].count(*key) > 0) return ngrams[n][*key] / skiptokencount[n];
    }
    return 0;
}

double EncGramModel::relfreq(EncSkipGram* key) {    
    for (int n = 1; n <= MAXLENGTH; n++) {
        if (skipgrams[n].count( *( (EncSkipGram*) key)) > 0) return skipgrams[n][ *( (EncSkipGram*) key)].count / tokencount[n];
    }
    return 0;
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

