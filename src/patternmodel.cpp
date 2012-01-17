#include "patternmodel.h"
#include "algorithms.h"
#include <limits>

using namespace std;


/*enum ModelType { 
    INDEXEDPATTERNMODEL = 10,
    GRAPHPATTERNMODEL = 20,
}*/


CorpusReference::CorpusReference(uint32_t sentence, unsigned char token) {
    this->sentence = sentence;
    this->token = token;
}

CorpusReference::CorpusReference(istream * in) {
    in->read( (char*) &sentence, sizeof(uint32_t)); 
    in->read( (char*) &token, sizeof(unsigned char)); 
}

void CorpusReference::writeasbinary(ostream * out) const {  
    out->write( (char*) &sentence, sizeof(uint32_t) );
    out->write( (char*) &token, sizeof(unsigned char) );
}


void NGramData::writeasbinary(ostream * out) const {
    const uint32_t s = refs.size();
    out->write( (char*) &s, sizeof(uint32_t) );
    for (set<CorpusReference>::iterator iter = refs.begin(); iter != refs.end(); iter++) {
        iter->writeasbinary(out);
    }
}



void ModelReader::readfile(const string & filename) {
    const bool DEBUG = false;
    

    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    

    uint64_t id = this->id();
    f.read( (char*) &id, sizeof(uint64_t));        
    uint64_t totaltokens;
    f.read( (char*) &totaltokens, sizeof(uint64_t));        
    uint64_t totaltypes;
    f.read( (char*) &totaltypes, sizeof(uint64_t)); 
    
    readheader(&f, totaltokens, totaltypes);
    
    for (int i = 0; i < totaltypes; i++) {           
        char gapcount;
        if (DEBUG) cerr << "\t@" << i;
        f.read(&gapcount, sizeof(char));
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            const EncNGram ngram = EncNGram(&f); //read from file            
            readngram(&f, ngram);          
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            const EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            readskipgram(&f, skipgram);
        }
        if (DEBUG)  cerr << endl;      //DEBUG  
    }
    readfooter(&f);    
    f.close();
}

void ModelWriter::writefile(const string & filename) {
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    

    uint64_t _id = this->id();
    f.write( (char*) &_id, sizeof(uint64_t));        
    uint64_t _totaltokens = tokens();
    f.write( (char*) &_totaltokens, sizeof(uint64_t));        
    uint64_t _totaltypes = types();
    f.write( (char*) &_totaltypes, sizeof(uint64_t)); 
        
    writeheader(&f);
    writengrams(&f);
    writeskipgrams(&f);
    writefooter(&f);
    f.close();
}


IndexedPatternModel::IndexedPatternModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  int MINSKIPTYPES, bool DOREVERSEINDEX, bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = true;
    this->DOINITIALONLYSKIP = DOINITIALONLYSKIP;
    this->DOFINALONLYSKIP = DOFINALONLYSKIP;
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;

    const int BUFFERSIZE = 65536;
    unsigned char line[BUFFERSIZE];


    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;            
        
        uint32_t sentence = 0;
    
        tokencount[n] = 0;
        skiptokencount[n] = 0;
        typecount[n] = 0;
        skiptypecount[n] = 0;
                
        vector< vector< pair<int,int> > > gaps;
        compute_multi_skips(gaps, vector<pair<int,int> >(), n);    
        
        
        ifstream *IN =  new ifstream( corpusfile.c_str() );    
        vector<unsigned int> words;
        while (IN->good()) {
            const int linesize = readline(IN, line, BUFFERSIZE );            
                    
            sentence++;

            if (sentence % 10000 == 0) {
                cerr << "\t@" << sentence << endl;
            }
                            
            
            const int l = countwords(line, linesize);            
            if (l >= 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping!" << endl;
                continue;
            }
            
            if (linesize > 0) //no { on purpose! applies to next for loop
            for (unsigned char i = 0; ((i < l - n + 1) && (i < 256)); i++) {                
                EncNGram * ngram = getencngram(i,n, line, linesize, sentence);  
                                                

                //cout << "NGRAM("<<ngram->n()<<","<<(int)ngram->size() << ")" << endl;
                              
                if (n > 2) {                    
                    EncNGram * subngram1 = ngram->slice(0, n - 1);
                    if (!(ngrams.count(*subngram1))) {
                        delete subngram1;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram1;

                                                             
                    EncNGram * subngram2 = ngram->slice(1, n - 1);
                    if (!(ngrams.count(*subngram2))) {
                        delete subngram2;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram2;                    
                }
                
                
                
                CorpusReference ref = CorpusReference(sentence,i);
                if (ngrams[*ngram].refs.empty()) typecount[n]++;
                ngrams[*ngram].refs.insert(ref);                        
                tokencount[n]++;            

                if (DOSKIPGRAMS) {
                	//Iterate over all possible gap configurations
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
                        
                    	//Iterate over all gaps in this configuration    
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
                                oc = ngrams.count(*subngram);
                                if (oc) oc = ngrams[*subngram].count();
                                if ((((oc == 0)) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) ))    {
                                    docount = false;
                                    break;
                                }
                            }
                            cursor = begin + length;
                        }   
                        if (cursor < n) {
                            EncNGram * subngram = ngram->slice(cursor,n-cursor);
                            subngrams.push_back(subngram);
                            oc = ngrams.count(*subngram);
                            if (oc) oc = ngrams[*subngram].count();
                            if (((oc == 0)) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) )  {
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
                                oc = ngrams.count(*subskip);
                                if (oc) oc = ngrams[*subskip].count();
                                if ((oc == 0) || (oc < MINTOKENS) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS) ) )  {
                                	//skip content doesn't occur often enough, abort count
                                	docount = false;
                                	break;
                                }
                                if (cursor > 0) skipcontent_skipref.push_back(begin - cursor);
                                cursor = begin+length;
                            }   
                            if (docount) { //second check, we might have aborted
		                        EncSkipGram skipcontent = EncSkipGram(skipcontent_subngrams, skipcontent_skipref, false, false);                                                        
		                        if (skipgrams[skipgram].count() == 0) skiptypecount[n]++;                            
		                        skipgrams[skipgram]._count++;
		                        skipgrams[skipgram].skipcontent[skipcontent].refs.insert(ref);
		                        skiptokencount[n]++;                            
                            }
                            for (size_t k = 0; k < skipcontent_subngrams.size(); k++) {       
		                        delete skipcontent_subngrams[k];
		                    }             
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

       cerr << "Found " << typecount[n] << " " << n << "-grams (" << tokencount[n] << " tokens)";
       if (DOSKIPGRAMS) {
        cerr << " and " << skiptypecount[n] << " skipgrams (" << skiptokencount[n] << " tokens)" << endl;
       } else {
        cerr << endl;
       }
    

       //prune n-grams
       int pruned = 0;
       for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
            if (iter->first.n() == n) {
                if (iter->second.count() < MINTOKENS) {
                    tokencount[n] -= iter->second.count();
                    typecount[n]--;
                    pruned++;
                    ngrams.erase(iter->first);                        
                } else {
                    ngramtokencount += iter->second.count();
                }
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << typecount[n] <<  " left (" << tokencount[n] << " tokens)" << endl;
    
       
       if (DOSKIPGRAMS) {       
           //prune skipgrams
           pruned = 0;
           for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {     
               if (iter->first.n() == n) {                          
                    bool pruneskipgram = false;
                    if ((iter->second.count() < MINTOKENS) || ((DOSKIPCONTENT && iter->second.skipcontent.size() < MINSKIPTYPES)))  {
                        pruneskipgram = true;
                    } else {             // if (DOSKIPCONTENT)
                        int prunedskiptokens = 0;
                        for(std::unordered_map<EncSkipGram,NGramData>::iterator iter2 = iter->second.skipcontent.begin(); iter2 != iter->second.skipcontent.end(); iter2++ ) {
                            if (iter2->second.count() < MINSKIPTOKENS) {
                                //prune skip
                                iter->second._count -= iter2->second.count();
                                prunedskiptokens += iter2->second.count();
                                iter->second.skipcontent.erase(iter2->first); 
                            }
                        }
                        if ( (prunedskiptokens > 0) && ( (iter->second.skipcontent.size() < MINSKIPTYPES) || (iter->second.count() - prunedskiptokens < MINTOKENS) ) ) { //reevaluate
                            pruneskipgram = true;
                            skiptokencount[n] -= prunedskiptokens;
                        } else {
                            //skipgramtokencount += iter->second.count;                            
                        }
                    }
                    if (pruneskipgram) {
                        skiptokencount[n] -= iter->second.count();
                        skiptypecount[n]--;
                        pruned++;
                        skipgrams.erase(iter->first);
                    }
               }
           }
           cerr << "Pruned " << pruned << " skipgrams, " << skiptypecount[n] <<  " left (" << skiptokencount[n] << " tokens)" << endl;
           
        }
        ngramtokencount += tokencount[n];
        ngramtypecount += typecount[n];
        skipgramtokencount += skiptokencount[n];
        skipgramtypecount += skiptypecount[n];
    }

    if (DOREVERSEINDEX) {
        //TODO: Compute reverse index        
    }
        
}


IndexedPatternModel::IndexedPatternModel(const string & filename, bool DOREVERSEINDEX) {    
    const bool DEBUG = true;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;    
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    if (!filename.empty()) readfile(filename);
}

void IndexedPatternModel::readngram(std::istream * f, const EncNGram & ngram) {
    ngramtypecount++;
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    for (int j = 0; j < count; j++) {
        CorpusReference ref = CorpusReference(f); //read from file
        ngrams[ngram].refs.insert(ref);
        /*if (DOREVERSEINDEX) {
            bool found = false;
            for (int k = 0; k < ngram_reverse_index[index].size(); k++) if (ngram_reverse_index[index][k] == ngram) { found = true; break; };
            if (!found) ngram_reverse_index[index].push_back(ngram);
        }*/
               
    }
    ngramtokencount += ngrams[ngram].count();    
    if (ngram.n() > MAXLENGTH) MAXLENGTH = ngram.n();
    tokencount[ngram.n()] += ngrams[ngram].count();
}



void IndexedPatternModel::readskipgram(std::istream * f, const EncSkipGram & skipgram) {
    skipgramtypecount++;
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count            
    skipgrams[skipgram]._count = count; //assign
    skipgramtokencount += count;
    skiptokencount[skipgram.n()] += count;            
    uint32_t skipcontentcount;
    f->read((char*) &skipcontentcount, sizeof(uint32_t));   
    for (int j = 0; j < skipcontentcount; j++) {                                
        EncSkipGram skipcontent = EncSkipGram(f);  //also when !DOSKIPCONTENT, bytes have to be read
        f->read((char*) &count, sizeof(uint32_t)); //read occurrence count                
        for (int k = 0; k < count; k++) {
            CorpusReference ref = CorpusReference(f); //read from file
            skipgrams[skipgram].skipcontent[skipcontent].refs.insert(ref);
        }        
    }    
}

void IndexedPatternModel::writengram(std::ostream * f, const EncNGram & ngram) {
    const uint32_t c = ngrams[ngram].count();
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                                     
    for (set<CorpusReference>::iterator iter = ngrams[ngram].refs.begin(); iter != ngrams[ngram].refs.end(); iter++) {                    
        iter->writeasbinary(f);
    }                    
}


void IndexedPatternModel::writengrams(std::ostream * f) {    
    const char czero = 0;
    for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengram(f, iter->first);       
    }   
}


void IndexedPatternModel::writeskipgram(std::ostream * f, const EncSkipGram & skipgram) {
    const uint32_t c  = skipgrams[skipgram].count();
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                         
    const uint32_t skipcontentcount = (uint32_t) skipgrams[skipgram].skipcontent.size();                        
    f->write( (char*) &skipcontentcount, sizeof(uint32_t) );
    for(unordered_map<EncSkipGram,NGramData>::iterator iter = skipgrams[skipgram].skipcontent.begin(); iter != skipgrams[skipgram].skipcontent.end(); iter++ ) {                    
        iter->first.writeasbinary(f);
        iter->second.writeasbinary(f);
    }
}


void IndexedPatternModel::writeskipgrams(std::ostream * f) {                 
    for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {                                
        iter->first.writeasbinary(f);        
        writeskipgram(f, iter->first);
    }     
}


bool IndexedPatternModel::exists(const EncAnyGram* key) const {    
    if (key->gapcount() == 0) {
        return (ngrams.count(*( (EncNGram*) key) ) > 0);
    } else {
        return (skipgrams.count(*( (EncSkipGram*) key) ) > 0);
    }
    return false;
}

const EncAnyGram* IndexedPatternModel::getkey(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<EncNGram,NGramData >::iterator iter = ngrams.find(*( (EncNGram*) key) );
        if (iter != ngrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<EncSkipGram,SkipGramData >::iterator iter = skipgrams.find(*( (EncSkipGram*) key) );
        if (iter != skipgrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }        
    }
}


const AnyGramData* IndexedPatternModel::getdata(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<EncNGram,NGramData >::iterator iter = ngrams.find(*( (EncNGram*) key) );
        if (iter != ngrams.end()) {
            return &iter->second;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<EncSkipGram,SkipGramData >::iterator iter = skipgrams.find(*( (EncSkipGram*) key) );
        if (iter != skipgrams.end()) {
            return &iter->second;
        } else {
            return NULL;
        }        
    }
}



int IndexedPatternModel::count(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ].count();
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return skipgrams[*( (EncSkipGram*) key) ].count();   
    }
    return 0;
}


double IndexedPatternModel::freq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ].count() / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)].count() / tokens();
    }
    return 0;
}


double IndexedPatternModel::relfreq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ].count() / tokencount[key->n()];
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)].count() / skiptokencount[key->n()];
    }
    return 0;
}


/*set<int> * IndexedPatternModel::index(const EncAnyGram* key) {
    if (key->gapcount() == 0) {        
        if (ngram_index.count(*( (EncNGram*) key) ) > 0) return &ngram_index[*( (EncNGram*) key) ];
    } else {
        if (skipgram_index.count( *( (EncSkipGram*) key)) > 0) return &skipgram_index[ *( (EncSkipGram*) key)];
    }
}*/




std::set<int> IndexedPatternModel::reverse_index_keys() {
    set<int> keys;
    for (unordered_map<int,vector<EncNGram> >::iterator iter = ngram_reverse_index.begin(); iter != ngram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }   
    for (unordered_map<int,vector<EncSkipGram> >::iterator iter = skipgram_reverse_index.begin(); iter != skipgram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }    
    return keys;
}


int IndexedPatternModel::reverse_index_size(const int i) {
    int s = 0;
    if (ngram_reverse_index.count(i)) s += ngram_reverse_index[i].size();
    if (skipgram_reverse_index.count(i)) s += skipgram_reverse_index[i].size();
    return s;
    
}

int IndexedPatternModel::reverse_index_size() {
    return ngram_reverse_index.size() + skipgram_reverse_index.size();
}

bool IndexedPatternModel::reverse_index_haskey(const int i) const {
    return ((ngram_reverse_index.count(i) > 0) || (skipgram_reverse_index.count(i) > 0));
}


vector<EncAnyGram*> IndexedPatternModel::reverse_index(const int i) {
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


EncAnyGram* IndexedPatternModel::get_reverse_index_item(const int key, const int i) {
    const int s = ngram_reverse_index[key].size();
    if (i < s) {                
        vector<EncNGram>::iterator iter = ngram_reverse_index[key].begin() + i;
        return &(*iter);
    } else {
        vector<EncSkipGram>::iterator iter = skipgram_reverse_index[key].begin() + (i - s);
        return &(*iter);
    }    
}

void IndexedPatternModel::decode(ClassDecoder & classdecoder, ostream *NGRAMSOUT, ostream *SKIPGRAMSOUT) {
    const int grandtotal = ngramtokencount + skipgramtokencount;   

    for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
       const double freq1 = (double) iter->second.count() / tokencount[iter->first.n()];
       const double freq2 = (double) iter->second.count() / ngramtokencount;
       const double freq3 = (double) iter->second.count() / grandtotal;
       const EncNGram ngram = iter->first;
        *NGRAMSOUT << (int) ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second.count() << '\t' << freq1 << '\t' << freq2 << '\t' << freq3;
        *NGRAMSOUT << '\t';
        for (set<CorpusReference>::iterator iter2 = iter->second.refs.begin() ; iter2 != iter->second.refs.end(); iter2++) {
            *NGRAMSOUT << iter2->sentence << ':' << (int) iter2->token << ' ';
        }                
        *NGRAMSOUT << endl;
    }
   

   if (SKIPGRAMSOUT != NULL) {
       for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
           const double freq1 = (double) iter->second.count() / skiptokencount[iter->first.n()]; 
           const double freq2 = (double) iter->second.count() / skipgramtokencount;           
           const double freq3 = (double) iter->second.count() / grandtotal;                          
           const EncSkipGram skipgram = iter->first;                              
           *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count() << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';
           const int skiptypes = iter->second.skipcontent.size();               
           const double entropy = iter->second.entropy();
           *SKIPGRAMSOUT << skiptypes << '\t' << iter->second.count() << '\t' << entropy << '\t';
            for(unordered_map<EncSkipGram,NGramData>::iterator iter2 = iter->second.skipcontent.begin(); iter2 != iter->second.skipcontent.end(); iter2++ ) {
                *SKIPGRAMSOUT << iter2->first.decode(classdecoder) << '|' << iter2->second.count() << '|';
                for (set<CorpusReference>::iterator iter3 = iter2->second.refs.begin() ; iter3 != iter2->second.refs.end(); iter3++) {
                    *SKIPGRAMSOUT << iter3->sentence << ':' << (int) iter3->token;        
                    *SKIPGRAMSOUT << ',';
                    //if (iter3 != iter2->second.refs.end() - 1) 
                }
                //MAYBE TODO: output references?
            }
           *SKIPGRAMSOUT << endl;
       }
    }

}


UnindexedPatternModel::UnindexedPatternModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = true;
    this->DOINITIALONLYSKIP = DOINITIALONLYSKIP;
    this->DOFINALONLYSKIP = DOFINALONLYSKIP;
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;

    const int BUFFERSIZE = 65536;
    unsigned char line[BUFFERSIZE];


    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;            
        
        uint32_t sentence = 0;
    
        tokencount[n] = 0;
        skiptokencount[n] = 0;
        typecount[n] = 0;
        skiptypecount[n] = 0;
                
        vector< vector< pair<int,int> > > gaps;
        compute_multi_skips(gaps, vector<pair<int,int> >(), n);    
        
        unordered_map<const EncSkipGram,size_t> lastskipcontenthash; //0 if not seen yet , 1 if already verified, hash key otherwise
        
        ifstream *IN =  new ifstream( corpusfile.c_str() );    
        vector<unsigned int> words;
        while (IN->good()) {
            const int linesize = readline(IN, line, BUFFERSIZE );            
                    
            sentence++;

            if (sentence % 10000 == 0) {
                cerr << "\t@" << sentence << endl;
            }
                            
            
            const int l = countwords(line, linesize);            
            if (l >= 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping!" << endl;
                continue;
            }
            
            if (linesize > 0) //no { on purpose! applies to next for loop
            for (unsigned char i = 0; ((i < l - n + 1) && (i < 256)); i++) {                
                EncNGram * ngram = getencngram(i,n, line, linesize, sentence);  
                                                

                //cout << "NGRAM("<<ngram->n()<<","<<(int)ngram->size() << ")" << endl;
                              
                if (n > 2) {                    
                    EncNGram * subngram1 = ngram->slice(0, n - 1);
                    if (!(ngrams.count(*subngram1))) {
                        delete subngram1;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram1;

                                                             
                    EncNGram * subngram2 = ngram->slice(1, n - 1);
                    if (!(ngrams.count(*subngram2))) {
                        delete subngram2;
                        delete ngram;
                        continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                    }
                    delete subngram2;                    
                }
                
                
                if (ngrams[*ngram] == 0) typecount[n]++;
                ngrams[*ngram]++; //increase occurence count                     
                tokencount[n]++;            

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
                                oc = ngrams.count(*subngram);
                                if (oc) oc = ngrams[*subngram];
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
                            oc = ngrams.count(*subngram);
                            if (oc) oc = ngrams[*subngram];
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
                            if (skipgrams[skipgram] == 0) skiptypecount[n]++;                            
                            skipgrams[skipgram]++;
                            skiptokencount[n]++;                            
                        }
                        if (docount) {
                            EncSkipGram skipgram = EncSkipGram(subngrams, skipref, initialskip, finalskip);                                                    
                            vector<EncNGram*> skipcontent_subngrams;
                            vector<int> skipcontent_skipref;
                            cursor = 0;
                            //iterate over all gaps in this configuration (to check their occurrence count)
                            for (size_t k = 0; k < gaps[j].size(); k++) {
                                const int begin = gaps[j][k].first;  
                                const int length = gaps[j][k].second;
                                EncNGram * subskip = ngram->slice(begin,length);                                                               
                                oc = ngrams.count(*subskip);
                                if (oc) oc = ngrams[*subskip];                                
                                if (lastskipcontenthash[skipgram] != 1) {
                                	const size_t newskipcontenthash = subskip->hash(); 
									if (newskipcontenthash != lastskipcontenthash[skipgram]) {
        								lastskipcontenthash[skipgram] = 1;
        							} else {
        								lastskipcontenthash[skipgram] = newskipcontenthash;
        							}          	
                                }                             
                                delete subskip;                                     
                                if ((oc == 0) || (oc < MINTOKENS) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS) ) )  {
                                	//skip content doesn't occur often enough, abort count
                                	docount = false;
                                	break;
                                }
                                cursor = begin+length;                                
                            }   
                            if (docount) { //second check, we might have aborted                                                        
		                        if (skipgrams[skipgram] == 0) skiptypecount[n]++;                            
		                        skipgrams[skipgram]++; //increase count
		                        skiptokencount[n]++;                            		                        
                            }
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

       cerr << "Found " << typecount[n] << " " << n << "-grams (" << tokencount[n] << " tokens)";
       if (DOSKIPGRAMS) {
        cerr << " and " << skiptypecount[n] << " skipgrams (" << skiptokencount[n] << " tokens)" << endl;
       } else {
        cerr << endl;
       }
    

       //prune n-grams
       int pruned = 0;
       for(unordered_map<EncNGram,uint32_t>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
            if (iter->first.n() == n) {
                if (iter->second < MINTOKENS) {
                    tokencount[n] -= iter->second;
                    typecount[n]--;
                    pruned++;
                    ngrams.erase(iter->first);                        
                } else {
                    ngramtokencount += iter->second;
                }
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << typecount[n] <<  " left (" << tokencount[n] << " tokens)" << endl;
    
       
       if (DOSKIPGRAMS) {       
           //prune skipgrams
           pruned = 0;
           for(unordered_map<EncSkipGram,uint32_t>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {     
               if (iter->first.n() == n) {                          
                    bool pruneskipgram = false;
                    if ( (iter->second < MINTOKENS) || (lastskipcontenthash[iter->first] != 1) ) {
                        pruneskipgram = true;
                        skiptokencount[n] -= iter->second;
                        skiptypecount[n]--;
                        pruned++;
                        skipgrams.erase(iter->first);
                    }
               }
           }
           cerr << "Pruned " << pruned << " skipgrams, " << skiptypecount[n] <<  " left (" << skiptokencount[n] << " tokens)" << endl;
           
        }
        ngramtokencount += tokencount[n];
        ngramtypecount += typecount[n];
        skipgramtokencount += skiptokencount[n];
        skipgramtypecount += skiptypecount[n];
    }

        
}


UnindexedPatternModel::UnindexedPatternModel(const string & filename, bool DOREVERSEINDEX) {    
    const bool DEBUG = true;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;    
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    readfile(filename);
}

void UnindexedPatternModel::readngram(std::istream * f, const EncNGram & ngram) {
    ngramtypecount++;
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    ngrams[ngram] = count;
    ngramtokencount += count;    
    if (ngram.n() > MAXLENGTH) MAXLENGTH = ngram.n();
    tokencount[ngram.n()] += count;
}



void UnindexedPatternModel::readskipgram(std::istream * f, const EncSkipGram & skipgram) {
    skipgramtypecount++;
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count            
    skipgrams[skipgram] = count; //assign
    skipgramtokencount += count;
    skiptokencount[skipgram.n()] += count;            
}

void UnindexedPatternModel::writengram(std::ostream * f, const EncNGram & ngram) {
    const uint32_t c = ngrams[ngram];
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                                       
}


void UnindexedPatternModel::writengrams(std::ostream * f) {    
    const char czero = 0;
    for(unordered_map<EncNGram,uint32_t>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengram(f, iter->first);       
    }   
}


void UnindexedPatternModel::writeskipgram(std::ostream * f, const EncSkipGram & skipgram) {
    const uint32_t c  = skipgrams[skipgram];
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                         
}


void UnindexedPatternModel::writeskipgrams(std::ostream * f) {             
    for(unordered_map<EncSkipGram,uint32_t>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {                               
        iter->first.writeasbinary(f);        
        writeskipgram(f, iter->first);
    }     
}


bool UnindexedPatternModel::exists(const EncAnyGram* key) const {    
    if (key->gapcount() == 0) {
        return (ngrams.count(*( (EncNGram*) key) ) > 0);
    } else {
        return (skipgrams.count(*( (EncSkipGram*) key) ) > 0);
    }
    return false;
}

const EncAnyGram* UnindexedPatternModel::getkey(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<EncNGram,uint32_t >::iterator iter = ngrams.find(*( (EncNGram*) key) );
        if (iter != ngrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<EncSkipGram,uint32_t >::iterator iter = skipgrams.find(*( (EncSkipGram*) key) );
        if (iter != skipgrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }        
    }
}




int UnindexedPatternModel::count(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ];
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return skipgrams[*( (EncSkipGram*) key) ];   
    }
    return 0;
}


double UnindexedPatternModel::freq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ] / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)] / tokens();
    }
    return 0;
}


double UnindexedPatternModel::relfreq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ] / tokencount[key->n()];
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)] / skiptokencount[key->n()];
    }
    return 0;
}


void UnindexedPatternModel::decode(ClassDecoder & classdecoder, ostream *NGRAMSOUT, ostream *SKIPGRAMSOUT) {
    const int grandtotal = ngramtokencount + skipgramtokencount;   

    for(unordered_map<EncNGram,uint32_t>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
       const double freq1 = (double) iter->second / tokencount[iter->first.n()];
       const double freq2 = (double) iter->second / ngramtokencount;
       const double freq3 = (double) iter->second / grandtotal;
       const EncNGram ngram = iter->first;
        *NGRAMSOUT << (int) ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq1 << '\t' << freq2 << '\t' << freq3;
        *NGRAMSOUT << endl;
    }
   

   if (SKIPGRAMSOUT != NULL) {
       for(unordered_map<EncSkipGram,uint32_t>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
           const double freq1 = (double) iter->second / skiptokencount[iter->first.n()]; 
           const double freq2 = (double) iter->second / skipgramtokencount;           
           const double freq3 = (double) iter->second / grandtotal;                          
           const EncSkipGram skipgram = iter->first;                              
           *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';          
           *SKIPGRAMSOUT << endl;
       }
    }

}


double SkipGramData::entropy() const {
    double entropy = 0;
    for(unordered_map<EncSkipGram,NGramData>::const_iterator iter = skipcontent.begin(); iter != skipcontent.end(); iter++ ) {
      double p = iter->second.count() / (double) _count;
      entropy += p * log2(p);
    }    
    return -1 * entropy;
}


std::set<CorpusReference> SkipGramData::get_refs() const {
    std::set<CorpusReference> s;
    for (std::unordered_map<EncSkipGram,NGramData>::const_iterator iter = skipcontent.begin(); iter != skipcontent.end(); iter++) {
        for (std::set<CorpusReference>::const_iterator iter2 = iter->second.refs.begin(); iter2 != iter->second.refs.end(); iter2++) {
            s.insert(*iter2);
        }
    }
    return s;
}




int intersection( set<CorpusReference> & a, set<CorpusReference> & b) {    
    //Jaccard co-occurrence    
    int intersectioncount = 0;        
    set<CorpusReference>::iterator iter_a = a.begin();    
    set<CorpusReference>::iterator iter_b = b.begin();    
    while ((iter_a !=a.end()) && (iter_b!=b.end())) {
        if (*iter_a < *iter_b) { 
            iter_a++;
        } else if (*iter_b < *iter_a) {
            iter_b++;
        } else {  //equal
            intersectioncount++;
            iter_a++;
            iter_b++;
        }
    }
    return intersectioncount;
}



GraphPatternModel::GraphPatternModel(IndexedPatternModel * model, bool DOPARENTS, bool DOCHILDREN, bool DOXCOUNT) {
    this->model = model;
    this->DOPARENTS = DOPARENTS;
    this->DOCHILDREN = DOCHILDREN;   
    this->DOXCOUNT = DOXCOUNT;     
    DELETEMODEL = false;
    
    cerr << "Computing relations on n-grams" << endl;
    for(std::unordered_map<EncNGram,NGramData >::iterator iter = model->ngrams.begin(); iter != model->ngrams.end(); iter++ ) {
        const EncNGram * ngram = &(iter->first);
        vector<EncNGram*> subngrams;
        ngram->subngrams(subngrams);
        for (vector<EncNGram*>::iterator iter2 = subngrams.begin(); iter2 != subngrams.end(); iter2++) {                
            const EncAnyGram * subngram = model->getkey(*iter2);
            if (subngram != NULL) {
                //subgram exists, add relation:
                if (DOCHILDREN)
                 rel_subsumption_children[ngram].insert(subngram);
                                    
                //reverse:
                if ((DOPARENTS) || (DOXCOUNT))
                 rel_subsumption_parents[subngram].insert(ngram);
            }
        }                  
    }
    
    

    cerr << "Computing relations on skip-grams" << endl;
    for(std::unordered_map<EncSkipGram,SkipGramData >::iterator iter = model->skipgrams.begin(); iter != model->skipgrams.end(); iter++ ) {        
        const EncSkipGram * skipgram = &(iter->first);
        vector<EncNGram*> parts;
        skipgram->parts(parts);
        for (vector<EncNGram*>::iterator iter2 = parts.begin(); iter2 != parts.end(); iter2++) {
            const EncNGram * ngram = *iter2;
            vector<EncNGram*> subngrams;
            ngram->subngrams(subngrams);
            for (vector<EncNGram*>::iterator iter3 = subngrams.begin(); iter3 != subngrams.end(); iter3++) {
                //subgram exists, add relation:
                const EncAnyGram * subngram = model->getkey(*iter3);
                if (subngram != NULL) {
                    if (DOCHILDREN)
                     rel_subsumption_children[skipgram].insert(subngram);
                     
                    if ((DOPARENTS) || (DOXCOUNT))
                     rel_subsumption_parents[subngram].insert(skipgram);
                }
            }
            delete ngram;
        }          
    }

    if (DOXCOUNT) {
        cerr << "Computing exclusive count" << endl;
        for (std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> >::const_iterator iter = rel_subsumption_parents.begin(); iter != rel_subsumption_parents.end(); iter++) {
            const EncAnyGram* anygram = iter->first;
            data_xcount[anygram] = xcount(anygram);
        }        
        if (!DOPARENTS) rel_subsumption_parents.clear();
    }
 
    if (DOCHILDREN) cerr << "Child subsumption relations: " << rel_subsumption_children.size() << endl;
    if (DOPARENTS) cerr << "Parent subsumption relations: " << rel_subsumption_parents.size() << endl;
    if (DOXCOUNT) cerr << "Exclusive count: " << data_xcount.size() << endl;
    
}


int GraphPatternModel::xcount(const EncAnyGram* anygram) {
    const AnyGramData* data = model->getdata(anygram);
    if (data == NULL) throw "xcount: No such anygram";
    set<CorpusReference> allrefs = data->get_refs();
        
    //compute union of all parent references
    set<CorpusReference> parentrefs;
    for (std::unordered_set<const EncAnyGram*>::iterator iter = rel_subsumption_parents[anygram].begin(); iter != rel_subsumption_parents[anygram].end(); iter++) {
        const AnyGramData * parentdata = model->getdata(*iter);
        if (parentdata != NULL) {
            parentrefs = parentdata->get_refs();
            for (set<CorpusReference>::iterator iter2 = parentrefs.begin(); iter2 != parentrefs.end(); iter2++) {
                CorpusReference r = *iter2; 
                //const CorpusReference * r = &(*(iter2)); //ugly, I know
                parentrefs.insert(r);
            }
        }        
    }
        
    //compute: union - intersection
    return allrefs.size() - intersection(allrefs, parentrefs);
}


void GraphPatternModel::readrelations(std::istream * in, const EncAnyGram * anygram, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relationhash, bool ignore) {
    uint16_t count;
    in->read((char*) &count,  sizeof(uint16_t));
    char gapcount;
    for (int i = 0; i < count; i++) {                        
       in->read(&gapcount, sizeof(char));
       if (gapcount == 0) {
        EncNGram ngram = EncNGram(in);
        if (!ignore) relationhash[model->getkey(anygram)].insert(model->getkey((EncAnyGram*) &ngram));
       } else {
        EncSkipGram skipgram = EncSkipGram( in, gapcount);
        if (!ignore) relationhash[model->getkey(anygram)].insert(model->getkey((EncAnyGram*) &skipgram));
       }           
    }    
}

void GraphPatternModel::writerelations(std::ostream * out,const EncAnyGram * anygram, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relationhash) {
    uint16_t count = relationhash[model->getkey(anygram)].size();
    out->write((char*) &count, sizeof(uint16_t));
    char gapcount;
    for (int i = 0; i < count; i++) {                        
        for (unordered_set<const EncAnyGram*>::iterator iter = relationhash[model->getkey(anygram)].begin(); iter != relationhash[model->getkey(anygram)].end(); iter++) {
            const EncAnyGram * anygram2 = *iter;
            anygram2->writeasbinary(out);
        }
    }    
}


void GraphPatternModel::readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes) {
    in->read((char*) &HASPARENTS,  sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASCHILDREN, sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASXCOUNT, sizeof(bool)); //1 byte, not 1 bit
}

void GraphPatternModel::writeheader(std::ostream * out) {
    out->write((char*) &DOPARENTS,  sizeof(bool)); //1 byte, not 1 bit
    out->write((char*) &DOCHILDREN, sizeof(bool)); //1 byte, not 1 bit
    out->write((char*) &DOXCOUNT, sizeof(bool)); //1 byte, not 1 bit
}

void GraphPatternModel::readngram(std::istream * in, const EncNGram & ngram) {
	model->readngram(in, ngram);   
    if (HASXCOUNT) {
        uint32_t _xcount;
        in->read((char*) &_xcount, sizeof(uint32_t));
        data_xcount[model->getkey((const EncAnyGram*) &ngram)] = _xcount;
    }
    if (HASPARENTS) readrelations(in, (const EncAnyGram*) &ngram, rel_subsumption_parents, !DOPARENTS);

    if (HASCHILDREN) 
        readrelations(in, (const EncAnyGram*) &ngram, rel_subsumption_children, !DOCHILDREN);    
}

void GraphPatternModel::readskipgram(std::istream * in, const EncSkipGram & skipgram) {
	model->readskipgram(in,skipgram);
    if (HASXCOUNT) {
        uint32_t _xcount;
        in->read((char*) &_xcount, sizeof(uint32_t));
        data_xcount[model->getkey((const EncAnyGram*) &skipgram)] = _xcount;
    }
    if (HASPARENTS) {
        readrelations(in, (const EncAnyGram*) &skipgram, rel_subsumption_parents, !DOPARENTS);
    }
    if (HASCHILDREN) 
        readrelations(in, (const EncAnyGram*) &skipgram, rel_subsumption_children, !DOCHILDREN);        
}


void GraphPatternModel::writengram(std::ostream * out, const EncNGram & ngram) {
	model->writengram(out,ngram);
    if (DOXCOUNT) {
        uint32_t _xcount = xcount((const EncAnyGram*) &ngram);
        out->write( (char*) &_xcount, sizeof(uint32_t));
    }
    
    if (DOPARENTS)
          writerelations(out, (const EncAnyGram*) &ngram, rel_subsumption_parents);
        
    if (DOCHILDREN)  
        writerelations(out, (const EncAnyGram*) &ngram, rel_subsumption_children);
        
            
}

void GraphPatternModel::writeskipgram(std::ostream * out, const EncSkipGram & skipgram) {
	model->writeskipgram(out,skipgram);    
    if (DOXCOUNT) {
        uint32_t _xcount = xcount((const EncAnyGram*) &skipgram);
        out->write( (char*) &_xcount, sizeof(uint32_t));
    }
    if (DOPARENTS)  writerelations(out, (const EncAnyGram*) &skipgram, rel_subsumption_parents);
        
    if (DOCHILDREN)  
        writerelations(out, (const EncAnyGram*) &skipgram, rel_subsumption_children);        
}



void GraphPatternModel::writengrams(std::ostream * f) {    
    const char czero = 0;
    for(unordered_map<EncNGram,NGramData>::iterator iter = model->ngrams.begin(); iter !=  model->ngrams.end(); iter++ ) {        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengram(f, iter->first);       
    }   
}


void GraphPatternModel::writeskipgrams(std::ostream * f) {                 
    for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = model->skipgrams.begin(); iter !=  model->skipgrams.end(); iter++ ) {                                
        iter->first.writeasbinary(f);        
        writeskipgram(f, iter->first);
    }     
}

GraphPatternModel::~GraphPatternModel() {
    if (DELETEMODEL) delete model;
}


void GraphPatternModel::decode(ClassDecoder & classdecoder, ostream *NGRAMSOUT, ostream *SKIPGRAMSOUT) {
    const int grandtotal = model->tokens();

    cerr << "Outputting n-grams" << endl;    
    for(unordered_map<const EncNGram,NGramData>::const_iterator iter = model->ngrams.begin(); iter != model->ngrams.end(); iter++ ) {
       const double freq1 = (double) iter->second.count() / model->tokencount[iter->first.n()];
       const double freq2 = (double) iter->second.count() / model->ngramtokencount;
       const double freq3 = (double) iter->second.count() / grandtotal;
       const EncAnyGram * ngram = &iter->first;       
        *NGRAMSOUT << (int) ngram->n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram->decode(classdecoder) << '\t' << iter->second.count() << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';
        if (DOXCOUNT) {            
            if (data_xcount.count(ngram) ) {
                int xc = data_xcount[ngram];
                double xratio = data_xcount[ngram] / (double) iter->second.count() ;
                *NGRAMSOUT << xc << '\t' << xratio << '\t';                
            } else {            
                *NGRAMSOUT << iter->second.count() << '\t' << 1.0 << '\t';
            }
        }
        for (set<CorpusReference>::iterator iter2 = iter->second.refs.begin() ; iter2 != iter->second.refs.end(); iter2++) {
            *NGRAMSOUT << iter2->sentence << ':' << (int) iter2->token << ' ';
        }                
        *NGRAMSOUT << endl;
    }
   

   if (SKIPGRAMSOUT != NULL) {
       cerr << "Outputting skip-grams" << endl;
       for(unordered_map<const EncSkipGram,SkipGramData>::const_iterator iter =  model->skipgrams.begin(); iter !=  model->skipgrams.end(); iter++ ) {
           const double freq1 = (double) iter->second.count() / model->skiptokencount[iter->first.n()]; 
           const double freq2 = (double) iter->second.count() / model->skipgramtokencount;           
           const double freq3 = (double) iter->second.count() / grandtotal;                          
           const EncAnyGram * skipgram = &iter->first;                              
           *SKIPGRAMSOUT << (int) skipgram->n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram->decode(classdecoder) << '\t' << iter->second.count() << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';
            if (DOXCOUNT) {            
                if (data_xcount.count( skipgram ) ) {
                    int xc = data_xcount[skipgram];
                    double xratio = data_xcount[skipgram ] / (double) iter->second.count() ;
                    *SKIPGRAMSOUT << xc << '\t' << xratio << '\t';                
                } else {            
                    *SKIPGRAMSOUT << iter->second.count() << '\t' << 1.0 << '\t';
                }
            }           	
           const int skiptypes = iter->second.skipcontent.size();               
           const double entropy = iter->second.entropy();
           *SKIPGRAMSOUT << skiptypes << '\t' << iter->second.count() << '\t' << entropy << '\t';
            for(unordered_map<const EncSkipGram,NGramData>::const_iterator iter2 = iter->second.skipcontent.begin(); iter2 != iter->second.skipcontent.end(); iter2++ ) {
                *SKIPGRAMSOUT << iter2->first.decode(classdecoder) << '|' << iter2->second.count() << '|';
                for (set<CorpusReference>::iterator iter3 = iter2->second.refs.begin() ; iter3 != iter2->second.refs.end(); iter3++) {
                    *SKIPGRAMSOUT << iter3->sentence << ':' << (int) iter3->token;        
                    *SKIPGRAMSOUT << ',';
                    //if (iter3 != iter2->second.refs.end() - 1) 
                }
                //MAYBE TODO: output references?
            }
           *SKIPGRAMSOUT << endl;
       }
    }

}

DoubleIndexedGraphPatternModel::DoubleIndexedGraphPatternModel(const std::string & filename) { //read a normal graph pattern model in another way optimised for Cooc alignment
	readfile(filename);
}


void DoubleIndexedGraphPatternModel::readheader(std::istream * in, uint64_t & totaltokens, uint64_t & totaltypes) {
    in->read((char*) &HASPARENTS,  sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASCHILDREN, sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASXCOUNT, sizeof(bool)); //1 byte, not 1 bit
    if (!HASXCOUNT) {
    	cerr << "WARNING: No Xcount data in graph model! Aligners may rely on this and not function without!" << endl;
    }
}

void DoubleIndexedGraphPatternModel::readrelations(std::istream * in, const EncAnyGram * anygram) {
	//Read and ignore relations (we don't care about them but we do need to read them)
    uint16_t count;
    in->read((char*) &count,  sizeof(uint16_t));
    char gapcount;
    for (int i = 0; i < count; i++) {   
       in->read(&gapcount, sizeof(char));
       if (gapcount == 0) {
        EncNGram ngram = EncNGram(in);
       } else {
        EncSkipGram skipgram = EncSkipGram( in, gapcount);
       }           
    }    
}


void DoubleIndexedGraphPatternModel::readngram(std::istream * in, const EncNGram & ngram) {
	ngrams[ngram]; //will create the ngram if it does not exist yet in the hash
	std::unordered_map<EncNGram,IndexCountData>::iterator iter = ngrams.find(ngram); //pointer to the ngram in the hash
	const EncAnyGram * anygram = &iter->first;

	ngramtypecount++;
    uint32_t count;
    in->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    ngramtokencount += count;    
    for (int j = 0; j < count; j++) {
        CorpusReference ref = CorpusReference(in); //read from file                
        reverseindex[ref.sentence].push_back(anygram);
    }	    
    if (HASXCOUNT) {
        uint32_t xcount;
        in->read((char*) &xcount, sizeof(uint32_t));
        ngrams[ngram].xcount = xcount;                   
    }
    if (HASPARENTS) readrelations(in, anygram); //read and ignore
    if (HASCHILDREN) readrelations(in, anygram);  //read and ignore
    //NOTE MAYBE TODO: make sure to update when GraphModel updates!    
}

void DoubleIndexedGraphPatternModel::readskipgram(std::istream * in, const EncSkipGram & skipgram) {
	skipgrams[skipgram]; //will create the ngram if it does not exist yet in the hash
	std::unordered_map<EncSkipGram,IndexCountData>::iterator iter = skipgrams.find(skipgram); //pointer to the skipgram in the hash
	const EncAnyGram * anygram = &iter->first;

    	  
    skipgramtypecount++;
    uint32_t count;
    in->read((char*) &count, sizeof(uint32_t)); //read occurrence count            
    skipgramtokencount += count;               
    uint32_t skipcontentcount;
    in->read((char*) &skipcontentcount, sizeof(uint32_t));   
    for (int j = 0; j < skipcontentcount; j++) {                                
        EncSkipGram skipcontent = EncSkipGram(in);  //also when !DOSKIPCONTENT, bytes have to be read
        in->read((char*) &count, sizeof(uint32_t)); //read occurrence count                
        for (int k = 0; k < count; k++) {
            CorpusReference ref = CorpusReference(in); //read from file                
        	reverseindex[ref.sentence].push_back(anygram);
        }        
    }    
    if (HASXCOUNT) {
        uint32_t xcount;
        in->read((char*) &xcount, sizeof(uint32_t));
        skipgrams[skipgram].xcount = xcount;    
    }
    if (HASPARENTS) readrelations(in, anygram);  //read and ignore
    if (HASCHILDREN) readrelations(in, anygram);  //read and ignore
    //NOTE MAYBE TODO: make sure to update when GraphModel updates!    
}



/*
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

*/
