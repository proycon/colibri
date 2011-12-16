#include "patternmodel.h"
#include "algorithms.h"

using namespace std;

CorpusReference::CorpusReference(uint32_t sentence, unsigned char token) {
    this->sentence = sentence;
    this->token = token;
}

void CorpusReference::writeasbinary(ostream * out) const {  
    out->write( (char*) &sentence, sizeof(uint32_t) );
    out->write( (char*) &token, sizeof(unsigned char) );
}


void NGramData::writeasbinary(ostream * out) const {
    for (set<CorpusReference>::iterator iter = refs.begin(); iter != refs.end(); iter++) {
        iter->writeasbinary(out);
    }
}




EncGramIndexedModel::EncGramIndexedModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  int MINSKIPTYPES, bool DOREVERSEINDEX, bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
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
            if (l > 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping!";
                continue;
            }
            
            if (linesize > 0) //no { on purpose! applies to next for loop
            for (unsigned char i = 0; ((i < l - n + 1) && (i < 256)); i++) {                
                EncNGram * ngram = getencngram(i,n, line, linesize);  
                                                

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
                                if (oc) oc = ngrams[*subngram].count();
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
                            if (oc) oc = ngrams[*subngram].count();
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
                            if (skipgrams[skipgram].count() == 0) skiptypecount[n]++;                            
                            skipgrams[skipgram]._count++;
                            skipgrams[skipgram].skipcontent[skipcontent].refs.insert(ref);
                            skiptokencount[n]++;                            
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


EncGramIndexedModel::EncGramIndexedModel(const string & filename, bool DOREVERSEINDEX) {
    const bool DEBUG = false;
    this->DOREVERSEINDEX = DOREVERSEINDEX;
    this->DOSKIPCONTENT = DOSKIPCONTENT;    
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    
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
            int count;
            f.read((char*) &count, sizeof(uint32_t)); //read occurrence count
            for (int j = 0; j < count; j++) {
                int index;
                CorpusReference ref = CorpusReference(&f); //read from file
                ngrams[ngram].refs.insert(ref);
                /*if (DOREVERSEINDEX) {
                    bool found = false;
                    for (int k = 0; k < ngram_reverse_index[index].size(); k++) if (ngram_reverse_index[index][k] == ngram) { found = true; break; };
                    if (!found) ngram_reverse_index[index].push_back(ngram);
                }*/
                       
            }
            ngramtokencount += ngrams[ngram].count();
            tokencount[ngram.n()] += ngrams[ngram].count();
        } else {
            //SKIPGRAM
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            skipgramtypecount++;            
            EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            int count;
            f.read((char*) &count, sizeof(uint32_t)); //read occurrence count            
            skipgrams[skipgram]._count = count; //assign
            skipgramtokencount += count;
            skiptokencount[skipgram.n()] += count;            
            int skipcontentcount;
            f.read((char*) &skipcontentcount, sizeof(int));   
            if (DEBUG) cerr << "\tcontentcount=" << (int) skipcontentcount;
            for (int j = 0; j < skipcontentcount; j++) {                                
                EncSkipGram skipcontent = EncSkipGram(&f);  //also when !DOSKIPCONTENT, bytes have to be read
                f.read((char*) &count, sizeof(int)); //read occurrence count                
                for (int k = 0; k < count; k++) {
                    CorpusReference ref = CorpusReference(&f); //read from file
                    skipgrams[skipgram].skipcontent[skipcontent].insert(ref);
                }
            }
            /*
                if (DOREVERSEINDEX) {
                    bool found = false;
                    for (int k = 0; k < skipgram_reverse_index[index].size(); k++) if (skipgram_reverse_index[index][k] == skipgram) { found = true; break; };
                    if (!found) skipgram_reverse_index[index].push_back(skipgram);                    
                }
            */
        }
        if (DEBUG)  cerr << endl;      //DEBUG  
    }
    f.close();
}



void EncGramIndexedModel::save(const std::string & filename) {
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    
    const char czero = 0;
    const int zero = 0;
    //const unsigned char check = 255;
    const unsigned long totaltokens = tokens();
    const unsigned long totaltypes = types();
    
    f.write( (char*) &totaltokens, sizeof(unsigned long) );
    f.write( (char*) &totaltypes, sizeof(unsigned long) );
    for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {        
        //f.write((char*) &check, sizeof(char)); 
        f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(&f);
        const uint32_t c = iter->second.count();
        f.write( (char*) &c, sizeof(uint32_t) ); //occurrence count                                     
        for (set<CorpusReference>::iterator iter2 = iter->second.refs.begin(); iter2 != iter->second.refs.end(); iter2++) {                    
            iter2->writeasbinary(&f);
        }                
    }    

    for (int n = 1; n <= MAXLENGTH; n++) {                
        for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {                            
            //f.write((char*)  &check, sizeof(char)); 
            iter->first.writeasbinary(&f);
            const uint32_t c  = iter->second.count();
            f.write( (char*) &c, sizeof(uint32_t) ); //occurrence count                         
            const uint32_t nrofskips = (uint32_t) iter->second.skipcontent.size();                
            f.write( (char*) &nrofskips, sizeof(uint32_t) );
            for(unordered_map<EncSkipGram,NGramData>::iterator iter2 = iter->second.skipcontent.begin(); iter2 != iter->second.skipcontent.end(); iter2++ ) {                    
                iter2->first.writeasbinary(&f);
                iter2->second.writeasbinary(&f);
            }
        }            
    }
    f.close();    
}


bool EncGramIndexedModel::exists(const EncAnyGram* key) const {    
    if (key->gapcount() == 0) {
        return (ngrams.count(*( (EncNGram*) key) ) > 0);
    } else {
        return (skipgrams.count(*( (EncSkipGram*) key) ) > 0);
    }
    return false;
}



int EncGramIndexedModel::count(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ];
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return skipgrams[*( (EncSkipGram*) key) ].count;   
    }
    return 0;
}


double EncGramIndexedModel::freq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ] / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)].count / tokens();
    }
    return 0;
}


double EncGramIndexedModel::relfreq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return ngrams[*( (EncNGram*) key) ] / tokencount[n];
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return skipgrams[ *( (EncSkipGram*) key)].count / skiptokencount[n];
    }
    return 0;
}


/*set<int> * EncGramIndexedModel::index(const EncAnyGram* key) {
    if (key->gapcount() == 0) {        
        if (ngram_index.count(*( (EncNGram*) key) ) > 0) return &ngram_index[*( (EncNGram*) key) ];
    } else {
        if (skipgram_index.count( *( (EncSkipGram*) key)) > 0) return &skipgram_index[ *( (EncSkipGram*) key)];
    }
}*/




std::set<int> EncGramIndexedModel::reverse_index_keys() {
    set<int> keys;
    for (unordered_map<int,vector<EncNGram> >::iterator iter = ngram_reverse_index.begin(); iter != ngram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }   
    for (unordered_map<int,vector<EncSkipGram> >::iterator iter = skipgram_reverse_index.begin(); iter != skipgram_reverse_index.end(); iter++) {
        keys.insert(iter->first);
    }    
    return keys;
}


int EncGramIndexedModel::reverse_index_size(const int i) {
    int s = 0;
    if (ngram_reverse_index.count(i)) s += ngram_reverse_index[i].size();
    if (skipgram_reverse_index.count(i)) s += skipgram_reverse_index[i].size();
    return s;
    
}

int EncGramIndexedModel::reverse_index_size() {
    return ngram_reverse_index.size() + skipgram_reverse_index.size();
}

bool EncGramIndexedModel::reverse_index_haskey(const int i) const {
    return ((ngram_reverse_index.count(i) > 0) || (skipgram_reverse_index.count(i) > 0));
}


vector<EncAnyGram*> EncGramIndexedModel::reverse_index(const int i) {
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


EncAnyGram* EncGramIndexedModel::get_reverse_index_item(const int key, const int i) {
    const int s = ngram_reverse_index[key].size();
    if (i < s) {                
        vector<EncNGram>::iterator iter = ngram_reverse_index[key].begin() + i;
        return &(*iter);
    } else {
        vector<EncSkipGram>::iterator iter = skipgram_reverse_index[key].begin() + (i - s);
        return &(*iter);
    }    
}

void EncGramIndexedModel::decode(ClassDecoder & classdecoder, ostream *NGRAMSOUT, ostream *SKIPGRAMSOUT) {
    const int grandtotal = ngramtokencount + skipgramtokencount;   

    for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
       const double freq1 = (double) iter->second / tokencount[iter->first.n()];
       const double freq2 = (double) iter->second / ngramtokencount;
       const double freq3 = (double) iter->second / grandtotal;
       const EncNGram ngram = iter->first;
        *NGRAMSOUT << (int) ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq1 << '\t' << freq2 << '\t' << freq3;
        *NGRAMSOUT << '\t';
        for (set<CorpusReference>::iterator iter2 = iter->second.refs.begin() ; iter2 != iter->second.refs.end(); iter2++) {
            *NGRAMSOUT << iter2->sentence << ':' << iter2->token << ' ';
        }                
        *NGRAMSOUT << endl;
    }
   

   if (SKIPGRAMSOUT != NULL) {
       for(unordered_map<EncNGram,SkipGramData>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
           const double freq1 = (double) iter->second.count() / skiptokencount[iter->first.n()]; 
           const double freq2 = (double) iter->second.count() / skipgramtokencount;           
           const double freq3 = (double) iter->second.count() / grandtotal;                          
           const EncSkipGram skipgram = iter->first;                              
           *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t';
           const int skiptypes = iter->second.skips.size();               
           const double entropy = iter->second.entropy();
           *SKIPGRAMSOUT << skiptypes << '\t' << iter->second.count << '\t' << entropy << '\t';
            for(unordered_map<EncSkipGram,NGramData>::iterator iter2 = iter->second.skipcontent.begin(); iter2 != iter->second.skipcontent.end(); iter2++ ) {
                *SKIPGRAMSOUT << iter2->first.decode(classdecoder) << '|' << iter2->second.count() << '|';
                for (set<CorpusReference>::iterator iter3 = iter2->second.refs.begin() ; iter3 != iter2->second.refs.end(); iter3++) {
                    *SKIPGRAMSOUT << iter3->sentence << ':' << iter3->token;        
                    if (iter3 != iter2->second.refs.ends() - 1) *SKIPGRAMSOUT << ',';
                }
                //MAYBE TODO: output references?
            }
           *SKIPGRAMSOUT << endl;
       }
    }

}




double SkipGramData::entropy() {
    double entropy = 0;
    for(unordered_map<EncSkipGram,NGramData>::iterator iter = skipcontent.begin(); iter != skipcontent.end(); iter++ ) {
      double p = iter->second.count() / (double) _count;
      entropy += p * log2(p);
    }    
    return -1 * entropy;
}




/*
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






EMAlignmentModel::EMAlignmentModel(EncGramIndexedModel & sourcemodel, EncGramIndexedModel & targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD) {
    int round = 0;    
    unsigned long c;
    double totaldivergence = 0;
    double prevavdivergence = 0;
    bool converged = false;
    
    set<int> reverseindexkeys = sourcemodel.reverse_index_keys();
            
    do {       
        round++; 
        c = 0;
        cerr << "  EM Round " << round << "... ";
        //use reverse index to iterate over all sentences
        for (set<int>::iterator iter = reverseindexkeys.begin(); iter != reverseindexkeys.end(); iter++) {
            const int key = *iter;        
            const int sourcegrams_size =  sourcemodel.reverse_index_size(key);
            //vector<EncAnyGram*> sourcegrams = sourcemodel.reverse_index(key);
            const int targetgrams_size =  targetmodel.reverse_index_size(key);
            
            if (targetgrams_size > 0 ) { //target model contains sentence?               
                //vector<EncAnyGram*> targetgrams = targetmodel.reverse_index(key);  
                cerr << key << ":" << sourcegrams_size << "x" << targetgrams_size << " ";

                //compute sentencetotal for normalisation later in count step, sum_s(p(t|s))                                      
                unordered_map<const EncAnyGram*, double> sentencetotal;                
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        sentencetotal[targetgram] += alignprob[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration
                    }
                }               
                
                //collect counts (for evidence that a targetgram is aligned to a sourcegram)
                std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
                unordered_map<const EncAnyGram*, double> total;
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        const double countvalue = alignprob[sourcegram][targetgram] / sentencetotal[targetgram];
                        count[sourcegram][targetgram] += countvalue;
                        total[targetgram] += countvalue;
                    }                    
                }
                
                 
                double prevtransprob;
                //estimate new probabilities (normalised count is the new estimated probability)
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);
                        
                        prevtransprob = alignprob[sourcegram][targetgram];
                        const double newtransprob = (double) count[sourcegram][targetgram] / total[targetgram];
                        alignprob[sourcegram][targetgram] = newtransprob;                        
                        
                        //for computation of convergence
                        const double divergence = abs(newtransprob - prevtransprob);
                        cerr << " prevtransprob=" << prevtransprob << " ";
                        cerr << " newtransprob=" << newtransprob << " ";
                        cerr << " div=" << divergence << " ";
                    
                        totaldivergence += divergence;
                        c++;
                    }
                }
            }
        }
        const double avdivergence = (double) totaldivergence / c;
        converged = (((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence)) <= CONVERGEDTHRESHOLD);       
        prevavdivergence = avdivergence;
        cerr << " average divergence = " << avdivergence << ", alignprob size = " << alignprob.size() << endl;
    } while (!converged);    
}


double CoocAlignmentModel::cooc( set<CorpusReference> & sourceindex, set<CorpusReference> & targetindex) {    
    //Jaccard co-occurrence    
    int intersectioncount = 0;    
    
    set<int>::iterator sourceiter = sourceindex.begin();    
    set<int>::iterator targetiter = targetindex.begin();
    
    while ((sourceiter !=sourceindex.end()) && (targetiter!=targetindex.end())) {
        if (sourceiter->sentence < targetiter->sentence) { 
            sourceiter++;
        } else if (targetiter->sentence < sourceiter->sentence) {
            targetiter++;
        } else {  //equal
            intersectioncount++;
            sourceiter++;
            targetiter++;
        }
    }
    const int unioncount = (sourceindex.size() + targetindex.size()) - intersectioncount; 
    //cerr << "union=" << unioncount << " intersection=" << intersectioncount << " ";
    return (double) intersectioncount / unioncount;
}


int CoocAlignmentModel::compute(const EncAnyGram * sourcegram, set<int> & sourceindex, EncGramIndexedModel & targetmodel) {        
    int c = 0;
    double bestcooc = 0;
    for (set<int>::iterator iter2 = sourceindex.begin(); iter2 != sourceindex.end(); iter2++) {
        const int sentencenumber = *iter2;        
        const int targetgrams_size =  targetmodel.reverse_index_size(sentencenumber);    
        c += targetgrams_size;
        for (int i = 0; i < targetgrams_size; i++) {    
            const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(sentencenumber,i);  
            set<int> * targetindices;
            if (targetgram->gapcount() == 0) {
               targetindices = &targetmodel.ngram_index[*( (EncNGram*) targetgram)];
            } else {
               targetindices = &targetmodel.skipgram_index[*( (EncSkipGram*) targetgram)];
            }
            
            const double coocvalue = cooc(sourceindex, *targetindices);            
            if ((relthreshold) && (coocvalue > bestcooc)) bestcooc = coocvalue;            
            if (coocvalue >= absthreshold) {                
                alignprob[sourcegram][targetgram] = coocvalue;
            }
        }
        if (relthreshold) {
            //TODO: prune based on relative threshold
        }   
    }    
    return c;
}

CoocAlignmentModel::CoocAlignmentModel(EncGramModel & sourcemodel, EncGramModel & targetmodel, const double absthreshold, const double relthreshold) {
    this->absthreshold = absthreshold;
    this->relthreshold = relthreshold;
    int c = 0;
    for (unordered_map<EncNGram,NGramData >::iterator iter = sourcemodel.ngrams.begin();  iter != sourcemodel.ngrams.end(); iter++) {
        cerr << ++c << " ";
        compute(&iter->first, iter->second.refs, targetmodel);
    }    
    for (unordered_map<EncSkipGram,SkipGramData >::iterator iter = sourcemodel.skipgrams.begin();  iter != sourcemodel.skipgrams.end(); iter++) {
        cerr << ++c << " ";
        compute(&iter->first, iter->second.refs, targetmodel);
    }            
}
    
void AlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT) {
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignprob.begin(); iter != alignprob.end(); iter++) {
        const EncAnyGram* sourcegram = iter->first;
        *OUT << sourcegram->decode(sourceclassdecoder) << "\t";
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram* targetgram = iter2->first;            
            *OUT << targetgram->decode(targetclassdecoder) << "\t" << iter2->second << "\t";
        }
        *OUT << endl;
    }
}
*/
