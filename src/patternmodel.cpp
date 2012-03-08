#include "patternmodel.h"
#include "algorithms.h"
#include <limits>

using namespace std;




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
	int last = 0;
	//EncNGram lastngram;
	//EncSkipGram lastskipgram;
	
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
        
    f.read( (char*) &model_id, sizeof(uint64_t));        
    f.read( (char*) &totaltokens, sizeof(uint64_t));        
    f.read( (char*) &totaltypes, sizeof(uint64_t)); 
    
    readheader(&f);
    unsigned char check;    
    for (int i = 0; i < totaltypes; i++) {           
        char gapcount;
        if (DEBUG) cerr << "\t@" << i;
        f.read((char*) &check, sizeof(char));
        if (check != 0xff) {
        	cerr << "ERROR processing " + filename + " at construction " << i << " of " << totaltypes << ". Expected check-byte, got " << (int) check << endl;
        	f.read(&gapcount, sizeof(char));
        	cerr << "DEBUG: next byte should be gapcount, value=" << (int) gapcount << endl;
        	if (last == 1) {
        		cerr << "DEBUG: previous construction was a n-gram" << endl;
        	} else if (last == 2) {
        		cerr << "DEBUG: previous construction was a skip-gram" << endl;
        	} 
        	exit(13);        	
        }
        f.read(&gapcount, sizeof(char));
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            const EncNGram ngram = EncNGram(&f); //read from file            
            readngramdata(&f, ngram);      
            last = 1;
            //lastngram = ngram;    
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            const EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            readskipgramdata(&f, skipgram);
            last = 2;
            //lastskipgram = skipgram;
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


ModelQuerier::ModelQuerier() {
	for (int n = 1; n <= MAXN; n++) {
    	compute_multi_skips(gapconf[n], vector<pair<int,int> >(), n);
    }
}

std::vector<pair<const EncAnyGram*, CorpusReference> > ModelQuerier::getpatterns(const unsigned char * data, const unsigned char datasize, bool doskipgrams, uint32_t linenum) {
	
	std::vector<pair<const EncAnyGram*, CorpusReference> > patterns;

	//extract all patterns in an input string
	if (maxlength() > MAXN) {
       	cerr << "FATAL ERROR: Maximum n-gram size " << maxlength() << " exceeds the internal maximum MAXN=" << MAXN << endl;
       	exit(14);
    }    
    
	const int l = countwords(data, datasize);
	for (int begin = 0; begin <= l; begin++) {
		for (int length = 1; (length <= maxlength()) && (begin+length <= l);  length++) {
			EncNGram * ngram = getencngram(begin,length, data, datasize);
			const EncAnyGram * anygram =  ngram;			
			if (count(anygram) > 0) {
				patterns.push_back( make_pair<const EncAnyGram*,CorpusReference>(getkey(anygram), CorpusReference(linenum, (char) begin) ) ); //stores the actual pointer used by the model
				if (doskipgrams) {
					//TODO: make more efficient for complete models that are guaranteed not to prune sub-parts
				
					//iterate over all gap configurations
					const int n = ngram->n();
					for (size_t j = 0; j < gapconf[n].size(); j++) {
						//iterate over all gaps

                        vector<EncNGram*> subngrams;
                        vector<int> skipref;
                        bool initialskip = false;
                        bool finalskip = false;
                        int cursor = 0;
                        bool docount = true;    
                        int oc = 0;             

                        
                    	//Iterate over all gaps in this configuration    
                        for (size_t k = 0; k < gapconf[n][j].size(); k++) {                                                        
                            const int begin = gapconf[n][j][k].first;  
                            const int length = gapconf[n][j][k].second;                        
                            //cerr << begin << ';' << length << ';' << n << endl;
                            skipref.push_back( length); 
                            if (k == 0) {
                                initialskip = (begin == 0);
                            }                            
                            if (begin > cursor) {
                                EncNGram * subngram = ngram->slice(cursor,begin-cursor);
                                subngrams.push_back(subngram);                                                               
                                /*oc = count((const EncAnyGram *) subngram);
                                if (oc) oc = ngrams[*subngram].count();
                                if ((((oc == 0)) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) ))    {
                                    docount = false;
                                    break;
                                }*/
                            }
                            cursor = begin + length;
                        }
                        if (cursor < n) {
                            EncNGram * subngram = ngram->slice(cursor,n-cursor);
                            subngrams.push_back(subngram);
                            /*oc = ngrams.count(*subngram);
                            if (oc) oc = ngrams[*subngram].count();
                            if (((oc == 0)) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) )  {
                                docount = false;
                                break;
                            }*/
                        } else {
                            finalskip = true;
                        }
                        if (initialskip && finalskip && skipref.size() <= 1) docount = false; //the whole n-gram is a skip, discard
                        if (docount) {
                            EncSkipGram skipgram = EncSkipGram(subngrams, skipref, initialskip, finalskip);
                            const EncAnyGram * anygram2 = &skipgram;			
							if (count(anygram2) > 0) {
								patterns.push_back( make_pair<const EncAnyGram*,CorpusReference>(getkey(anygram2), CorpusReference(linenum, (char) begin)) ); //stores the actual pointer used by the model
							}			
						}
						for (size_t k = 0; k < subngrams.size(); k++) {
                            delete subngrams[k];
                        }
					}
				}
			} else {
				//if this pattern doesn't occur, we could break because longer patterns won't occur either in most models, but this is not always the case, so we don't break
				//TODO: make more efficient for complete models that are guaranteed not to prune sub-parts								
			}
			delete ngram;
		}
	}
	return patterns;
}

void ModelQuerier::querier(ClassEncoder & encoder, ClassDecoder & decoder, bool exact, bool repeat) {
	unsigned char buffer[65536];
	uint32_t linenum = 0;
    std::string line;
    do {
    	linenum++;
    	cout << linenum << ">> "; 
    	getline(cin,line);    	
    	if (!line.empty()) {
			int buffersize = encoder.encodestring(line, buffer);
			if (exact) {
				//TODO
			} else {    	
				vector<pair<const EncAnyGram*, CorpusReference> > patterns = getpatterns(buffer,buffersize, true, linenum);
				for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = patterns.begin(); iter != patterns.end(); iter++) {
					const EncAnyGram * anygram = iter->first;
					const CorpusReference ref = iter->second;
					outputinstance(anygram, ref, decoder);
				} 
			}
		}
    } while (!cin.eof() && (repeat));		
}


IndexedPatternModel::IndexedPatternModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  int MINSKIPTYPES, bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
    this->DOINITIALONLYSKIP = DOINITIALONLYSKIP;
    this->DOFINALONLYSKIP = DOFINALONLYSKIP;
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;

    const int BUFFERSIZE = 65536;
    unsigned char line[BUFFERSIZE];


    if (MAXLENGTH > MAXN) {
       	cerr << "FATAL ERROR: Maximum n-gram size " << MAXLENGTH << " exceeds the internal maximum MAXN=" << MAXN << endl;
       	exit(14);
    }

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
        if (!IN->good()) {
        	cerr << "ERROR: Unable to open file " << corpusfile << endl;
        	exit(5);
        }    
        vector<unsigned int> words;
        while (IN->good()) {
            const int linesize = readline(IN, line, BUFFERSIZE );            
                    
            sentence++;

            if (sentence % 10000 == 0) {
                cerr << "\t@" << sentence << endl;
            }
                            
            
            const int l = countwords(line, linesize);            
            if (l >= 256) {
                cerr << "WARNING: Sentence " << sentence << " exceeds maximum word-length 256, skipping! (" << linesize << " bytes)" << endl;
                continue;
           } else if (l == 0) {
            	cerr << "WARNING: Sentence " << sentence << " contains no words, skipping! (" << linesize << " bytes)" << endl;
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
                    if ((iter->second.count() < MINTOKENS) || ((iter->second.skipcontent.size() < MINSKIPTYPES)))  {
                        pruneskipgram = true;
                    } else {            
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


        
}


IndexedPatternModel::IndexedPatternModel(const string & filename) {    
    const bool DEBUG = true;
      
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    if (!filename.empty()) readfile(filename);
}

void IndexedPatternModel::readngramdata(std::istream * f, const EncNGram & ngram, bool ignore) {
    if (!ignore) ngramtypecount++;
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    for (int j = 0; j < count; j++) {
        CorpusReference ref = CorpusReference(f); //read from file
        if (!ignore) ngrams[ngram].refs.insert(ref);
        /*if (DOREVERSEINDEX) {
            bool found = false;
            for (int k = 0; k < ngram_reverse_index[index].size(); k++) if (ngram_reverse_index[index][k] == ngram) { found = true; break; };
            if (!found) ngram_reverse_index[index].push_back(ngram);
        }*/
               
    }
    if (!ignore)  {
    	ngramtokencount += ngrams[ngram].count();    
    	const char n = ngram.n();
    	if (n > MAXN) {
      		cerr << "FATAL ERROR: N-gram size " << n << " exceeds the internal maximum MAXN=" << MAXN << endl;
	       	exit(14);
    	}
    	if (n > MAXLENGTH) MAXLENGTH = n;    	
    	tokencount[n] += ngrams[ngram].count();
    }    
}



void IndexedPatternModel::readskipgramdata(std::istream * f, const EncSkipGram & skipgram, bool ignore) {
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count            
    if (!ignore) {
     	skipgramtypecount++;
    	skipgrams[skipgram]._count = count; //assign
    	skipgramtokencount += count;
    	const char n = skipgram.n();
    	if (n > MAXN) {
      		cerr << "FATAL ERROR: N-gram size " << n << " exceeds the internal maximum MAXN=" << MAXN << endl;
	       	exit(14);
    	}
    	skiptokencount[n] += count;
    }            
    uint32_t skipcontentcount;
    f->read((char*) &skipcontentcount, sizeof(uint32_t));   
    for (int j = 0; j < skipcontentcount; j++) {                                
        EncSkipGram skipcontent = EncSkipGram(f);  
        f->read((char*) &count, sizeof(uint32_t)); //read occurrence count                
        for (int k = 0; k < count; k++) {
            CorpusReference ref = CorpusReference(f); //read from file
            if (!ignore) skipgrams[skipgram].skipcontent[skipcontent].refs.insert(ref);
        }        
    }    
}

void IndexedPatternModel::writengramdata(std::ostream * f, const EncNGram & ngram) {
    const uint32_t c = ngrams[ngram].count();
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                                     
    for (set<CorpusReference>::iterator iter = ngrams[ngram].refs.begin(); iter != ngrams[ngram].refs.end(); iter++) {                    
        iter->writeasbinary(f);
    }                    
}


void IndexedPatternModel::writengrams(std::ostream * f) {
	const unsigned char check = 0xff;    
    const char czero = 0;
    for(unordered_map<EncNGram,NGramData>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
    	f->write((char*) &check, sizeof(char));        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengramdata(f, iter->first);       
    }   
}


void IndexedPatternModel::writeskipgramdata(std::ostream * f, const EncSkipGram & skipgram) {
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
	const unsigned char check = 0xff;                 
    for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
    	f->write((char*) &check, sizeof(char));                                
        iter->first.writeasbinary(f);        
        writeskipgramdata(f, iter->first);
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
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return (double) ngrams[*( (EncNGram*) key) ].count() / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return (double) skipgrams[ *( (EncSkipGram*) key)].count() / tokens();
    }
    return 0;
}


double IndexedPatternModel::relfreq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return (double) ngrams[*( (EncNGram*) key) ].count() / tokencount[key->n()];
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return (double) skipgrams[ *( (EncSkipGram*) key)].count() / skiptokencount[key->n()];
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

void IndexedPatternModel::outputinstance(const EncAnyGram * anygram, CorpusReference ref, ClassDecoder & decoder) {
	cout << ref.sentence << ':' << (int) ref.token << "\t" << anygram->decode(decoder) << "\t" << count(anygram) << "\t" << setprecision(numeric_limits<double>::digits10 + 1) << freq(anygram) << endl; 
}



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
        *NGRAMSOUT << "\t0\t-\t";
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
           *SKIPGRAMSOUT << skiptypes << '\t' << entropy << '\t';
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

void IndexedPatternModel::writeanygram(const EncAnyGram * anygram,std::ostream * out) {
	/*helper function to write an ngram or skipgram, provided it occurs in this model */
	const char zero = 0;
    const EncAnyGram * anygram2 = getkey(anygram); //find the anygram pointer as it appear in the model
    if (!anygram2->isskipgram()) {
    	//ngram
    	out->write( (char*) &zero, sizeof(char)); //for ngrams
    	//((EncNGram*) anygram2)->writeasbinary(out);
    	const EncNGram * tmp = (EncNGram*) anygram2;
    	tmp->writeasbinary(out);    	
    } else {
    	//skipgram
    	const EncSkipGram * tmp = (EncSkipGram*) anygram2;
    	tmp->writeasbinary(out);
    	//((EncSkipGram*) anygram2)->writeasbinary(out);
    }        
}         
        


UnindexedPatternModel::UnindexedPatternModel(const string & corpusfile, int MAXLENGTH, int MINTOKENS, bool DOSKIPGRAMS, int MINSKIPTOKENS,  bool DOINITIALONLYSKIP, bool DOFINALONLYSKIP) {
    
    this->MAXLENGTH = MAXLENGTH;
    this->MINTOKENS = MINTOKENS;
    this->DOSKIPGRAMS = DOSKIPGRAMS;
    this->MINSKIPTOKENS = MINSKIPTOKENS;
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
        if (!IN->good()) {
        	cerr << "ERROR: Unable to open file " << corpusfile << endl;
        	exit(5);
        }        
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
            } else if (l == 0) {
            	cerr << "WARNING: Sentence " << sentence << " contains no words, skipping!" << endl;
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


UnindexedPatternModel::UnindexedPatternModel(const string & filename) {    
    const bool DEBUG = true;    
    
    ngramtokencount = 0;
    skipgramtokencount = 0; 
    ngramtypecount = 0;
    skipgramtypecount = 0;
    MAXLENGTH = 0;
    
    readfile(filename);
}

void UnindexedPatternModel::readngramdata(std::istream * f, const EncNGram & ngram, bool ignore ) {    
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    if (!ignore) {
    	ngramtypecount++;
    	ngrams[ngram] = count;
    	ngramtokencount += count;
        if (ngram.n() > MAXLENGTH) MAXLENGTH = ngram.n();
    	tokencount[ngram.n()] += count;
    }
}



void UnindexedPatternModel::readskipgramdata(std::istream * f, const EncSkipGram & skipgram, bool ignore) {    
    uint32_t count;
    f->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    if (!ignore) {
        skipgramtypecount++;    
    	skipgrams[skipgram] = count; //assign
	    skipgramtokencount += count;
	    skiptokencount[skipgram.n()] += count;
	}            
}

void UnindexedPatternModel::writengramdata(std::ostream * f, const EncNGram & ngram) {
    const uint32_t c = ngrams[ngram];
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                                       
}


void UnindexedPatternModel::writengrams(std::ostream * f) {    
    const unsigned char check = 0xff;
    const char czero = 0;
    for(unordered_map<EncNGram,uint32_t>::iterator iter = ngrams.begin(); iter != ngrams.end(); iter++ ) {
    	f->write((char*) &check, sizeof(char)); //check        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengramdata(f, iter->first);       
    }   
}


void UnindexedPatternModel::writeskipgramdata(std::ostream * f, const EncSkipGram & skipgram) {
    const uint32_t c  = skipgrams[skipgram];
    f->write( (char*) &c, sizeof(uint32_t) ); //occurrence count                         
}


void UnindexedPatternModel::writeskipgrams(std::ostream * f) {
	const unsigned char check = 0xff;             
    for(unordered_map<EncSkipGram,uint32_t>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
    	f->write((char*) &check, sizeof(char)); //check                               
        iter->first.writeasbinary(f);        
        writeskipgramdata(f, iter->first);
    }     
}

void UnindexedPatternModel::outputinstance(const EncAnyGram * anygram, CorpusReference ref, ClassDecoder & decoder) {
	cout << ref.sentence << ':' << (int) ref.token << "\t" << anygram->decode(decoder) << "\t" << count(anygram) << "\t" << setprecision(numeric_limits<double>::digits10 + 1) << freq(anygram) << endl; 
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
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return (double) ngrams[*( (EncNGram*) key) ] / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return (double) skipgrams[ *( (EncSkipGram*) key)] / tokens();
    }
    return 0;
}


double UnindexedPatternModel::relfreq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return (double) ngrams[*( (EncNGram*) key) ] / tokencount[key->n()];
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return (double) skipgrams[ *( (EncSkipGram*) key)] / skiptokencount[key->n()];
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

int transitivereduction(std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relations ) {
	int pruned = 0;
	//compute transitive closure by pruning
	for (std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> >::iterator iter = relations.begin(); iter != relations.end(); iter++) {
		const EncAnyGram * node = iter->first;		
		if (relations.count(node) > 0) {
			for (std::unordered_set<const EncAnyGram*>::iterator parentiter = relations[node].begin(); parentiter != relations[node].end(); parentiter++) {
				const EncAnyGram * parent = *parentiter;
				if (relations.count(parent) > 0) {
					for (std::unordered_set<const EncAnyGram*>::iterator grandparentiter = relations[parent].begin(); grandparentiter != relations[parent].end(); grandparentiter++) {
						const EncAnyGram * grandparent = *grandparentiter;
						if (relations[node].count(grandparent) > 0) {
							relations[node].erase(grandparent);
							pruned++;
						}						
					}
				} 
			}
		}
	}
	return pruned;
}

int GraphPatternModel::transitivereduction() {
	int pruned = 0;
	if (DOCHILDREN) pruned += ::transitivereduction(rel_subsumption_children);
    if (DOPARENTS) pruned += ::transitivereduction(rel_subsumption_parents);
    return pruned;	
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
        if ((!ignore) && (secondpass)) relationhash[model->getkey(anygram)].insert(model->getkey((EncAnyGram*) &ngram));
       } else {
        EncSkipGram skipgram = EncSkipGram( in, gapcount);
        if ((!ignore) && (secondpass)) relationhash[model->getkey(anygram)].insert(model->getkey((EncAnyGram*) &skipgram));
       }           
    }    
}

void GraphPatternModel::writerelations(std::ostream * out,const EncAnyGram * anygram, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > & relationhash) {
	const char zero = 0;
	unordered_set<const EncAnyGram*> * relations = &relationhash[model->getkey(anygram)];
	
    uint16_t count = relations->size();  // size() doesn't correspond to actual iterations???
    /*for (unordered_set<const EncAnyGram*>::iterator iter = relationhash[model->getkey(anygram)].begin(); iter != relationhash[model->getkey(anygram)].end(); iter++) {
    	count++;    
    }*/
    out->write((char*) &count, sizeof(uint16_t));
    char gapcount;
    int c = 0;                           
    for (unordered_set<const EncAnyGram*>::iterator iter = relations->begin(); iter != relations->end(); iter++) {
    	c++;
		model->writeanygram(*iter, out);
        /*const EncAnyGram * anygram2 = model->getkey(*iter);
        if (!anygram2->isskipgram()) out->write( (char*) &zero, sizeof(char)); //for ngrams         
        ((EncSkipGram*) anygram2)->writeasbinary(out);*/
    }
    //sanity check:
    if (c != count) {
    	cerr << "INTERNAL ERROR: GraphPatternModel::writerelations: Sanity check failed, wrote " << c << " constructions instead of expected " << count << endl;
    	cerr < "DEBUG: " << relations->size() << endl;        	
    	exit(13);
    }    
}


void GraphPatternModel::readheader(std::istream * in, bool ignore) {
    in->read((char*) &HASPARENTS,  sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASCHILDREN, sizeof(bool)); //1 byte, not 1 bit
    in->read((char*) &HASXCOUNT, sizeof(bool)); //1 byte, not 1 bit
}

void GraphPatternModel::writeheader(std::ostream * out) {
    out->write((char*) &DOPARENTS,  sizeof(bool)); //1 byte, not 1 bit
    out->write((char*) &DOCHILDREN, sizeof(bool)); //1 byte, not 1 bit
    out->write((char*) &DOXCOUNT, sizeof(bool)); //1 byte, not 1 bit
}

void GraphPatternModel::readngramdata(std::istream * in, const EncNGram & ngram, bool ignore) {
	model->readngramdata(in, ngram, ignore || secondpass);
    if (HASXCOUNT) {
        uint32_t _xcount;
        in->read((char*) &_xcount, sizeof(uint32_t));
        if ((DOXCOUNT) && (secondpass) && (!ignore)) {
        	data_xcount[model->getkey((const EncAnyGram*) &ngram)] = _xcount;
        	
        }
    }
    if (HASPARENTS) readrelations(in, (const EncAnyGram*) &ngram, rel_subsumption_parents, (ignore || !secondpass || !DOPARENTS));

    if (HASCHILDREN) 
        readrelations(in, (const EncAnyGram*) &ngram, rel_subsumption_children, (ignore || !secondpass || !DOCHILDREN));
}

void GraphPatternModel::readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore) {
	model->readskipgramdata(in,skipgram, ignore || secondpass);
    if (HASXCOUNT) {
        uint32_t _xcount;
        in->read((char*) &_xcount, sizeof(uint32_t));
        if ((DOXCOUNT) && (secondpass) && (!ignore)) data_xcount[model->getkey((const EncAnyGram*) &skipgram)] = _xcount;
    }
    if (HASPARENTS) {
        readrelations(in, (const EncAnyGram*) &skipgram, rel_subsumption_parents, (ignore || !secondpass || !DOPARENTS));
    }
    if (HASCHILDREN) 
        readrelations(in, (const EncAnyGram*) &skipgram, rel_subsumption_children, (ignore || !secondpass || !DOCHILDREN));        
}


void GraphPatternModel::writengramdata(std::ostream * out, const EncNGram & ngram) {
	model->writengramdata(out,ngram);
    if (DOXCOUNT) {
        uint32_t _xcount = xcount((const EncAnyGram*) &ngram);
        out->write( (char*) &_xcount, sizeof(uint32_t));
    }
    
    if (DOPARENTS)
          writerelations(out, (const EncAnyGram*) &ngram, rel_subsumption_parents);
        
    if (DOCHILDREN)  
        writerelations(out, (const EncAnyGram*) &ngram, rel_subsumption_children);
        
            
}

void GraphPatternModel::writeskipgramdata(std::ostream * out, const EncSkipGram & skipgram) {
	model->writeskipgramdata(out,skipgram);    
    if (DOXCOUNT) {
        uint32_t _xcount = xcount((const EncAnyGram*) &skipgram);
        out->write( (char*) &_xcount, sizeof(uint32_t));
    }
    if (DOPARENTS)  writerelations(out, (const EncAnyGram*) &skipgram, rel_subsumption_parents);
        
    if (DOCHILDREN)  
        writerelations(out, (const EncAnyGram*) &skipgram, rel_subsumption_children);        
}



void GraphPatternModel::writengrams(std::ostream * f) {
	const unsigned char check = 0xff;    
    const char czero = 0;
    for(unordered_map<EncNGram,NGramData>::iterator iter = model->ngrams.begin(); iter !=  model->ngrams.end(); iter++ ) {
    	f->write((char*) &check, sizeof(char)); //check        
        f->write(&czero, sizeof(char)); //gapcount, always zero for ngrams
        iter->first.writeasbinary(f);
        writengramdata(f, iter->first);       
    }   
}


void GraphPatternModel::writeskipgrams(std::ostream * f) {
	const char check = 0xff;                 
    for(unordered_map<EncSkipGram,SkipGramData>::iterator iter = model->skipgrams.begin(); iter !=  model->skipgrams.end(); iter++ ) {
    	f->write(&check, sizeof(char)); //gapcount, always zero for ngrams                                
        iter->first.writeasbinary(f);        
        writeskipgramdata(f, iter->first);
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


/****************************************************************************************************************************************/

SelectivePatternModel::SelectivePatternModel(const std::string & filename,  bool DOFORWARDINDEX,  bool DOREVERSEINDEX, bool DOXCOUNT, int COUNTTHRESHOLD, double FREQTHRESHOLD, double XCOUNTRATIOTHRESHOLD, int XCOUNTTHRESHOLD, bool DOSKIPGRAMS,  int MINLENGTH, int MAXLENGTH, bool DOPARENTS, bool DOCHILDREN, AlignConstraintInterface * alignconstrain, bool alignconstrainsource) { //read a normal graph pattern model in another way optimised for Cooc alignment
	this->DOFORWARDINDEX = DOFORWARDINDEX;
	this->DOREVERSEINDEX = DOREVERSEINDEX;
	this->DOXCOUNT = DOXCOUNT;
	this->COUNTTHRESHOLD = COUNTTHRESHOLD;
	this->FREQTHRESHOLD = FREQTHRESHOLD;
	this->XCOUNTRATIOTHRESHOLD = XCOUNTRATIOTHRESHOLD;
	this->XCOUNTTHRESHOLD = XCOUNTTHRESHOLD;
	this->DOSKIPGRAMS = DOSKIPGRAMS;
	this->MINLENGTH = MINLENGTH;
	this->MAXLENGTH = MAXLENGTH;
	this->DOPARENTS = DOPARENTS;
	this->alignconstrain = alignconstrain;
	this->alignconstrainsource = alignconstrainsource;
	

	ngramtokencount = 0;
	skipgramtokencount = 0;
	ngramtypecount = 0;
	skipgramtypecount = 0;
	ignoredtypes = 0;
	ignoredtokens = 0;
	        
    secondpass = false;       
	readfile(filename);	
	secondpass = DOPARENTS;
	if (secondpass) {
		readfile(filename);
	}
}


void SelectivePatternModel::readheader(std::istream * in, bool ignore) {
	if (model_id == GRAPHPATTERNMODEL) { 
    	in->read((char*) &HASPARENTS,  sizeof(bool)); //1 byte, not 1 bit
	    in->read((char*) &HASCHILDREN, sizeof(bool)); //1 byte, not 1 bit
	    in->read((char*) &HASXCOUNT, sizeof(bool)); //1 byte, not 1 bit
	} else {
		HASPARENTS = false;
		HASCHILDREN = false;
		HASXCOUNT = false;
	}	
	if ((model_id == UNINDEXEDPATTERNMODEL) && (DOFORWARDINDEX || DOREVERSEINDEX)) 
	  cerr << "WARNING!!! You opted to load indexes but you the model you are loading is unindexed!" << endl;
	if (!HASXCOUNT && DOXCOUNT) 
	  cerr << "WARNING!!! You opted to load exclusive count data but the model you are loading does not contain this data! (Make sure to load a graphmodel generated with exclusive count data)" << endl;
	if (!HASPARENTS && DOPARENTS) 
	  cerr << "WARNING!!! You opted to load parent relation data but the model you are loading does not contain this data! (Make sure to load a graphmodel generated with parent relation data)" << endl;
	if (!HASCHILDREN && DOCHILDREN) 
	  cerr << "WARNING!!! You opted to load children relation data but the model you are loading does not contain this data! (Make sure to load a graphmodel generated with children relation data)" << endl;
}

/*void SelectivePatternModel::readrelations(std::istream * in) {
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
}*/

int SelectivePatternModel::transitivereduction() {
	int pruned = 0;
	if ((HASCHILDREN) && (DOCHILDREN)) ::transitivereduction(rel_subsumption_children);
    if ((HASPARENTS) && (DOPARENTS)) ::transitivereduction(rel_subsumption_parents);
    return pruned;	
}

void SelectivePatternModel::readrelations(std::istream * in, const EncAnyGram * anygram, std::unordered_map<const EncAnyGram*,std::unordered_set<const EncAnyGram*> > * relationhash, bool ignore) {
    uint16_t count;
    in->read((char*) &count,  sizeof(uint16_t));
    char gapcount;
    for (int i = 0; i < count; i++) {                        
       in->read(&gapcount, sizeof(char));
       if (gapcount == 0) {
        EncNGram ngram = EncNGram(in);
        if ((!ignore) && (secondpass) && (anygram != NULL) && (relationhash != NULL)) (*relationhash)[getkey(anygram)].insert(getkey((EncAnyGram*) &ngram));
       } else {
        EncSkipGram skipgram = EncSkipGram( in, gapcount);
        if ((!ignore) && (secondpass) && (anygram != NULL) && (relationhash != NULL)) (*relationhash)[getkey(anygram)].insert(getkey((EncAnyGram*) &skipgram));
       }           
    }    
}


void SelectivePatternModel::readngramdata(std::istream * in, const EncNGram & ngram, bool ignore) {
	//NOTE MAYBE TODO: make sure to update when GraphModel updates!

	if ((ngram.n() < MINLENGTH) || (ngram.n() > MAXLENGTH)) ignore = true;

	//READ STAGE -- nothing permanently stored yet

    uint32_t count;
    vector<uint32_t> index;
    uint32_t xcount = 0;     
    
    in->read((char*) &count, sizeof(uint32_t)); //read occurrence count		
	if (model_id != UNINDEXEDPATTERNMODEL) {	
		for (int j = 0; j < count; j++) {
		    CorpusReference ref = CorpusReference(in); //read from file
		    if (!ignore) index.push_back(ref.sentence);         
		}
	}    
    if (HASXCOUNT) in->read((char*) &xcount, sizeof(uint32_t)); //read, process later
    
    if (secondpass) {    	
    	const EncAnyGram * anygram = getkey(&ngram);	
    	if (HASPARENTS) readrelations(in, anygram, &rel_subsumption_parents);
    	if (HASCHILDREN) readrelations(in, anygram, &rel_subsumption_children);
    } else {
    	if (HASPARENTS) readrelations(in); //read and ignore
    	if (HASCHILDREN) readrelations(in);  //read and ignore
    }
    
    //THRESHOLD CHECK STAGE - deciding whether to ignore based on unreached thresholds
    
	if ((count < COUNTTHRESHOLD) || ((double) count / totaltypes < FREQTHRESHOLD)) ignore = true;
	if ((DOXCOUNT) && (HASXCOUNT) && ((xcount < XCOUNTTHRESHOLD) || (( (double) xcount / count) < XCOUNTRATIOTHRESHOLD) ) ) ignore = true;


	//ALIGNMENT MODEL CONSTRAINTS
	if (alignconstrain != NULL) {
		if (alignconstrainsource) {
			if (alignconstrain->getsourcekey((const EncAnyGram*) &ngram) == NULL) {
				ignore = true;
			}
		} else {
			if (alignconstrain->gettargetkey((const EncAnyGram*) &ngram) == NULL) {
				ignore = true;
			}
		}
	}

    //STORAGE STAGE
    if (!ignore) {
		ngrams[ngram]; //will create the ngram if it does not exist yet in the hash
		std::unordered_map<EncNGram,IndexCountData>::iterator iter = ngrams.find(ngram); //pointer to the ngram in the hash
		const EncAnyGram * anygram = &iter->first;	
		ngramtypecount++;
		ngramtokencount += count;
		ngrams[ngram].count = count; // 0 if n/a
		ngrams[ngram].xcount = xcount; // 0 if n/a
		if (!index.empty()) {			
			for (vector<uint32_t>::iterator iter = index.begin(); iter != index.end(); iter++) {
				ngrams[ngram].sentences.insert(*iter);
		    	reverseindex[*iter].push_back(anygram);
			}
		}	
	} else {
		ignoredtypes++;
		ignoredtokens += count;
	}        	        
}



void SelectivePatternModel::readskipgramdata(std::istream * in, const EncSkipGram & skipgram, bool ignore) {
	//NOTE MAYBE TODO: make sure to update when GraphModel updates!

	if ( (!DOSKIPGRAMS) || (skipgram.n() < MINLENGTH) || (skipgram.n() > MAXLENGTH) ) ignore = true;
	 
	//READ STAGE -- nothing permanently stored yet
	
    uint32_t count;
    uint32_t skipcontentcount;
    uint32_t xcount = 0;
    vector<uint32_t> index;
    
    in->read((char*) &count, sizeof(uint32_t)); //read occurrence count
    
    if (model_id != UNINDEXEDPATTERNMODEL) {
		in->read((char*) &skipcontentcount, sizeof(uint32_t));   
		for (int j = 0; j < skipcontentcount; j++) {                                
		    EncSkipGram skipcontent = EncSkipGram(in);  
		    in->read((char*) &count, sizeof(uint32_t)); //read occurrence count                
		    for (int k = 0; k < count; k++) {
		        CorpusReference ref = CorpusReference(in); //read from file       
		        if (!ignore) index.push_back(ref.sentence);            
		    }        
		}    
	}
    if (HASXCOUNT) in->read((char*) &xcount, sizeof(uint32_t));     

    if (secondpass) {
    	const EncAnyGram * anygram = getkey(&skipgram);	
    	if (HASPARENTS) readrelations(in, anygram, &rel_subsumption_parents);
    	if (HASCHILDREN) readrelations(in, anygram, &rel_subsumption_children);
    } else {
    	if (HASPARENTS) readrelations(in); //read and ignore
    	if (HASCHILDREN) readrelations(in);  //read and ignore
    }
            
    //THRESHOLD CHECK STAGE - deciding whether to ignore based on unreached thresholds
    
	if ((count < COUNTTHRESHOLD) || ((double) count / totaltypes < FREQTHRESHOLD)) ignore = true;
	if ((DOXCOUNT) && (HASXCOUNT) && ((xcount < XCOUNTTHRESHOLD) || ((double) (xcount / count) < XCOUNTRATIOTHRESHOLD) ) ) ignore = true;
	
	//ALIGNMENT MODEL CONSTRAINTS
	if (alignconstrain != NULL) {
		if (alignconstrainsource) {
			if (alignconstrain->getsourcekey((const EncAnyGram*) &skipgram) == NULL) {
				ignore = true;
			}
		} else {
			if (alignconstrain->gettargetkey((const EncAnyGram*) &skipgram) == NULL) {
				ignore = true;
			}
		}
	}	
	
     //STORAGE STAGE
    if (!ignore) {
    	skipgrams[skipgram]; //will create the ngram if it does not exist yet in the hash
		std::unordered_map<EncSkipGram,IndexCountData>::iterator iter = skipgrams.find(skipgram); //pointer to the skipgram in the hash
		const EncAnyGram * anygram = &iter->first;
    	skipgramtypecount++;            
    	skipgramtokencount += count;
    	skipgrams[skipgram].count = count; // 0 if n/a
    	skipgrams[skipgram].xcount = xcount;  // 0 if n/a
		if (!index.empty()) {			
			for (vector<uint32_t>::iterator iter = index.begin(); iter != index.end(); iter++) {
		    	skipgrams[skipgram].sentences.insert(*iter);
		    	reverseindex[*iter].push_back(anygram);
			}
		}	    	
    } else {
		ignoredtypes++;
		ignoredtokens += count;
	}   
    
        
}


int SelectivePatternModel::count(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key)) > 0) return ngrams[*( (EncNGram*) key) ].count;
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return skipgrams[*( (EncSkipGram*) key) ].count;   
    }
    return 0;
}

int SelectivePatternModel::xcount(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key)) > 0) return ngrams[*( (EncNGram*) key) ].xcount;
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return skipgrams[*( (EncSkipGram*) key) ].xcount;   
    }
    return 0;
}

double SelectivePatternModel::xcountratio(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key)) > 0) return (double) ngrams[*( (EncNGram*) key) ].xcount / ngrams[*( (EncNGram*) key) ].count;
    } else {
        if (skipgrams.count(*( (EncSkipGram*) key) ) > 0) return (double) skipgrams[*( (EncSkipGram*) key) ].xcount / skipgrams[*( (EncSkipGram*) key) ].count;   
    }
    return 0;
}



double SelectivePatternModel::freq(const EncAnyGram* key) {    
    if (key->gapcount() == 0) {        
        if (ngrams.count(*( (EncNGram*) key) ) > 0) return (double) ngrams[*( (EncNGram*) key) ].count / tokens();
    } else {
        if (skipgrams.count( *( (EncSkipGram*) key)) > 0) return (double) skipgrams[ *( (EncSkipGram*) key)].count / tokens();
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

void SelectivePatternModel::outputinstance(const EncAnyGram * anygram, CorpusReference ref, ClassDecoder & decoder) {
	cout << ref.sentence << ':' << (int) ref.token << "\t" << anygram->decode(decoder) << "\t" << count(anygram) << "\t";	
	cout << "\t" << setprecision(numeric_limits<double>::digits10 + 1) << freq(anygram);
	if (DOXCOUNT) cout << "\t" << xcount(anygram) << "\t" << xcountratio(anygram) << endl;
	cout << endl;	
}


const EncAnyGram* SelectivePatternModel::getkey(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<const EncNGram,IndexCountData >::iterator iter = ngrams.find(*( (EncNGram*) key) );
        if (iter != ngrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<const EncSkipGram,IndexCountData >::iterator iter = skipgrams.find(*( (EncSkipGram*) key) );
        if (iter != skipgrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }        
    }
}
