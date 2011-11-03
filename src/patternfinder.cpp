#include <ngram.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <limits>
#include <iomanip> // contains setprecision()
#include <map>
#include <unistd.h>
#include <cmath>
#include <set>
#include <algorithms.h>

using namespace std;

typedef unordered_map<EncNGram,int> freqlist;

class skipgramdata {
   public:
    int count;
    unordered_map<char,freqlist> skips;
    skipgramdata() {
        count = 0;
    }
};

typedef unordered_map<EncSkipGram,skipgramdata> skipgrammap;

long total(freqlist& f) { 
    long s = 0;
    for(freqlist::iterator iter = f.begin(); iter != f.end(); iter++ ) {
      s += iter->second;
    }
    return s;
}


double compute_entropy(unordered_map<char,freqlist> & data, const int total) {
    double entropy = 0;
    for(unordered_map<char,freqlist>::iterator iter = data.begin(); iter != data.end(); iter++ ) {
        for(freqlist::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++ ) {
          double p = iter2->second / (double) total;
          //cout << setprecision(numeric_limits<double>::digits10 + 1) << iter->second << " / " << total << " = " << p << endl;
          entropy += p * log2(p);
        }    
    }
    return -1 * entropy;
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



void usage() {
    cerr << "Syntax: patternfinder -c classfile -f encoded-corpus" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-t <number>      Token threshold: n-grams and skipgrams occuring less than this will be pruned (default: 2)" << endl;
    cerr << "\t-l <number>      Maximum n-gram/skipgram length (in words, default: 9)" << endl;
    cerr << "\t-s               Compute skip-grams" << endl;    
    cerr << "\t-T <number>      Skip threshold: only skip content that occurs at least x times will be considered (default: 1) " << endl;
    cerr << "\t-S <number>      Skip type threshold: only skipgrams with x possible types for the skip will be considered, otherwise the skipgram will be pruned  (default: 2)" << endl;
    cerr << "\t-i               Compute and output index" << endl;
    //cerr << "\t-C               Compute and output compositionality" << endl;
    cerr << "\t-L               Output invididual skips" << endl;
    cerr << "\t-o <string>      Output prefix" << endl;
        
}


int main( int argc, char *argv[] ) {
    
    string classfile = "";
    string corpusfile = "";
    string outputprefix = "";
    
    int MINTOKENS = 2;
    int MINSKIPTOKENS = 1;
    unsigned int MINSKIPTYPES = 2;
    int MAXLENGTH = 8;
    bool DOSKIPGRAMS = false;
    bool DOINDEX = false;
    bool DOSKIPOUTPUT = false;
    //bool DOCOMPOSITIONALITY = false;
    
    char c;    
    while ((c = getopt(argc, argv, "c:f:t:T:S:l:o:siLh")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
            break;
        case 'f':
            corpusfile = optarg;
            break;        
        case 't':
            MINTOKENS = atoi(optarg);
            break;
        case 'T':
            MINSKIPTOKENS = atoi(optarg);            
            break;
        case 'S':
            MINSKIPTYPES = atoi(optarg);            
            break;
        case 'l':
            MAXLENGTH = atoi(optarg);            
            break;
        case 's':
            DOSKIPGRAMS = true;
            break;
        case 'i':
            DOINDEX = true;
            break;
        case 'L':
            DOSKIPOUTPUT = true;
            break;
        case 'o': 
            outputprefix = optarg;
            break;
        case '?':
            if (optopt == 'c') {
                cerr <<  "Option -" << optopt << " requires an argument." << endl;
            } else {
                cerr << "Unknown option: -" <<  optopt << endl;
            }
            
            return 1;
        default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
    
    if (classfile.empty() || corpusfile.empty()) {
        cerr << "ERROR: Need to specify -c classfile and -f corpusfile" << endl;
        usage();
        exit(2);
    }
    
    if (outputprefix.empty()) {
        outputprefix = corpusfile;
    }

    const string ngramoutputfile = outputprefix + ".ngrams";
    ofstream *NGRAMSOUT =  new ofstream( ngramoutputfile.c_str() );      
    const string skipgramoutputfile = outputprefix + ".skipgrams";
    const string ngramindexfile = outputprefix + ".ngrams.index";
    const string skipgramindexfile = outputprefix + ".skipgrams.index";
    ofstream *SKIPGRAMSOUT = NULL;
    ofstream *NGRAMINDEX = NULL;
    ofstream *SKIPGRAMINDEX = NULL;
    if (DOINDEX) NGRAMINDEX = new ofstream( ngramindexfile.c_str() );      
    if (DOSKIPGRAMS) {    
         SKIPGRAMSOUT = new ofstream( skipgramoutputfile.c_str() );      
         if (DOINDEX) SKIPGRAMINDEX = new ofstream( skipgramindexfile.c_str() );      
    }
    
    ClassDecoder classdecoder = ClassDecoder(classfile);

    
    cerr << "Processing " << classfile << endl;
    //int n = atoi(argv[2]);
    
    
    unsigned char line[65536];    
    vector<freqlist> ngrams;
    vector<skipgrammap> skipgrams;
    
    
    ngrams.push_back(freqlist()); //for 1-based indexing
    skipgrams.push_back(skipgrammap()); //for 1-based indexing
    
    int tokencount[MAXLENGTH+1];
    int skiptokencount[MAXLENGTH+1];
    
    unsigned long ngramtotal = 0;
    unsigned long skipgramtotal = 0;

    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;
        ngrams.push_back(freqlist());
        skipgrams.push_back(skipgrammap());
        
        unordered_map<EncNGram,set<int>> ngram_index;
        unordered_map<EncSkipGram,set<int> > skipgram_index;
        
        int linenum = 0;
        tokencount[n] = 0;
        skiptokencount[n] = 0;
        
                
        vector< vector< pair<int,int> > > gaps;
        compute_multi_skips(gaps, vector<pair<int,int> >(), n);
        //const vector< pair<int,int> > gaps = get_consecutive_gaps(n);

        
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
                
                //cout << "NGRAM("<<ngram.n()<<","<<(int)ngram.size() << "): <" << ngram.decode(classdecoder) << ">" << endl;
                              
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
                        
                        vector<EncNGram*> subngrams;
                        vector<int> skipref;
                        bool initialskip = false;
                        bool finalskip = false;
                        int cursor = 0;
                        bool docount = true;                    
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
                                const int oc = ngrams[subngram->n()].count(*subngram);
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
                            const int oc = ngrams[subngram->n()].count(*subngram);
                            if ((oc == 0) || ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS) && (oc < MINSKIPTOKENS)) )    {
                                docount = false;
                                break;
                            }
                        } else {
                            finalskip = true;
                        }
                        if (initialskip && finalskip && skipref.size() <= 1) docount = false; //the whole n-gram is a skip, discard
                        if (docount) {
                            /*if (n == 4) { //DEBUG
                                cerr << "SKIPREF=" << skipref.size() << " INITIAL=" << initialskip  << " FINAL=" << finalskip << " SUBNGRAMS=" << subngrams.size() << endl;
                                cerr << "-- INITIAL: " << initialskip << endl;
                                cerr << "-- FINAL: " << finalskip << endl;
                                for (char x = 0; x < (char) subngrams.size(); x++) {
                                    //_size += (subngrams[i])->size();
                                    //if (i < (char) subngrams.size() - 1) _size += 2; //double 0 byte delimiting subngrams
                                    cerr << "--SUBNGRAM-- n=" << (int) (subngrams[x])->n() << ",size=" << (int) (subngrams[x])->size()  << endl;        
                                    (subngrams[x])->out();
                                    cout << endl;
                                }
                                    
                            }*/
                                
                            EncSkipGram skipgram = EncSkipGram(subngrams, skipref, initialskip, finalskip);
                            for (size_t k = 0; k < gaps[j].size(); k++) {
                                const int begin = gaps[j][k].first;  
                                const int length = gaps[j][k].second;
                                EncNGram * skip = ngram->slice(begin,length);
                                skipgrams[n][skipgram].count++;
                                skipgrams[n][skipgram].skips[(char) k][*skip] += 1;
                                skiptokencount[n]++;
                                delete skip;
                            }
                            if (DOINDEX) skipgram_index[skipgram].insert(linenum);
                        }
                        //cleanup
                        for (size_t k = 0; k < subngrams.size(); k++) {
                            delete subngrams[k];
                        }
                    
                        /* OLD: 
                        const int begin = gaps[j].first;  
                        const int length = gaps[j].second;
                        
                        
                        
                        //Don't count skipgram if its consecutive subparts are not in the ngram lists
                        EncNGram * skipgram_preskip = ngram->slice(0,begin);
                        const int preskip_n = skipgram_preskip->n();
                        if (!(ngrams[preskip_n].count(*skipgram_preskip))) {
                            delete skipgram_preskip;
                            continue;
                        }
                        
                        EncNGram * skipgram_postskip = ngram->slice(begin+length,ngram->n() - begin - length);
                        const int postskip_n = skipgram_postskip->n();
                        if (!(ngrams[postskip_n].count(*skipgram_postskip))) {
                            delete skipgram_preskip;
                            delete skipgram_postskip;
                            continue;

                        }
                        
                        EncSkipGram skipgram = EncSkipGram(*skipgram_preskip, *skipgram_postskip, n);
                        
                        int skipcount = skipgram.skipcount;
                        
                        EncNGram * skip = ngram->slice(begin,length);
                        
                        
                        bool docount = true;
                        
                        if ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS)) {
                            const int skip_n = skip->n();
                            if (ngrams[skip_n].count(*skip) == 0) {
                                docount = false;
                            } else {
                                docount = (ngrams[skip_n][*skip] >= MINSKIPTOKENS);
                            } 
                        }
                        
                        if (docount) {
                            skipgrams[n][skipgram].count++;
                            skipgrams[n][skipgram].skips[skipcount][*skip] += 1;
                            skiptokencount[n]++;
                            
                            if (DOINDEX) skipgram_index[skipgram].insert(linenum);

                        }
                        
                        delete skip;
                        delete skipgram_preskip;
                        delete skipgram_postskip;
                        */
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
       //unordered_map<EncNGram,set<int>>::iterator iter2 = ngram_index.begin();
       for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
            //if (DOINDEX) iter2++;
            if (iter->second < MINTOKENS) {
                if (DOINDEX) ngram_index.erase( iter->first);        
                tokencount[n] -= iter->second;
                pruned++;
                ngrams[n].erase(iter->first);                        
            } else {
                ngramtotal += iter->second;
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << ngrams[n].size() <<  " left (" << tokencount[n] << " tokens)" << endl;
    
       
       if (DOSKIPGRAMS) {       
           //prune skipgrams
           pruned = 0;
           for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {                
                bool pruneskipgram = false;
                if ((iter->second.count < MINTOKENS) || ((char) iter->second.skips.size() < iter->first.skipcount))  {
                    pruneskipgram = true;
                } else {
                    for(unordered_map<char,freqlist>::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                        if ((iter2->second.size() < MINSKIPTYPES) || (iter->second.count < MINTOKENS))  {
                            //not enough types or overall tokens, prune entire skipgram
                            pruneskipgram = true;                        
                        } else  {               
                            bool prunedskip = false;
                            for(freqlist::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++ ) {
                                if (iter3->second < MINSKIPTOKENS) {
                                    //prune skip
                                    prunedskip = true;
                                    iter->second.count -= iter3->second;
                                    iter2->second.erase(iter3->first); 
                                }
                            }
                            if ( (prunedskip) && ( (iter2->second.size() < MINSKIPTYPES) || (iter->second.count < MINTOKENS) ) ) { //reevaluate
                                pruneskipgram = true;
                            } else {
                                skipgramtotal += iter->second.count;
                            }
                        }
                        if (pruneskipgram) break;
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
        
        if (DOINDEX) {
            cerr << "Writing ngram index" << endl;

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
                
            }
        }
    }


    
   
   
    cerr << "Writing results to file: " << ngramoutputfile;
    if (DOSKIPGRAMS) {
        cerr << " , " << skipgramoutputfile << endl;
    } else {
        cerr << endl;
    }
    const int grandtotal = ngramtotal + skipgramtotal;
   
    for (int n = 1; n <= MAXLENGTH; n++) {   
        for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
           const double freq1 = (double) iter->second / tokencount[n];
           const double freq2 = (double) iter->second / ngramtotal;
           const double freq3 = (double) iter->second / grandtotal;
           const EncNGram ngram = iter->first;
           //cout << "NGRAM " << "N=" << ngram.n() << " SIZE=" << (int) ngram.size() << " DECODED=" << ngram.decode(classdecoder) << endl;
           *NGRAMSOUT << (int) ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << endl;
        }
       
        if (DOSKIPGRAMS) {           
           
           for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {
               const double freq1 = (double) iter->second.count / skiptokencount[n]; 
               const double freq2 = (double) iter->second.count / skipgramtotal;           
               const double freq3 = (double) iter->second.count / grandtotal;                          
               int skiptypes = 0;
               for(unordered_map<char,freqlist>::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                   skiptypes += iter2->second.size();
               }               
               const double entropy = compute_entropy(iter->second.skips, iter->second.count);
               const EncSkipGram skipgram = iter->first;                              
               *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t' << skiptypes << '\t' << iter->second.count << '\t' << entropy << '\t';
               if (DOSKIPOUTPUT) {
                   for (char c = 0; c < 5; c++) {
                        if (iter->second.skips.count(c)) {
                            for(freqlist::iterator iter3 = iter->second.skips[c].begin(); iter3 != iter->second.skips[c].end(); iter3++ ) {
                                const EncNGram skipcontent = iter3->first;
                                *SKIPGRAMSOUT << skipcontent.decode(classdecoder) << '|' << iter3->second << '|';
                            }
                            *SKIPGRAMSOUT << '|'; //double pipe at end and delimiter for multiple skips
                        }
                   }                   
               }
               *SKIPGRAMSOUT << endl;
           }
           
        }
    }
    
    NGRAMSOUT->close();
    if (DOINDEX) NGRAMINDEX->close();
    if (DOSKIPGRAMS) {
        SKIPGRAMSOUT->close();
        if (DOINDEX) SKIPGRAMINDEX->close();
    }
    
      

    

}

