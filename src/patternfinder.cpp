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

using namespace std;

typedef unordered_map<EncNGram,int> freqlist;

class skipgramdata {
   public:
    int count;
    freqlist skips;
    skipgramdata() {
        count = 0;
    }
};

typedef unordered_map<EncSingleSkipGram,skipgramdata> skipgrammap;

long total(freqlist& f) { 
    long s = 0;
    for(freqlist::iterator iter = f.begin(); iter != f.end(); iter++ ) {
      s += iter->second;
    }
    return s;
}


double compute_entropy(freqlist & data) {
    double entropy = 0;
    for(freqlist::iterator iter = data.begin(); iter != data.end(); iter++ ) {
      entropy += iter->second * -1 * log2(iter->second);
    }    
    return -1 * entropy;
}


vector< pair<int,int> > get_consecutive_gaps(const int n) {
    vector< pair<int,int> > gaps;
    int begin = 1;
    while (begin < n) {
        int length = (n - 1) - begin;
        while (length > 0) {
            pair<int,int> gap = make_pair(begin, length);
            gaps.push_back(gap);
            length--;
        }
        begin++;
    }      
    return gaps;
}


void usage() {
    cerr << "Syntax: patternfinder -c classfile -f encoded-corpus" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-t <number>      Token threshold: n-grams and skipgrams occuring less than this will be pruned (default: 2)" << endl;
    cerr << "\t-l <number>      Maximum n-gram/skipgram length (in words, default: 9)" << endl;
    cerr << "\t-s               Compute skip-grams" << endl;    
    cerr << "\t-T <number>      Skip threshold: only skip content that occurs at least x times will be considered (default: 1) " << endl;
    cerr << "\t-S <number>      Skip type threshold: only skipgrams with x possible types for the skip will be considered, otherwise the skipgram will be pruned  (default: 2)" << endl;
    cerr << "\t-i               Compute index" << endl;
        
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
    
    char c;    
    while ((c = getopt(argc, argv, "c:f:t:T:S:l:o:si")) != -1)
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
        cerr << "Need to specify -c classfile and -f corpusfile" << endl;
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
    
    
    unsigned char line[1024];    
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
        unordered_map<EncSingleSkipGram,set<int> > skipgram_index;
        
        int linenum = 0;
        tokencount[n] = 0;
        skiptokencount[n] = 0;

        const vector< pair<int,int> > gaps = get_consecutive_gaps(n);
        
        
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

            
                for (size_t j = 0; j < gaps.size(); j++) {
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
                    
                    EncSingleSkipGram skipgram = EncSingleSkipGram(*skipgram_preskip, *skipgram_postskip, n);
                    EncNGram * skip = ngram->slice(begin,length);
                    
                    
                    bool docount = true;
                    
                    if ((MINSKIPTOKENS > 1) && (MINSKIPTOKENS >= MINTOKENS)) {
                        const int skip_n = skipgram_postskip->n();
                        if (!(ngrams[skip_n].count(*skip))) {
                            docount = false;
                        } else {
                            docount = (ngrams[skip_n-1][*skip] >= MINSKIPTOKENS);
                        } 
                    }
                    
                    if (docount) {
                        skipgrams[n][skipgram].count++;
                        skipgrams[n][skipgram].skips[*skip] += 1;
                        skiptokencount[n]++;
                        
                        if (DOINDEX) skipgram_index[skipgram].insert(linenum);

                    }
                    
                    delete skip;
                    delete skipgram_preskip;
                    delete skipgram_postskip;
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
            if (iter->second < MINTOKENS) {
                tokencount[n] -= iter->second;
                pruned++;
                ngrams[n].erase(iter->first);        
                if (DOINDEX) ngram_index.erase(iter->first);        
            } else {
                ngramtotal += iter->second;
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << ngrams[n].size() <<  " left (" << tokencount[n] << " tokens)" << endl;
       
    
       
       if (DOSKIPGRAMS) {       
           //prune skipgrams
           pruned = 0;
           for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {
                if ((iter->second.count < MINTOKENS) || (iter->second.skips.size() < MINSKIPTYPES)) {
                    //prune skipgram
                    skiptokencount[n] -= iter->second.count;
                    pruned++;
                    skipgrams[n].erase(iter->first);                    
                    if (DOINDEX) skipgram_index.erase(iter->first);
                } else {
                    for(freqlist::iterator iter2 = iter->second.skips.begin(); iter2 != iter->second.skips.end(); iter2++ ) {
                        if (iter2->second < MINSKIPTOKENS) {
                            //prune skip
                            iter->second.count -= iter2->second;
                            iter->second.skips.erase(iter2->first); 
                        }
                    }
                    skipgramtotal += iter->second.count;
                }
           }
           cerr << "Pruned " << pruned << " skipgrams, " << skipgrams[n].size() <<  " left (" << skiptokencount[n] << " tokens)" << endl;
           
        }
        
        if (DOINDEX) {
            cerr << "Writing ngram index";

            for(unordered_map<EncNGram,set<int>>::iterator iter = ngram_index.begin(); iter != ngram_index.end(); iter++ ) {
                const EncNGram ngram = iter->first;
                *NGRAMINDEX << ngram.decode(classdecoder);
                for (set<int>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                    *NGRAMINDEX << ' ' << *iter2;
                }
            }
            *NGRAMINDEX << endl;
            
            if (DOSKIPGRAMS) {
                cerr << "Writing skipgram index";

                for(unordered_map<EncSingleSkipGram,set<int>>::iterator iter = skipgram_index.begin(); iter != skipgram_index.end(); iter++ ) {
                    const EncSingleSkipGram skipgram = iter->first;
                    *SKIPGRAMINDEX << skipgram.decode(classdecoder);
                    for (set<int>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                        *SKIPGRAMINDEX << ' ' << *iter2;
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
               const int skiptypes = iter->second.skips.size();
               const double entropy = compute_entropy(iter->second.skips);
               const EncSingleSkipGram skipgram = iter->first;                              
               *SKIPGRAMSOUT << (int) skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count << '\t' << freq1 << '\t' << freq2 << '\t' << freq3 << '\t' << skiptypes << '\t' << iter->second.count << '\t' << entropy << endl;
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

