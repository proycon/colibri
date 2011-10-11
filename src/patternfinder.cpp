#include <ngram.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <limits>
#include <iomanip> // contains setprecision()
#include <map>

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



int main( int argc, char *argv[] ) {
    
    if (argc != 3) {
        cerr << "Usage: patternfinder classfile encoded-corpus" << endl;
        exit(2);
    }
    
    const string classfile = argv[1];
    
    ClassDecoder classdecoder = ClassDecoder(classfile);
    
    string corpusfile = argv[2];
    
    cerr << "Processing " << classfile << endl;
    //int n = atoi(argv[2]);
    
    
    unsigned char line[1024];    
    vector<freqlist> ngrams;
    vector<skipgrammap> skipgrams;
    ngrams.push_back(freqlist()); //for 1-based indexing
    skipgrams.push_back(skipgrammap()); //for 1-based indexing
    const int MINTOKENS = 2;
    const int MAXLENGTH = 6;
    
    int tokencount[MAXLENGTH+1];
    int skiptokencount[MAXLENGTH+1];

    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;
        ngrams.push_back(freqlist());
        skipgrams.push_back(skipgrammap());
        
        
        
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

            
                for (int j = 0; j < gaps.size(); j++) {
                    const int begin = gaps[j].first;  
                    const int length = gaps[j].second;
                    
                    
                    
                    //Don't count skipgram if its consecutive subparts are not in the ngram lists
                    EncNGram * skipgram_preskip = ngram->slice(0,begin);
                    const int preskip_n = skipgram_preskip->n();
                    if ((preskip_n > 1) && (!(ngrams[preskip_n-1].count(*skipgram_preskip)))) {
                        delete skipgram_preskip;
                        continue;
                    }
                    
                    EncNGram * skipgram_postskip = ngram->slice(begin+length,ngram->n() - begin - length);
                    const int postskip_n = skipgram_postskip->n();
                    if ((postskip_n > 1) && (!(ngrams[postskip_n-1].count(*skipgram_postskip)))) {
                        delete skipgram_preskip;
                        delete skipgram_postskip;
                        continue;

                    }
                    
                    EncSingleSkipGram skipgram = EncSingleSkipGram(*skipgram_preskip, *skipgram_postskip);
                    EncNGram * skip = ngram->slice(begin,length);
                    
                    skipgrams[n][skipgram].count++;
                    skipgrams[n][skipgram].skips[*skip] += 1;
                    skiptokencount[n]++;
                    
                    delete skip;
                    delete skipgram_preskip;
                    delete skipgram_postskip;
                }
                
                delete ngram;                 
            }            
            
        };

       cerr << "Found " << tokencount[n] << " " << n << "-grams and " << skiptokencount[n] << " skipgrams" << endl;

       //prune n-grams
       int pruned = 0;
       int ngramtotal = 0;
       for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
            if (iter->second <= MINTOKENS) {
                pruned += ngrams[n].erase(iter->first);
            } else {
                ngramtotal += iter->second;
            }
       }
       cerr << "Pruned " << pruned << " " << n << "-grams, " << (tokencount[n] - pruned) <<  " left" << endl;
       
       //prune skipgrams
       pruned = 0;
       int skipgramtotal = 0; //total tokens
       for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {
            if (iter->second.count <= MINTOKENS) {
                pruned += skipgrams[n].erase(iter->first);
            } else {
                skipgramtotal += iter->second.count;
            }
       }
       cerr << "Pruned " << pruned << " " << "skipgrams, " << (skiptokencount[n] - pruned) <<  " left" << endl;
       
       
       for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
           const double freq = (double) iter->second / ngramtotal;
           const EncNGram ngram = iter->first;
           //cout << "NGRAM " << "N=" << ngram.n() << " SIZE=" << (int) ngram.size() << " DECODED=" << ngram.decode(classdecoder) << endl;
           cout << ngram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq << endl;
       }
       
       for(skipgrammap::iterator iter = skipgrams[n].begin(); iter != skipgrams[n].end(); iter++ ) {
           const double freq = (double) iter->second.count / skipgramtotal;           
           const EncSingleSkipGram skipgram = iter->first;
           cout << skipgram.n() << '\t' << setprecision(numeric_limits<double>::digits10 + 1) << skipgram.decode(classdecoder) << '\t' << iter->second.count << '\t' << freq << endl;
       }
 
    }   

}

