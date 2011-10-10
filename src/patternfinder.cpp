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
    ngrams.push_back(freqlist()); //for 1-based indexing
    vector<freqlist> skipgrams;
    const int MINTOKENS = 2;
    const int MAXLENGTH = 6;
    
    int tokencount[MAXLENGTH+1];
    

    for (int n = 1; n <= MAXLENGTH; n++) {
        cerr << "Counting " << n << "-grams" << endl;
        ngrams.push_back(freqlist());
        
        
        
        int linenum = 0;
        tokencount[n] = 0;
        int skiptokens = 0;
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
                EncNGram ngram = getencngram(i,n, line, linesize);  
                
                //cout << "NGRAM("<<ngram.n()<<","<<(int)ngram.size() << "): <" << ngram.decode(classdecoder) << ">" << endl;
                              
                if (n > 2) {                    
                    EncNGram subngram1 = ngram.slice(0, n - 1);
                    if (!(ngrams[n-1].count(subngram1))) continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                                                             
                    EncNGram subngram2 = ngram.slice(1, n - 1);
                    if (!(ngrams[n-1].count(subngram2))) continue; //if subngram does not exist (count==0), no need to count ngram, skip to next
                }
                                    
                ngrams[n][ngram] += 1;
                tokencount[n]++;
                //cout << ngram_s << endl;
            
                /*for (int j = 0; j < gaps.size(); j++) {
                    const int begin = gaps[j].first;  
                    const int length = gaps[j].second;
                    
                    
                    
                    const vector<string> skipgram_pregap = vector<string>(words.begin() + i, words.begin() + i + begin);
                    const string skipgram_pregap_s = join(skipgram_pregap);
                    
                    if (!(*ngrams).count(skipgram_pregap_s)) continue;

                    const vector<string> skipgram_postgap = vector<string>(words.begin() + i + begin + length, words.begin() + i + n);                    
                    const string skipgram_postgap_s = join(skipgram_postgap);
                    if (!(*ngrams)[n].count(skipgram_postgap_s)) continue;
                    
                    const string skipgram = skipgram_pregap_s + " | " + skipgram_postgap_s;
                    skipgrams[skipgram] += 1;
                    skiptokens++;
                    
                }*/
                
                 
            }            
            
        };

       cerr << "Found " << tokencount[n] << " " << n << "-grams" << endl; //"and " << skiptokens << " skipgrams" << endl;

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
       /*pruned = 0;
       int skipgramtotal = 0;
       for(unordered_map<string,int>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
            if (iter->second <= MINTOKENS) {
                pruned += skipgrams.erase(iter->first);
            } else {
                skipgramtotal += iter->second;
            }
       }
       cerr << "Pruned " << pruned << " " << "skipgrams, " << (skiptokens - pruned) <<  " left" << endl;*/
       
       
       for(freqlist::iterator iter = ngrams[n].begin(); iter != ngrams[n].end(); iter++ ) {
           const double freq = (double) iter->second / ngramtotal;
           const EncNGram ngram = iter->first;
           ngram.decode(classdecoder) ;
           //cout << "NGRAM " << "N=" << ngram.n() << " SIZE=" << (int) ngram.size() << " DECODED=" << ngram.decode(classdecoder) << endl;
           cout << setprecision(numeric_limits<double>::digits10 + 1) << ngram.decode(classdecoder) << '\t' << iter->second << '\t' << freq << endl;
       }
       
       /*for(unordered_map<string,int>::iterator iter = skipgrams.begin(); iter != skipgrams.end(); iter++ ) {
           const double freq = (double) iter->second / skipgramtotal;
           cout << setprecision(numeric_limits<double>::digits10 + 1) << iter->first << '\t' << iter->second << '\t' << freq << endl;
       } */      
 
    }   

}

