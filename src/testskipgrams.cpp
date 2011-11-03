#include <algorithms.h>
#include <ngram.h>
#include <vector>
#include <iostream>
#include <utility>

using namespace std;

int main( int argc, char *argv[] ) {
    
    const unsigned char data[4] = {2,4,0,5};         
    EncNGram bigram = EncNGram(data, (char) 4);
    const unsigned char data2[5] = {6,0,7,0,8};  
    EncNGram trigram = EncNGram(data2, (char) 5);
    
    vector<EncNGram*> subngrams;
    subngrams.push_back(&bigram);
    subngrams.push_back(&trigram);
    
    
    
    vector<int> skipref;
    skipref.push_back(1);
    //skipref.push_back(3)
    cerr << "SKIPS: 1 3 2" << endl;
    cerr << "SUBNGRAM #1 N=" << (int) bigram.n() << endl;
    bigram.out();
    cout << endl;
    
    cerr << "SUBNGRAM #2 N=" << (int) trigram.n() << endl;
    trigram.out(); 
    cout << endl;
    

         
    EncSkipGram singleskipgram = EncSkipGram(subngrams,skipref,false,false);
    cerr << "SINGLESKIPGRAM N=" << (int) singleskipgram.n() << endl;
    singleskipgram.out();
    cout << endl;
    
    
    skipref.push_back(3);
    
    EncSkipGram initial2skipgram = EncSkipGram(subngrams,skipref,true,false);
    cerr << "INITIAL2SKIPGRAM N=" <<  (int)  initial2skipgram.n() << endl;
    initial2skipgram.out();
    cout << endl;
    
    
    EncSkipGram final2skipgram = EncSkipGram(subngrams,skipref,false,true);
    cerr << "FINAL2SKIPGRAM N=" <<  (int)  final2skipgram.n() << endl;
    final2skipgram.out();
    cout << endl;
    
    skipref.push_back(2);
    EncSkipGram threeskipgram = EncSkipGram(subngrams,skipref,true,true);
    cerr << "SKIPGRAM N=" <<  (int) threeskipgram.n() << endl;
    threeskipgram.out(); 
    cout << endl;
    
    
    const unsigned char data3[15] = {2,0,3,0,4,0,5,0,6,0,7,0,8,0,9};    
    EncNGram ngram = EncNGram(data3, (char) 15);
    
    cerr << "NGRAM N=" << (int) ngram.n() << endl;
    ngram.out();
    cout << endl;
        
    
    vector< vector< pair<int,int> > > gaps;
    compute_multi_skips(gaps, vector<pair<int,int> >(), ngram.n());
    
    for (int i = 0; i < gaps.size(); i++)  {
        for (int j = 0; j < gaps[i].size(); j++)  {
            cerr << "b"<<gaps[i][j].first<<"l"<<gaps[i][j].second << " ";
        }   
        cerr << endl;
    }
        
    //void compute_multi_skips(, vector<pair<int,int> > path, const int n, const int maxskips, const int skipnum, const int leftmargin) {    
}
