#include "classdecoder.h"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

unsigned int bytestoint(unsigned char* a, const int l) {
    int result = 0;
    for (int i = 0; i < l+1; i++) {
        result += *(a + i) * pow(256,i);
    }
    return result;
}


ClassDecoder::ClassDecoder(const string filename) {
       ifstream *IN =  new ifstream( filename.c_str() );    

        while (IN->good()) {
          string line;
          getline(*IN, line);              
          for (int i = 0; i < line.size(); i++) {
              if (line[i] == '\t') {
                  const string cls_s = string(line.begin(), line.begin() + i);
                  unsigned int cls = (unsigned int) atoi(cls_s.c_str());
                  const string word = string(line.begin() + i + 1, line.end());
                  classes[cls] = word;
                  //cerr << "CLASS=" << cls << " WORD=" << word << endl;
              }
              
          }
        }
        IN->close();  
}

        
vector<string> ClassDecoder::decodeseq(vector<int> seq) {
    vector<string> result;
    const int l = seq.size();
    for (int i; i < l; i++) 
        result.push_back( classes[seq[i]] ); 
    return result;
}



void ClassDecoder::decodefile(const string filename) {
    ifstream *IN = new ifstream(filename.c_str()); //, ios::in | ios::binary);
    unsigned char buffer[10];
    int n = 0;
    while (IN->good()) {
        char bufchar;
        IN->get(bufchar);
        
        unsigned char c = (unsigned char) bufchar;
        buffer[n] = c;
        //cout << "READ: " << ((int) c) << endl;
        if (c == 0) {
            //cout << "N: " << n << endl;
            const unsigned int cls = bytestoint(buffer, n - 1);  
            if (cls == 1) {
                cout << endl;
            } else {
                //cout << cls << ' ';
                cout << classes[cls] << ' ';
            }
            n = 0;
        } else {
            n++;
        }
    }
    IN->close();                    
} 



int readline(istream* IN, unsigned char* buffer) {
    int n = 0;
    short eolsequence = 0; //3 stages: 0 1 0 , when all three are found, we have a sentence
    while (IN->good()) {
        char bufchar;
        IN->get(bufchar);
        unsigned char c = (unsigned char) bufchar;
        buffer[n] = c;        
        if (c == 0) {
            eolsequence++;
            if (eolsequence == 3) return n - 3; //minus last three 0 1 0 bytes            
        } else if (c == 1) {
            if (eolsequence == 1) {
                eolsequence++;
            } else {
                eolsequence = 0;
            }
        } else {
            eolsequence = 0;
        }
    }    
    return 0;
}
