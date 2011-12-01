#include "classdecoder.h"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

unsigned int bytestoint(const unsigned char* a, const int l) {
    int result = 0;
    for (int i = 0; i < l; i++) {
        result += *(a + i) * pow(256,i);
    }
    return result;
}


ClassDecoder::ClassDecoder(const string filename) {
        
       ifstream *IN =  new ifstream( filename.c_str() );    
       if (!(*IN)) {
           cerr << "File does not exist: " << filename << endl;
           exit(3);
       }
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
    for (int i = 0; i < l; i++) 
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
            const unsigned int cls = bytestoint(buffer, n);  
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



int readline(istream* IN, unsigned char* buffer, const int MAXBUFFERSIZE) {
    int n = 0;
    short eolsequence = 0; //3 stages: 0 1 0 , when all three are found, we have a sentence
    while (IN->good()) {
        char bufchar;
        IN->get(bufchar);
        unsigned char c = (unsigned char) bufchar;
        if (n >= MAXBUFFERSIZE) {
            cerr << "Buffer overflow in classdecoder readline(): " << n << ". This indicates a sentence in the data that is far longer than anything reasonable, probably an error in sentence segmentation or tokenisation." << endl;
            exit(666);
        }
        buffer[n] = c;        
        if (c == 0) {
            eolsequence++;
            if (eolsequence == 3) return n - 2; //minus last two bytes of the 0 1 0 bytes (retaining final \0)
        } else if (c == 1) {
            if (eolsequence == 1) {
                eolsequence++;
            } else {
                eolsequence = 0;
            }
        } else {
            eolsequence = 0;
        }
        n++;
    }    
    return 0;
}

const int countwords(const unsigned char* data, const int l) {    
    int words = 1;
    for (int i = 0; i < l; i++) {
        if ( (data[i] == 0) && (i > 0) && (i < l - 1) ) words++;
    }
    return words;
}
