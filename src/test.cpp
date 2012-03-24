#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>
#include <patternmodel.h>
#include <common.h>


using namespace std;

int main( int argc, char *argv[] ) {
	//string model = argv[1];
	string classfile = argv[1];
	string querystring = argv[2];

	ClassDecoder classdecoder = ClassDecoder(classfile);
    ClassEncoder encoder = ClassEncoder(classfile);
	
	unsigned char buffer[65536];
	char buffersize = encoder.encodestring(querystring, buffer);
    EncNGram ngram = EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte
    
    vector<EncNGram*> subngrams;
    ngram.subngrams(subngrams);
    for (vector<EncNGram*>::iterator iter2 = subngrams.begin(); iter2 != subngrams.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }
    
    
	
}
