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
	//string classfile = argv[1];
	
	string classfile = "/tmp/colibritest";
	
    
    ofstream f;
    f.open(classfile.c_str(), ios::out);
    f << "2\tbe\n3\tTo\n4\tto\n5\tor\n6\tnot\n";            
    f.close();

    
    
	
	
	
	ClassDecoder classdecoder = ClassDecoder(classfile);
    ClassEncoder encoder = ClassEncoder(classfile);
	
 	
 	cerr << "Encoding n-gram from string input" << endl;
 	string querystring = "To be or not to be";
	EncNGram ngram = encoder.input2ngram(querystring); 	

	cerr << "Ngram: " << ngram.decode(classdecoder) << endl;
	cerr << "N: " << (int) ngram.n() << endl;
	cerr << "Size: " << (int) ngram.size() << endl;
 	cerr << "Subgrams: " << endl;
 	
    vector<EncNGram*> subngrams;
    ngram.subngrams(subngrams);
    for (vector<EncNGram*>::iterator iter2 = subngrams.begin(); iter2 != subngrams.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }
    cerr << "----------------------------------------------------" << endl;
    cerr << "Encoding skip-gram from string input" << endl;
	string querystring2 = "To {*1*} or {*1*} to be";
	
	EncSkipGram skipgram = encoder.input2skipgram(querystring2);
	
	cerr << "Skipgram: " << skipgram.decode(classdecoder) << endl;
	cerr << "N: " << (int) skipgram.n() << endl;
	cerr << "Size: " << (int) skipgram.size() << endl;
 
    cerr << "Parts: " << endl;
    vector<EncNGram*> parts;
    skipgram.parts(parts);
    for (vector<EncNGram*>::iterator iter2 = parts.begin(); iter2 != parts.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }
	cerr << "----------------------------------------------------" << endl;
	vector<EncNGram*> parts2;
/*
    cerr << "Encoding skip-gram from string input" << endl;
	string querystring3 = "{*1*} be {*1*} not {*2*}";
	
	EncSkipGram skipgraminv = encoder.input2skipgram(querystring3);
	
	cerr << "Skipgram: " << skipgraminv.decode(classdecoder) << endl;
	cerr << "N: " << (int) skipgraminv.n() << endl;
	cerr << "Size: " << (int) skipgraminv.size() << endl;	
	
	cerr << "Parts: " << endl;
    
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }
  */  
    cerr << "----------------------------------------------------" << endl;

    cerr << "Encoding skip-gram from string input" << endl;
	string querystring4 = "be {*1*} not";
	
	EncSkipGram skipgraminv2 = encoder.input2skipgram(querystring4);
	
	cerr << "Skipgram: " << skipgraminv2.decode(classdecoder) << endl;
	cerr << "N: " << (int) skipgraminv2.n() << endl;
	cerr << "Size: " << (int) skipgraminv2.size() << endl;	
	
	cerr << "Parts: " << endl;
    parts2.clear();
    skipgraminv2.parts(parts2);
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }    
    cerr << "----------------------------------------------------" << endl;
    cerr << "Re-instantiating skipgram with skipcontent" << endl;
    EncNGram rengram = skipgram.instantiate(&skipgraminv2, parts, parts2);
    
    cerr << "Ngram: " << rengram.decode(classdecoder) << endl;
	cerr << "N: " << (int) rengram.n() << endl;
	cerr << "Size: " << (int) rengram.size() << endl;
	 
 		  	
}
