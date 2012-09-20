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
    f << "2\tbe\n3\tTo\n4\tto\n5\tor\n6\tnot\n73477272\tblah\n";            
    f.close();

    
    
	
	
	
	ClassDecoder classdecoder = ClassDecoder(classfile);
    ClassEncoder encoder = ClassEncoder(classfile);
	
 	
 	cerr << "Encoding n-gram from string input" << endl;
 	string querystring = "To be or not to be";
	EncNGram ngram = encoder.input2ngram(querystring, true); 	

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
	
	EncSkipGram skipgram = encoder.input2skipgram(querystring2, true);
	
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

    cerr << "Extracting skip content based on skip gram and full instance" << endl;
    EncSkipGram skipcontent = skipgram.extractskipcontent(ngram);
    
    cerr << "Skipcontent: " << skipcontent.decode(classdecoder) << endl;
	cerr << "N: " << (int) skipcontent.n() << endl;
	cerr << "Size: " << (int) skipcontent.size() << endl;


	cerr << "----------------------------------------------------" << endl;

    cerr << "Getwords" << endl;
    pair<int,int> wordspos = getwords(ngram.data, ngram.size(), 2, 2);
    EncNGram ngram2 = EncNGram(ngram.data + wordspos.first, wordspos.second);
    
    cerr << "Ngram: " << ngram2.decode(classdecoder) << endl;
	cerr << "N: " << (int) ngram2.n() << endl;
	cerr << "Size: " << (int) ngram2.size() << endl;
	
	
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
	
	EncSkipGram skipgraminv2 = encoder.input2skipgram(querystring4, true);
	
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
    EncNGram rengram = skipgram.instantiate(&skipgraminv2, parts2);
    
    cerr << "Ngram: " << rengram.decode(classdecoder) << endl;
	cerr << "N: " << (int) rengram.n() << endl;
	cerr << "Size: " << (int) rengram.size() << endl;
	 
	cerr << "----------------------------------------------------" << endl;
	string querystring5 = "be {*1*} not {*2*} be";
	EncSkipGram skipgram5 = encoder.input2skipgram(querystring5, true);
	cout << skipgram5.decode(classdecoder) << endl;
	
	cerr << "Parts: " << endl;
    parts2.clear();
    skipgram5.parts(parts2);
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }    	 
	cerr << "getgaps: " << endl;
	vector<pair<int,int> > gaps;
	skipgram5.getgaps(gaps);
	for (vector<pair<int,int >>::iterator iter2 = gaps.begin(); iter2 != gaps.end(); iter2++) {
	    cout << iter2->first << ':' << iter2->second << endl; 
	}
	cerr << "getparts: " << endl;
	vector<pair<int,int> > p;
	skipgram5.getparts(p);
	for (vector<pair<int,int >>::iterator iter2 = p.begin(); iter2 != p.end(); iter2++) {
        cout << iter2->first << ':' << iter2->second << endl; 
	}	
	cerr << "mask: " << endl;
	vector<bool> m;
	skipgram5.mask(m);
	for (vector<bool>::iterator iter2 = m.begin(); iter2 != m.end(); iter2++) {
        if (*iter2)  {
            cout << "1";
        } else {
            cout << "0";
        }
	}		
	cout << endl;
	
    cerr << "gettoken(5): " << endl;
    EncNGram * token = skipgram5.gettoken(5);
    cout << token->decode(classdecoder) << endl;
	
	cerr << "----------------------------------------------------" << endl;
	string querystring6 = "be {*1*} not";
	EncSkipGram skipgram6 = encoder.input2skipgram(querystring6, true);
	cout << skipgram6.decode(classdecoder) << endl;
	skipgram6.out();
	
	cerr << "Parts: " << endl;
    parts2.clear();
    skipgram6.parts(parts2);
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }    	 
	cerr << "getgaps: " << endl;
	vector<pair<int,int> > gaps6;
	skipgram6.getgaps(gaps6);
	for (vector<pair<int,int >>::iterator iter2 = gaps6.begin(); iter2 != gaps6.end(); iter2++) {
	    cout << iter2->first << ':' << iter2->second << endl; 
	}
	cerr << "getparts: " << endl;
	vector<pair<int,int> > p6;
	skipgram6.getparts(p6);
	for (vector<pair<int,int >>::iterator iter2 = p6.begin(); iter2 != p6.end(); iter2++) {
        cout << iter2->first << ':' << iter2->second << endl; 
	}	



    cerr << "----------------------------------------------------" << endl;
    string querystring7 = "blah {*1*} or {*2*} blah";
	EncSkipGram skipgram7 = encoder.input2skipgram(querystring7, true);
	cout << skipgram7.decode(classdecoder) << endl;
	
	cerr << "Parts: " << endl;
    parts2.clear();
    skipgram7.parts(parts2);
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }    	 
	cerr << "getgaps: " << endl;
	gaps.clear();
	skipgram7.getgaps(gaps);
	for (vector<pair<int,int >>::iterator iter2 = gaps.begin(); iter2 != gaps.end(); iter2++) {
	    cout << iter2->first << ':' << iter2->second << endl; 
	}
	cerr << "getparts: " << endl;
	p.clear();
	skipgram7.getparts(p);
	for (vector<pair<int,int >>::iterator iter2 = p.begin(); iter2 != p.end(); iter2++) {
        cout << iter2->first << ':' << iter2->second << endl; 
	}	
    cerr << "mask: " << endl;
	m.clear();
	skipgram7.mask(m);
	for (vector<bool>::iterator iter2 = m.begin(); iter2 != m.end(); iter2++) {
        if (*iter2)  {
            cout << "1";
        } else {
            cout << "0";
        }
	}		
	cout << endl;
    
    cerr << "gettoken(5): " << endl;
    token = skipgram7.gettoken(5);
    cout << token->decode(classdecoder) << endl;
    
    cerr << "----------------------------------------------------" << endl;
    string querystring8 = "{*1*} or blah {*2*}";
	EncSkipGram skipgram8 = encoder.input2skipgram(querystring8, true);
	cout << skipgram8.decode(classdecoder) << endl;
	
		
	cerr << "Parts: " << endl;
    parts2.clear();
    skipgram8.parts(parts2);
    for (vector<EncNGram*>::iterator iter2 = parts2.begin(); iter2 != parts2.end(); iter2++) {                
        const EncAnyGram * subngram = *iter2;
    	cout << "'" << subngram->decode(classdecoder) << "'" << endl;
    }    	 
	cerr << "getgaps: " << endl;
	gaps.clear();
	skipgram8.getgaps(gaps);
	for (vector<pair<int,int >>::iterator iter2 = gaps.begin(); iter2 != gaps.end(); iter2++) {
	    cout << iter2->first << ':' << iter2->second << endl; 
	}
	cerr << "getparts: " << endl;
	p.clear();
	skipgram8.getparts(p);
	for (vector<pair<int,int >>::iterator iter2 = p.begin(); iter2 != p.end(); iter2++) {
        cout << iter2->first << ':' << iter2->second << endl; 
	}	
    cerr << "mask: " << endl;
	m.clear();
	skipgram8.mask(m);
	for (vector<bool>::iterator iter2 = m.begin(); iter2 != m.end(); iter2++) {
        if (*iter2)  {
            cout << "1";
        } else {
            cout << "0";
        }
	}		
	cout << endl;
    
    cerr << "gettoken(2): " << endl;
    token = skipgram8.gettoken(2);
    cout << token->decode(classdecoder) << endl;
	    
     		  	
}
