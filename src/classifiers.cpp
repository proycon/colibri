#include <classifiers.h>

using namespace std;
using namespace Timbl;

BuildClassifier::BuildClassifier(const std::string & id, bool appendmode = false, bool exemplarweights = false) {
    ID = id;
    trainfile = string(id + ".train");
    opened = false;        
    featurevectorsize = 0;    
    this->appendmode = appendmode;
    this->exemplarweights = exemplarweights;
}        

BuildClassifier::~BuildClassifier() {
    if (opened) outputfile.close();    
}

void BuildClassifier::addinstance(vector<const string> featurevector, const string & label, double exemplarweight) {
    if (!opened) {
        if (appendmode) {
            outputfile.open(trainfile, ios::app);
        } else {
            outputfile.open(trainfile);
        }
    }

    if (featurevectorsize == 0) {
        featurevectorsize = featurevector.size();
    } else if (featurevector.size() != featurevectorsize) {
        cerr << "INTERNAL ERROR: Passed feature vector of size " << featurevector.size() << ", expected: " << featurevectorsize << endl;
        exit(6); 
    } 
    
    for (vector<const string>::const_iterator iter = featurevector.begin(); iter != featurevector.end(); iter++) {    
        const string feature = *iter;
        outputfile << feature << ' ';  
    } 
    outputfile << label;
    if (examplarweights) {
        outputfile << ' ' << exemplarweight;
    }
    outputfile << endl;
}

void BuildClassifier::train(const string & timbloptions) {
    ibasefile = string(id + ".ibase");
    TimblAPI * timblexp = new TimblAPI( timbloptions , ID );
    timblexp->Learn(trainfile);   
    timblexp->WriteInstanceBase( ibasefile );
    delete timblexp;    
}

NClassifierArray::NClassifierArray(const string & id, int maxn, int leftcontextsize, int rightcontextsize) {
    this->leftcontextsize = leftcontextsize;
    this->rightcontextsize = rightcontextsize;
    this->maxn = maxn;
}

NClassifierArray::build(const string & enctraincorpusfile, const TranslationTable * ttable, ClassDecoder * sourceclassdecoder, ClassDecoder * targetclassdecoder, bool exemplarweights) {
    vector<Classifier*> buildarray;
    ifstream IN;
    IN.open(enctraincorpusfile);
    if (!IN.good()) {
        cerr << "Training corpus not found: " << enctraincorpusfile << endl;
        exit(2);
    }
    buildarray.push_back(NULL);    
    for (int n = 1; n <= 9; n++) {
        buildarray.push_back( new BuildClassifier(String(id + ".n." + n), false, exemplarweights) ); 
    }
    const int BUFFERSIZE = 65536;
    unsigned char line[BUFFERSIZE];
	uint32_t linenum = 0;
    while (IN.good()) {
        const int linesize = readline(IN, line, BUFFERSIZE);
		vector<pair<const EncAnyGram*, CorpusReference> > patterns = getpatterns(buffer,buffersize, true, linenum,1,maxn);
		for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = patterns.begin(); iter != patterns.end(); iter++) {
			const EncAnyGram * anygram = iter->first;
			const CorpusReference ref = iter->second;
			
			//check if pattern occurs in phrasetable
			if (ttable->getsourcekey(anygram) != NULL) {
			     //yes
			     
			     vector<const string> featurevector;
			     //obtain left context
			     
			     
			     //obtain focus part			    
                 const int n = anygram->n();			     			        
			     for (int i = 0; i < n; i++) {
			        const EncNGram * unigram = anygram->gettoken(i);
			        featurevector.push_back(unigram->decode(sourceclassdecoder) )
			        delete unigram;
			     }			      
			     
			     //obtain right context
			     
			     //compose feature vector
			     const string label = 
			}
			
		}         
    }
    
    
    
    
    do {
    	linenum++; 
    	getline(&fin,line);    	
    	if (!line.empty()) {
			int buffersize = encoder.encodestring(line, buffer, allowunknown); //last bool is
			if (exact) {
				//TODO
			} else {    	
				vector<pair<const EncAnyGram*, CorpusReference> > patterns = getpatterns(buffer,buffersize, true, linenum,minn,maxn);
				for (vector<pair<const EncAnyGram*, CorpusReference> >::iterator iter = patterns.begin(); iter != patterns.end(); iter++) {
					const EncAnyGram * anygram = iter->first;
					const CorpusReference ref = iter->second;
					outputinstance(anygram, ref, decoder);
				} 
			}
		}
    } while (!cin.eof() && (repeat));	  
    
    for (int n = 1; n <= 9; n++) {
        delete buildarray[n];
    }    
}



