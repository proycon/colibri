#include "classencoder.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include "common.h"

using namespace std;
using namespace folia;


bool validclass(unsigned int cls) {
	/* checks if the integer value is a valid class that has no null bytes in its representation */
    do {
    	int r = cls % 256;
    	if (r == 0) return false;
  		cls = cls / 256;
    } while (cls >= 256);
    return true;
}

unsigned char * inttobytes(unsigned int cls, int & length) {	
	//compute length of byte array
	unsigned int cls2 = cls;
	length = 0;
	do {
		cls2 = cls2 / 256;
		length++;
	} while (cls2 > 0);
	unsigned char * byterep = new unsigned char[length];
	int i = 0;
    do {
    	int r = cls % 256;
    	byterep[i++] = (unsigned char) r;
    	cls = cls / 256;
    } while (cls > 0);
	return byterep;    
}

ClassEncoder::ClassEncoder() {
    unknownclass = 2;
    bosclass = 3;
    eosclass = 4;
    highestclass = 5; //5 and lower are reserved
}

ClassEncoder::ClassEncoder(const string & filename) {
       unknownclass = 2;
       highestclass = 0; 
       bosclass = 3;
       eosclass = 4;
       
	   ifstream IN;
	   IN.open( filename.c_str() );    
       if (!(IN)) {
           cerr << "File does not exist: " << filename << endl;
           exit(3);
       }
        while (IN.good()) {
          string line;
          getline(IN, line);              
          const int s = line.size();
          for (int i = 0; i < s; i++) {
              if (line[i] == '\t') {
                  const string cls_s = string(line.begin(), line.begin() + i);
                  unsigned int cls = (unsigned int) atoi(cls_s.c_str());
                  const string word = string(line.begin() + i + 1, line.end());                  
                  classes[word] = cls;
                  if (cls == 2) {
                    unknownclass = 0;
                  }
                  if (cls > (unsigned int) highestclass) highestclass = cls;
                  //cerr << "CLASS=" << cls << " WORD=" << word << endl;
              }
              
          }
        }        
        IN.close();  
        
        if (unknownclass == 0) {
            highestclass++;
            unknownclass = highestclass;
            classes["{UNKNOWN}"] = unknownclass;
        } else {        
            classes["{UNKNOWN}"] = unknownclass;
            classes["{BEGIN}"] = bosclass;
            classes["{END}"] = eosclass;
        }
}


void ClassEncoder::processcorpus(const string & filename, unordered_map<string,int> & freqlist) {
	   //compute frequency list of all words        
       ifstream IN;
       IN.open( filename.c_str() );    
       if (!(IN)) {
           cerr << "File does not exist: " << filename << endl;
           exit(3);
       }       
        while (IN.good()) {
          string line;
          getline(IN, line);         
          int begin = 0;
          const int s = line.size();
          for (int i = 0; i < s; i++) {
              if ((line[i] == ' ') || (i == s - 1)) {
              	  int offset = 0;
              	  if (i == s - 1) offset = 1;              	  
              	  string word = string(line.begin() + begin, line.begin() + i + offset);              	  
              	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ")) {
              	    word = trim(word, " \t\n\r"); //trim whitespace, control characters
              	  	freqlist[word]++;
              	  }
              	  begin = i+ 1; 
              }
              
          }
        }        
        IN.close();       
}


void ClassEncoder::processfoliacorpus(const string & filename, unordered_map<string,int> & freqlist) {
    Document doc;
    doc.readFromFile(filename);
    
    vector<Word*> words = doc.words();
    for (vector<Word*>::iterator iterw = words.begin(); iterw != words.end(); iterw++) {
        Word * word = *iterw;
        const string wordtext = word->str();
        freqlist[wordtext]++;
    }
    
}

void ClassEncoder::buildclasses(unordered_map<string,int> freqlist) {

        //sort by occurrence count  using intermediate representation
        multimap<const int, const string> revfreqlist;
        for (unordered_map<string,int>::iterator iter = freqlist.begin(); iter != freqlist.end(); iter++) {
        	revfreqlist.insert( pair<const int,const string>(-1 * iter->second, iter->first) );
        } 
        
        freqlist.clear();
        
        int cls = highestclass;
        for (multimap<const int,const string>::iterator iter = revfreqlist.begin(); iter != revfreqlist.end(); iter++) {
            if (!classes.count(iter->second)) { //check if it doesn't already exist, in case we are expanding on existing classes 
        	    cls++;
        	    while (!validclass(cls)) cls++;        	
        	    classes[iter->second] = cls;
            }
        }
}

void ClassEncoder::build(const string & filename) {
	    unordered_map<string,int> freqlist;
	    if (filename.rfind(".xml") != string::npos) {
	        processfoliacorpus(filename, freqlist);
	    } else {
	        processcorpus(filename, freqlist);
	    }
        buildclasses(freqlist);
}


void ClassEncoder::build(vector<string> & files) {
	    unordered_map<string,int> freqlist;
	    	    
	    for (vector<string>::iterator iter = files.begin(); iter != files.end(); iter++) {
	        const string filename = *iter;
	        cerr << "Processing " << filename << endl;
	        if (filename.rfind(".xml") != string::npos) {
	            processfoliacorpus(filename, freqlist);
	        } else {
	            processcorpus(filename, freqlist);
	        }
	        
	    } 	    	    
        buildclasses(freqlist);
}

void ClassEncoder::save(const string & filename) {
	ofstream OUT;
	OUT.open( filename.c_str() );
	for (std::unordered_map<std::string,unsigned int>::iterator iter = classes.begin(); iter != classes.end(); iter++) {
	    if (iter->second != unknownclass) OUT << iter->second << '\t' << iter->first << endl;
	}
	OUT.close();
}

        
vector<unsigned int> ClassEncoder::encodeseq(const vector<string> & seq) {
    vector<unsigned int> result;
    const int l = seq.size();
    for (int i = 0; i < l; i++) 
        result.push_back( classes[seq[i]] ); 
    return result;
}

int ClassEncoder::encodestring(const string & line, unsigned char * outputbuffer, bool allowunknown, bool autoaddunknown) {
	  int outputcursor = 0;
      int begin = 0;      
      const int l = line.length();
      for (int i = 0; i < l; i++) {
      	  if ((line[i] == ' ') || (i == l - 1)) {
          	  string word;
          	  if (line[i] == ' ') {
          	  	word  = string(line.begin() + begin, line.begin() + i);
          	  } else {
			   	word  = string(line.begin() + begin, line.begin() + i + 1);
          	  }
          	  begin = i+1;
          	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ")) {
          	    unsigned int cls;
          	    if (classes.count(word) == 0) {
                    if (autoaddunknown) {
                        cls = ++highestclass;
                        classes[word] = cls;  
                        added[cls] = word;
          	    	} else if (!allowunknown) {	
  	        			cerr << "ERROR: Unknown word '" << word << "', does not occur in model" << endl;
  	        			return 0;         
	  	        	} else {
	  	        		cerr << "WARNING: Unknown word '" << word << "', does not occur in model. Replacing with placeholder" << endl;
	  	        		cls = unknownclass;	
	  	        	}    	
          	    } else {
          	  		cls = classes[word];
          	  	}
          	  	int length = 0;
  	        	const unsigned char * byterep = inttobytes(cls, length);
  	        	if (length == 0) {
  	        		cerr << "INTERNAL ERROR: Error whilst encoding '" << word << "' (class " << cls << "), length==0, not possible!" << endl;
  	        		exit(13);
  	        	}  	        		
  	        	//cerr << "writing " << word << " as " << cls << " in " << length << " bytes" << endl;
  	        	for (int j = 0; j < length; j++) {
  	        		outputbuffer[outputcursor++] = byterep[j];
  	        	}  	        	
  	        	//OUT.write((const char *) byterep, length);
  	        	delete [] byterep;
  	        	outputbuffer[outputcursor++] = 0; //write separator
  	        	//OUT.write(&zero, sizeof(char)); //write separator 
          	  }			 
          }
      }
      return outputcursor; //doing -1 to break of possible trailing zero bytes breaks stuff
}

int ClassEncoder::encodestring(const string & line, unsigned char * outputbuffer, unsigned char * skipconf, char * skipcount,  bool allowunknown,  bool autoaddunknown) {
	  int outputcursor = 0;
	  *skipcount = 0;
      int begin = 0;      
      bool finalskip = false;
      const int l = line.length();
      for (int i = 0; i < l; i++) {
      	  if ((line[i] == ' ') || (i == l - 1)) {
          	  string word;
          	  if (line[i] == ' ') {
          	  	word  = string(line.begin() + begin, line.begin() + i);
          	  } else {
			   	word  = string(line.begin() + begin, line.begin() + i + 1);
          	  }
          	  begin = i+1;
          	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ")) {
          	    unsigned int cls = 0;
          	    if (word == "{*1*}") { //not very elegant, but gets the job done for now
          	    	skipconf[(int) (*skipcount)++] = 1; 
          	    } else if (word == "{*2*}") {
          	    	skipconf[(int) (*skipcount)++] = 2;
          	    } else if (word == "{*3*}") {
          	    	skipconf[(int) (*skipcount)++] = 3;
          	    } else if (word == "{*4*}") {
          	    	skipconf[(int) (*skipcount)++] = 4;
          	    } else if (word == "{*5*}") {
          	    	skipconf[(int) (*skipcount)++] = 5;
          	    } else if (word == "{*6*}") {
          	    	skipconf[(int) (*skipcount)++] = 6;
          	    } else if (word == "{*7*}") {
          	    	skipconf[(int) (*skipcount)++] = 7;
          	    } else if (word == "{*8*}") {
          	    	skipconf[(int) (*skipcount)++] = 8;
          	    } else if (word == "{*9*}") {
          	    	skipconf[(int) (*skipcount)++] = 9;
          	    } else if (classes.count(word) == 0) {
          	        if (autoaddunknown) {
                        cls = ++highestclass;
                        classes[word] = cls;
                        added[cls] = word;
          	    	} else if (!allowunknown) { 	
  	        			cerr << "ERROR: Unknown word '" << word << "', does not occur in model" << endl;
  	        			return 0;         
	  	        	} else {
	  	        		cerr << "WARNING: Unknown word '" << word << "', does not occur in model. Replacing with placeholder" << endl;
	  	        		cls = unknownclass;	
	  	        	}    	
          	    } else {
          	  		cls = classes[word];
          	  	}
          	  	
	      	  	
		      	if (cls > 0) {
		      	    //word, no skip
		      	    finalskip = false;
		      		int length = 0;
	  	        	unsigned char * byterep = inttobytes(cls, length);
   	   	        	if (length == 0) {
	  	        		cerr << "INTERNAL ERROR: Error whilst encoding '" << word << "' (class " << cls << "), length==0, not possible!" << endl;
	  	        		exit(13);
	  	        	}  	        		
	  	        	//cerr << "writing " << word << " as " << cls << " in " << length << " bytes" << endl;
	  	        	for (int j = 0; j < length; j++) {
	  	        		outputbuffer[outputcursor++] = byterep[j];
	  	        	}  	        	
	  	        	//OUT.write((const char *) byterep, length);	  	        
	  	        	delete [] byterep;
				} else {
				    //skip
				    if (outputcursor == 0) {
				        //initial skip
				        outputbuffer[outputcursor++] = 0;
				    }
				    finalskip = true;
				}
  	        	outputbuffer[outputcursor++] = 0; //write separator
  	        	//OUT.write(&zero, sizeof(char)); //write separator 
          	  }			 
          }
      }
      if (finalskip) outputbuffer[outputcursor++] = 0;
      return outputcursor;
}

EncSkipGram ClassEncoder::input2skipgram(const std::string & querystring, bool allowunknown, bool autoaddunknown) {
	unsigned char buffer[65536];
	unsigned char skipref[4];
	char skipcount = 0;
	char buffersize = encodestring(querystring, buffer, skipref, &skipcount, allowunknown, autoaddunknown);
	EncSkipGram skipgram = EncSkipGram(buffer,buffersize-1,skipref,skipcount); //-1 to strip last \0 byte
	return skipgram;
}

EncNGram ClassEncoder::input2ngram(const std::string & querystring, bool allowunknown,  bool autoaddunknown) {
	unsigned char buffer[65536];
	char buffersize = encodestring(querystring, buffer, allowunknown, autoaddunknown);
	EncNGram ngram = EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte
	return ngram;
}

EncAnyGram * ClassEncoder::input2anygram(const std::string & querystring, bool allowunknown,  bool autoaddunknown) {
	unsigned char buffer[65536];
	unsigned char skipref[4];
	char skipcount = 0;
	char buffersize = encodestring(querystring, buffer, skipref, &skipcount, allowunknown, autoaddunknown);
	if (skipcount == 0) {
		return new EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte;
	} else {
		return new EncSkipGram(buffer,buffersize-1,skipref,skipcount); //-1 to strip last \0 byte;
	}
}

void ClassEncoder::add(std::string s, unsigned int cls) {
    classes[s] = cls;
    if (cls > highestclass) highestclass = cls;
}

void ClassEncoder::encodefile(const std::string & inputfilename, const std::string & outputfilename, bool allowunknown, bool autoaddunknown, bool append) {
    const char zero = 0;
    const char one = 1;
	    
    if ((inputfilename.rfind(".xml") != string::npos) ||  (inputfilename.rfind(".bz2") != string::npos) ||  (inputfilename.rfind(".gz") != string::npos)) {
        //FoLiA
        Document doc;
        doc.readFromFile(inputfilename);
        
	    ofstream OUT;
	    if (append) {
	        OUT.open(outputfilename.c_str(), ios::app | ios::binary);
	        if (OUT.tellp() > 0) {
	            OUT.write(&one, sizeof(char)); //newline          
          	    OUT.write(&zero, sizeof(char)); //write separator
	        }
	    } else {
	        OUT.open(outputfilename.c_str(), ios::out | ios::binary);
	    }
	    unsigned char outputbuffer[65536];
	    int outputsize = 0;
	    unsigned int linenum = 1;
	    vector<Word*> words = doc.words();
	    const size_t wl = words.size();
	    FoliaElement * prevparent = NULL;
	    string line = "";
	    for (int i = 0; i < wl; i++) {
	        Word * word = words[i];
	        if ((!line.empty()) && (word->parent() != prevparent) && (i< wl -1)) {
	            outputsize = encodestring(line, outputbuffer, allowunknown, autoaddunknown);     
	            if (outputsize > 0) OUT.write((const char *) outputbuffer, outputsize);
	        	OUT.write(&one, sizeof(char)); //newline          
          	    OUT.write(&zero, sizeof(char)); //write separator
          	    linenum++;        
          	    line = "";
	        }
            prevparent = word->parent();
        	if (line.empty()) {
                line += word->str(); 
            } else {
                line += " " + word->str();
            }
        } 
        if (!line.empty()) {
            outputsize = encodestring(line, outputbuffer, allowunknown, autoaddunknown);     
	        if (outputsize > 0) OUT.write((const char *) outputbuffer, outputsize);
        }
	    cerr << "Encoded " << linenum << " lines" << endl;
	    OUT.close();
	            
    } else {
	    ofstream OUT;
	    ifstream IN;
	    if (append) {
	        OUT.open(outputfilename.c_str(), ios::app | ios::binary);
	        if (OUT.tellp() > 0) {
	            OUT.write(&one, sizeof(char)); //newline          
          	    OUT.write(&zero, sizeof(char)); //write separator
	        }
	    } else {
	        OUT.open(outputfilename.c_str(), ios::out | ios::binary);
	    }	
	    IN.open(inputfilename.c_str());
	    unsigned char outputbuffer[65536];
	    int outputsize = 0;
	    unsigned int linenum = 1;
	    while (IN.good()) {	
          string line = "";
          getline(IN, line);
          if ((outputsize > 0) && (!IN.eof())) {      
          	 OUT.write(&one, sizeof(char)); //newline          
          	 OUT.write(&zero, sizeof(char)); //write separator
          	 linenum++;      	 
          }           
          outputsize = encodestring(line, outputbuffer, allowunknown, autoaddunknown);
          if (outputsize > 0) OUT.write((const char *) outputbuffer, outputsize);                        
        }
        cerr << "Encoded " << linenum << " lines" << endl;
	    IN.close();
	    OUT.close();
	}
}
