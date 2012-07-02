#include "classencoder.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include "common.h"

using namespace std;



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
    highestclass = 5; //5 and lower are reserved
}

ClassEncoder::ClassEncoder(const string & filename) {
       unknownclass = 2;
       highestclass = 0; 
    
	   ifstream IN;
	   IN.open( filename.c_str() );    
       if (!(IN)) {
           cerr << "File does not exist: " << filename << endl;
           exit(3);
       }
        while (IN.good()) {
          string line;
          getline(IN, line);              
          for (int i = 0; i < line.size(); i++) {
              if (line[i] == '\t') {
                  const string cls_s = string(line.begin(), line.begin() + i);
                  unsigned int cls = (unsigned int) atoi(cls_s.c_str());
                  const string word = string(line.begin() + i + 1, line.end());                  
                  classes[word] = cls;
                  if (cls == 2) {
                    unknownclass = 0;
                  }
                  if (cls > highestclass) highestclass = cls;
                  //cerr << "CLASS=" << cls << " WORD=" << word << endl;
              }
              
          }
        }        
        IN.close();  
        
        if (unknownclass == 0) {
            highestclass++;
            unknownclass = highestclass;
            classes["{UNKNOWN}"] = unknownclass;
        }
}

void ClassEncoder::build(const string & filename) {
	   unordered_map<string,int> freqlist;
	    
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
          for (int i = 0; i < line.size(); i++) {
              if ((line[i] == ' ') || (i == line.size() - 1)) {
              	  int offset = 0;
              	  if (i == line.size() - 1) offset = 1;              	  
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
        
        //sort by occurrence count  using intermediate representation
        multimap<const int, const string> revfreqlist;
        for (unordered_map<const string,int>::iterator iter = freqlist.begin(); iter != freqlist.end(); iter++) {
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

int ClassEncoder::encodestring(const string & line, unsigned char * outputbuffer, bool allowunknown) {
	  int outputcursor = 0;
      int begin = 0;      
      for (int i = 0; i < line.length(); i++) {
      	  if ((line[i] == ' ') || (i == line.length() - 1)) {
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
          	    	if (!allowunknown) {	
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
  	        	delete byterep;
  	        	outputbuffer[outputcursor++] = 0; //write separator
  	        	//OUT.write(&zero, sizeof(char)); //write separator 
          	  }			 
          }
      }
      return outputcursor; //doing -1 to break of possible trailing zero bytes breaks stuff
}

int ClassEncoder::encodestring(const string & line, unsigned char * outputbuffer, unsigned char * skipconf, char * skipcount,  bool allowunknown) {
	  int outputcursor = 0;
	  *skipcount = 0;
      int begin = 0;      
      for (int i = 0; i < line.length(); i++) {
      	  if ((line[i] == ' ') || (i == line.length() - 1)) {
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
          	    	skipconf[(*skipcount)++] = 1; 
          	    } else if (word == "{*2*}") {
          	    	skipconf[(*skipcount)++] = 2;
          	    } else if (word == "{*3*}") {
          	    	skipconf[(*skipcount)++] = 3;
          	    } else if (word == "{*4*}") {
          	    	skipconf[(*skipcount)++] = 4;
          	    } else if (word == "{*5*}") {
          	    	skipconf[(*skipcount)++] = 5;
          	    } else if (word == "{*6*}") {
          	    	skipconf[(*skipcount)++] = 6;
          	    } else if (word == "{*7*}") {
          	    	skipconf[(*skipcount)++] = 7;
          	    } else if (word == "{*8*}") {
          	    	skipconf[(*skipcount)++] = 8;
          	    } else if (word == "{*9*}") {
          	    	skipconf[(*skipcount)++] = 9;
          	    } else if (classes.count(word) == 0) {
          	    	if (!allowunknown) { 	
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
	  	        	delete byterep;
				}
  	        	outputbuffer[outputcursor++] = 0; //write separator
  	        	//OUT.write(&zero, sizeof(char)); //write separator 
          	  }			 
          }
      }
      return outputcursor;
}

EncSkipGram ClassEncoder::input2skipgram(const std::string & querystring, bool allowunknown) {
	unsigned char buffer[65536];
	unsigned char skipref[4];
	char skipcount = 0;
	char buffersize = encodestring(querystring, buffer, skipref, &skipcount, allowunknown);
	EncSkipGram skipgram = EncSkipGram(buffer,buffersize-1,skipref,skipcount); //-1 to strip last \0 byte
	return skipgram;
}

EncNGram ClassEncoder::input2ngram(const std::string & querystring, bool allowunknown) {
	unsigned char buffer[65536];
	char buffersize = encodestring(querystring, buffer, allowunknown);
	EncNGram ngram = EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte
	return ngram;
}

EncAnyGram * ClassEncoder::input2anygram(const std::string & querystring, bool allowunknown) {
	unsigned char buffer[65536];
	unsigned char skipref[4];
	char skipcount = 0;
	char buffersize = encodestring(querystring, buffer, skipref, &skipcount, allowunknown);
	if (skipcount == 0) {
		return new EncNGram(buffer,buffersize-1); //-1 to strip last \0 byte;
	} else {
		return new EncSkipGram(buffer,buffersize-1,skipref,skipcount); //-1 to strip last \0 byte;
	}
}

void ClassEncoder::encodefile(const std::string & inputfilename, const std::string & outputfilename, bool allowunknown) {
	const char zero = 0;
	const char one = 1;
	ofstream OUT;
	ifstream IN;
	OUT.open(outputfilename.c_str(), ios::out | ios::binary);	
	IN.open(inputfilename.c_str());
	unsigned char outputbuffer[65536];
	int outputsize = 0;
	unsigned int linenum = 1;
	bool empty = true;
	while (IN.good()) {	
      string line = "";
      getline(IN, line);
      if ((outputsize > 0) && (!IN.eof())) {      
      	 OUT.write(&one, sizeof(char)); //newline          
      	 OUT.write(&zero, sizeof(char)); //write separator
      	 linenum++;      	 
      }           
      outputsize = encodestring(line, outputbuffer, allowunknown);
      if (outputsize > 0) OUT.write((const char *) outputbuffer, outputsize);                        
    }
      /*
      empty = true;                   
      int begin = 0;      
      for (int i = 0; i < line.length(); i++) {
      	  if ((line[i] == ' ') || (i == line.length() - 1)) {
          	  string word;
          	  if (line[i] == ' ') {
          	  	word  = string(line.begin() + begin, line.begin() + i);
          	  } else {
			   	word  = string(line.begin() + begin, line.begin() + i + 1);
          	  }
          	  begin = i+1;
          	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ")) {
          	    empty = false;
          	    if (classes.count(word) == 0) {
  	        		cerr << "ERROR: Unknown word '" << word << "', does not occur in model" << endl;
  	        		exit(6);          	    	
          	    }
          	  	unsigned int cls = classes[word];
          	  	int length = 0;
  	        	const unsigned char * byterep = inttobytes(cls, length);
  	        	if (length == 0) {
  	        		cerr << "INTERNAL ERROR: Error whilst encoding '" << word << "' (class " << cls << "), length==0, not possible!" << endl;
  	        		exit(13);
  	        	}  	        		
  	        	//cerr << "writing " << word << " as " << cls << " in " << length << " bytes" << endl;
  	        	OUT.write((const char *) byterep, length);
  	        	delete byterep;
  	        	OUT.write(&zero, sizeof(char)); //write separator 
          	  }			 
          }
      }
     */     
        
    
    
        
    cerr << "Encoded " << linenum << " lines" << endl;
	IN.close();
	OUT.close();
}
