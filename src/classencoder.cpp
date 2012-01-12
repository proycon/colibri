#include "classencoder.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>

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
}

ClassEncoder::ClassEncoder(const string & filename) {
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
                  //cerr << "CLASS=" << cls << " WORD=" << word << endl;
              }
              
          }
        }        
        IN.close();  
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
              	  string word = string(line.begin() + begin, line.begin() + i);
              	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " "))
              	  	freqlist[word]++;
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
        
        int cls = 1; //one is reserved
        for (multimap<const int,const string>::iterator iter = revfreqlist.begin(); iter != revfreqlist.end(); iter++) {
        	cls++;
        	while (!validclass(cls)) cls++;
        	classes[iter->second] = cls;
        }
}

void ClassEncoder::save(const string & filename) {
	ofstream OUT;
	OUT.open( filename.c_str() );
	for (std::unordered_map<std::string,unsigned int>::iterator iter = classes.begin(); iter != classes.end(); iter++) {
		OUT << iter->second << '\t' << iter->first << endl;
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


void ClassEncoder::encodefile(const std::string & inputfilename, const std::string & outputfilename) {
	const char zero = 0;
	const char one = 1;
	ofstream OUT;
	ifstream IN;
	OUT.open(outputfilename.c_str(), ios::out | ios::binary);	
	IN.open(inputfilename.c_str());
	while (IN.good()) {	
      string line;
      getline(IN, line);              
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
          	  	unsigned int cls = classes[word];
          	  	int length = 0;
  	        	const unsigned char * byterep = inttobytes(cls, length);
  	        	//cerr << "writing " << word << " as " << cls << " in " << length << " bytes" << endl;
  	        	OUT.write((const char *) byterep, length);
  	        	delete byterep; 
          	  }
			 OUT.write(&zero, sizeof(char)); //write separator
          }
      }
      OUT.write(&one, sizeof(char)); //newline          
      OUT.write(&zero, sizeof(char)); //write separator          
    }        
	IN.close();
	OUT.close();
}
