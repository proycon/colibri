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
	length = ceil(cls / 256);	
	unsigned char * byterep = new unsigned char[length];
	int i = 0;
    do {
    	int r = cls % 256;
    	byterep[i++] = (unsigned char) r;
    	cls = cls / 256;
    } while (cls >= 256);
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
              	  freqlist[word]++;
              }
              
          }
        }        
        IN.close();
        
        //sort by occurrence count  using intermediate representation
        multimap<const int, const string> revfreqlist;
        for (unordered_map<const string,int>::iterator iter = freqlist.begin(); iter != freqlist.end(); iter++) {
        	revfreqlist.insert( pair<const int,const string>(iter->second, iter->first) );
        } 
        
        freqlist.clear();
        
        int cls = 2; //one is reserved
        for (multimap<const int,const string>::iterator iter = revfreqlist.begin(); iter != revfreqlist.end(); iter++) {
        	while (!validclass(cls)) cls++;
        	classes[iter->second] = iter->first;
        }
}

void ClassEncoder::save(const string & filename) {
	ofstream OUT;
	OUT.open( filename.c_str() );
	for (std::unordered_map<std::string,unsigned int>::iterator iter = classes.begin(); iter != classes.end(); iter++) {
		OUT << iter->first << '\t' << iter->second << endl;
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
	const char one = 0;
	ofstream OUT;
	ifstream IN;
	OUT.open(outputfilename.c_str(), ios::out | ios::binary);	
	IN.open(inputfilename.c_str());
	while (IN.good()) {	
      string line;
      getline(IN, line);              
      int begin = 0;
      for (int i = 0; i < line.size(); i++) {
          if ((line[i] == ' ') || (i == line.size() - 1)) {              
          	  const string word = string(line.begin() + begin, line.begin() + i);
          	  begin = i;
          	  if (word == "\n") {
          	  	if (i < line.size() - 1) OUT.write(&one, sizeof(char)); //write newline
          	  } else {
          	  	unsigned int cls = classes[word];
          	  	int length = 0;
  	        	const unsigned char * byterep = inttobytes(cls, length);
  	        	OUT.write((const char *) byterep, length);
  	        	delete byterep; 
          	  }
          }
          if (i < line.size() - 1) OUT.write(&zero, sizeof(char)); //write separator
      }
    }        
	IN.close();
	OUT.close();
}
