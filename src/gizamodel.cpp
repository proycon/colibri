#include "gizamodel.h"

using namespace std;


GizaModel::GizaModel(const string & filename, ClassEncoder * classencoder, ClassDecoder * classdecoder) {
    IN = new ifstream(filename.c_str(), ios::in);
    this->classencoder = classencoder;
    this->classdecoder = classdecoder;    
    sentenceindex = 0;
}

GizaModel::~GizaModel() {
    IN->close();
    delete IN;
}

GizaSentenceAlignment::GizaSentenceAlignment(const string & sourceline,const string & targetline, ClassEncoder * classencoder, const int index) {
    parsesource(sourceline, classencoder);
    parsetarget(targetline, classencoder);
    this->index = index;
}

/*GizaSentenceAlignment::GizaSentenceAlignment(const EncNGram * source, const EncNGram * target, const int index) {
    this->source = new EncNGram(*source);
    this->target = new EncNGram(*target);
    this->index = index;    
} */
     
GizaSentenceAlignment::GizaSentenceAlignment(const GizaSentenceAlignment& ref) {
    index = ref.index;
    source = new EncNGram(*(ref.source));
    target = new EncNGram(*(ref.target));
    alignment = ref.alignment; //not sure if this will work?
}

GizaSentenceAlignment::~GizaSentenceAlignment() {
    if (source != NULL) delete source;
    if (target != NULL) delete target;
}

GizaSentenceAlignment GizaModel::readsentence() {
      sentenceindex++;
      string line;
      
      getline(*IN, line); 
      if (line[0] != '#') {
        cerr << "Error parsing GIZA++ Alignment, excepted new fragment, found: " << line << endl;
      }       
      
      string targetline;
      string sourceline;
      getline(*IN, targetline);
      getline(*IN, sourceline);
      return GizaSentenceAlignment(sourceline, targetline, classencoder, sentenceindex);    
}

void GizaSentenceAlignment::parsetarget(const string & line, ClassEncoder * classencoder) {
  unsigned char buffer[65536];
    
  int size = classencoder->encodestring(line, buffer, true);
  target = new EncNGram(buffer,size);
}
  
void GizaSentenceAlignment::parsesource(const string & line, ClassEncoder * classencoder) {
  string cleanline = ""; 
  unsigned char buffer[65536];
  bool inalignment = false;
   
  int begin = 0;
  unsigned char sourceindex = 1; //first is 1
  for (int i = 0; i < line.size(); i++) {
      if ((line[i] == ' ') || (i == line.size() - 1)) {
      	  int offset = 0;      	        	  
      	  if (i == line.size() - 1) offset = 1;
      	   
      	  string word =  string(line.begin() + begin, line.begin() + i + offset); 
      	  if (word == "})") {
      	    inalignment = false;
      	    begin = i + 1;
      	    continue;
      	  } 
      	  if (word == "({") {
      	    inalignment = true;
      	    begin = i + 1;
      	    continue;   
      	  } 
      	        	         	  
      	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ")) {
      	    if (!inalignment) {
      	        if (!cleanline.empty()) cleanline += " ";
      	  	    cleanline += word;
      	  	    sourceindex++;
      	  	} else {
      	  	    int targetindex = atoi(word.c_str());
      	  	    //istringstream(word) >> targetindex;
      	  	    alignment[sourceindex].push_back((unsigned char) targetindex);      	  	          	  	    
      	  	}
      	  }
      	  begin = i+1; 
      }
  }
  
  int size = classencoder->encodestring(cleanline, buffer, true);
}
