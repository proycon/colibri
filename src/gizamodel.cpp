#include "gizamodel.h"

using namespace std;


GizaModel::GizaModel(const string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder) {
    IN = new ifstream(filename.c_str(), ios::in);
    this->sourceencoder = sourceencoder;
    this->targetencoder = targetencoder;    
    sentenceindex = 0;
}

GizaModel::~GizaModel() {
    IN->close();
    delete IN;
}

GizaSentenceAlignment::GizaSentenceAlignment(const string & sourceline,const string & targetline,ClassEncoder * sourceencoder, ClassEncoder * targetencoder, const int index) {
    parsesource(sourceline, sourceencoder);
    parsetarget(targetline, targetencoder);
    this->index = index;
}

/*GizaSentenceAlignment::GizaSentenceAlignment(const EncNGram * source, const EncNGram * target, const int index) {
    this->source = new EncNGram(*source);
    this->target = new EncNGram(*target);
    this->index = index;    
}*/
     
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
      return GizaSentenceAlignment(sourceline, targetline, sourceencoder, targetencoder, sentenceindex);    
}

void GizaSentenceAlignment::parsetarget(const string & line, ClassEncoder * classencoder) {
  unsigned char buffer[65536];
  
  int size = classencoder->encodestring(line, buffer, true);
  target = new EncNGram(buffer,size-1);  
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
      	        	         	  
      	  if ((word.length() > 0) && (word != "\r") && (word != "\t") && (word != " ") && (word != "NULL")) {
      	    if (!inalignment) {
      	        if (!cleanline.empty()) cleanline += " ";
      	  	    cleanline += word;
      	  	    sourceindex++;
      	  	} else {
      	  	    int targetindex = atoi(word.c_str());
      	  	    //istringstream(word) >> targetindex;
      	  	    alignment.insert( pair<const unsigned char, const unsigned char>(sourceindex, (unsigned char) targetindex) );
      	  	    //alignment[sourceindex].push_back((unsigned char) targetindex);      	  	          	  	    
      	  	}
      	  }
      	  begin = i+1; 
      }
  }
  
  //cerr << endl << "SOURCE=[" << cleanline << "]" << endl;
  int size = classencoder->encodestring(cleanline, buffer, true);
  source = new EncNGram(buffer,size-1);  
}


GizaSentenceAlignment GizaSentenceAlignment::intersect(const GizaSentenceAlignment & other) {
    GizaSentenceAlignment intersection = *this;
     
    for (multimap<const unsigned char, const unsigned char>::const_iterator iter = alignment.begin(); iter != alignment.end(); iter++) {
        unsigned char sourceindex = iter->first;
        unsigned char targetindex = iter->second;
        if (other.alignment.count(targetindex)) {
            for (multimap<const unsigned char, const unsigned char>::const_iterator iter2 = other.alignment.lower_bound(targetindex); iter2 != other.alignment.upper_bound(targetindex); iter2++) {
                if (iter2->second == sourceindex) {
                    intersection.alignment.insert( pair<const unsigned char, const unsigned char>(sourceindex, (unsigned char) targetindex) );
                }
            } 
        }
    }
    return intersection;
}  


GizaSentenceAlignment GizaSentenceAlignment::unify(const GizaSentenceAlignment & other) {
    GizaSentenceAlignment unified = *this;     
    for (multimap<const unsigned char, const unsigned char>::const_iterator iter = alignment.begin(); iter != alignment.end(); iter++) {
        unsigned char sourceindex = iter->first;
        unsigned char targetindex = iter->second;
        unified.alignment.insert( pair<const unsigned char, const unsigned char>(sourceindex, (unsigned char) targetindex) );        
    }

    for (multimap<const unsigned char, const unsigned char>::const_iterator iter = other.alignment.begin(); iter != other.alignment.end(); iter++) {
        unsigned char sourceindex = iter->first;
        unsigned char targetindex = iter->second;        
        unified.alignment.insert( pair<const unsigned char, const unsigned char>(targetindex, (unsigned char) sourceindex) );
    }    
    return unified;
}  

void GizaSentenceAlignment::out(std::ostream* OUT, ClassDecoder & sourcedecoder, ClassDecoder & targetdecoder) {
    //OUT << "<html><head><title>Word Alignments</title></head><body>" << endl;
    *OUT << "<table>";
    *OUT << "<tr><td></td>";
        
    for (int i = 0; i < source->n(); i++) {
        *OUT << "<td>";         
        EncNGram * unigram = source->slice(i, 1);
        *OUT << unigram->decode(sourcedecoder);
        delete unigram;
        *OUT << "</td>";
    }
    
    *OUT << "</tr>" << endl;
    
    //init
     
    map<unsigned char, map<unsigned char, bool>> matrix;     
    for (int i = 0; i < source->n(); i++) {
        for (int j = 0; j < target->n(); j++) {
            matrix[i][j] = false;
        }
    } 
       

    
    for (multimap<const unsigned char, const unsigned char>::const_iterator iter = alignment.begin(); iter != alignment.end(); iter++) {
        unsigned char sourceindex = iter->first;
        unsigned char targetindex = iter->second;
        matrix[sourceindex][targetindex] = true;   
    }

    
    
    
    for (int i = 0; i < target->n(); i++) {    
        *OUT << "<tr><td>";   
        EncNGram * unigram = target->slice(i, 1);
        *OUT << unigram->decode(targetdecoder);
        delete unigram;   
        *OUT << "</td>";        
        for (int j = 0; j < source->n(); j++) {        
            if (matrix[j][i]) {
              *OUT << "<td style=\"background: black\"></td>" << endl;
            } else {
                *OUT << "<td></td>" << endl;
            } 
        }
        *OUT << "</tr>";
    }        
    *OUT << "</table>" << endl;
}
