#include <lm.h>

using namespace std;

LanguageModel::LanguageModel(const std::string & filename, ClassEncoder & encoder) { 
    order = 0;
    bool hasunk = false;
    ifstream f;    
    f.open(filename.c_str(), ios::in);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }    
    while (!f.eof()) {               
        string line;
        getline(f, line);                
        if (line == "\\data\\") {
            order = 0;
        } else if (line == "\\1-grams:") { //bit inelegant, but simplest
            order = 1;
        } else if (line == "\\2-grams:") {
            order = 2;
        } else if (line == "\\3-grams:") {
            order = 3;            
        } else if (line == "\\4-grams:") {
            order = 4;
        } else if (line == "\\5-grams:") {
            order = 5;            
        } else if (line == "\\6-grams:") {
            order = 6;            
        } else if (line == "\\7-grams:") {
            order = 7;            
        } else if (line == "\\8-grams:") {
            order = 8;            
        } else if (line == "\\9-grams:") {
            order = 9;                        
        } else if (!line.empty()) {
            if (order == 0) {
              if (line.substr(0,5) == "ngram") {
                string n_s = line.substr(6,1);
                string v_s = line.substr(8);
                int n = atoi(n_s.c_str());
                int v = atoi(v_s.c_str());
                total[n] = v;
              }   
            } else if (order > 0) {
                string logprob_s = "";
                string backofflogprob_s = "";
                string ngramcontent = "";
                int fields = 0;
                int begin = 0;
                for (int i = 0; i  < line.length(); i++) {
                    if (line[i] == '\t') {
                        if (fields == 0) {
                            logprob_s = line.substr(begin, i - begin);
                        } else if (fields == 1) {
                            ngramcontent = line.substr(begin, i - begin);
                        } else if (fields == 2) {
                            backofflogprob_s = line.substr(begin, i - begin);
                        }
                        begin = i + 1;
                        fields++;
                    }
                }
                if ((!logprob_s.empty()) && (!ngramcontent.empty())) {
                    if (ngramcontent == "<unk>") {
                        unsigned char unknownclass = 2;
                        EncNGram ngram = EncNGram(&unknownclass, 1);
                        ngrams[ngram] = atof(logprob_s.c_str());
                        hasunk = true;
                    } else {
                        EncNGram ngram = encoder.input2ngram(ngramcontent,  true);
                        if (!ngram.unknown()) {
                            ngrams[ngram] = atof(logprob_s.c_str());
                            if (!backofflogprob_s.empty()) {
                                backoff[ngram] = atof(backofflogprob_s.c_str());
                            }
                        }
                    }
                }
            }
        }
        
    }
    f.close();
    
    if (!hasunk) {
        cerr << "ERROR: Language Model has no value <unk>, make sure to generate SRILM model with -unk parameter" << endl;
        exit(3);
    }
}


double LanguageModel::score(EncData & data, bool fullsentence) {
    //TODO: handle fullsentence: add <s></s> wrappers
    if (data.length() <= order) {
        return score(EncNGram(data));
    } else {
        double result = 0;
        for (int begin = 0; begin < data.length() - order; begin++) {
            EncNGram * ngram = data.slice(begin, order);
            result += score(*ngram);
            delete ngram;        
        }
        return result;
    } 
}


double LanguageModel::score(EncNGram ngram) {
    double result;
    const int n = ngram.n();
    if (n > order) {
        result = 0;
        for (int begin = 0; begin < n - order; begin++) {
            EncNGram * subngram = ngram.slice(begin, order);
            result += score(*subngram);
            delete subngram;        
        }    
        return result;
    } else {
        //if n-gram exists, return probability    
        for (unordered_map<EncNGram, double>::iterator iter = ngrams.find(ngram); iter != ngrams.end(); iter++) {
            return iter->second;
        }

        //ngram not found: back-off: alpha(history) * p(n-1)
        EncNGram * history = ngram.slice(0,n-1);
        EncNGram * head = ngram.slice(1,n-1);

        double backoffscore = 0;
        for (unordered_map<EncNGram, double>::iterator iter = backoff.find(*history); iter != ngrams.end(); iter++) {
            backoffscore = iter->second;
        }
        if (backoffscore != 0) {
            //use backoffscore from history
            result = backoffscore + score(*head);
        } else {
            //not found, backoff some more
            result = score(*head); 
        }

        delete history;
        delete head;
        return result;    
    }
}
