#include <lm.h>


using namespace std;

LanguageModel::LanguageModel(const std::string & filename, ClassEncoder & encoder, ClassDecoder * classdecoder, bool debug) {
    this->DEBUG = debug; 
    this->classdecoder = classdecoder;
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
                for (unsigned int i = 0; i  <= line.length(); i++) {
                    if ((line[i] == '\t') || (line[i] == '\n') || (i == line.length())) {
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
                        ngrams[ngram] = atof(logprob_s.c_str()) * log(10); //* log(10) does log10 to log_e conversion
                        hasunk = true;
                        if (DEBUG) {
                            cerr << " Adding UNKNOWN to LM: " << (int) ngram.n() << "\t" <<  ngramcontent << "\t" << ngrams[ngram] << endl;
                        }
                    } else {
                        EncNGram ngram = encoder.input2ngram(ngramcontent,  true);
                        if (!ngram.unknown()) {
                            ngrams[ngram] = atof(logprob_s.c_str()) * log(10); //* log(10) does log10 to log_e conversion
                            if (!backofflogprob_s.empty()) {
                                backoff[ngram] = atof(backofflogprob_s.c_str()) * log(10); //* log(10) does log10 to log_e conversion
                                if (DEBUG) cerr << " Adding to LM: " << (int) ngram.n() << "\t" <<  ngramcontent << "\t" << ngrams[ngram] << "\t" << backoff[ngram] << endl;
                            } else {
                                if (DEBUG) cerr << " Adding to LM: " << (int) ngram.n() << "\t" << ngramcontent << "\t" << ngrams[ngram] << endl;
                            }
                        }
                    }
                } else {
                    cerr << "WARNING: Ignoring line: " << line << endl;
                }
            } else {
                cerr << "WARNING: Don't know what to do with line: " << line << endl;
            }
        }
        
    }
    f.close();
    
    if (!hasunk) {
        cerr << "ERROR: Language Model has no value <unk>, make sure to generate SRILM model with -unk parameter" << endl;
        exit(3);
    }
}


double LanguageModel::score(const EncNGram * ngram, const EncNGram * history) { //returns logprob (base 10)
    
    if (DEBUG) {
        if (history == NULL) {
            cerr << "\t\t\tLM DEBUG: score(): history=NULL ngram='" << ngram->decode(*classdecoder) << "'" << endl;
        } else {
            cerr << "\t\t\tLM DEBUG: score(): history='" << history->decode(*classdecoder) << "' ngram='" << ngram->decode(*classdecoder) << "'" << endl;
        }
    }  
    double result = 0;
    const int n = ngram->n();       
    for (int i = 0; i < n; i++) {
        EncNGram * word = ngram->slice(i,1);
        EncNGram * newhistory = NULL;
        if ((i >= order-1) || ((i > 0) && (history == NULL))) {
            //new history can be fetched from ngram alone
            int begin = i - (order - 1);
            if (begin < 0) begin = 0;              
            newhistory = ngram->slice(begin,i - begin);            
        } else if (history != NULL) {
            //new history has to be fetched from old history and ngram
            EncNGram * slice = NULL;
            if (i > 0) slice = ngram->slice(0,i);
            const int leftover = order - 1 - i;
            if (leftover > 0) {
                EncNGram * historyslice = history->slice(history->n() - leftover, leftover);
                if (slice == NULL) {
                    newhistory = historyslice;
                } else {   
                    newhistory = new EncNGram(*historyslice + *slice);
                    delete historyslice;
                }            
            } 
            if (slice != NULL) delete slice;
        }
        
        result += scoreword(word, newhistory);
         
        delete word;
        if (newhistory != NULL) delete newhistory;
    }
    if (DEBUG) cerr << "\t\t\tLM DEBUG: score() = " << result << endl;
    return result; 
}

double LanguageModel::scoreword(const EncNGram * word, const EncNGram * history) {

    if (DEBUG) {
        if (history == NULL) {
            cerr << "\t\t\tLM DEBUG: scoreword(): history=NULL word='" << word->decode(*classdecoder) << "'" << endl;
        } else {
            cerr << "\t\t\tLM DEBUG: scoreword(): history='" << history->decode(*classdecoder) << "' word='" << word->decode(*classdecoder) << "'" << endl;
        }
    }  

    const EncNGram * lookup;
    
    if (history != NULL) {
        lookup = new EncNGram( *history + *word);
    } else {
        lookup = word;
    }
    const int n = lookup->n();
    unordered_map<EncNGram, double>::iterator iter = ngrams.find(*lookup);
         
    if (iter != ngrams.end()) {
        if (DEBUG) cerr << "\t\t\tLM DEBUG: scoreword(): Found " << n << "-gram, score=" << iter->second << endl;
        if (history != NULL) delete lookup;        
        return iter->second;        
    } else {
        if (DEBUG) cerr << "\t\t\tLM DEBUG: scoreword(): " << n << "-gram not found. Backing off..." << endl;        
    }
    
    if (history != NULL) delete lookup;
    
         
    //not found, back-off    
    double result = 0;

    double backoffweight = 0; //backoff weight will be 0 if not found 
    /*
    Not all N-grams in the model file have backoff weights. The highest order N-grams do not need a backoff weight. For lower order N-grams backoff weights are only recorded for those that appear as the prefix of a longer N-gram included in the model. For other lower order N-grams the backoff weight is implicitly 1 (or 0, in log representation).
    */    

    if (history != NULL) {
        //ngram not found: back-off: alpha(history) * p(n-1)    


        

        const int history_n = history->n();
                    
        EncNGram * newhistory;        
    
        
        if (history_n > 1) { 
            newhistory = history->slice(1,history_n-1);    
        } else {
            newhistory = NULL;
        }
        unordered_map<EncNGram, double>::iterator iter = backoff.find(*history);
        if (iter != backoff.end()) {
            if (DEBUG) cerr << "\t\t\tLM DEBUG: backoffpart=" << iter->first.decode(*classdecoder) << " = " << iter->second << endl;
            backoffweight = iter->second;        
        }
        //if (DEBUG) cerr << "LM DEBUG: scoreword(): Backoffpart=" << backoffpart->decode(*classdecoder) << endl;
              
        
        if (DEBUG) cerr << "\t\t\tLM DEBUG: scoreword(): Backoffweight=" << backoffweight << endl;
        result = backoffweight + scoreword(word, newhistory);
        if (DEBUG) cerr << "\t\t\tLM DEBUG: scoreword() =" << result << endl;
        if (newhistory != NULL) delete newhistory;
        
    } else {
        cerr << "INTERNAL ERROR: LanguageModel::scoreword() ... unigram not found, and no history.. this should not happen" << endl;
        exit(6);
    }    
    return result;
}



/*double LanguageModel::score(EncData & data, bool fullsentence) {
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
            if (DEBUG) cerr << "LM DEBUG: Found " << (int) ngram.n() << "-gram, score=" << iter->second << endl; 
            return iter->second;
        }
        
        if (n == 1) {
            cerr << "INTERNAL ERROR: Unexpected unigram passed to LM for scoring. This should not happen" << endl;
            exit(6);
        }

        //ngram not found: back-off: alpha(history) * p(n-1)
        EncNGram * history = ngram.slice(0,n-1);
        EncNGram * head = ngram.slice(1,n-1);

        double backoffweight = 0; //backoff weight will be 0 if not found 

        
        for (unordered_map<EncNGram, double>::iterator iter = backoff.find(*history); iter != backoff.end(); iter++) {
            backoffweight = iter->second;
        }
        if (DEBUG) cerr << "LM DEBUG: Backoffweight " << backoffweight <<", backing off.." << endl;
        result = backoffweight + score(*head);
        if (DEBUG) cerr << "LM DEBUG: Result=" << result << endl;
        
        delete history;
        delete head;
        return result;    
    }
}
*/
