#include <alignmodel.h>

using namespace std;

double CoocAlignmentModel::cooc( const multiset<uint32_t> & sourceindex, const multiset<uint32_t> & targetindex, const double threshold) {    
    //Jaccard co-occurrence    
    int intersectioncount = 0;    
  	
  	if ((threshold > 0) && (mode == JACCARD)) {
  		//pre-computation, computing best possible cooc score based on set sizes only, return 0 when threshold not attained, saving processing time (especially on large sets)
  		if (sourceindex.size() <= targetindex.size()) {  		
  			if (((double) sourceindex.size() /  targetindex.size() ) < threshold) return 0;
  		} else {
  			if (((double) targetindex.size() / sourceindex.size() ) < threshold) return 0;
  		} 	
  	}  
  
    multiset<uint32_t>::const_iterator sourceiter = sourceindex.begin();    
    multiset<uint32_t>::const_iterator targetiter = targetindex.begin();
    
    while ((sourceiter !=sourceindex.end()) && (targetiter!=targetindex.end())) {
        if (*sourceiter < *targetiter) { 
            sourceiter++;
        } else if (*targetiter < *sourceiter) {
            targetiter++;
        } else {  //equal
            intersectioncount++;
            sourceiter++;
            targetiter++;
        }
    }
    
    if (mode == JACCARD) {
    	//cerr << "union=" << unioncount << " intersection=" << intersectioncount << " ";
    	const int unioncount = (sourceindex.size() + targetindex.size()) - intersectioncount;
    	return (double) intersectioncount / unioncount;
    } else if (mode == DICE) {
    	return (double) ((2*intersectioncount) / (sourceindex.size() + targetindex.size()));
    } else {
    	cerr << "ERROR: No valid co-occurence metric selected!  Unable to compute!" << endl;
     	return 0;
    }
}


unsigned int CoocAlignmentModel::compute(const EncAnyGram * sourcegram, const multiset<uint32_t> & sourceindex, SelectivePatternModel & targetmodel) {        
   
    unsigned int found = 0;
    unsigned int prunedabs = 0;
    unsigned int prunedprob = 0;
    double totalcooc = 0;
    uint32_t prevsentencenumber = 0;
	unordered_set<const EncAnyGram *> targetpatterns;
    //cerr << "Processing new construction" << endl;
    if (DEBUG) cerr << "\t\tForward index yields " << sourceindex.size() << " sentence references" << endl;
    for (multiset<uint32_t>::const_iterator iter = sourceindex.begin(); iter != sourceindex.end(); iter++) {
        const uint32_t sentencenumber = *iter;        
        if (sentencenumber == prevsentencenumber) continue;
		if (targetmodel.reverseindex.count(sentencenumber) > 0) {
			if (DEBUG) cerr << "\t\t\tReverseindex for sentence " << sentencenumber << " yields " << targetmodel.reverseindex[sentencenumber].size() << " target-side patterns" << endl;
			for (vector<const EncAnyGram*>::const_iterator reviter = targetmodel.reverseindex[sentencenumber].begin(); reviter != targetmodel.reverseindex[sentencenumber].end(); reviter++) {
				const EncAnyGram* targetgram = *reviter;
				targetpatterns.insert(targetgram);
			}
		}
		prevsentencenumber = sentencenumber;
    }
	if (DEBUG) cerr << "\t\tGathered " << targetpatterns.size() << " target-side patterns for given source pattern, computing co-occurence..." << endl;
	for (unordered_set<const EncAnyGram *>::const_iterator iter = targetpatterns.begin(); iter != targetpatterns.end(); iter++) {
			const EncAnyGram* targetgram = *iter;
	        multiset<uint32_t> * targetindex;
		    if (targetgram->gapcount() == 0) {
		       targetindex = &targetmodel.ngrams[*( (EncNGram*) targetgram)].sentences;
		    } else {
		       targetindex = &targetmodel.skipgrams[*( (EncSkipGram*) targetgram)].sentences;
		    }				    
		    const double coocvalue = cooc(sourceindex, *targetindex, absthreshold);                    
		    if (coocvalue >= absthreshold) {
		    	//prune based on absolute co-occurrence value
		    	found++;
		        alignmatrix[sourcegram][targetgram] = coocvalue;
		        totalcooc += coocvalue;				       
		    } else {
		    	prunedabs++;
		    }
	}				
    if ((totalcooc > 0) && (probthreshold > 0)) {
    	//prune based on probability threshold
    	for (std::unordered_map<const EncAnyGram*, double>::const_iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
    		const double alignprob = (double) iter->second / totalcooc;
    		if (alignprob < probthreshold) {
    			alignmatrix[sourcegram].erase(iter->first);
    			prunedprob++;
    		}    		
    	}   
    	if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram); 
    }   
    if (DEBUG) cerr << "\t\t" << (found - prunedabs) << " alignments found (after pruning " << prunedabs << " on co-occurence value and " << prunedprob << " on alignment probability)" << endl;
    return found - prunedabs;
}

CoocAlignmentModel::CoocAlignmentModel(CoocMode mode, SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, const double absthreshold, const double probthreshold, bool DEBUG) {
    this->mode = mode;
    this->absthreshold = absthreshold;
    this->probthreshold = probthreshold;
    this->DEBUG = DEBUG;
    unsigned int c = 0;
    unsigned int found = 0;
    for (unordered_map<EncNGram,IndexCountData >::const_iterator iter = sourcemodel.ngrams.begin();  iter != sourcemodel.ngrams.end(); iter++) {
    	c++;
        if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (ngram) -- " << found << " alignment possibilities thus-far" << endl;
        found += compute(&iter->first, iter->second.sentences, targetmodel);
    }    
    for (unordered_map<EncSkipGram,IndexCountData >::const_iterator iter = sourcemodel.skipgrams.begin();  iter != sourcemodel.skipgrams.end(); iter++) {
    	c++;
    	if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (skipgram) -- " << found << " alignment possibilities thus-far" << endl;
        found += compute(&iter->first, iter->second.sentences, targetmodel);
    }            
}
    
void AlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT) {
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
        const EncAnyGram* sourcegram = iter->first;
        *OUT << sourcegram->decode(sourceclassdecoder) << "\t";
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram* targetgram = iter2->first;            
            *OUT << targetgram->decode(targetclassdecoder) << "\t" << iter2->second << "\t";
        }
        *OUT << endl;
    }
}



/*
EMAlignmentModel::EMAlignmentModel(IndexedPatternModel & sourcemodel, IndexedPatternModel & targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD) {
    int round = 0;    
    unsigned long c;
    double totaldivergence = 0;
    double prevavdivergence = 0;
    bool converged = false;
    
    set<int> reverseindexkeys = sourcemodel.reverse_index_keys();
            
    do {       
        round++; 
        c = 0;
        cerr << "  EM Round " << round << "... ";
        //use reverse index to iterate over all sentences
        for (set<int>::iterator iter = reverseindexkeys.begin(); iter != reverseindexkeys.end(); iter++) {
            const int key = *iter;        
            const int sourcegrams_size =  sourcemodel.reverse_index_size(key);
            //vector<EncAnyGram*> sourcegrams = sourcemodel.reverse_index(key);
            const int targetgrams_size =  targetmodel.reverse_index_size(key);
            
            if (targetgrams_size > 0 ) { //target model contains sentence?               
                //vector<EncAnyGram*> targetgrams = targetmodel.reverse_index(key);  
                cerr << key << ":" << sourcegrams_size << "x" << targetgrams_size << " ";

                //compute sentencetotal for normalisation later in count step, sum_s(p(t|s))                                      
                unordered_map<const EncAnyGram*, double> sentencetotal;                
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        sentencetotal[targetgram] += alignprob[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration
                    }
                }               
                
                //collect counts (for evidence that a targetgram is aligned to a sourcegram)
                std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
                unordered_map<const EncAnyGram*, double> total;
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);                        
                        const double countvalue = alignprob[sourcegram][targetgram] / sentencetotal[targetgram];
                        count[sourcegram][targetgram] += countvalue;
                        total[targetgram] += countvalue;
                    }                    
                }
                
                 
                double prevtransprob;
                //estimate new probabilities (normalised count is the new estimated probability)
                for (int i = 0; i < targetgrams_size; i++) {    
                    const EncAnyGram* targetgram = targetmodel.get_reverse_index_item(key,i);   
                    for (int j = 0; j < sourcegrams_size; j++) {                                            
                        const EncAnyGram* sourcegram = sourcemodel.get_reverse_index_item(key,j);
                        
                        prevtransprob = alignprob[sourcegram][targetgram];
                        const double newtransprob = (double) count[sourcegram][targetgram] / total[targetgram];
                        alignprob[sourcegram][targetgram] = newtransprob;                        
                        
                        //for computation of convergence
                        const double divergence = abs(newtransprob - prevtransprob);
                        cerr << " prevtransprob=" << prevtransprob << " ";
                        cerr << " newtransprob=" << newtransprob << " ";
                        cerr << " div=" << divergence << " ";
                    
                        totaldivergence += divergence;
                        c++;
                    }
                }
            }
        }
        const double avdivergence = (double) totaldivergence / c;
        converged = (((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence)) <= CONVERGEDTHRESHOLD);       
        prevavdivergence = avdivergence;
        cerr << " average divergence = " << avdivergence << ", alignprob size = " << alignprob.size() << endl;
    } while (!converged);    
}
*/

