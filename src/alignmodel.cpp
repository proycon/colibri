#include <alignmodel.h>

using namespace std;

double CoocAlignmentModel::cooc( multiset<uint32_t> & sourceindex, multiset<uint32_t> & targetindex) {    
    //Jaccard co-occurrence    
    int intersectioncount = 0;    
    
    multiset<uint32_t>::iterator sourceiter = sourceindex.begin();    
    multiset<uint32_t>::iterator targetiter = targetindex.begin();
    
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
    const int unioncount = (sourceindex.size() + targetindex.size()) - intersectioncount; 
    //cerr << "union=" << unioncount << " intersection=" << intersectioncount << " ";
    return (double) intersectioncount / unioncount;
}


int CoocAlignmentModel::compute(const EncAnyGram * sourcegram, multiset<uint32_t> & sourceindex, DoubleIndexedGraphPatternModel & targetmodel) {        
    int c = 0;
    double bestcooc = 0;
    //cerr << "Processing new construction" << endl;
    for (multiset<uint32_t>::iterator iter = sourceindex.begin(); iter != sourceindex.end(); iter++) {
        const uint32_t sentencenumber = *iter;
      	//cerr << "Reverseindex lookup: " << sentencenumber << endl;        
		if (targetmodel.reverseindex.count(sentencenumber) > 0) {
			//cerr << "\tFound" << endl;
			c += targetmodel.reverseindex[sentencenumber].size();
			//cerr << "\tPatterns: " << targetmodel.reverseindex[sentencenumber].size() << endl;
			for (vector<const EncAnyGram*>::iterator reviter = targetmodel.reverseindex[sentencenumber].begin(); reviter != targetmodel.reverseindex[sentencenumber].end(); reviter++) {
					
					const EncAnyGram* targetgram = *reviter;
			        multiset<uint32_t> * targetindex;
				    if (targetgram->gapcount() == 0) {
				       targetindex = &targetmodel.ngrams[*( (EncNGram*) targetgram)].sentences;
				    } else {
				       targetindex = &targetmodel.skipgrams[*( (EncSkipGram*) targetgram)].sentences;
				    }				    
				    const double coocvalue = cooc(sourceindex, *targetindex);        
				    if ((relthreshold) && (coocvalue > bestcooc)) bestcooc = coocvalue;            
				    if (coocvalue >= absthreshold) {
				    	//cerr << "\t\tRegistered cooc: " << coocvalue << endl;                
				        alignprob[sourcegram][targetgram] = coocvalue;				       
				    }
			}				
	        if (relthreshold) {
            //TODO: prune based on relative threshold
	        }   
	    }			        
    }    
    return c;
}

CoocAlignmentModel::CoocAlignmentModel(DoubleIndexedGraphPatternModel & sourcemodel, DoubleIndexedGraphPatternModel & targetmodel, const double absthreshold, const double relthreshold) {
    this->absthreshold = absthreshold;
    this->relthreshold = relthreshold;
    unsigned int c = 0;
    for (unordered_map<EncNGram,IndexCountData >::iterator iter = sourcemodel.ngrams.begin();  iter != sourcemodel.ngrams.end(); iter++) {
    	c++;
        if (c % 10000 == 0) cerr << "\t@" << c << endl;
        compute(&iter->first, iter->second.sentences, targetmodel);
    }    
    for (unordered_map<EncSkipGram,IndexCountData >::iterator iter = sourcemodel.skipgrams.begin();  iter != sourcemodel.skipgrams.end(); iter++) {
    	c++;
    	if (c % 10000 == 0) cerr << "\t@" << c << endl;
        compute(&iter->first, iter->second.sentences, targetmodel);
    }            
}
    
void AlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT) {
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignprob.begin(); iter != alignprob.end(); iter++) {
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

