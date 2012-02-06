#include <alignmodel.h>

using namespace std;

double CoocAlignmentModel::cooc( const multiset<uint32_t> & sourceindex, const multiset<uint32_t> & targetindex, const double heurthreshold) {    
    //Jaccard co-occurrence    
    int intersectioncount = 0;    
  	
  	if ((heurthreshold > 0) && (mode == JACCARD)) {
  		//pre-computation, computing best possible cooc score based on set sizes only, when heuristic threshold not attained, saving processing time on computation (especially on large sets) and simply return the highest possibly attainable score (never above the given threshold)
  		if (sourceindex.size() <= targetindex.size()) {  		
  			const double coocheur = (double) sourceindex.size() /  targetindex.size();
  			if (coocheur < heurthreshold) return coocheur;
  		} else {
  			const double coocheur = (double) targetindex.size() /  sourceindex.size();
  			if (coocheur < heurthreshold) return coocheur;
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
    double bestcooc = 0;
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
		    } else {
		    	prunedabs++;
		    }
		    totalcooc += coocvalue;				
		    if (coocvalue > bestcooc) bestcooc = coocvalue;
		    
	}				
    if ((totalcooc > 0) && (normalize || bestonly || probthreshold > 0)) {
    	//normalisation and pruning step (based on probability threshold)
    	for (std::unordered_map<const EncAnyGram*, double>::const_iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
    		const double alignprob = (double) iter->second / totalcooc;
			if ((alignprob < probthreshold) || (bestonly && iter->second < bestcooc)) {
				//prune
    			alignmatrix[sourcegram].erase(iter->first);
    			prunedprob++;    			
    		} else if (normalize) {
    			//normalise
    			alignmatrix[sourcegram][iter->first] = alignprob;
    		}    		    		 
    	}   
    	if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram); 
    }   
    if (DEBUG) cerr << "\t\t" << (found - prunedprob) << " alignments found (after pruning " << prunedabs << " on co-occurence value and " << prunedprob << " on alignment probability)" << endl;
    return found - prunedprob;
}

CoocAlignmentModel::CoocAlignmentModel(CoocMode mode, SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, const double absthreshold, const double probthreshold, bool bestonly, bool normalize, bool DEBUG) {
    this->mode = mode;
    this->absthreshold = absthreshold;
    this->probthreshold = probthreshold;
    this->normalize = normalize;
    this->bestonly = bestonly;
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
        map<double, const EncAnyGram*> sorted;        
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	sorted[-1 * iter2->second] = iter2->first;            
        }
        for (map<double, const EncAnyGram*>::iterator iter2 = sorted.begin(); iter2 != sorted.end(); iter2++) {
			const EncAnyGram* targetgram = iter2->second;            
            *OUT << targetgram->decode(targetclassdecoder) << "\t" << (-1 * iter2->first) << "\t";
        }            
        *OUT << endl;
    }
}



EMAlignmentModel::EMAlignmentModel(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, bool bestonly, bool DEBUG) {
    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
    //initialise uniformly
    cerr << "  Initialisation step" << endl; 
    for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel.reverseindex.begin(); reviter_source != sourcemodel.reverseindex.end(); reviter_source++) {
    	uint32_t sentence = reviter_source->first;
		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
		if (targetmodel.reverseindex.count(sentence) > 0) {
			vector<const EncAnyGram*> * targetpatterns = &targetmodel.reverseindex[sentence];
			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
			double v = (double) 1 / targetpatterns->size();
			for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
    			const EncAnyGram * targetgram = *targetiter;
                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
    			    const EncAnyGram * sourcegram = *sourceiter;
					alignmatrix[sourcegram][targetgram] = v;
    			}
    		}
		}
    }
          
            
            
    do {       
        round++; 
        c = 0;        
        cerr << "  EM Round " << round << "... ";
        
		std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
		unordered_map<const EncAnyGram*, double> total;
        
        
        //ESTIMATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel.reverseindex.begin(); reviter_source != sourcemodel.reverseindex.end(); reviter_source++) {        		
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel.reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel.reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal;      
        			for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
        				const EncAnyGram * targetgram = *targetiter;
        				for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
        					const EncAnyGram * sourcegram = *sourceiter;
        					sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration 
        				} 
        			}
        			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		                    const double countvalue = alignmatrix[sourcegram][targetgram] / sentencetotal[targetgram];
		                    count[sourcegram][targetgram] += countvalue;
		                    total[sourcegram] += countvalue;
		                }
		            }
		       
        		}	
		} //end loop over corpus
		
        double prevtransprob;                
        double totaldivergence = 0;
        //MAXIMISATION STEP: improved model  update probability estimates (Maximum Likelihood Estimation) (normalised count is the new estimated probability)
        for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
        	const EncAnyGram * sourcegram = sourceiter->first;
        	for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
        		const EncAnyGram * targetgram = targetiter->first;
        		
				prevtransprob = targetiter->second;
                const double newtransprob = (double) count[sourcegram][targetgram] / total[sourcegram];
                alignmatrix[sourcegram][targetgram] = newtransprob;
                
                //for computation of convergence
                const double divergence = abs(newtransprob - prevtransprob);
                if (DEBUG) {/*
	                cerr << " prevtransprob=" << prevtransprob << " ";
	                cerr << " newtransprob=" << newtransprob << " ";
	                cerr << " div=" << divergence << " " << endl;*/
                }
                totaldivergence += divergence;
                c++;        		
        	}
        }
		
        const double avdivergence = (double) totaldivergence / c;
        converged = ((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence) <= CONVERGEDTHRESHOLD);               
        cerr << "   average divergence = " << avdivergence << ", delta with prev divergence = " << abs(avdivergence - prevavdivergence) << " > " << CONVERGEDTHRESHOLD << ", alignprob size = " << alignmatrix.size() << endl;
        prevavdivergence = avdivergence;
    } while (!converged);    
    
    if ((probthreshold > 0) || (bestonly)) {
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			double best = 0;
			for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
				const EncAnyGram * targetgram = targetiter->first;
				if (targetiter->second < probthreshold) {
					alignmatrix[sourcegram].erase(targetgram);
				} else if (bestonly && targetiter->second > best) best = targetiter->second;
			}
			if (bestonly) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < best) alignmatrix[sourcegram].erase(targetgram);
				}
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}


void AlignmentModel::intersect(AlignmentModel * reversemodel, double probthreshold) {
	 //Compute intersection with reverse model
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;
		for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			if ((reversemodel->alignmatrix.count(targetgram) > 0) && (reversemodel->alignmatrix[targetgram].count(sourcegram) > 0) && (reversemodel->alignmatrix[targetgram][sourcegram] * targetiter->second >= probthreshold)) {
				alignmatrix[sourcegram][targetgram] = reversemodel->alignmatrix[targetgram][sourcegram] * targetiter->second;
			} else {
				//prune
				alignmatrix[sourcegram].erase(targetgram);
			} 			
		}		
		if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);	
	}	 
}	 
