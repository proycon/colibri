#include <alignmodel.h>

using namespace std;

const EncNullGram * NULLGRAM = new EncNullGram();


AlignmentModel::AlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, bool DEBUG) {    
    this->DEBUG = DEBUG;
    this->sourcemodel = sourcemodel;
    this->targetmodel = targetmodel;     
}


void AlignmentModel::intersect(AlignmentModel * reversemodel, double probthreshold, int bestn) {
	 //Compute intersection with reverse model
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		if (sourceiter->first == NULLGRAM) continue;
		const EncAnyGram * sourcegram = sourceiter->first;
		
		double lowerbound = 0.0;
        list<double> bestq;
		
		for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			if ((reversemodel->alignmatrix.count(targetgram) > 0) && (reversemodel->alignmatrix[targetgram].count(sourcegram) > 0) && (reversemodel->alignmatrix[targetgram][sourcegram] * targetiter->second >= probthreshold)) {
				const double p = reversemodel->alignmatrix[targetgram][sourcegram] * targetiter->second;
				alignmatrix[sourcegram][targetgram] = p;
				if ((bestn) && ((p > lowerbound) || (bestq.size() < bestn))) {
		    		orderedinsert(bestq, p);
		    		while (bestq.size() > bestn) bestq.pop_front();
		    		lowerbound = *bestq.begin();		    			
		    	}
			} else {
				alignmatrix[sourcegram].erase(targetgram); //no intersection, prune
			} 			
		}		
		if (bestn) {
			//pruning stage
			for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
				const EncAnyGram * targetgram = targetiter->first;
				if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
			}
		}
		if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);	
	}	 
}	 


const EncAnyGram * AlignmentModel::getsourcekey(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<EncNGram,bool>::iterator iter = sourcengrams.find(*( (EncNGram*) key) );
        if (iter != sourcengrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<EncSkipGram,bool>::iterator iter = sourceskipgrams.find(*( (EncSkipGram*) key) );
        if (iter != sourceskipgrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }        
    }
}


const EncAnyGram * AlignmentModel::gettargetkey(const EncAnyGram* key) {
    if (key->gapcount() == 0) {
        std::unordered_map<EncNGram,bool>::iterator iter = targetngrams.find(*( (EncNGram*) key) );
        if (iter != targetngrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }
    } else {
        std::unordered_map<EncSkipGram,bool>::iterator iter = targetskipgrams.find(*( (EncSkipGram*) key) );
        if (iter != targetskipgrams.end()) {
            return &iter->first;
        } else {
            return NULL;
        }        
    }
}


AlignmentModel::AlignmentModel(const string & filename, const int bestn) {
	DEBUG = false;
	unsigned char check;
	
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    uint64_t model_id;    
    uint64_t sourcecount = 0;    
    f.read( (char*) &model_id, sizeof(uint64_t));        
    f.read( (char*) &sourcecount, sizeof(uint64_t));        
     
    char gapcount;    
    for (int i = 0; i < sourcecount; i++) {	    
	    if (DEBUG) cerr << "\t@" << i << endl;
        f.read((char*) &check, sizeof(char));
        if (check != 0xff) {
        	cerr << "ERROR processing " + filename + " at construction " << i << " of " << sourcecount << ". Expected check-byte, got " << (int) check << endl;
        	f.read(&gapcount, sizeof(char));
        	cerr << "DEBUG: next byte should be gapcount, value=" << (int) gapcount << endl; 
        	exit(13);        	
        }
        f.read(&gapcount, sizeof(char));	 
        const EncAnyGram * sourcegram;   
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            EncNGram ngram = EncNGram(&f); //read from file            
            if (!getsourcekey((EncAnyGram*) &ngram)) {
            	sourcengrams[ngram] = true;            	
            }   
            sourcegram = getsourcekey((EncAnyGram*) &ngram);                                           
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            if (!getsourcekey((EncAnyGram*) &skipgram)) {
            	sourceskipgrams[skipgram] = true;            	
            }   
            sourcegram = getsourcekey((EncAnyGram*) &skipgram);                     
        }        
        uint64_t targetcount;
        f.read( (char*) &targetcount, sizeof(uint64_t));
        const EncAnyGram * besttargetgram = NULL;
        double lowerbound = 0.0;
        list<double> bestq;
        for (int j = 0; j < targetcount; j++) {
        	const EncAnyGram * targetgram = NULL;   
            f.read(&gapcount, sizeof(char));	    
		    if (gapcount == 0) {
		        if (DEBUG)  cerr << "\tNGRAM";
		        EncNGram ngram = EncNGram(&f); //read from file
		        if (!gettargetkey((EncAnyGram*) &ngram)) {
		        	targetngrams[ngram] = true;		        	
		        }   
		        targetgram = gettargetkey((EncAnyGram*) &ngram);                                           
		    } else {
		        if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
		        EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              		        
		        if (!gettargetkey((EncAnyGram*) &skipgram)) {
		        	targetskipgrams[skipgram] = true;		        	
		        }   
		        targetgram = gettargetkey((EncAnyGram*) &skipgram);                      
		    }		    
		    double p;
		    f.read((char*) &p, sizeof(double));
		    if (sourcegram != NULL and targetgram != NULL) {
		    	alignmatrix[sourcegram][targetgram] = p;
		    	if ((bestn) && ((p > lowerbound) || (bestq.size() < bestn))) {
		    		orderedinsert(bestq, p);
		    		while (bestq.size() > bestn) bestq.pop_front();
		    		lowerbound = *bestq.begin();		    			
		    	}
		    } else {
		    	cerr << "SOURCEGRAM or TARGETGRAM is NULL";
		    	exit(6);
		    }
        }
		if (bestn) {
			//pruning stage
			for (unordered_map<const EncAnyGram *, double>::iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
				if (iter->second < lowerbound) {
					alignmatrix[sourcegram].erase(iter->first);
				}				
			}		
		}            
	}
    f.close();
}

double AlignmentModel::cooc( CoocMode mode, const multiset<uint32_t> & sourceindex, const multiset<uint32_t> & targetindex, const double heurthreshold) {    
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



unsigned int AlignmentModel::trainCooc(CoocMode mode, const EncAnyGram * sourcegram, const multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel, const int bestn , const double absthreshold,  const double probthreshold, bool normalize) {           
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
		if (targetmodel->reverseindex.count(sentencenumber) > 0) {
			if (DEBUG) cerr << "\t\t\tReverseindex for sentence " << sentencenumber << " yields " << targetmodel->reverseindex[sentencenumber].size() << " target-side patterns" << endl;
			for (vector<const EncAnyGram*>::const_iterator reviter = targetmodel->reverseindex[sentencenumber].begin(); reviter != targetmodel->reverseindex[sentencenumber].end(); reviter++) {
				const EncAnyGram* targetgram = *reviter;
				targetpatterns.insert(targetgram);
			}
		}
		prevsentencenumber = sentencenumber;
    }
	if (DEBUG) cerr << "\t\tGathered " << targetpatterns.size() << " target-side patterns for given source pattern, computing co-occurence..." << endl;
	
	double lowerbound = 0.0;
	list<double> bestq;
	for (unordered_set<const EncAnyGram *>::const_iterator iter = targetpatterns.begin(); iter != targetpatterns.end(); iter++) {
			const EncAnyGram* targetgram = *iter;
	        multiset<uint32_t> * targetindex;
		    if (targetgram->gapcount() == 0) {
		       targetindex = &targetmodel->ngrams[*( (EncNGram*) targetgram)].sentences;
		    } else {
		       targetindex = &targetmodel->skipgrams[*( (EncSkipGram*) targetgram)].sentences;
		    }				    
		    const double coocvalue = cooc(mode, sourceindex, *targetindex, absthreshold);                    
		    if (coocvalue >= absthreshold) {
		    	//prune based on absolute co-occurrence value
		    	found++;
		        alignmatrix[sourcegram][targetgram] = coocvalue;
		    } else {
		    	prunedabs++;
		    }
		    totalcooc += coocvalue;	
		    if ((bestn) && ((coocvalue > lowerbound) || (bestq.size() < bestn)))  {
		    	orderedinsert(bestq, coocvalue);		    	
				while (bestq.size() > bestn) bestq.pop_front();
				lowerbound = *bestq.begin();
		    } 		    
	}
	if ((DEBUG) && (bestn))  cerr << "\t\tbest-n lowerbound=" << lowerbound << endl;				
    if ((totalcooc > 0) && (normalize || bestn || probthreshold > 0)) {
    	//normalisation and pruning step (based on probability threshold)
    	for (std::unordered_map<const EncAnyGram*, double>::const_iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
    		const double alignprob = (double) iter->second / totalcooc;
    		if ((alignprob < probthreshold) || ((bestn) && (iter->second < lowerbound)))  {
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
    if (DEBUG) cerr << "\t\t" << (found - prunedprob) << " alignments found (after pruning " << prunedabs << " on co-occurence value and " << prunedprob << " on alignment probability, including best-n)" << endl;
    return found - prunedprob;
}


unsigned int AlignmentModel::trainCooc(CoocMode mode,const int bestn ,  const double absthreshold, const double probthreshold) {
    unsigned int c = 0;
    unsigned int found = 0;
    for (unordered_map<EncNGram,IndexCountData >::const_iterator iter = sourcemodel->ngrams.begin();  iter != sourcemodel->ngrams.end(); iter++) {
    	c++;
        if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (ngram) -- " << found << " alignment possibilities thus-far" << endl;
        found += trainCooc(mode, &iter->first, iter->second.sentences, targetmodel, bestn, absthreshold, probthreshold);
    }    
    for (unordered_map<EncSkipGram,IndexCountData >::const_iterator iter = sourcemodel->skipgrams.begin();  iter != sourcemodel->skipgrams.end(); iter++) {
    	c++;
    	if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (skipgram) -- " << found << " alignment possibilities thus-far" << endl;
        found += trainCooc(mode, &iter->first, iter->second.sentences, targetmodel, bestn, absthreshold, probthreshold);
    }  
    return found;          
}


void AlignmentModel::trainEM(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn, const bool DONULL, const bool INIT) {
	// Compute p(target|source)      alignmatrix[source][target]
	/* 
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all patterns t in t_s
            for all patterns s in s_s
              total_s(t) += t(t|s)
         for all patterns t in t_s
             for all patterns s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)
	*/



    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    unsigned long cells = 0;
    
    if (INIT) {
		cerr << "  Initialisation step" << endl;
		double v = (double) 1 / targetmodel->types();     
		for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
			uint32_t sentence = reviter_source->first;
			const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
			if (targetmodel->reverseindex.count(sentence) > 0) {
				vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
				if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
				for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
					const EncAnyGram * targetgram = *targetiter;
					if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram] = v;
						cells++;
					}
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
        
        
        //EXPECTATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences    		
        
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal; 
        			for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			const EncAnyGram * sourcegram = *sourceiter;
		    			if (alignmatrix.count(sourcegram)) {
							for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {        				  
									const EncAnyGram * targetgram = *targetiter;
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
		    				const double countvalue_null = alignmatrix[NULLGRAM][targetgram] / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}
						
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {		                	
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		    			    if ((alignmatrix.count(sourcegram) && alignmatrix[sourcegram].count(targetgram))) {
		                    	const double countvalue = alignmatrix[sourcegram][targetgram] / sentencetotal[targetgram];
		                    	count[sourcegram][targetgram] += countvalue;
		                    	total[sourcegram] += countvalue;
		                    }
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
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}

void AlignmentModel::normalize() {
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;
		//compute sum
		double sum = 0;
		for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			sum += targetiter->second;			
		}
		//normalize
		for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			alignmatrix[sourcegram][targetgram] = alignmatrix[sourcegram][targetgram] / sum; 
		}
	}
}

void AlignmentModel::save(const string & filename) {
	//TODO: Merge with CoocAlignMentModel::save()

	const unsigned char check = 0xff;
	const char czero = 0;
		
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    uint64_t _id = 101;
    f.write( (char*) &_id, sizeof(uint64_t));
            
    uint64_t sourcecount = alignmatrix.size();
    if (alignmatrix.count(NULLGRAM) > 0) sourcecount--;
    f.write( (char*) &sourcecount, sizeof(uint64_t));         

    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;

    	
		f.write((char*) &check, sizeof(char));
			  
    	const EncAnyGram * sourcegram = iter->first;
    	if (sourcegram->isskipgram()) {
    		const EncSkipGram * skipgram = (const EncSkipGram*) sourcemodel->getkey(sourcegram);
    		skipgram->writeasbinary(&f);
    	} else {
    	    const EncNGram * ngram = (const EncNGram*) sourcemodel->getkey(sourcegram);
    	    f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
    		ngram->writeasbinary(&f);    		
    	}                
    	uint64_t targetcount = iter->second.size();
    	f.write( (char*) &targetcount, sizeof(uint64_t));
    	            
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	const EncAnyGram* targetgram = iter2->first;
        	const double p = iter2->second;
        	if (targetgram->isskipgram()) {
    			const EncSkipGram * skipgram = (const EncSkipGram*) targetmodel->getkey(targetgram);
    			if (skipgram == NULL) {
    				cerr << "TARGET-SIDE SKIPGRAM NOT FOUND!\n";
    				exit(3);
    			}  else {
	    			skipgram->writeasbinary(&f);
	    		}
    		} else {
    	    	const EncNGram * ngram = (const EncNGram*) targetmodel->getkey(targetgram);     	
	    		if (ngram == NULL) {
    				cerr << "TARGET-SIDE NGRAM NOT FOUND!";
    				exit(3);
    			}  else {
    				f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
	    			ngram->writeasbinary(&f);
	    		}
    		}                
        	f.write( (char*) &p, sizeof(double));
        }

    }    
    f.close();	
}



/*********************************/
/*
unsigned int CoocAlignmentModel::compute(const EncAnyGram * sourcegram, const multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel) {        
   
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
		if (targetmodel->reverseindex.count(sentencenumber) > 0) {
			if (DEBUG) cerr << "\t\t\tReverseindex for sentence " << sentencenumber << " yields " << targetmodel->reverseindex[sentencenumber].size() << " target-side patterns" << endl;
			for (vector<const EncAnyGram*>::const_iterator reviter = targetmodel->reverseindex[sentencenumber].begin(); reviter != targetmodel->reverseindex[sentencenumber].end(); reviter++) {
				const EncAnyGram* targetgram = *reviter;
				targetpatterns.insert(targetgram);
			}
		}
		prevsentencenumber = sentencenumber;
    }
	if (DEBUG) cerr << "\t\tGathered " << targetpatterns.size() << " target-side patterns for given source pattern, computing co-occurence..." << endl;
	
	double lowerbound = 0.0;
	list<double> bestq;
	for (unordered_set<const EncAnyGram *>::const_iterator iter = targetpatterns.begin(); iter != targetpatterns.end(); iter++) {
			const EncAnyGram* targetgram = *iter;
	        multiset<uint32_t> * targetindex;
		    if (targetgram->gapcount() == 0) {
		       targetindex = &targetmodel->ngrams[*( (EncNGram*) targetgram)].sentences;
		    } else {
		       targetindex = &targetmodel->skipgrams[*( (EncSkipGram*) targetgram)].sentences;
		    }				    
		    const double coocvalue = cooc(mode, sourceindex, *targetindex, absthreshold);                    
		    if (coocvalue >= absthreshold) {
		    	//prune based on absolute co-occurrence value
		    	found++;
		        alignmatrix[sourcegram][targetgram] = coocvalue;
		    } else {
		    	prunedabs++;
		    }
		    totalcooc += coocvalue;	
		    if ((bestn) && ((coocvalue > lowerbound) || (bestq.size() < bestn)))  {
		    	orderedinsert(bestq, coocvalue);		    	
				while (bestq.size() > bestn) bestq.pop_front();
				lowerbound = *bestq.begin();
		    } 		    
	}
	if ((DEBUG) && (bestn))  cerr << "\t\tbest-n lowerbound=" << lowerbound << endl;				
    if ((totalcooc > 0) && (normalize || bestn || probthreshold > 0)) {
    	//normalisation and pruning step (based on probability threshold)
    	for (std::unordered_map<const EncAnyGram*, double>::const_iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
    		const double alignprob = (double) iter->second / totalcooc;
    		if ((alignprob < probthreshold) || ((bestn) && (iter->second < lowerbound)))  {
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
    if (DEBUG) cerr << "\t\t" << (found - prunedprob) << " alignments found (after pruning " << prunedabs << " on co-occurence value and " << prunedprob << " on alignment probability, including best-n)" << endl;
    return found - prunedprob;
}

CoocAlignmentModel::CoocAlignmentModel(CoocMode mode,SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int bestn ,  const double absthreshold, const double probthreshold, bool normalize, bool DEBUG) {
    this->mode = mode;
    this->absthreshold = absthreshold;
    this->probthreshold = probthreshold;
    this->normalize = normalize;
    this->bestn = bestn;
    this->DEBUG = DEBUG;
    this->sourcemodel = sourcemodel;
    this->targetmodel = targetmodel;
    unsigned int c = 0;
    unsigned int found = 0;
    for (unordered_map<EncNGram,IndexCountData >::const_iterator iter = sourcemodel->ngrams.begin();  iter != sourcemodel->ngrams.end(); iter++) {
    	c++;
        if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (ngram) -- " << found << " alignment possibilities thus-far" << endl;
        found += compute(&iter->first, iter->second.sentences, targetmodel);
    }    
    for (unordered_map<EncSkipGram,IndexCountData >::const_iterator iter = sourcemodel->skipgrams.begin();  iter != sourcemodel->skipgrams.end(); iter++) {
    	c++;
    	if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (skipgram) -- " << found << " alignment possibilities thus-far" << endl;
        found += compute(&iter->first, iter->second.sentences, targetmodel);
    }            
}
    */
void AlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT) {
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;
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

void BiAlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT) {
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;
        const EncAnyGram* sourcegram = iter->first;
        *OUT << sourcegram->decode(sourceclassdecoder) << "\t";
        map<double, const EncAnyGram*> sorted;        
        for (unordered_map<const EncAnyGram*, double >::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	sorted[-1 * iter2->second] = iter2->first;            
        }
        for (map<double, const EncAnyGram*>::iterator iter2 = sorted.begin(); iter2 != sorted.end(); iter2++) {
			const EncAnyGram* targetgram = iter2->second;            
            *OUT << targetgram->decode(targetclassdecoder) << "\t" << alignmatrix[sourcegram][targetgram] << "\t" << alignmatrixrev[targetgram][sourcegram] << "\t";
        }            
        *OUT << endl;
    }
}


void AlignmentModel::simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT, bool targetfirst,  bool wordbased, bool mosesformat) {
	/* output a simple word-based lexicon, similar to the one used in moses (s2t, t2s) */
	string delimiter;
	if (mosesformat) {
		delimiter = " ||| ";
	} else {
		delimiter = " ";
	}
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;    
        const EncAnyGram* sourcegram = iter->first;
        if (wordbased && sourcegram->n() > 1) continue;         
        map<double, const EncAnyGram*> sorted;        
        double total = 0;
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	if (wordbased && iter2->first->n() > 1) continue; 
        	sorted[-1 * iter2->second] = iter2->first;
        	total += iter2->second;            
        }
        for (map<double, const EncAnyGram*>::iterator iter2 = sorted.begin(); iter2 != sorted.end(); iter2++) {
			const EncAnyGram* targetgram = iter2->second;
			if (targetfirst) {
				*OUT << targetgram->decode(targetclassdecoder) << delimiter << sourcegram->decode(sourceclassdecoder);				 
			} else {			
				*OUT << sourcegram->decode(sourceclassdecoder) << delimiter << targetgram->decode(targetclassdecoder);
			}
			*OUT << delimiter << ((-1 * iter2->first) / total) << endl;
        }            
    }	
}


void BiAlignmentModel::simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT, bool targetfirst,  bool wordbased, bool mosesformat) {
	/* output a simple word-based lexicon, similar to the one used in moses (s2t, t2s) */
	string delimiter;
	if (mosesformat) {
		delimiter = " ||| ";
	} else {
		delimiter = " ";
	}
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double > >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;    
        const EncAnyGram* sourcegram = iter->first;
        if (wordbased && sourcegram->n() > 1) continue;         
        map<double, const EncAnyGram*> sorted;        
        double sourcetotal = 0;        
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	if (wordbased && iter2->first->n() > 1) continue; 
        	sorted[-1 * iter2->second] = iter2->first;
        	sourcetotal += iter2->second;            
        }
        for (map<double, const EncAnyGram*>::iterator iter2 = sorted.begin(); iter2 != sorted.end(); iter2++) {
			const EncAnyGram* targetgram = iter2->second;
			double targettotal = 0;
			if (alignmatrixrev.count(targetgram) > 0) { 
				for (unordered_map<const EncAnyGram*, double>::iterator iter3 = alignmatrixrev[targetgram].begin(); iter3 != alignmatrixrev[targetgram].end(); iter3++) {
					if (wordbased && iter3->first->n() > 1) continue; 
					targettotal += iter3->second;            
				}			
				if (targetfirst) {
					*OUT << targetgram->decode(targetclassdecoder) << delimiter << sourcegram->decode(sourceclassdecoder);				 
				} else {			
					*OUT << sourcegram->decode(sourceclassdecoder) << delimiter << targetgram->decode(targetclassdecoder);
				}
				*OUT << delimiter << ((alignmatrix[sourcegram][targetgram]) / sourcetotal) << ' ' << ((alignmatrixrev[targetgram][sourcegram]) / targettotal) << endl;
			}
        }            
    }	
}

/***************************** BEGIN EM ************************************/
/*
EMAlignmentModel::EMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, bool INIT, bool DONULL, bool DEBUG) {
	this->sourcemodel = sourcemodel;
	this->targetmodel = targetmodel;
	this->INIT = INIT;
	this->DEBUG = DEBUG;
	this->DONULL = DONULL;
	//initialise uniformly
    if (INIT) {
		cerr << "  Initialisation step" << endl;
		double v = (double) 1 / targetmodel->types();     
		for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
			uint32_t sentence = reviter_source->first;
			const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
			if (targetmodel->reverseindex.count(sentence) > 0) {
				vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
				if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
				for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
					const EncAnyGram * targetgram = *targetiter;
					if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram] = v;
					}
				}
			}
		}
	}    
}


void EMAlignmentModel::train(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn) {
	// Compute p(target|source)      alignmatrix[source][target]
	 
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all words t in t_s
            for all words s in s_s
              total_s(t) += t(t|s)
         for all words t in t_s
             for all words s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)
	



    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
      
            
            
    do {       
        round++; 
        c = 0;        
        cerr << "  EM Round " << round << "... ";
        
		std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
		unordered_map<const EncAnyGram*, double> total;
        
        
        //EXPECTATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences    		
        
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal; 
        			for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			const EncAnyGram * sourcegram = *sourceiter;
		    			if (alignmatrix.count(sourcegram)) {
							for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {        				  
									const EncAnyGram * targetgram = *targetiter;
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
		    				const double countvalue_null = alignmatrix[NULLGRAM][targetgram] / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}
						
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {		                	
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		    			    if ((alignmatrix.count(sourcegram) && alignmatrix[sourcegram].count(targetgram))) {
		                    	const double countvalue = alignmatrix[sourcegram][targetgram] / sentencetotal[targetgram];
		                    	count[sourcegram][targetgram] += countvalue;
		                    	total[sourcegram] += countvalue;
		                    }
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
                totaldivergence += divergence;
                c++;        		
        	}
        }
        
        
		
        const double avdivergence = (double) totaldivergence / c;
        converged = ((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence) <= CONVERGEDTHRESHOLD);               
        cerr << "   average divergence = " << avdivergence << ", delta with prev divergence = " << abs(avdivergence - prevavdivergence) << " > " << CONVERGEDTHRESHOLD << ", alignprob size = " << alignmatrix.size() << endl;
        prevavdivergence = avdivergence;
    } while (!converged);    
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}
*/




/************************************************* END EM ***************/
/*
ItEMAlignmentModel::ItEMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn, bool DONULL, bool DEBUG) {
	// Compute p(target|source)      alignmatrix[source][target]
	
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all words t in t_s
            for all words s in s_s
              total_s(t) += t(t|s)
         for all words t in t_s
             for all words s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)


	this->sourcemodel = sourcemodel;
	this->targetmodel = targetmodel;

    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
    int maxn = 0;
    
    //initialise uniformly
    cerr << "  Initialisation step" << endl;
    double v = (double) 1 / targetmodel->types();     
    for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
    	uint32_t sentence = reviter_source->first;
		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;		
		if (targetmodel->reverseindex.count(sentence) > 0) {
			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
			for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
    			const EncAnyGram * targetgram = *targetiter;
    			if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
    			    const EncAnyGram * sourcegram = *sourceiter;
					alignmatrix[sourcegram][targetgram] = v;
					if (sourcegram->n() > maxn) maxn = sourcegram->n();
    			}
    		}
		}
    }
          
            
            
    do {       
        round++; 
        c = 0;        
        cerr << "  EM Round " << round << "... " << endl;
        
        std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
		unordered_map<const EncAnyGram*, double> total;
			
        for (int n = 1; n <= maxn; n++) {
        	cerr << "      n=" << n << endl;
		    
		    //EXPECTATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
		    for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences    		
		    		uint32_t sentence = reviter_source->first;
		    		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
		    		if (targetmodel->reverseindex.count(sentence) > 0) {
		    			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
		    			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
		    			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
		    			unordered_map<const EncAnyGram*, double> sentencetotal;
		    			for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    				const EncAnyGram * sourcegram = *sourceiter;      
		    				if (sourcegram->n() == n) {
								for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
									const EncAnyGram * targetgram = *targetiter;
									if (DONULL) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram];
									sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration 
								} 
		    				}
		    			}
		    			
							
				        //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
				        //for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {
				        
				        for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sentencetotal.begin(); targetiter != sentencetotal.end(); targetiter++) {  
							const EncAnyGram * targetgram = targetiter->first;
							
							//the null condition:
							if (DONULL) {
								const double countvalue_null = alignmatrix[NULLGRAM][targetgram] / targetiter->second;
				            	count[NULLGRAM][targetgram] += countvalue_null;
								total[NULLGRAM] += countvalue_null;
							}
						
				            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
							    const EncAnyGram * sourcegram = *sourceiter;
							    if (sourcegram->n() == n) {                                                                 
				                	const double countvalue = alignmatrix[sourcegram][targetgram] / targetiter->second;
				                	count[sourcegram][targetgram] += countvalue;
				                	total[sourcegram] += countvalue;
				                }
				            }
				        }
				   
		    		}	
			} //end loop over corpus
		
		}
		
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
                totaldivergence += divergence;
                c++;        		
        	}
        }
        
        
		
        const double avdivergence = (double) totaldivergence / c;
        converged = ((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence) <= CONVERGEDTHRESHOLD);               
        cerr << "   average divergence = " << avdivergence << ", delta with prev divergence = " << abs(avdivergence - prevavdivergence) << " > " << CONVERGEDTHRESHOLD << ", alignprob size = " << alignmatrix.size() << endl;
        prevavdivergence = avdivergence;
    } while (!converged);    
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}

EMAlignmentModel3::EMAlignmentModel3(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn, bool DONULL, bool DEBUG) {
	// Compute p(target|source)      alignmatrix[source][target]
	
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all words t in t_s
            for all words s in s_s
              total_s(t) += t(t|s)
         for all words t in t_s
             for all words s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)
	

	this->sourcemodel = sourcemodel;
	this->targetmodel = targetmodel;

    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
    
    
    
    //initialise uniformly
    cerr << "  Initialisation step" << endl;
    double v = (double) 1 / targetmodel->types();     
    for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
    	uint32_t sentence = reviter_source->first;
		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
		if (targetmodel->reverseindex.count(sentence) > 0) {
			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
			for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
    			const EncAnyGram * targetgram = *targetiter;
    			if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
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
		unsigned int found = 0;

		//EXPECTATION STEP: collect counts to estimate improved model
		for (unordered_map<EncNGram,IndexCountData >::const_iterator iter = sourcemodel->ngrams.begin();  iter != sourcemodel->ngrams.end(); iter++) {
			c++;
		    if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (ngram) -- " << found << " alignment possibilities thus-far" << endl;
		    found += expectation(&iter->first, iter->second.sentences, targetmodel, count, total);
		}    
		for (unordered_map<EncSkipGram,IndexCountData >::const_iterator iter = sourcemodel->skipgrams.begin();  iter != sourcemodel->skipgrams.end(); iter++) {
			c++;
			if ((c % 1000 == 0) || (DEBUG)) cerr << "\t@" << c << " (skipgram) -- " << found << " alignment possibilities thus-far" << endl;
		    found += expectation(&iter->first, iter->second.sentences, targetmodel, count, total);
		}            

        		
		
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
                totaldivergence += divergence;
                c++;        		
        	}
        }
        
        
		
        const double avdivergence = (double) totaldivergence / c;
        converged = ((round >= MAXROUNDS) || abs(avdivergence - prevavdivergence) <= CONVERGEDTHRESHOLD);               
        cerr << "   average divergence = " << avdivergence << ", delta with prev divergence = " << abs(avdivergence - prevavdivergence) << " > " << CONVERGEDTHRESHOLD << ", alignprob size = " << alignmatrix.size() << endl;
        prevavdivergence = avdivergence;
    } while (!converged);    
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}


unsigned int EMAlignmentModel3::expectation(const EncAnyGram * sourcegram, const multiset<uint32_t> & sourceindex, SelectivePatternModel * targetmodel, std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > & count, unordered_map<const EncAnyGram*, double> & total) {        
   
    uint32_t prevsentencenumber = 0;
    
    unordered_map<const EncAnyGram*, double> sourcegramtotal;   
	unordered_set<const EncAnyGram *> targetpatterns;	
    //cerr << "Processing new construction" << endl;
    
    //compute sourcegramtotal for normalisation later in count step, sum_s(p(t|s))
    if (DEBUG) cerr << "\t\tForward index yields " << sourceindex.size() << " sentence references" << endl;
    for (multiset<uint32_t>::const_iterator iter = sourceindex.begin(); iter != sourceindex.end(); iter++) {
        const uint32_t sentencenumber = *iter;        
        if (sentencenumber == prevsentencenumber) continue;
		if (targetmodel->reverseindex.count(sentencenumber) > 0) {
			if (DEBUG) cerr << "\t\t\tReverseindex for sentence " << sentencenumber << " yields " << targetmodel->reverseindex[sentencenumber].size() << " target-side patterns" << endl;
			for (vector<const EncAnyGram*>::const_iterator reviter = targetmodel->reverseindex[sentencenumber].begin(); reviter != targetmodel->reverseindex[sentencenumber].end(); reviter++) {
				const EncAnyGram* targetgram = *reviter;
				targetpatterns.insert(targetgram);
				sourcegramtotal[targetgram] += alignmatrix[sourcegram][targetgram];
			}
		}
		prevsentencenumber = sentencenumber;
    }
	if (DEBUG) cerr << "\t\tGathered " << targetpatterns.size() << " target-side patterns for given source pattern, computing expectation..." << endl;


	unsigned int found = 0;
	
    //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
    for (unordered_set<const EncAnyGram *>::const_iterator targetiter = targetpatterns.begin(); targetiter != targetpatterns.end(); targetiter++) {  
		const EncAnyGram * targetgram = *targetiter;		

        const double countvalue = alignmatrix[sourcegram][targetgram] / sourcegramtotal[targetgram];
        count[sourcegram][targetgram] += countvalue;
        total[sourcegram] += countvalue;        
        found++;
    }
	return found;
}

void CoocAlignmentModel::save(const string & filename) {
	const unsigned char check = 0xff;
	const char czero = 0;
		
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    uint64_t _id = 101;
    f.write( (char*) &_id, sizeof(uint64_t));
            
    uint64_t sourcecount = alignmatrix.size();
    f.write( (char*) &sourcecount, sizeof(uint64_t));         

    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
		f.write((char*) &check, sizeof(char));
			  
    	const EncAnyGram * sourcegram = iter->first;
    	if (sourcegram->isskipgram()) {
    		const EncSkipGram * skipgram = (const EncSkipGram*) sourcemodel->getkey(sourcegram);
    		skipgram->writeasbinary(&f);
    	} else {
    	    const EncNGram * ngram = (const EncNGram*) sourcemodel->getkey(sourcegram);
    	    f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
    		ngram->writeasbinary(&f);    		
    	}                
    	uint64_t targetcount = iter->second.size();
    	f.write( (char*) &targetcount, sizeof(uint64_t));
    	            
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	const EncAnyGram* targetgram = iter2->first;
        	const double p = iter2->second;
        	if (targetgram->isskipgram()) {
    			const EncSkipGram * skipgram = (const EncSkipGram*) targetmodel->getkey(targetgram);
    			if (skipgram == NULL) {
    				cerr << "TARGET-SIDE SKIPGRAM NOT FOUND!\n";
    				exit(3);
    			}  else {
	    			skipgram->writeasbinary(&f);
	    		}
    		} else {
    	    	const EncNGram * ngram = (const EncNGram*) targetmodel->getkey(targetgram);     	
	    		if (ngram == NULL) {
    				cerr << "TARGET-SIDE NGRAM NOT FOUND!";
    				exit(3);
    			}  else {
    				f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
	    			ngram->writeasbinary(&f);
	    		}
    		}                
        	f.write( (char*) &p, sizeof(double));
        }

    }    
    f.close();	
}
*/







BiAlignmentModel::BiAlignmentModel(const string & s2tfilename, const string & t2sfilename): AlignmentModel(s2tfilename) {
	//TODO: Code duplication, merge with AlignmentModel() constructor
	DEBUG = false;
	unsigned char check;
	
    ifstream f;
    f.open(t2sfilename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << t2sfilename << endl;
       exit(3);
    }
    
    uint64_t model_id;    
    uint64_t targetcount = 0;    
    f.read( (char*) &model_id, sizeof(uint64_t));        
    f.read( (char*) &targetcount, sizeof(uint64_t));        
     
    char gapcount;    
    for (int i = 0; i < targetcount; i++) {	    
	    if (DEBUG) cerr << "\t@" << i << endl;
        f.read((char*) &check, sizeof(char));
        if (check != 0xff) {
        	cerr << "ERROR processing " << t2sfilename << " at construction " << i << " of " << targetcount << ". Expected check-byte, got " << (int) check << endl;
        	f.read(&gapcount, sizeof(char));
        	cerr << "DEBUG: next byte should be gapcount, value=" << (int) gapcount << endl; 
        	exit(13);        	
        }
        f.read(&gapcount, sizeof(char));	 
        const EncAnyGram * targetgram;   
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            EncNGram ngram = EncNGram(&f); //read from file            
            if (!gettargetkey((EncAnyGram*) &ngram)) {
            	targetngrams[ngram] = true;            	
            }   
            targetgram = gettargetkey((EncAnyGram*) &ngram);                                           
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              
            if (!gettargetkey((EncAnyGram*) &skipgram)) {
            	targetskipgrams[skipgram] = true;            	
            }   
            targetgram = gettargetkey((EncAnyGram*) &skipgram);                     
        }        
        uint64_t sourcecount;
        f.read( (char*) &sourcecount, sizeof(uint64_t));
        for (int j = 0; j < sourcecount; j++) {
        	const EncAnyGram * sourcegram = NULL;   
            f.read(&gapcount, sizeof(char));	    
		    if (gapcount == 0) {
		        if (DEBUG)  cerr << "\tNGRAM";
		        EncNGram ngram = EncNGram(&f); //read from file
		        if (!getsourcekey((EncAnyGram*) &ngram)) {
		        	sourcengrams[ngram] = true;		        	
		        }   
		        sourcegram = getsourcekey((EncAnyGram*) &ngram);                                           
		    } else {
		        if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
		        EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file              		        
		        if (!getsourcekey((EncAnyGram*) &skipgram)) {
		        	sourceskipgrams[skipgram] = true;		        	
		        }   
		        sourcegram = getsourcekey((EncAnyGram*) &skipgram);                      
		    }
		    double p;
		    f.read((char*) &p, sizeof(double));
		    if (sourcegram != NULL and targetgram != NULL) {
		    	alignmatrixrev[targetgram][sourcegram] = p;
		    } else {
		    	cerr << "SOURCEGRAM or TARGETGRAM is NULL";
		    	exit(6);
		    }
        }        
	}
    f.close();
}


int AlignmentModel::graphalign(SelectivePatternModel & sourcemodel, SelectivePatternModel & targetmodel, double impactfactor) {
	int adjustments = 0;

	if (!sourcemodel.has_parents()) {
		cerr << "ERROR: No parent relations available for source model; required for graphalign" << endl;
		exit(6);
	}
	if (!targetmodel.has_parents()) {
		cerr << "ERROR: No parent relations available for target model; required for graphalign" << endl;
		exit(6);
	}

	std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > weightmatrix;  

	//for sourcegram in source
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;	
		//if sourcegram has parent relations;
		if (sourcemodel.rel_subsumption_parents.count(sourcegram) > 0) {
			//for sourceparentgram in parents_sourcegram:
			for (unordered_set<const EncAnyGram*>::const_iterator parentiterS = sourcemodel.rel_subsumption_parents[sourcegram].begin(); parentiterS != sourcemodel.rel_subsumption_parents[sourcegram].end(); parentiterS++) {
				const EncAnyGram * sourceparentgram = *parentiterS;
				//for t in targets_aligned_to_sourcegram:
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					//if targetgram has parent relations;
					if (targetmodel.rel_subsumption_parents.count(targetgram) > 0) {
						for (unordered_set<const EncAnyGram*>::const_iterator parentiterT = targetmodel.rel_subsumption_parents[targetgram].begin(); parentiterT != targetmodel.rel_subsumption_parents[targetgram].end(); parentiterT++) {
								const EncAnyGram * targetparentgram = *parentiterT;
								//check if parents are aligned too
								if (alignmatrix.count(sourceparentgram) > 0 && (alignmatrix[sourceparentgram].count(targetparentgram) > 0)) {
									if (DEBUG) cerr << "graphalign: parallel alignments found for parents, strenghtening ties (+1)" << endl;
									//yes!
									//strenghten relation between aligned parents
									weightmatrix[sourceparentgram][targetparentgram]++;
									if (DEBUG) cerr << "\t" << weightmatrix[sourceparentgram][targetparentgram] << endl;
									//strenghten relation between aligned pair
									weightmatrix[sourcegram][targetgram]++;									
									if (DEBUG) cerr << "\t" << weightmatrix[sourcegram][targetgram] << endl;
								} else {
									//check for crossed-alignments: sourceparent with targetchild
									if (alignmatrix.count(sourceparentgram) > 0 && (alignmatrix[sourceparentgram].count(targetgram) > 0)) {
										if (DEBUG) cerr << "graphalign: cross-alignment 1 found, weakening ties (-0.5)" << endl;
										weightmatrix[sourceparentgram][targetgram] -= 0.5;
										if (DEBUG) cerr << "\t" << weightmatrix[sourceparentgram][targetgram] << endl;
										weightmatrix[sourcegram][targetgram] -= 0.5;												
									}
									//check for crossed-alignments: source with targetparent
									if (alignmatrix.count(sourcegram) > 0 && (alignmatrix[sourcegram].count(targetparentgram) > 0)) {
										if (DEBUG) cerr << "graphalign: cross-alignment 2 found, weakening ties (-0.5)" << endl;
										weightmatrix[sourcegram][targetparentgram] -= 0.5;
										if (DEBUG) cerr << "\t" << weightmatrix[sourcegram][targetparentgram] << endl;
										weightmatrix[sourcegram][targetgram] -= 0.5;											
									}
								}
						}	
					}				
				}								
			}
		}	
	}
	
	//apply weights	
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = weightmatrix.begin(); iter != weightmatrix.end(); iter++) {
		const EncAnyGram * sourcegram = iter->first;	
		for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = weightmatrix[sourcegram].begin(); targetiter != weightmatrix[sourcegram].end(); targetiter++) {		
			const EncAnyGram * targetgram = targetiter->first;
			//convert the weights from to: -inf -- 0 -- inf  --> 0 -- 1 -- inf , so they can be applied directly
			const double weight = weightmatrix[sourcegram][targetgram];
			
			const double a = alignmatrix[sourcegram][targetgram];
			if (weight < 0) {					 
				alignmatrix[sourcegram][targetgram] = pow(a, weight+1);
			} else if (weight > 0) {
				alignmatrix[sourcegram][targetgram] = pow(a, -1 * weight + 1);
			}
			adjustments++;							
		}
	}
	
	return adjustments;	
}			


void orderedinsert(list<double> & l, double value) {
	size_t index = 0;
	for (list<double>::iterator iter = l.begin(); iter != l.end(); iter++) {
		if (value < *iter) {
			l.insert(iter, value);
			return;
		}
	} 
	l.push_back(value);
}



/**************************** EXPERIMENTAL EM MODEL **********************************/






EMAlignmentModel2::EMAlignmentModel2(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, bool INIT, bool DONULL, bool DEBUG) {
	this->sourcemodel = sourcemodel;
	this->targetmodel = targetmodel;
	this->INIT = INIT;
	this->DEBUG = DEBUG;
	this->DONULL = DONULL;
	//initialise uniformly
    if (INIT) {
		cerr << "  Initialisation step" << endl;
		double v = (double) 1 / targetmodel->types();     
		for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
			uint32_t sentence = reviter_source->first;
			const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
			if (targetmodel->reverseindex.count(sentence) > 0) {
				vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
				if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
				for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
					const EncAnyGram * targetgram = *targetiter;
					if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram] = v;
					}
				}
			}
		}
	}    
}

void EMAlignmentModel2::train(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn) {
	// Compute p(target|source)      alignmatrix[source][target]
	/* 
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all words t in t_s
            for all words s in s_s
              total_s(t) += t(t|s)
         for all words t in t_s
             for all words s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)
	*/



    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
      
            
            
    do {       
        round++; 
        c = 0;        
        cerr << "  EM Round " << round << "... ";
        
		std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> > count;                
		unordered_map<const EncAnyGram*, double> total;
        
        
        //EXPECTATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences    		
        
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal; 
        			for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			const EncAnyGram * sourcegram = *sourceiter;
		    			if (alignmatrix.count(sourcegram)) {
							for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {        				  
									const EncAnyGram * targetgram = *targetiter;
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
		    				const double countvalue_null = alignmatrix[NULLGRAM][targetgram] / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}
						
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {		                	
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		    			    if ((alignmatrix.count(sourcegram) && alignmatrix[sourcegram].count(targetgram))) {
		                    	const double countvalue = alignmatrix[sourcegram][targetgram] / sentencetotal[targetgram];
		                    	count[sourcegram][targetgram] += countvalue;
		                    	total[sourcegram] += countvalue;
		                    }
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
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}


void EMAlignmentModel2::save(const string & filename) {
	//TODO: Merge with CoocAlignMentModel::save()

	const unsigned char check = 0xff;
	const char czero = 0;
		
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    uint64_t _id = 101;
    f.write( (char*) &_id, sizeof(uint64_t));
            
    uint64_t sourcecount = alignmatrix.size();
    if (alignmatrix.count(NULLGRAM) > 0) sourcecount--;
    f.write( (char*) &sourcecount, sizeof(uint64_t));         

    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;

    	
		f.write((char*) &check, sizeof(char));
			  
    	const EncAnyGram * sourcegram = iter->first;
    	if (sourcegram->isskipgram()) {
    		const EncSkipGram * skipgram = (const EncSkipGram*) sourcemodel->getkey(sourcegram);
    		skipgram->writeasbinary(&f);
    	} else {
    	    const EncNGram * ngram = (const EncNGram*) sourcemodel->getkey(sourcegram);
    	    f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
    		ngram->writeasbinary(&f);    		
    	}                
    	uint64_t targetcount = iter->second.size();
    	f.write( (char*) &targetcount, sizeof(uint64_t));
    	            
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	const EncAnyGram* targetgram = iter2->first;
        	const double p = iter2->second;
        	if (targetgram->isskipgram()) {
    			const EncSkipGram * skipgram = (const EncSkipGram*) targetmodel->getkey(targetgram);
    			if (skipgram == NULL) {
    				cerr << "TARGET-SIDE SKIPGRAM NOT FOUND!\n";
    				exit(3);
    			}  else {
	    			skipgram->writeasbinary(&f);
	    		}
    		} else {
    	    	const EncNGram * ngram = (const EncNGram*) targetmodel->getkey(targetgram);     	
	    		if (ngram == NULL) {
    				cerr << "TARGET-SIDE NGRAM NOT FOUND!";
    				exit(3);
    			}  else {
    				f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
	    			ngram->writeasbinary(&f);
	    		}
    		}                
        	f.write( (char*) &p, sizeof(double));
        }

    }    
    f.close();	
}





void AlignmentModel::trainEM2(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn, const bool DONULL, const bool INIT) {
    /* VERSION that considers multiple occurrences of the same type in a phrase pair */
    

	// Compute p(target|source)      alignmatrix[source][target]
	/* 
	initialize t(t|s) uniformly
   do until convergence
   	  set count(t|s) to 0 for all t,s
  	  set total(s) to 0 for all s
      for all sentence pairs (t_s,s_s)
         set total_s(t) = 0 for all t
         for all words t in t_s
            for all words s in s_s
              total_s(t) += t(t|s)
         for all words t in t_s
             for all words s in s_s
                count(t|s) += t(t|s) / total_s(t)
                total(s)   += t(t|s) / total_s(t)
      for all s
     	for all t
           t(t|s) = count(t|s) / total(s)
	*/



    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
    if (INIT) {
		cerr << "  Initialisation step" << endl;
		double v = (double) 1 / targetmodel->types();     
		for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
			uint32_t sentence = reviter_source->first;
			const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
			if (targetmodel->reverseindex.count(sentence) > 0) {
				vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
				if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;			
				for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
					const EncAnyGram * targetgram = *targetiter;
					if (DONULL) alignmatrix[NULLGRAM][targetgram] = v;
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram] = v;
					}
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
        
        
        //EXPECTATION STEP: collect counts to estimate improved model -- use reverse index to iterate over all sentences in training data 
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences    		
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;

        		
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal; 
        			for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			const EncAnyGram * sourcegram = *sourceiter;		    			
		    			const int sourcemultioc = sourcemodel->countforsentence(sourcegram, sentence);
		    			if (alignmatrix.count(sourcegram)) {
							for (vector<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {        				  
									const EncAnyGram * targetgram = *targetiter;
									const int targetmultioc = targetmodel->countforsentence(targetgram, sentence);
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += sourcemultioc * targetmultioc * alignmatrix[sourcegram][targetgram]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			const int targetmultioc = targetmodel->countforsentence(targetgram, sentence);
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
		    				const double countvalue_null = (targetmultioc * alignmatrix[NULLGRAM][targetgram]) / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}
						
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {		                	
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		    			    if ((alignmatrix.count(sourcegram) && alignmatrix[sourcegram].count(targetgram))) {
		    			        const int sourcemultioc = sourcemodel->countforsentence(sourcegram, sentence);
		                    	const double countvalue = (sourcemultioc * targetmultioc * alignmatrix[sourcegram][targetgram]) / sentencetotal[targetgram];
		                    	count[sourcegram][targetgram] += countvalue;
		                    	total[sourcegram] += countvalue;
		                    }
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
    
    
    
    
    if ((probthreshold > 0) || (bestn > 0)) {
    
		cerr << "  Pruning stage... ";
		for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if ((targetiter->second > lowerbound) || (bestq.size() < bestn)) {
						orderedinsert(bestq, targetiter->second);
						while (bestq.size() > bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}

/*void AlignmentModel::growdiag(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s) {
    
}*/


int AlignmentModel::prune(const double prunethreshold) {
    int pruned = 0;
	for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;
	    for (unordered_map<const EncAnyGram*, double>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
		    const EncAnyGram * targetgram = targetiter->first;
		    if (targetiter->second < prunethreshold) {
		        pruned++;
		        alignmatrix[sourcegram].erase(targetgram);
		    }
		}
		if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
    }
    return pruned;
}

void AlignmentModel::extractgizapatterns(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold,ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
        GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
        GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);
        
        int found = 0;
        cerr << "@" << sentenceindex << endl;
        
        //extract patterns in sentence pair
     
        /*
                    
                            
           2)  
            patterns_t = all patterns in sentence_t           
            for word_s in sentence_s:
                patterns_s = find patterns BEGINNING WITH word_s
                for pattern_s in patterns_s:
                    bestscore = 0
                    for pattern_t in patterns_t:
                        score = 0
                        for patternword_s in pattern_s:
                            if patternword_s aligned with anything outside pattern_t (according to intersection)
                                score = -2; break  
                            elseif patternword_s aligned with something in pattern_t (according to intersection)
                                score += 2
                                mark index in pattern_t as aligned                                
                            elseif patternword_s aligned with anything outside pattern_t (according to union)
                                score = score - 1
                            elseif patternword_s not aligned:
                                pass 
                        for patternword_t in pattern_t:
                            if not marked as aligned: (in above loop)
                                   if patternword_t aligned with anything outside pattern_s (according to intersection)
                                        score = -2; break                
                                    elseif patternword_t aligned with anything outside pattern_s (according to union)
                                        score = score - 1
                                    elseif patternword_t not aligned:
                                        pass                                  
                       if score > 0
                            bestscore = score
                            bestpattern_t = pattern_t
                            
      3)  
            patterns_t = all patterns in sentence_t           
            for word_s in sentence_s:
                patterns_s = find patterns BEGINNING WITH word_s
                for pattern_s in patterns_s:
                    bestscore = 0
                    for pattern_t in patterns_t:
                        aligned = 0
                        halfaligned = 0
                        unaligned = 0
                        firstsourcealigned = false
                        lastsourcealigned = false
                        firsttargetaligned = false
                        lasttargetaligned = false
                        for alignedsourceindex, alignedtargetindex in intersection:
                            if alignedsourceindex not in pattern_s or alignedtargetindex not in pattern_t:
                                aligned--; break;
                            else:
                                aligned++;
                                if alignedsourceindex == sourceindex: firstsourcealigned = true
                                if alignedsourceindex == sourceindex + patternsize_s: lastsourcealigned = true
                                if alignedtargetindex == targetindex: firstsourcealigned = true
                                if alignedtargetindex == targetindex + patternsize_t: lastsourcealigned = true
                            else:
                                unaligned++;                                
                        if ((aligned < 0) || (!firstaligned) || (!lastaligned)) break;
                        for alignedsourceindex, alignedtargetindex in union:                            
                            if (alignedsourceindex in pattern_s and alignedtargetindex not in pattern_t) or (alignedsourceindex not in pattern_s and alignedtargetindex in pattern_t):
                                halfaligned++;
                            else if not (alignedsourceindex in pattern_s and alignedtargetindex not in pattern_t):
                                halfaligned--;
                                        
                       if score > 0
                            bestscore = score
                            bestpattern_t = pattern_t                            
                        
        
        */

        //get all patterns in the target sentence
        vector<const EncAnyGram*> * sourcepatterns = &sourcemodel->reverseindex[sentenceindex];
        vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentenceindex];
        
        //reconstruct token-offset in reverse index (making up for not having full index anymore in SelectivePatternModel)
        unordered_map<const EncAnyGram *, vector<int> > sourcetokenfwindex;
        unordered_map<int, vector<const EncAnyGram *> > sourcetokenrevindex;
        recompute_token_index(sourcetokenfwindex, sourcetokenrevindex, sentence_s2t.source, sourcepatterns);
        unordered_map<const EncAnyGram *, vector<int> > targettokenfwindex;
        unordered_map<int, vector<const EncAnyGram *> > targettokenrevindex;
        recompute_token_index(targettokenfwindex, targettokenrevindex, sentence_s2t.target, targetpatterns);
        
        
        
        //iterate over source sentence word by word
        for (int sourceindex = 0; sourceindex < sentence_s2t.source->length(); sourceindex++) {  
            //find source patterns valid at this position  
            if (sourcetokenrevindex.count(sourceindex)) {
                for (vector<const EncAnyGram*>::iterator iter_s = sourcetokenrevindex[sourceindex].begin(); iter_s !=  sourcetokenrevindex[sourceindex].end(); iter_s++) {
                    //now find what target patterns are aligned, and how well the aligment is (expressed through a score)
                    const EncAnyGram * sourcepattern = *iter_s;
                    const unsigned char sourcepatternsize = sourcepattern->n(); 
                    double bestscore = 0;             
                    const EncAnyGram * besttargetpattern = NULL;                           
                    for (vector<const EncAnyGram*>::iterator iter_t = targetpatterns->begin(); iter_t != targetpatterns->end(); iter_t++) {
                         const EncAnyGram * targetpattern = *iter_t;
                         const unsigned char targetpatternsize = sourcepattern->n();
                         const unsigned char maxpatternsize = sourcepatternsize ? (sourcepatternsize > targetpatternsize) : targetpatternsize;
                         
                         if (targettokenfwindex.count(targetpattern)) {
                             for (vector<int>::iterator iter2 = targettokenfwindex[targetpattern].begin(); iter2 != targettokenfwindex[targetpattern].end(); iter2++) { //loops over all occurences of the target pattern in the target sentence                         
                                 const int targetindex = *iter2; //begin index of target pattern
                                
                                 int aligned = 0;                                 
                                 bool firstsourcealigned = false;
                                 bool firsttargetaligned = false;
                                 bool lastsourcealigned = false;
                                 bool lasttargetaligned = false;
                                 //check alignment points in intersection                                                                      
                                 for (multimap<const unsigned char, const unsigned char>::const_iterator alignmentiter = sentence_i.alignment.begin(); alignmentiter != sentence_i.alignment.end(); alignmentiter++) {
                                    const unsigned char alignedsourceindex = alignmentiter->first -1; //-1 because of 1-based-indexing
                                    const unsigned char alignedtargetindex = alignmentiter->second -1; //-1 because of 1-based-indexing
                                 
                                    const bool insource = (alignedsourceindex >= sourceindex) && (alignedsourceindex < sourceindex + sourcepattern->n());
                                    const bool intarget = (alignedtargetindex >= targetindex) && (alignedtargetindex < targetindex + targetpattern->n());
                                     
                                     if ((!insource) || (!intarget)) {
                                        aligned--; break;
                                     } else {
                                        aligned++;
                                        if (alignedsourceindex == sourceindex) firstsourcealigned = true;
                                        if (alignedsourceindex == sourceindex + (sourcepatternsize - 1) ) lastsourcealigned = true;
                                        if (alignedtargetindex == targetindex) firsttargetaligned = true;
                                        if (alignedtargetindex == targetindex + (targetpatternsize - 1) ) lasttargetaligned = true;
                                     }                                                                        
                                 }
                                 if ((aligned < 0) || (!firstsourcealigned)  || (!lastsourcealigned) || (!firsttargetaligned)  || (!lasttargetaligned)) break;
                                 
                                 if (coocthreshold > 0) {
                                    //TODO: Implement cooc check
                                 }                                 
                                 
                                 double score = (double) aligned / maxpatternsize;
                                 if (score < 1) {
                                     //check alignment points in union 
                                     int halfaligned = -aligned; //start negative so intersection points are not counted twice
                                     for (multimap<const unsigned char, const unsigned char>::const_iterator alignmentiter = sentence_u.alignment.begin(); alignmentiter != sentence_u.alignment.end(); alignmentiter++) {
                                        const unsigned char alignedsourceindex = alignmentiter->first -1; //-1 because of 1-based-indexing
                                        const unsigned char alignedtargetindex = alignmentiter->second -1; //-1 because of 1-based-indexing
                                     
                                        const bool insource = (alignedsourceindex >= sourceindex) && (alignedsourceindex < sourceindex + sourcepattern->n());
                                        const bool intarget = (alignedtargetindex >= targetindex) && (alignedtargetindex < targetindex + targetpattern->n());
                                         
                                         if ((insource) && (intarget)) {
                                            halfaligned++;
                                         //} else if ((insource) || (intarget)) {
                                            //halfaligned--;
                                         }                                   
                                     }
                                     if (halfaligned > 0) {
                                        //adjust score favourably according to union point (each point weighing 4 times less than intersection alignment points)
                                        score = score + (double) halfaligned / (maxpatternsize * 4);
                                        if (score > 1) score = 1;
                                     }
                                }          
                                if (score > bestscore) { 
                                    //retain only the best target pattern given an occurrence of a source pattern
                                    bestscore = score;
                                    besttargetpattern = targetpattern;
                                }                       
                            }
                         }
                    } //iteration over all target patterns    
                    
                    
                    if ((besttargetpattern != NULL) && (bestscore >= alignscorethreshold)) {
                        //add alignment
                        alignmatrix[sourcepattern][besttargetpattern] += 1;
                        found++;
                        if ((sourcedecoder != NULL) && (targetdecoder != NULL)) {
                            cout << sourcepattern->decode(*sourcedecoder) << " ||| " << besttargetpattern->decode(*targetdecoder) << endl;                             
                        }
                    }
                                                    
                }
             }
          }  //sourceindex iterator
          cerr << " Found " << found << endl;
}



void AlignmentModel::extractgizapatterns(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold,ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {

    while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {         
        GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
        GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();    
        
        extractgizapatterns(sentence_s2t, sentence_t2s, gizamodel_s2t.index(), pairoccurrencethreshold, coocthreshold, alignscorethreshold, sourcedecoder,targetdecoder);
      } //alignment read        
      
      
      if (pairoccurrencethreshold > 0) prune(pairoccurrencethreshold);
            
      //normalize alignment matrix
      normalize();      
}

void recompute_token_index(unordered_map<const EncAnyGram *, vector<int> > & tokenfwindex, unordered_map<int, vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, const vector<const EncAnyGram*> * patterns ) {
    for (int i = 0; i < sentence->length(); i++) {
        for (vector<const EncAnyGram*>::const_iterator iter = patterns->begin(); iter !=  patterns->end(); iter++) {
            const EncAnyGram * anygram = *iter;
            bool match;
            if (anygram->isskipgram()) {
                match = sentence->match((EncSkipGram*) anygram, i);  
            } else {
                match = sentence->match((EncNGram*) anygram, i);
            }
            if (match) {
                //TODO: make sure anygram pointer remains valid (may need getkey?)
                tokenfwindex[anygram].push_back(i);
                tokenrevindex[i].push_back(anygram);
            }
        }
    }       
}
