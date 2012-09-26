#include <alignmodel.h>

using namespace std;

const EncNullGram * NULLGRAM = new EncNullGram();
const unsigned char BOSCLASS = 3;
const unsigned char EOSCLASS = 4;

AlignmentModel::AlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, unsigned char leftsourcecontext, unsigned char rightsourcecontext, bool DEBUG) {    
    this->DEBUG = DEBUG;
    this->sourcemodel = sourcemodel;
    this->targetmodel = targetmodel;
    this->leftsourcecontext = leftsourcecontext;
    this->rightsourcecontext = rightsourcecontext;     
}


void AlignmentModel::intersect(AlignmentModel * reversemodel, double probthreshold, int bestn) {
	 //Compute intersection with reverse model
	for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		if (sourceiter->first == NULLGRAM) continue;
		const EncAnyGram * sourcegram = sourceiter->first;
		const EncAnyGram * revsourcegram = reversemodel->gettargetkey(sourcegram);
		if (revsourcegram == NULL) continue;
		double lowerbound = 0.0;
		list<double> bestq;

		
		for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			const EncAnyGram * revtargetgram = reversemodel->getsourcekey(targetgram);
			
			if ((revtargetgram != NULL) &&  (reversemodel->alignmatrix.count(revtargetgram) > 0) && (reversemodel->alignmatrix[revtargetgram].count(revsourcegram) > 0) && (listproduct(reversemodel->alignmatrix[revtargetgram][revsourcegram]) * listproduct(targetiter->second) >= probthreshold)) {
			    if (reversemodel->alignmatrix[revtargetgram][revsourcegram].size() != targetiter->second.size()) {
			        cerr << "AlignmentModel::intersect: Unable to compute intersection with reverse model, score vectors are of different length " << targetiter->second.size() << "vs " << reversemodel->alignmatrix[revtargetgram][revsourcegram].size();
			        exit(6);
			    }			    
			    //score vector [a,b] * [x,y]  == [a*x,b*y]
			    const int scores = targetiter->second.size();
			    for (int i = 0; i < scores; i++) {
			        alignmatrix[sourcegram][targetgram][i] = alignmatrix[sourcegram][targetgram][i] * reversemodel->alignmatrix[revtargetgram][revsourcegram][i];
			    } 			 
				const double p = listproduct(alignmatrix[sourcegram][targetgram]);
				if ((bestn) && ((p > lowerbound) || (bestq.size() < (size_t) bestn))) {
		    		orderedinsert(bestq, p);
		    		while (bestq.size() > (size_t) bestn) bestq.pop_front();
		    		lowerbound = *bestq.begin();		    			
		    	}
			} else {
				alignmatrix[sourcegram].erase(targetgram); //no intersection, prune
			} 			
		}		
		if (bestn) {
			//pruning stage
			for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
				const EncAnyGram * targetgram = targetiter->first;
				const double score = listproduct(targetiter->second);
				if (score < lowerbound) alignmatrix[sourcegram].erase(targetgram);
			}
		}
		if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);	
	}	 
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
		        alignmatrix[sourcegram][targetgram].push_back(coocvalue);
		    } else {
		    	prunedabs++;
		    }
		    totalcooc += coocvalue;	
		    if ((bestn) && ((coocvalue > lowerbound) || (bestq.size() < (size_t) bestn)))  {
		    	orderedinsert(bestq, coocvalue);		    	
				while (bestq.size() > (size_t) bestn) bestq.pop_front();
				lowerbound = *bestq.begin();
		    } 		    
	}
	if ((DEBUG) && (bestn))  cerr << "\t\tbest-n lowerbound=" << lowerbound << endl;				
    if ((totalcooc > 0) && (normalize || bestn || probthreshold > 0)) {
    	//normalisation and pruning step (based on probability threshold)
    	for (t_aligntargets::const_iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
    	    const double score = listproduct(iter->second); 
    		const double alignprob = (double) score / totalcooc;
    		if ((alignprob < probthreshold) || ((bestn) && (score < lowerbound)))  {
    			//prune
    			alignmatrix[sourcegram].erase(iter->first);
    			prunedprob++;
    		} else if (normalize) {
    			//normalise
    			alignmatrix[sourcegram][iter->first][0] = alignprob;
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
					if ((DONULL) && (alignmatrix[NULLGRAM][targetgram].empty()))  alignmatrix[NULLGRAM][targetgram].push_back(v);
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram].push_back(v);
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
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram][0]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram][0]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
		    				const double countvalue_null = alignmatrix[NULLGRAM][targetgram][0] / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}
						
		                for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {		                	
		    			    const EncAnyGram * sourcegram = *sourceiter;                                                                 
		    			    if ((alignmatrix.count(sourcegram) && alignmatrix[sourcegram].count(targetgram))) {
		                    	const double countvalue = alignmatrix[sourcegram][targetgram][0] / sentencetotal[targetgram];
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
        for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
        	const EncAnyGram * sourcegram = sourceiter->first;
        	for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
        		const EncAnyGram * targetgram = targetiter->first;
        		
				prevtransprob = targetiter->second[0];
                const double newtransprob = (double) count[sourcegram][targetgram] / total[sourcegram];
                if (alignmatrix[sourcegram][targetgram].empty()) {
                    alignmatrix[sourcegram][targetgram].push_back( newtransprob );
                } else {
                    alignmatrix[sourcegram][targetgram][0] = newtransprob;
                }
                
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
		for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
			const EncAnyGram * sourcegram = sourceiter->first;
			
			if (!bestn) {
				for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second[0] < probthreshold) {
						alignmatrix[sourcegram].erase(targetgram);
					} 
				}
			} else {
				double lowerbound = 0.0;
				list<double> bestq;
				for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					if ((targetiter->second[0] > lowerbound) || (bestq.size() < (size_t) bestn)) {
						orderedinsert(bestq, targetiter->second[0]);
						while (bestq.size() > (size_t) bestn) bestq.pop_front();
						lowerbound = *bestq.begin();
					}
				}
				for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
					const EncAnyGram * targetgram = targetiter->first;
					if (targetiter->second[0] < lowerbound) alignmatrix[sourcegram].erase(targetgram);
				}				
			}
			if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
		}    
    }
}

void AlignmentModel::normalize() {
	for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;
		//compute sum
		vector<double> sum;		
		for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			while (sum.size() < targetiter->second.size()) sum.push_back(0); //init
			for (unsigned int i = 0; i < targetiter->second.size(); i++) {
			    sum[i] += targetiter->second[i]; 
			}			
		}
		//normalize
		for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
			const EncAnyGram * targetgram = targetiter->first;
			for (unsigned int i = 0; i < targetiter->second.size(); i++) {
			    alignmatrix[sourcegram][targetgram][i] = alignmatrix[sourcegram][targetgram][i] / sum[i];
			 } 
		}
	}
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
	for (t_alignmatrix::iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;	
		//if sourcegram has parent relations;
		if (sourcemodel.rel_subsumption_parents.count(sourcegram) > 0) {
			//for sourceparentgram in parents_sourcegram:
			for (unordered_set<const EncAnyGram*>::const_iterator parentiterS = sourcemodel.rel_subsumption_parents[sourcegram].begin(); parentiterS != sourcemodel.rel_subsumption_parents[sourcegram].end(); parentiterS++) {
				const EncAnyGram * sourceparentgram = *parentiterS;
				//for t in targets_aligned_to_sourcegram:
				for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
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
	for (std::unordered_map<const EncAnyGram*,std::unordered_map<const EncAnyGram*, double> >::iterator iter = weightmatrix.begin(); iter != weightmatrix.end(); iter++) {
		const EncAnyGram * sourcegram = iter->first;	
		for (std::unordered_map<const EncAnyGram*, double>::const_iterator targetiter = weightmatrix[sourcegram].begin(); targetiter != weightmatrix[sourcegram].end(); targetiter++) {		
			const EncAnyGram * targetgram = targetiter->first;
			//convert the weights from to: -inf -- 0 -- inf  --> 0 -- 1 -- inf , so they can be applied directly
			
			const double weight = weightmatrix[sourcegram][targetgram];
            for (unsigned int i = 0; i < alignmatrix[sourcegram][targetgram].size(); i++) {									
			    const double a = alignmatrix[sourcegram][targetgram][i];
			    if (weight < 0) {					 
				    alignmatrix[sourcegram][targetgram][i] = pow(a, weight+1);
			    } else if (weight > 0) {
				    alignmatrix[sourcegram][targetgram][i] = pow(a, -1 * weight + 1);
			    }
			}
			adjustments++;							
		}
	}
	
	return adjustments;	
}			


/*
void AlignmentModel::trainEM2(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn, const bool DONULL, const bool INIT) {
    VERSION that considers multiple occurrences of the same type in a phrase pair
    

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


/*void AlignmentModel::growdiag(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s) {
    
}*/


int AlignmentModel::prune(const double prunethreshold) {
    int pruned = 0;
	for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourceiter->first;
	    for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
		    const EncAnyGram * targetgram = targetiter->first;
		    if (listproduct(targetiter->second) < prunethreshold) {
		        pruned++;
		        alignmatrix[sourcegram].erase(targetgram);
		    }
		}
		if (alignmatrix[sourcegram].size() == 0) alignmatrix.erase(sourcegram);
    }
    return pruned;
}

EncAnyGram * AlignmentModel::addcontext(const EncData * sentence, const EncAnyGram * focus, int sourceindex) {
    

    EncNGram * leftcontext;
    {
        EncNGram * leftcontext_body = NULL;
        EncNGram * leftcontext_dummies = NULL;
        int begin = sourceindex - leftsourcecontext;
        int length = leftsourcecontext;
        int leftdummies = 0;
        if (begin < 0) {
            leftdummies = -1 * begin;
            length = length - leftdummies;
            begin = 0;
        }                         
        if (length > 0) {
            leftcontext_body = sentence->slice(begin, length);                                                        
        } 
        if (leftdummies > 0) {
            unsigned char * buffer = new unsigned char[1024];
            char buffersize = 0;
            for (int i = 0; i < leftdummies; i++) { buffer[buffersize++] = BOSCLASS; buffer[buffersize++] = 0; } 
            leftcontext_dummies = new EncNGram(buffer, buffersize-1);
            delete [] buffer;
        } else {
            leftcontext = leftcontext_body;
        }                 
        
        if ((leftcontext_body != NULL) && (leftcontext_dummies != NULL)) {
            leftcontext = new EncNGram(*leftcontext_dummies + *leftcontext_body);
        } else if (leftcontext_body == NULL) {
            leftcontext = leftcontext_dummies;
        } else if (leftcontext_dummies == NULL) {
            leftcontext = leftcontext_body;
        }
                    
        if ((leftcontext_body != NULL) && (leftcontext != leftcontext_body)) delete leftcontext_body;
        if ((leftcontext_dummies != NULL) && (leftcontext != leftcontext_dummies)) delete leftcontext_dummies;
    }
    
    EncNGram * rightcontext;
    {
        EncNGram * rightcontext_body = NULL;
        EncNGram * rightcontext_dummies = NULL;
        int begin = sourceindex + focus->n();
        int length = rightsourcecontext;
        int rightdummies = 0;
        if (begin + length > sentence->length()) {
            rightdummies = sentence->length() - (begin + length);
            length = length - rightdummies;
            begin = 0;
        }                         
        if (length > 0) {
            rightcontext_body = sentence->slice(begin, length);                                                        
        } 
        if (rightdummies > 0) {
            unsigned char * buffer = new unsigned char[1024];
            char buffersize = 0;
            for (int i = 0; i < rightdummies; i++) { buffer[buffersize++] = EOSCLASS; buffer[buffersize++] = 0; } 
            rightcontext_dummies = new EncNGram(buffer, buffersize-1);
            delete [] buffer;
        } else {
            rightcontext = rightcontext_body;
        }                 
        
        if ((rightcontext_body != NULL) && (rightcontext_dummies != NULL)) {
            rightcontext = new EncNGram(*rightcontext_body + *rightcontext_dummies);
        } else if (rightcontext_body == NULL) {
            rightcontext = rightcontext_dummies;
        } else if (rightcontext_dummies == NULL) {
            rightcontext = rightcontext_body;
        }
                    
        if ((rightcontext_body != NULL) && (rightcontext != rightcontext_body)) delete rightcontext_body;
        if ((rightcontext_dummies != NULL) && (rightcontext != rightcontext_dummies)) delete rightcontext_dummies;
    }
    return focus->addcontext(leftcontext, rightcontext);
}


int AlignmentModel::extractgizapatterns(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
        GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
        GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);
        
        int found = 0;

        
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
                            
      3)  (implemented) 
      
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
        
        if (DEBUG) { 
            cerr << "DEBUG source token rev index size: " << sourcetokenrevindex.size() << endl;
            cerr << "DEBUG source token fw index size: " << sourcetokenfwindex.size() << endl;
            cerr << "DEBUG target token rev index size: " << targettokenrevindex.size() << endl;
            cerr << "DEBUG target token fw index size: " << targettokenfwindex.size() << endl;
        }
        
        //iterate over source sentence word by word
        for (int sourceindex = 0; sourceindex < sentence_s2t.source->length(); sourceindex++) {
            int count_s = 0;  
            const EncAnyGram * prevsourcepatternwithcontext = NULL;
            //find source patterns valid at this position  
            if (sourcetokenrevindex.count(sourceindex)) {
                for (vector<const EncAnyGram*>::iterator iter_s = sourcetokenrevindex[sourceindex].begin(); iter_s !=  sourcetokenrevindex[sourceindex].end(); iter_s++) {
                    //now find what target patterns are aligned, and how well the aligment is (expressed through a score)
                    count_s++;                    
                    
                    const EncAnyGram * sourcepattern = *iter_s;
                    unsigned char sourcepatternsize = sourcepattern->n();
                    
                    const EncAnyGram * sourcepatternwithcontext = NULL;
                    bool sourcepatternused = false; 
                    if (leftsourcecontext || rightsourcecontext) {
                        //extract context
                        
                        //make sure we re-use any previously used sourcepatternwithcontext rather than using this newly generated one
                        cout << endl; //DEBUG
                        sourcepatternwithcontext = addcontext(sentence_s2t.source, sourcepattern, sourceindex);
                        const EncAnyGram * key = getsourcekey(sourcepatternwithcontext, false); // allowfallback=false (do not fall back to sourcemodel, which doesn't know contexts)
                        if (key != NULL) {
                            delete sourcepatternwithcontext; //delete generated one
                            sourcepatternwithcontext = key; //use existing key
                            //cout << "CONTEXT EXISTS @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        //} else {
                            //cout << "CONTEXT NEW @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        }
                        
                    }                    
                     
                    double bestscore = 0;             
                    const EncAnyGram * besttargetpattern = NULL;
                    int count_t = 0;                           
                    
                    multiset<uint32_t> * sourcesentenceindex = NULL; //used only if coocthreshold > 0                         
                    if (coocthreshold > 0) {
                        if (sourcepattern->isskipgram()) {
                            sourcesentenceindex = &sourcemodel->skipgrams[*( (EncSkipGram*) sourcepattern)].sentences;    
                        } else {
                            sourcesentenceindex = &sourcemodel->ngrams[*( (EncNGram*) sourcepattern)].sentences;
                        }
                    }                    
                    
                    for (vector<const EncAnyGram*>::iterator iter_t = targetpatterns->begin(); iter_t != targetpatterns->end(); iter_t++) {
                         count_t++;                         
                         const EncAnyGram * targetpattern = *iter_t;
                         if (DEBUG) cerr << " @" << sentenceindex << "-" << count_s << "-" << count_t << endl;
                         if (DEBUG) cerr << "   DEBUG: [" << sourcepattern->decode(*sourcedecoder) << "] vs [" << targetpattern->decode(*targetdecoder) << "]" << endl;
                         
                         const unsigned char targetpatternsize = targetpattern->n();
                         const unsigned char maxpatternsize = (sourcepatternsize > targetpatternsize) ? sourcepatternsize : targetpatternsize;
                         
                          
                         if (targettokenfwindex.count(targetpattern)) {
                             if (coocthreshold > 0) {
                                //filter by co-occurrence threshold
                                multiset<uint32_t> * targetsentenceindex; //used only if coocthreshold > 0                         
                                if (targetpattern->isskipgram()) {
                                    targetsentenceindex = &targetmodel->skipgrams[*( (EncSkipGram*) targetpattern)].sentences;    
                                } else {
                                    targetsentenceindex = &targetmodel->ngrams[*( (EncNGram*) targetpattern)].sentences;
                                }                                                                                        			
                                double coocvalue = cooc(JACCARD, *sourcesentenceindex, *targetsentenceindex, coocthreshold);
                                if (DEBUG) cerr << "     DEBUG COOC=" << coocvalue << endl;
                                if (coocvalue < coocthreshold) {
                                        continue; //threshold not reached, skipping
                                }                                    
                             }                               
                             
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
                                     
                                     if ((insource) && (intarget)) {
                                        aligned++;
                                        if (alignedsourceindex == sourceindex) firstsourcealigned = true;
                                        if (alignedsourceindex == sourceindex + (sourcepatternsize - 1) ) lastsourcealigned = true;
                                        if (alignedtargetindex == targetindex) firsttargetaligned = true;
                                        if (alignedtargetindex == targetindex + (targetpatternsize - 1) ) lasttargetaligned = true;
                                     } else if ((insource) || (intarget)) {
                                        //alignment point outside pattern, discard
                                        aligned = -1; break;                                        
                                     }
                                 }
                                 if (DEBUG) cerr << "     DEBUG ALIGNED=" << aligned << "  ENDS: " << (int) firstsourcealigned << " " << (int) lastsourcealigned << " " << (int) firsttargetaligned << " " << (int) lasttargetaligned << endl;
                                 if ((aligned < 0) || (!firstsourcealigned)  || (!lastsourcealigned) || (!firsttargetaligned)  || (!lasttargetaligned)) break;
                                 if (DEBUG) cerr << "     DEBUG FOUND" << endl;
                                                        
                                 
                                 double score = (double) aligned / maxpatternsize;

                                 if (DEBUG) cerr << "     DEBUG score(1): " << aligned << " / " << (int) maxpatternsize << " = " << score << endl;
                                 
                                 
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
                                        if (DEBUG) cerr << "     DEBUG score(2): " << score << endl;
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
                        
                        //use pattern-in-context or without context as key in alignmatrix?
                        const EncAnyGram * key;
                        if (sourcepatternwithcontext != NULL) {
                            key = sourcepatternwithcontext;
                            //cout << "KEY=CONTEXT @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;

                        } else {
                            key = sourcepattern;
                            //cout << "KEY=NORMAL @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        }
                                        
                        //add alignment
                        if (alignmatrix[key][besttargetpattern].empty()) {
                            alignmatrix[key][besttargetpattern].push_back(1);
                        } else {
                            alignmatrix[key][besttargetpattern][0] += 1;
                        }
                        
                        /* DEBUG: if ((sourcepatternwithcontext != NULL) && (getsourcekey(key) != key)) {
                                cerr << "ERROR: CAN NOT RETRIEVE KEY" << endl;
                        } */
                        
                        sourcepatternused = true;
                        found++;
                        if ((sourcedecoder != NULL) && (targetdecoder != NULL)) {
                            cout << key->decode(*sourcedecoder) << " ||| " << besttargetpattern->decode(*targetdecoder) << " ||| " << bestscore << endl;
                        }
                    }
                    
                    if (sourcepatternwithcontext != NULL) {
                        if (sourcepatternused) {
                            sourcecontexts[sourcepattern].insert(sourcepatternwithcontext);
                        } else {
                            //delete sourcepatternwithcontext; //TODO: Reenable: causes segfault
                        }
                    }                                   
                }
             }
          }  //sourceindex iterator
         return found;
}



int AlignmentModel::extractgizapatterns(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold,ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
    unsigned int totalfound = 0;
    while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {         
        GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
        GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();    
        const int i = gizamodel_s2t.index();
        cerr << " @" << i << endl;
        int found = extractgizapatterns(sentence_s2t, sentence_t2s, i, pairoccurrencethreshold, coocthreshold, alignscorethreshold, sourcedecoder,targetdecoder);
        totalfound += found;
        cerr << "   found " << found << endl;
      } //alignment read        
      
      
      if (pairoccurrencethreshold > 0) prune(pairoccurrencethreshold);
            
      //normalize alignment matrix
      normalize();
      return totalfound;      
}

int AlignmentModel::extractskipgrams(const int absolutecoocthreshold) {
    unsigned int found = 0;
    //extract skipgram alignments on the basis of n-gram alignments found, and a graphmodel containing template information (irregardless of whether stored as transitive reduction or not)
    
    //TODO: what about valid skipgram to n-gram alignment?
    
    /*

    for ngram_s in alignmatrix:
        if |alignmatrix[ngram_s]| > 1: (nothing to generalise otherwise)
            skipgrams_s = get_templates_recursively(rel_templates_s[ngram_s])
            if skipgrams_s:            
                for ngram_t in alignmatrix[ngram_s]
                    skipgrams_t = get_templates_recursively(rel_templates_t[ngram_t])
                    if skipgrams_t:
                        for skipgram_s in skipgrams_s:
                            for skipgram_t in skipgrams_t:
                                submatrix[skipgram_s][skipgram_t]++
    prune 1-values from submatrix                     
       
    for skipgram_s in submatrix: //conflict resolution
        if |submatrix[skipgram_s]| == 1:
            skipgram_t = first and only
            alignmatrix[skipgram_s][skipgram_t] = submatrix[skipgram_s][skipgram_t]
        else:                                
            clusters = get_independent_cluster(skipgrams_t) //complete subgraphs
            for cluster in clusters:
              maxcooc = 0
              for skipgram_t in cluster:
                if cooc(skipgram_s, skipgram_t) > maxcooc:
                    maxcooc = _
                    best = skipgram_t
                elseif   == maxcooc:
                   if skipgram_t more abstract than best:
                        best = skipgram_t
              if best:
                 alignmatrix[skipgram_s][best] = submatrix[skipgram_s][best]                      
    */
    
    //check if both models contain template relations
    if ((!sourcemodel->HASTEMPLATES) || (!targetmodel->HASTEMPLATES)) {
        cerr << "WARNING: No template relations available in source and/or target model, unable to extract skipgrams" << endl;
        return 0;
    }
    
    unordered_map<const EncSkipGram *, unordered_map<const EncSkipGram *, uint16_t> > prealignmatrix; //temporary matrix
    
    cerr << "Extracting skipgrams";
    
    //harvest possible skipgram prealignments from (ngram) align matrix
    for (t_alignmatrix::const_iterator sourceiter = alignmatrix.begin(); sourceiter != alignmatrix.end(); sourceiter++) {
		const EncAnyGram * sourcegram = sourcemodel->getkey(sourceiter->first);
		if (sourcegram == NULL) {
		    cerr << "WARNING: Source-side pattern in align matrix not found in model! Skipping.." << endl;
		    continue; 
		}
		if (sourceiter->second.size() > 1) { //are there multiple candidates?		    
		    unordered_set<const EncSkipGram *> sourceskipgrams;
		    if (get_templates(sourcegram, sourcemodel, sourceskipgrams)) { //get all templates for this pattern, recursively
		        //cerr << ".";
                for (t_aligntargets::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
		            const EncAnyGram * targetgram = targetmodel->getkey(targetiter->first);
                    if (targetgram == NULL) {
                        cerr << "WARNING: Target-side pattern in align matrix not found in model! Skipping.." << endl;
                        continue; 
                    }		            
		            unordered_set<const EncSkipGram *> targetskipgrams;
		            if (get_templates(targetgram, targetmodel, targetskipgrams)) {
		                //cerr << ":";
		                //register all possible combinations of skipgrams as pre-alignments
		                for (unordered_set<const EncSkipGram *>::const_iterator sourceiter2 = sourceskipgrams.begin(); sourceiter2 != sourceskipgrams.end(); sourceiter2++) {
		                    const EncSkipGram * sourcegram2 = *sourceiter2;
		                    for (unordered_set<const EncSkipGram *>::const_iterator targetiter2 = targetskipgrams.begin(); targetiter2 != targetskipgrams.end(); targetiter2++) {
		                        cerr << ".";
		                        const EncSkipGram * targetgram2 = *targetiter2;
		                        prealignmatrix[sourcegram2][targetgram2]++;
		                    }		                    
		                } 		                    
		            }
		        }
		    } 
		}
	}
	cerr << endl;
	cerr << "Prealignments found for " << prealignmatrix.size() << " source-side skipgrams" << endl;
		
	//prune skipgram pre-alignments that occur only once
	for (unordered_map<const EncSkipGram*,unordered_map<const EncSkipGram*, uint16_t> >::const_iterator sourceiter = prealignmatrix.begin(); sourceiter != prealignmatrix.end(); sourceiter++) {
		const EncSkipGram * sourcegram = sourceiter->first;
	    for (unordered_map<const EncSkipGram*, uint16_t>::const_iterator targetiter = sourceiter->second.begin(); targetiter != sourceiter->second.end(); targetiter++) {
		    const EncSkipGram * targetgram = targetiter->first;
		    if (targetiter->second < absolutecoocthreshold) {
		        prealignmatrix[sourcegram].erase(targetgram);
		    }
		}
		if (prealignmatrix[sourcegram].size() == 0) prealignmatrix.erase(sourcegram);
    }
    cerr << "Prealignments left after pruning: " << prealignmatrix.size() << " source-side skipgrams" << endl;		

	 
	//assign final alignments to alignmatrix, do conflict resolution if multiple candidates exist
	for (unordered_map<const EncSkipGram*,unordered_map<const EncSkipGram*, uint16_t> >::const_iterator sourceiter = prealignmatrix.begin(); sourceiter != prealignmatrix.end(); sourceiter++) {
	 	const EncSkipGram * sourcegram = sourceiter->first;
	 	if (sourceiter->second.size() == 1) {
	 	    //no conflict
	 	    const EncSkipGram * targetgram = sourceiter->second.begin()->first;	 	    
	 	    alignmatrix[(const EncAnyGram *) sourcegram][(const EncAnyGram *) targetgram].push_back( prealignmatrix[sourcegram][targetgram] ) ;
	 	    found++;
	 	} else {
	 	    //possible conflict, multiple candidates
	 	    cerr << "MULTIPLE CANDIDATES: " << sourceiter->second.size() << endl;
            //find clusters of related skipgrams (complete subgraphs), relation through template/instantiation
	 	    vector<unordered_set<const EncSkipGram *>> clusters;
	 	    find_clusters(prealignmatrix[sourcegram], clusters, targetmodel); 
	 	    cerr << " clusters=" << clusters.size() << endl;
	 	    
	 	    const std::multiset<uint32_t> * sourcesentences = &sourcemodel->skipgrams[*( (EncSkipGram*) sourcegram)].sentences;
	 	    
	 	    for (vector<unordered_set<const EncSkipGram *>>::iterator clusteriter = clusters.begin(); clusteriter != clusters.end(); clusteriter++) {
	 	        double bestcooc = 0; 
	 	        unordered_set<const EncSkipGram *> cluster = *clusteriter;
	 	        unordered_set<const EncSkipGram *>::iterator best_iter = cluster.begin();
                for (unordered_set<const EncSkipGram *>::iterator targetiter = cluster.begin(); targetiter != cluster.end(); targetiter++) {
                    const EncAnyGram * targetgram = *targetiter;
		            const std::multiset<uint32_t> * targetsentences = &targetmodel->skipgrams[*( (EncSkipGram*) targetgram)].sentences;		                        
                    double coocvalue = cooc(JACCARD, *sourcesentences, *targetsentences);
                    if (coocvalue >= bestcooc) {
                        bestcooc = coocvalue;
                        best_iter = targetiter;
                    } else if (coocvalue == bestcooc) {
                        //multiple candidates with the same cooc value, prefer the most abstracting one
                        const  EncSkipGram * aligncandidate = *best_iter; 
                        if (((const EncSkipGram *) targetgram)->gapratio() > aligncandidate->gapratio()) {
                            best_iter = targetiter;
                        }
                    }
                }	 	        
                const EncSkipGram * bestaligntarget = *best_iter;
                alignmatrix[(const EncAnyGram *) sourcegram][(const EncAnyGram *) bestaligntarget].clear();
                alignmatrix[(const EncAnyGram *) sourcegram][(const EncAnyGram *) bestaligntarget].push_back(  prealignmatrix[sourcegram][bestaligntarget] );
                found++;
	 	    }
	 	    
	 	    
	 	}
	 	
	 }
    
     cerr << "Alignments found: " << found << endl;
     
     
     if (found) normalize();
     
     return found;    
}


void recompute_token_index(unordered_map<const EncAnyGram *, vector<int> > & tokenfwindex, unordered_map<int, vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, const vector<const EncAnyGram*> * patterns, bool includeskipgrams ) {
    for (int i = 0; i < sentence->length(); i++) {
        for (vector<const EncAnyGram*>::const_iterator iter = patterns->begin(); iter !=  patterns->end(); iter++) {
            const EncAnyGram * anygram = *iter;
            if ((!includeskipgrams) && (anygram->isskipgram())) continue;
            bool match;
            if (anygram->isskipgram()) {
                match = sentence->match((EncSkipGram*) anygram, i);  
            } else {
                match = sentence->match((EncNGram*) anygram, i);
            }
            if (match) {
                tokenfwindex[anygram].push_back(i);
                tokenrevindex[i].push_back(anygram);
            }
        }
    }       
}


size_t get_templates(const EncAnyGram * anygram, SelectivePatternModel * model, unordered_set<const EncSkipGram *> & container) {    
    if (model->rel_templates.count(anygram)) {
        for (unordered_set<const EncAnyGram*>::iterator iter = model->rel_templates[anygram].begin(); iter != model->rel_templates[anygram].end(); iter++) {
               const EncSkipGram * tmplate = (const EncSkipGram*) *iter;   
               container.insert(tmplate);
               get_templates(tmplate, model, container);
        }   
    }
    return container.size();
}  

void find_clusters(unordered_map<const EncSkipGram*,uint16_t> skipgrams, vector<unordered_set<const EncSkipGram*> > & clusters, SelectivePatternModel * model ) {
    /*
        for skipgram in skipgrams:
            for cluster in clusters:
                if related(skipgram, cluster[0]):
                    targetcluster = cluster
                if targetcluster:
                    clusters.add( [skipgram] )        
    */
    for (unordered_map<const EncSkipGram*,uint16_t>::iterator sgiter = skipgrams.begin(); sgiter != skipgrams.end(); sgiter++ ) {
        const EncAnyGram * skipgram = model->getkey(sgiter->first);
        if (skipgram == NULL) {
            cerr << "WARNING: [in find_clusters()] Skipgram not found in model!" << endl;
            continue;
        }
        unordered_set<const EncAnyGram *> relationcandidates;
        model->getrelations( model->rel_templates, skipgram, relationcandidates); //forward search
        //model->getrelations( model->rel_instances, skipgram, relationcandidates); //backward search
        cerr << "DEBUG: " << relationcandidates.size() << endl;
                
        bool found = false;                
        for (vector<unordered_set<const EncSkipGram*> >::iterator clusteriter = clusters.begin(); clusteriter != clusters.end(); clusteriter++) {
            const EncSkipGram * refskipgram = *(clusteriter->begin());
            //check if the two are related
            bool related = false;
            
            //cerr << endl << "REF: ";
            //refskipgram->out();
            //cerr << endl;
            for (unordered_set<const EncAnyGram *>::iterator iter2 = relationcandidates.begin(); iter2 != relationcandidates.end(); iter2++) {
                const EncSkipGram * candidate = (const EncSkipGram*) *iter2;
                //cerr << endl << "CDD: " << endl;
                //candidate->out();
                //cerr << endl;
                if (*refskipgram == *candidate) {
                    //cerr << "MATCH!!!" << endl;
                    related = true; 
                    break;
                } 
            } 
        
            if (related) {
                //cerr << "FOUND RELATED CLUSTER" << endl;
                found = true;
                clusteriter->insert((const EncSkipGram *) skipgram);
                break;
            }            
        }
        if (!found) {
            unordered_set<const EncSkipGram*> newcluster;
            newcluster.insert((const EncSkipGram *) skipgram);
            clusters.push_back(newcluster);
        }
    }  
}


unsigned int AlignmentModel::prunepatternmodel(IndexedPatternModel & patternmodel, double threshold) {
    unsigned int pruned = 0; 
    cerr << "SIZE=" << alignmatrix.size() << endl;
    for (unordered_map<const EncNGram,NGramData >::iterator iter = patternmodel.ngrams.begin(); iter != patternmodel.ngrams.end(); iter++) {
        const EncNGram ngram = iter->first;
        const EncAnyGram * anygram = getsourcekey(&ngram);
        if (anygram != NULL) {
            double maxscore = 0;
            for (t_aligntargets::iterator iter2 = alignmatrix[anygram].begin(); iter2 != alignmatrix[anygram].end(); iter2++) {
                const double score = listproduct(iter2->second);
                if (score > maxscore) maxscore = score;
            }
            if (maxscore < threshold) {
                patternmodel.ngrams.erase(ngram);
                pruned++;
            }        
        } else {
            patternmodel.ngrams.erase(ngram);
            pruned++;        
        }
    }
   for (unordered_map<const EncSkipGram,SkipGramData >::iterator iter = patternmodel.skipgrams.begin(); iter != patternmodel.skipgrams.end(); iter++) {
        const EncSkipGram skipgram = iter->first;
        const EncAnyGram * anygram = getsourcekey(&skipgram);
        if (anygram != NULL) {
            double maxscore = 0;
            for (t_aligntargets::iterator iter2 = alignmatrix[anygram].begin(); iter2 != alignmatrix[anygram].end(); iter2++) {
                const double score = listproduct(iter2->second);
                if (score > maxscore) maxscore = score;
            }
            if (maxscore < threshold) {
                patternmodel.skipgrams.erase(skipgram);
                pruned++;
            }        
        } else {
            patternmodel.skipgrams.erase(skipgram);
            pruned++;        
        }
    }
 
    return pruned;
}



AlignmentModel::AlignmentModel(const string & s2tfilename, const string & t2sfilename, const double s2tthreshold, const double t2sthreshold, const double productthreshold, bool DEBUG) {
    this->DEBUG = DEBUG;
    //TODO: optimize reading directly into target matrix?
    AlignmentModel s2tmodel = AlignmentModel(s2tfilename);
    AlignmentModel t2smodel = AlignmentModel(t2sfilename);
    sourcemodel = NULL;
    targetmodel = NULL;
    this->leftsourcecontext = s2tmodel.leftsourcecontext;
    this->rightsourcecontext = s2tmodel.rightsourcecontext;     
    load(s2tmodel, t2smodel, s2tthreshold, t2sthreshold, productthreshold);    
}


AlignmentModel::AlignmentModel(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold, const double t2sthreshold, const double productthreshold, bool DEBUG) {
    this->DEBUG = DEBUG;
    sourcemodel = NULL;
    targetmodel = NULL;
    load(s2tmodel, t2smodel, s2tthreshold, t2sthreshold, productthreshold);
}

void AlignmentModel::load(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold, const double t2sthreshold, const double productthreshold) {
    for (t_alignmatrix::iterator iter = s2tmodel.alignmatrix.begin(); iter != s2tmodel.alignmatrix.end(); iter++) {
        const EncAnyGram * copysource_s2t = s2tmodel.getsourcekey(iter->first);
        const EncAnyGram * copysource_t2s = t2smodel.gettargetkey(iter->first);
        if ((copysource_s2t != NULL) && (copysource_t2s != NULL)) {
            for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                const EncAnyGram * copytarget_s2t = s2tmodel.gettargetkey(iter2->first);
                const EncAnyGram * copytarget_t2s = t2smodel.getsourcekey(iter2->first);                
                if ((copytarget_s2t != NULL) && (copytarget_t2s != NULL)  && (listproduct(iter2->second) >= s2tthreshold)) {                    
                    if (t2smodel.alignmatrix.count(copytarget_t2s) && t2smodel.alignmatrix[copytarget_t2s].count(copysource_t2s)) {                                                        
                        if ((t2sthreshold == 0) || (listproduct(t2smodel.alignmatrix[copytarget_t2s][copysource_t2s]) >= t2sthreshold)) {                            
                            if ((productthreshold == 0) || (listproduct(iter2->second) * listproduct(t2smodel.alignmatrix[copytarget_t2s][copysource_t2s]) >= productthreshold)) {                                
                                if (copysource_t2s->isskipgram()) {                                    
                                    alignmatrix[copysource_t2s]; //should be enough to insert the pointer
                                    //EncSkipGram skipgram = *((const EncSkipGram *) copysource_t2s);
                                    //sourceskipgrams.insert(skipgram); 
                                } else {                                    
                                    alignmatrix[copysource_t2s]; //should be enough to insert the pointer
                                    //EncNGram ngram = *((const EncNGram *) copysource_t2s);
                                    //sourcengrams.insert(ngram);
                                }
                                if (copytarget_t2s->isskipgram()) {
                                    EncSkipGram skipgram = *((const EncSkipGram *) copytarget_t2s);
                                    targetskipgrams.insert(skipgram); 
                                } else {
                                    EncNGram ngram = *((const EncNGram *) copytarget_t2s);
                                    targetngrams.insert(ngram);
                                }   
                                for (unsigned int i = 0; i < iter2->second.size(); i++) {
                                    alignmatrix[getsourcekey(copysource_t2s)][gettargetkey(copytarget_t2s)].push_back(iter2->second[0]);
                                }
                                for (unsigned int i = 0; i < t2smodel.alignmatrix[copytarget_t2s][copysource_t2s].size(); i++) {
                                    alignmatrix[getsourcekey(copysource_t2s)][gettargetkey(copytarget_t2s)].push_back(t2smodel.alignmatrix[copytarget_t2s][copysource_t2s][i]);
                                }   
                            }                  
                        }
                    }        
                }                            
            }   
        } 
    }    

}


AlignmentModel::AlignmentModel(const string & filename, bool logprobs, bool allowskipgrams, const int bestn, bool DEBUG) {
    this->DEBUG = DEBUG;
    load(filename,logprobs, allowskipgrams, bestn);
}

void AlignmentModel::load(const string & filename, bool logprobs, bool allowskipgrams, const int bestn) {
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
    bool multiscore;
    if ((model_id < (unsigned int) ALIGNMENTMODEL) || (model_id > (unsigned int) ALIGNMENTMODEL + 99))  {
        cerr << "File '" << filename << "' is not an alignment model" << endl;
        exit(6);
    }
    if ((model_id >= (unsigned int) ALIGNMENTMODEL) && (model_id <= 101))  { 
        multiscore = false; //backward compatibility with old models
    } else {
        multiscore = true;
    }
    if (model_id >= (unsigned int) ALIGNMENTMODEL + 4) { 
        f.read( (char*) &leftsourcecontext, sizeof(unsigned char));
        f.read( (char*) &rightsourcecontext, sizeof(unsigned char));
    }        
    f.read( (char*) &sourcecount, sizeof(uint64_t));  
    

    char gapcount;    
    for (unsigned int i = 0; i < sourcecount; i++) {	    
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
        bool sourceisskipgram = false;
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            const EncNGram * ngram = new EncNGram(&f); //read from file
            sourcegram = getsourcekey((EncAnyGram*) ngram); //does the key already exist?
            if (sourcegram == NULL) {            
                //no
                alignmatrix[(const EncAnyGram*) ngram];
                sourcegram = getsourcekey((const EncAnyGram*) ngram);
                if (sourcegram == NULL) { cerr << "INTERNAL ERROR: sourcegram still not found after insertion! Should never happen!"; exit(6); }
            	//sourcengrams.insert(ngram);            	            
            } else {
                //yes
                delete ngram; //ngram not necessary
            }               
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            const EncSkipGram * skipgram = new EncSkipGram( &f, gapcount); //read from file
            if (allowskipgrams) {              
                sourcegram = getsourcekey((EncAnyGram*) skipgram);  //does the key already exist?
                if (sourcegram == NULL) {
                    alignmatrix[(const EncAnyGram*) skipgram];
                	//sourceskipgrams.insert(skipgram);
                	sourcegram = getsourcekey((const EncAnyGram*) skipgram);
                    if (sourcegram == NULL) { cerr << "INTERNAL ERROR: sourcegram still not found after insertion! Should never happen!"; exit(6); }             	
                } else {
                    delete skipgram;
                }                 
            } else {
                delete skipgram;
            }
            sourceisskipgram = true;                  
        }        
        if ( (leftsourcecontext || rightsourcecontext) && ( !(sourceisskipgram && !allowskipgrams) ) ) {
            //deal with source-side context
            const EncAnyGram * sourcegramfocus = sourcegram->slice(leftsourcecontext, sourcegram->n() - leftsourcecontext - rightsourcecontext);
            sourcecontexts[sourcegramfocus].insert(sourcegram);
        }                                                   
        uint64_t targetcount;
        
        f.read( (char*) &targetcount, sizeof(uint64_t));
        if (DEBUG)  cerr << "\t--" << targetcount << ":";
        double lowerbound = 0.0;
        list<double> bestq;
        for (unsigned int j = 0; j < targetcount; j++) {
        	const EncAnyGram * targetgram = NULL;   
            f.read(&gapcount, sizeof(char));
            bool targetisskipgram = false;	    
		    if (gapcount == 0) {
		        if (DEBUG)  cerr << "\tNGRAM";
		        EncNGram ngram = EncNGram(&f); //read from file
		        if (DEBUG)  cerr << "\tread";
		        if (!gettargetkey((EncAnyGram*) &ngram)) {
		        	targetngrams.insert(ngram);		        	
		        }   
		        targetgram = gettargetkey((EncAnyGram*) &ngram);                                           
		    } else {
		        if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
		        EncSkipGram skipgram = EncSkipGram( &f, gapcount); //read from file
		        if (allowskipgrams) {		                          		        
		            if (!gettargetkey((EncAnyGram*) &skipgram)) {
		            	targetskipgrams.insert(skipgram);		        	
		            }   
		            targetgram = gettargetkey((EncAnyGram*) &skipgram);
		        }
		        targetisskipgram = true;                    
		    }		    
		    char scores;
		    if (multiscore) {
		        f.read((char*) &scores, sizeof(char));
		    } else {
		        scores = 1;
		    }		    
		
            double p;
		    for (int i = 0; i < scores; i++) {            
		        f.read((char*) &p, sizeof(double));
		        if ((p > 0) && (logprobs)) p = log(p); //base e		        
		        if ((allowskipgrams) || ((!sourceisskipgram) && (!targetisskipgram))) {  
           		    if ((sourcegram == NULL) || (targetgram == NULL)) {
		             	cerr << "SOURCEGRAM or TARGETGRAM is NULL";
		            	exit(6);
		            }		 		      
		            alignmatrix[sourcegram][targetgram].push_back(p);
		        }		    
		    }		  			    
		    if (bestn) {
		        if (logprobs) {
		           p = listsum(alignmatrix[sourcegram][targetgram]);   
		        } else {
		           p = listproduct(alignmatrix[sourcegram][targetgram]);
		        }
		        if ((p > lowerbound) || (bestq.size() < (size_t) bestn)) {
    	    		orderedinsert(bestq, p);
	    		    while (bestq.size() > (size_t) bestn) bestq.pop_front();
	    		    lowerbound = *bestq.begin();		    			
		        }
		    }
        }
        if (bestn) {
			//pruning stage
			for (t_aligntargets::iterator iter = alignmatrix[sourcegram].begin(); iter != alignmatrix[sourcegram].end(); iter++) {
			    double score;
			    if (logprobs) {
			        score = listsum(iter->second);
			    } else {
			        score = listproduct(iter->second);
			    }			    
				if (score < lowerbound) {
					alignmatrix[sourcegram].erase(iter->first);
				}				
			}		
			//TODO: prune orphaned sourcecontexts
		}  
	}
    f.close();
}

AlignmentModel::AlignmentModel(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs, bool DEBUG) {
    //load from moses-style phrasetable file
    this->DEBUG = DEBUG;
    sourcemodel = NULL;
    targetmodel = NULL;    
    leftsourcecontext = rightsourcecontext = 0;
		
    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    while (!f.eof()) {
        string line;
        getline(f, line);
        int mode = 0;
        string source = "";
        string target = "";
        vector<double> scores;
        int begin = 0;        
        for (unsigned int i = 0; i < line.size(); i++) {
            if (line.substr(i,5) == " ||| ") {
                if (mode == 0) {
                    source = line.substr(begin, i - begin);
                } else if (mode == 1) {
                    target = line.substr(begin, i - begin);
                }
                begin = i+5;
                mode++;
            } else if ((mode == 2) && (line[i] == ' ')) {
                double score = atof(line.substr(begin, i - begin).c_str());
                if ((score > 0) && (logprobs)) score = log(score); //base e
                scores.push_back(score);
                begin = i + 1;
            }
        }
        if ((mode == 2) && ((size_t) begin < line.size())) {
            double score = atof(line.substr(begin, line.size() - begin).c_str());
            if ((score > 0) && (logprobs)) score = log(score); //base e
            scores.push_back(score);            
        }
        if ((!source.empty()) && (!target.empty()) && (!scores.empty())) {
            //add to phrasetable
            EncAnyGram * sourcegram = sourceencoder->input2anygram(source,true,true);
            
            /*
            if (sourcegram->isskipgram()) {
                EncSkipGram skipgram = *((EncSkipGram *) sourcegram);
                if (!getsourcekey(sourcegram)) sourceskipgrams.insert(skipgram);            	
            } else {
                EncNGram ngram = *((EncNGram *) sourcegram);
                if (!getsourcekey(sourcegram)) sourcengrams.insert(ngram);
            }*/      
            const EncAnyGram * sourcegramkey = getsourcekey(sourcegram);
            if (sourcegramkey == NULL) {
                alignmatrix[sourcegram]; //does the insert
                sourcegramkey = sourcegram; //getsourcekey(sourcegram);
            }
            
                        
            EncAnyGram * targetgram = targetencoder->input2anygram(target,true,true);
            if (targetgram->isskipgram()) {
                EncSkipGram skipgram = *((EncSkipGram *) targetgram);
                if (!gettargetkey(targetgram)) targetskipgrams.insert(skipgram);            	
            } else {
                EncNGram ngram = *((EncNGram *) targetgram);
                if (!gettargetkey(targetgram)) targetngrams.insert(ngram);
            }
            const EncAnyGram * targetgramkey = gettargetkey(targetgram);            
            
            alignmatrix[sourcegramkey][targetgramkey] = scores;
        }
    }    
}


void AlignmentModel::save(const string & filename) {
	const unsigned char check = 0xff;
	const char czero = 0;
		
    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "File does not exist: " << filename << endl;
       exit(3);
    }
    
    uint64_t _id = ALIGNMENTMODEL + ALIGNMENTMODELVERSION;
    f.write( (char*) &_id, sizeof(uint64_t));
    
    f.write( (char*) &leftsourcecontext, sizeof(unsigned char));
    f.write( (char*) &rightsourcecontext, sizeof(unsigned char));
            
    uint64_t sourcecount = alignmatrix.size();
    if (alignmatrix.count(NULLGRAM) > 0) sourcecount--;
    f.write( (char*) &sourcecount, sizeof(uint64_t));   
    

    for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;

    	
		f.write((char*) &check, sizeof(char));
			  
    	const EncAnyGram * sourcegram = iter->first; // getsourcekey(iter->first);
    	/*if (sourcegram == NULL) { //check not needed? pointer should be right anyhow? getsourcekey == NULL will happen with source-context!
    	    cerr << "AlignmentModel::save(): Source key not found! This should not happen!"; 
    	    exit(6); 
    	}*/
    	if (sourcegram->isskipgram()) {
    		const EncSkipGram * skipgram = (const EncSkipGram*) sourcegram;
    		skipgram->writeasbinary(&f);
    	} else {
    	    const EncNGram * ngram = (const EncNGram*) sourcegram;
    	    f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
    		ngram->writeasbinary(&f);    		
    	}                
    	uint64_t targetcount = iter->second.size();
    	f.write( (char*) &targetcount, sizeof(uint64_t));
    	            
        for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	const EncAnyGram* targetgram = gettargetkey(iter2->first);
        	if (targetgram == NULL) {
        	    cerr << "AlignmentModel::save(): Target key not found! This should not happen!";
        	    exit(6); 
        	}        	        	
        	        	
        	if (targetgram->isskipgram()) {
    			const EncSkipGram * skipgram = (const EncSkipGram*) targetgram;
	    		skipgram->writeasbinary(&f);
	    	
    		} else {
    	    	const EncNGram * ngram = (const EncNGram*) targetgram;     	
    			f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
	    		ngram->writeasbinary(&f);	    		
    		}
    		
    		const char scores = (char) iter2->second.size();
    		f.write( (char*) &scores, sizeof(char));       
    		for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {         
    		    const double p = *iter3;
        	    f.write( (char*) &p, sizeof(double));
        	}
        }

    }    
    f.close();	
}


void AlignmentModel::decode(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT, bool mosesformat) {
    string delimiter = "\t";
    if (mosesformat) delimiter = " ||| ";
    for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
    	if (iter->first == NULLGRAM) continue;
        const EncAnyGram* sourcegram = iter->first;
        //cerr << "DEBUG: SIZE=" << iter->second.size() << " @" << (size_t) iter->first << " #" << iter->first->hash() << endl;         
        for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram* targetgram = iter2->first;
            *OUT << sourcegram->decode(sourceclassdecoder) << delimiter;
            *OUT << targetgram->decode(targetclassdecoder) << delimiter;            
            for (vector<double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                *OUT << *iter3 << ' ';
            }
            *OUT << endl;
        }
    }
}

const EncAnyGram * AlignmentModel::getfocuskey(const EncAnyGram * key) {
    t_contexts::iterator keyiter = sourcecontexts.find(key);
    if (keyiter != sourcecontexts.end()) {
        return keyiter->first;    
    } else {
        return NULL;
    }    
}
    

const EncAnyGram * AlignmentModel::getsourcekey(const EncAnyGram * key,  bool allowfallback) {
    //cout << "SEARCHING #" << (size_t) key->hash() << endl;
    t_alignmatrix::iterator keyiter = alignmatrix.find(key);
    if (keyiter != alignmatrix.end()) {
        //cout << "FOUND @" << (size_t) keyiter->first << endl;
        return keyiter->first;
    } else if ((sourcemodel != NULL) && (allowfallback)) {    
        return sourcemodel->getkey(key);
    } else {
        //cout << "NOT FOUND" << endl;
        return NULL;
    }
    /*            
    if (sourcemodel != NULL) return sourcemodel->getkey(key);
    if (key->gapcount() == 0) {
        std::unordered_set<EncNGram>::iterator iter = sourcengrams.find(*( (EncNGram*) key) );
        if (iter != sourcengrams.end()) {
            return &(*iter);
        } else {
            return NULL;
        }
    } else {
        std::unordered_set<EncSkipGram>::iterator iter = sourceskipgrams.find(*( (EncSkipGram*) key) );
        if (iter != sourceskipgrams.end()) {
            return &(*iter);
        } else {
            return NULL;
        }        
    }
    */
}


const EncAnyGram * AlignmentModel::gettargetkey(const EncAnyGram* key) {
    if (targetmodel != NULL) return targetmodel->getkey(key);
    if (key->gapcount() == 0) {
        std::unordered_set<EncNGram>::iterator iter = targetngrams.find(*( (EncNGram*) key) );
        if (iter != targetngrams.end()) {
            return &(*iter);
        } else {
            return NULL;
        }
    } else {
        std::unordered_set<EncSkipGram>::iterator iter = targetskipgrams.find(*( (EncSkipGram*) key) );
        if (iter != targetskipgrams.end()) {
            return &(*iter);
        } else {
            return NULL;
        }        
    }
}






/**************************** EXPERIMENTAL EM MODEL **********************************/



/*


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
					if (DONULL) alignmatrix[NULLGRAM][targetgram][0] = v;
		            for (vector<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
					    const EncAnyGram * sourcegram = *sourceiter;
						alignmatrix[sourcegram][targetgram][0] = v;
					}
				}
			}
		}
	}    
}

void EMAlignmentModel2::train(const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, const int bestn) {
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
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram][0]; //compute sum over all source conditions for a targetgram under consideration																	 
							}
						}
        			}
		    			
		    			
		    			
		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (vector<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {  
		    			const EncAnyGram * targetgram = *targetiter;
		    			
		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram][0]; //belongs to previous step technically, but moved into this loop for efficieny
		    			
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
*/


