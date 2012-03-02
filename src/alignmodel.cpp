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
	for (unordered_set<const EncAnyGram *>::const_iterator iter = targetpatterns.begin(); iter != targetpatterns.end(); iter++) {
			const EncAnyGram* targetgram = *iter;
	        multiset<uint32_t> * targetindex;
		    if (targetgram->gapcount() == 0) {
		       targetindex = &targetmodel->ngrams[*( (EncNGram*) targetgram)].sentences;
		    } else {
		       targetindex = &targetmodel->skipgrams[*( (EncSkipGram*) targetgram)].sentences;
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

CoocAlignmentModel::CoocAlignmentModel(CoocMode mode,SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const double absthreshold, const double probthreshold, bool bestonly, bool normalize, bool DEBUG) {
    this->mode = mode;
    this->absthreshold = absthreshold;
    this->probthreshold = probthreshold;
    this->normalize = normalize;
    this->bestonly = bestonly;
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

void AlignmentModel::simpletableoutput(ClassDecoder & sourceclassdecoder, ClassDecoder & targetclassdecoder, ostream * OUT, bool wordbased) {
	/* output a simple word-based lexicon, similar to the one used in moses (s2t, t2s) */
    for (unordered_map<const EncAnyGram*,unordered_map<const EncAnyGram*, double> >::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
        const EncAnyGram* sourcegram = iter->first;        
        map<double, const EncAnyGram*> sorted;        
        double total = 0;
        for (unordered_map<const EncAnyGram*, double>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        	if (wordbased && iter2->first->n() > 1) continue; 
        	sorted[-1 * iter2->second] = iter2->first;
        	total += iter2->second;            
        }
        for (map<double, const EncAnyGram*>::iterator iter2 = sorted.begin(); iter2 != sorted.end(); iter2++) {
			const EncAnyGram* targetgram = iter2->second;
			*OUT << sourcegram->decode(sourceclassdecoder) << " " << targetgram->decode(targetclassdecoder) << " " << ((-1 * iter2->first) / total) << endl;      
        }            
        *OUT << endl;
    }	
}

EMAlignmentModel::EMAlignmentModel(SelectivePatternModel * sourcemodel, SelectivePatternModel * targetmodel, const int MAXROUNDS, const double CONVERGEDTHRESHOLD, double probthreshold, bool bestonly, bool DEBUG) {
    int round = 0;    
    unsigned long c;
    double prevavdivergence = 0;
    bool converged = false;
    
    //initialise uniformly
    cerr << "  Initialisation step" << endl; 
    for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
    	uint32_t sentence = reviter_source->first;
		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
		if (targetmodel->reverseindex.count(sentence) > 0) {
			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
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
        for ( unordered_map<uint32_t,std::vector<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {        		
        		uint32_t sentence = reviter_source->first;
        		const vector<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			vector<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
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


AlignmentModel::AlignmentModel(const string & filename) {
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

