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
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
}



AlignmentModel::AlignmentModel(unsigned char leftsourcecontext, unsigned char rightsourcecontext, int ptsfield, bool DEBUG) {
    this->DEBUG = DEBUG;
    this->leftsourcecontext = leftsourcecontext;
    this->rightsourcecontext = rightsourcecontext;
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
    this->ptsfield = ptsfield;
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
			        throw InternalError();
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
			for (unordered_set<const EncAnyGram*>::const_iterator reviter = targetmodel->reverseindex[sentencenumber].begin(); reviter != targetmodel->reverseindex[sentencenumber].end(); reviter++) {
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
		for ( unordered_map<uint32_t,unordered_set<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {
			uint32_t sentence = reviter_source->first;
			const unordered_set<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
			if (targetmodel->reverseindex.count(sentence) > 0) {
				unordered_set<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
				if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
				for (unordered_set<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {
					const EncAnyGram * targetgram = *targetiter;
					if ((DONULL) && (alignmatrix[NULLGRAM][targetgram].empty()))  alignmatrix[NULLGRAM][targetgram].push_back(v);
		            for (unordered_set<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
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
        for ( unordered_map<uint32_t,std::unordered_set<const EncAnyGram*> >::const_iterator reviter_source = sourcemodel->reverseindex.begin(); reviter_source != sourcemodel->reverseindex.end(); reviter_source++) {   //iterate over sentences

        		uint32_t sentence = reviter_source->first;
        		const unordered_set<const EncAnyGram*> * sourcepatterns = &reviter_source->second;
        		if (targetmodel->reverseindex.count(sentence) > 0) {
        			unordered_set<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentence];
        			if ((DEBUG) || (sentence % 1000 == 0)) cerr << "@" << sentence << " (" << sourcepatterns->size() << "x" << targetpatterns->size() << ")" << endl;
        			//compute sentencetotal for normalisation later in count step, sum_s(p(t|s))
        			unordered_map<const EncAnyGram*, double> sentencetotal;
        			for (unordered_set<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
		    			const EncAnyGram * sourcegram = *sourceiter;
		    			if (alignmatrix.count(sourcegram)) {
							for (unordered_set<const EncAnyGram*>::iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {
									const EncAnyGram * targetgram = *targetiter;
									if (alignmatrix[sourcegram].count(targetgram)) sentencetotal[targetgram] += alignmatrix[sourcegram][targetgram][0]; //compute sum over all source conditions for a targetgram under consideration
							}
						}
        			}



		            //collect counts to estimate improved model   (for evidence that a targetgram is aligned to a sourcegram)
		            for (unordered_set<const EncAnyGram*>::const_iterator targetiter = targetpatterns->begin(); targetiter != targetpatterns->end(); targetiter++) {
		    			const EncAnyGram * targetgram = *targetiter;

		    			//the null condition:
		    			if (DONULL) {
		    				if (alignmatrix[NULLGRAM].count(targetgram)) sentencetotal[targetgram] += alignmatrix[NULLGRAM][targetgram][0]; //belongs to previous step technically, but moved into this loop for efficieny

		    				const double countvalue_null = alignmatrix[NULLGRAM][targetgram][0] / sentencetotal[targetgram];
		                	count[NULLGRAM][targetgram] += countvalue_null;
							total[NULLGRAM] += countvalue_null;
						}

		                for (unordered_set<const EncAnyGram*>::const_iterator sourceiter = sourcepatterns->begin(); sourceiter != sourcepatterns->end(); sourceiter++) {
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

void AlignmentModel::normalize(t_alignmatrix * matrix) {
    if (matrix == NULL) matrix = &alignmatrix;
	for (t_alignmatrix::const_iterator sourceiter = matrix->begin(); sourceiter != matrix->end(); sourceiter++) {
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
			    (*matrix)[sourcegram][targetgram][i] = (*matrix)[sourcegram][targetgram][i] / sum[i];
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
        int begin = sourceindex + focus->n(); //begin of the rightcontext
        int length = rightsourcecontext;
        int rightdummies = 0;
        if (begin + length >= sentence->length()) {
            rightdummies = (begin + length) - sentence->length();
            length = length - rightdummies;
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



    EncAnyGram * result = focus->addcontext(leftcontext, rightcontext);

    delete leftcontext;
    delete rightcontext;

    return result;
}



EncAnyGram * AlignmentModel::addcontext(const EncData * sentence, const EncAnyGram * focus, int sourceindex, int leftsourcecontext, int rightsourcecontext) {

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
        int begin = sourceindex + focus->n(); //begin of the rightcontext
        int length = rightsourcecontext;
        int rightdummies = 0;
        if (begin + length >= sentence->length()) {
            rightdummies = (begin + length) - sentence->length();
            length = length - rightdummies;
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



    EncAnyGram * result = focus->addcontext(leftcontext, rightcontext);

    delete leftcontext;
    delete rightcontext;

    return result;
}


int AlignmentModel::extractgizapatterns2(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold, int computereverse, bool bestonly, bool weighbyalignmentscore, bool expandunaligned, bool combine, double unionweight, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {

        /*

        Aligned pattern extraction on the basis of pattern models and GIZA++ Word Alignments

            - alignscorethreshold   - The threshold for the score measuring alignment strength, only patterns passing the threshold can be added
            - bestonly              - Consider only the strongest targe alignment for a given source pattern occurrence
            - weighbyalignmentscore - Also include alignment scores rather than mere frequency in the phrase table scores
            - expandunaligned       - Unaligned words at the end or beginning of a pattern are considered part of the pattern (only sensible without bestonly) //TODO: not implemented yet
            - combine               - Post-processing step, if patterns X and Z are found in a sequence X Y Z where Y contains only unaligned words of max length min(|X|,|Z|), then add X Y Z as a new pattern (provided it passes the threshold)  //TODO: not implemented yet
            - double unionweight    - Value > 1 determining the increase in score per extra union point (set to 1 for intersectiononly)




        basic algorithm: (without expansion and recombination)

        for all patterns s in source sentence:
            for all patterns t in target sentence:
                maxpatternsize = max(|s|,|t|)

                - Is s consistent with t according to the GIZA word alignment?
                - count intersection points, compute intersectionscore: intersectionpoints / maxpatternsize

                if intersectionscore < 1 and not intersectiononly :
                    - count union points - intersectionpoints,  compute unionscore: unionpoints / maxpatternsize

                score = max(intersectionscore + (1 - intersectionscore) * unionscore, 1)

                if not bestonly and score > alignscorethreshold:
                    add pattern with score or with fixed value 1 if not weighbyalignmentscore

            if bestonly and bestscore > alignscorethreshold:
                add best pattern with score or with fixed value 1 if not weighbyalignmentscore

        */



        GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
        GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);

        int found = 0;


        //extract patterns in sentence pair


        //get all patterns in the target sentence
        unordered_set<const EncAnyGram*> * sourcepatterns = &sourcemodel->reverseindex[sentenceindex];
        unordered_set<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentenceindex];

        //reconstruct token-offset in reverse index (making up for not having full index anymore in SelectivePatternModel)
        unordered_map<const EncAnyGram *, vector<int> > sourcetokenfwindex;
        unordered_map<int, vector<const EncAnyGram *> > sourcetokenrevindex;
        recompute_token_index(sourcetokenfwindex, sourcetokenrevindex, sentence_s2t.source, sourcepatterns);
        unordered_map<const EncAnyGram *, vector<int> > targettokenfwindex;
        unordered_map<int, vector<const EncAnyGram *> > targettokenrevindex;
        recompute_token_index(targettokenfwindex, targettokenrevindex, sentence_s2t.target, targetpatterns);


        struct extractedpair {
            const EncAnyGram * source;
            unsigned char sourceindex;
            const EncAnyGram * target;
            unsigned char targetindex;
            int intersectionpoints;
            int unionpoints;
            int unaligned;

            extractedpair(const EncAnyGram *  source, unsigned char sourceindex, const EncAnyGram *  target, unsigned char targetindex, int intersectionpoints, int unionpoints, int unaligned = 0) {
                this->source = source;
                this->sourceindex = sourceindex;
                this->target = target;
                this->targetindex = targetindex;
                this->intersectionpoints = intersectionpoints;
                this->unionpoints = unionpoints;
                this->unaligned = unaligned;
            }
            extractedpair(const extractedpair &other) {
                this->source = other.source;
                this->sourceindex = other.sourceindex;
                this->target = other.target;
                this->targetindex = other.targetindex;
                this->intersectionpoints = other.intersectionpoints;
                this->unionpoints = other.unionpoints;
                this->unaligned = other.unaligned;
            }
            bool operator==(const extractedpair &other) const {
                return ((source == other.source) && (target == other.target) && (sourceindex == other.sourceindex) && (targetindex == other.targetindex));
            }
        };

        vector<extractedpair> extractedpairs;
        //unordered_map<int, vector<extractedpair> > sourcetokenextracted;
        //unordered_map<int, vector<extractedpair> > targettokenextracted;


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

                    count_s++;

                    const EncAnyGram * sourcepattern = *iter_s;
                    unsigned char sourcepatternsize = sourcepattern->n();

                    //grab context if applicable
                    const EncAnyGram * sourcepatternwithcontext = NULL;
                    bool sourcepatternused = false;
                    if (leftsourcecontext || rightsourcecontext) {
                        //extract context
                        sourcepatternwithcontext = addcontext(sentence_s2t.source, sourcepattern, sourceindex);

                        //make sure we re-use any previously used sourcepatternwithcontext rather than using this newly generated one
                        const EncAnyGram * key = getsourcekey(sourcepatternwithcontext, false); // allowfallback=false (do not fall back to sourcemodel, which doesn't know contexts)
                        if (key != NULL) {
                            delete sourcepatternwithcontext; //delete generated one
                            sourcepatternwithcontext = key; //use existing key
                            //cout << "CONTEXT EXISTS @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        //} else {
                            //cout << "CONTEXT NEW @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        }
                    }


                    //now find what target patterns are aligned, and how well the alignment is (expressed through a score)
                    double bestscore = 0;
                    const EncAnyGram * besttargetpattern = NULL;
                    int besttargetindex = 0;
                    int bestintersectionpoints = 0;
                    int bestunionpoints = 0;
                    int count_t = 0;

                    multiset<uint32_t> * sourcesentenceindex = NULL; //used only if coocthreshold > 0
                    if (coocthreshold > 0) {
                        if (sourcepattern->isskipgram()) {
                            sourcesentenceindex = &sourcemodel->skipgrams[*( (EncSkipGram*) sourcepattern)].sentences;
                        } else {
                            sourcesentenceindex = &sourcemodel->ngrams[*( (EncNGram*) sourcepattern)].sentences;
                        }
                    }

                    for (unordered_set<const EncAnyGram*>::iterator iter_t = targetpatterns->begin(); iter_t != targetpatterns->end(); iter_t++) { //iterate over all target patterns IN THIS SENTENCE
                         count_t++;
                         const EncAnyGram * targetpattern = *iter_t;
                         if (DEBUG) cerr << " @" << sentenceindex << "-" << count_s << "-" << count_t << endl;
                         if ((DEBUG) && (sourcedecoder != NULL) && (targetdecoder != NULL)) cerr << "   DEBUG: [" << sourcepattern->decode(*sourcedecoder) << "] vs [" << targetpattern->decode(*targetdecoder) << "]" << endl;

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

                                 int intersectionpoints = 0;
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
                                        intersectionpoints++;
                                        if (alignedsourceindex == sourceindex) firstsourcealigned = true;
                                        if (alignedsourceindex == sourceindex + (sourcepatternsize - 1) ) lastsourcealigned = true;
                                        if (alignedtargetindex == targetindex) firsttargetaligned = true;
                                        if (alignedtargetindex == targetindex + (targetpatternsize - 1) ) lasttargetaligned = true;
                                     } else if ((insource) || (intarget)) {
                                        //alignment point outside pattern, invalid alignment, discard
                                        intersectionpoints = 0; break;
                                     }
                                 }
                                 if (DEBUG) cerr << "     DEBUG INTERSECTIONPOINTS=" << intersectionpoints << "  ENDS: " << (int) firstsourcealigned << " " << (int) lastsourcealigned << " " << (int) firsttargetaligned << " " << (int) lasttargetaligned << endl;
                                 if ((intersectionpoints == 0) || (!firstsourcealigned)  || (!lastsourcealigned) || (!firsttargetaligned)  || (!lasttargetaligned)) break;
                                 //if (DEBUG) cerr << "     DEBUG FOUND" << endl;


                                 const double intersectionscore = (double) intersectionpoints / maxpatternsize; //1.0 if all points are aligned

                                 if (DEBUG) cerr << "     DEBUG intersectionscore(1): " << intersectionpoints << " / " << (int) targetpatternsize << " = " << intersectionscore << endl;



                                 int unionpoints = 0;
                                 if ((intersectionscore < 1) && (unionweight != 1)) {
                                     //check alignment points in union

                                     for (multimap<const unsigned char, const unsigned char>::const_iterator alignmentiter = sentence_u.alignment.begin(); alignmentiter != sentence_u.alignment.end(); alignmentiter++) {
                                        const unsigned char alignedsourceindex = alignmentiter->first -1; //-1 because of 1-based-indexing
                                        const unsigned char alignedtargetindex = alignmentiter->second -1; //-1 because of 1-based-indexing

                                        const bool insource = (alignedsourceindex >= sourceindex) && (alignedsourceindex < sourceindex + sourcepattern->n());
                                        const bool intarget = (alignedtargetindex >= targetindex) && (alignedtargetindex < targetindex + targetpattern->n());

                                         if ((insource) && (intarget)) {
                                            unionpoints++;
                                         }
                                     }
                                     unionpoints = unionpoints - intersectionpoints; //substract intersection points, we want only points in the union and not in the intersection
                                     if (DEBUG) cerr << "     DEBUG UNIONPOINTS=" << unionpoints << endl;
                                     /*if (unionpoints > 0) {
                                        unionscore =  (double) unionpoints / maxpatternsize;

                                        if (DEBUG) cerr << "     DEBUG unionscore(2): " << unionpoints << " / " << (int) maxpatternsize << " = " << unionscore << endl;
                                     }*/
                                }

                                double score = intersectionscore;
                                if (unionweight != 1) score = score * pow(unionweight,unionpoints);
                                if (score > 1) score = 1;
                                if (DEBUG) cerr << "     DEBUG score(2): " << intersectionscore << " * " << unionweight << "^" << unionpoints << " (max 1) = " << score << endl;

                                if (bestonly) {
                                    if (score > bestscore) {
                                        //retain only the best target pattern given an occurrence of a source pattern
                                        bestscore = score;
                                        besttargetpattern = targetpattern;
                                        besttargetindex = targetindex;
                                        bestintersectionpoints = intersectionpoints;
                                        bestunionpoints = unionpoints;
                                    }
                                } else if (score >= alignscorethreshold) {
                                    addextractedpattern(sourcepattern, targetpattern, (weighbyalignmentscore ? score : 1), computereverse, sourcepatternwithcontext);
                                    if (combine) {
                                        extractedpairs.push_back(extractedpair(sourcepattern, sourceindex, targetpattern, targetindex, intersectionpoints, unionpoints));
                                    }
                                    sourcepatternused = true;
                                    found++;
                                    if ((sourcedecoder != NULL) && (targetdecoder != NULL)) {
                                        cout << sourcepattern->decode(*sourcedecoder) << " ||| " << targetpattern->decode(*targetdecoder) << " ||| " << score << endl;
                                    }

                                }

                                if ((!bestonly) && (expandunaligned)) {
                                    //try and expand the pattern with unaligned points to the target pattern, on both ends

                                    //target indices
                                    int lastleftaligned = -1;
                                    int firstrightaligned = sentence_s2t.target->size();

                                    //(using union instead of intersection, strictest sense of what it means to be 'unaligned')
                                    for (multimap<const unsigned char, const unsigned char>::const_iterator alignmentiter = sentence_u.alignment.begin(); alignmentiter != sentence_u.alignment.end(); alignmentiter++) {
                                        const unsigned char alignedtargetindex = alignmentiter->second -1; //-1 because of 1-based-indexing
                                        if ((alignedtargetindex > lastleftaligned) && (alignedtargetindex < targetindex)) {
                                            lastleftaligned = alignedtargetindex;
                                        }
                                        if ((alignedtargetindex < firstrightaligned) && (alignedtargetindex > targetindex + (targetpatternsize - 1))) {
                                            firstrightaligned = alignedtargetindex;
                                        }
                                    }
                                    if (DEBUG) cerr << "     DEBUG lastleftaligned=" << lastleftaligned << " firstrightaligned=" << firstrightaligned << " targetindex=" << targetindex << " +targetpatternsize-1=" << targetindex + (targetpatternsize - 1) << endl;

                                    for (int begin = lastleftaligned+1; begin <= targetindex; begin++) {
                                        for (int end = targetindex + targetpatternsize; end <= firstrightaligned-1; end++) {
                                            if (!((begin == targetindex) && (end == targetindex + (targetpatternsize - 1)))) { //exclude the non-expanded pattern we already have anyway
                                                if (DEBUG) cerr << "     DEBUG Found expansion using unaligned points from " << begin << " to " << end << endl;
                                                const int length = end-begin + 1;
                                                EncNGram * newtargetpattern = sentence_s2t.target->slice(begin, length);
                                                const EncAnyGram * targetkey = gettargetkey(newtargetpattern);
                                                if (targetkey != NULL) { //check if not found in pattern model (most will be new patterns)
                                                    //found in pattern model anyway, use that one
                                                    delete newtargetpattern;
                                                    newtargetpattern = (EncNGram *) targetkey;
                                                }
                                                const unsigned char maxpatternsize2 = ( length > maxpatternsize) ? length : maxpatternsize;
                                                const int unaligned = length - targetpatternsize;

                                                double score2 = (double) intersectionpoints / maxpatternsize2;
                                                if ((intersectionscore < 1) && (unionweight != 1)) {
                                                    score2 = score2 * pow(unionweight,unionpoints);
                                                }
                                                if (score2 > 1) score2 = 1;

                                                //substract penalty for unaligned points. The size of the penalty depends on the length of the SOURCE patterns, the longer the source, the lower the penalty for unaligned points, the shorter the source, the higher the penalty for adding unaligned points (up to -50% per unaligned point).. score2 drops to 0 quite fast if more points are added
                                                const double penaltyweight = (double) 1 / (sourcepatternsize+1);
                                                score2 = score2 * (1 - (unaligned * penaltyweight));
                                                if (DEBUG) cerr << "     DEBUG score2=" << score2 << " penaltyweight=" << penaltyweight << endl;
                                                if ((score2 >= alignscorethreshold) && (score2 > 0)) {
                                                    addextractedpattern(sourcepattern, newtargetpattern, (weighbyalignmentscore ? score2 : 1), computereverse, sourcepatternwithcontext);
                                                    sourcepatternused = true;
                                                    found++;
                                                } else if (targetkey == NULL) {
                                                    delete newtargetpattern;
                                                }

                                            }
                                        }
                                    }


                                }

                            }
                         }
                    } //iteration over all target patterns

                    if ((besttargetpattern != NULL) && (bestscore >= alignscorethreshold)) {
                        addextractedpattern(sourcepattern, besttargetpattern, (weighbyalignmentscore ? bestscore : 1), computereverse, sourcepatternwithcontext);
                        if (combine) {
                            extractedpairs.push_back(extractedpair(sourcepattern, sourceindex, besttargetpattern, besttargetindex, bestintersectionpoints, bestunionpoints));
                        }

                        sourcepatternused = true;
                        found++;
                        if ((sourcedecoder != NULL) && (targetdecoder != NULL)) {
                            cout << sourcepattern->decode(*sourcedecoder) << " ||| " << besttargetpattern->decode(*targetdecoder) << " ||| " << bestscore << endl;
                        }
                    }

                    if (sourcepatternwithcontext != NULL) {
                        if (sourcepatternused) {
                            sourcecontexts[sourcepattern].insert(sourcepatternwithcontext);
                        } else {
                            //delete sourcepatternwithcontext; //TODO: Reenable: causes segfault... I think we already deleted this earlier, then this is not needed
                        }
                    }
                }
             }
          }  //sourceindex iterator


         //------ RECOMBINATION STEP -----------
         // if two patterns can be combined into a larger one, on both source and target side, then add the larger one too
         //      patterns may be combined if they follow eachother and have none or only unaligned points in between (for both sides)  (swap order is allowed)
         //      the maximum size of a gap between patterns is limited by the smallest pattern


         if ((combine) && (leftsourcecontext || rightsourcecontext)) {
            // TODO: context not yet supported with recombination!
            cerr << "ERROR: Context and recombination can not be mixed (not implemented yet)" << endl;
            throw InternalError();
         }

         if (combine) {
              map<int,bool> sourcecoverage;
              map<int,bool> targetcoverage;
              for (vector<extractedpair>::iterator iter = extractedpairs.begin(); iter != extractedpairs.end(); iter++) {
                 for (int i = iter->sourceindex; i < iter->sourceindex + iter->source->n(); i++) sourcecoverage[i] = true;
                 for (int i = iter->targetindex; i < iter->targetindex + iter->source->n(); i++) targetcoverage[i] = true;
              }
              for (multimap<const unsigned char, const unsigned char>::const_iterator alignmentiter = sentence_u.alignment.begin(); alignmentiter != sentence_u.alignment.end(); alignmentiter++) {
                const unsigned char alignedsourceindex = alignmentiter->second -1; //-1 because of 1-based-indexing
                const unsigned char alignedtargetindex = alignmentiter->second -1; //-1 because of 1-based-indexing
                sourcecoverage[(int) alignedsourceindex] = true;
                targetcoverage[(int) alignedtargetindex] = true;
              }

              for (vector<extractedpair>::iterator iter = extractedpairs.begin(); iter != extractedpairs.end(); iter++) {
                    const int n1 = iter->source->n();
                    for (vector<extractedpair>::iterator iter2 = extractedpairs.begin(); iter2 != extractedpairs.end(); iter2++) {
                        if (iter2->sourceindex >= iter->sourceindex + n1) { //does iter2 occur after iter1?
                            const int n2 = iter2->source->n();
                            const int minsize = (n1 > n2) ? n1 : n2;
                            const int gapsize = iter2->sourceindex - (iter->sourceindex + n1);
                            if (gapsize > minsize) {
                                //gap is too big, skip
                                continue;
                            } else if (gapsize > 0) {
                                //if there is a gap between iter and iter2, make sure it is not already covered by a third extraction OR an alignment
                                bool gapalreadycovered = false;
                                for (int i = iter->sourceindex + n1; i < iter2->sourceindex; i++) {
                                    if (sourcecoverage.count(i)) {
                                        gapalreadycovered = true;
                                        break;
                                    }
                                }
                                if (gapalreadycovered) continue;
                            }

                            const unsigned char s_begin = iter->sourceindex;
                            const unsigned char s_length = (iter2->sourceindex + n2) - s_begin;
                            const EncAnyGram * recombinedsource = sentence_s2t.source->slice(s_begin, s_length);
                            const EncAnyGram * recombinedsourcekey = getsourcekey(recombinedsource);

                            if (recombinedsourcekey != NULL) {
                                delete recombinedsource;
                                recombinedsource = recombinedsourcekey;

                                //check if the pattern was already extracted:
                                double duplicate = false;
                                for (vector<extractedpair>::iterator iter3 = extractedpairs.begin(); iter3 != extractedpairs.end(); iter3++) {
                                    if (iter3->source->hash() == iter->source->hash()) {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (duplicate) {
                                    if (DEBUG) cerr << "     DEBUG recombination yields source that has already been extracted, skipping" << endl;
                                    continue; //if so, skip
                                }
                            }



                            //ok, we have two good source patterns, not extracted yet, now check if the target side holds up too
                            const int t_n1 = iter->target->n();
                            const int t_n2 = iter2->target->n();
                            const int t_minsize = (t_n1 > t_n2) ? t_n1 : t_n2;
                            unsigned char t_begin = 0;
                            unsigned char t_length = 0;
                            int t_gapsize = 0;
                            if ((iter2->targetindex >= iter->targetindex + t_n1)) {
                                //t_iter2 after t_iter
                                t_gapsize = iter2->targetindex - (iter->targetindex + t_n1);
                                bool gapalreadycovered = false;
                                for (int i = iter->targetindex + n1; i < iter2->targetindex; i++) {
                                    if (targetcoverage.count(i)) {
                                        gapalreadycovered = true;
                                        break;
                                    }
                                }
                                if (gapalreadycovered) continue;
                                t_begin = iter->targetindex;
                                t_length = (iter2->targetindex + t_n2) - t_begin;
                            } else if (iter->targetindex >= iter2->targetindex + t_n2) {
                                //t_iter after t_iter2 (swapped order)
                                t_gapsize = iter->targetindex - (iter2->targetindex + t_n2);
                                bool gapalreadycovered = false;
                                for (int i = iter2->targetindex + n2; i < iter->targetindex; i++) {
                                    if (targetcoverage.count(i)) {
                                        gapalreadycovered = true;
                                        break;
                                    }
                                }
                                if (gapalreadycovered) continue;
                                t_begin = iter2->targetindex;
                                t_length = (iter->targetindex + t_n1) - t_begin;
                            }

                            if (t_length > 0) {
                                //we have a recombination!
                                const EncAnyGram * recombinedtarget = sentence_s2t.target->slice(t_begin, t_length);
                                const unsigned char maxpatternsize2 = ( s_length > t_length) ? s_length : t_length;

                                //check if it happens to already occur as pattern
                                const EncAnyGram * recombinedtargetkey = gettargetkey(recombinedtarget);
                                if (recombinedtargetkey != NULL) {
                                    delete recombinedtarget;
                                    recombinedtarget = recombinedtargetkey;
                                }

                                //compute new score
                                double score2 = (double) (iter->intersectionpoints + iter2->intersectionpoints) / maxpatternsize2;
                                if ((score2 < 1) && (unionweight != 1)) {
                                    score2 = score2 * pow(unionweight,iter->unionpoints + iter2->unionpoints);
                                }
                                if (score2 > 1) score2 = 1;

                                //substract penalty for unaligned points. The size of the penalty depends on the length of the SOURCE patterns, the longer the source, the lower the penalty for unaligned points, the shorter the source, the higher the penalty for adding unaligned points (up to -50% per unaligned point).. score2 drops to 0 quite fast if more points are added
                                const int unaligned = iter->unaligned + iter2->unaligned + t_gapsize;
                                const double penaltyweight = (double) 1 / (s_length+1);
                                score2 = score2 * (1 - (unaligned * penaltyweight));
                                if (DEBUG) cerr << "     DEBUG recombination found, score=" << score2 << " penaltyweight=" << penaltyweight << endl;
                                if ((score2 >= alignscorethreshold) && (score2 > 0)) {
                                    addextractedpattern(recombinedsource, recombinedtarget, (weighbyalignmentscore ? score2 : 1), computereverse);
                                    found++;
                                } else {
                                    if (recombinedsourcekey == NULL) delete recombinedsource;
                                    if (recombinedtargetkey == NULL) delete recombinedtarget;
                                }


                            }
                        }
                    }
              }

         }

         return found;
}


void AlignmentModel::addextractedpattern(const EncAnyGram * sourcepattern, const EncAnyGram * targetpattern, double score, int computereverse, const EncAnyGram * sourcepatternwithcontext) {

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
            if (alignmatrix[key][targetpattern].empty()) {
                alignmatrix[key][targetpattern].push_back(score);
            } else {
                alignmatrix[key][targetpattern][0] +=score;
            }

            if (computereverse) {
                if (reversealignmatrix[targetpattern][key].empty()) {
                    reversealignmatrix[targetpattern][key].push_back(score);
                } else {
                    reversealignmatrix[targetpattern][key][0] += score;
                }
            }


}


int AlignmentModel::extractgizapatterns(GizaSentenceAlignment & sentence_s2t, GizaSentenceAlignment & sentence_t2s, int sentenceindex, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold, int computereverse, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
        //SUPERVISED
        GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
        GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);

        int found = 0;


        //extract patterns in sentence pair

        /*


           2)  (OLD IDEA)
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

                       maxpatternsize = max(|pattern_s|,|pattern_t|)


                       if score > 0
                            bestscore = score
                            bestpattern_t = pattern_t


        */

        //get all patterns in the target sentence
        unordered_set<const EncAnyGram*> * sourcepatterns = &sourcemodel->reverseindex[sentenceindex];
        unordered_set<const EncAnyGram*> * targetpatterns = &targetmodel->reverseindex[sentenceindex];

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

                    count_s++;

                    const EncAnyGram * sourcepattern = *iter_s;
                    unsigned char sourcepatternsize = sourcepattern->n();

                    //grab context if applicable
                    const EncAnyGram * sourcepatternwithcontext = NULL;
                    bool sourcepatternused = false;
                    if (leftsourcecontext || rightsourcecontext) {
                        //extract context
                        sourcepatternwithcontext = addcontext(sentence_s2t.source, sourcepattern, sourceindex);

                        //make sure we re-use any previously used sourcepatternwithcontext rather than using this newly generated one
                        const EncAnyGram * key = getsourcekey(sourcepatternwithcontext, false); // allowfallback=false (do not fall back to sourcemodel, which doesn't know contexts)
                        if (key != NULL) {
                            delete sourcepatternwithcontext; //delete generated one
                            sourcepatternwithcontext = key; //use existing key
                            //cout << "CONTEXT EXISTS @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        //} else {
                            //cout << "CONTEXT NEW @" << (size_t) sourcepatternwithcontext << " #" << sourcepatternwithcontext->hash() << endl;
                        }
                    }


                    //now find what target patterns are aligned, and how well the alignment is (expressed through a score)
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

                    for (unordered_set<const EncAnyGram*>::iterator iter_t = targetpatterns->begin(); iter_t != targetpatterns->end(); iter_t++) { //iterate over all target patterns IN THIS SENTENCE
                         count_t++;
                         const EncAnyGram * targetpattern = *iter_t;
                         if (DEBUG) cerr << " @" << sentenceindex << "-" << count_s << "-" << count_t << endl;
                         if ((DEBUG) && (sourcedecoder != NULL) && (targetdecoder != NULL)) cerr << "   DEBUG: [" << sourcepattern->decode(*sourcedecoder) << "] vs [" << targetpattern->decode(*targetdecoder) << "]" << endl;

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
                                        //alignment point outside pattern, invalid alignment, discard
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
                        addextractedpattern(sourcepattern, besttargetpattern, 1, computereverse, sourcepatternwithcontext);
                        sourcepatternused = true;
                        found++;
                        if ((sourcedecoder != NULL) && (targetdecoder != NULL)) {
                            cout << sourcepattern->decode(*sourcedecoder) << " ||| " << besttargetpattern->decode(*targetdecoder) << " ||| " << bestscore << endl;
                        }
                    }

                    if (sourcepatternwithcontext != NULL) {
                        if (sourcepatternused) {
                            sourcecontexts[sourcepattern].insert(sourcepatternwithcontext);
                        } else {
                            //delete sourcepatternwithcontext; //TODO: Reenable: causes segfault... I think we already deleted this earlier, then this is not needed
                        }
                    }
                }
             }
          }  //sourceindex iterator
         return found;
}



int AlignmentModel::extractgizapatterns2(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold, int computereverse, bool bestonly, bool weighbyalignmentscore, bool expandunaligned, bool combine, double unionweight, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
    unsigned int totalfound = 0;
    while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {
        GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
        GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();
        const int i = gizamodel_s2t.index();
        cerr << " @" << i << endl;
        int found = extractgizapatterns2(sentence_s2t, sentence_t2s, i, pairoccurrencethreshold, coocthreshold, alignscorethreshold, computereverse, bestonly, weighbyalignmentscore, expandunaligned, combine, unionweight, sourcedecoder,targetdecoder);
        totalfound += found;
        cerr << "   found " << found << endl;
      } //alignment read


      if (pairoccurrencethreshold > 1) {
        cerr << "Pruning according to pair occurrence threshold (" << pairoccurrencethreshold << ")" << endl;
        prune(pairoccurrencethreshold);
      }

      //normalize alignment matrix
      cerr << "Normalising" << endl;
      normalize();
      if (computereverse) {
            cerr << "Integrating reverse model" << endl;
            integratereverse(computereverse);
      }

      return totalfound;
}


int AlignmentModel::extractgizapatterns(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, int pairoccurrencethreshold, const double coocthreshold, const double alignscorethreshold, int computereverse, ClassDecoder * sourcedecoder, ClassDecoder * targetdecoder) {
    unsigned int totalfound = 0;
    while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {
        GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
        GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();
        const int i = gizamodel_s2t.index();
        cerr << " @" << i << endl;
        int found = extractgizapatterns(sentence_s2t, sentence_t2s, i, pairoccurrencethreshold, coocthreshold, alignscorethreshold, computereverse, sourcedecoder,targetdecoder);
        totalfound += found;
        cerr << "   found " << found << endl;
      } //alignment read


      if (pairoccurrencethreshold > 0) prune(pairoccurrencethreshold);

      //normalize alignment matrix
      normalize();
      if (computereverse) {
            cerr << "Integrating reverse model" << endl;
            integratereverse(computereverse);
      }

      return totalfound;
}



void AlignmentModel::integratereverse(int computereverse) {
    normalize(&reversealignmatrix);

    //add reversealignmodel to normal model and delete reversemodel
    for (t_alignmatrix::iterator iter = reversealignmatrix.begin(); iter != reversealignmatrix.end(); iter++) {
        const EncAnyGram * target = iter->first;
        for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram * source = iter2->first;
            if (iter2->second.size() != 1) {
                cerr << "INTERNAL ERROR: Expected one score for pair in reverse matrix, got " << iter2->second.size() << endl;
                throw InternalError();
            }
            if ( alignmatrix.count(source) == 0) {
                cerr << "INTERNAL ERROR: [integratereverse] Source pattern not found, this should not happen" << endl;
                throw InternalError();
            }
            if ( alignmatrix[source].count(target) == 0) {
                cerr << "INTERNAL ERROR: [integratereverse] Target pattern not found, this should not happen" << endl;
                throw InternalError();
            }
            if ( alignmatrix[source][target].size() != 1) {
                cerr << "INTERNAL ERROR: Expected one score for pair in alignmatrix, got " <<alignmatrix[source][target].size() << endl;
                throw InternalError();
            }
            if (computereverse == 1) {
                alignmatrix[source][target][0] = alignmatrix[source][target][0] * iter2->second[0];
            } else if (computereverse == 2) {
                alignmatrix[source][target].push_back(iter2->second[0]);
            }
        }
        //reversealignmatrix.erase(target); //make space whilst processing
    }
}


int AlignmentModel::extractgizapatterns_heur(GizaSentenceAlignment & sentence_a, int sentenceindex, int computereverse) {
    //extracts phrases based on alignment point dataset



    //From statmt.org: We collect all aligned phrase pairs that are consistent with the word alignment: The words in a legal phrase pair are only aligned to each other, and not to words outside.
    // See also Statistical Machine Translation (Philip Koehn), pag 133

    int found = 0;

    const int t_length = sentence_a.target->length();
    const int s_length = sentence_a.source->length();


    //phrases will be stored here ((s_start,s_end),(t_start,t_end))
    vector< pair< pair<int,int>,pair<int,int> > > phrases;

    for (int t_start = 0; t_start < t_length; t_start++) {
        for (int t_end = t_start; t_end < t_length; t_end++) {


            //find the minimally matching source phrase given t_tstart,t_end
            int s_start = s_length;
            int s_end = 0;
            for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_a.alignment.begin(); iter != sentence_a.alignment.end(); iter++) {
                const unsigned char s = iter->first;
                const unsigned char t = iter->second;
                if ((t_start <= t) && (t <= t_end)) {
                    if (s < s_start) s_start = s;
                    if (s > s_end) s_end = s;
                }

            }

            if ((s_end == 0) || (s_start > s_end)) continue; //nothing found

            //extract phrase


            //check if alignment points violate consistency
            bool valid = true;
            for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_a.alignment.begin(); iter != sentence_a.alignment.end(); iter++) {
                const unsigned char s = iter->first;
                const unsigned char t = iter->second;
                if ( (s_start <= s) && (s_end >= s) && ((t < t_start) || (t > t_end))  ) {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue; //not valid!

            //good, valid

            //add phrase pairs (incl. additional unaligned s; unaligned words on edges)
            bool s_startaligned = false;
            bool s_endaligned = false;

            int s_start2 = s_start;
            do {
                int s_end2 = s_end;

                do {
                    //add phrasepair
                    if (DEBUG) cerr << "\t[extractgizapatterns_heur] found phrase (" << s_start2 << "," << s_end2 << ") -> (" << t_start << "," << t_end << ")" << endl;
                    phrases.push_back(pair< pair<int,int>,pair<int,int> >( pair<int,int>(s_start2,s_end2), pair<int,int>(t_start,t_end) ) );

                    //attempt to add free unaligned points (widen source phrase to the right):
                    s_end2++;
                    if (s_end2 >= s_length) break;
                    //is s_end2 aligned?
                    s_endaligned= (sentence_a.alignment.count((unsigned char) s_end2) > 0);
                } while (!s_endaligned);

                //attempt to add free unaligned points (widen source phrase to the left):
                s_start2--;
                if (s_start2 < 0) break;
                s_startaligned= (sentence_a.alignment.count((unsigned char) s_start2) > 0);
            } while (!s_startaligned);



        }
    }




    //add actual phrases
    for (vector< pair< pair<int,int>,pair<int,int> > >::iterator iter = phrases.begin(); iter != phrases.end(); iter++) {
        int s_start = iter->first.first;
        int s_end = iter->first.second;
        int t_start = iter->first.first;
        int t_end = iter->first.second;


        EncNGram * sourcegramnocontext = sentence_a.source->slice(s_start,s_end - s_start + 1);
        EncNGram * sourcegram;

        bool context = true;
        s_start = s_start - leftsourcecontext;
        s_end = s_end + rightsourcecontext;
        if ((s_start < 0) && (s_end > s_length)) {
            //add begin and end markers
            EncNGram bos = EncNGram(&BOSCLASS, 1);
            EncNGram eos = EncNGram(&EOSCLASS, 1);
            EncNGram tmp2 = EncNGram(bos + *sourcegramnocontext);
            sourcegram = new EncNGram(tmp2 + eos);
        } else if (s_start < 0) {
            //add begin marker
            EncNGram bos = EncNGram(&BOSCLASS, 1);
            sourcegram = new EncNGram(bos + *sourcegramnocontext);
        } else if (s_end > s_length) {
            //add end marker
            EncNGram eos = EncNGram(&EOSCLASS, 1);
            sourcegram = new EncNGram(*sourcegramnocontext + eos);
        } else {
            sourcegram = sourcegramnocontext;
            context = false;
        }


        const EncAnyGram * sourcegramkey = getsourcekey(sourcegram);
        if (sourcegramkey != NULL) {
            if (DEBUG) cerr << "\t[extractgizapatterns_heur] sourcegram has been used earlier, reusing" << endl;
            delete sourcegram;
            sourcegram = (EncNGram *) sourcegramkey;
        }


        if (context) {
            if (DEBUG) cerr << "\t[extractgizapatterns_heur] context" << endl;
            sourcecontexts[(const EncAnyGram *) sourcegramnocontext].insert((const EncAnyGram *) sourcegram);
        }


        const EncAnyGram * targetgram = (const EncAnyGram *) sentence_a.target->slice(t_start,t_end - t_start + 1);
        const EncAnyGram * targetgramkey = gettargetkey(targetgram);
        if (targetgramkey == NULL) {
            if (targetgram->isskipgram()) {
                //ok, never happens, but perhaps for future
                cerr << "targetgram is skipgram? not possible!" << endl;
                throw InternalError();
            } else {
                if (DEBUG) cerr << "\t[extractgizapatterns_heur] targetgram is new, adding" << endl;
                pair<unordered_set<EncNGram>::iterator, bool> returnvalue;
                returnvalue = targetngrams.insert(  *( (const EncNGram *) targetgram) );
                cerr << "\t[extractgizapatterns_heur] h=" << targetgram->hash() << " gapcount=" << (int) targetgram->gapcount() << endl;
                targetgram = gettargetkey(targetgram);
                if (targetgram == NULL) {
                    cerr << "\t[extractgizapatterns_heur] targetgram == NULL _after_ insertion! should not happen!" << endl;
                    const EncNGram * t = &(*(returnvalue.first));
                    cerr << "\t[extractgizapatterns_heur] h2=" << t->hash() << endl;
                    cerr << targetgram << " vs " << t << endl;
                    throw InternalError();
                }
            }
        } else {
            if (DEBUG) cerr << "\t[extractgizapatterns_heur] targetgram has been used earlier, reusing" << endl;
            delete targetgram;
            targetgram = targetgramkey;
        }

        //sanity check
        sourcegram->hash();
        targetgram->hash();

        if (DEBUG) cerr << "\t[extractgizapatterns_heur] adding alignment" << endl;

        //add alignment
        if (alignmatrix[(const EncAnyGram *)sourcegram][targetgram].empty()) {
            alignmatrix[(const EncAnyGram *)sourcegram][targetgram].push_back(1);
            found++;
        } else {
            alignmatrix[(const EncAnyGram *)sourcegram][targetgram][0] += 1;
        }

        if (computereverse) {
            if (reversealignmatrix[targetgram][(const EncAnyGram *)sourcegram].empty()) {
                reversealignmatrix[targetgram][(const EncAnyGram *)sourcegram].push_back(1);
            } else {
                reversealignmatrix[targetgram][(const EncAnyGram *)sourcegram][0] += 1;
            }
        }
    }

    return found;

}

int AlignmentModel::extractgizapatterns_heur(GizaModel & gizamodel_s2t, GizaModel & gizamodel_t2s, PhraseAlignHeuristic phrasealignheuristic, int computereverse) {
    unsigned int totalfound = 0;
    while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {
        GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
        GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();
        const int i = gizamodel_s2t.index();
        cerr << " @" << i << endl;


        int found = 0;
        if ((phrasealignheuristic == PAH_GROWDIAG) || (phrasealignheuristic == PAH_GROWDIAGFINAL)) {
            GizaSentenceAlignment sentence_a = extractgiza_growdiag(sentence_s2t, sentence_t2s);
            if (phrasealignheuristic == PAH_GROWDIAGFINAL) {
                extractgiza_final(sentence_a, sentence_s2t, sentence_t2s);
            }
            found = extractgizapatterns_heur(sentence_a, i, computereverse);
        } else if (phrasealignheuristic == PAH_S2T) {
            found = extractgizapatterns_heur(sentence_s2t, i, computereverse );
        } else if (phrasealignheuristic == PAH_INTERSECTION) {
            GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
            found = extractgizapatterns_heur(sentence_i, i, computereverse);
        } else if (phrasealignheuristic == PAH_UNION) {
            GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);
            found = extractgizapatterns_heur(sentence_u, i, computereverse);
        }

        //int found = extractgizapatterns(sentence_s2t, sentence_t2s, i, pairoccurrencethreshold, coocthreshold, alignscorethreshold, computereverse, sourcedecoder,targetdecoder);
        totalfound += found;
        cerr << "   found " << found << endl;
      } //alignment read


      normalize();
      if (computereverse) {
            cerr << "Integrating reverse model" << endl;
            integratereverse(computereverse);
      }

      return totalfound;
}




GizaSentenceAlignment AlignmentModel::extractgiza_growdiag(GizaSentenceAlignment & sentence_s2t ,GizaSentenceAlignment & sentence_t2s) {
    /*
     GROW-DIAG():
      iterate until no new points added
        for english word e = 0 ... en
          for foreign word f = 0 ... fn
            if ( e aligned with f )
              for each neighboring point ( e-new, f-new ):
                if ( ( e-new not aligned or f-new not aligned ) and
                     ( e-new, f-new ) in union( e2f, f2e ) )
                  add alignment point ( e-new, f-new )
    */

    int added;
    GizaSentenceAlignment sentence_a = sentence_s2t.intersect(sentence_t2s); //alignment starts with intersection
    GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);



    do {
        if (DEBUG) cerr << "\t[growdiag] Next iteration" << endl;
        added = 0;

        for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_a.alignment.begin(); iter != sentence_a.alignment.end(); iter++) {
            //for each source word e
            //for each target word f
            //where e aligned with f

            //for each neighbouring point
            for (int x = iter->first - 1; x <= iter->first + 1; x++) {
                for (int y = iter->second - 1; y <= iter->second + 1; y++) {
                    if ((x < 0) || (x >= sentence_s2t.source->length())) continue; //impossible point
                    if ((y < 0) || (y >= sentence_s2t.target->length())) continue; //impossible point

                    if (!(( x == (int) iter->first) && (y == (int) iter->second))) { //if not the same point (we only want neigbours)
                        //check if e-new (y) not aligned or f-new (x) not aligned


                        //is this neighbour aligned already?
                        if (sentence_a.alignment.count( (unsigned char) x)) { //check x
                            continue; //yes, aligned already, break
                        } else {
                            //check y
                            bool found = false;
                            for (multimap<const unsigned char,const unsigned char>::iterator iter2 = sentence_a.alignment.lower_bound((unsigned char)  x); iter2 != sentence_a.alignment.upper_bound((unsigned char)  x); iter2++) {
                                if (iter2->second == (unsigned char)  y) {
                                    found = true;
                                    break;
                                }
                            }
                            if (found) {
                                continue; //yes, aligned already, break
                            } else {
                                //not aligned yet.. is it in the union?
                                if (sentence_u.alignment.count( (unsigned char) x) == 0) { //check x
                                    continue; //no, therefore no candidate, break
                                } else {
                                    //check y in union
                                    found = false;
                                    for (multimap<const unsigned char,const unsigned char>::iterator iter2 = sentence_u.alignment.lower_bound((unsigned char)  x); iter2 != sentence_u.alignment.upper_bound((unsigned char)  x); iter2++) {
                                        if (iter2->second == (unsigned char)  y) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if (found) {
                                        //yes, (x,y) found in union, add alignment point:
                                        sentence_a.alignment.insert(pair<const unsigned char, const unsigned char>((unsigned char)  x, (unsigned char) y ));
                                        added += 1;
                                        if (DEBUG) cerr << "\t[growdiag] added alignment point (" << x << "," << y << ")" << endl;
                                    } else {
                                        //no, break
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
        if (DEBUG) cerr << "\t[growdiag] " << added << " alignment points added this iteration" << endl;
    } while (added > 0);
    return sentence_a;
}

void AlignmentModel::extractgiza_final(GizaSentenceAlignment & sentence_a ,GizaSentenceAlignment & sentence_s2t , GizaSentenceAlignment & sentence_t2s ) {
    /*FINAL(a):
        for english word e-new = 0 ... en
            for foreign word f-new = 0 ... fn
              if ( ( e-new not aligned or f-new not aligned ) and
                   ( e-new, f-new ) in alignment a )
                add alignment point ( e-new, f-new )
    */

    //does final step both ways


    //temporary map to easily find targets
    set<unsigned char> targets_a;
    for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_a.alignment.begin(); iter != sentence_a.alignment.end(); iter++) {
        targets_a.insert(iter->second);
    }


    for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_s2t.alignment.begin(); iter != sentence_s2t.alignment.end(); iter++) {
        const unsigned char s = iter->first;
        const unsigned char t = iter->second;

        const bool s_aligned = sentence_a.alignment.count(s);
        const bool t_aligned = targets_a.count(t);
        if (!s_aligned || !t_aligned) {
            sentence_a.alignment.insert(pair<const unsigned char, const unsigned char>(s,t));
            if (DEBUG) cerr << "\t[final] added alignment point (" << (int) s << "," << (int) t << ")" << endl;
            targets_a.insert(t);
        }
    }


    //inverse direction:
    for (multimap<const unsigned char,const unsigned char>::iterator iter = sentence_t2s.alignment.begin(); iter != sentence_t2s.alignment.end(); iter++) {
        const unsigned char s = iter->second;
        const unsigned char t = iter->first;

        const bool t_aligned = targets_a.count(t);;
        const bool s_aligned = sentence_a.alignment.count(s);
        if (!s_aligned || !t_aligned) {
            sentence_a.alignment.insert(pair<const unsigned char, const unsigned char>(s,t));
            if (DEBUG) cerr << "\t[final inverse] added alignment point (" << (int) s << "," << (int) t << ")" << endl;
        }
    }
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


void recompute_token_index(unordered_map<const EncAnyGram *, vector<int> > & tokenfwindex, unordered_map<int, vector<const EncAnyGram *> > & tokenrevindex, EncData * sentence, unordered_set<const EncAnyGram*> * patterns, bool includeskipgrams ) {
    for (int i = 0; i < sentence->length(); i++) {
        for (unordered_set<const EncAnyGram*>::const_iterator iter = patterns->begin(); iter !=  patterns->end(); iter++) {
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
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
    load(s2tmodel, t2smodel, s2tthreshold, t2sthreshold, productthreshold);
}


AlignmentModel::AlignmentModel(AlignmentModel & s2tmodel, AlignmentModel & t2smodel,  const double s2tthreshold, const double t2sthreshold, const double productthreshold, bool DEBUG) {
    this->DEBUG = DEBUG;
    sourcemodel = NULL;
    targetmodel = NULL;
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
    load(s2tmodel, t2smodel, s2tthreshold, t2sthreshold, productthreshold);
}

AlignmentModel::~AlignmentModel() {
    //TODO: More and better cleanup
/*    for (unordered_map<const EncAnyGram*, unordered_map<const EncAnyGram*, unordered_map<const EncAnyGram*, double> > >::iterator iter = keywords.begin(); iter != keywords.end(); iter++) {
        for (unordered_map<const EncAnyGram*, unordered_map<const EncAnyGram*, double> >::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            for (unordered_map<const EncAnyGram*, double>::iterator iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                const EncAnyGram * keygram = iter2->first;
                if (getsourcekey(keygram) == NULL) {
                    delete keygram;
                }
            }
        }
    }
    */

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


AlignmentModel::AlignmentModel(const string & filename, bool logprobs, int ptsfield, bool allowskipgrams, const int bestn, bool DEBUG, const int bestnkeywords, const double keywordprobthreshold) {
    this->DEBUG = DEBUG;
    sourcemodel = NULL;
    targetmodel = NULL;
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
    this->ptsfield = ptsfield;
    load(filename,logprobs, allowskipgrams, bestn, bestnkeywords, keywordprobthreshold);
}

void AlignmentModel::load(const string & filename, bool logprobs, bool allowskipgrams, const int bestn, const int bestnkeywords, const double keywordprobthreshold) {
	unsigned char check;

    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "Aligment model file does not exist: " << filename << endl;
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
    if (DEBUG) cerr << "model_id=" << model_id << endl;
    if ((model_id >= (unsigned int) ALIGNMENTMODEL) && (model_id <= 101))  {
        multiscore = false; //backward compatibility with old models
    } else {
        multiscore = true;
    }
    if (model_id >= (unsigned int) ALIGNMENTMODEL + 4) {
        f.read( (char*) &leftsourcecontext, sizeof(unsigned char));
        f.read( (char*) &rightsourcecontext, sizeof(unsigned char));
    } else {
        leftsourcecontext = 0;
        rightsourcecontext = 0;
    }
    int ngramversion = 1;
    if (model_id < (unsigned int) ALIGNMENTMODEL + 5) {
        ngramversion = 0;
    }
    bool dokeywords = true;
    if (model_id < (unsigned int) ALIGNMENTMODEL + 6) {
        dokeywords = false;
    }

    if (DEBUG) {
        cerr << "leftsourcecontext=" << (int) leftsourcecontext << endl;
        cerr << "rightsourcecontext=" << (int) rightsourcecontext << endl;
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
        	throw InternalError();
        }
        f.read(&gapcount, sizeof(char));

        const EncAnyGram * sourcegram;
        bool sourceisskipgram = false;
        if (gapcount == 0) {
            if (DEBUG)  cerr << "\tNGRAM";
            const EncNGram * ngram = new EncNGram(&f, ngramversion); //read from file
            if (DEBUG)  cerr << " n=" << (int) ngram->n() << " size=" << (int) ngram->size();
            sourcegram = getsourcekey((EncAnyGram*) ngram); //does the key already exist?
            if (sourcegram == NULL) {
                //no
                alignmatrix[(const EncAnyGram*) ngram];
                sourcegram = getsourcekey((const EncAnyGram*) ngram);
                if (sourcegram == NULL) { cerr << "INTERNAL ERROR: sourcegram still not found after insertion! Should never happen!";throw InternalError(); }
            	//sourcengrams.insert(ngram);
            } else {
                //yes
                delete ngram; //ngram not necessary
            }
        } else {
            if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
            const EncSkipGram * skipgram = new EncSkipGram( &f, gapcount, ngramversion); //read from file
            if (allowskipgrams) {
                sourcegram = getsourcekey((EncAnyGram*) skipgram);  //does the key already exist?
                if (sourcegram == NULL) {
                    alignmatrix[(const EncAnyGram*) skipgram];
                	//sourceskipgrams.insert(skipgram);
                	sourcegram = getsourcekey((const EncAnyGram*) skipgram);
                    if (sourcegram == NULL) { cerr << "INTERNAL ERROR: sourcegram still not found after insertion! Should never happen!"; throw InternalError(); }
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
            if (sourcegram->n() - leftsourcecontext - rightsourcecontext <= 0) {
                cerr << "INTERNAL ERROR: unable to remove context from n-gram, nothing left! nwithcontext=" << (int) sourcegram->n() << " leftcontextsize=" << (int) leftsourcecontext << " rightcontextsize=" << (int) rightsourcecontext <<  endl;
                cerr << sourcegram->out() << endl;
                throw InternalError();
            }
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
		        EncNGram ngram = EncNGram(&f, ngramversion); //read from file
		        if (DEBUG)  cerr << " n=" << (int) ngram.n() << " size=" << (int) ngram.size();
		        if (!gettargetkey((EncAnyGram*) &ngram)) {
		        	targetngrams.insert(ngram);
		        }
		        targetgram = gettargetkey((EncAnyGram*) &ngram);
		    } else {
		        if (DEBUG)  cerr << "\tSKIPGRAM, " << (int) gapcount << " gaps";
		        EncSkipGram skipgram = EncSkipGram( &f, gapcount, ngramversion); //read from file
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
		        if (DEBUG) cerr << " scores=" << (int) scores;
		    } else {
		        scores = 1;
		    }

            double p;
		    for (int i = 0; i < scores; i++) {
		        f.read((char*) &p, sizeof(double));
		        if (DEBUG) cerr << " score: " << p;
		        if ((p > 0) && (logprobs)) p = log(p); //base e
		        if ((allowskipgrams) || ((!sourceisskipgram) && (!targetisskipgram))) {
           		    if ((sourcegram == NULL) || (targetgram == NULL)) {
		             	cerr << "SOURCEGRAM or TARGETGRAM is NULL";
		            	throw InternalError();
		            }
		            alignmatrix[sourcegram][targetgram].push_back(p);
		        }
		    }

		    if (dokeywords) {

                const EncAnyGram * sourcegramfocus;
                if ((leftsourcecontext != 0) || (rightsourcecontext != 0)) {
                    const EncAnyGram * tmp = sourcegram->slice(leftsourcecontext, sourcegram->n() - leftsourcecontext - rightsourcecontext);
                    sourcegramfocus = getsourcekey(tmp);
                    delete tmp;
                } else {
                    sourcegramfocus = sourcegram;
                }                               

		        uint32_t keywordcount;
		        f.read((char*) &keywordcount, sizeof(uint32_t));
		        for (int i = 0; i < keywordcount; i++) {
                    f.read(&gapcount, sizeof(char));
                    if (gapcount == 0) {
                        if (DEBUG)  cerr << "\tKEYWORD-NGRAM";
                        EncNGram * ngram = new EncNGram(&f, ngramversion); //read from file
                        if (DEBUG)  cerr << " n=" << (int) ngram->n() << " size=" << (int) ngram->size();
                        const EncAnyGram * sourcekey = getsourcekey((EncAnyGram*) ngram); //for keyword
                        double keywordprob;
                        f.read((char*) &keywordprob, sizeof(double));
                        if (sourcekey != NULL) {
                            if ((i <= bestnkeywords) && (keywordprob >= keywordprobthreshold))  keywords[sourcegramfocus][targetgram][sourcekey] = keywordprob;
                            delete ngram;
                        } else {
                            if ((i <= bestnkeywords) && (keywordprob >= keywordprobthreshold)) {
                                keywords[sourcegramfocus][targetgram][ngram] = keywordprob;
                            } else {
                                delete ngram;
                            }
                        }
                    } else {
                        if (DEBUG)  cerr << "\tKEYWORD-SKIPGRAM, " << (int) gapcount << " gaps";
                        EncSkipGram * skipgram = new EncSkipGram( &f, gapcount, ngramversion); //read from file
                        double keywordprob;
                        f.read((char*) &keywordprob, sizeof(double));
                        if (allowskipgrams) {
                            const EncAnyGram * sourcekey = getsourcekey((EncAnyGram*) skipgram); //for keyword
                            if (sourcekey != NULL) {
                                if ((i <= bestnkeywords) && (keywordprob >= keywordprobthreshold)) keywords[sourcegramfocus][targetgram][sourcekey] = keywordprob;
                                delete skipgram;
                            } else {
                                if ((i <= bestnkeywords) && (keywordprob >= keywordprobthreshold)) {
                                    keywords[sourcegramfocus][targetgram][skipgram] = keywordprob;
                                } else {
                                    delete skipgram;
                                }
                            }
                        } else {
                            delete skipgram;
                        }
                    }
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
		if (DEBUG)  cerr << endl;
	}
    f.close();

}

AlignmentModel::AlignmentModel(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs, int ptsfield, bool DEBUG) {
    //load from moses-style phrasetable file
    this->DEBUG = DEBUG;
    sourcemodel = NULL;
    targetmodel = NULL;
    this->debug_sourceclassdecoder = NULL;
    this->debug_targetclassdecoder = NULL;
    this->ptsfield = ptsfield;
    leftsourcecontext = rightsourcecontext = 0;
    load(filename, sourceencoder, targetencoder, logprobs, ptsfield);
}

void AlignmentModel::load(const std::string & filename, ClassEncoder * sourceencoder, ClassEncoder * targetencoder, bool logprobs, int ptsfield) {
    //load from moses-style phrasetable file
    this->ptsfield = ptsfield;
    leftsourcecontext = rightsourcecontext = 0;

    ifstream f;
    f.open(filename.c_str(), ios::in | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "Aligment model file does not exist: " << filename << endl;
       exit(3);
    }

    while (!f.eof()) {
        string line;
        getline(f, line);
        int mode = 0;
        string source = "";
        string target = "";
        string scores_s;
        vector<double> scores;
        int begin = 0;
        for (unsigned int i = 0; i < line.size(); i++) {
            if (line.substr(i,5) == " ||| ") {
                if (mode == 0) {
                    source = line.substr(begin, i - begin);
                } else if (mode == 1) {
                    target = line.substr(begin, i - begin);
                } else if (mode == 2) {
                    scores_s = line.substr(begin, i - begin);
                }
                begin = i+5;
                mode++;
            }
        }
        scores_s = scores_s + " ";
        begin = 0;
        //cerr << "DEBUG: scores_s=" << scores_s << endl;
        for (unsigned int i = 0; i < scores_s.size(); i++) {
            if ((scores_s[i] == ' ')  && (i > begin)) {
                double score = atof(scores_s.substr(begin, i - begin).c_str());
                //cerr << scores_s.substr(begin, i - begin) << " -> " << score << endl;
                if ((score > 0) && (logprobs)) score = log(score); //base e
                scores.push_back(score);
                begin = i + 1;
            }
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


void AlignmentModel::save(const string & filename, const int bestnkeywords) {
	const unsigned char check = 0xff;
	const char czero = 0;
    unordered_set<const EncAnyGram *> processedkws; //temporary map will store processed keywords, if context is presents, keywords will only be stored at the first occurrence of the focus

    ofstream f;
    f.open(filename.c_str(), ios::out | ios::binary);
    if ((!f) || (!f.good())) {
       cerr << "Aligment model file unable to save: " << filename << endl;
       exit(3);
    }

    uint64_t _id = ALIGNMENTMODEL + ALIGNMENTMODELVERSION;
    f.write( (char*) &_id, sizeof(uint64_t));

    f.write( (char*) &leftsourcecontext, sizeof(unsigned char));
    f.write( (char*) &rightsourcecontext, sizeof(unsigned char));

    uint64_t sourcecount = alignmatrix.size();
    if (alignmatrix.count(NULLGRAM) > 0) sourcecount--;
    f.write( (char*) &sourcecount, sizeof(uint64_t));

    unsigned int i = 0;

    for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
        i++;
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
        	const EncAnyGram* targetgram = gettargetkey(iter2->first, true);
        	if (targetgram == NULL) {
        	    cerr << "AlignmentModel::save(): Target key not found! This should not happen! Error whilst processing source-pattern " << i << endl;
        	    throw InternalError();
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



            const EncAnyGram * sourcegramfocus;
            if ((leftsourcecontext != 0) || (rightsourcecontext != 0)) {
                const EncAnyGram * tmp = sourcegram->slice(leftsourcecontext, sourcegram->n() - leftsourcecontext - rightsourcecontext);
                sourcegramfocus = getsourcekey(tmp);
                delete tmp;
            } else {
                sourcegramfocus = sourcegram;
            }                               
                    

        	if ((keywords.count(sourcegramfocus)) && (keywords[sourcegramfocus].count(targetgram))) {

                if (processedkws.count(sourcegramfocus)) {
                    f.write(&czero, sizeof(char));
                } else {
                    processedkws.insert(sourcegramfocus);

                    uint32_t keywordcount = keywords[sourcegramfocus][targetgram].size();
                    if (keywordcount > bestnkeywords) keywordcount = bestnkeywords;
                    f.write( (char*) &keywordcount, sizeof(uint32_t));
                    //sort before saving
                    multimap<double, const EncAnyGram *> sortedkeywords;
                    for (unordered_map<const EncAnyGram*, double>::iterator iter3 = keywords[sourcegramfocus][targetgram].begin(); iter3 != keywords[sourcegramfocus][targetgram].end(); iter3++) {
                        sortedkeywords.insert(pair<double, const EncAnyGram *>(iter3->second, iter3->first));
                    }
                    int kwcount = 0;
                    for (multimap<double, const EncAnyGram*>::iterator iter3 = sortedkeywords.begin(); iter3 != sortedkeywords.end(); iter3++) {
                        kwcount++;
                        if (kwcount > bestnkeywords) break;
                        const EncAnyGram * keyword = iter3->second;
                        if (keyword->isskipgram()) {
                            const EncSkipGram * skipgram = (const EncSkipGram*) keyword;
                            skipgram->writeasbinary(&f);
                        } else {
                            const EncNGram * ngram = (const EncNGram*) keyword;
                            f.write(&czero, sizeof(char)); //gapcount, always zero for ngrams
                            ngram->writeasbinary(&f);
                        }
                        const double p = iter3->first;
                        f.write( (char*) &p, sizeof(double));
                    }

                }
        	} else {
        	    const uint32_t keywordcount = 0;
        	    f.write( (char*) &keywordcount, sizeof(uint32_t));
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
            if ( (keywords.count(sourcegram)) && (keywords[sourcegram].count(targetgram)) && (!mosesformat) ) {
                if (DEBUG) cerr << "keywordcount=" << keywords[sourcegram][targetgram].size() << endl;
                *OUT << endl << "#KEYWORDS:";
                multimap<double,const EncAnyGram *> sortedkw;
                for (unordered_map<const EncAnyGram *, double>::iterator iter3 = keywords[sourcegram][targetgram].begin(); iter3 != keywords[sourcegram][targetgram].end(); iter3++) sortedkw.insert(pair<double,const EncAnyGram *>(-1* iter3->second, iter3->first));
                for (multimap<double,const EncAnyGram *>::iterator iter3 = sortedkw.begin(); iter3 != sortedkw.end(); iter3++) *OUT << "\t" << iter3->second->decode(sourceclassdecoder) << ' ' << (iter3->first * -1);
            }
            *OUT << endl;
        }
    }
}

const EncAnyGram * AlignmentModel::getfocuskey(const EncAnyGram * key) {
    if (!leftsourcecontext && !rightsourcecontext) return getsourcekey(key); //no context
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
        //cout << "Trying sourcemodel" << endl;
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


const EncAnyGram * AlignmentModel::gettargetkey(const EncAnyGram* key, bool returnselfifnotfound, bool forcemodel) {
    if (targetmodel != NULL) {
        const EncAnyGram * modelkey = targetmodel->getkey(key);
        if (modelkey != NULL) {
            return modelkey;
        } else if (forcemodel) {
            return returnselfifnotfound ? key : NULL;
        }
    }
    if (returnselfifnotfound && targetngrams.empty() && targetskipgrams.empty()) return key; //no targetngrams/targetskipgrams, we'll just have to assume key is a valid pointer
    if (key->gapcount() == 0) {
        std::unordered_set<EncNGram>::iterator iter = targetngrams.find(*( (EncNGram*) key) );
        if (iter != targetngrams.end()) {
            return &(*iter);
        } else {
            return returnselfifnotfound ? key : NULL;
        }
    } else {
        std::unordered_set<EncSkipGram>::iterator iter = targetskipgrams.find(*( (EncSkipGram*) key) );
        if (iter != targetskipgrams.end()) {
            return &(*iter);
        } else {
            return returnselfifnotfound ? key : NULL;
        }
    }
}


void AlignmentModel::computereverse() {
    reversealignmatrix.clear();
    for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
        const EncAnyGram * source = iter->first;
        for (t_aligntargets::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram * target = iter2->first;
            reversealignmatrix[target][source] = alignmatrix[source][target];
        }
    }
}


t_aligntargets AlignmentModel::sumtranslationoptions(const EncAnyGram * sourcefocus, bool debug) {

        //compute translate options, aggregating context-data into non-context based scores
        int scorevectorsize = 0;

        t_aligntargets translationoptions;
        for (unordered_set<const EncAnyGram*>::const_iterator iter = sourcecontexts[sourcefocus].begin(); iter != sourcecontexts[sourcefocus].end(); iter++) {
            const EncAnyGram * sourcekey = *iter;
            for (t_aligntargets::iterator iter2 = alignmatrix[sourcekey].begin(); iter2 != alignmatrix[sourcekey].end(); iter2++) {
                if (scorevectorsize == 0) {
                    scorevectorsize = iter2->second.size();
                     if (scorevectorsize > 2) {
                        cerr << "ERROR: Score vector contains more than two scores. Unknown how to interpret these in sumtranslationoptions() (due to model having context)" << endl;
                        throw InternalError();
                     }
                }
                const EncAnyGram * targetgram = iter2->first;
                if (translationoptions.count(targetgram) == 0) {
                    //targetgram does not exist yet
                    for (int i = 0; i < iter2->second.size(); i++) {
                        const double p =  pow(exp(1), iter2->second[i]);
                        translationoptions[targetgram].push_back(p);
                    }
                } else {
                    //targetgram exists, sum
                    for (int i = 0; i < iter2->second.size(); i++) {
                        const double p =  pow(exp(1), iter2->second[i]);
                        translationoptions[targetgram][i] += p;
                    }
                }
            }
        }

        //compute total
        double total = 0;
        for (t_aligntargets::iterator iter = translationoptions.begin(); iter != translationoptions.end(); iter++) {
            total += iter->second[0];
        }


        if ((scorevectorsize == 2) && (reversealignmatrix.empty())) computereverse();

        //convert computed scores back to logprobs
        for (t_aligntargets::iterator iter = translationoptions.begin(); iter != translationoptions.end(); iter++) {
            const EncAnyGram * targetgram = iter->first;

            translationoptions[targetgram][0] = log(translationoptions[targetgram][0] / total);


            if (scorevectorsize == 2) {
                double revtotal = 0;
                //normalize reverse probability
                for (t_aligntargets::iterator iter2 = reversealignmatrix[targetgram].begin(); iter2 != reversealignmatrix[targetgram].end(); iter2++) {
                    revtotal += pow(exp(1),iter2->second[1]);
                }
                translationoptions[targetgram][1] = log(translationoptions[targetgram][1] / revtotal);

            }


        }

        return translationoptions;
}


AlignmentModel * AlignmentModel::removecontext() {
    if (!this->leftsourcecontext && !this->rightsourcecontext) {
        cerr << "ERROR: Model has no context information" << endl;
        throw InternalError();
    }
    AlignmentModel * newmodel = new AlignmentModel();
    for (t_contexts::iterator iter = sourcecontexts.begin(); iter != sourcecontexts.end();  iter++) {
        const EncAnyGram * sourcekey = iter->first;
        newmodel->alignmatrix[sourcekey] = sumtranslationoptions(sourcekey);
    }
    newmodel->targetngrams = targetngrams;
    newmodel->targetskipgrams = targetskipgrams;
    return newmodel;
}


void AlignmentModel::stats() {
    unsigned int sourcecount = alignmatrix.size();
    unsigned int totalcount = 0;
    int scorecount = 0;
    if (leftsourcecontext || rightsourcecontext) {
        sourcecount = sourcecontexts.size();
        for (t_contexts::iterator iter = sourcecontexts.begin(); iter != sourcecontexts.end();  iter++) {
            for (unordered_set<const EncAnyGram *>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                totalcount += alignmatrix[*iter2].size();
                scorecount = alignmatrix[*iter2].begin()->second.size();
            }
        }
    } else {
        for (t_alignmatrix::iterator iter = alignmatrix.begin(); iter != alignmatrix.end(); iter++) {
            totalcount += iter->second.size();
            if (iter->second.begin() != iter->second.end()) scorecount = iter->second.begin()->second.size();
        }
    }
    cout << "sources aligned: " << sourcecount << endl;
    cout << "total alignments: " << totalcount << endl;
    cout << "average alignments per source: " << (double) totalcount / sourcecount << endl;
    cout << "left context: " << (int) leftsourcecontext << endl;
    cout << "right context: " << (int) rightsourcecontext << endl;
    if (leftsourcecontext || rightsourcecontext) cout << "source contexts: " << alignmatrix.size() << endl;
    cout << "score vector size: " << scorecount << endl;
}


int AlignmentModel::computekeywords(IndexedPatternModel & sourcepatternmodel, IndexedPatternModel & targetpatternmodel, int include_threshold, int absolute_threshold, double probability_threshold, int filter_threshold, int bestnkeywords) {
    int keywordsfound = 0;
    int c = 0;
    int total = sourcepatternmodel.ngrams.size() + sourcepatternmodel.skipgrams.size();
    for (unordered_map<const EncNGram,NGramData >::iterator iter = sourcepatternmodel.ngrams.begin(); iter != sourcepatternmodel.ngrams.end(); iter++) {
        c++;
        const EncAnyGram * sourcegram = (const EncAnyGram *) &(iter->first); //without context
        if (iter->second.count() >= include_threshold) keywordsfound += computekeywords(sourcepatternmodel, targetpatternmodel, sourcegram, absolute_threshold, probability_threshold, filter_threshold, bestnkeywords);
        if ((DEBUG) || (c % 10000 == 0)) cerr << " Computekeywords @" << c << "/" << total << " -- " << keywordsfound << " keywords found in total thus far" << endl;
    }
    for (unordered_map<const EncSkipGram,SkipGramData >::iterator iter = sourcepatternmodel.skipgrams.begin(); iter != sourcepatternmodel.skipgrams.end(); iter++) {
        c++;
        const EncAnyGram * sourcegram = (const EncAnyGram *) &(iter->first); //without context
        if (iter->second.count() >= include_threshold) keywordsfound += computekeywords(sourcepatternmodel, targetpatternmodel, sourcegram, absolute_threshold, probability_threshold, filter_threshold, bestnkeywords);
        if ((DEBUG) || (c % 10000 == 0)) cerr << " Computekeywords @" << c << "/" << total << " -- " << keywordsfound << " keywords found in total thus far" << endl;
    }
    return keywordsfound;
}


int AlignmentModel::computekeywords(SelectivePatternModel & sourcepatternmodel, SelectivePatternModel & targetpatternmodel, int include_threshold, int absolute_threshold, double probability_threshold, int filter_threshold, int bestnkeywords) {
    int keywordsfound = 0;
    int c = 0;
    int total = sourcepatternmodel.ngrams.size() + sourcepatternmodel.skipgrams.size();
    for (unordered_map<const EncNGram,IndexCountData >::iterator iter = sourcepatternmodel.ngrams.begin(); iter != sourcepatternmodel.ngrams.end(); iter++) {
        c++;
        const EncAnyGram * sourcegram = (const EncAnyGram *) &(iter->first); //without context
        if (iter->second.count >= include_threshold) keywordsfound += computekeywords(sourcepatternmodel, targetpatternmodel, sourcegram, absolute_threshold, probability_threshold, filter_threshold, bestnkeywords);
        if ((DEBUG) || (c % 10000 == 0)) cerr << " Computekeywords @" << c << "/" << total << " -- " << keywordsfound << " keywords found in total thus far" << endl;
    }
    for (unordered_map<const EncSkipGram,IndexCountData >::iterator iter = sourcepatternmodel.skipgrams.begin(); iter != sourcepatternmodel.skipgrams.end(); iter++) {
        c++;
        const EncAnyGram * sourcegram = (const EncAnyGram *) &(iter->first); //without context
        if (iter->second.count >= include_threshold) keywordsfound +=computekeywords(sourcepatternmodel, targetpatternmodel, sourcegram, absolute_threshold, probability_threshold, filter_threshold, bestnkeywords);
        if ((DEBUG) || (c % 10000 == 0)) cerr << " Computekeywords @" << c << "/" << total << " -- " << keywordsfound << " keywords found in total thus far" << endl;
    }
    return keywordsfound;
}





//Bit ugly, code duplication! Fix sometime... Make sure I edit the right one!
int AlignmentModel::computekeywords(IndexedPatternModel & sourcepatternmodel, IndexedPatternModel & targetpatternmodel, const EncAnyGram * sourcegram, int absolute_threshold, double probability_threshold , int filter_threshold, int bestnkeywords) {
    int keywordsfound = 0;

    //const EncAnyGram * sourcegram    is without context
    //keywords are stored in keywords[] without context

    const EncAnyGram * sourcefocuskey = getfocuskey(sourcegram); //without context

    if ((leftsourcecontext == 0) && (rightsourcecontext == 0)) {
        //cheat, temporarily add to sourcecontexts
        const EncAnyGram * sourcekey = getsourcekey(sourcegram); 
        if (sourcekey == NULL) return 0;
        sourcecontexts[sourcefocuskey].insert(sourcekey);
    }

    for (std::unordered_set<const EncAnyGram *>::iterator ctiter = sourcecontexts[sourcefocuskey].begin(); ctiter != sourcecontexts[sourcefocuskey].end(); ctiter++) {
        const EncAnyGram * sourcekey = *ctiter; //with context

        if ( (alignmatrix.count(sourcekey)) && (alignmatrix[sourcekey].size() > 1)  && (sourcepatternmodel.exists(sourcegram))) { //don't bother searching keywords if there is only one translation for a sourcegram

            unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> > countmap; // targetgram -> key -> count

            for (t_aligntargets::iterator iter = alignmatrix[sourcekey].begin(); iter != alignmatrix[sourcekey].end(); iter++) {
                    const EncAnyGram * targetgram = iter->first;

                    if (targetpatternmodel.exists(targetgram)) {
                        set<int> sentenceconstraints = targetpatternmodel.getsentences(targetgram);
                        countmap[targetgram] = sourcepatternmodel.getcooccurrences(sourcegram, NULL, &sentenceconstraints); ////returns map for counting key -> counts
                    }
            }



            keywordsfound = computekeywords(&sourcepatternmodel, sourcekey,sourcegram,countmap, absolute_threshold, probability_threshold , filter_threshold, bestnkeywords);
        }

    }

    if ((leftsourcecontext == 0) && (rightsourcecontext == 0)) sourcecontexts.clear(); //cleanup cheat
    

    return keywordsfound;

}



int AlignmentModel::computekeywords(ModelQuerier * sourcepatternmodel, const EncAnyGram * sourcekey, const EncAnyGram * sourcegram, unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> > & countmap, int absolute_threshold, double probability_threshold , int filter_threshold, int bestnkeywords) {

    //sourcegram and sourcekey are without context, and are stored as such in
    //keywords[]. No need to concern about context here

    int keywordsfound = 0;

    if (DEBUG) {
        cerr << "\tTemporary map for " << countmap.size() << " target patterns" << endl;
    }


    double Nkloc = 0;
    for (unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> >::iterator iter = countmap.begin(); iter != countmap.end(); iter++) {
        unordered_map<const EncAnyGram *, int>::iterator keyiter = iter->second.begin();
        while (keyiter != iter->second.end()) {
            const EncAnyGram * keygram = keyiter->first;
            //exclude keywords that are already subsumed by the sourcegram  (does not work for skipgrams yet)
            if ((keygram->gapcount() == 0) && (sourcegram->gapcount() == 0) && (((const EncNGram *) sourcegram)->contains((const EncNGram *) keygram))) {
                if (DEBUG) cerr << "\tKey is subpart, deleting key " << iter->second.size();
                keyiter = countmap[iter->first].erase(keyiter);
                if (DEBUG) cerr << " " << iter->second.size() << endl;
            } else {
                Nkloc += keyiter->second;
                keyiter++;
            }
        }
    }


    multimap<double, const EncAnyGram *> tmp; //temporary map holding all keywords

    for (unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> >::iterator iter = countmap.begin(); iter != countmap.end(); iter++) {
        const EncAnyGram * targetgram = iter->first;
        for (unordered_map<const EncAnyGram *, int>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
            const EncAnyGram * keygram = iter2->first;
            //compute probability of translation given keyword

            if (DEBUG) cerr << "\tConsidering keyword - ";

            const double Ns_kloc = iter2->second;
            const double Nkcorp = sourcepatternmodel->occurrencecount(keygram);
            if( (Ns_kloc >= absolute_threshold) && (Nkcorp >= filter_threshold)) {
                const double p = (Ns_kloc / Nkloc) * (1/Nkcorp);
                if (p >= probability_threshold) {
                    //add to keywords
                    keywords[sourcekey][targetgram][keygram] = p;
                    tmp.insert(pair<double, const EncAnyGram *>(p, keygram));
                    keywordsfound++;
                    if (DEBUG) cerr << "\tACCEPTED p=" << p << " Ns_kloc=" << Ns_kloc << " Nkloc=" << Nkloc << " Nkcorp=" << Nkcorp << endl;
                } else if (DEBUG) {
                    cerr << "\trejected p=" << p << endl;
                }
            } else if (DEBUG) {
                 cerr << "\trejected Ns_kloc=" << Ns_kloc << " Nkcorp=" << Nkcorp << endl;
            }

        }
    }

    //remove keywords that have no disambiguative power, i.e keywords that are present regardless for all targets
    //also remove keywords above bestnkeywords threshold!
    int kwcount = 0;
    for (multimap<double,const EncAnyGram *>::iterator tmpiter = tmp.begin(); tmpiter != tmp.end(); tmpiter++) {
        kwcount++;
        const EncAnyGram * keygram = tmpiter->second;
        bool omnipresent = true;
        for (unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> >::iterator iter = countmap.begin(); iter != countmap.end(); iter++) {
            if (iter->second.count(keygram)) {
                omnipresent = false;
                break;
            }
        }
        if ((omnipresent) || (kwcount > bestnkeywords)) {
            unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> >::iterator iter = countmap.begin();
            while (iter != countmap.end()) {
                const EncAnyGram * targetgram = iter->first;
                iter->second.erase(keygram);
                if (iter->second.empty()) {
                    iter = countmap.erase(iter);
                } else {
                    iter++;
                }
            }
        }
    }

   return keywordsfound;
}

int AlignmentModel::computekeywords(SelectivePatternModel & sourcepatternmodel, SelectivePatternModel & targetpatternmodel, const EncAnyGram * sourcegram, int absolute_threshold, double probability_threshold , int filter_threshold, int bestnkeywords) {

    //const EncAnyGram * sourcegram    is without context
    //keywords are stored in keywords[] without context
    
    int keywordsfound = 0;
    const EncAnyGram * sourcefocuskey = getfocuskey(sourcegram);

    if ((leftsourcecontext == 0) && (rightsourcecontext == 0)) {
        //cheat, temporarily add to sourcecontexts
        const EncAnyGram * sourcekey = getsourcekey(sourcegram); 
        if (sourcekey == NULL) return 0;
        sourcecontexts[sourcefocuskey].insert(sourcekey);
    }

    for (unordered_set<const EncAnyGram *>::iterator ctiter = sourcecontexts[sourcefocuskey].begin(); ctiter != sourcecontexts[sourcefocuskey].end(); ctiter++) {
        const EncAnyGram * sourcekey = *ctiter; //with context

        if (alignmatrix.count(sourcekey)  && (alignmatrix[sourcekey].size() > 1) && (sourcepatternmodel.exists(sourcekey))) {

        unordered_map<const EncAnyGram *, unordered_map<const EncAnyGram *, int> > countmap; // targetgram -> key -> count

        for (t_aligntargets::iterator iter = alignmatrix[sourcekey].begin(); iter != alignmatrix[sourcekey].end(); iter++) {
                const EncAnyGram * targetgram = iter->first;

                if (targetpatternmodel.exists(targetgram)) {
                    set<int> sentenceconstraints = targetpatternmodel.getsentences(targetgram);
                    if (sentenceconstraints.empty()) cerr << "\tWARNING: sentenceconstraints is empty set" << endl;
                    countmap[targetgram] = sourcepatternmodel.getcooccurrences(sourcegram, NULL, &sentenceconstraints); ////returns map for counting key -> counts
                    if (countmap[targetgram].empty()) cerr << "\tWARNING: No co-occurrences found for a source pattern" << endl;
                }
        }

        keywordsfound = computekeywords(&sourcepatternmodel, sourcekey,sourcegram,countmap, absolute_threshold, probability_threshold , filter_threshold, bestnkeywords);

        }

    }

    if ((leftsourcecontext == 0) && (rightsourcecontext == 0)) sourcecontexts.clear(); //cleanup cheat

    return keywordsfound;
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


