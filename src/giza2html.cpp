#include <getopt.h>
#include "alignmodel.h"


using namespace std;

void usage() {
    cerr << "Usage: giza2html -W gizamodel [-R reverse-giza-model] -S source-class-file -T target-class-file" << endl;
    cerr << " Input:" << endl;
    cerr << "\t-W gizamodel              GIZA++ word-alignment model (*.A3.final)" << endl;    
    cerr << "\t-R gizamodel              Reverse GIZA++ word-alignment model (*.A3.final)"  << endl;
    cerr << "\t-S sourceclassfile        Source class file (for decoding)" << endl;
    cerr << "\t-T targetclassfile        Target class file (for decoding)" << endl;
    cerr << "\t-i                        Compute intersection" << endl;
    cerr << "\t-u                        Compute union" << endl;
    cerr << "\t-b [sentencenum]          Start output beginning at this sentence number" << endl;
    cerr << "\t-e [sentencenum]          End output beginning at this sentence number" << endl;
    cerr << "\t-s [sourcepatternmodel]   Source pattern model (for use with -x)" << endl;
    cerr << "\t-t [targetpatternmodel]   Target pattern model (for use with -x)" << endl;
    cerr << "\t-x [threshold]            Extract phrases (using given alignmentthreshold, requires -s and -t)" << endl;    
}




int main( int argc, char *argv[] ) {
    string gizast = "";
    string gizats = "";
    string sourceclassfile = "";
    string targetclassfile = "";
    string sourcemodelfile = "";
    string targetmodelfile = "";

    bool INTERSECTION = false;
    bool UNION = false;

    unsigned int begin = 0;
    unsigned int end = 0; 

    bool DOEXTRACT = false;
    double alignthreshold = 0.5;

    
    
    char c;    
    while ((c = getopt(argc, argv, "W:R:S:T:iuhb:e:x:s:t:")) != -1)
        switch (c)
        {       
            case 0:
                break;
            case 'i':
            	INTERSECTION = true;
            	break;
            case 'u':
            	UNION = true;
            	break;        	
            case 'h':
            	usage();
            	exit(0);
            case 'W':
                gizast = optarg;
                break;
            case 'R':
                gizats = optarg;
                break;
            case 'S':
                sourceclassfile = optarg;
                break;
            case 'T':
                targetclassfile = optarg;
                break; 
            case 'b':
                begin = atoi(optarg);
                break;
            case 'e':
                end = atoi(optarg);
                break;
            case 's':
                sourcemodelfile = optarg;
                break;
            case 't':
                targetmodelfile = optarg;
                break;                
            case 'x':
                DOEXTRACT = true;
                alignthreshold = atof(optarg);
                break;
            default:
                cerr << "Unknown option: -" <<  optopt << endl;
                abort ();
        }
        

        if (sourceclassfile.empty() || targetclassfile.empty() || gizast.empty()) {
            usage();
            exit(2);
        }

        

        if (DOEXTRACT && (sourcemodelfile.empty() || targetmodelfile.empty())) {
            cerr << "ERROR: Using -x requires -s and -t" << endl;
            usage();
            exit(2);                
        }
        

        

		cerr << "Loading source class encoder " << sourceclassfile << endl;
	    ClassDecoder sourceclassdecoder = ClassDecoder(sourceclassfile);
	    ClassEncoder sourceclassencoder = ClassEncoder(sourceclassfile);

	    cerr << "Loading target class encoder " << targetclassfile << endl;
		ClassDecoder targetclassdecoder = ClassDecoder(targetclassfile);    
		ClassEncoder targetclassencoder = ClassEncoder(targetclassfile);
			
	    cerr << "Initialising GIZA++ Word Alignments" << endl;
	    GizaModel gizamodel_s2t = GizaModel(gizast, &sourceclassencoder, &targetclassencoder);	   

        cout << "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" /><title>Colibri: Word Alignments</title><style type=\"text/css\">div { background: #eee; margin: 10px; margin-left: auto; margin-right: auto; width: 100%; text-align: center }\ntable { margin: 10px; border-collapse: collapse; } \ntable td { border: solid 1px #999; }</style></head><body>" << endl;

        unsigned int sentence = 0;
	    if (!gizats.empty()) {

	    
	        if (DOEXTRACT) {
                cerr << "Loading source model " << sourcemodelfile << endl;
	            SelectivePatternModel sourcemodel = SelectivePatternModel(sourcemodelfile);
	            cerr << "Loading target model " << targetmodelfile << endl;
	            SelectivePatternModel targetmodel = SelectivePatternModel(targetmodelfile);
	        
	            AlignmentModel alignmentmodel = AlignmentModel(&sourcemodel,&targetmodel);
	        
	            GizaModel gizamodel_t2s = GizaModel(gizats, &targetclassencoder, &sourceclassencoder);
	            while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {
	                sentence++;
	                GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
	                GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();
	                
	                if ((sentence >= begin) && ((sentence <= end) || (end == 0))) {
	                    cout << "<div>" << endl;	
	                    cout << "<h3>" << sentence << "</h3>" << endl;         
                        sentence_s2t.out( (ostream*) &cout,  sourceclassdecoder, targetclassdecoder);
                        sentence_t2s.out( (ostream*) &cout,  targetclassdecoder, sourceclassdecoder);
                        if (INTERSECTION) {
                            GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
                            sentence_i.out( (ostream*) &cout, sourceclassdecoder, targetclassdecoder);
                        }
                        if (UNION) {
                            GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);
                            sentence_u.out( (ostream*) &cout, sourceclassdecoder, targetclassdecoder);
                        }
                        cout << "</div>" << endl;
                        cout << "<pre>" << endl;
                        alignmentmodel.extractgizapatterns(sentence_s2t, sentence_t2s, sentence, 0,0,alignthreshold, &sourceclassdecoder, &targetclassdecoder);
                        alignmentmodel.alignmatrix.clear(); //no need to retain matrix in memory
                        cout << "</pre>" << endl;  
                        cout << "<hr />" << endl;
                    }
                }	        
	        } else {
	            GizaModel gizamodel_t2s = GizaModel(gizats, &targetclassencoder, &sourceclassencoder);
	            while (!gizamodel_s2t.eof() && !gizamodel_t2s.eof()) {
	                sentence++;
	                GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
	                GizaSentenceAlignment sentence_t2s = gizamodel_t2s.readsentence();
	                
	                if ((sentence >= begin) && ((sentence <= end) || (end == 0))) {
	                    cout << "<div>" << endl;	            
                        sentence_s2t.out( (ostream*) &cout,  sourceclassdecoder, targetclassdecoder);
                        sentence_t2s.out( (ostream*) &cout,  targetclassdecoder, sourceclassdecoder);
                        if (INTERSECTION) {
                            GizaSentenceAlignment sentence_i = sentence_s2t.intersect(sentence_t2s);
                            sentence_i.out( (ostream*) &cout, sourceclassdecoder, targetclassdecoder);
                        }
                        if (UNION) {
                            GizaSentenceAlignment sentence_u = sentence_s2t.unify(sentence_t2s);
                            sentence_u.out( (ostream*) &cout, sourceclassdecoder, targetclassdecoder);
                        }
                        cout << "</div><hr />" << endl;
                    }                    	       
	            }
	        }         
	    } else {
	       while (!gizamodel_s2t.eof()) {
	            sentence++;        
                GizaSentenceAlignment sentence_s2t = gizamodel_s2t.readsentence();
                if ((sentence >= begin) && ((sentence <= end) || (end == 0))) {
                    cout << "<div>" << endl;
                    sentence_s2t.out( (ostream*) &cout,  sourceclassdecoder, targetclassdecoder);
                    cout << "</div><hr />" << endl;
                }
                       
            }
	    }
	    
	   cout << "</body></html>" << endl;
	
	
}
