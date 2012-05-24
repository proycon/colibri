#include <string>
#include <iostream>
#include "classencoder.h"
#include "getopt.h"
#include <common.h>

using namespace std;

void usage() {
    cerr << "Syntax: classencoder -f corpus [ -c classmodel ]" << endl;
    cerr << "Description: Encodes a corpus. If used with -c, encodes a corpus according to the specified pre-existing class model" << endl;
}

int main( int argc, char *argv[] ) {    
    string classfile = "";
    string corpusfile = "";
    string outputprefix = "";
       
    char c;    
    while ((c = getopt(argc, argv, "f:h:c:")) != -1)
        switch (c)
        {
        case 'f':
            corpusfile = optarg;
            break;
        case 'c':
            classfile = optarg;
            break;                
        case 'h':
            usage();
            exit(0);
		default:
            cerr << "Unknown option: -" <<  optopt << endl;
            abort ();
        }
        
    if (corpusfile.empty()) {
    	usage();
    	exit(2);
    }
    
    if (outputprefix.empty()) {
        outputprefix = corpusfile; 
        strip_extension(outputprefix,"txt");    
    }


    ClassEncoder classencoder;
    
    bool allowunknown = false;
    
    if (!classfile.empty()) {
        cerr << "Loading classes from file" << endl;
        classencoder = ClassEncoder(classfile);
        allowunknown = true;
        cerr << "Building classes from corpus (extending existing classes)" << endl;
        classencoder.build(corpusfile);
        classencoder.save(outputprefix + ".cls");          
    } else {
        cerr << "Building classes from corpus" << endl;
        classencoder = ClassEncoder();
        classencoder.build(corpusfile);
        classencoder.save(outputprefix + ".cls");        
    }        
    classencoder.encodefile(corpusfile, outputprefix + ".clsenc", allowunknown);
}
