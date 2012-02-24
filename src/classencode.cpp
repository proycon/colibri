#include <string>
#include <iostream>
#include "classencoder.h"
#include "getopt.h"
#include <common.h>

using namespace std;

void usage() {
    cerr << "Syntax: classencoder -f corpus" << endl;
    cerr << "Descriptions: Encodes a corpus" << endl;
}

int main( int argc, char *argv[] ) {    
    string classfile = "";
    string corpusfile = "";
    string outputprefix = "";
       
    char c;    
    while ((c = getopt(argc, argv, "f:h")) != -1)
        switch (c)
        {
        case 'f':
            corpusfile = optarg;
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

    ClassEncoder classencoder = ClassEncoder();
    classencoder.build(corpusfile);
    classencoder.save(outputprefix + ".cls");
    classencoder.encodefile(corpusfile, outputprefix + ".clsenc");
}
