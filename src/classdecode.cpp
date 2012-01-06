#include <string>
#include <iostream>
#include "classdecoder.h"
#include "getopt.h"

using namespace std;

void usage() {
    cerr << "Syntax: classdecoder -f encoded-corpus -c class-file" << endl;
    cerr << "Descriptions: Decodes an encoded corpus" << endl;
}

int main( int argc, char *argv[] ) {    
    string classfile = "";
    string corpusfile = "";
       
    char c;    
    while ((c = getopt(argc, argv, "c:f:h")) != -1)
        switch (c)
        {
        case 'c':
            classfile = optarg;
            break;
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
        
    if (classfile.empty() || corpusfile.empty()) {
    	usage();
    	exit(2);
    }
    
    ClassDecoder classdecoder = ClassDecoder(classfile);
    classdecoder.decodefile(corpusfile);   
}
