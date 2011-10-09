#include <string>
#include "classdecoder.h"

using namespace std;

int main( int argc, char *argv[] ) {
    const string cls = argv[1];
    const string clsdec = argv[2];
    ClassDecoder classdecoder = ClassDecoder(cls);
    classdecoder.decodefile(clsdec);   
}
