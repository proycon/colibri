#include <common.h>

std::string get_extension(const std::string& filename) {
    if(filename.find_last_of(".") != std::string::npos)
        return filename.substr(filename.find_last_of(".")+1);
    return "";
}

bool strip_extension(std::string& filename, const std::string extension) {
    if(filename.find_last_of(".") != std::string::npos)
        if (filename.substr(filename.find_last_of(".")+1) == extension) {
            filename = filename.substr(0,filename.find_last_of("."));
            return true;
        }
    return false;
}
