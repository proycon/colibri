#ifndef COLIBRICOMMON_H
#define COLIBRICOMMON_H

#include <string>
#include <list>
#include <vector>

std::string trim(const std::string &t, const std::string &ws);
std::string get_extension(const std::string& filename);
bool strip_extension(std::string& filename, const std::string extension);
double listproduct(const std::vector<double> & l);
double listsum(const std::vector<double> & l);
void orderedinsert(std::list<double> & l, double value);

#endif
