#ifndef GETFILES_H
#define GETFILES_H

#include <vector>
#include <string>

std::vector<std::string> GetListOfFiles(const std::string &dir);
std::vector<std::string> GetListOfFiles(const std::string &dir, const std::string &pattern);



#endif // BETACOPE_UTILITIES_DIR_H
