#include "utilities/getFiles.hpp"
#include "utilities/logger.hpp"

#include <filesystem>
#include <string>

namespace fs = std::filesystem;

std::vector<std::string> GetListOfFiles(const std::string &dir)
{
  std::vector<std::string> files = {};

  // looping through all of the files under directory
  for(const auto &entry : fs::directory_iterator(dir.c_str()))
  {
    // removing the directory part
    std::string filename = entry.path();
    filename.erase(0, filename.find(dir) + dir.length());
    files.push_back(filename);
  }
  return files;
}

std::vector<std::string> GetListOfFiles(
    const std::string &dir,
    const std::string &pattern)
{
  std::vector<std::string> files = GetListOfFiles(dir);
  std::vector<std::string> filtered_files = {};

  for(const auto &filename : files){
    // spliting pattern based on wildcard symbol *
    std::vector<std::string> pattern_fragment = {};
    std::string pattern_copy = pattern;
    std::string wildcard_delimiter = "*";
    std::size_t pos = 0;
    while((pos = pattern_copy.find(wildcard_delimiter)) != std::string::npos)
    {
      std::string pattern_fragment = pattern_copy.substr(0, pos);
      pattern_copy.erase(0, pos + wildcard_delimiter.length());
    }

    pattern_fragment.push_back(pattern_copy);

    bool matched = false;
    for(const auto &frag : pattern_fragment) {
      if (filename.find(frag) != std::string::npos)
        matched = true;
      else {
        matched = false;
        break;
      }
    }

    if (matched) {
      std::string file_name = dir + "/" + filename;
      filtered_files.push_back(file_name);
    }
  }

  return filtered_files;
}
