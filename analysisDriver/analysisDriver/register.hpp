#ifndef ANA_REGISTER_H
#define ANA_REGISTER_H

#include "baseAna/baseAna.hpp"

#include <string>
#include <vector>
#include <map>

//==============================================================================
template<typename T>
BaseAna *CreateUserAna(){return new T;};

//==============================================================================
struct AnalysisFactory{
  typedef std::map<std::string, BaseAna*(*)()> ana_map;

  static BaseAna *SelectAnalysis(const std::string &name){
    auto &my_map = GetMap();
    auto it = my_map.find(name);
    if(it == my_map.end()) return 0;
    return it->second();
  }

  static bool CheckAnalysis(const std::string &name){
    auto &my_map = GetMap();
    auto it = my_map.find(name);
    if(it == my_map.end()) return false;
    return true;
  }

  static ana_map &GetMap(){
    static ana_map *my_map = new ana_map;
    return *my_map;
  }
};

//==============================================================================
// template<typename T>
// struct AnalysisRegister:AnalysisFactory{
//   AnalysisRegister(const std::string &name){
//     auto &my_map = GetMap();
//     my_map.insert(std::make_pair(name, &CreateUserAna<T>));
//   }
// };

struct AnalysisRegister:AnalysisFactory{
  AnalysisRegister(){};

  template<typename T>
  static bool Register(const std::string &name){
    auto &my_map = GetMap();
    my_map.insert(std::make_pair(name, &CreateUserAna<T>));
    return true;
  }
};

#endif
