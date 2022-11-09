#ifndef BASEANA_H
#define BASEANA_H

#include "configMgr/betaConfigMgr.hpp"

#include <string>
#include <vector>


struct BaseAna {

  BaseAna(){};
  virtual ~BaseAna(){};

  virtual void initialize(BetaConfigMgr *configMgr) = 0;
  virtual void execute(BetaConfigMgr *configMgr) = 0;
  virtual void finalize(BetaConfigMgr *configMgr) = 0;
};

//==============================================================================
template<typename T>
BaseAna *CreateUserAna(){return new T;}

//==============================================================================
struct AnalysisFactory{
  typedef std::map<std::string, BaseAna*(*)()> ana_map;

  static BaseAna *SelectAnalysis(const std::string &name){
    auto *my_map = GetMap();
    auto it = my_map->find(name);
    if(it == my_map->end()) return 0;
    return it->second();
  }

  static bool CheckAnalysis(const std::string &name){
    auto *my_map = GetMap();
    auto it = my_map->find(name);
    if(it == my_map->end()) return false;
    return true;
  }

  static ana_map *GetMap(){
    static ana_map my_map;
    return &my_map;
  }

private:
  inline static ana_map my_map;
};

//==============================================================================
template<typename T>
struct AnalysisRegister:AnalysisFactory{
  AnalysisRegister(const std::string &name){
    auto *my_map = GetMap();
    my_map->insert(std::make_pair(name, &CreateUserAna<T>));
  }
};

#endif
