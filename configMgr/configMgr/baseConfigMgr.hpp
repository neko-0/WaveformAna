#ifndef BASECONFIGMGR_H
#define BASECONFIGMGR_H

#include <vector>
#include <string>
#include <map>

class BaseConfigMgr {
  std::string name_;
public:
  BaseConfigMgr(){};
  virtual ~BaseConfigMgr(){};

  virtual bool Initialize() = 0;
  virtual bool NextEvent() = 0;
  virtual bool Fill() = 0;
  virtual bool Finalize() = 0;

  virtual std::string GetName() = 0;
};

template<typename T>
BaseConfigMgr *CreateT(){return new T;}

template<typename T>
T *_InitUserConfigMgr(BaseConfigMgr *base){return dynamic_cast<T*>(base);}

//==============================================================================
struct ConfigMgrFactory{
  typedef std::map<std::string, BaseConfigMgr*(*)()> config_map;

  static BaseConfigMgr *CreateConfigMgr(const std::string &name){
    auto *my_map = GetMap();
    auto it = my_map->find(name);
    if(it == my_map->end()) return 0;
    return it->second();
  };

  template <typename T>
  static T *CreateConfigMgrT(const std::string &name){
    auto *base = CreateConfigMgr(name);
    return _InitUserConfigMgr<T>(base);
  }

  template <typename T>
  static T *InitUserConfigMgr(BaseConfigMgr *base){
    return _InitUserConfigMgr<T>(base);
  }

protected:
  static config_map *GetMap(){
    static config_map my_map;
    return &my_map;
  };

private:
  static config_map my_map;
};


template<typename T>
struct ConfigMgrRegister:ConfigMgrFactory{

  ConfigMgrRegister(const std::string &name){
    auto *my_map = GetMap();
    my_map->insert(std::make_pair(name, &CreateT<T>));
  }

};

#endif
