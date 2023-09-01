#ifndef BASE_ANA_H
#define BASE_ANA_H

#include "configMgr/betaConfigMgr.hpp"

#include <string>
#include <vector>


struct BaseAna {

  BaseAna();
  virtual ~BaseAna();

  virtual void setup(BetaConfigMgr* const configMgr) = 0;
  virtual void initialize(BetaConfigMgr* const configMgr) = 0;
  virtual void finalize(BetaConfigMgr* const configMgr) = 0;
  virtual bool execute(BetaConfigMgr* const configMgr) = 0;
};


#endif
