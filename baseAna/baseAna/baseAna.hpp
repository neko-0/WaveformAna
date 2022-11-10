#ifndef BASE_ANA_H
#define BASE_ANA_H

#include "configMgr/betaConfigMgr.hpp"

#include <string>
#include <vector>


struct BaseAna {

  BaseAna();
  virtual ~BaseAna();

  virtual void initialize(BetaConfigMgr *configMgr) = 0;
  virtual void execute(BetaConfigMgr *configMgr) = 0;
  virtual void finalize(BetaConfigMgr *configMgr) = 0;
};


#endif
