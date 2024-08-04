#ifndef ANA_HIGH_BW_H
#define ANA_HIGH_BW_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaHighBW : BaseAna {

  AnaHighBW(){};
  ~AnaHighBW(){};

  virtual void setup(BetaConfigMgr* const configMgr);
  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);
  virtual bool execute(BetaConfigMgr* const configMgr);

private:
    
    // ===========================================================================
    // input waveform channels
    static const int num_ch_ = 4;
    std::vector<double> *w[num_ch_];
    std::vector<double> *t[num_ch_];

};

#endif
