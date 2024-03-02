#ifndef ANA_QUICK_SCAN_H
#define ANA_QUICK_SCAN_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <TProfile2D.h>

#include <vector>

struct AnaQuickScan : BaseAna {

  AnaQuickScan(){};
  ~AnaQuickScan(){};

  virtual void setup(BetaConfigMgr* const configMgr);
  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);
  virtual bool execute(BetaConfigMgr* const configMgr);

private:
  const int ch_start_ = 0;
  static const int num_ch_ = 16;

  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t;
  std::vector<double> *pos;

  // output variables
  double *output_pmax[num_ch_];
  double *output_x, *output_y, *output_z;

  TProfile2D *heatmap[num_ch_];
};

#endif
