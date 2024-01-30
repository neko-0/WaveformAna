#ifndef ANA_QUICK_SCAN_H
#define ANA_QUICK_SCAN_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

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
  double *posx = 0;
  double *posy = 0;
 
  double **posx_addr;

  // output variables
  double *output_pmax[num_ch_];
  double *output_tmax[num_ch_];
  double *output_rise[num_ch_];
  double *output_area[num_ch_];
  double *output_fwhm[num_ch_];
  double *output_rms[num_ch_];
  double *output_x, *output_y, *output_z;

  std::vector<double> *output_w[num_ch_];
  std::vector<double> *output_t;
  // std::vector<double> *output_cfd[num_ch_];
};

#endif
