#ifndef ANA_SSRL_H
#define ANA_SSRL_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaSSRL : BaseAna {

  AnaSSRL(){};
  ~AnaSSRL(){};

  virtual void initialize(BetaConfigMgr *configMgr);
  virtual void execute(BetaConfigMgr *configMgr);
  virtual void finalize(BetaConfigMgr *configMgr);

private:
  const int ch_start_ = 1;
  static const int num_ch_ = 4;
  std::vector<int> active_ch_ = {};
  bool store_waveform = true;
  bool use_single_t_trace = true;
  bool found_single_t_trace = false;

  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t[num_ch_];

  // output variables
  bool *output_basecorr[num_ch_];
  int *output_nsignal[num_ch_];
  std::vector<double> *output_rise[num_ch_];
  std::vector<double> *output_area[num_ch_];
  std::vector<double> *output_fwhm[num_ch_];
  std::vector<double> *output_20cfd[num_ch_];
  std::vector<double> *output_50cfd[num_ch_];
  std::vector<double> *output_pmax[num_ch_];
  std::vector<double> *output_tmax[num_ch_];


  std::vector<double> *output_w[num_ch_];
  std::vector<double> *output_t[num_ch_];
};

#endif
