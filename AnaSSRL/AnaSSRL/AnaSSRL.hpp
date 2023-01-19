#ifndef ANA_SSRL_H
#define ANA_SSRL_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaSSRL : BaseAna {

  AnaSSRL(){};
  ~AnaSSRL(){};

  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void execute(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);

private:
  const int ch_start_ = 1;
  static const int num_ch_ = 4;
  std::vector<int> active_ch_ = {};
  bool store_waveform = true;
  bool use_single_t_trace = true;
  bool found_single_t_trace = false;

  const double fixed_win_min = -0.5e-9;
  const double fixed_win_max = 2.5e-9;

  // ===========================================================================
  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t[num_ch_];

  // ===========================================================================
  // output variables
  bool *output_basecorr[num_ch_];
  int *output_nsignal[num_ch_];
  float *output_rms[num_ch_];

  std::vector<float> *output_rise[num_ch_];
  std::vector<float> *output_area[num_ch_];
  std::vector<float> *output_fwhm[num_ch_];
  std::vector<float> *output_20cfd[num_ch_];
  std::vector<float> *output_50cfd[num_ch_];
  std::vector<float> *output_pmax[num_ch_];
  std::vector<float> *output_tmax[num_ch_];
  std::vector<float> *output_tmax_diff[num_ch_];
  std::vector<float> *output_raw_pmax[num_ch_];

  std::vector<float> *output_w[num_ch_];
  std::vector<float> *output_corr_w[num_ch_];
  std::vector<float> *output_t[num_ch_];
};

#endif
