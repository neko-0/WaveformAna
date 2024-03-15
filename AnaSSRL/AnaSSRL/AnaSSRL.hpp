#ifndef ANA_SSRL_H
#define ANA_SSRL_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaSSRL : BaseAna {

  AnaSSRL(){};
  ~AnaSSRL(){};

  virtual void setup(BetaConfigMgr* const configMgr);
  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);
  virtual bool execute(BetaConfigMgr* const configMgr);

private:
  void bucket_time_difference(
    const int &ch,
    const double &bucket_start,
    const double &bucket_step);

  void prepare_fix_window_branches(BetaConfigMgr* const configMgr);
  void fill_fix_window_branches(
    int ch,
    const std::vector<double> &v_trace,
    const std::vector<double> &t_trace,
    double t_min,
    double t_max,
    bool fill_previous=false,
    double threshold=0.0);

  void find_max_ch(
    const std::vector<int> &chlist,
    std::vector<int> &buffer,
    std::vector<double> &output,
    const std::vector<int> &sumCh,
    int targetCh = -1,
    int targetCh2 = -1,
    double scale = 0.8,
    int max_1st = -1);

  void find_max_ch(
    const std::vector<int> &chlist,
    std::vector<int> &buffer,
    std::vector<double> &output);

  void regular_routine(std::vector<double> &corr_w, int ch);
  void simple_routine(std::vector<double> &corr_w, int ch);
  void trigger_routine(std::vector<double> &corr_w, int ch);

private:
  const int ch_start_ = 0;
  static const int num_ch_ = 16;
  std::vector<int> active_ch_ = {};
  bool store_waveform = true;
  bool use_single_t_trace = true;
  bool use_single_input_t_trace = true;
  bool found_single_t_trace = false;

  bool skip_it = false;

  double threshold = 0.0;

  double bucket_t_start_ = -44e-9;
  double bucket_t_end_ = 2.08e-9;
  int nbuckets_= 28;

  bool fill_fix_window = true;
  double fix_win_start_ = 560.0;
  double fix_win_step_size_ = 11;
  int fix_win_nstep_ = 40;
  double fix_win_edge_dist_ = 0.0;

  bool do_max_ch_ = false;

  int baseline_opt = 0;
  int run_type = 0;

  std::vector<int> invert_ch;
  std::vector<int> simple_ana_ch;


  int trigger_ch = -1;

  // ===========================================================================
  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t[num_ch_];

  std::vector<double> *trig;

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

  std::vector<float> *output_bucket_pmax[num_ch_];
  std::vector<float> *output_bucket_area[num_ch_];
  std::vector<float> *output_bucket_tmax[num_ch_];
  std::vector<float> *output_bucket_cfd20[num_ch_];
  std::vector<float> *output_bucket_cfd50[num_ch_];
  std::vector<int> *output_bucket_index[num_ch_];
  std::vector<float> *output_bucket_tmax_diff[num_ch_];
  std::vector<float> *output_bucket_cfd20_diff[num_ch_];
  std::vector<float> *output_bucket_cfd50_diff[num_ch_];
  float *output_bucket_corr[num_ch_];

  std::vector<int> *output_max_ch;
  std::vector<int> *output_2nd_max_ch;
  std::vector<int> *output_small_pad_max_ch;

  std::vector<double> *output_sum_large;
  std::vector<double> *output_sum_small;
  std::vector<double> *output_sum_strip_set1; // [12, 5, 6]
  std::vector<double> *output_sum_strip_set2; // [11, 7, 7]
  std::vector<double> *output_sum_strip_set3; // [11, 7, 7]

  std::vector<float> *output_w[num_ch_];
  std::vector<float> *output_corr_w[num_ch_];
  std::vector<float> *output_t[num_ch_];

  std::vector<float> *output_fix_pmax[num_ch_];
  std::vector<float> *output_fix_tmax[num_ch_];
  std::vector<float> *output_fix_area[num_ch_];

  double *thresholdTime[num_ch_];

  double *trig_time;

};

#endif
