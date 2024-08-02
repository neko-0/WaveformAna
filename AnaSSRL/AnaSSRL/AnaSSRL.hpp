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

  void prepare_bunch_window_branches(BetaConfigMgr* const configMgr);
  void fill_bunch_window_branches(
    int ch,
    const std::vector<double> &v_trace,
    const std::vector<double> &t_trace,
    double t_min,
    double t_max,
    bool fill_previous=false,
    double threshold=0.0);
  void fill_bunch_window_wp(int ch);

  void prepare_leading_signal_branches(BetaConfigMgr* const configMgr);
  void fill_leading_signal_branches(
    int ch,
    const std::vector<double> &v_trace,
    const std::vector<double> &t_trace,
    double t_min,
    double t_max);
  
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
  void scan_routinue(std::vector<double> &corr_w, int ch);

private:
  int ch_start_ = 0;
  static const int num_ch_ = 16;
  std::vector<int> active_ch_ = {};
  bool store_waveform = true;
  bool use_single_t_trace = true;
  bool use_single_input_t_trace = true;
  bool found_single_t_trace = false;

  int routine_ = 1;

  bool skip_it = false;

  double threshold = 0.0;

  int bunch_nstep_ = 40;
  double bunch_start_ = 560.0;
  double bunch_step_size_ = 11;
  double bunch_edge_dist_ = 0.0;

  double leading_tmin_ = 0.0;
  double leading_tmax_ = 0.0;

  double rms_start_ = 0.0;
  double rms_end_ = 0.0;

  bool do_max_ch_ = false;

  int baseline_opt = 0;
  int run_type = 0;

  std::vector<int> invert_ch;
  std::vector<int> simple_ana_ch;

  std::vector<double> corr_common_time_;

  int trigger_ch = -1;

  // ===========================================================================
  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t[num_ch_];

  std::vector<double> *trg0;
  std::vector<double> *trg1;
  
  std::vector<double> *pos;

  // ===========================================================================
  // output variables
  bool *output_basecorr[num_ch_];
  float *output_rms[num_ch_];

  std::vector<float> *output_pmax[num_ch_];
  std::vector<float> *output_tmax[num_ch_];
  std::vector<float> *output_rise[num_ch_];
  std::vector<float> *output_fall[num_ch_];
  std::vector<float> *output_area[num_ch_];
  std::vector<float> *output_20cfd[num_ch_];
  std::vector<float> *output_50cfd[num_ch_];
  
  std::vector<bool> *output_rms_wp_tight[num_ch_];
  std::vector<bool> *output_rms_wp_loose[num_ch_];
  std::vector<bool> *output_bunch_wp_tight[num_ch_];
  std::vector<bool> *output_bunch_wp_loose[num_ch_];
  std::vector<bool> *output_fall_wp_tight[num_ch_];
  std::vector<bool> *output_fall_wp_loose[num_ch_];
  std::vector<bool> *output_wp_tight[num_ch_];
  std::vector<bool> *output_wp_loose[num_ch_];

  
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

  float *output_leading_fall[num_ch_];
  float *output_leading_pmax[num_ch_];
  float *output_leading_tmax[num_ch_];
  float *output_leading_area[num_ch_];
  float *output_leading_rise[num_ch_];
  float *output_leading_20cfd[num_ch_];
  float *output_leading_50cfd[num_ch_];

  double *thresholdTime[num_ch_];

  double *trig_time;

  double trg_threshold_time_;

  double *output_x, *output_y;

};

#endif
