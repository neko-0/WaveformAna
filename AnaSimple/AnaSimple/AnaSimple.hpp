#ifndef ANA_SIMPLE_H
#define ANA_SIMPLE_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaSimple : BaseAna {

  AnaSimple(){};
  ~AnaSimple(){};

  virtual void setup(BetaConfigMgr* const configMgr);
  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);
  virtual bool execute(BetaConfigMgr* const configMgr);

private:
    std::vector<int> invert_channels_;
    double search_tmin_ = -2e-9;
    double search_tmax_ = 2e-9;
    bool store_waveform_ = true;
    bool use_single_input_time_ = false;
    int ch_start_ = 0;
    std::vector<int> active_ch_;
    bool found_single_t_trace_ = false;
    
    // ===========================================================================
    // input waveform channels
    static const int num_ch_ = 4;
    std::vector<double> *w[num_ch_];
    std::vector<double> *t[num_ch_];

    // ===========================================================================
    // output variables
    double *output_pmax[num_ch_];
    double *output_tmax[num_ch_];
    double *output_rms[num_ch_];
    double *output_area[num_ch_];
    
    std::vector<float> *output_w_corr[num_ch_];
    std::vector<float> *output_time[num_ch_];
};

#endif
