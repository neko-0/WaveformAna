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

    int ch_start_ = 1;
    int dut_ch_ = 1;
    std::vector<int> enbale_ch_;
    std::vector<int> active_ch_;
    double winsize_ = 0.25e-6;
    double peak_finding_threshold_ = 0.0015;
    
    // ===========================================================================
    // input waveform channels
    static const int num_ch_ = 4;
    std::vector<double> *i_w[num_ch_];
    std::vector<double> *i_t[num_ch_];

    std::vector<double> *o_w[num_ch_];
    std::vector<double> *o_t[num_ch_];
    
    std::vector<double> *dut_process_time;
    std::vector<double> *dut_process_voltage;
    
    std::vector<double> *w_step;
    std::vector<double> *t_step;
};

#endif
