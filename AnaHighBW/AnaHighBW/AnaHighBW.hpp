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
    int clock_ch_ = 2;
    std::vector<int> enbale_ch_;
    std::vector<int> active_ch_;
    double winsize_ = 0.25e-6;
    double peak_finding_threshold_ = 0.0015;

    int cap_t_start_;
    int cap_t_win_;
    int cap_cycle_size_;
    int cap_ncycles_;
    
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

    void offset_correction(std::vector<double> &v_trace);
    void find_capacitors(const int dut_ch); // find capacitor based on signal shape and features.
    void find_capacitors(const int dut_ch, const int clock_ch); // use the clock channel to find capacitors.
    void find_capacitors(
      const int dut_ch, 
      const int t_start,
      const int t_win,
      const int cycle_size,
      const int ncycles); // easiest approach assuming clock is relatively statble. 
};

#endif
