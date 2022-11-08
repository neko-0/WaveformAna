#ifndef ANA_SSRL_H
#define ANA_SSRL_H

#include "analysisDriver/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaSSRL : BaseAna {

  AnaSSRL(){};
  ~AnaSSRL(){};

  virtual void initialize(BetaConfigMgr *configMgr);
  virtual void execute(BetaConfigMgr *configMgr);
  virtual void finalize(BetaConfigMgr *configMgr);

private:
  const int ch_start_ = 2;
  static const int num_ch_ = 2;

  // input variables
  std::vector<double> *w[num_ch_];
  std::vector<double> *t[num_ch_];

  // output variables
  double *output_pmax[num_ch_];
  double *output_tmax[num_ch_];
  double *output_rise[num_ch_];
  double *output_area[num_ch_];
  double *output_fwhm[num_ch_];

  std::vector<double> *output_w[num_ch_];
  std::vector<double> *output_t[num_ch_];
  std::vector<double> *output_cfd[num_ch_];

  static AnalysisRegister<AnaSSRL> reg;
};


AnalysisRegister<AnaSSRL> AnaSSRL::reg("SSRL");

#endif
