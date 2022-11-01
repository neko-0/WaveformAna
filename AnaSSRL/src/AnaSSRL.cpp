#include "AnaSSRL/AnaSSRL.hpp"
#include "waveformMethods/waveformMethods.hpp"

#include <vector>
#include <string>

namespace wm = waveform_methods;

void AnaSSRL::initialize(BetaConfigMgr *configMgr){

  for(int i = 0; i < this->num_ch_; i++){
    // input branches
    std::string current_ch = std::to_string(ch_start_ + i);
    w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
    t[i] = configMgr->SetInputBranch<std::vector<double>>("t" + current_ch);

    // output branches
    output_pmax[i] = configMgr->SetOutputBranch<double>("pmax" + current_ch);
    output_tmax[i] = configMgr->SetOutputBranch<double>("tmax" + current_ch);
    output_rise[i] = configMgr->SetOutputBranch<double>("rise" + current_ch);
    output_area[i] = configMgr->SetOutputBranch<double>("area" + current_ch);

    output_w[i] = configMgr->SetOutputBranch<std::vector<double>>("w" + current_ch);
    output_t[i] = configMgr->SetOutputBranch<std::vector<double>>("t" + current_ch);
  }
}

void AnaSSRL::execute(BetaConfigMgr *configMgr){

  for(int ch = 0; ch < this->num_ch_; ch++){
    // inverting signal
    for(int i=0; i < w[ch]->size(); i++){
      w[ch]->at(i) *= 1.0;
      output_w[ch]->push_back(w[ch]->at(i));
      output_t[ch]->push_back(t[ch]->at(i));
    }

    auto wave_pt = wm::FindSignalMax(*w[ch], *t[ch]);
    auto rise = wm::CalcRiseTime(*w[ch], *t[ch], wave_pt.index);
    auto area = wm::CalcPulseArea(*w[ch], *t[ch], wave_pt.index);

    *output_pmax[ch] = wave_pt.v;
    *output_tmax[ch] = wave_pt.t;
    *output_rise[ch] = rise;
    *output_area[ch] = area;
  }

}

void AnaSSRL::finalize(BetaConfigMgr *configMgr){
  // pass
}
