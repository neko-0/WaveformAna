#include "AnaTCT/AnaTCT.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace wm = waveform_methods;

void AnaTCT::initialize(BetaConfigMgr* const configMgr){

  t = configMgr->SetInputBranch<std::vector<double>>("t0");

  posx = configMgr->SetInputBranch<double>("x");
  posy = configMgr->SetInputBranch<double>("y");
  posz = configMgr->SetInputBranch<double>("z");

  output_t = configMgr->SetOutputBranch<std::vector<double>>("t");
  output_x = configMgr->SetOutputBranch<double>("x");
  output_y = configMgr->SetOutputBranch<double>("y");
  output_z = configMgr->SetOutputBranch<double>("z");

  for(int i = 0; i < this->num_ch_; i++){
    // input branches
    std::string current_ch = std::to_string(ch_start_ + i);
    w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);

    // output branches
    output_pmax[i] = configMgr->SetOutputBranch<double>("pmax" + current_ch);
    output_tmax[i] = configMgr->SetOutputBranch<double>("tmax" + current_ch);
    output_rise[i] = configMgr->SetOutputBranch<double>("rise" + current_ch);
    output_area[i] = configMgr->SetOutputBranch<double>("area" + current_ch);
    output_fwhm[i] = configMgr->SetOutputBranch<double>("fwhm" + current_ch);
    output_rms[i] = configMgr->SetOutputBranch<double>("rms" + current_ch);

    // output_cfd[i] = configMgr->SetOutputBranch<std::vector<double>>("cfd" + current_ch);
    output_w[i] = configMgr->SetOutputBranch<std::vector<double>>("w" + current_ch);
  }
}

void AnaTCT::execute(BetaConfigMgr* const configMgr){

  const double v_scale = 1.0; // V to mV
  const double t_scale = 1.0; // s to ps

  std::move(t->begin(), t->end(), std::back_inserter(*output_t));
  // std::cout << "posx " << *posx << "\n";
  *output_x = *posx;
  *output_y = *posy;
  *output_z = *posz;

  for(int ch = 0; ch < this->num_ch_; ch++){
    int trace_size = w[ch]->size();
    double baseline = wm::Baseline::CalcBaseline(*w[ch], 0, (int)(0.25*trace_size));
    // inverting signal
    for(int i=0; i < trace_size; i++){
      w[ch]->at(i) -= baseline;
      w[ch]->at(i) *= -1.0*v_scale;
      // t[ch]->at(i) *= t_scale;
      output_w[ch]->push_back(w[ch]->at(i));
      // output_t[ch]->push_back(t[ch]->at(i));
    }

    auto wave_pt = wm::FindSignalMax(*w[ch], *t);
    auto rise = wm::CalcRiseTime(*w[ch], *t, wave_pt.index);
    auto area = wm::CalcPulseArea(*w[ch], *t, wave_pt.index);
    auto fwhm = wm::CalcFWHM(*w[ch], *t, wave_pt.index);
    auto rms = wm::CalcNoise(*w[ch], 0.25);
    // auto cfd_times = wm::CalcCFDTime(*w[ch], *t, wave_pt.index, 0.1, 0.1);

    *output_pmax[ch] = wave_pt.v;
    *output_tmax[ch] = wave_pt.t;
    *output_rise[ch] = rise;
    *output_area[ch] = area;
    *output_fwhm[ch] = fwhm;
    *output_rms[ch] = rms;

    // std::move(cfd_times.begin(), cfd_times.end(), std::back_inserter(*output_cfd[ch]));
  }

}

void AnaTCT::finalize(BetaConfigMgr* const configMgr){
  // pass
}
