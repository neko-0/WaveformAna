#include "AnaSSRL/AnaSSRL.hpp"
#include "waveformMethods/waveformMethods.hpp"

#include <vector>
#include <string>
#include <omp.h>

namespace wm = waveform_methods;

void AnaSSRL::initialize(BetaConfigMgr *configMgr){

  for(int i = 0; i < this->num_ch_; i++){
    // input branches
    std::string current_ch = std::to_string(ch_start_ + i);
    w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
    if(!w[i]) continue;
    t[i] = configMgr->SetInputBranch<std::vector<double>>("t" + current_ch);

    active_ch_.push_back(i);

    // output branches
    output_basecorr[i] = configMgr->SetOutputBranch<bool>("basecorr" + current_ch);
    output_nsignal[i] = configMgr->SetOutputBranch<int>("nsignal" + current_ch);
    output_rms[i] = configMgr->SetOutputBranch<double>("rms" + current_ch);
    output_pmax[i] = configMgr->SetOutputBranch<std::vector<double>>("pmax" + current_ch);
    output_tmax[i] = configMgr->SetOutputBranch<std::vector<double>>("tmax" + current_ch);
    output_rise[i] = configMgr->SetOutputBranch<std::vector<double>>("rise" + current_ch);
    output_area[i] = configMgr->SetOutputBranch<std::vector<double>>("area" + current_ch);
    output_fwhm[i] = configMgr->SetOutputBranch<std::vector<double>>("fwhm" + current_ch);
    output_20cfd[i] = configMgr->SetOutputBranch<std::vector<double>>("20cfd" + current_ch);
    output_50cfd[i] = configMgr->SetOutputBranch<std::vector<double>>("20cfd" + current_ch);
    output_tmax_diff[i] = configMgr->SetOutputBranch<std::vector<double>>("tmaxdiff" + current_ch);

    if(store_waveform){
      output_w[i] = configMgr->SetOutputBranch<std::vector<double>>("w" + current_ch);
      if(use_single_t_trace){
        if(found_single_t_trace){
          output_t[i] = nullptr;
          continue;
        }
        output_t[i] = configMgr->SetOutputBranch<std::vector<double>>("t");
        found_single_t_trace = true;
        continue;
      }
      output_t[i] = configMgr->SetOutputBranch<std::vector<double>>("t" + current_ch);
    }
  }

  LOG_INFO("number of active chanenls: " + std::to_string(active_ch_.size()));

  LOG_INFO("external config: " + configMgr->ext_config_name());
}

void AnaSSRL::execute(BetaConfigMgr *configMgr){

  #pragma omp parallel for
  for(auto &ch : active_ch_){
    if(w[ch]->size() == 0){
      LOG_WARNING("Trace size 0");
      continue;
    }

    auto mix_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);

    // check polarity of the signal. determine threshold for pmax finder
    double polarity;
    // double abs_max;
    if(mix_params.pos_max > mix_params.neg_max*-1.0){
      polarity = 1.0;
      // abs_max = mix_params.pos_max;
    } else {
      polarity = -1.0;
      // abs_max = -1.0*mix_params.neg_max;
    }
    double threshold = 5.0 * mix_params.rms;

    for(int i=0; i < w[ch]->size(); i++){
      w[ch]->at(i) = (w[ch]->at(i) - mix_params.baseline)*polarity;
      // output_w[ch]->emplace_back(w[ch]->at(i));
      // output_t[ch]->emplace_back(t[ch]->at(i));
    }

    auto n_wave_pts = wm::FindMultipleSignalMax(*w[ch], *t[ch], threshold);
    // auto n_wave_pts = wm::FindMultipleSignalMaxAlt1(*w[ch], *t[ch], threshold, 5);

    *output_basecorr[ch] = wm::MultiSignalBaselineCorrection(
      n_wave_pts, *w[ch], *t[ch], 0.5, 5e-9, 5e-9);

    // std::move(cfd_times.begin(), cfd_times.end(), std::back_inserter(*output_cfd[ch]));
    // #pragma omp parallel for
    // for(int pt_i=0; pt_i < n_wave_pts.size(); pt_i++){
    //   auto pt = n_wave_pts.at(pt_i);
    for(auto &pt : n_wave_pts){
      if(pt.index < 0) continue;
      if(!output_tmax[ch]->empty()){
        output_tmax_diff[ch]->push_back(pt.t - output_tmax[ch]->back());
      } else {
        output_tmax_diff[ch]->push_back(0.0);
      }
      output_pmax[ch]->push_back(pt.v);
      output_tmax[ch]->push_back(pt.t);
      output_rise[ch]->push_back(wm::CalcRiseTime(*w[ch], *t[ch], pt.index));
      output_area[ch]->push_back(wm::CalcPulseArea(*w[ch], *t[ch], pt.index));
      auto cfd = wm::CalcCFDTime(*w[ch], *t[ch], pt.index, 0.2, 0.5, 0.3);
      output_20cfd[ch]->push_back(cfd.at(0));
      output_50cfd[ch]->push_back(cfd.at(1));
      output_fwhm[ch]->push_back(wm::CalcFWHM(*w[ch], *t[ch], pt.index));
    }
    *output_nsignal[ch] = output_pmax[ch]->size();
    *output_rms[ch] = mix_params.rms;

    if(store_waveform){
      std::move(w[ch]->begin(), w[ch]->end(), std::back_inserter(*output_w[ch]));
      if(!output_t[ch]) continue;
      std::move(t[ch]->begin(), t[ch]->end(), std::back_inserter(*output_t[ch]));
    }
  }

}

void AnaSSRL::finalize(BetaConfigMgr *configMgr){
  // pass
}
