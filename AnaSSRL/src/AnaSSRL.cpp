#include "AnaSSRL/AnaSSRL.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

namespace wm = waveform_methods;

void AnaSSRL::initialize(BetaConfigMgr* const configMgr){

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
    output_rms[i] = configMgr->SetOutputBranch<float>("rms" + current_ch);
    output_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("pmax" + current_ch);
    output_tmax[i] = configMgr->SetOutputBranch<std::vector<float>>("tmax" + current_ch);
    output_rise[i] = configMgr->SetOutputBranch<std::vector<float>>("rise" + current_ch);
    output_area[i] = configMgr->SetOutputBranch<std::vector<float>>("area" + current_ch);
    output_fwhm[i] = configMgr->SetOutputBranch<std::vector<float>>("fwhm" + current_ch);
    output_20cfd[i] = configMgr->SetOutputBranch<std::vector<float>>("cfd20_" + current_ch);
    output_50cfd[i] = configMgr->SetOutputBranch<std::vector<float>>("cfd50_" + current_ch);
    output_tmax_diff[i] = configMgr->SetOutputBranch<std::vector<float>>("tmaxdiff" + current_ch);
    output_raw_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("raw_pmax" + current_ch);

    // output_bucket_corr[i] = configMgr->SetOutputBranch<float>("bucket_corr" + current_ch);
    output_bucket_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_pmax" + current_ch);
    output_bucket_tmax[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_tmax" + current_ch);
    output_bucket_area[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_area" + current_ch);
    output_bucket_cfd20[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_cfd20_" + current_ch);
    output_bucket_cfd50[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_cfd50_" + current_ch);
    output_bucket_index[i] = configMgr->SetOutputBranch<std::vector<int>>("bucket_index" + current_ch);
    output_bucket_tmax_diff[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_tmax_diff" + current_ch);
    output_bucket_cfd20_diff[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_cfd20_diff" + current_ch);
    output_bucket_cfd50_diff[i] = configMgr->SetOutputBranch<std::vector<float>>("bucket_cfd50_diff" + current_ch);
    output_bucket_pmax[i]->reserve(nbuckets_);
    output_bucket_tmax[i]->reserve(nbuckets_);
    output_bucket_area[i]->reserve(nbuckets_);
    output_bucket_cfd20[i]->reserve(nbuckets_);
    output_bucket_cfd50[i]->reserve(nbuckets_);
    output_bucket_index[i]->reserve(nbuckets_);
    output_bucket_tmax_diff[i]->reserve(nbuckets_);
    output_bucket_cfd20_diff[i]->reserve(nbuckets_);
    output_bucket_cfd50_diff[i]->reserve(nbuckets_);

    if(store_waveform){
      output_w[i] = configMgr->SetOutputBranch<std::vector<float>>("w" + current_ch);
      output_corr_w[i] = configMgr->SetOutputBranch<std::vector<float>>("w_corr" + current_ch);
      if(use_single_t_trace){
        if(found_single_t_trace){
          output_t[i] = nullptr;
          continue;
        }
        output_t[i] = configMgr->SetOutputBranch<std::vector<float>>("t");
        found_single_t_trace = true;
        continue;
      }
      output_t[i] = configMgr->SetOutputBranch<std::vector<float>>("t" + current_ch);
    }
  }

  LOG_INFO("number of active chanenls: " + std::to_string(active_ch_.size()));

  LOG_INFO("external config: " + configMgr->ext_config_name());
}

void AnaSSRL::execute(BetaConfigMgr* const configMgr){
  // #pragma omp parallel for
  for(auto &ch : active_ch_){
    if(w[ch]->size() == 0){
      LOG_WARNING("Trace size 0");
      continue;
    }

    // auto corr_w = wm::Baseline::ARPLS_PLS(*w[ch], 1.0e13);
    int win_size = w[ch]->size() / 12;
    auto corr_w = wm::Baseline::NoiseMedian(*w[ch], win_size);

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
      corr_w.at(i) *= polarity;
      // corr_w.at(i) = (corr_w.at(i)+w[ch]->at(i)) / 2.0; // look good with average
    }

    auto n_wave_pts = wm::FindMultipleSignalMax(corr_w, *t[ch], threshold);
    auto n_raw_wave_pts = wm::FindMultipleSignalMax(*w[ch], *t[ch], threshold);

    *output_basecorr[ch] = wm::Baseline::MultiSignalBaselineCorrection(
      n_raw_wave_pts, *w[ch], *t[ch], 0.5, 5e-9, 5e-9);

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
      output_rise[ch]->push_back(wm::CalcRiseTime(corr_w, *t[ch], pt.index));
      output_area[ch]->push_back(wm::CalcPulseArea(corr_w, *t[ch], pt.index));
      auto cfd = wm::CalcCFDTime(corr_w, *t[ch], pt.index, 0.2, 0.5, 0.3);
      output_20cfd[ch]->push_back(cfd.at(0));
      output_50cfd[ch]->push_back(cfd.at(1));
      output_fwhm[ch]->push_back(wm::CalcFWHM(corr_w, *t[ch], pt.index));
    }
    *output_nsignal[ch] = output_pmax[ch]->size();
    *output_rms[ch] = mix_params.rms;

    for(auto &pt : n_raw_wave_pts){
      if(pt.index < 0) continue;
      output_raw_pmax[ch]->push_back(pt.v);
    }

    bucket_time_difference(ch, -44e-9, 2.08e-9);

    if(store_waveform){
      std::move(corr_w.begin(), corr_w.end(), std::back_inserter(*output_corr_w[ch]));
      std::move(w[ch]->begin(), w[ch]->end(), std::back_inserter(*output_w[ch]));
      if(!output_t[ch]) continue;
      std::move(t[ch]->begin(), t[ch]->end(), std::back_inserter(*output_t[ch]));
    }
  }

}

void AnaSSRL::finalize(BetaConfigMgr* const configMgr){
  // pass
}

void AnaSSRL::bucket_time_difference(
  const int &ch,
  const double &bucket_start,
  const double &bucket_step)
{
  int num_filled = output_tmax[ch]->size();
  for(int i = 0; i < nbuckets_; i++){
    double lbound = bucket_start + bucket_step * i;
    double hbound = lbound + bucket_step;
    int nfilled = 0;
    int fill_index = -1;
    for(int j = 0; j < num_filled; j++){
      auto tmax = output_tmax[ch]->at(j);
      // auto tmax = output_50cfd[ch]->at(j);
      if(tmax < lbound || tmax > hbound) continue;
      if(fill_index == -1){
        fill_index = j;
      } else {
        fill_index = -1;
        break;
      }
    }
    if(fill_index == -1){
      output_bucket_pmax[ch]->emplace_back(0.0);
      output_bucket_tmax[ch]->emplace_back(0.0);
      output_bucket_area[ch]->emplace_back(0.0);
      output_bucket_cfd20[ch]->emplace_back(0.0);
      output_bucket_cfd50[ch]->emplace_back(0.0);
      output_bucket_index[ch]->emplace_back(-1.0);
    } else{
      output_bucket_tmax[ch]->emplace_back(output_tmax[ch]->at(fill_index));
      output_bucket_pmax[ch]->emplace_back(output_pmax[ch]->at(fill_index));
      output_bucket_area[ch]->emplace_back(output_area[ch]->at(fill_index));
      output_bucket_cfd20[ch]->emplace_back(output_20cfd[ch]->at(fill_index));
      output_bucket_cfd50[ch]->emplace_back(output_50cfd[ch]->at(fill_index));
      output_bucket_index[ch]->emplace_back(i+1);
    }
  }

  // calculate difference
  for(int i = 0; i < nbuckets_; i++){
    if(i==0){
      output_bucket_tmax_diff[ch]->emplace_back(-1.0);
      output_bucket_cfd20_diff[ch]->emplace_back(-1.0);
      output_bucket_cfd50_diff[ch]->emplace_back(-1.0);
      continue;
    }
    auto bucket_i_1 = output_bucket_index[ch]->at(i);
    auto bucket_i_2 = output_bucket_index[ch]->at(i-1);
    if(bucket_i_1 == -1 || bucket_i_2 == -1){
      output_bucket_tmax_diff[ch]->emplace_back(-1.0);
      output_bucket_cfd20_diff[ch]->emplace_back(-1.0);
      output_bucket_cfd50_diff[ch]->emplace_back(-1.0);
    } else {
      int ndiff = bucket_i_1 - bucket_i_2;
      float tmax_diff = (output_bucket_tmax[ch]->at(i) - output_bucket_tmax[ch]->at(i-1)) / ndiff;
      float cfd20_diff = (output_bucket_cfd20[ch]->at(i) - output_bucket_cfd20[ch]->at(i-1)) / ndiff;
      float cfd50_diff = (output_bucket_cfd50[ch]->at(i) - output_bucket_cfd50[ch]->at(i-1)) / ndiff;
      output_bucket_tmax_diff[ch]->emplace_back(tmax_diff);
      output_bucket_cfd20_diff[ch]->emplace_back(cfd20_diff);
      output_bucket_cfd50_diff[ch]->emplace_back(cfd50_diff);
    }
  }
}
