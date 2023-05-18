#include "AnaSSRL/AnaSSRL.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaSSRL::initialize(BetaConfigMgr* const configMgr){

  LOG_INFO("external config: " + configMgr->ext_config_name());
  YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

  const auto &general = yaml_config["general"];
  this->store_waveform = general["store_waveform"].as<bool>();
  this->use_single_t_trace = general["use_single_t_trace"].as<bool>();
  this->use_single_input_t_trace = general["use_single_input_t_trace"].as<bool>();
  int num_ch = this->num_ch_;
  if(general["nchannels"].as<int>() > 0){
    num_ch = general["nchannels"].as<int>();
  }

  const auto &buckets = yaml_config["buckets"];
  this->bucket_t_start_ = buckets["bucket_t_start"].as<double>();
  this->bucket_t_end_ = buckets["bucket_t_end"].as<double>();
  this->nbuckets_ = buckets["nbuckets"].as<int>();

  const auto &fix_win = yaml_config["fix_window"];
  this->fill_fix_window = fix_win["fill_fix_window"].as<bool>();
  this->fix_win_start_ = fix_win["fix_win_start"].as<double>();
  this->fix_win_step_size_ = fix_win["fix_win_step_size"].as<double>();
  this->fix_win_nstep_ = fix_win["fix_win_nstep"].as<int>();

  for(int i = 0; i < num_ch_; i++){
    // input branches
    std::string current_ch = std::to_string(ch_start_ + i);
    w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
    if(!w[i]) continue;
    if(!use_single_input_t_trace){
      t[i] = configMgr->SetInputBranch<std::vector<double>>("t" + current_ch);
    } else{
        if(i==0) t[i] = configMgr->SetInputBranch<std::vector<double>>("time");
        else t[i] = t[0];
    }

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

  if(fill_fix_window) prepare_fix_window_branches(configMgr);

  // trig = configMgr->SetInputBranch<std::vector<double>>("trg0");
  // trig_time = configMgr->SetOutputBranch<double>("trig_time");

  LOG_INFO("number of active chanenls: " + std::to_string(active_ch_.size()));
}

void AnaSSRL::execute(BetaConfigMgr* const configMgr){
  // auto timestamp = wm::FindTimeAtThreshold(*trig, *t[0], 3000);
  // if(timestamp.size()>0) *trig_time = timestamp.at(0);
  // else *trig_time = -1.0;

  // #pragma omp parallel for if(active_ch_.size() > 5)
  for(auto &ch : active_ch_){
    if(w[ch]->size() == 0){
      LOG_WARNING("Trace size 0");
      continue;
    }

    // auto corr_w = wm::Baseline::ARPLS_PLS(*w[ch], 1.0e13);
    auto corr_w = wm::Baseline::NoiseMedian(*w[ch], w[ch]->size() / 12);

    auto mix_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);
    // auto range_rms = wm::CalcNoise(*w[ch], *t[ch], 480, 550);
    // mix_params.rms = range_rms;

    // check polarity of the signal. determine threshold for pmax finder
    double polarity = -1.0; // forcing signal inverstion
    // if(std::abs(mix_params.pos_max) >= std::abs(mix_params.neg_max)) {
    //   polarity = 1.0;
    // } else {
    //   polarity = -1.0;
    // }
    // double threshold = 5.0 * mix_params.rms;
    double threshold = 5.0 * 6.0; // AC Digitizer Run

    for(int i=0; i < w[ch]->size(); i++){
      w[ch]->at(i) = (w[ch]->at(i) - mix_params.baseline)*polarity;
      corr_w.at(i) *= polarity;
    }

    auto n_wave_pts = wm::FindMultipleSignalMax(corr_w, *t[ch], threshold);
    auto n_raw_wave_pts = wm::FindMultipleSignalMax(w[ch], t[ch], threshold);

    *output_basecorr[ch] = wm::Baseline::MultiSignalBaselineCorrection(
      n_raw_wave_pts, w[ch], t[ch], 0.5, 5e-9, 5e-9
    );

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

    bucket_time_difference(ch, bucket_t_start_, bucket_t_end_); // 0.5, 1.0 for Digitizer run

    if(fill_fix_window){
      for(std::size_t _step = 0; _step < fix_win_nstep_; _step++){
        auto _begin = fix_win_start_ + _step * fix_win_step_size_;
        auto _end = _begin + fix_win_step_size_;
        fill_fix_window_branches(ch, corr_w, *t[ch], _begin, _end);
      }
    }

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

// =============================================================================
void AnaSSRL::prepare_fix_window_branches(BetaConfigMgr* const configMgr){
  for(auto &i : active_ch_){
    std::string current_ch = std::to_string(i+1);
    output_fix_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_pmax" + current_ch);
    output_fix_tmax[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_tmax" + current_ch);
    output_fix_area[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_area" + current_ch);
  }
}

// =============================================================================
void AnaSSRL::fill_fix_window_branches(
  int ch,
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  double t_min,
  double t_max)
{
  double threshold = 0.0;
  auto s_max = wm::FindSignalMax(v_trace, t_trace, t_min, t_max);
  if(s_max.v > threshold){
    output_fix_pmax[ch]->push_back(s_max.v);
    output_fix_tmax[ch]->push_back(s_max.t);
    output_fix_area[ch]->push_back(wm::CalcPulseArea(v_trace, t_trace, s_max.index));
  }
  else{
    output_fix_pmax[ch]->push_back(0.0);
    output_fix_tmax[ch]->push_back(-1.0);
    output_fix_area[ch]->push_back(0.0);
  }
}
