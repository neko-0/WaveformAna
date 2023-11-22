#include "AnaSSRL/AnaSSRL.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaSSRL::setup(BetaConfigMgr* const configMgr) { }

void AnaSSRL::initialize(BetaConfigMgr* const configMgr){

  LOG_INFO("external config: " + configMgr->ext_config_name());
  YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

  const auto &general = yaml_config["general"];
  this->store_waveform = general["store_waveform"].as<bool>();
  this->use_single_t_trace = general["use_single_t_trace"].as<bool>();
  this->use_single_input_t_trace = general["use_single_input_t_trace"].as<bool>();
  this->baseline_opt = general["baseline_opt"].as<int>();
  this->run_type = general["run_type"].as<int>();
  this->do_max_ch_ = general["do_max_ch"].as<bool>();
  this->threshold = general["threshold"].as<double>();
  int num_ch = this->num_ch_;
  if(general["nchannels"].as<int>() > 0){
    num_ch = general["nchannels"].as<int>();
  }
  this->invert_ch = general["invert_ch"].as<std::vector<int>>();
  this->trigger_ch = general["trigger_ch"].as<int>();
  this->simple_ana_ch = general["simple_ana_ch"].as<std::vector<int>>();

  const auto &buckets = yaml_config["buckets"];
  this->bucket_t_start_ = buckets["bucket_t_start"].as<double>();
  this->bucket_t_end_ = buckets["bucket_t_end"].as<double>();
  this->nbuckets_ = buckets["nbuckets"].as<int>();

  const auto &fix_win = yaml_config["fix_window"];
  this->fill_fix_window = fix_win["fill_fix_window"].as<bool>();
  this->fix_win_start_ = fix_win["fix_win_start"].as<double>();
  this->fix_win_step_size_ = fix_win["fix_win_step_size"].as<double>();
  this->fix_win_nstep_ = fix_win["fix_win_nstep"].as<int>();
  this->fix_win_edge_dist_ = fix_win["fix_win_edge_dist"].as<double>();

  for(int i = 0; i < num_ch; i++){
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

    thresholdTime[i] = configMgr->SetOutputBranch<double>("thresholdTime" + current_ch);

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

  LOG_INFO("number of active chanenls: " + std::to_string(active_ch_.size()));
}

//==============================================================================
void AnaSSRL::trigger_routine(std::vector<double> &corr_w, int ch){

  double polarity = 1.0;
  // check user specified invertion
  if(std::find(this->invert_ch.begin(), this->invert_ch.end(), ch) != this->invert_ch.end() ) {
    polarity = -1.0;
  }

  auto mix_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);

  corr_w.reserve(w[ch]->size());
  for(int i=0; i < w[ch]->size(); i++){
      corr_w.push_back((w[ch]->at(i) - mix_params.baseline) * polarity);
  }

  if(fill_fix_window){
    for(std::size_t _step = 0; _step < fix_win_nstep_; _step++){
      auto _begin = *t[ch]->begin();
      auto _end = *(t[ch]->end()-2);
      fill_fix_window_branches(ch, corr_w, *t[ch], _begin, _end, true);
    }
  }
}

//==============================================================================
void AnaSSRL::regular_routine(std::vector<double> &corr_w, int ch){
  double polarity = 1.0;
  double threshold;
  // check user specified invertion
  if(std::find(this->invert_ch.begin(), this->invert_ch.end(), ch) != this->invert_ch.end() ) {
    polarity = -1.0;
  }

  if(this->baseline_opt == 1) {
    std::vector<double> inv_w;
    inv_w.reserve(w[ch]->size());
    for(int i=0; i < w[ch]->size(); i++){
      inv_w.emplace_back(w[ch]->at(i) * -1.0); // need to handle positive signal
    }
    corr_w = wm::Baseline::ARPLS_PLS(inv_w, 1.0e5);
  } else {
    corr_w = wm::Baseline::NoiseMedian(*w[ch], w[ch]->size() / 12);
  }

  auto mix_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);
  // auto range_rms = wm::CalcNoise(*w[ch], *t[ch], 480, 550);
  // mix_params.rms = range_rms;

  if(this->run_type == 0){
    polarity = -1.0;
    // threshold = 30.0; // 5.0 * 6.0 AC Digitizer Run
    threshold = 5.0 * mix_params.rms;
  } else if (this->run_type == 1) {
    polarity = 1.0;
    threshold = 5.0 * mix_params.rms;
  } else {
    // check polarity of the signal. determine threshold for pmax finder
    if(std::abs(mix_params.pos_max) >= std::abs(mix_params.neg_max)) {
      polarity = 1.0;
    } else {
      polarity = -1.0;
    }
    threshold = 5.0 * mix_params.rms;
  }

  if(this->threshold > 0.0) threshold = this->threshold;

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

  if(*output_nsignal[ch] == 0){
    output_tmax_diff[ch]->push_back(0.0);
    output_pmax[ch]->push_back(0.0);
    output_tmax[ch]->push_back(0.0);
    output_rise[ch]->push_back(0.0);
    output_area[ch]->push_back(0.0);
    output_20cfd[ch]->push_back(0.0);
    output_50cfd[ch]->push_back(0.0);
    output_fwhm[ch]->push_back(0.0);
  }

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
}

//==============================================================================
void AnaSSRL::simple_routine(std::vector<double> &corr_w, int ch){
  double polarity = 1.0;
  // check user specified invertion
  if(std::find(this->invert_ch.begin(), this->invert_ch.end(), ch) != this->invert_ch.end() ) {
    polarity = -1.0;
  }
  corr_w.reserve(w[ch]->size());
  for(int i=0; i < w[ch]->size(); i++){
      corr_w.push_back(w[ch]->at(i) * polarity);
  }

  auto mix_params = wm::CalcMaxNoiseBase(corr_w, 0.25);
  *output_rms[ch] = mix_params.rms;

  if(fill_fix_window){
    for(std::size_t _step = 0; _step < fix_win_nstep_; _step++){
      auto _begin = *t[ch]->begin();
      auto _end = *(t[ch]->end()-2);
      fill_fix_window_branches(ch, corr_w, *t[ch], _begin, _end, true);
    }
  }
}

//==============================================================================
=======
>>>>>>> master
bool AnaSSRL::execute(BetaConfigMgr* const configMgr){
  // auto timestamp = wm::FindTimeAtThreshold(*trig, *t[0], 3000);
  // if(timestamp.size()>0) *trig_time = timestamp.at(0);
  // else *trig_time = -1.0;
  // #pragma omp parallel for if(active_ch_.size() > 5)
  for(auto &ch : active_ch_){
    if(w[ch]->size() == 0){
      LOG_WARNING("Trace size 0");
      continue;
    }

    std::vector<double> corr_w;

    if(std::find(this->simple_ana_ch.begin(), this->simple_ana_ch.end(), ch) != this->simple_ana_ch.end() ) {
      simple_routine(corr_w, ch);
    } else if(this->trigger_ch == ch){
      trigger_routine(corr_w, ch);
    } else {
      regular_routine(corr_w, ch);
    }

    auto th_time = wm::FindTimeAtThreshold(corr_w, *t[ch], 0.5);
    if(th_time.size() != 0){
      *thresholdTime[ch] = th_time.at(0);
    } else {
      *thresholdTime[ch] = -999;
    }

    if(store_waveform){
      std::move(corr_w.begin(), corr_w.end(), std::back_inserter(*output_corr_w[ch]));
      std::move(w[ch]->begin(), w[ch]->end(), std::back_inserter(*output_w[ch]));
      if(!output_t[ch]) continue;
      std::move(t[ch]->begin(), t[ch]->end(), std::back_inserter(*output_t[ch]));
    }
  }

  if(do_max_ch_) {
    std::vector<int> large_pad = {4, 5, 6, 7, 8, 9, 10, 11};
    std::vector<int> small_pad = {0, 1, 2, 3, 12, 13, 14, 15};
    std::vector<int> strip_set1 = {11, 4, 5};
    std::vector<int> strip_set2 = {10, 5, 11};
    std::vector<int> strip_set3 = {5, 4, 11, 10, 6};
    find_max_ch(large_pad, *output_max_ch, *(this->output_sum_large), large_pad);
    find_max_ch(large_pad, *output_2nd_max_ch, *(this->output_sum_large), large_pad, -1, -1, 0.8, 11);
    find_max_ch(small_pad, *output_small_pad_max_ch, *(this->output_sum_small));
    std::vector<int> tmp;
    find_max_ch(large_pad, tmp, *(this->output_sum_strip_set1), strip_set1, 11);
    find_max_ch(large_pad, tmp, *(this->output_sum_strip_set2), strip_set2, 5);
    find_max_ch(large_pad, tmp, *(this->output_sum_strip_set3), strip_set3, 5);
  }

  return true;
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
    std::string current_ch = std::to_string(ch_start_ + i);
    output_fix_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_pmax" + current_ch);
    output_fix_tmax[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_tmax" + current_ch);
    output_fix_area[i] = configMgr->SetOutputBranch<std::vector<float>>("fix_area" + current_ch);
  }

  if(!do_max_ch_) return;

  output_max_ch = configMgr->SetOutputBranch<std::vector<int>>("max_ch");
  output_2nd_max_ch = configMgr->SetOutputBranch<std::vector<int>>("max2_ch");
  output_small_pad_max_ch = configMgr->SetOutputBranch<std::vector<int>>("small_max_ch");
  output_sum_large = configMgr->SetOutputBranch<std::vector<double>>("sum_large");
  output_sum_small = configMgr->SetOutputBranch<std::vector<double>>("sum_small");
  output_sum_strip_set1 = configMgr->SetOutputBranch<std::vector<double>>("sum_strip_set1");
  output_sum_strip_set2 = configMgr->SetOutputBranch<std::vector<double>>("sum_strip_set2");
  output_sum_strip_set3 = configMgr->SetOutputBranch<std::vector<double>>("sum_strip_set3");
}

// =============================================================================
void AnaSSRL::fill_fix_window_branches(
  int ch,
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  double t_min,
  double t_max,
  bool fill_previous)
{
  if(output_fix_pmax[ch]->size()!=0 &&fill_previous){
      output_fix_pmax[ch]->push_back(output_fix_pmax[ch]->back());
      output_fix_tmax[ch]->push_back(output_fix_tmax[ch]->back());
      output_fix_area[ch]->push_back(output_fix_area[ch]->back());
  }

  double threshold = 0.0;
  auto s_max = wm::FindSignalMax(v_trace, t_trace, t_min, t_max);
  if(s_max.v < threshold
    || abs(s_max.t - t_min) <= this->fix_win_edge_dist_
    || abs(s_max.t - t_max) <= this->fix_win_edge_dist_) {
      output_fix_pmax[ch]->push_back(0.0);
      output_fix_tmax[ch]->push_back(-1.0);
      output_fix_area[ch]->push_back(0.0);
  }
  else{
    output_fix_pmax[ch]->push_back(s_max.v);
    output_fix_tmax[ch]->push_back(s_max.t);
    output_fix_area[ch]->push_back(wm::CalcPulseArea(v_trace, t_trace, s_max.index));
  }
}

// =============================================================================
void AnaSSRL::find_max_ch(
  const std::vector<int> &chlist,
  std::vector<int> &buffer,
  std::vector<double> &output,
  const std::vector<int> &sumCh,
  int targetCh,
  int targetCh2,
  double scale,
  int max_1st)
{
  int npmax = output_fix_pmax[active_ch_.at(0)]->size();
  for(std::size_t x = 0; x < npmax; x++){
    double v_max = -1.0;
    int ch_max = -1;
    for(auto &i : chlist) {
      if(i == max_1st) continue;
      if(output_fix_pmax[i]->at(x) > v_max){
        v_max = output_fix_pmax[i]->at(x);
        ch_max = i;
      }
    }
    buffer.push_back(ch_max);
    if( targetCh != -1 || targetCh2 != -1) {
      if( (ch_max != targetCh) && (ch_max != targetCh2) ) {
        output.push_back(-1.0);
        continue;
      }
    }
    if(ch_max != -1){
      bool do_sum = true;
      double summed = 0.0;
      for(auto &i : sumCh) {
        if(i == ch_max){
          summed+=v_max;
          continue;
        }
        if(output_fix_pmax[i]->at(x) > v_max*scale){
          do_sum = false;
          break;
        }
        if(abs(output_fix_tmax[ch_max]->at(x) - output_fix_tmax[i]->at(x)) > 3){
          do_sum = false;
          break;
        }
        summed += output_fix_pmax[i]->at(x);
      }
      if(do_sum) output.push_back(summed);
      else output.push_back(-1.0);
    }
  }
}

void AnaSSRL::find_max_ch(
  const std::vector<int> &chlist,
  std::vector<int> &buffer,
  std::vector<double> &output)
{
  AnaSSRL::find_max_ch(chlist, buffer, output, chlist, -1, -1);
}
