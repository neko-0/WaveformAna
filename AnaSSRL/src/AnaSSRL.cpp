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
  this->ch_start_ = general["ch_start"].as<int>();
  this->store_waveform = general["store_waveform"].as<bool>();
  this->use_single_t_trace = general["use_single_t_trace"].as<bool>();
  this->use_single_input_t_trace = general["use_single_input_t_trace"].as<bool>();
  this->baseline_opt = general["baseline_opt"].as<int>();
  this->run_type = general["run_type"].as<int>();
  this->do_max_ch_ = general["do_max_ch"].as<bool>();
  this->threshold_ = general["threshold"].as<double>();
  int num_ch = this->num_ch_;
  if(general["nchannels"].as<int>() > 0){
    num_ch = general["nchannels"].as<int>();
  }
  this->invert_ch = general["invert_ch"].as<std::vector<int>>();
  this->trigger_ch = general["trigger_ch"].as<int>();
  this->simple_ana_ch = general["simple_ana_ch"].as<std::vector<int>>();
  this->routine_ = general["routine"].as<int>();
  this->rms_start_ = general["rms_start"].as<double>();
  this->rms_end_ = general["rms_end"].as<double>();
  
  const auto &bunch_win = yaml_config["bunch_window"];
  this->bunch_start_ = bunch_win["bunch_start"].as<double>();
  this->bunch_step_size_ = bunch_win["bunch_step_size"].as<double>();
  this->bunch_nstep_ = bunch_win["bunch_nstep"].as<int>();
  this->bunch_edge_dist_ = bunch_win["bunch_edge_dist"].as<double>();

  const auto &leading = yaml_config["leading_signal"];
  this->leading_tmin_ = leading["tmin"].as<double>();
  this->leading_tmax_ = leading["tmax"].as<double>();

  // CAEN trigger channels
  trg0 = configMgr->SetInputBranch<std::vector<double>>("trg0");
  trg1 = configMgr->SetInputBranch<std::vector<double>>("trg1");
  
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
    thresholdTime[i] = configMgr->SetOutputBranch<double>("thresholdTime" + current_ch);
    output_rms[i] = configMgr->SetOutputBranch<float>("rms" + current_ch);

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

  prepare_bunch_window_branches(configMgr);
  prepare_leading_signal_branches(configMgr);

  // branch for scan routine
  if(routine_ == 2){
    pos = configMgr->SetInputBranch<std::vector<double>>("pos");
    output_x = configMgr->SetOutputBranch<double>("posx");
    output_y = configMgr->SetOutputBranch<double>("posy");
  }
  

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

  for(std::size_t _step = 0; _step < bunch_nstep_; _step++){
    auto _begin = *t[ch]->begin();
    auto _end = *(t[ch]->end()-2);
    fill_bunch_window_branches(ch, corr_w, *t[ch], _begin, _end, true);
  }
  
}

//==============================================================================
void AnaSSRL::regular_routine(std::vector<double> &corr_w, int ch){
  double polarity = 1.0;
  double threshold = 0.0;

  auto mix_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);
  *output_rms[ch] = wm::CalcNoise(*w[ch], *t[ch], this->rms_start_, this->rms_end_);
  
  if(this->run_type == 0) { // auto polarity detection
    // check polarity of the signal. determine threshold for pmax finder
    if(std::abs(mix_params.pos_max) >= std::abs(mix_params.neg_max)) polarity = 1.0;
    else polarity = -1.0;
  } 
  
  // overwrite by user threshold
  if(this->threshold_ > 0.0) threshold = this->threshold_;
  else threshold = *output_rms[ch] * 5.0;
  
  std::shared_ptr<std::vector<double>> inv_w = nullptr;
  // check user specified invertion
  if(std::find(this->invert_ch.begin(), this->invert_ch.end(), ch) != this->invert_ch.end() ) {
    inv_w = std::make_shared<std::vector<double>>();
    inv_w->reserve(w[ch]->size());
    for(int i=0; i < w[ch]->size(); i++) {
      inv_w->emplace_back(w[ch]->at(i) * -1.0);
    }
    polarity = 1.0; // reset to 1.0
  } else {
    inv_w = std::shared_ptr<std::vector<double>>(w[ch], [](std::vector<double>*) {});
  }

  if(this->baseline_opt == 1) {
    corr_w = wm::Baseline::ARPLS_PLS(*inv_w, 1.0e5);
  } else {
    corr_w = wm::Baseline::NoiseMedian(*inv_w, inv_w->size() / 12);
  }

  for(int i=0; i < w[ch]->size(); i++) {
    w[ch]->at(i) = (w[ch]->at(i) - mix_params.baseline) * polarity;
    corr_w.at(i) *= polarity;
  }

  // search for all possible Pmax
  // auto n_wave_pts = wm::FindMultipleSignalMax(corr_w, *t[ch], threshold);
  // auto n_raw_wave_pts = wm::FindMultipleSignalMax(w[ch], t[ch], threshold);

  // *output_basecorr[ch] = wm::Baseline::MultiSignalBaselineCorrection(
  //   n_raw_wave_pts, w[ch], t[ch], 0.5, 5e-9, 5e-9
  // );

  // std::move(cfd_times.begin(), cfd_times.end(), std::back_inserter(*output_cfd[ch]));
  // #pragma omp parallel for
  // for(int pt_i=0; pt_i < n_wave_pts.size(); pt_i++){
  //   auto pt = n_wave_pts.at(pt_i);

  fill_leading_signal_branches(ch, corr_w, *t[ch], this->leading_tmin_, this->leading_tmax_);
  
  for(std::size_t _step = 0; _step < bunch_nstep_; _step++){
    auto _begin = bunch_start_ + _step * bunch_step_size_;
    auto _end = _begin + bunch_step_size_;
    fill_bunch_window_branches(ch, corr_w, *t[ch], _begin, _end);
  }
  fill_bunch_window_wp(ch);

}

//==============================================================================
/* 
Simple routine to with baseline correction using only 25% points.
This routine does NOT distinguish bunches, and treat entire waveform as single signal waveform.
*/
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

  fill_leading_signal_branches(ch, corr_w, *t[ch], this->leading_tmin_, this->leading_tmax_);

  // baseline correction with first 25% of the points.
  auto mix_params = wm::CalcMaxNoiseBase(corr_w, 0.25);
  *output_rms[ch] = mix_params.rms;


  for(std::size_t _step = 0; _step < bunch_nstep_; _step++){
    auto _begin = *t[ch]->begin();
    auto _end = *(t[ch]->end()-2);
    // The true flag just repeatly fill with the first entry. good for traces without bunches.
    fill_bunch_window_branches(ch, corr_w, *t[ch], _begin, _end, true);
  }
  
}

//==============================================================================
void AnaSSRL::scan_routinue(std::vector<double> &corr_w, int ch) {
  // perform time shift correction with the trigger.
  if(output_t[ch] || corr_common_time_.empty()){
    corr_common_time_.clear();
    corr_common_time_.reserve(t[ch]->size());
    for(int i = 0; i < t[ch]->size(); i++) {
      corr_common_time_.push_back(t[ch]->at(i) - trg_threshold_time_);
    }
  }

  // polarity invertion
  double polarity = 1.0;
  // check user specified invertion
  if(std::find(this->invert_ch.begin(), this->invert_ch.end(), ch) != this->invert_ch.end() ) {
    polarity = -1.0;
  }

  // Make sure the signal is positive for ARPLS baseline correction.
  std::vector<double> inv_w(w[ch]->size(), 0.0);
  for(int i=0; i < w[ch]->size(); i++){
    inv_w[i] = w[ch]->at(i) * polarity;
  }
  // corr_w = wm::Baseline::ARPLS_PLS(inv_w, 1.0e5);
  corr_w = wm::Baseline::NoiseMedian(inv_w, inv_w.size() / 12);

  // need to estimate the noise. 
  // use the gap between bunch leading signal and sensor response
  // currently hard code from [300, 500] * 200sp
  double rms = wm::CalcNoise(corr_w, corr_common_time_, this->rms_start_, this->rms_end_);
  *output_rms[ch] = rms;

  // search for all possible Pmax
  double threshold = this->threshold_;
  if(threshold <= 0.0) threshold = rms * 5;

  fill_leading_signal_branches(ch, corr_w, corr_common_time_, this->leading_tmin_, this->leading_tmax_);
  
  for(std::size_t _step = 0; _step < bunch_nstep_; _step++){
    auto _begin = bunch_start_ + _step * bunch_step_size_;
    auto _end = _begin + bunch_step_size_;
    bool _fill_previous = false;
    if(_end >= corr_common_time_.size() ){
      _fill_previous = true;
    }
    fill_bunch_window_branches(ch, corr_w, corr_common_time_, _begin, _end, _fill_previous);
  }
  fill_bunch_window_wp(ch);

  if(output_t[ch]) {
    std::move(corr_common_time_.begin(), corr_common_time_.end(), std::back_inserter(*output_t[ch]));
  }
  
}

//==============================================================================
bool AnaSSRL::execute(BetaConfigMgr* const configMgr){

  if( trg0 ) {
    trg_threshold_time_ = wm::FindTimeAtThreshold(*trg0, *t[0], this->threshold_, true).at(0);
  }
  // auto timestamp = wm::FindTimeAtThreshold(*trig, *t[0], 3000);
  // if(timestamp.size()>0) *trig_time = timestamp.at(0);
  // else *trig_time = -1.0;
  // #pragma omp parallel for if(active_ch_.size() > 5)
  for(auto &ch : active_ch_) {
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
      if(routine_ == 2) {
        scan_routinue(corr_w, ch);
      }
      else {
        regular_routine(corr_w, ch);
      }
    }

    *thresholdTime[ch] = wm::FindTimeAtThreshold(corr_w, *t[ch], this->threshold_, true).at(0);
   
    if(store_waveform){
      std::move(corr_w.begin(), corr_w.end(), std::back_inserter(*output_corr_w[ch]));
      std::move(w[ch]->begin(), w[ch]->end(), std::back_inserter(*output_w[ch]));
      if(routine_ == 2) continue; // skip time for scan routine here
      if(!output_t[ch]) continue;
      std::move(t[ch]->begin(), t[ch]->end(), std::back_inserter(*output_t[ch]));
    }
  }

  // fill position for scan routine
  if(routine_ == 2) {
    *output_x = pos->at(0);
    *output_y = pos->at(1);
  }

  // you probably don't need this except for previous old AC-LGAD runs.
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

// =============================================================================
void AnaSSRL::prepare_leading_signal_branches(BetaConfigMgr* const configMgr){
  for(auto &i : active_ch_){
    std::string current_ch = std::to_string(ch_start_ + i);
    output_leading_pmax[i] = configMgr->SetOutputBranch<float>("leading_pmax" + current_ch);
    output_leading_tmax[i] = configMgr->SetOutputBranch<float>("leading_tmax" + current_ch);
    output_leading_area[i] = configMgr->SetOutputBranch<float>("leading_area" + current_ch);
    output_leading_rise[i] = configMgr->SetOutputBranch<float>("leading_risetime" + current_ch);
    output_leading_fall[i] = configMgr->SetOutputBranch<float>("leading_falltime" + current_ch);
    output_leading_20cfd[i] = configMgr->SetOutputBranch<float>("leading_cfd20_" + current_ch);
    output_leading_50cfd[i] = configMgr->SetOutputBranch<float>("leading_cfd50_" + current_ch);
  }
}

// =============================================================================
void AnaSSRL::prepare_bunch_window_branches(BetaConfigMgr* const configMgr){
  for(auto &i : active_ch_){
    std::string current_ch = std::to_string(ch_start_ + i);
    output_pmax[i] = configMgr->SetOutputBranch<std::vector<float>>("pmax" + current_ch);
    output_tmax[i] = configMgr->SetOutputBranch<std::vector<float>>("tmax" + current_ch);
    output_area[i] = configMgr->SetOutputBranch<std::vector<float>>("area" + current_ch);
    output_rise[i] = configMgr->SetOutputBranch<std::vector<float>>("risetime" + current_ch);
    output_fall[i] = configMgr->SetOutputBranch<std::vector<float>>("falltime" + current_ch);
    output_20cfd[i] = configMgr->SetOutputBranch<std::vector<float>>("cfd20_" + current_ch);
    output_50cfd[i] = configMgr->SetOutputBranch<std::vector<float>>("cfd50_" + current_ch);
    
    output_rms_wp_loose[i] = configMgr->SetOutputBranch<std::vector<bool>>("rms_wp_loose" + current_ch);
    output_rms_wp_tight[i] = configMgr->SetOutputBranch<std::vector<bool>>("rms_wp_tight" + current_ch);
    output_bunch_wp_loose[i] = configMgr->SetOutputBranch<std::vector<bool>>("bunch_wp_loose" + current_ch);
    output_bunch_wp_tight[i] = configMgr->SetOutputBranch<std::vector<bool>>("bunch_wp_tight" + current_ch);
    output_fall_wp_loose[i] = configMgr->SetOutputBranch<std::vector<bool>>("fall_wp_loose" + current_ch);
    output_fall_wp_tight[i] = configMgr->SetOutputBranch<std::vector<bool>>("fall_wp_tight" + current_ch);
    output_wp_loose[i] = configMgr->SetOutputBranch<std::vector<bool>>("wp_loose" + current_ch);
    output_wp_tight[i] = configMgr->SetOutputBranch<std::vector<bool>>("wp_tight" + current_ch);

    output_pmax[i]->reserve(bunch_nstep_);
    output_tmax[i]->reserve(bunch_nstep_);
    output_area[i]->reserve(bunch_nstep_);
    output_rise[i]->reserve(bunch_nstep_);
    output_fall[i]->reserve(bunch_nstep_);
    output_20cfd[i]->reserve(bunch_nstep_);
    output_50cfd[i]->reserve(bunch_nstep_);
    output_rms_wp_loose[i]->reserve(bunch_nstep_);
    output_rms_wp_tight[i]->reserve(bunch_nstep_);
    output_bunch_wp_loose[i]->reserve(bunch_nstep_);
    output_bunch_wp_tight[i]->reserve(bunch_nstep_);
    output_fall_wp_loose[i]->reserve(bunch_nstep_);
    output_fall_wp_tight[i]->reserve(bunch_nstep_);
    output_wp_loose[i]->reserve(bunch_nstep_);
    output_wp_tight[i]->reserve(bunch_nstep_);
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
void AnaSSRL::fill_leading_signal_branches(
  int ch,
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  double t_min,
  double t_max)
{
  auto s_max = wm::FindSignalMax(v_trace, t_trace, t_min, t_max);
  *output_leading_pmax[ch] = s_max.v;
  *output_leading_tmax[ch] = s_max.t;
  *output_leading_area[ch] = wm::CalcPulseArea(v_trace, t_trace, s_max.index);
  *output_leading_fall[ch] = wm::CalcFallTime(v_trace, t_trace, s_max.index, 0.1, 0.9);
  *output_leading_rise[ch] = wm::CalcRiseTime(v_trace, t_trace, s_max.index, 0.1, 0.9);
  
  auto cfd = wm::CalcCFDTime(v_trace, t_trace, s_max.index, 0.2, 0.5, 0.3);
  *output_leading_20cfd[ch] = cfd.at(0);
  *output_leading_50cfd[ch] = cfd.at(1);
}

// =============================================================================
/*
Find signal maximum in a given time interval.
The default threshold for accepting is 0.0
*/
void AnaSSRL::fill_bunch_window_branches(
  int ch,
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  double t_min,
  double t_max,
  bool fill_previous,
  double threshold)
{
  // just repeatly fill with the first entry if fill_previous = true.
  if(output_pmax[ch]->size()!=0 && fill_previous){
      output_pmax[ch]->push_back(output_pmax[ch]->back());
      output_tmax[ch]->push_back(output_tmax[ch]->back());
      output_area[ch]->push_back(output_area[ch]->back());
      output_rise[ch]->push_back(output_rise[ch]->back());
      output_fall[ch]->push_back(output_fall[ch]->back());
      output_20cfd[ch]->push_back(output_20cfd[ch]->back());
      output_50cfd[ch]->push_back(output_50cfd[ch]->back());
  } else {
    auto s_max = wm::FindSignalMax(v_trace, t_trace, t_min, t_max);
    // do check on Pmax value, and distance of Tmax from the left and right time interval bounds. 
    if(s_max.v < threshold
      || abs(s_max.t - t_min) <= this->bunch_edge_dist_
      || abs(s_max.t - t_max) <= this->bunch_edge_dist_) {
        output_pmax[ch]->push_back(-1.0);
        output_tmax[ch]->push_back(-1.0);
        output_area[ch]->push_back(0.0);
        output_rise[ch]->push_back(0.0);
        output_fall[ch]->push_back(0.0);
        output_20cfd[ch]->push_back(0.0);
        output_50cfd[ch]->push_back(0.0);
    }
    else{
      output_pmax[ch]->push_back(s_max.v);
      output_tmax[ch]->push_back(s_max.t);
      output_area[ch]->push_back(wm::CalcPulseArea(v_trace, t_trace, s_max.index));
      output_rise[ch]->push_back(wm::CalcRiseTime(v_trace, t_trace, s_max.index, 0.1, 0.9));
      output_fall[ch]->push_back(wm::CalcFallTime(v_trace, t_trace, s_max.index, 0.1, 0.9));
      
      auto cfd = wm::CalcCFDTime(v_trace, t_trace, s_max.index, 0.2, 0.5, 0.3);
      output_20cfd[ch]->push_back(cfd.at(0));
      output_50cfd[ch]->push_back(cfd.at(1));
    }
  }
}

// ===============================================
void AnaSSRL::fill_bunch_window_wp(int ch)
{
  // check for working points
  // loose working point requires pmax > 3*rms, and tmax difference > 0.5*fix_win_step_size_
  // tight working point requires pmax > 5*rms, and tmax difference > 0.8*fix_win_step_size_
  
  std::vector<int> bad_rms_loose;
  std::vector<int> bad_rms_tight;
  std::vector<int> bad_bunch_loose;
  std::vector<int> bad_bunch_tight;
  std::vector<int> bad_fall_loose;
  std::vector<int> bad_fall_tight;
  
  for(auto i = 0; i < output_pmax[ch]->size()-1; i++) {
    double pmax1 = output_pmax[ch]->at(i);
    double tmax1 = output_tmax[ch]->at(i);
    double pmax2 = output_pmax[ch]->at(i+1);
    double tmax2 = output_tmax[ch]->at(i+1);
    
    // loose
    if( pmax1 < (*output_rms[ch] * 3)) bad_rms_loose.push_back(i);
    if( pmax2 < (*output_rms[ch] * 3)) bad_rms_loose.push_back(i+1);
  
    if(tmax2 - tmax1 < 0.2 * bunch_step_size_) {
      if(pmax1 < pmax2) bad_bunch_loose.push_back(i);
      else bad_bunch_loose.push_back(i+1);
    } 
    
    if(tmax2 < tmax1 + 0.5 * output_fall[ch]->at(i)) bad_fall_loose.push_back(i+1);
    
    // tight
    if( pmax1 < (*output_rms[ch] * 5)) bad_rms_tight.push_back(i);
    if( pmax2 < (*output_rms[ch] * 5)) bad_rms_tight.push_back(i+1);
    
    if(tmax2 - tmax1 < 0.5 * bunch_step_size_) {
      if(pmax1 < pmax2) bad_bunch_tight.push_back(i);
      else bad_bunch_tight.push_back(i+1);
    } 
    
    if(tmax2 < tmax1 + 0.8 * output_fall[ch]->at(i)) bad_fall_tight.push_back(i+1);

  }

  // filling wp
  for(auto i = 0; i < output_pmax[ch]->size(); i++) {
    // filling loose
    if(std::find(bad_rms_loose.begin(), bad_rms_loose.end(), i) != bad_rms_loose.end()){
      output_rms_wp_loose[ch]->push_back(false);
    } else {
      output_rms_wp_loose[ch]->push_back(true);
    }

    if(std::find(bad_bunch_loose.begin(), bad_bunch_loose.end(), i) != bad_bunch_loose.end()){
      output_bunch_wp_loose[ch]->push_back(false);
    } else {
      output_bunch_wp_loose[ch]->push_back(true);
    }

    if(std::find(bad_fall_loose.begin(), bad_fall_loose.end(), i) != bad_fall_loose.end()){
      output_fall_wp_loose[ch]->push_back(false);
    } else {
      output_fall_wp_loose[ch]->push_back(true);
    }

    if(output_rms_wp_loose[ch]->back() && output_bunch_wp_loose[ch]->back() && output_fall_wp_loose[ch]->back()) {
      output_wp_loose[ch]->push_back(true);
    } else {
      output_wp_loose[ch]->push_back(false);
    }

    // filling tight
    if(std::find(bad_rms_tight.begin(), bad_rms_tight.end(), i) != bad_rms_tight.end()){
      output_rms_wp_tight[ch]->push_back(false);
    } else {
      output_rms_wp_tight[ch]->push_back(true);
    }

    if(std::find(bad_bunch_tight.begin(), bad_bunch_tight.end(), i) != bad_bunch_tight.end()){
      output_bunch_wp_tight[ch]->push_back(false);
    } else {
      output_bunch_wp_tight[ch]->push_back(true);
    }

    if(std::find(bad_fall_tight.begin(), bad_fall_tight.end(), i) != bad_fall_tight.end()){
      output_fall_wp_tight[ch]->push_back(false);
    } else {
      output_fall_wp_tight[ch]->push_back(true);
    }

    if(output_rms_wp_tight[ch]->back() && output_bunch_wp_tight[ch]->back() && output_fall_wp_tight[ch]->back()) {
      output_wp_tight[ch]->push_back(true);
    } else {
      output_wp_tight[ch]->push_back(false);
    }

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
  int npmax = output_pmax[active_ch_.at(0)]->size();
  for(std::size_t x = 0; x < npmax; x++){
    double v_max = -1.0;
    int ch_max = -1;
    for(auto &i : chlist) {
      if(i == max_1st) continue;
      if(output_pmax[i]->at(x) > v_max){
        v_max = output_pmax[i]->at(x);
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
        if(output_pmax[i]->at(x) > v_max*scale){
          do_sum = false;
          break;
        }
        if(abs(output_tmax[ch_max]->at(x) - output_tmax[i]->at(x)) > 3){
          do_sum = false;
          break;
        }
        summed += output_pmax[i]->at(x);
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
