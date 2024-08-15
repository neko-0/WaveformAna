#include "AnaHighBW/AnaHighBW.hpp"
#include "waveformMethods/core.hpp"
#include "waveformMethods/filter.hpp"

#include <math.h>
#include <numeric>
#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaHighBW::setup(BetaConfigMgr* const configMgr) { }

void AnaHighBW::initialize(BetaConfigMgr* const configMgr) {

    LOG_INFO("external config: " + configMgr->ext_config_name());
    
    // parameters are all passed by the YAML file.
    YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

    const auto &general = yaml_config["general"];
    enbale_ch_ = general["enable_channels"].as<std::vector<int>>();
    dut_ch_ = general["dut_channel"].as<int>();
    clock_ch_ = general["clock_channel"].as<int>();
    winsize_ = general["winsize"].as<double>();
    peak_finding_threshold_ = general["peak_finding_threshold"].as<double>();

    const auto &cap_param = yaml_config["capacitor_param"];
    cap_t_start_ = cap_param["time_start"].as<int>();
    cap_t_win_ = cap_param["time_win"].as<int>();
    cap_cycle_size_ = cap_param["cycle_size"].as<int>();
    cap_ncycles_ = cap_param["ncycles"].as<int>();
    
    // prepare input and output voltage and time branches.
    for(int i = 0; i < enbale_ch_.size(); i++) {
        std::string current_ch = std::to_string(ch_start_ + enbale_ch_.at(i));
        i_w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
        i_t[i] = configMgr->SetInputBranch<std::vector<double>>("t" + current_ch);

        o_w[i] = configMgr->SetOutputBranch<std::vector<double>>("w" + current_ch);
        o_t[i] = configMgr->SetOutputBranch<std::vector<double>>("t" + current_ch);

        active_ch_.push_back(i);
    }

    // storing processed voltage-time traces.
    // dut_process_time = configMgr->SetOutputBranch<std::vector<double>>("dut_p_t");
    // dut_process_voltage = configMgr->SetOutputBranch<std::vector<double>>("dut_p_w");

    // estimate of the voltage and time at each capacitor step.
    w_step = configMgr->SetOutputBranch<std::vector<double>>("w_step");
    t_step = configMgr->SetOutputBranch<std::vector<double>>("t_step");
}

bool AnaHighBW::execute(BetaConfigMgr* const configMgr){
    for(auto &ch : active_ch_) {
        if(i_w[ch]->size() == 0){
            LOG_WARNING("Trace size 0");
            continue;
        }
    }
    
    // find_capacitors(dut_ch_);
    // find_capacitors(dut_ch_, clock_ch_);
    find_capacitors(dut_ch_, cap_t_start_, cap_t_win_, cap_cycle_size_, cap_ncycles_);
    
    // dumping traces to output.
    for(auto &ch : active_ch_) {
        std::move(i_t[ch]->begin(), i_t[ch]->end(), std::back_inserter(*o_t[ch]));
        std::move(i_w[ch]->begin(), i_w[ch]->end(), std::back_inserter(*o_w[ch]));
    }

    return true;
}

void AnaHighBW::finalize(BetaConfigMgr* const configMgr){
  // pass
}

// User functions
// ===============================================================================================================
void AnaHighBW::find_capacitors(const int dut_ch) {
    // calculate the window for average filter.
    double low_b = *(i_t[dut_ch_]->begin()) + winsize_;
    auto _end = std::lower_bound(i_t[dut_ch_]->begin(), i_t[dut_ch_]->end(), low_b);
    int winsize = std::distance(i_t[dut_ch_]->begin(), _end);
    
    // smoothing the voltage-time trace with fixed window moving averaging.
    *dut_process_voltage = wm::Filter::WindowMean(*i_w[dut_ch_], *i_t[dut_ch_], winsize_);
    
    // compute the divation as score for the given window size.
    *dut_process_voltage = wm::Filter::StdScore(*dut_process_voltage, winsize);

    // this is for baseline correction. we probably don't need it.
    // *dut_process_voltage = wm::Baseline::NoiseMedian(*dut_process_voltage, dut_process_voltage->size()/100);
    // dut_process_time = i_t[dut_ch_];

    // finding the maximum of the deviation score.
    // and use the move the maximum time to the left by the rise time of 50% to 90%.
    auto max_pts = wm::FindMultipleSignalMax(*dut_process_voltage, *i_t[dut_ch_], peak_finding_threshold_);
    for(auto &pt : max_pts) {
        auto fall_time =  wm::CalcRiseTime(*dut_process_voltage, *i_t[dut_ch_], pt.index, 0.5, 0.9);
        auto t_low_iter = std::lower_bound(i_t[dut_ch_]->begin(), i_t[dut_ch_]->end(), i_t[dut_ch_]->at(pt.index) - fall_time);
        auto fall_time_loc = std::distance(i_t[dut_ch_]->begin(), t_low_iter);
        w_step->push_back(i_w[dut_ch_]->at(fall_time_loc));
        t_step->push_back(i_t[dut_ch_]->at(fall_time_loc));
        // w_step->push_back(i_w[dut_ch_]->at(pt.index));
        // t_step->push_back(i_t[dut_ch_]->at(pt.index));
    }
}

// ===============================================================================================================
// ===============================================================================================================
void AnaHighBW::find_capacitors(const int dut_ch, const int clock_ch) {

    // calculate the window for average filter.
    double low_b = *(i_t[clock_ch]->begin()) + winsize_;
    auto _end = std::lower_bound(i_t[clock_ch]->begin(), i_t[clock_ch]->end(), low_b);
    int winsize = std::distance(i_t[clock_ch]->begin(), _end);
    
    // smoothing the voltage-time trace with fixed window moving averaging.
    offset_correction(*i_w[clock_ch]);
    auto th_time = wm::FindTimeAtThreshold(*i_w[clock_ch], *i_t[clock_ch], 0.0);

    double win_tolerance = 0.45e-6;
    for(int i = 0; i < th_time.size() - 1; i++) {
        if(th_time[i+1]-th_time[i] < win_tolerance) continue;
        int _t_start = std::distance(i_t[dut_ch]->begin(), std::lower_bound(i_t[dut_ch]->begin(), i_t[dut_ch]->end(), th_time[i]));
        int _t_end = std::distance(i_t[dut_ch]->begin(), std::lower_bound(i_t[dut_ch]->begin(), i_t[dut_ch]->end(), th_time[i+1]));
        auto _iter_start = i_w[dut_ch]->begin();
        auto _iter_end = i_w[dut_ch]->begin();
        std::advance(_iter_start, _t_start);
        std::advance(_iter_end, _t_end);
        auto mean = std::accumulate(_iter_start, _iter_end, 0.0) / std::distance(_iter_start, _iter_end);
        w_step->push_back(mean);
        int mid_time = std::distance(i_t[dut_ch]->begin(), std::lower_bound(i_t[dut_ch]->begin(), i_t[dut_ch]->end(), (th_time[i+1]+th_time[i])*0.5));
        t_step->push_back(i_t[dut_ch_]->at(mid_time));
    }
}

// ===============================================================================================================
// ===============================================================================================================
void AnaHighBW::find_capacitors(
    const int dut_ch, 
    const int t_start,
    const int t_win,
    const int cycle_size,
    const int ncycles)
{   
    // int cycle_size = 1300;

    int trace_size = i_w[dut_ch]->size();

    for(int i = 0; i < ncycles; i++) { 
        int _start_pt = t_start + i * cycle_size;
        int _end_pt = t_start + i * cycle_size + t_win;
        int _mid_pt = t_start + i * cycle_size + t_win * 0.5;
        
        if(_end_pt >= trace_size) break;

        auto _start = i_w[dut_ch]->begin();
        auto _end = i_w[dut_ch]->begin();
        std::advance(_start, _start_pt);
        std::advance(_end, _end_pt);
        auto sum = std::accumulate(_start, _end, 0.0);
        int dN = std::distance(_start, _end);

        w_step->push_back(sum/dN);
        t_step->push_back(_mid_pt);
    }
}


// ===============================================================================================================
// ===============================================================================================================
void AnaHighBW::offset_correction(std::vector<double> &v_trace) {
    double mean = std::accumulate(v_trace.begin(), v_trace.end(), 0.0) / v_trace.size();
    auto _iter = v_trace.begin();
    while(++_iter != v_trace.end()) {
        *_iter = std::move(*_iter) - mean;
    }
}