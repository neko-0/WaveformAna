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
    winsize_ = general["winsize"].as<double>();
    peak_finding_threshold_ = general["peak_finding_threshold"].as<double>();
    
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
    dut_process_time = configMgr->SetOutputBranch<std::vector<double>>("dut_p_t");
    dut_process_voltage = configMgr->SetOutputBranch<std::vector<double>>("dut_p_w");

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
