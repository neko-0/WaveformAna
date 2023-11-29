#include "AnaSimple/AnaSimple.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaSimple::setup(BetaConfigMgr* const configMgr) { }

void AnaSimple::initialize(BetaConfigMgr* const configMgr) {

    LOG_INFO("external config: " + configMgr->ext_config_name());
    // YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

    YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

    const auto &general = yaml_config["general"];
    this->store_waveform_ = general["store_waveform"].as<bool>();
    this->invert_channels_ = general["invert_ch"].as<std::vector<int>>();
    this->use_single_input_time_ = general["use_single_input_t_trace"].as<bool>();
    this->ch_start_ = general["start_ch_index"].as<int>();
    this->search_tmax_ = general["search_tmax"].as<double>();
    this->search_tmin_ = general["search_tmin"].as<double>();
    

    for(int i = 0; i < num_ch_; i++) {
        // input branches
        std::string current_ch = std::to_string(ch_start_ + i);
        w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
        if(!w[i]) continue;
        if(!use_single_input_time_) {
            t[i] = configMgr->SetInputBranch<std::vector<double>>("t" + current_ch);
        } else{
            if(i==0) t[i] = configMgr->SetInputBranch<std::vector<double>>("time");
            else t[i] = t[0];
        }

        active_ch_.push_back(i);

        output_pmax[i] = configMgr->SetOutputBranch<double>("pmax"+current_ch);
        output_tmax[i] = configMgr->SetOutputBranch<double>("tmax"+current_ch);
        output_rms[i] = configMgr->SetOutputBranch<double>("rms"+current_ch);
        output_area[i] = configMgr->SetOutputBranch<double>("area"+current_ch);

        if(this->store_waveform_) {
            output_w_corr[i] = configMgr->SetOutputBranch<std::vector<float>>("w_corr" + current_ch);
            if(use_single_input_time_){
                if(found_single_t_trace_){
                    output_time[i] = nullptr;
                    continue;
                }
                output_time[i] = configMgr->SetOutputBranch<std::vector<float>>("t");
                found_single_t_trace_ = true;
                continue;
            }
            output_time[i] = configMgr->SetOutputBranch<std::vector<float>>("t" + current_ch);
        }
    }
}

bool AnaSimple::execute(BetaConfigMgr* const configMgr){
    for(auto &ch : active_ch_){
        if(w[ch]->size() == 0){
            LOG_WARNING("Trace size 0");
            continue;
        }
    
        auto baseline_params = wm::CalcMaxNoiseBase(*w[ch], 0.25);

        std::vector<double> corr_w;
        corr_w.reserve(w[ch]->size());
        double polarity = 1.0;
        if(std::find(this->invert_channels_.begin(), this->invert_channels_.end(), ch+ch_start_) != this->invert_channels_.end() ) {
            polarity = -1.0;
        }
        for(int i = 0; i < w[ch]->size(); i++){
            corr_w.push_back((w[ch]->at(i) - baseline_params.baseline) * polarity);
        }

        *output_rms[ch] = baseline_params.rms;

        auto wave_pt = wm::FindSignalMax(corr_w, *t[ch], search_tmin_, search_tmax_);
        
        *output_pmax[ch] = wave_pt.v;
        *output_tmax[ch] = wave_pt.t;

        *output_area[ch] = wm::CalcPulseArea(corr_w, *t[ch], wave_pt.index);

        if(store_waveform_){
            std::move(corr_w.begin(), corr_w.end(), std::back_inserter(*output_w_corr[ch]));
            if(!output_time[ch]) continue;
            std::move(t[ch]->begin(), t[ch]->end(), std::back_inserter(*output_time[ch]));
        }
    
    }

    return true;
}

void AnaSimple::finalize(BetaConfigMgr* const configMgr){
  // pass
}
