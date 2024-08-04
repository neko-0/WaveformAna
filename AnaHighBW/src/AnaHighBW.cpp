#include "AnaHighBW/AnaHighBW.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaHighBW::setup(BetaConfigMgr* const configMgr) { }

void AnaHighBW::initialize(BetaConfigMgr* const configMgr) {

    LOG_INFO("external config: " + configMgr->ext_config_name());
    // YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

    // YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());
}

bool AnaHighBW::execute(BetaConfigMgr* const configMgr){
    // for(auto &ch : active_ch_){
    //     if(w[ch]->size() == 0){
    //         LOG_WARNING("Trace size 0");
    //         continue;
    //     }
    // }

    return true;
}

void AnaHighBW::finalize(BetaConfigMgr* const configMgr){
  // pass
}
