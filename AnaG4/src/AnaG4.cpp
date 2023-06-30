#include "AnaG4/AnaG4.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaG4::initialize(BetaConfigMgr* const configMgr){

  LOG_INFO("external config: " + configMgr->ext_config_name());
  YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());
}

void AnaG4::execute(BetaConfigMgr* const configMgr){

}

void AnaG4::finalize(BetaConfigMgr* const configMgr){
  // pass
}
