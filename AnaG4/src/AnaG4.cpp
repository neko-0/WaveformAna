#include "AnaG4/AnaG4.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <omp.h>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaG4::setup(BetaConfigMgr* const configMgr){
  LOG_INFO("Starting setup for configMgr");
  configMgr->input_treename("lgad");
  configMgr->output_treename("lgad");
}

void AnaG4::initialize(BetaConfigMgr* const configMgr){

  LOG_INFO("external config: " + configMgr->ext_config_name());
  // YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());

  pion_x = configMgr->SetInputBranch<std::vector<double>>("pion_x");
  pion_y = configMgr->SetInputBranch<std::vector<double>>("pion_y");
  pion_z = configMgr->SetInputBranch<std::vector<double>>("pion_z");
  pion_edep = configMgr->SetInputBranch<std::vector<double>>("pion_edep");
  pion_ke = configMgr->SetInputBranch<double>("pion_ke");
  pion_region = configMgr->SetOutputBranch<std::vector<int>>("pion_region");


  muon_x = configMgr->SetInputBranch<std::vector<double>>("muon_x");
  muon_y = configMgr->SetInputBranch<std::vector<double>>("muon_y");
  muon_z = configMgr->SetInputBranch<std::vector<double>>("muon_z");
  muon_edep = configMgr->SetInputBranch<std::vector<double>>("muon_edep");
  muon_ke = configMgr->SetInputBranch<double>("muon_ke");
  muon_region = configMgr->SetInputBranch<std::vector<int>>("muon_region");


  output_pion_edep = configMgr->SetOutputBranch<std::vector<double>>("pion_edep");
  output_pion_z = configMgr->SetOutputBranch<std::vector<double>>("pion_z");
  output_pion_last_x = configMgr->SetOutputBranch<double>("pion_last_x");
  output_pion_last_y = configMgr->SetOutputBranch<double>("pion_last_y");
  output_pion_last_z = configMgr->SetOutputBranch<double>("pion_last_z");
  output_pion_last_edep = configMgr->SetOutputBranch<double>("pion_last_edep");
  // output_pion_region = configMgr->SetOutputBranch<std::vector<int>>("pion_region");

  output_muon_edep = configMgr->SetOutputBranch<std::vector<double>>("muon_edep");
  output_muon_z = configMgr->SetOutputBranch<std::vector<double>>("muon_z");
  output_muon_last_x = configMgr->SetOutputBranch<double>("muon_last_x");
  output_muon_last_y = configMgr->SetOutputBranch<double>("muon_last_y");
  output_muon_last_z = configMgr->SetOutputBranch<double>("muon_last_z");
  output_muon_last_edep = configMgr->SetOutputBranch<double>("muon_last_edep");
  // output_muon_region = configMgr->SetOutputBranch<std::vector<int>>("muon_region");

}

bool AnaG4::execute(BetaConfigMgr* const configMgr){
  // only considering the stop pion
  if(*pion_ke > 0.0) return false;

  *output_pion_last_x = pion_x->end()[-1];
  *output_pion_last_y = pion_y->end()[-1];
  *output_pion_last_z = pion_z->end()[-1];
  *output_pion_last_edep = pion_edep->end()[-2];

  if(*muon_ke <= 0.0){
    *output_muon_last_x = muon_x->end()[-1];
    *output_muon_last_y = muon_y->end()[-1];
    *output_muon_last_z = muon_z->end()[-1];
    *output_muon_last_edep = muon_edep->end()[-2];
  } else {
    *output_muon_last_x = -999.0;
    *output_muon_last_y = -999.0;
    *output_muon_last_z = -999.0;
    *output_muon_last_edep = -999.0;
  }

  for(std::size_t i = 0; i < pion_edep->size(); i++){
    output_pion_edep->push_back(pion_edep->at(i));
    output_pion_z->push_back(pion_z->at(i));
  }
  for(std::size_t i = 0; i < muon_edep->size(); i++){
    output_muon_edep->push_back(muon_edep->at(i));
    output_muon_z->push_back(muon_z->at(i));
  }

  // std::move(pion_edep->begin(), pion_edep->end(), std::back_inserter(*output_pion_edep));
  // std::move(muon_edep->begin(), muon_edep->end(), std::back_inserter(*output_muon_edep));

  // std::move(pion_z->begin(), pion_z->end(), std::back_inserter(*output_pion_z));
  // std::move(muon_z->begin(), muon_z->end(), std::back_inserter(*output_muon_z));

  return true;
}

void AnaG4::finalize(BetaConfigMgr* const configMgr){
  // pass
}
