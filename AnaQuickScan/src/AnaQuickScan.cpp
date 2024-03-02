#include "AnaQuickScan/AnaQuickScan.hpp"
#include "waveformMethods/core.hpp"

#include <vector>
#include <string>
#include <iostream>

#include "yaml-cpp/yaml.h"

namespace wm = waveform_methods;

void AnaQuickScan::setup(BetaConfigMgr* const configMgr) { }

void AnaQuickScan::initialize(BetaConfigMgr* const configMgr){  

  int x_bins = 1000;
  double x_min = -1e4;
  double x_max = 1e4;

  int y_bins = 1000;
  double y_min = -1e4;
  double y_max = 1e4;

  double z_min = -1e4;
  double z_max = 1e4;
  
  if( !configMgr->ext_config_name().empty() ){
    LOG_INFO("external config: " + configMgr->ext_config_name());
    YAML::Node yaml_config = YAML::LoadFile(configMgr->ext_config_name());
    const auto ranges = yaml_config["ranges"];
    x_bins = ranges["x_bins"].as<int>();
    x_min = ranges["x_min"].as<double>();
    x_max = ranges["x_max"].as<double>();
    y_bins = ranges["y_bins"].as<int>();
    y_min = ranges["y_min"].as<double>();
    y_max = ranges["y_max"].as<double>();
    z_min = ranges["z_min"].as<double>();
    z_max = ranges["z_max"].as<double>();
  }

  pos = configMgr->SetInputBranch<std::vector<double>>("pos");
  t = configMgr->SetInputBranch<std::vector<double>>("time");

  for(int i = 0; i < this->num_ch_; i++){
    // input branches
    std::string current_ch = std::to_string(ch_start_ + i);
    w[i] = configMgr->SetInputBranch<std::vector<double>>("w" + current_ch);
    std::string pname = "prof" + current_ch;
    heatmap[i] = new TProfile2D(pname.c_str(), pname.c_str(), x_bins, x_min, x_max, y_bins, y_min, y_max, z_min, z_max);
  }

}

bool AnaQuickScan::execute(BetaConfigMgr* const configMgr){

  const double v_scale = 1.0; // V to mV
  const double t_scale = 1.0; // s to ps

  for(int ch = 0; ch < this->num_ch_; ch++){
    int trace_size = w[ch]->size();
    double baseline = wm::Baseline::CalcBaseline(*w[ch], *t, 240, 340);
    // inverting signal
    for(int i=0; i < trace_size; i++){
      w[ch]->at(i) -= baseline;
      w[ch]->at(i) *= -1.0*v_scale;
      // t[ch]->at(i) *= t_scale;
    }

    auto wave_pt = wm::FindSignalMax(*w[ch], *t, 600, 800);
    // std::cout << "max : " << wave_pt.v << ", " << baseline << "\n";
    auto f = heatmap[ch]->Fill(pos->at(0), pos->at(1), wave_pt.v, 1);
    // auto f = heatmap[ch]->Fill(pos->at(0), pos->at(1), 1e3, 1);
    // std::cout << ch << "ch, " << "max : " << wave_pt.v << ", " << pos->at(0) << ", " << pos->at(1) << ", " << f << "\n";
  }

  return false;
}

void AnaQuickScan::finalize(BetaConfigMgr* const configMgr){

  configMgr->GetOutputFile()->cd();
  for(int ch = 0; ch < this->num_ch_; ch++) {
    heatmap[ch]->Write();
    delete heatmap[ch];
  }
}
