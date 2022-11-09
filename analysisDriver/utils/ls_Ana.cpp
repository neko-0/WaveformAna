#include "analysisDriver/driver.hpp"
#include "utilities/getFiles.hpp"
#include "utilities/logger.hpp"


int main(int argc, char **argv){
  auto &ana_map = AnalysisFactory::GetMap();
  LOG_INFO("Number of registered analysis: " + std::to_string(ana_map.size()));
  for(auto &ana : ana_map){
    LOG_INFO("Registered analysis: " + ana.first);
  }

  return 0;
}
