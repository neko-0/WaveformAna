#include "analysisDriver/driver.hpp"
#include "utilities/getFiles.hpp"
#include "utilities/logger.hpp"


int main(int argc, char **argv){
  auto selector = std::string(argv[1]);
  auto ana_map = AnalysisFactory::GetMap();
  if(AnalysisFactory::CheckAnalysis(selector)) {
    LOG_INFO("Found analysis " + selector);
  } else {
    LOG_ERROR("Cannot find analysis " + selector);
  }

  return 0;
}
