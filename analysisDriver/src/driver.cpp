#include "analysisDriver/driver.hpp"
#include "utilities/logger.hpp"

#include <string>
#include <stdexcept>

//==============================================================================
void AnalysisDriver::AnalysisSelector(const std::string &name){
  this->user_ana = AnalysisFactory::SelectAnalysis(name);
  if(this->user_ana){
    LOG_INFO("Selected user analysis: " + name);
  }
  else{
    auto ana_map = AnalysisFactory::GetMap();
    LOG_ERROR("Number of registered analysis: " + std::to_string(ana_map->size()));
    throw std::runtime_error("Unable to select analysis: " + name);
  }
}

//==============================================================================
void AnalysisDriver::Initialize(const std::string &fname){
  LOG_INFO("Start initialization with file: " + fname);
  this->configMgr->input_filename(fname);
  this->configMgr->Initialize();
  this->user_ana->initialize(this->configMgr);
  LOG_INFO("Initialization finisehd.");
}

//==============================================================================
void AnalysisDriver::EventLoop(){
  LOG_INFO("Starting EventLoop.");
  while(this->configMgr->NextEvent()){
    if( counter_ % 1000 == 0 ||
        (counter_ % 10 ==0 && counter_ <=100) ||
        (counter_ == this->configMgr->GetInputEntries()-1) ){
      LOG_INFO("Proccesed number of events:" + std::to_string(counter_));
    }
    DoAnalysis();
    this->configMgr->Fill();
    counter_++;
  }
  LOG_INFO("EventLoop finished.");
}

//==============================================================================
void AnalysisDriver::DoAnalysis(){
  user_ana->execute(configMgr);
}

//==============================================================================
void AnalysisDriver::Finalize(){
  this->user_ana->finalize(this->configMgr);
  this->configMgr->Finalize();
}
