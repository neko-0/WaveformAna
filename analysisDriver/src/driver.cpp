#include "analysisDriver/driver.hpp"
#include "analysisDriver/register.hpp"
#include "utilities/logger.hpp"

#include <string>
#include <stdexcept>

#include <chrono>
#include <stdio.h>

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<double> Second;

// stackoverflow
// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int get_ram(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while(fgets(line, 128, file) != NULL){
        if(strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

//==============================================================================
void AnalysisDriver::AnalysisSelector(const std::string &name){
  this->user_ana = AnalysisFactory::SelectAnalysis(name);
  if(this->user_ana){
    LOG_INFO("Selected user analysis: " + name);
  }
  else{
    auto &ana_map = AnalysisFactory::GetMap();
    LOG_ERROR("Number of registered analysis: " + std::to_string(ana_map.size()));
    throw std::runtime_error("Unable to select analysis: " + name);
  }
}

//==============================================================================
void AnalysisDriver::ReportSatus(){
  Second dt = Time::now() - t0_;
  int dN = (counter_ - previous_count_);
  float evt_s =  dN / dt.count();
  int eta = dN ? (total_entries_ - counter_) / evt_s / 60 : -1;
  int rss = get_ram();
  LOG_INFO("Proccesed number of events: "
    + std::to_string(counter_) + "/" + total_entries_ + ", "
    + std::to_string(evt_s) + " evt/s, "
    + std::to_string(rss) + " KB, "
    + std::to_string(eta) + " ETA(min)");
  previous_count_ = counter_;
  t0_ = Time::now();
}

//==============================================================================
void AnalysisDriver::AddExternalConfig(const std::string &name){
  this->configMgr->ext_config_name(name);
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
  t0_ = Time::now();
  total_entries_ = this->configMgr->GetInputEntries();
  while(this->configMgr->NextEvent()){
    if( counter_ % 1000 == 0 || (counter_ % 10 ==0 && counter_ <=100) ||
        (counter_ == this->configMgr->GetInputEntries()-1) ) {
          ReportSatus();
    }
    if(DoAnalysis()){
      this->configMgr->Fill();
    }
    counter_++;
  }
  LOG_INFO("EventLoop finished.");
}

//==============================================================================
bool AnalysisDriver::DoAnalysis(){
  return this->user_ana->execute(configMgr);
}

//==============================================================================
void AnalysisDriver::Finalize(){
  this->user_ana->finalize(this->configMgr);
  this->configMgr->Finalize();
}
