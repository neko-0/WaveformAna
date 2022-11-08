#ifndef ANALYSIS_DRIVER_H
#define ANALYSIS_DRIVER_H

#include "configMgr/betaConfigMgr.hpp"
#include "analysisDriver/baseAna.hpp"

#include <vector>
#include <string>
#include <chrono>

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<double> Second;
typedef std::chrono::time_point<Time> TimePoint;

class AnalysisDriver {
  int total_entries_ = 0;
  int counter_ = 0;
  int previous_count_ = 0;
  TimePoint t0_ = Time::now();
  std::vector<std::string> file_list_;

protected:
  BaseAna *user_ana;
  BetaConfigMgr *configMgr = new BetaConfigMgr();

public:

  virtual ~AnalysisDriver(){delete configMgr;}

  virtual void AnalysisSelector(const std::string &name);
  virtual void Initialize(const std::string &fname);
  virtual void EventLoop();
  virtual void DoAnalysis();
  virtual void Finalize();
  virtual void ReportSatus();
};


#endif
