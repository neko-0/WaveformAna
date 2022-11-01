#ifndef ANALYSIS_DRIVER_H
#define ANALYSIS_DRIVER_H

#include "configMgr/betaConfigMgr.hpp"
#include "analysisDriver/baseAna.hpp"

#include <vector>
#include <string>

class AnalysisDriver {
  int counter_;
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
};


#endif
