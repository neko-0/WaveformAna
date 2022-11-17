#include "analysisDriver/register.hpp"

// include analysis header
#include "AnaSSRL/AnaSSRL.hpp"
#include "AnaTCT/AnaTCT.hpp"


bool AnalysisRegister::Run(){
  auto _ssrl_reg = AnalysisRegister::Register<AnaSSRL>("SSRL");
  auto _tct_reg = AnalysisRegister::Register<AnaTCT>("TCT");

  return true;
}
