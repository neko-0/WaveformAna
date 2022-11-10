#include "analysisDriver/register.hpp"

// include analysis header
#include "AnaSSRL/AnaSSRL.hpp"


bool AnalysisRegister::Run(){
  auto _ssrl_reg = AnalysisRegister::Register<AnaSSRL>("SSRL");

  return true;
}
