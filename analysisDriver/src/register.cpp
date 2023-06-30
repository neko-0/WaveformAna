#include "analysisDriver/register.hpp"

// include analysis header
#include "AnaSSRL/AnaSSRL.hpp"
#include "AnaTCT/AnaTCT.hpp"
#include "AnaG4/AnaG4.hpp"


bool AnalysisRegister::Run(){
  auto _ssrl_reg = AnalysisRegister::Register<AnaSSRL>("SSRL");
  auto _tct_reg = AnalysisRegister::Register<AnaTCT>("TCT");
  auto _g4_reg = AnalysisRegister::Register<AnaTCT>("G4");

  return true;
}
