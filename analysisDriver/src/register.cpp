#include "analysisDriver/register.hpp"

// include analysis header
#include "AnaSSRL/AnaSSRL.hpp"
#include "AnaTCT/AnaTCT.hpp"
#include "AnaG4/AnaG4.hpp"
#include "AnaSimple/AnaSimple.hpp"
#include "AnaQuickScan/AnaQuickScan.hpp"
#include "AnaHighBW/AnaHighBW.hpp"


bool AnalysisRegister::Run(){
    auto _ssrl_reg = AnalysisRegister::Register<AnaSSRL>("SSRL");
    auto _tct_reg = AnalysisRegister::Register<AnaTCT>("TCT");
    auto _g4_reg = AnalysisRegister::Register<AnaG4>("G4");
    auto _simple_reg = AnalysisRegister::Register<AnaSimple>("Simple");
    auto _quick_scan_reg = AnalysisRegister::Register<AnaQuickScan>("QuickScan");
    auto _high_bw_reg = AnalysisRegister::Register<AnaHighBW>("HighBW");

    return true;
}
