#include "configMgr/betaConfigMgr.hpp"
#include "utilities/getFiles.hpp"
#include "waveformMethods/waveformMethods.hpp"

#include <iostream>

int main(int argc, char **argv){
  auto baseConfig = ConfigMgrFactory::CreateConfigMgr("BetaConfigMgr");

  auto my_betaconfig = *ConfigMgrFactory::InitUserConfigMgr<BetaConfigMgr>(baseConfig);

  std::cout << my_betaconfig.GetName() << "\n";

  my_betaconfig.input_filename("test.root");
  my_betaconfig.Initialize();

  auto *w2 = my_betaconfig.SetInputBranch<std::vector<double>>("w2");
  auto *t2 = my_betaconfig.SetInputBranch<std::vector<double>>("t2");

  auto files = GetListOfFiles(".", "*.root");
  for(auto &f : files){
    std::cout << f << "\n";
  }

  my_betaconfig.NextEvent();
  my_betaconfig.NextEvent();

  std::cout << "size "<< t2->size() << "\n";
  auto area = waveform_methods::CalcPulseArea(*w2, *t2);

  std::cout << "Area " << area << "\n";

  auto rise = waveform_methods::CalcRiseTime(*w2, *t2);
  auto fall = waveform_methods::CalcFallTime(*w2, *t2);

  std::cout << "Rise " << rise << " Fall " << fall << "\n";

  for(auto &i : *w2){
    i *= -1.0;
  }
  auto pmaxes = waveform_methods::FindIdenticalSignalMax(*w2, *t2);
  std::cout << "N Max " << pmaxes.size() << "\n";

  auto output_area = my_betaconfig.SetOutputBranch<double>("area");

  auto dy = my_betaconfig.GetOutputBranchValue<double>("double_type");

  *output_area = area;

  dy = &w2->at(10);

  if(w2){
    std::cout << "\nGOOD\n" << *dy;
  }
  my_betaconfig.SetOutputBranch<std::vector<double>>("double_v_type");

  my_betaconfig.Fill();
  my_betaconfig.Finalize();
}
