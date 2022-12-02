#include "analysisDriver/driver.hpp"
#include "analysisDriver/register.hpp"
#include "utilities/getFiles.hpp"
#include "utilities/logger.hpp"


#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

void RunAnalysis(
    const std::string &ana,
    const std::string &fname,
    const std::string &ext_config)
{
  auto ana_driver = AnalysisDriver();
  ana_driver.AnalysisSelector(ana);
  ana_driver.AddExternalConfig(ext_config);
  ana_driver.Initialize(fname);
  ana_driver.EventLoop();
  ana_driver.Finalize();
}

//==============================================================================
int main(int argc, char **argv){
  bpo::options_description desc("Analysis options.");
  bpo::variables_map vm;
  bpo::command_line_style::style_t style = bpo::command_line_style::style_t(
    bpo::command_line_style::unix_style|
    bpo::command_line_style::allow_long_disguise
  );
  desc.add_options()
  ("help,h", "help message.")
  ("directory,d", bpo::value<std::string>()->required(), "directory for the input files")
  ("selector,s", bpo::value<std::string>()->required(), "analysis selector")
  ("config", bpo::value<std::string>()->default_value(""), "configuration file")
  // ("skipWaveform", bpo::bool_switch()->default_value(false), "skipping waveform in output file.")
  // ("skim", bpo::bool_switch()->default_value(false), "skim the output file.")
  // ("mp", bpo::bool_switch()->default_value(false), "internal mp")
  // ("thread,j", bpo::value<unsigned>()->default_value(std::thread::hardware_concurrency()), "using number of thread.")
  ;
  bpo::store(bpo::parse_command_line(argc, argv, desc, style), vm);
  bpo::notify(vm);

  if(vm.count("help")){
    std::cout << desc << std::endl;
    return 0;
  }

  // verifying analysis
  std::string selector = vm["selector"].as<std::string>();
  if(AnalysisRegister::Run() && !AnalysisFactory::CheckAnalysis(selector)) {
    LOG_INFO("No analysis " + selector);
    return 1;
  }
  LOG_INFO("Getting list of input files ");

  auto files = GetListOfFiles(vm["directory"].as<std::string>(), "*.root");

  if(files.size() == 0){
    LOG_ERROR("Cannot find any files.")
    return 1;
  }

  LOG_INFO("Here is list of input files ");
  for(auto &fname : files){
    LOG_INFO(fname);
  }

  // do analysis
  for(auto &fname : files){
    LOG_INFO("Staring analysis with input file: " + fname);
    RunAnalysis(selector, fname, vm["config"].as<std::string>());
  }

  return 0;
}
