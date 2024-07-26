#include "analysisDriver/driver.hpp"
#include "analysisDriver/register.hpp"
#include "utilities/getFiles.hpp"
#include "utilities/logger.hpp"

#include<unistd.h>
#include<sys/types.h>
#include<sys/wait.h>
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
  ("directory,d", bpo::value<std::string>()->default_value("."), "directory for the input files")
  ("selector,s", bpo::value<std::string>()->required(), "analysis selector")
  ("config,c", bpo::value<std::string>()->default_value(""), "configuration file")
  ("mp", bpo::bool_switch()->default_value(false), "internal mp")
  ("nfile,n", bpo::value<int>()->default_value(5), "number of files for mp")
  ("input,i", bpo::value<std::string>()->default_value(""), "input file.")
  // ("skipWaveform", bpo::bool_switch()->default_value(false), "skipping waveform in output file.")
  // ("skim", bpo::bool_switch()->default_value(false), "skim the output file.")
  // ("mp", bpo::bool_switch()->default_value(false), "internal mp")
  // ("thread,j", bpo::value<unsigned>()->default_value(std::thread::hardware_concurrency()), "using number of thread.")
  ;
  bpo::store(bpo::parse_command_line(argc, argv, desc, style), vm);
  
  if(vm.count("help")){
    std::cout << desc << std::endl;
    return 0;
  }

  bpo::notify(vm);

  // verifying analysis
  std::string selector = vm["selector"].as<std::string>();
  if(AnalysisRegister::Run() && !AnalysisFactory::CheckAnalysis(selector)) {
    LOG_INFO("No analysis " + selector);
    return 1;
  }

  std::vector<std::string> files;
  if(!vm["input"].as<std::string>().empty()){
    files.push_back(vm["input"].as<std::string>());
  }
  else{
    LOG_INFO("Getting list of input files ");
    files = GetListOfFiles(vm["directory"].as<std::string>(), "*.root");
  }

  if(files.size() == 0){
    LOG_ERROR("Cannot find any files.")
    return 1;
  }

  LOG_INFO("Here is list of input files ");
  for(auto &fname : files){
    LOG_INFO(fname);
  }

  // do analysis uising fork
  if(vm["mp"].as<bool>()){
    LOG_INFO("Using internal MP.")
    int fcount = 0;
    int maxf = vm["nfile"].as<int>();
    int status = 0;
    int total_jobs = 0;
    for(auto &fname : files){
      fcount++;
      total_jobs++;
      pid_t pid = fork();
      if(pid==0){ // only child run the job
        LOG_INFO("Staring analysis with input file: " + fname);
        RunAnalysis(selector, fname, vm["config"].as<std::string>());
        exit(0);
      }
      for(; fcount >= maxf; fcount--){
        wait(&status);
        total_jobs--;
      }
    }
    for(; total_jobs!=0; total_jobs--){ // last waiting call
      wait(&status);
    }
    return 0;
  }

  // do analysis
  for(auto &fname : files){
    LOG_INFO("Staring analysis with input file: " + fname);
    RunAnalysis(selector, fname, vm["config"].as<std::string>());
  }

  return 0;
}
