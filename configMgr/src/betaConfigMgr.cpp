#include "configMgr/betaConfigMgr.hpp"
#include "utilities/logger.hpp"

BetaConfigMgr::BetaConfigMgr(const std::string &filename)
  : input_filename_(filename)
{
  // pass
}

bool BetaConfigMgr::Initialize()
{
  LOG_INFO("Initializing with input file: " + this->input_filename_ );
  this->input_file = TFile::Open(this->input_filename_.c_str());
  this->is_open_ = !this->input_file->IsZombie();
  if(!this->is_open_)
  {
    LOG_ERROR(input_filename_ + ": it's Zombie file");
    return false;
  }
  this->input_tree = (TTree *)this->input_file->Get(input_treename_.c_str());
  this->input_entries_ = input_tree->GetEntries();

  LOG_INFO("Input file looks fine. continue...");

  std::string delimiter = "/";
  std::string ofile_name = this->input_filename_;
  while (int(ofile_name.find(delimiter)) != -1)
  {
    ofile_name.erase(0, ofile_name.find(delimiter) + delimiter.length());
  }
  this->output_filename_ = this->output_file_prefix_ += ofile_name;
  this->output_file = TFile::Open(this->output_filename_.c_str(), "RECREATE");
  this->output_tree = new TTree(this->output_treename_.c_str(), "");
  this->output_tree->SetDirectory(this->output_file);
  LOG_INFO("Created output file: " + output_filename_);
  LOG_INFO("Created output tree: " + output_treename_);

  return true;
}

void BetaConfigMgr::Clear()
{
  for(auto &i : this->output_vector_branch_index_){
    this->output_branches_.at(i)->clear();
  }
}

bool BetaConfigMgr::Finalize()
{
    return true;
}

bool BetaConfigMgr::NextEvent(){
  BetaConfigMgr::Clear();
  auto evt_check =  this->input_tree->GetEntry(current_entry_);
  current_entry_++;
  return evt_check;
}

bool BetaConfigMgr::Fill(){
  return this->output_tree->Fill();
}

BaseBranch *BetaConfigMgr::GetInputBranch(const std::string &name){
  auto index = this->input_branch_map_.find(name);
  if(index == this->input_branch_map_.end() ){
    return nullptr;
  }
  return this->input_branches_.at(index->second);
}

BaseBranch *BetaConfigMgr::GetOutputBranch(const std::string &name){
  auto index = this->output_branch_map_.find(name);
  if(index == this->output_branch_map_.end() ){
    return nullptr;
  }
  return this->output_branches_.at(index->second);
}


ConfigMgrRegister<BetaConfigMgr> reg("BetaConfigMgr");
