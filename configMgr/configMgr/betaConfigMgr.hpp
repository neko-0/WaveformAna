#ifndef BETACONFIGMGR_H
#define BETACONFIGMGR_H

#include "configMgr/baseConfigMgr.hpp"
#include "utilities/logger.hpp"

#include <string>
#include <vector>
#include <unordered_map>

#include <TFile.h>
#include <TTree.h>

//==============================================================================

// helper function from stack overflow

template <class N> struct is_vector
{
  static const bool value = false;
  using dtype = N;
};


template <class N, class A> struct is_vector<std::vector<N, A>>
{
  static const bool value = true;
  using dtype = std::vector<N, A>;
};


//==============================================================================
struct BaseBranch
{
  BaseBranch(){};
  virtual ~BaseBranch(){};

  virtual void setup() = 0;
  virtual void clear() = 0;
  virtual void **branch_addr() = 0;
};

//==============================================================================

template <typename dtype>
struct BranchWrapper : public BaseBranch
{
private:
  std::string name_;
  using btype = dtype;
  bool delete_ = false;
  dtype *branch_ = nullptr;
public:
  BranchWrapper(){};
  BranchWrapper(const std::string &name):name_(name){};
  ~BranchWrapper(){
    if(delete_ && branch_){
      // LOG_INFO("clean up branch: " + name_)
      delete branch_;
    }
  }

  void setup(){ branch_ = new dtype; delete_ = true;}

  void **branch_addr(){ return (void**)&branch_; }
  btype **branch_dtype_addr(){ return &branch_; }

  dtype *branch(){ return branch_; }
  void branch(dtype *value){ branch_ = value; }
  void branch(dtype &value){ branch_ = &value; }

  void clear(){
    if constexpr(is_vector<dtype>::value){
      this->branch_->clear();
    }
  }

  void SetAutoDelete(bool value){delete_ = value;}
};

//==============================================================================
class BetaConfigMgr : public BaseConfigMgr
{
  std::string name_ = "BetaConfigMgr";
  std::string input_filename_ = "input_dummy_name";
  std::string input_treename_ = "wfm";
  std::string output_filename_ = "output_dummy_name";
  std::string output_treename_ = "wfm";
  std::string output_file_prefix_ = "stats_";

  bool is_open_ = false;

  int input_branch_counter_ = 0;
  int output_branch_counter_ = 0;

  int input_entries_;
  int current_entry_ = 0;

  std::vector<BaseBranch*> input_branches_ = {};
  std::vector<BaseBranch*> output_branches_ = {};
  std::map<std::string, int> input_branch_map_ = {};
  std::map<std::string, int> output_branch_map_ = {};
  std::vector<int> output_vector_branch_index_ = {};

protected:
  TFile *input_file = nullptr;
  TTree *input_tree = nullptr;

  TFile *output_file = nullptr;
  TTree *output_tree = nullptr;

public:

  BetaConfigMgr(){};
  BetaConfigMgr(const std::string &filename);

  ~BetaConfigMgr(){
    for(auto &branch : input_branches_){
      if(branch) delete branch;
    }
    for(auto &branch : output_branches_){
      if(branch) delete branch;
    }
    if(output_tree) output_tree->AutoSave();
    if(input_file) input_file->Close();
    if(output_file) output_file->Close();
  }

  virtual bool Initialize();
  virtual bool NextEvent();
  virtual bool Fill();
  virtual bool Finalize();
  void Clear();

  virtual std::string GetName(){return name_;};

  template <typename dtype>
  dtype *SetInputBranch(const std::string &name);

  template <typename dtype>
  dtype *SetOutputBranch(const std::string &name);

  template <typename dtype>
  dtype *GetInputBranchValue(const std::string &name);

  template <typename dtype>
  dtype *GetOutputBranchValue(const std::string &name);

  BaseBranch *GetInputBranch(const std::string &name);
  BaseBranch *GetOutputBranch(const std::string &name);

  void input_filename(const std::string &value){input_filename_ = value;}
  std::string input_filename(){return input_filename_;}

  void input_treename(const std::string &value){input_treename_ = value;}
  std::string input_treename(){return input_treename_;}

  void output_treename(const std::string &value){output_treename_ = value;}
  std::string output_treename(){return output_treename_;}

  void output_file_prefix(const std::string &value){output_file_prefix_ = value;}
  std::string output_file_prefix(){return output_file_prefix_;}

  TFile *GetInputFile(){return input_file;}
  TTree *GetInputTree(){return input_tree;}
  TFile *GetOutputFile(){return output_file;}
  TTree *GetOutputTree(){return output_tree;}

  int GetInputEntries(){return input_entries_;}
};

//==============================================================================
template <typename dtype>
dtype *BetaConfigMgr::SetInputBranch(const std::string &name)
{
  if(!this->input_tree) return nullptr;

  auto *my_branch = new BranchWrapper<dtype>(name);
  my_branch->setup();
  this->input_branches_.push_back(my_branch);

  auto branch_index = std::pair<std::string, int>(name, this->input_branch_counter_);
  this->input_branch_map_.insert(branch_index);

  // this->input_tree->SetBranchAddress(name.c_str(), (dtype**)my_branch->get());

  if constexpr(is_vector<dtype>::value){
    LOG_INFO(name + ": dtype=<stl::vector>, buffer/clear will be handle internally");
    this->input_tree->SetBranchAddress(name.c_str(), my_branch->branch_dtype_addr());
  } else {
    LOG_INFO(name + ": dtype=" + std::string(typeid(dtype).name()));
    this->input_tree->SetBranchAddress(name.c_str(), my_branch->branch());
  }

  this->input_branch_counter_++;
  return my_branch->branch();
}

template <typename dtype>
dtype *BetaConfigMgr::SetOutputBranch(const std::string &name)
{
  if(!this->output_tree) return nullptr;

  auto *my_branch = new BranchWrapper<dtype>(name);
  my_branch->setup();
  this->output_branches_.push_back(my_branch);

  auto branch_index = std::pair<std::string, int>(name, this->output_branch_counter_);
  this->output_branch_map_.insert(branch_index);

  this->output_tree->Branch(name.c_str(), my_branch->branch());
  // this->output_tree->Branch(name.c_str(), my_branch->branch_dtype_addr());

  if constexpr(is_vector<dtype>::value){
    LOG_INFO(name + ": dtype=<stl::vector>, buffer/clear will be handle internally");
    this->output_vector_branch_index_.push_back(this->output_branch_counter_);
    // my_branch->branch()->reserve(2000);
  } else {
    LOG_INFO(name + ": dtype=" + std::string(typeid(dtype).name()));
  }

  this->output_branch_counter_++;
  return my_branch->branch();
}

template <typename dtype>
dtype *BetaConfigMgr::GetInputBranchValue(const std::string &name){
  auto my_branch = BetaConfigMgr::GetInputBranch(name);
  if(my_branch) return static_cast<BranchWrapper<dtype>*>(my_branch)->branch();
  return nullptr;
}

template <typename dtype>
dtype *BetaConfigMgr::GetOutputBranchValue(const std::string &name){
  auto my_branch = BetaConfigMgr::GetOutputBranch(name);
  if(my_branch) return static_cast<BranchWrapper<dtype>*>(my_branch)->branch();
  return nullptr;
}

#endif
