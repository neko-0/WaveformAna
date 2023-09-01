#ifndef ANA_G4_H
#define ANA_G4_H

#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

#include <vector>

struct AnaG4 : BaseAna {

  AnaG4(){};
  ~AnaG4(){};

  virtual void setup(BetaConfigMgr* const configMgr);
  virtual void initialize(BetaConfigMgr* const configMgr);
  virtual void finalize(BetaConfigMgr* const configMgr);
  virtual bool execute(BetaConfigMgr* const configMgr);

private:
  // ===========================================================================
  // input variables
  std::vector<double> *pion_x;
  std::vector<double> *pion_y;
  std::vector<double> *pion_z;
  std::vector<double> *pion_edep;
  std::vector<int> *pion_region;
  double *pion_ke;

  std::vector<double> *muon_x;
  std::vector<double> *muon_y;
  std::vector<double> *muon_z;
  std::vector<double> *muon_edep;
  std::vector<int> *muon_region;
  double *muon_ke;

  // ===========================================================================
  // output variables
  double *output_pion_last_x;
  double *output_pion_last_y;
  double *output_pion_last_z;
  double *output_pion_last_edep;
  std::vector<double> *output_pion_z;
  std::vector<double> *output_pion_edep;
  // std::vector<int> *output_pion_region;

  double *output_muon_last_x;
  double *output_muon_last_y;
  double *output_muon_last_z;
  double *output_muon_last_edep;
  std::vector<double> *output_muon_z;
  std::vector<double> *output_muon_edep;
  // std::vector<int> *output_muon_region;

};

#endif
