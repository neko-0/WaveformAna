#include "waveformMethods/waveformMethods.hpp"

double waveform_methods::CalcCFDTime(
  const double &frac,
  const int &max_index,
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace)
{
  double frac_v_value = v_trace.at(max_index) * frac;
  double frac_t_value = 0.0;
  double cfd_time = 0.0;
  int frac_i = 0;

  bool success = false;

  for(int i = max_index; i > -1; i--){
    double v_value = v_trace.at(i);
    if(v_value <= frac_v_value){
      frac_i = i;
      frac_t_value = t_trace.at(i);
      success = true;
      break;
    }
  }

  if(frac_i == v_trace.size() - 1 ) frac_i--;

  if(success){
    frac_t_value = waveform_methods::LinearInterpolationX(
      frac_t_value, v_trace.at(frac_i),
      t_trace.at(frac_i+1), v_trace.at(frac_i),
      frac_v_value);
  }

  return frac_t_value;
}
