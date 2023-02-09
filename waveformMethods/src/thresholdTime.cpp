#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

namespace waveform_methods {

std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start,
  const int &end,
  const double &threshold)
{
  std::vector<double> threshold_time = {};
  for(int i = start; i < end; i++){
    double v_value = v_trace.at(i);
    if(v_value >= threshold && v_trace.at(i-1) < threshold){
      t_value = Utils::LinearInterpolationX(
        t_trace.at(i-1), v_trace.at(i-1),
        t_trace.at(i), v_trace.at(i), threshold);
      threhosld_time.push_back(t_value);
    }
  }

  return threshold_time;
}

//==============================================================================
std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &threshold)
{
  return FindTimeAtThreshold(v_trace, t_trace, 0, v_trace.size(), threshold);
}

//==============================================================================
std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax,
  const double &threshold)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int start = std::distance(t_trace.begin(), lower);
  int end = std::distance(t_trace.begin(), upper);

  return FindTimeAtThreshold(v_trace, t_trace, start, end, threshold);
}

}
