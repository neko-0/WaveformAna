#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

namespace waveform_methods {

std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start,
  const int &end,
  const double &threshold,
  const bool &early_stop)
{
  std::vector<double> threshold_time = {};

  for(std::size_t i = start > 0 ? start : 1; i < end; i++){
    if(v_trace.at(i) < threshold) continue;
    if(v_trace.at(i-1) >= threshold) continue;
    double t_value = Utils::LinearInterpolationX(
      t_trace.at(i-1), v_trace.at(i-1),
      t_trace.at(i), v_trace.at(i), threshold);
    threshold_time.push_back(t_value);
    if( early_stop ) break;
  }

  if(threshold_time.size() == 0) threshold_time.push_back(std::numeric_limits<double>::min());
  
  return threshold_time;
}

//==============================================================================
std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &threshold,
  const bool &early_stop)
{
  return FindTimeAtThreshold(v_trace, t_trace, 1, v_trace.size(), threshold, early_stop);
}

//==============================================================================
std::vector<double> FindTimeAtThreshold(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax,
  const double &threshold,
  const bool &early_stop)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int start = std::distance(t_trace.begin(), lower);
  int end = std::distance(t_trace.begin(), upper);

  return FindTimeAtThreshold(v_trace, t_trace, start, end, threshold, early_stop);
}

}
