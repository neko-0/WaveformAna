#include "waveformMethods/waveformMethods.hpp"

using namespace waveform_methods;

//==============================================================================
WavePoint waveform_methods::FindSignalMax(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start,
  const int &end)
{
  double vmax = v_trace.at(start);
  double tmax = t_trace.at(start);
  int index = -1;

  int trace_size = t_trace.size();
  for(int i = start; i < (end <= trace_size ? end : trace_size); i++)
  {
    if(v_trace.at(i) < vmax) continue;
    vmax = v_trace.at(i);
    tmax = t_trace.at(i);
    index = i;
  }

  return waveform_methods::WavePoint{index, vmax, tmax};
}

//==============================================================================
WavePoint waveform_methods::FindSignalMax(
  const TraceD &v_trace,
  const TraceD &t_trace)
{
  return waveform_methods::FindSignalMax(v_trace, t_trace, 0, v_trace.size());
}

//==============================================================================
WavePoint waveform_methods::FindSignalMax(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &min,
  const double &max)
{

  // find the closest starting point
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), min);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), max);

  int start_i = 0;
  int end_i = t_trace.size();
  // fall back to entire time window if no lower bound found.
  if( lower != t_trace.end() ) {
    start_i =  std::distance(t_trace.begin(), lower);
  }
  if( upper != t_trace.end() ) {
    end_i = std::distance(t_trace.begin(), upper);
  }

  return waveform_methods::FindSignalMax(v_trace, t_trace, start_i, end_i);
}
