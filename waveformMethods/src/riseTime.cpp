#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

namespace waveform_methods {

double CalcRiseTime(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start_index,
  const double &low,
  const double &high,
  const EdgeType &type)
{
  int trace_size = v_trace.size();

  int max_index = start_index;

  if(max_index < 0) {
      auto wave_max = FindSignalMax(v_trace, t_trace);
  }

  if(max_index > trace_size - 2 ) {
    max_index = trace_size - 2;
  }

  if(max_index < 0) return -1.0;

  int top_i = trace_size - 2, bot_i = 0;
  bool found_top = false, found_bot = false;
  double low_b = 0.0, high_b = 0.0;

  // find the lower and upper bound for rise time. (10% to 90% default)
  double vmax = v_trace.at(max_index);
  if(type == EdgeType::Fall){
    high_b = vmax * low;
    low_b = vmax * high;
  } else {
    low_b = vmax * low;
    high_b = vmax * high;
  }

  int i = max_index;
  while(true) {
    if(!found_top && v_trace.at(i) <= high_b){
      top_i = i;
      found_top = true;
    }
    if(!found_bot && v_trace.at(i) <= low_b){
      bot_i = i;
      found_bot = true;
    }
    if(found_top && found_bot){
      // make sure it's not he last point
      if(top_i >= trace_size-2) top_i--;
      if(bot_i >= trace_size-2) bot_i--;
      break;
    }
    if(type == EdgeType::Fall) {
      i++;
      if(i >= trace_size) break;
    } else {
      i--;
      if(i < 0) break;
    }
  }

  if(bot_i > trace_size - 2 ) bot_i = trace_size - 2;
  if(top_i > trace_size - 2 ) top_i = trace_size - 2;
  
  double t1 = Utils::LinearInterpolationX(
      t_trace.at(bot_i), v_trace.at(bot_i),
      t_trace.at(bot_i+1), v_trace.at(bot_i+1), low_b);

  double t2 = Utils::LinearInterpolationX(
      t_trace.at(top_i), v_trace.at(top_i),
      t_trace.at(top_i+1), v_trace.at(top_i+1), high_b);

  return t2 - t1;
}

//==============================================================================
double CalcFallTime(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start_index,
  const double &low,
  const double &high)
{
  return CalcRiseTime(v_trace, t_trace, start_index, low, high, EdgeType::Fall);
}

}
