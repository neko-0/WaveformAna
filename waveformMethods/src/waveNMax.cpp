#include "waveformMethods/waveformMethods.hpp"

#include <math.h>

using namespace waveform_methods;

WavePoints waveform_methods::_FindIdenticalSignalMax(
  const TraceDType &v_trace,
  const TraceDType&t_trace,
  const int &start, const int &end, const double &baseline)
{
  WavePoints wave_list = {};
  std::vector<int> max_indexes = {};

  double vmax = baseline;
  for(int i = start; i < end; i++) {
    double current_value = v_trace.at(i);
    if(current_value < vmax) continue;
    else if(current_value == vmax){
      max_indexes.push_back(i);
    }
    else {
      vmax = current_value;
      max_indexes.clear();
      max_indexes.push_back(i);
    }
  }

  for(const auto &i : max_indexes) {
    auto wave_max = WavePoint{i, v_trace.at(i), t_trace.at(i)};
    wave_list.push_back(wave_max);
  }

  return wave_list;
}

//==============================================================================
WavePoints waveform_methods::FindIdenticalSignalMax(
  const TraceDType &v_trace,
  const TraceDType &t_trace,
  const int &start, const int &end)
{
  int trace_size = v_trace.size();
  int start_i = std::min(start, trace_size-1);
  int end_i = std::min(end, trace_size-2);
  return waveform_methods::_FindIdenticalSignalMax(
    v_trace, t_trace, start_i, end_i);
}

//==============================================================================
WavePoints waveform_methods::FindIdenticalSignalMax(
  const TraceDType &v_trace,
  const TraceDType &t_trace)
{
  return waveform_methods::_FindIdenticalSignalMax(
    v_trace, t_trace, 0, v_trace.size());
}

//==============================================================================
WavePoints waveform_methods::_FindMultipleSignalMax(
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  const int &start_i,
  const int &end_i,
  const double &threshold,
  const double &scale /* = 2.0 */)
{
    WavePoints max_pts = {};

    double pmax = 0.0;
    int pmax_index = 0;
    bool is_signal = false;

    for(int i = start_i; i < end_i; i++)
    {
      double v_value = v_trace.at(i);
      if(!is_signal) // trying to find the 1st pmax
      {
        if(v_value < threshold) continue; // need to be at least greater than the threshold
        if(v_value < pmax) continue;
        pmax = v_value; // update current pmax value
        pmax_index = i;
        is_signal = true;
      }
      else // continue to update pmax
      {
        if(v_value > pmax)
        {
            pmax = v_value;
            pmax_index = i;
        }
        else if( (v_value < threshold) &&
                 (threshold - abs(v_value)) <= (threshold / scale) )
        { // if the v_value go well below the threshold, store the previous wave point
            max_pts.push_back( WavePoint{pmax_index, pmax, t_trace.at(pmax_index)} );
            // reset the pmax value to current point.
            pmax = v_value;
            pmax_index = i;
            is_signal = false;
        }
        else if(i == end_i - 1) // sotre the last point if no others are found.
        {
            max_pts.push_back( WavePoint{pmax_index, pmax, t_trace.at(pmax_index)} );
        }
      }
    }

    // feed a default value if nothing is found.
    if(max_pts.size() == 0) max_pts.push_back( WavePoint{-1, 10e11, 10e11} );

    return max_pts;
}

//==============================================================================
WavePoints waveform_methods::FindMultipleSignalMax(
  const std::vector<double> &v_trace,
  const std::vector<double> &t_trace,
  const double &threshold,
  const double &scale /* = 2.0 */)
{
  return _FindMultipleSignalMax(
    v_trace, t_trace, 0, v_trace.size(), threshold, scale);
}
