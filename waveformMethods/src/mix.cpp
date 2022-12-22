#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>

using namespace waveform_methods;

/**
 * Calculate the noise, baseline, and both negatve & positive max index.
 *
 * @param v_trace Containers for waveform trace.
 * @param start Staring index for baseline and noise
 * @param end Ending index for baseline and noise
 * @return sum of `values`, or 0.0 if `values` is empty.
 */
MixParamsSet1 waveform_methods::CalcMaxNoiseBase(
  const TraceD &v_trace,
  const int &start,
  const int &end)
{
  double avg = 0, sq = 0;
  double neg_max = 0, pos_max = 0;
  int neg_max_i = 0, pos_max_i = 0;

  #pragma omp simd reduction(+:avg,sq)
  for(int i = 0; i < v_trace.size(); i++){
    auto value = v_trace.at(i);
    if(i >= start && i <= end){
      avg += value;
      sq += value*value;
    }
    if(neg_max > value){
      neg_max = value;
      neg_max_i = i;
    }
    if(pos_max < value){
      pos_max = value;
      pos_max_i = i;
    }
  }

  double ndiff = end - start;
  avg /= ndiff;
  sq /= ndiff;
  double rms = pow(sq - avg*avg, 0.5);

  return MixParamsSet1{avg, rms, neg_max_i, pos_max_i, neg_max, pos_max};
}

//==============================================================================

MixParamsSet1 waveform_methods::CalcMaxNoiseBase(
  const TraceD &v_trace,
  const double &frac)
{
  int npts = v_trace.size()*frac;
  return waveform_methods::CalcMaxNoiseBase(v_trace, 0, npts);
}
