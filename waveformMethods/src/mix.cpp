#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>
#include <stdio.h>
#include <iostream>

namespace waveform_methods {

/**
 * Calculate the noise, baseline, and both negatve & positive max index.
 *
 * @param v_trace Containers for waveform trace.
 * @param start Staring index for baseline and noise
 * @param end Ending index for baseline and noise
 * @return sum of `values`, or 0.0 if `values` is empty.
 */
MixParamsSet1 CalcMaxNoiseBase(
  const TraceD &v_trace,
  const int &start,
  const int &end)
{
  double avg = 0.0, sq = 0.0;
  double neg_max = 0.0, pos_max = 0.0;
  int neg_max_i = 0, pos_max_i = 0;

  #pragma omp simd reduction(+:avg,sq)
  for(std::size_t i = 0; i < v_trace.size(); i++){
    double value = v_trace.at(i);
    if(i >= start && i <= end){
      avg += value;
      sq += value*value;
    }
    if(neg_max < value){
      neg_max = value;
      neg_max_i = i;
    }
    if(pos_max > value){
      pos_max = value;
      pos_max_i = i;
    }
  }

  double ndiff = std::max(end - start + 1, 1); // make sure not dividing zero.
  avg /= ndiff;
  sq /= ndiff;
  double rms = pow(sq - avg*avg, 0.5);

  return MixParamsSet1{avg, rms, neg_max_i, pos_max_i, neg_max, pos_max};
}

//==============================================================================

MixParamsSet1 CalcMaxNoiseBase(
  const TraceD &v_trace,
  const double &frac)
{
  int npts = v_trace.size()*frac;
  if(npts == 0) npts = v_trace.size(); // fall back to entire trace.

  return CalcMaxNoiseBase(v_trace, 0, npts);
}

}
