#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>

namespace waveform_methods {

double CalcNoise(
  const TraceD &v_trace,
  const int &imin,
  const int &imax)
{
  // const auto &begin = v_trace.begin() + imin;
  // const auto &end = v_trace.begin() + imax;
  // double avg = 0, sq = 0;
  // for(auto i=begin; i!=end; i++){
  //   avg += *i;
  //   sq += (*i)*(*i);
  // }

  double avg = 0, sq = 0;
  #pragma omp simd reduction(+:avg,sq)
  for(int i=imin; i!=imax; i++){
    avg += v_trace.at(i);
    sq += v_trace.at(i)*v_trace.at(i);
  }

  double ndiff = (imax - imin);
  avg /= ndiff;
  sq /= ndiff;

  return pow(sq - avg*avg, 0.5);
}

//==============================================================================
double CalcNoise(
  const TraceD &v_trace,
  const int &npts)
{
  return CalcNoise(v_trace, 0, npts);
}

//==============================================================================
double CalcNoise(
  const TraceD &v_trace,
  const double &frac)
{
  int npts = v_trace.size()*frac;
  return CalcNoise(v_trace, 0, npts);
}

//==============================================================================
double CalcNoise(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int lower_i = 0;
  int upper_i = t_trace.size();
  if(lower != t_trace.end()) {
    lower_i = std::distance(t_trace.begin(), lower);
  }
  if(upper != t_trace.end()) {
    upper_i = std::distance(t_trace.begin(), upper);
  }
  return CalcNoise(v_trace, lower_i, upper_i);
}

}
