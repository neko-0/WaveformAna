#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>

using namespace waveform_methods;

double waveform_methods::CalcNoise(
  const TraceD &v_trace,
  const int &imin,
  const int &imax)
{
  const auto &begin = v_trace.begin() + imin;
  const auto &end = v_trace.begin() + imax;

  double avg = 0, sq = 0;
  for(auto i=begin; i!=end; i++){
    avg += *i;
    sq += (*i)*(*i);
  }

  avg /= (imax - imin);
  sq /= (imax - imin);
  return pow(sq - avg*avg, 0.5);
}

//==============================================================================
double waveform_methods::CalcNoise(
  const TraceD &v_trace,
  const int &npts)
{
  return waveform_methods::CalcNoise(v_trace, 0, npts);
}

//==============================================================================
double waveform_methods::CalcNoise(
  const TraceD &v_trace,
  const double &frac)
{
  int npts = v_trace.size()*frac;
  return waveform_methods::CalcNoise(v_trace, 0, npts);
}

//==============================================================================
double waveform_methods::CalcNoise(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int lower_i = std::distance(t_trace.begin(), lower);
  int upper_i = std::distance(t_trace.begin(), upper);
  return waveform_methods::CalcNoise(v_trace, lower_i, upper_i);
}
