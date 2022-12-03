#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>
#include <iostream>
//==============================================================================
double waveform_methods::CalcBaseline(
  const TraceD &v_trace,
  const int &start,
  const int &end)
{
  auto iter0 = v_trace.begin();
  auto iter1 = v_trace.begin();
  std::advance(iter0, start);
  std::advance(iter1, end);
  double sum = std::accumulate(iter0, iter1, 0.0);

  return sum / (end-start);
}

//==============================================================================
double waveform_methods::CalcBaseline(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int dN = std::distance(lower, upper);
  double sum = std::accumulate(lower, upper, 0);

  return sum / dN;
}

//==============================================================================
bool waveform_methods::MultiSignalBaselineCorrection(
  WavePoints &signals,
  TraceD &v_trace,
  const double &frac_npts,
  const int &npts_forward,
  const int &npts_backward)
{
  int nsignal = signals.size();
  if(nsignal < 2) return false;

  bool corrected = false;
  for(int i = 1; i < nsignal; i++){
    int sig1_i = signals.at(i-1).index;
    int sig2_i = signals.at(i).index;
    if((sig2_i - (sig1_i+npts_backward)) > npts_forward){
      int start = sig2_i-npts_forward;
      int end = sig2_i+npts_backward;
      int end_base = npts_forward * frac_npts;
      if(start < 0) continue;
      if(end >= v_trace.size()) continue;
      if(start+end_base >= v_trace.size()) continue;
      double baseline = waveform_methods::CalcBaseline(v_trace, start, start+end_base);
      for(int j = start; j<end; j++){
        v_trace.at(j) -= baseline;
      }
      signals.at(i).v -= baseline;
      corrected = true;
    }
  }

  return corrected;
}

//==============================================================================
bool waveform_methods::MultiSignalBaselineCorrection(
  WavePoints &signals,
  TraceD &v_trace,
  const TraceD &t_trace,
  const double &frac_npts,
  const double &forward_time,
  const double &backward_time)
{
  double dt = t_trace.at(1) - t_trace.at(0);
  int npt_f = forward_time / dt;
  int npt_b = backward_time / dt;
  return waveform_methods::MultiSignalBaselineCorrection(
    signals, v_trace, frac_npts, npt_f, npt_b);
}
