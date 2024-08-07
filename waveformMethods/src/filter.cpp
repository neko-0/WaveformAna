#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/filter.hpp"
#include "waveformMethods/utils.hpp"

#include <math.h>
#include <numeric>
#include <iostream>
#include <stdio.h>
// #include <unistd.h>

namespace waveform_methods::Filter {

//==============================================================================
std::vector<double> WindowMean(
  std::vector<double> &v_trace, 
  std::vector<double> &t_trace,
  const double winsize)
{
  std::vector<double> output(t_trace.size());

  _WindowMean(
    v_trace.begin(),
    v_trace.end(),
    t_trace.begin(),
    t_trace.end(),
    output.begin(),
    winsize);

  return output;
}

//==============================================================================
std::vector<double> StdScore(
  std::vector<double> &v_trace,
  const int winsize)
{
  int tot_pts = v_trace.size();

  std::vector<double> output(tot_pts);
  auto o_iter = output.begin();

  for(int i = 0; i < tot_pts; i++) {
    auto _start = v_trace.begin() + i;
    auto _end = i + winsize < tot_pts ? _start + winsize : v_trace.end();
    int dN = std::distance(_start, _end) - 1;
    if (dN <= 0) dN = 1;
    *o_iter = std::abs(*_start - std::move(std::accumulate(_start, _end, 0.0) / dN));
    o_iter++;
  }

  return output;
}

}
