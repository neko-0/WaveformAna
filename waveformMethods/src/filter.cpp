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

  WindowMean(v_trace.begin(), v_trace.end(), t_trace.begin(), t_trace.end(), output.begin(), winsize);

  return output;
}

}
