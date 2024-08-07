#ifndef WAVEFOMR_METHODS_FILTER_H
#define WAVEFOMR_METHODS_FILTER_H

#include "waveformMethods/waveformMethods.hpp"

#include <math.h>
#include <numeric>

namespace waveform_methods::Filter
{
  
  template<class InputIt, class OutputIt> constexpr
  OutputIt WindowMean(
    InputIt v_iter_first,
    InputIt v_iter_last, 
    InputIt t_iter_first,
    InputIt t_iter_last,
    OutputIt o_iter,
    const double winsize) 
    {
      // typedef typename std::iterator_traits<InputIt>::value_type value_t;
      // value_t mean = *v_iter_first;

      // lower bound.
      double low_b = *t_iter_first + winsize;

      // calculate number of points for the given window size.
      int npts = std::distance(t_iter_first, std::lower_bound(t_iter_first, t_iter_last, low_b));
      int tot_pts = std::distance(t_iter_first, t_iter_last);

      if(v_iter_first == v_iter_last) return o_iter;

      for(int i = 0; i < tot_pts; i++) {
        auto _start = v_iter_first + i;
        auto _end = i + npts < tot_pts ? _start + npts : v_iter_last;
        int dN = std::distance(_start, _end) - 1;
        if (dN <= 0) dN = 1;
        *o_iter = std::accumulate(_start, _end, 0.0) / dN;
        o_iter++;
      }

      return o_iter;
    }

  std::vector<double> WindowMean(
    std::vector<double> &v_trace, 
    std::vector<double> &t_trace,
    const double winsize);


};

#endif
