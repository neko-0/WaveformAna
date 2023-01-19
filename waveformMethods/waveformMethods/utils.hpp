#ifndef WAVEFOMR_METHODS_UTILS_H
#define WAVEFOMR_METHODS_UTILS_H

#include "waveformMethods/waveformMethods.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace waveform_methods::Utils {

  double LinearInterpolationX(
    const double &x1, const double &y1,
    const double &x2, const double &y2,
    const double &y);

  std::vector<double> GaussianKernel(
    const int &window_size,
    const double &height,
    const double &mean,
    const double &sigma);

  template<class dtype>
  std::vector<dtype> MedianFilter(
    const std::vector<dtype> &data,
    const int &window_size);


}

#endif
