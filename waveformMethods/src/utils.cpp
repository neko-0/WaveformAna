#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

namespace waveform_methods::Utils {

double LinearInterpolationX(
  const double &x1, const double &y1,
  const double &x2, const double &y2,
  const double &y)
{
  return x1 + (y - y1) * (x2 - x1) / (y2 - y1);
}

template<class dtype>
dtype MedienFilter(
  const std::vector<dtype> &data,
  const int &window_size)
{
  std::vector<dtype> hi;
  return hi;
}

}
