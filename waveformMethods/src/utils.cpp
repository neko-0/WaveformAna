#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

#include <cmath>
#include <iostream>

namespace waveform_methods::Utils {

double LinearInterpolationX(
  const double &x1, const double &y1,
  const double &x2, const double &y2,
  const double &y)
{
  return x1 + (y - y1) * (x2 - x1) / (y2 - y1);
}

template<class dtype>
std::vector<dtype> MedianFilter(
  const std::vector<dtype> &data,
  const int &window_size)
{
  int size = data.size();
  int half_lw = window_size / 2;
  int half_rw = window_size - half_lw;
  std::vector<dtype> output(size);

  // storing the left boundary as unchanged from data
  for(int i=0; i<half_lw; i++){
    output[i] = data[i];
  }

  for(int i=half_lw; i < size-half_rw; i++){
    auto start = data.begin()+i-half_lw;
    auto end = start + window_size;
    std::vector<dtype> temp(start, end);
    std::nth_element(temp.begin(), temp.begin()+half_lw, temp.end());
    output[i] = temp[half_lw];
  }

  // right boundary
  for(int i=size-1; i>size-half_rw; i--){
    output[i] = data[i];
  }

  return output;
}
template std::vector<int> MedianFilter<int>(const std::vector<int>&, const int&);
template std::vector<float> MedianFilter<float>(const std::vector<float>&, const int&);
template std::vector<double> MedianFilter<double>(const std::vector<double>&, const int&);

std::vector<double> GaussianKernel(
  const int &window_size,
  const double &height,
  const double &mean,
  const double &sigma)
{
  int half_win = window_size / 2;
  std::vector<double> gaus;
  gaus.reserve(window_size);
  double sum = 0.0;
  double sig_sq = -0.5 / (sigma*sigma);
  for(int i = -half_win; i < half_win; i++){
    double value = height*exp(sig_sq*(i-mean)*(i-mean));
    gaus.emplace_back(value);
    sum += value;
}
  for(auto &x : gaus){x /= sum;}

  return gaus;
}

}
