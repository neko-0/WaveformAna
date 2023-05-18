#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/utils.hpp"

#include <cmath>
#include <iostream>
#include <stdio.h>

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

  std::vector<dtype> lbound(half_lw + window_size);
  std::vector<dtype> rbound(half_rw + window_size);

  // extending left boundaries by mirroring.
  for(std::size_t i = 0; i < half_lw + window_size; i++){
    std::size_t j;
    if(i < half_lw){
      j = half_lw - i;
    } else {
      j = i - half_lw;
    }
    lbound[i] = data[j];
  }
  for(std::size_t i = 0; i < half_lw; i++){
    auto start = lbound.begin() + i;
    auto end = start + window_size;
    std::vector<dtype> temp(start, end);
    std::nth_element(temp.begin(), temp.begin() + half_lw, temp.end());
    output[i] = temp[half_lw];
  }

  // central region
  for(std::size_t i=half_lw; i < size-half_rw; i++){
    auto start = data.begin()+i-half_lw;
    auto end = start + window_size;
    std::vector<dtype> temp(start, end);
    std::nth_element(temp.begin(), temp.begin()+half_lw, temp.end());
    output[i] = temp[half_lw];
  }

  // extending right boundary by mirroring
  for(std::size_t i = 0; i < half_rw + window_size; i++){
    std::size_t j;
    if(i < window_size){
      j = i + size - window_size - 1;
    } else {
      j = size - (i - window_size) - 1;
    }
    rbound[i] = data[j];
  }
  for(std::size_t i = 0; i < half_rw; i++){
    auto start = rbound.begin() + i;
    auto end = start + window_size;
    std::vector<dtype> temp(start, end);
    std::nth_element(temp.begin(), temp.begin() + half_lw, temp.end());
    output[i + size - half_rw] = temp[half_lw];
  }

  // for(auto x : output) std::cout << x << "\n";
  // for(auto x : lbound) std::cout << x << "\n";
  // getchar();

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
  for(int i = -half_win; i <= half_win; i++){
    double dx = i-mean;
    double value = height*exp(sig_sq*dx*dx);
    gaus.emplace_back(value);
    sum += value;
  }
  for(auto &x : gaus){x /= sum;}

  return gaus;
}

}
