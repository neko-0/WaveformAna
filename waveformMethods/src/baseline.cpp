#include "waveformMethods/waveformMethods.hpp"
#include "waveformMethods/baseline.hpp"
#include "waveformMethods/utils.hpp"

#include <math.h>
#include <numeric>
#include <iostream>
#include <stdio.h>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
// #include <Eigen/IterativeLinearSolvers>
#include <Eigen/Geometry>
// #include <unistd.h>

namespace waveform_methods::Baseline {

//==============================================================================
double CalcBaseline(
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
double CalcBaseline(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const double &tmin,
  const double &tmax)
{
  auto lower = std::lower_bound(t_trace.begin(), t_trace.end(), tmin);
  auto upper = std::lower_bound(t_trace.begin(), t_trace.end(), tmax);
  int dN = std::distance(lower, upper);
  double sum = std::accumulate(lower, upper, 0.0);

  return sum / dN;
}

//==============================================================================
bool MultiSignalBaselineCorrection(
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
      double baseline = CalcBaseline(v_trace, start, start+end_base);
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
bool MultiSignalBaselineCorrection(
  WavePoints &signals,
  TraceD * const v_trace,
  const double &frac_npts,
  const int &npts_forward,
  const int &npts_backward)
{
  return MultiSignalBaselineCorrection(
    signals, *v_trace, frac_npts, npts_forward, npts_backward);
}

//==============================================================================
bool MultiSignalBaselineCorrection(
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
  return MultiSignalBaselineCorrection(
    signals, v_trace, frac_npts, npt_f, npt_b);
}

//==============================================================================
bool MultiSignalBaselineCorrection(
  WavePoints &signals,
  TraceD * const v_trace,
  const TraceD * const t_trace,
  const double &frac_npts,
  const double &forward_time,
  const double &backward_time)
{
  return MultiSignalBaselineCorrection(
    signals, *v_trace, *t_trace, frac_npts, forward_time, backward_time);
}

//==============================================================================
/*
  ARPLS, see
  Baseline correction using asymmetrically reweighted penalized least squares smoothing.
    Baek, S.J., et al.
    Analyst, 2015, 140, 250-257.
*/
TraceD ARPLS(
  TraceD &v_trace,
  const double &lam,
  const int &max_iter,
  const double &tol,
  const bool &do_correction)
{
  int size = v_trace.size();
  // Eigen::setNbThreads(5);
  // Assuming second order difference matrix diff_D
  Eigen::SparseMatrix<double> m_H = Get_2nd_H_Matrix(size, lam);
  // initial weight vector
  Eigen::VectorXd v_w = Eigen::VectorXd::Ones(size);
  // measured and smoothed vectors
  Eigen::VectorXd v_y = Eigen::Map<Eigen::VectorXd>(v_trace.data(),size);
  Eigen::VectorXd v_z;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  for(int i = 0; i < max_iter; i++){
    std::cout << "Iteration " << i << "\n";
    auto m_W = Eigen::SparseMatrix<double>(v_w.asDiagonal());
    solver.compute(m_W + m_H);
    v_z = solver.solve(m_W * v_y);
    if(i == max_iter-1) break;
    Eigen::VectorXd v_d = v_y - v_z;
    double mean = 0.0;
    int count = 0;
    std::vector<int> d_minus_i;
    for(int x = 0; x < size; x++){
      if(v_d.coeff(x) >= 0) continue;
      mean += v_d.coeff(x);
      d_minus_i.push_back(count);
      count++;
    }
    mean /= (count-1);
    double var = 0.0;
    for(auto &x : d_minus_i){
      var += (v_d.coeff(x) - mean) * (v_d.coeff(x) - mean);
    }
    double std = sqrt(var/(count-1));
    Eigen::VectorXd v_w_new(size);
    #pragma omp parallel for num_threads(2)
    for(int x = 0; x < size; x++){
      auto coeff = v_d.coeff(x);
      v_w_new.coeffRef(x) = 1.0 / (1.0 + exp(2.0*(coeff-(-mean+2.0*std))/std));
    }
    // checking weight tolerance
    bool early_exit = true;
    for(int x = 0; x < size; x++){
      if(abs(v_w_new.coeffRef(x)-v_w.coeffRef(x))/v_w.coeffRef(x) > tol){
        break;
        early_exit = false;
      }
    }
    v_w = v_w_new;
    if(early_exit) break;
  }

  std::vector<double> output(v_z.data(), v_z.data()+v_z.size());

  if(do_correction){
    #pragma omp parallel for simd schedule(static, 800)
    for(int i = 0; i < size; i++){
      output[i] = v_trace[i] - output[i];
    }
  }

  return output;
}

TraceD ARPLS_Batch(
  TraceD &v_trace,
  const double &lam,
  const int &max_iter,
  const bool &do_correction)
{
  int max_size = v_trace.size();
  int batch_size = 1000;
  int start = 0;
  std::vector<double> output;
  output.reserve(max_size);
  bool iter = true;
  while(iter){
    int stop = start + batch_size;
    if(stop>max_size){
      iter = false;
      stop = max_size;
    }
    // std::cout << "at " << start << " " << stop << "\n";
    std::vector<double> temp_v;
    temp_v.reserve(stop-start);
    std::move(v_trace.begin()+start, v_trace.begin()+stop, std::back_inserter(temp_v));
    auto z = ARPLS(temp_v, lam, max_iter, do_correction);
    output.insert(output.end(), z.begin(), z.end());
    start = stop;
  }

  return output;
}

//==============================================================================
/*
  see https://arxiv.org/pdf/1409.4802.pdf
*/
Eigen::VectorXd PLS_PTRANS_I(
  const Eigen::SparseMatrix<double> &A,
  const Eigen::VectorXd &y)
{
  int size = A.rows();

  double mu[size];
  double alpha[size];
  double beta[size];
  double zeta[size];
  double gamma[size];

  // i = 1
  mu[0] = A.coeff(0,0);
  alpha[0] = A.coeff(0, 1) / mu[0];
  beta[0] = A.coeff(0, 2) / mu[0];
  zeta[0] = y.coeff(0) / mu[0];

  // i = 2
  gamma[1] = A.coeff(1,0);
  mu[1] = A.coeff(1,1) - alpha[0]*gamma[1];
  alpha[1] = (A.coeff(1,2) - beta[0]*gamma[1]) / mu[1];
  beta[1] = A.coeff(1,3) / mu[1];
  zeta[1] = (y.coeff(1) - zeta[0]*gamma[1]) / mu[1];

  /*
  γi = ci − αi−2ei,
  μi = di − βi−2ei − αi−1γi,
  αi = ai−βi−1γi / μi ,
  βi = bi / μi ,
  zi = yi−zi−2 ei−zi−1 γi / μi
  */
  // from 3,4,..., n-2
  for(int i=2; i<size-3; i++){
    gamma[i] = A.coeff(i,i+1) - alpha[i-2]*A.coeff(i,i-2);
    mu[i] = A.coeff(i,i) - beta[i-2]*A.coeff(i,i-2) - alpha[i-1]*gamma[i];
    alpha[i] = (A.coeff(i,i+1) - beta[i-1]*gamma[i])/mu[i];
    beta[i] = A.coeff(i, i+2) / mu[i];
    zeta[i] = (y.coeff(i) - zeta[i-2]*A.coeff(i,i-2) - zeta[i-1]*gamma[i]) / mu[i];
  }

  gamma[size-2] = A.coeff(size-2,size-3) - alpha[size-4]*A.coeff(size-2,size-5);
  mu[size-2] = A.coeff(size-2,size-2) - beta[size-4]*A.coeff(size-2,size-5) - alpha[size-3]*gamma[size-2];
  alpha[size-2] = (A.coeff(size-2,size-1) - beta[size-3]*gamma[size-2]) / mu[size-2];
  gamma[size-1] = A.coeff(size-1,size-2) - alpha[size-3] * A.coeff(size-1,size-3);
  mu[size-1] = A.coeff(size-1,size-1) - beta[size-3] * A.coeff(size-1,size-3) - alpha[size-2]*gamma[size-1];
  zeta[size-2] = (y.coeff(size-2) - zeta[size-3]*(A.coeff(size-2,size-4)-gamma[size-2])) / mu[size-2];
  zeta[size-1] = (y.coeff(size-1) - zeta[size-2]*(A.coeff(size-1,size-3)-gamma[size-1])) / mu[size-1];

  // solving X
  Eigen::VectorXd x(size);
  x.coeffRef(size-1) = zeta[size-1];
  x.coeffRef(size-2) = zeta[size-2] - alpha[size-2]*x.coeffRef(size-1);
  for(int i = size-3; i > 0; i--){
    x.coeffRef(i) = zeta[i] - alpha[i]*x.coeffRef(i+1) - beta[i]*x.coeffRef(i+2);
  }

  return x;
}

TraceD ARPLS_PLS(
  TraceD &v_trace,
  const double &lam,
  const int &max_iter,
  const double &tol,
  const bool &do_correction)
{
  int size = v_trace.size();
  // Eigen::setNbThreads(5);
  // Assuming second order difference matrix diff_D
  Eigen::SparseMatrix<double> m_H = Get_2nd_H_Matrix(size, lam);
  // initial weight vector
  Eigen::VectorXd v_w = Eigen::VectorXd::Ones(size);
  // measured and smoothed vectors
  Eigen::VectorXd v_y = Eigen::Map<Eigen::VectorXd>(v_trace.data(),size);
  Eigen::VectorXd v_z;
  for(int i = 0; i < max_iter; i++){
    // std::cout << "Iteration " << i << "\n";
    auto m_W = Eigen::SparseMatrix<double>(v_w.asDiagonal());
    // std::cout << "before eval \n";
    auto m_A = (m_W + m_H).eval();
    auto v_b = (m_W * v_y).eval();
    // std::cout << "after eval \n";
    // std::cout << "before PLS \n";
    v_z = PLS_PTRANS_I(m_A, v_b);
    // std::cout << "after PLS \n";
    if(i == max_iter-1) break;
    Eigen::VectorXd v_d = v_y - v_z;
    double mean = 0.0;
    int count = 0;
    std::vector<int> d_minus_i;
    for(int x = 0; x < size; x++){
      if(v_d.coeff(x) >= 0) continue;
      mean += v_d.coeff(x);
      d_minus_i.push_back(count);
      count++;
    }
    mean /= (count-1);
    double var = 0.0;
    for(auto &x : d_minus_i){
      var += (v_d.coeff(x) - mean) * (v_d.coeff(x) - mean);
    }
    double std = sqrt(var/(count-1));
    double mean_sig = 2.0*std - mean;
    double std_inv = 2.0/std;
    Eigen::VectorXd v_w_new(size);
    // #pragma omp parallel for simd schedule(static, 512)
    for(int x = 0; x < size; x++){
      v_w_new.coeffRef(x) = 1.0 / (1.0 + exp(std_inv*(v_d.coeff(x)-mean_sig)));
    }
    // checking weight tolerance
    bool early_exit = true;
    for(int x = 0; x < size; x++){
      if(abs(v_w_new.coeffRef(x)-v_w.coeffRef(x))/v_w.coeffRef(x) > tol){
        break;
        early_exit = false;
      }
    }
    v_w = v_w_new;
    if(early_exit) break;
  }

  std::vector<double> output(v_z.data(), v_z.data()+v_z.size());

  if(do_correction){
    for(int i = 0; i < size; i++){
      output[i] = v_trace[i] - output[i];
    }
  }

  return output;
}

//==============================================================================
TraceD NoiseMedian(
  const TraceD &v_trace,
  const int &window_size,
  const double &sigma,
  const bool &do_correction)
{
  double gaus_sigma = sigma;
  if(gaus_sigma < 0){
    gaus_sigma = sqrt((window_size - 1.0)/12.0);
  }
  TraceD median = Utils::MedianFilter<double>(v_trace, window_size);
  TraceD gaus = Utils::GaussianKernel(window_size, 1, 0, gaus_sigma);

  int size = v_trace.size();
  int half_lw = window_size / 2;
  int half_rw = window_size - half_lw;

  int lw_stop = half_lw;
  int rw_stop = size-half_rw;

  std::vector<double> baseline(size);

  std::vector<double> lbound(half_lw + window_size);
  std::vector<double> rbound(half_rw + window_size);

  // extending left boundaries by mirroring.
  for(std::size_t i = 0; i < half_lw + window_size; i++){
    std::size_t j;
    if(i < half_lw){
      j = half_lw - i;
    } else {
      j = i - half_lw;
    }
    lbound[i] = median[j];
  }
  for(std::size_t i = 0; i < half_lw; i++){
    double sum = 0.0;
    for(std::size_t j = 0; j<window_size; j++){
      sum += lbound[i+j] * gaus[j];
    }
    baseline[i] = sum;
  }

  // convoluting central region
  for(std::size_t i = half_lw; i < rw_stop; i++){
    double sum = 0.0;
    // for(int j=i-half_lw; j<i+half_rw; j++){
    //   // if((i-j)<0) continue;
    //   sum += median[j] * gaus[i-j];
    // }
    std::size_t median_start = i-half_lw;
    for(std::size_t j = 0; j < window_size; j++){
      // if((i-j)<0) continue;
      sum += median[median_start+j] * gaus[j];
    }
    baseline[i] = sum;
  }

  // extending right boundaries by mirroring.
  for(std::size_t i = 0; i < half_rw + window_size; i++){
    std::size_t j;
    if(i < window_size){
      j = i + size - window_size - 1;
    } else {
      j = size - (i - window_size) - 1;
    }
    rbound[i] = median[j];
  }
  for(std::size_t i = 0; i < half_rw; i++){
    double sum = 0.0;
    for(std::size_t j = 0; j<window_size; j++){
      sum += rbound[i+j] * gaus[j];
    }
    baseline[i + size - half_rw] = sum;
  }


  // perform correction and return corrected waveform.
  if(do_correction){
    for(std::size_t i = 0; i < size; i++){
      baseline[i] = v_trace[i] - baseline[i];
    }
  }

  return baseline;
}

}
