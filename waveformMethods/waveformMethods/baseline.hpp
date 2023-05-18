#ifndef WAVEFOMR_METHODS_BASELINE_H
#define WAVEFOMR_METHODS_BASELINE_H

#include "waveformMethods/waveformMethods.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace waveform_methods::Baseline
{

  const Eigen::SparseMatrix<double> &Get_2nd_D_Matrix(const int &size)
  {
    static Eigen::SparseMatrix<double> m_2nd_D = [](const int &size){
      Eigen::SparseMatrix<double> d(size-2, size);
      for(int i = 0; i < size-2; i++){
        d.coeffRef(i, i) = 1.0;
        d.coeffRef(i, i+1) = -2.0;
        d.coeffRef(i, i+2) = 1.0;
      }
      return d;
    }(size);
    return m_2nd_D;
  };

  const Eigen::SparseMatrix<double> &Get_2nd_H_Matrix(
      const int &size,
      const double &lam)
  {
    static Eigen::SparseMatrix<double> m_2nd_H = [](int size, double lam){
      Eigen::SparseMatrix<double> d(size-2, size);
      for(int i = 0; i < size-2; i++){
        d.coeffRef(i, i) = 1.0;
        d.coeffRef(i, i+1) = -2.0;
        d.coeffRef(i, i+2) = 1.0;
      }
      return (d.transpose() * d * lam).eval();
    }(size, lam);
    return m_2nd_H;
  };

  //============================================================================
  double CalcBaseline(
    const TraceD &v_trace,
    const int &start,
    const int &end);

  double CalcBaseline(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &tmin,
    const double &tmax);

  bool MultiSignalBaselineCorrection(
    WavePoints &signals,
    TraceD &v_trace,
    const double &frac_npts,
    const int &npts_forward,
    const int &npts_backward);

  bool MultiSignalBaselineCorrection(
    WavePoints &signals,
    TraceD * const v_trace,
    const double &frac_npts,
    const int &npts_forward,
    const int &npts_backward);

  bool MultiSignalBaselineCorrection(
    WavePoints &signals,
    TraceD &v_trace,
    const TraceD &t_trace,
    const double &frac_npts,
    const double &forward_time,
    const double &backward_time);

  bool MultiSignalBaselineCorrection(
    WavePoints &signals,
    TraceD * const v_trace,
    const TraceD * const t_trace,
    const double &frac_npts,
    const double &forward_time,
    const double &backward_time);

  TraceD ARPLS(
    TraceD &v_trace,
    const double &lam = 1.0e5,
    const int &max_iter = 50,
    const double &tol = 1.0e-3,
    const bool &do_correction=true);

  TraceD ARPLS_PLS(
    TraceD &v_trace,
    const double &lam = 1.0e5,
    const int &max_iter = 50,
    const double &tol = 1.0e-3,
    const bool &do_correction=true);

  TraceD ARPLS_Batch(
    TraceD &v_trace,
    const double &lam = 1.0e5,
    const int &max_iter = 10,
    const bool &do_correction=true);

  TraceD NoiseMedian(
    const TraceD &v_trace,
    const int &window_size,
    const double &sigma=-1.0,
    const bool &do_correction=true);
};

#endif
