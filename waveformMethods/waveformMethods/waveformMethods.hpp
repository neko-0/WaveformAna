#ifndef WAVEFORMMETHODS_H
#define WAVEFORMMETHODS_H

#include <vector>
#include <omp.h>

namespace waveform_methods
{

  struct WavePoint{
      int index;
      double v; // v(oltage) value
      double t; // t(ime) value
  };

  struct MixParamsSet1{
    double baseline;
    double rms;
    int neg_max_i;
    int pos_max_i;
    double neg_max;
    double pos_max;
  };

  enum EdgeType{ Rise, Fall };

  typedef std::vector<WavePoint> WavePoints;
  typedef std::vector<double> TraceD;

  //============================================================================
  WavePoint FindSignalMax(const TraceD &v_trace, const TraceD &t_trace);
  WavePoint FindSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &min, const double &max);
  WavePoint FindSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start, const int &end);

  //============================================================================
  WavePoints FindIdenticalSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace);
  WavePoints FindIdenticalSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start, const int &end);
  WavePoints _FindIdenticalSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start, const int &end,
    const double &baseline=0.0);
  WavePoints FindMultipleSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &threshold,
    const double &scale = 2.0);
  WavePoints _FindMultipleSignalMax(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start_i,
    const int &end_i,
    const double &threshold,
    const double &scale = 2.0);

  //============================================================================
  WavePoints _FindThresholdPoints(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &threshold);

  //============================================================================
  std::vector<double> FindTimeAtThreshold(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start,
    const int &end,
    const double &threshold);
  std::vector<double> FindTimeAtThreshold(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &threshold);
  std::vector<double> FindTimeAtThreshold(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &tmin,
    const double &tmax,
    const double &threshold);

  //============================================================================
  double CalcRiseTime(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start_index=-1,
    const double &low=0.1,
    const double &high=0.9,
    const EdgeType &type=EdgeType::Rise);
  double CalcFallTime(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start_index=-1,
    const double &low=0.1,
    const double &high=0.9);

  //============================================================================
  double LinearInterpolationX(
    const double &x1, const double &y1,
    const double &x2, const double &y2,
    const double &y);

  //============================================================================
  double CalcNoise(const TraceD &v_trace, const int &npts);
  double CalcNoise(const TraceD &v_trace, const double &frac);
  double CalcNoise(
    const TraceD &v_trace,
    const int &imin, const int &imax);
  double CalcNoise(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &tmin, const double &tmax);

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
    TraceD &v_trace,
    const TraceD &t_trace,
    const double &frac_npts,
    const double &forward_time,
    const double &backward_time);

  //============================================================================
  double CalcPulseArea(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &max_index, const double &baseline= 0.0);
  double CalcPulseArea(
      const TraceD &v_trace,
      const TraceD &t_trace);
  double CalcPulseArea(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start_index, const int &end_index);
  double CalcPulseArea(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const double &time_start, const double &time_end);
  double _CalcPulseArea(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &lstart,
    const int &rstart,
    const double &baseline= 0.0);

  //============================================================================
  std::vector<double> CalcCFDTime(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &max_i,
    const double start_frac,
    const double incr_size);
  double CalcCFDTime(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &max_i,
    const double &frac);
  WavePoint _CalcCFDTime(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &start_i,
    const double &frac_max,
    const int &trace_size,
    const EdgeType &type);

  //============================================================================
  double CalcFWHM(
    const TraceD &v_trace,
    const TraceD &t_trace,
    const int &max_i,
    const double &frac=0.5);

  //============================================================================
  MixParamsSet1 CalcMaxNoiseBase(
    const TraceD &v_trace,
    const int &start,
    const int &end);
  MixParamsSet1 CalcMaxNoiseBase(
    const TraceD &v_trace,
    const double &frac);
};

#endif
