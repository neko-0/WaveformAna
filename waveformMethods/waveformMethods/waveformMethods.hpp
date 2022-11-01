#ifndef WAVEFORMMETHODS_H
#define WAVEFORMMETHODS_H

#include <vector>

namespace waveform_methods
{

  struct WavePoint{
      int index;
      double v; // v(oltage) value
      double t; // t(ime) value
  };

  enum RiseTimeType{ Rise, Fall };

  typedef std::vector<WavePoint> WavePoints;
  typedef std::vector<double> TraceDType;

  //============================================================================
  WavePoint FindSignalMax(const TraceDType &v_trace, const TraceDType &t_trace);
  WavePoint FindSignalMax(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const double &min, const double &max);
  WavePoint FindSignalMax(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start, const int &end);

  //============================================================================
  WavePoints FindIdenticalSignalMax(
    const TraceDType &v_trace,
    const TraceDType &t_trace);
  WavePoints FindIdenticalSignalMax(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start, const int &end);
  WavePoints _FindIdenticalSignalMax(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start, const int &end,
    const double &baseline=0.0);
  WavePoints FindMultipleSignalMax(
    const std::vector<double> &v_trace,
    const std::vector<double> &t_trace,
    const double &threshold,
    const double &scale = 2.0);
  WavePoints _FindMultipleSignalMax(
    const std::vector<double> &v_trace,
    const std::vector<double> &t_trace,
    const int &start_i,
    const int &end_i,
    const double &threshold,
    const double &scale = 2.0);

  //============================================================================
  WavePoints _FindThresholdPoints(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const double &threshold);

  //============================================================================
  double CalcRiseTime(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start_index=-1,
    const double &low=0.1,
    const double &high=0.9,
    const RiseTimeType &type=RiseTimeType::Rise);
  double CalcFallTime(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start_index=-1,
    const double &low=0.1,
    const double &high=0.9);

  //============================================================================
  double LinearInterpolationX(
    const double &x1, const double &y1,
    const double &x2, const double &y2,
    const double &y);

  //============================================================================
  double CalcNoise(const TraceDType &v_trace, const int &npts);
  double CalcNoise(const TraceDType &v_trace, const double &frac);
  double CalcNoise(
    const TraceDType &v_trace,
    const int &imin, const int &imax);
  double CalcNoise(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const double &tmin, const double &tmax);

  //============================================================================
  double CalcPulseArea(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &max_index, const double &baseline= 0.0);
  double CalcPulseArea(
      const TraceDType &v_trace,
      const TraceDType &t_trace);
  double CalcPulseArea(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const int &start_index, const int &end_index);
  double CalcPulseArea(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    const double &time_start, const double &time_end);
  double _CalcPulseArea(
    const TraceDType &v_trace,
    const TraceDType &t_trace,
    int lstart_i, int rstart_i, const double &baseline= 0.0);

    //============================================================================
    double CalcCFDTime(
      const double &frac,
      const int &max_index,
      const std::vector<double> &v_trace,
      const std::vector<double> &t_trace);
};

#endif
