#include "waveformMethods/waveformMethods.hpp"

waveform_methods::WavePoint waveform_methods::_CalcCFDTime(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &start_i,
  const double &frac_max,
  const int &trace_size,
  const EdgeType &type)
{

  int start, end;
  if(type == EdgeType::Rise){
    end = -1;
    start = start_i > -1 ? start_i : 0;
    start = start == trace_size ? start : start - 1;
  }
  else if(type == EdgeType::Fall){
    end = trace_size - 1;
    start = start_i < trace_size - 1 ? start_i : trace_size - 1;
  }

  bool success = false;
  int frac_i = 0;
  double frac_t_value = 0.0;
  double frac_v_value = 0.0;
  while(true){
    if(start == end) break;
    double v_value = v_trace.at(start);
    if(v_value <= frac_max){
      frac_i = start;
      frac_v_value = v_value;
      success = true;
      break;
    }
    if(type == EdgeType::Rise) start--;
    else if(type == EdgeType::Fall) start++;
  }

  if(!success) return WavePoint{-1, frac_max, frac_t_value};

  if(type == EdgeType::Rise){
    frac_t_value = waveform_methods::LinearInterpolationX(
      t_trace.at(frac_i), v_trace.at(frac_i),
      t_trace.at(frac_i+1), v_trace.at(frac_i+1),
      frac_v_value);
      return WavePoint{frac_i, frac_max, frac_t_value};
  }
  else {
    frac_t_value = waveform_methods::LinearInterpolationX(
      t_trace.at(frac_i-1), v_trace.at(frac_i-1),
      t_trace.at(frac_i), v_trace.at(frac_i),
      frac_v_value);
      return WavePoint{frac_i, frac_max, frac_t_value};
  }
}

//==============================================================================
double waveform_methods::CalcCFDTime(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &max_i,
  const double &frac)
{
  double frac_value = v_trace.at(max_i) * frac;
  int trace_size = v_trace.size();
  auto wav_pt = _CalcCFDTime(
    v_trace, t_trace, max_i, frac_value, trace_size, EdgeType::Rise);

  return wav_pt.t;
}

//==============================================================================
std::vector<double> waveform_methods::CalcCFDTime(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &max_i,
  const double start_frac,
  const double incr_size)
{
  std::vector<double> cfd_times = {};
  int incr = 0;
  int trace_size = v_trace.size();
  double max_value = v_trace.at(max_i);
  double frac_0 = start_frac;
  while(frac_0 < 1.0){
    double frac_value = max_value * frac_0;
    auto wav_pt = _CalcCFDTime(
      v_trace, t_trace, max_i, frac_value, trace_size, EdgeType::Rise);
    cfd_times.push_back(wav_pt.t);
    frac_0 += incr_size;
  }

  return cfd_times;
}

//==============================================================================
double waveform_methods::CalcFWHM(
  const TraceD &v_trace,
  const TraceD &t_trace,
  const int &max_i,
  const double &frac)
{
  double frac_value = frac * v_trace.at(max_i);
  int trace_size = v_trace.size();
  auto wav_pt_rise = _CalcCFDTime(
    v_trace, t_trace, max_i, frac_value, trace_size, EdgeType::Rise);
  auto wav_pt_fall = _CalcCFDTime(
    v_trace, t_trace, max_i, frac_value, trace_size, EdgeType::Fall);
  return wav_pt_fall.t - wav_pt_rise.t;
}
