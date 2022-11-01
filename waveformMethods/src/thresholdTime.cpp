

/*==============================================================================
Find the time at where the voltage crosses the threshold level.

  param thresholdLevel := the threshold level
  param voltage        := the vector of voltage value of the waveform
  param time           := the vector of time value of the waveform
  param pmax_holder    := a std pair of the Pmax value and its index in the wave

  return : the time at the threshold.

/=============================================================================*/
double WaveformAnalysis::Find_Time_At_Threshold(
    const double &thresholdLevel,
    const std::vector<double> &voltageVec,
    const std::vector<double> &timeVec,
    const std::pair<double, unsigned int> &Pmax
)
{
    double timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

    unsigned int timeBelowThreshold_index = 0;

    unsigned int pmax_index = Pmax.second;
    double pmax = Pmax.first;
    std::size_t npoints = voltageVec.size();

    if (pmax_index == npoints - 1)
        pmax_index = pmax_index - 1; // preventing out of range

    if (pmax < thresholdLevel)
    {
        return 9999.0;
    }
    else
    {
        for (int i = pmax_index; i > -1; i--)
        {
            if (voltageVec.at(i) <= thresholdLevel)
            {
                timeBelowThreshold_index = i;

                timeBelowThreshold = timeVec.at(i);

                break;
            }
        }

        timeAtThreshold = xlinearInter( timeBelowThreshold, voltageVec.at(timeBelowThreshold_index), timeVec.at(timeBelowThreshold_index + 1), voltageVec.at(timeBelowThreshold_index + 1), thresholdLevel);

        return timeAtThreshold;
    }
}
