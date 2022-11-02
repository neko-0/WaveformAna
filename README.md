# Installation
`mkdir TestArea; cd TestArea`

`git clone https://github.com/neko-0/WaveformAna.git`

`mkdir -p build`

`cd WaveformAna; source setup.sh`

`build`

# Running Analysis

`run_Ana --selector SSRL --directory path_to_input_files`

where
  - `--selector, -s` is the name of the user analysis
  - `--directory, -d` is the path to the `ROOT` files.

# Creating your own analysis

New user analysis need to be registered in order to be used with `run_Ana`.

The user analysis need to be derived from class `BaseAna`, and user needs to
implement the following virtual member functions:
  - `void initialize(BetaConfigMgr *configMgr)`

    call only once for initialization purpose.

  - `void execute(BetaConfigMgr *configMgr)`

    call for every waveform event. Analysis on each event takes place within this method.

  - `void finialize(BetaConfigMgr *configMgr)`

    call only once for cleanup purpose.

To register the analysis, call the `AnalysisRegister` help function.
e.g. in the analysis header file, include a line

`AnalysisRegister<MyNewAna> reg("My_Analysis");`

then you should be able to selector the analysis with

`run_Ana --selector My_Analysis --directory path_to_input_files`
