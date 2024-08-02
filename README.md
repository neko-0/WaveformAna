# Installation

```bash
git clone https://github.com/neko-0/WaveformAna.git
source WaveformAna/setup.sh # for next login, just source the same setup.sh
build
```

NOTE: if you have problem compiling with the <std::filesystem>, try to update
gcc/g++ to version 11.

# Running Analysis

`run_Ana --selector SSRL --directory path_to_input_files`

where
  - `--selector, -s` is the name of the user analysis
  - `--directory, -d` is the path to the `ROOT` files.

optional:
  - `--config, -c` path to configuration file. User need to implement actual usage in the analysis routine.
  - `--nfile, -n` number of files for MP.
  - `--mp` turn on multiprocessing.


# Creating your own analysis

New user analysis need to be registered in order to be used with `run_Ana`.

First create a directory under `WaveformAna`, e.g `MyNewAna` with the following
structure:

```.
WaveformAna
 ├── CMakeLists.txt
 ├── MyNewAna
     ├── MyNewAna.hpp
     ├── src
         ├──MyNewAna.cpp
```

The user analysis need to be derived from class `BaseAna`, so in `MyNewAna.hpp`:

```cpp
#include "baseAna/baseAna.hpp"
#include "configMgr/betaConfigMgr.hpp"

struct MyNewAna : BaseAna {
  MyNewAna(){};
  ~MyNewAna(){};
  virtual void setup(BetaConfigMgr *configMgr);
  virtual void initialize(BetaConfigMgr *configMgr); // call only once for initialization purpose.
  virtual bool execute(BetaConfigMgr *configMgr); // call for every waveform event.
  virtual void finalize(BetaConfigMgr *configMgr); // call only once for cleanup purpose.
};
```

and user needs to implement the virtual member functions:

The `CMakelist.txt` should have something similar to this:

```cmake
cmake_minimum_required(VERSION 3.10)
project(MyNewAna VERSION 1.0 DESCRIPTION "MyNewAna")
add_library(MyNewAna SHARED src/MyNewAna.cpp)
target_include_directories(MyNewAna PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(MyNewAna PRIVATE BaseAna)
target_link_libraries(MyNewAna PRIVATE ConfigMgr)
```

To register the analysis, do the following in `analysisDriver/src/register.cpp`

```cpp
#include "analysisDriver/register.hpp"

// include analysis header
#include "AnaSSRL/AnaSSRL.hpp"
#include "MyNewAna/MyNewAna.hpp" // <---- here is the new analysis


bool AnalysisRegister::Run(){
  auto _ssrl_reg = AnalysisRegister::Register<AnaSSRL>("SSRL");
  auto _mynewana_reg = AnalysisRegister::Register<MyNewAna>("My_Analysis");

  return true;
}
```

and in the `analysisDriver/CMakeLists.txt`, add the packaged name that you have
named in your analysis

```cmake
# Need to include user analysis here as well.
set(ANA_LIST AnaSSRL MyNewAna)
```

then you should be able to selector the analysis with

`run_Ana --selector My_Analysis --directory path_to_input_files`
