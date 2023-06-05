OversamplingOscillators
Sam Pluta

uses the VariableOversampling class from Jatin Chowdhury's ChowDSP library. that source is included with the download.


run the following from this directory to build from source using cmake

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSC_PATH=<PATH TO SC SOURCE> 
cmake --build . --config Release

It should build TriBL, SquareBL, SawBL, SinOscOS, TriOS, SawOS, VarSawOS, PMOscOS, PM7OS, FM7OS. The OS plugins are Oversampled. The BL plugins are polynomial bandlimited.

After building ake sure this directory or all of the scx, sc, and schelp files are in the SC path, recompile the SC libary, and they should work. #OversamplingOscillators
