OversamplingOscillators
Sam Pluta

uses the VariableOversampling class from Jatin Chowdhury's ChowDSP library in the ported plugins library by Mads Kjeldgaard

to set this up to edit in .vscode, edit the c_cpp_properties.json (in .vscode hidded folder) to correctly point to these files on your system

run the following from this directory to build from source using cmake

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSC_PATH=<PATH TO SC SOURCE> -DPORTED_PATH=<PATH TO PORTED PLUGINS SOURCE>
cmake --build . --config Release

It should build OSaw2.scx and OVarSaw2.scx. Make sure these two files and OversamplingOscillators.sc are in the SC path, recompile the SC libary, and they should work.# OversamplingOscillators
