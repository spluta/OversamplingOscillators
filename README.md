OversamplingOscillators (and BandLimited as well)
Sam Pluta

Oversampling Oscillators use the VariableOversampling class from Jatin Chowdhury's ChowDSP library. That source is included with the download. Bandlimited oscillators are based on Julius Smith's SawN, etc from "Alias-Suppressed Oscillators based on Differentiated Polynomial Waveforms", Vesa Valimaki, Juhan Nam, Julius Smith, and Jonathan Abel.

These are luxurious oscillators that sound great, and have some inefficient code used to achieve these great sounds. They are not designed to be super efficient, so if you want that, use the standard SC oscillators.


run the following from this directory to build from source using cmake

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSC_PATH=<PATH TO SC SOURCE> 
cmake --build . --config Release
```

It should build SawBL, SquareBL, TriBL, ImpulseBL, SinOscOS, TriOS, SawOS, SquareOS, VarSawOS, PMOscOS, PM7OS, FM7OS. The OS plugins are Oversampled. The BL plugins are polynomial bandlimited.

After building ake sure this directory or all of the scx, sc, and schelp files are in the SC path, recompile the SC libary, and they should work. #OversamplingOscillators
