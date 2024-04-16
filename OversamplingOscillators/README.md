OversamplingOscillators (and BandLimited as well)
Sam Pluta

Optional Dependencies: BufFFT Plugins and SignalBox Quark

Oversampling Oscillators use the VariableOversampling class from Jatin Chowdhury's ChowDSP library. That source is included with the download. 

Standard Oversampling Oscillators include SinOscOS, TriOS, SawOS, SquareOS, VarSawOS, PMOscOS, PM7OS, FM7OS, FM7aOS, and FM7bOS.

Oversampled WAVETABLE Oscillators include OscOS and OscOS3. It is highly recommended to download the BufFFT library to be able to construct 3d wavetables for OscOS3 on the NRT server.

Oversampled Waveshapers include ShaperOS, ShaperOS2, BuchlaFoldOS, and SergeFoldOS.

Bandlimited oscillators are based on Julius Smith's SawN, etc from "Alias-Suppressed Oscillators based on Differentiated Polynomial Waveforms", Vesa Valimaki, Juhan Nam, Julius Smith, and Jonathan Abel.

Bandlimited oscillators include SawBL, SquareBL, TriBL, and ImpulseBL.

These are luxurious oscillators that sound great, and have some inefficient code used to achieve these great sounds. They are not designed to be super efficient, so if you want that, use the standard SC oscillators.


run the following from this directory to build from source using cmake

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSC_PATH=<PATH TO SC SOURCE> 
cmake --build . --config Release
```

It should build . The OS plugins are Oversampled. The BL plugins are polynomial bandlimited.

After building ake sure this directory or all of the scx, sc, and schelp files are in the SC path, recompile the SC libary, and they should work. 

If downloading a pre-build binary, mac users will need to unquarantine. In the terminal, run 'xattr -cr <the_oversampling_oscillators_directory>'

#OversamplingOscillators
