// OversamplingOscillators.hpp

#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include "dcblocker.h"

namespace SawOS {

class SawOSNext {
  public:
    SawOSNext();
    ~SawOSNext();
    void next(float freq, float *m_phase, float *m_freqMul, float *osBuffer, float overSamplingRatio);
  private:

};

class SawOS : public SCUnit {
public:
  SawOS();

  // Destructor
  ~SawOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);
  //float next(float freq);
  SawOSNext saw;

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{0.f};
  float m_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace SawOS

namespace TriOS {
class TriOS : public SCUnit {
public:
  TriOS();

  // Destructor
  ~TriOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{0.f};
  float m_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};
} //namespace TriOS

namespace VarSawOS {

class VarSawOS : public SCUnit {
public:
  VarSawOS();

  // Destructor
  ~VarSawOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  float next(float freq, float width);
  void next_aa(int nSamples);
  void next_ak(int nSamples);
  void next_ka(int nSamples);
  void next_kk(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, Width, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{in0(Freq)};
  float m_phase{(float)in0(Phase)};
  float m_width{in0(Width)};
  float m_invwidth{2.f/m_width};
  float m_inv1width{2.f/(1.f - m_width)};
  float m_freqMul{1.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace VarSawOS

namespace SquareOS {

class SquareOS : public SCUnit {
public:
  SquareOS();

  // Destructor
  ~SquareOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  float next(float freq, float width);
  void next_aa(int nSamples);
  void next_ak(int nSamples);
  void next_ka(int nSamples);
  void next_kk(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, Width, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{in0(Freq)};
  float m_phase{(float)in0(Phase)};
  float m_width{in0(Width)};
  float m_freqMul{1.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace SquareOS


namespace SawBL {

class Diff {       // The class

  public:
    Diff();
    ~Diff();
    
    float diff(float sample, float p0n);

  private:             // Access specifier
    float lastSample{0.f};        // Attribute (int variable)
    
};

class SawBLNext {
  public:
    SawBLNext();
    ~SawBLNext();
    float next(float freq, float* phase, int* counter, float sampleRate, float freqMul);

  private:
    // Diff diff6_1;
    // Diff diff6_2;
    // Diff diff6_3;
    // Diff diff6_4;
    // Diff diff6_5;

    Diff diff4_1;
    Diff diff4_2;
    Diff diff4_3;

    Diff diff2_1;
};

class SawBL : public SCUnit {
public:
  
  SawBL();

  // Destructor
  ~SawBL();

private:
  // Calc function
  //float next(float freq);
  void next_a(int nSamples);
  void next_k(int nSamples);
  //float diff(float x, int num);

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  //float m_p0n = sampleRate()/m_freq; 
  float m_freqMul{2.0f/(float)sampleRate()};

  //const float samplerate{(float)sampleRate()};

  int m_counter{0};

  float diffArray [3] = {m_phase, m_phase, m_phase};

  SawBLNext saw;

};

} // namespace SawBL

namespace SquareBL {

class SquareNext {
  public:
    SquareNext();
    ~SquareNext();
    float next(float freq, float duty);
    void setRatePhase(float rateIn, float phaseIn);

  //const float samplerate{(float)sampleRate()};

  private:
    //this is a bit luxurious. at 96K it can go down to 11.75hz
    SawBL::SawBLNext saw;
    int arraySize = 4096;
    float delArray [4096] = {1.f};
    int delArrayCounter{0};
    int delMax{arraySize-1};
    int m_counter{0};
    float m_sampleRate{96000.f};
    float m_phase{0.f};
    float m_freqMul{2.0f/m_sampleRate};
    float m_fmin{96000.f/float(2*delMax)};
    
};


class SquareBL : public SCUnit {
public:
  SquareBL();

  // Destructor
  ~SquareBL();

private:
  // Calc function
  void next_aa(int nSamples);
  void next_ak(int nSamples);
  void next_ka(int nSamples);
  void next_kk(int nSamples);

  
  
  //float diff(float x, int num);

  enum InputParams { Freq, Phase, Width, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  float m_sampleRate{(float)sampleRate()};

  SquareNext square;
};

} // namespace SquareBL

namespace TriBL {

  class TriBL : public SCUnit {
    public:
  TriBL();

  // Destructor
  ~TriBL();

private:
  // Calc function
  float next(float freq, float duty);
  void next_aa(int nSamples);
  void next_ak(int nSamples);
  void next_ka(int nSamples);
  void next_kk(int nSamples);
  
  //float diff(float x, int num);

  enum InputParams { Freq, Phase, Width, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  float m_sampleRate{(float)sampleRate()};

  float m_delay1{0.f};

  //dcblocker::Dcblocker blocker;

  SquareBL::SquareNext square;
  };

} // namespace TriBL

namespace ImpulseBL {

class Impulse {
  public:
    Impulse();
    ~Impulse();
    float next ();
    
  private:
    float impArray[2] = {0.f};
};

class ImpulseBL : public SCUnit {
public:
  ImpulseBL();
  // Destructor
  ~ImpulseBL();

private:
  // Calc function
  
  void next_a(int nSamples);
  void next_k(int nSamples);

  enum InputParams { Freq, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  int m_counter{0};
  float m_phase{1.f};
  float freqMul{2.0f/(float)sampleRate()};

  float m_delay1 = {0.f};

  Impulse impulse;

  SawBL::SawBLNext saw;
};
}