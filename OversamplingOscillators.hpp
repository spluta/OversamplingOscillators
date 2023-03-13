// OversamplingOscillators.hpp

#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include "../dcblocker.h"

namespace SawOS {

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

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{0.f};
  float m_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace SawOS

namespace VarSawOS {

class VarSawOS : public SCUnit {
public:
  VarSawOS();

  // Destructor
  ~VarSawOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_aa(int nSamples);
  void next_ak(int nSamples);
  void next_ka(int nSamples);
  void next_kk(int nSamples);

  enum InputParams { Freq, Phase, Width, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{in0(Freq)};
  double m_phase{(double)in0(Phase)};
  float m_width{in0(Width)};
  float m_invwidth{2.f/m_width};
  float m_inv1width{2.f/(1.f - m_width)};
  float m_freqMul{1.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace VarSawOS


namespace SawPn {

class Diff {       // The class

  public:
    Diff();
    ~Diff();
    
    float diff(float sample, float p0n);

  private:             // Access specifier
    float lastSample{0.f};        // Attribute (int variable)
    
};

class SawPn : public SCUnit {
public:
  SawPn();

  // Destructor
  ~SawPn();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);
  //float diff(float x, int num);

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  float m_poly2{in0(Phase)};
  float m_poly4{in0(Phase)};
  //float m_p0n = sampleRate()/m_freq; 
  float m_freqMul{2.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};

  Diff diff4_1;
  Diff diff4_2;
  Diff diff4_3;

  Diff diff2_1;

  int m_counter{0};

  float diffArray [3] = {m_phase, m_phase, m_phase};

};

} // namespace SawOS