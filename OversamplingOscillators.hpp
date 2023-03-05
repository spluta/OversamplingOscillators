// OversamplingOscillators.hpp

#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include "../dcblocker.h"

namespace OSaw2 {

class OSaw2 : public SCUnit {
public:
  OSaw2();

  // Destructor
  ~OSaw2();
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

} // namespace OSaw2

namespace OVarSaw2 {

class OVarSaw2 : public SCUnit {
public:
  OVarSaw2();

  // Destructor
  ~OVarSaw2();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);

  enum InputParams { Freq, Phase, Duty, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
  dcblocker::Dcblocker dcfilter;

  float m_freq_past{in0(Freq)};
  double m_phase{(double)in0(Phase)};
  float m_duty{in0(Duty)};
  float m_invduty{2.f/m_duty};
  float m_inv1duty{2.f/(1.f - m_duty)};
  float m_freqMul{1.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace OVarSaw2