#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include "SC_PlugIn.h"
#include "SC_PlugIn.hpp"
#include <array>

namespace ShaperOS {

class ShaperOS : public SCUnit {
public:
  ShaperOS(SCUnit *unit);

  // Destructor
  ~ShaperOS();
  VariableOversampling<> oversample;
  //std::array<double, 1000> sergeWavetable;
  float* m_inbuf;
  float m_prev_in;
  float m_fbufnum;
  SndBuf* m_buf;

private:
  // Calc function
  float next(float sig, float amp);
  float next_os(float sig, float amp);
  void next_aa(int nSamples);
  float force_inline ShaperOS::ShaperPerform(const float* table0, const float* table1, float in, float offset, float fmaxindex);

  enum InputParams { BufNum, Sig, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;

  int m_oversamplingIndex{0};
};

} // namespace ShaperOS