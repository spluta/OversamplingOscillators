#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include <array>

namespace BuchlaFoldOS {

class BuchlaFoldOS : public SCUnit {
public:
  BuchlaFoldOS();

  // Destructor
  ~BuchlaFoldOS();
  VariableOversampling<> oversample;
  float m_oversampling_ratio;

private:
  // Calc function
  float next(float sig, float amp);
  float next_os(float sig, float amp);
  void next_aa(int nSamples);
  float buchla_cell(float sig, float sign, float thresh, float sig_mul1, float sign_mul, float sig_mul2);

  enum InputParams { Sig, Amp, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;

  int m_oversamplingIndex{0};
};

} // namespace BuchlaFoldOS

namespace SergeFoldOS {

class SergeFoldOS : public SCUnit {
public:
  SergeFoldOS();

  // Destructor
  ~SergeFoldOS();
  VariableOversampling<> oversample;
  std::array<double, 1000> sergeWavetable;
  float m_oversampling_ratio;

private:
  // Calc function
  float next(float sig, float amp);
  float next_os(float sig, float amp);
  void next_aa(int nSamples);

  enum InputParams { Sig, Amp, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;

  int m_oversamplingIndex{0};
};

} // namespace SergeFoldOS

namespace BufUnit {
  class BufUnit : public SCUnit {
    public:

    BufUnit();

    // Destructor
    ~BufUnit();
      float m_fbufnum;
      SndBuf* m_buf;
    
      bool GetTable(World* world, float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData, int& tableSize);
    private:
};
} // namespace BufUnit

namespace ShaperOS {

class ShaperOS : public SCUnit {
public:
  ShaperOS();

  // Destructor
  ~ShaperOS();
  VariableOversampling<> oversample;
  //std::array<double, 1000> sergeWavetable;
  float* m_inbuf;
  float m_fbufnum;
  SndBuf* m_buf;
  BufUnit::BufUnit buf_unit;
  float m_oversampling_ratio;
  

private:
  // Calc function
  void next_aa(int nSamples);
  bool GetTable(float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData,
                               int& tableSize);
  //float force_inline ShaperOS::ShaperPerform(const float* table0, const float* table1, float in, float offset, float fmaxindex);

  float next_os(const float* table0, float in, float fmaxindex);

  float Perform(const float* table0, float in, float fmaxindex);

  enum InputParams { BufNum, Sig, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;

  int m_oversamplingIndex{0};
};

} // namespace ShaperOS

namespace OscOS {

class OscOS : public SCUnit {
public:
  OscOS();

  // Destructor
  ~OscOS();
  VariableOversampling<> oversample;
  VariableOversampling<> upsample_buf_loc;
  float* m_inbuf;
  float m_fbufnum;
  SndBuf* m_buf;
  BufUnit::BufUnit buf_unit;
  float m_last_phase;
  float m_last_buf_loc;
  float m_oversampling_ratio;
  

private:
  // Calc function
  void next_aa(int nSamples);
  bool GetTable(float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData,
                               int& tableSize);
  //float force_inline OscOS::ShaperPerform(const float* table0, const float* table1, float in, float offset, float fmaxindex);

  float next_os(const float* table0, float phase, float buf_divs, float buf_loc, int table_size, float fmaxindex);

  float Perform(const float* table0, float phase, float buf_divs, float buf_loc, int table_size, float fmaxindex);

  enum InputParams { BufNum, Phase, BufDivs, BufLoc, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;
  float *upsample_buf;

  int m_oversamplingIndex{0};
};

} // namespace OscOS