#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include <array>
#include "BufUnit.cpp"
#include "Extras.cpp"


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
  std::array<float, 1000> wavetable;
  float* sergeWavetable;
  float m_oversampling_ratio;
  Extras::ProcessFuncs process_funcs;

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

namespace ShaperOS {

class ShaperOS : public SCUnit {
public:
  ShaperOS();

  // Destructor
  ~ShaperOS();
  VariableOversampling<> oversample;
  float* m_inbuf;
  float m_fbufnum;
  SndBuf* m_buf;
  BufUnit::BufUnit buf_unit;
  float m_oversampling_ratio;
  Extras::ProcessFuncs process_funcs;
  

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

namespace ShaperOS2 {

class ShaperOS2 : public SCUnit {
public:
  ShaperOS2();

  // Destructor
  ~ShaperOS2();
  VariableOversampling<> oversample;
  VariableOversampling<> upsample_buf_loc;
  VariableOversampling<> upsample_input;
  float* m_inbuf;
  float m_fbufnum;
  SndBuf* m_buf;
  BufUnit::BufUnit buf_unit;
  float m_last_phase;
  float m_last_buf_loc;
  float m_oversampling_ratio;
  
  //float Perform(const float* table0, float input, float buf_divs, float buf_loc, int table_size, float fmaxindex);
  Extras::ProcessFuncs process_funcs;

private:
  // Calc function
  void next_aa(int nSamples);

  float next_os(const float* table0, float input, float buf_divs, float buf_loc, int table_size, float fmaxindex);

  
  enum InputParams { BufNum, In, BufDivs, BufLoc, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *osBuffer;
  float *upsample_buf_ptr;
  float *upsample_input_ptr;

  int m_oversamplingIndex{0};
};

} // namespace ShaperOS2


namespace OscOS3 {

class OscOS3 : public SCUnit {
public:
  OscOS3();

  // Destructor
  ~OscOS3();
  VariableOversampling<> oversample;
  VariableOversampling<> upsample_buf_loc;
  VariableOversampling<> upsample_chan_loc;
  VariableOversampling<> upsample_phase_buf_loc;

  float* m_inbuf;
  float m_fbufnum;

  SndBuf* m_buf;
  BufUnit::BufUnit buf_unit;
  BufUnit::BufUnit phase_buf_unit;

  float m_freq_past{0.f};
  float i_phase{in0(Phase)};
  double m_freqMul{2.0/sampleRate()};

  float m_last_phase;
  float m_last_buf_loc;
  double m_oversampling_ratio;
  float m_sync_trig;
  Extras::SawOSNext saw;
  Extras::ProcessFuncs process_funcs;

private:
  // Calc function
  void next_aa(int nSamples);
  float next_os(const float* buf_data, const float* phase_buf_data, const float freq, const float phase, float buf_divs, float buf_loc, float num_chans, float chan_loc, float phase_buf_divs, float phase_buf_loc, int each_table_size, float fmaxindex, int each_phase_table_size, float phase_fmaxindex);
  
  float Perform(const float* buf_data, const float* phase_buf_data, float phase, float buf_divs, float buf_loc, float num_chans, float chan_loc, float phase_buf_divs, float phase_buf_loc, int table_size, float fmaxindex, int phase_table_size, float phase_fmaxindex);

  enum InputParams { BufNum, PhaseBuf, Freq, Phase, SyncTrig, BufDivs, BufLoc, NumChans, ChanLoc, PhaseBufDivs, PhaseBufLoc, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float *os_buffer;
  float *os_buf_loc;
  float *os_chan_loc;
  float *os_phase_buf_loc;

  int m_oversampling_index{0};
};

} // namespace OscOS3

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
    int counter{0};
  
    SincExtras::SincFuncs process_funcs;  
    
  private:
    // Calc function
    void next_aa(int nSamples);
  
    float next_os(const float* table0, float phase, float buf_divs, float buf_loc, int table_size, float fmaxindex);
  
    enum InputParams { BufNum, Phase, BufDivs, BufLoc, OverSample, NumInputParams };
    enum Outputs { Out1, NumOutputParams };
  
    float *osBuffer;
    float *os_buf_loc;
    float *buf_loc_interp;
  
    int m_oversamplingIndex{0};
  };

} // namespace OscOS