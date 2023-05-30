// OversamplingOscillators.hpp

#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"

namespace SawOS {

class SawOSNext {
  public:
    SawOSNext();
    ~SawOSNext();
    float m_lastPhase {0.f};
    float m_phase {0.f};
    float next(float freq, float phaseIn, float m_freqMul);
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
  void next_aa(int nSamples);
  SawOSNext saw;

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq_past{0.f};
  float i_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};
  int m_oversamplingIndex{0};
};

} // namespace SawOS

namespace SinOscOS {
  class SinTable {
    public:
      float lookup(float x);
      SinTable();
      ~SinTable();
    private:
      float table [4097];
  };

  class SinOSNext {
    public:
      SinOSNext();
      SinOSNext(float startPhase);
      ~SinOSNext();
      float next(float freq, float phaseIn, float m_freqMul);
      float m_val{0};
      
    private:
      SawOS::SawOSNext saw;
      SinTable sinTable;

      float m_lastPhase;
  };

class SinOscOS : public SCUnit {
  public:
    SinOscOS();
    ~SinOscOS();

    SinOSNext sine;
    
    VariableOversampling<> oversample;

  private:
    // Calc function
    void next_aa(int nSamples);

    enum InputParams { Freq, Phase, OverSample, NumInputParams };
    enum Outputs { Out1, NumOutputParams };

    float sample_rate;
    float m_freqMul{2.0f/(float)sampleRate()};
    int m_oversamplingIndex{(int)in0(OverSample)};
  };
} //namespace SinOscOS

namespace PMOscOS {

class PMOscOS : public SCUnit {
  public:
    PMOscOS();
    ~PMOscOS();

    SinOscOS::SinOSNext sine0;
    SinOscOS::SinOSNext sine1;

    VariableOversampling<> oversample;

    float *osBuffer;

  private:
    // Calc function
    void next_aa(int nSamples);

//carfreq, modfreq, pmindex: 0.0, modphase
    enum InputParams { CarFreq, ModFreq, PMMul, PMModPhase, OverSample, NumInputParams };
    enum Outputs { Out1, NumOutputParams };

    float sample_rate{(float)sampleRate()};
    float m_modphase{in0(PMModPhase)};
    float m_freqMul{2.0f/(float)sampleRate()};
    int m_oversamplingIndex{(int)in0(OverSample)};
  };
} //namespace PMOscOS

namespace FM7OS {

class FM7OS : public SCUnit {
  public:
    FM7OS();
    ~FM7OS();

    SinOscOS::SinOSNext sines[6];

    VariableOversampling<> oversamples[6];
    float *osBuffers[6];

  private:
    // Calc function
    void next_aa(int nSamples);

    enum InputParams { ctl0, ctl1, ctl2, ctl3, ctl4, ctl5, ctl6, ctl7, ctl8, ctl9, ctl10, ctl11, ctl12, ctl13, ctl14, ctl15, ctl16, ctl17, 
    modNum0, modNum1, modNum2, modNum3, modNum4, modNum5, modNum6, modNum7, modNum8, modNum9, modNum10, modNum11, modNum12, modNum13, modNum14, modNum15, modNum16, modNum17, modNum18, modNum19, modNum20, modNum21, modNum22, modNum23, modNum24, modNum25, modNum26, modNum27, modNum28, modNum29, modNum30, modNum31, modNum32, modNum33, modNum34, modNum35,
    OverSample, NumInputParams };
    enum Outputs { Out1, Out2, Out3, Out4, Out5, Out6, NumOutputParams };

    float sample_rate{(float)sampleRate()};
    
    float m_phases[6] = {in0(ctl1),in0(ctl4),in0(ctl7),in0(ctl0),in0(ctl3),in0(ctl6)};

    float m_freqMul{2.0f/(float)sampleRate()};
    int m_oversamplingIndex{sc_clip((int)in0(OverSample), 0, 4)};
    int m_oversampleRatio{(int)pow(2, m_oversamplingIndex)};
  };
} //namespace FM7OS

namespace PM7OS {
class PM7OS : public SCUnit {
  public:
    PM7OS();
    ~PM7OS();

    SinOscOS::SinOSNext sines[6];

    VariableOversampling<> oversamples[6];
    float *osBuffers[6];

  private:
    // Calc function
    void next_aa(int nSamples);

    enum InputParams { ctl0, ctl1, ctl2, ctl3, ctl4, ctl5, ctl6, ctl7, ctl8, ctl9, ctl10, ctl11, ctl12, ctl13, ctl14, ctl15, ctl16, ctl17, 
    modNum0, modNum1, modNum2, modNum3, modNum4, modNum5, modNum6, modNum7, modNum8, modNum9, modNum10, modNum11, modNum12, modNum13, modNum14, modNum15, modNum16, modNum17, modNum18, modNum19, modNum20, modNum21, modNum22, modNum23, modNum24, modNum25, modNum26, modNum27, modNum28, modNum29, modNum30, modNum31, modNum32, modNum33, modNum34, modNum35,
    OverSample, NumInputParams };
    enum Outputs { Out1, Out2, Out3, Out4, Out5, Out6, NumOutputParams };

    float sample_rate{(float)sampleRate()};
    
    float m_phases[6] = {in0(ctl1),in0(ctl4),in0(ctl7),in0(ctl0),in0(ctl3),in0(ctl6)};

    float m_freqMul{2.0f/(float)sampleRate()};
    int m_oversamplingIndex{sc_clip((int)in0(OverSample), 0, 4)};
    int m_oversampleRatio{(int)pow(2, m_oversamplingIndex)};
  };
} //namespace PM7OS

namespace TriOS {
class TriOS : public SCUnit {
public:
  TriOS();

  // Destructor
  ~TriOS();
  VariableOversampling<> oversample;

private:
  // Calc function
  void next_aa(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq_past{0.f};
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
  float next(float freq, float phase, float width);
  void next_aa(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, Width, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

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
  float next(float freq, float phase, float width);
  void next_aa(int nSamples);

  SawOS::SawOSNext saw;

  enum InputParams { Freq, Phase, Width, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq_past{in0(Freq)};
  float m_phase{(float)in0(Phase)};
  float m_width{in0(Width)};
  float m_freqMul{2.0f/(float)sampleRate()};
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

    Diff diff4_1;
    Diff diff4_2;
    Diff diff4_3;

    Diff diff2_1;
};

class SawBL : public SCUnit {
public:
  SawBL();
  ~SawBL();

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);

  enum InputParams { Freq, Phase, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};

  int m_counter{0};

  float diffArray [3] = {m_phase, m_phase, m_phase};

  SawBLNext saw;
};

} // namespace SawBL


namespace SawH {

class SawH : public SCUnit {
public:
  SawH();
  ~SawH();

private:
  // Calc function
  void next_a(int nSamples);
  void next_k(int nSamples);

  enum InputParams { Freq, Phase, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

  float m_freq{abs(in0(Freq))};
  float m_phase{in0(Phase)};
  float m_freqMul{2.0f/(float)sampleRate()};

  int m_counter{0};

  float diffArray [3] = {m_phase, m_phase, m_phase};

  SawBL::SawBLNext saw;
  VariableOversampling<> oversample;
};

} // namespace SawH

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