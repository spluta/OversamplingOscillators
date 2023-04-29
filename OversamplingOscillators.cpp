// OversamplingOscillators.cpp
// Sam Pluta
// uses the VariableOversampling class from Jatin Chowdhurry's ChowDSP library in the ported plugins library by Mads Kjeldgaard
// edit the c_cpp_properties.json and CMakeLists.txt files to correctly point to these files on your system

#include "OversamplingOscillators.hpp"
#include "SC_PlugIn.hpp"
#include "SC_PlugIn.h"

static InterfaceTable *ft;

namespace SawOS
{

  SawOSNext::SawOSNext()
  {
  }

  SawOSNext::~SawOSNext() {}

  // float SawBLNext::next(float freq, float *phase, int *counter, float p0n, float freqMul)
  void SawOSNext::next(float freq, float *m_phase, float m_freqMul, float *osBuffer, int overSamplingRatio)
  {

    float phaseDiff = (phaseIn - m_lastPhase);
    //m_phase = sc_wrap(m_phase, -1.f, 1.f);
    m_lastPhase = phaseIn;

    for (int k = 0; k < overSamplingRatio; k++)
    {
      m_phase += (phaseDiff/overSamplingRatio)
      *m_phase += freq * m_freqMul / (float)overSamplingRatio;
      if (*m_phase >= 1.f)
        *m_phase -= 2.f;
      else if (*m_phase <= -1.f)
        *m_phase += 2.f;
      osBuffer[k] = z;
    }
  }

  SawOS::SawOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    mCalcFunc = make_calc_function<SawOS, &SawOS::next_aa>();
  }
  SawOS::~SawOS() {}

  void SawOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    float *outbuf = out(Out1);

    float *osBuffer = oversample.getOSBuffer();
    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      saw.next(freq[i], &m_phase, m_freqMul, osBuffer, oversample.getOversamplingRatio());
      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
  }
}

namespace SinOscOS
{
  SinTable::SinTable() 
  {
    for (int i = 0; i < 4096; ++i)
      table[i] = sin((float)i/2048.f*pi)*(-1);
    table[4096] = 0.f;
  }

  float SinTable::lookup(float index)
  {
    int low = (int)floor(index);
    float dec = index-(float)low;
    int high = (int)ceil(index);
    float val = table[low]*(1.f-dec)+(table[high]*dec);
    //Print("%f %f \n", index, val);
    return val;
  }
  SinTable::~SinTable() {}

  SinOSNext::SinOSNext() {}

  SinOSNext::SinOSNext(float startPhase) {
    m_lastPhase = startPhase;
    m_phase = startPhase;
  }

  void SinOSNext::next(float freq, float phaseIn, float m_freqMul, float *osBuffer, int overSamplingRatio)
  {
    m_phase += (phaseIn - m_lastPhase);
    m_phase = sc_wrap(m_phase, -1.f, 1.f);
    m_lastPhase = phaseIn;

    saw.next(freq, &m_phase, m_freqMul, osBuffer, overSamplingRatio);
    for (int i2 = 0; i2 < overSamplingRatio; i2++)
      osBuffer[i2] = sinTable.lookup((osBuffer[i2]+1.f)*2048.f);
  }

  SinOSNext::~SinOSNext() {}

  SinOscOS::SinOscOS()
  {
    sample_rate = sampleRate();
    mCalcFunc = make_calc_function<SinOscOS, &SinOscOS::next_aa>();
    next_aa(1);
  }
  SinOscOS::~SinOscOS() {}


  void SinOscOS::next_aa(int nSamples)
  {
    const float *freq = in(Freq);
    const float *phase = in(Phase);
    float *outbuf = out(Out1);

    osBuffer = oversample.getOSBuffer();
    oversample.reset(sample_rate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      float phaseIn = phase[i];

      sine.next(freq[i], phaseIn, m_freqMul, osBuffer, oversample.getOversamplingRatio());
      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
  }

  // void SinOscOS::next_k(int nSamples)
  // {
  //   SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);

  //   float *outbuf = out(Out1);

  //   int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
  //   if (osIndexIn != m_oversamplingIndex)
  //   {
  //     m_oversamplingIndex = osIndexIn;
  //     oversample.setOversamplingIndex(m_oversamplingIndex);
  //   }

  //   float *osBuffer = oversample.getOSBuffer();
  //   for (int i = 0; i < nSamples; ++i)
  //   {
  //     float out;
  //     saw.next(slopedFreq.consume(), &m_phase, &m_freqMul, osBuffer, oversample.getOversamplingRatio());
  //     // make the saw a sine
  //     for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
  //       //osBuffer[i2] = sin((osBuffer[i2]+1.f)*pi);
  //       osBuffer[i2] = sinTable.lookup((osBuffer[i2]+1.f)*2048.f);

  //     if (m_oversamplingIndex != 0)
  //       out = oversample.downsample();
  //     else
  //       out = osBuffer[0];
  //     outbuf[i] = out;
  //   }
  //   m_freq_past = slopedFreq.value;
  // }

  // void SinOscOS::next(int nSamples, AudioSignal<float> freq, float *outbuf)
  // {
  //   float out;
  //   int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);

  //   if (osIndexIn != m_oversamplingIndex)
  //   {
  //     m_oversamplingIndex = osIndexIn;
  //     oversample.setOversamplingIndex(m_oversamplingIndex);
  //   }
  //   float *osBuffer = oversample.getOSBuffer();
  //   for (int i = 0; i < nSamples; ++i)
  //   {
      
  //     saw.next(freq.consume(), &m_phase, &m_freqMul, osBuffer, oversample.getOversamplingRatio());
  //     // make the saw a sin
  //     for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
  //       //osBuffer[i2] = sin((osBuffer[i2]+1.f)*pi);
  //       osBuffer[i2] = sinTable.lookup((osBuffer[i2]+1.f)*2048.f);

  //     if (m_oversamplingIndex != 0)
  //       out = oversample.downsample();
  //     else
  //       out = osBuffer[0];
  //     outbuf[i] = out;
  //   }
  //   //return out;
  // }

  // void SinOscOS::next_a(int nSamples)
  // {

  //   AudioSignal<float> freq =  makeSignal(Freq);
  //   float *outbuf = out(Out1);

  //   next(nSamples, freq, outbuf);
    
  // }

  // void SinOscOS::next_k(int nSamples)
  // {
  //   SlopeSignal<float> freq = makeSlope(in0(Freq), m_freq_past);

  //   //AudioSignal<float> freq =  makeSignal(Freq);
  //   float *outbuf = out(Out1);

  //   next(nSamples, freq, outbuf);
  // }
}

namespace TriOS
{
  TriOS::TriOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<TriOS, &TriOS::next_a>();
    else
      mCalcFunc = make_calc_function<TriOS, &TriOS::next_k>();
    next_a(1);
  }
  TriOS::~TriOS() {}

  void TriOS::next_a(int nSamples)
  {

    const float *freq = in(Freq);
    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);

    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }
    float *osBuffer = oversample.getOSBuffer();
    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      saw.next(freq[i], &m_phase, m_freqMul, osBuffer, oversample.getOversamplingRatio());
      // make the saw a triangle
      for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
        osBuffer[i2] = abs(osBuffer[i2]) * 2 - 1;

      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
  }

  void TriOS::next_k(int nSamples)
  {
    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);

    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    float *osBuffer = oversample.getOSBuffer();
    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      saw.next(slopedFreq.consume(), &m_phase, m_freqMul, osBuffer, oversample.getOversamplingRatio());
      // make the saw a triangle
      for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
        osBuffer[i2] = abs(osBuffer[i2]) * 2 - 1;
      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
    m_freq_past = slopedFreq.value;
  }
}

namespace VarSawOS
{

  VarSawOS::VarSawOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    if (inRate(0) == 2)
      if (inRate(2) == 2)
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_aa>();
      else
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_ak>();
    else if (inRate(2) == 2)
      mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_ka>();
    else
      mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_kk>();
    next_aa(1);
  }

  VarSawOS::~VarSawOS() {}

  float VarSawOS::next(float freq, float width)
  {
    float out;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);
    float *osBuffer = oversample.getOSBuffer();

    saw.next(freq, &m_phase, m_freqMul, osBuffer, oversample.getOversamplingRatio());
    for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
    {
      float temp = osBuffer[i2] / 2 + 0.5;
      osBuffer[i2] = temp < width ? temp * invwidth : (1.f - temp) * inv1width;
    }
    if (m_oversamplingIndex != 0)
      out = oversample.downsample() - 1.f;
    else
      out = osBuffer[0];
    return out;
  }

  void VarSawOS::next_ak(int nSamples)
  {

    const float *freq = in(Freq);
    float *outbuf = out(Out1);
    float width = in0(Width);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);

    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], width);
    }
  }

  void VarSawOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], inWidth[i]);
    }
  }

  void VarSawOS::next_kk(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float width = in0(Width);

    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(slopedFreq.consume(), width);
    }
    m_freq_past = slopedFreq.value;
  }

  void VarSawOS::next_ka(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(slopedFreq.consume(), inWidth[i]);
    }
    m_freq_past = slopedFreq.value;
  }
}

namespace SquareOS
{

  SquareOS::SquareOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    if (inRate(0) == 2)
      if (inRate(2) == 2)
        mCalcFunc = make_calc_function<SquareOS, &SquareOS::next_aa>();
      else
        mCalcFunc = make_calc_function<SquareOS, &SquareOS::next_ak>();
    else if (inRate(2) == 2)
      mCalcFunc = make_calc_function<SquareOS, &SquareOS::next_ka>();
    else
      mCalcFunc = make_calc_function<SquareOS, &SquareOS::next_kk>();
    next_aa(1);
  }

  SquareOS::~SquareOS() {}

  float SquareOS::next(float freq, float width)
  {
    float out;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);
    float *osBuffer = oversample.getOSBuffer();

    saw.next(freq, &m_phase, m_freqMul, osBuffer, oversample.getOversamplingRatio());
    for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
    {
      float temp = osBuffer[i2] / 2 + 0.5;
      osBuffer[i2] = temp < width ? -1.f : 1.f;
    }
    if (m_oversamplingIndex != 0)
      out = oversample.downsample();
    else
      out = osBuffer[0];
    return out;
  }

  void SquareOS::next_ak(int nSamples)
  {

    const float *freq = in(Freq);
    float *outbuf = out(Out1);
    float width = in0(Width);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);

    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], width);
    }
  }

  void SquareOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], inWidth[i]);
    }
  }

  void SquareOS::next_kk(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float width = in0(Width);

    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(slopedFreq.consume(), width);
    }
    m_freq_past = slopedFreq.value;
  }

  void SquareOS::next_ka(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(slopedFreq.consume(), inWidth[i]);
    }
    m_freq_past = slopedFreq.value;
  }
}

namespace SawBL
{
  SawBLNext::SawBLNext()
  {
  }

  SawBLNext::~SawBLNext() {}

  Diff::Diff()
  {
    lastSample = 0.f;
  }

  Diff::~Diff() {}

  float Diff::diff(float sample, float p0n)
  {
    float diff = (sample - lastSample) / (2.0 / p0n);
    lastSample = sample;
    return diff;
  }

  // this code is adapted from Julius Smith's implementation in Faust

  float SawBLNext::next(float freq, float *phase, int *counter, float p0n, float freqMul)
  {
    *phase += freq * freqMul;
    if (*phase >= 1.f)
      *phase -= 2.f;
    else if (*phase <= -1.f)
      *phase += 2.f;

    float poly4 = *phase * *phase * (*phase * *phase - 2.0);
    poly4 = diff4_1.diff(poly4, p0n);
    poly4 = diff4_2.diff(poly4, p0n);
    poly4 = diff4_3.diff(poly4, p0n) / 24.f;

    float poly2 = *phase * *phase;
    poly2 = diff2_1.diff(poly2, p0n) / 2.f;

    float cross = 1.f;

    if (freq < 4000)
      cross = 0.f;
    else if (freq > 8000)
      cross = 1.f;
    else
      cross = (freq - 4000) / 4000;

    float cross2 = 1.f;

    if (*counter < 2)
      poly2 = 0.f;
    if (*counter < 4)
    {
      *counter = *counter + 1;
      poly4 = 0.f;
    }

    float out = (poly2 * (1 - cross)) + (poly4 * cross);

    if (freq < 1.f)
      out = *phase;

    return out;
  }

  SawBL::SawBL()
  {

    m_counter = 0;

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<SawBL, &SawBL::next_a>();
    else
      mCalcFunc = make_calc_function<SawBL, &SawBL::next_k>();
    next_a(1);
  }

  SawBL::~SawBL() {}

  void SawBL::next_a(int nSamples)
  {

    const float *freqIn = in(Freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float freq = abs(freqIn[i]);
      float p0n = sampleRate() / freq;
      float out = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);

      outbuf[i] = out;
    }
  }

  void SawBL::next_k(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(abs(in0(Freq)), m_freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float freq = slopedFreq.consume();

      float p0n = sampleRate() / freq;
      float out = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);

      outbuf[i] = out;
    }
    m_freq = slopedFreq.value;
  }
}

namespace SawH {

  SawH::SawH()
  {
    m_counter = 0;

    oversample.reset(sampleRate());
    oversample.setOversamplingIndex(2);

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<SawH, &SawH::next_a>();
    else
      mCalcFunc = make_calc_function<SawH, &SawH::next_k>();
    next_a(1);
  }

  SawH::~SawH() {}

  void SawH::next_a(int nSamples)
  {
    const float *freqIn = in(Freq);
    float *outbuf = out(Out1);

    oversample.setOversamplingIndex(2);

    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      float freq = abs(freqIn[i])/(float)oversample.getOversamplingRatio();
      float p0n = sampleRate() / (freq);
      for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
        osBuffer[i2] = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);

      outbuf[i] = oversample.downsample();
    }
  }

  void SawH::next_k(int nSamples)
  {
    SlopeSignal<float> slopedFreq = makeSlope(abs(in0(Freq)), m_freq);
    float *outbuf = out(Out1);

    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      float freq = slopedFreq.consume()/(float)oversample.getOversamplingRatio();
      float p0n = sampleRate() / (freq);
      for (int i2 = 0; i2 < oversample.getOversamplingRatio(); i2++)
        osBuffer[i2] = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);
      outbuf[i] = oversample.downsample();
    }
  }
}

namespace SquareBL
{

  SquareNext::SquareNext()
  {
  }

  SquareNext::~SquareNext() {}

  float SquareNext::next(float freq, float duty)
  {
    freq = abs(freq);
    float p0n = m_sampleRate / freq;
    float ddel = duty * p0n;
    float del = sc_max(0, sc_min(ddel, (float)delMax));

    float sawval = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);

    delArray[delArrayCounter] = sawval;

    // diffdel(x,del) = x-x@int(del)*(1-ma.frac(del))-x@(int(del)+1)*ma.frac(del);

    float dec = del - (long)del;
    float delaySig = (delArray[sc_wrap(delArrayCounter - (int)del, 0, delMax)] * (1 - dec)) +
                     (delArray[sc_wrap(delArrayCounter - ((int)del + 1), 0, delMax)] * dec);

    float out = sawval - delaySig;

    delArrayCounter++;
    if (delArrayCounter >= (delMax + 1))
      delArrayCounter = 0;

    float lf_square = sawval + (0.5f - duty) * 2;
    if (lf_square >= 0)
      lf_square = 1.f;
    else
      lf_square = -1.f;

    if (freq < m_fmin)
      out = lf_square - ((0.5f - duty) * 2.f);

    return out; // returns the signal before offset correction because TriBL needs the offset signal
  }

  void SquareNext::setRatePhase(float rateIn, float phaseIn)
  {
    m_sampleRate = rateIn;
    m_phase = phaseIn;
    m_freqMul = 2.0f / m_sampleRate;
    m_fmin = m_sampleRate / float(2 * delMax);
  }

  SquareBL::SquareBL()
  {
    square.setRatePhase(sampleRate(), m_phase);

    if (inRate(0) == 2)
      if (inRate(2) == 2)
        mCalcFunc = make_calc_function<SquareBL, &SquareBL::next_aa>();
      else
        mCalcFunc = make_calc_function<SquareBL, &SquareBL::next_ak>();
    else if (inRate(2) == 2)
      mCalcFunc = make_calc_function<SquareBL, &SquareBL::next_ka>();
    else
      mCalcFunc = make_calc_function<SquareBL, &SquareBL::next_kk>();
    next_aa(1);
  }

  SquareBL::~SquareBL() {}

  void SquareBL::next_aa(int nSamples)
  {
    const float *freqIn = in(Freq);
    const float *dutyIn = in(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float duty = dutyIn[i];
      outbuf[i] = square.next(freqIn[i], duty) + ((0.5f - duty) * 2.f);
    }
  }
  void SquareBL::next_ak(int nSamples)
  {
    const float *freqIn = in(Freq);
    const float duty = in0(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = square.next(freqIn[i], duty) + ((0.5f - duty) * 2.f);
    }
  }
  void SquareBL::next_ka(int nSamples)
  {

    const float *dutyIn = in(Width);
    float *outbuf = out(Out1);

    SlopeSignal<float> freq = makeSlope(in0(Freq), m_freq);

    for (int i = 0; i < nSamples; ++i)
    {
      float duty = dutyIn[i];
      outbuf[i] = square.next(freq.consume(), duty) + ((0.5f - duty) * 2.f);
    }
    m_freq = freq.value;
  }
  void SquareBL::next_kk(int nSamples)
  {

    SlopeSignal<float> freq = makeSlope(in0(Freq), m_freq);

    const float duty = in0(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = square.next(freq.consume(), duty) + ((0.5f - duty) * 2.f);
    }
    m_freq = freq.value;
  }
}

namespace TriBL
{

  TriBL::TriBL()
  {
    square.setRatePhase(sampleRate(), m_phase);

    if (inRate(0) == 2)
      if (inRate(2) == 2)
        mCalcFunc = make_calc_function<TriBL, &TriBL::next_aa>();
      else
        mCalcFunc = make_calc_function<TriBL, &TriBL::next_ak>();
    else if (inRate(2) == 2)
      mCalcFunc = make_calc_function<TriBL, &TriBL::next_ka>();
    else
      mCalcFunc = make_calc_function<TriBL, &TriBL::next_kk>();
    next_aa(1);
  }

  TriBL::~TriBL() {}

  float TriBL::next(float freq, float duty)
  {
    freq = abs(freq);
    float gain = 4 * freq / sampleRate();
    duty = sc_clip(duty, 0.f, 1.f);
    float out = (square.next(freq, duty) + m_delay1 * 0.999);
    m_delay1 = out;
    float gain2 = pow(abs(0.5f - duty) * 2, 5) * 4 + 1;
    return out * gain * gain2;
  }

  void TriBL::next_aa(int nSamples)
  {
    const float *freqIn = in(Freq);
    const float *dutyIn = in(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freqIn[i], dutyIn[i]);
    }
  }

  void TriBL::next_ak(int nSamples)
  {
    const float *freqIn = in(Freq);
    const float duty = in0(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      // float freq = freqIn[i];
      // float gain = 4*freq/sampleRate();
      // float out = square.next(freq, duty)+m_delay1*0.999;
      // outbuf[i] = out*gain;
      // m_delay1 = out;
      outbuf[i] = next(freqIn[i], duty);
    }
  }
  void TriBL::next_ka(int nSamples)
  {

    const float *dutyIn = in(Width);
    float *outbuf = out(Out1);
    SlopeSignal<float> freq = makeSlope(in0(Freq), m_freq);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq.consume(), dutyIn[i]);
    }
    m_freq = freq.value;
  }
  void TriBL::next_kk(int nSamples)
  {
    const float duty = in0(Width);

    SlopeSignal<float> freq = makeSlope(in0(Freq), m_freq);

    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq.consume(), duty);
    }
    m_freq = freq.value;
  }
}

namespace ImpulseBL
{

  Impulse::Impulse() {}
  Impulse::~Impulse() {}

  float Impulse::next()
  {
    impArray[0] = 1;
    float output = (1.f - impArray[1]);
    impArray[1] = impArray[0];
    return output;
  }

  ImpulseBL::ImpulseBL()
  {

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<ImpulseBL, &ImpulseBL::next_a>();
    else
      mCalcFunc = make_calc_function<ImpulseBL, &ImpulseBL::next_k>();

    next_a(1);
  }

  ImpulseBL::~ImpulseBL() {}

  void ImpulseBL::next_a(int nSamples)
  {
    const float *freqIn = in(Freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float freq = abs(freqIn[i]);
      float p0n = sampleRate() / freq;
      float z = saw.next(freqIn[i], &m_phase, &m_counter, p0n, freqMul);
      outbuf[i] = impulse.next() + (m_delay1 - z);
      m_delay1 = z;
    }
  }

  void ImpulseBL::next_k(int nSamples)
  {
    const float freq = abs(in0(Freq));
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float p0n = sampleRate() / freq;
      float z = saw.next(freq, &m_phase, &m_counter, p0n, freqMul);
      outbuf[i] = impulse.next() + (m_delay1 - z);
      m_delay1 = z;
    }
  }

}

PluginLoad(SawOSUGens)
{
  ft = inTable;
  registerUnit<SawOS::SawOS>(ft, "SawOS", false);
  registerUnit<SinOscOS::SinOscOS>(ft, "SinOscOS", false);
  registerUnit<TriOS::TriOS>(ft, "TriOS", false);
  registerUnit<VarSawOS::VarSawOS>(ft, "VarSawOS", false);
  registerUnit<SquareOS::SquareOS>(ft, "SquareOS", false);
  
  registerUnit<SawBL::SawBL>(ft, "SawBL", false);
  registerUnit<SquareBL::SquareBL>(ft, "SquareBL", false);
  registerUnit<TriBL::TriBL>(ft, "TriBL", false);
  registerUnit<ImpulseBL::ImpulseBL>(ft, "ImpulseBL", false);

  registerUnit<SawH::SawH>(ft, "SawH", false);
}
