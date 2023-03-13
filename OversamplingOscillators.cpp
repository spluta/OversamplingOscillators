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

  SawOS::SawOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<SawOS, &SawOS::next_a>();
    else
      mCalcFunc = make_calc_function<SawOS, &SawOS::next_k>();
    next_a(1);
  }

  SawOS::~SawOS() {}

  void SawOS::next_a(int nSamples)
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
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        float z = m_phase;
        m_phase += freq[i] * m_freqMul / (float)oversample.getOversamplingRatio();
        if (m_phase >= 1.f)
          m_phase -= 2.f;
        else if (m_phase <= -1.f)
          m_phase += 2.f;
        osBuffer[k] = z;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
    }
  }

  void SawOS::next_k(int nSamples)
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
      const float freq = slopedFreq.consume();
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        float z = m_phase;
        m_phase += freq * m_freqMul / (float)oversample.getOversamplingRatio();
        if (m_phase >= 1.f)
          m_phase -= 2.f;
        else if (m_phase <= -1.f)
          m_phase += 2.f;
        osBuffer[k] = z;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
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
      if (inRate(2)==2)
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_aa>();
      else
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_ak>();
    else
      if (inRate(2)==2)
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_ka>();
      else
        mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_kk>();
    next_aa(1);
  }

  VarSawOS::~VarSawOS() {}

  void VarSawOS::next_ak(int nSamples)
  {

    const float *freq = in(Freq);

    float *outbuf = out(Out1);
    float nextWidth = sc_clip(in0(Width), 0.001, 0.999);
    float width = m_width;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    // float y;
    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      float phaseInc = freq[i] * m_freqMul / (float)oversample.getOversamplingRatio();

      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        if (m_phase >= 1.f)
        {
          m_phase -= 1.f;
          width = m_width = sc_clip(nextWidth, 0.001, 0.999);
          invwidth = m_invwidth = 2.f / width;
          inv1width = m_inv1width = 2.f / (1.f - width);
        }

        float z = m_phase < width ? m_phase * invwidth : (1.f - m_phase) * inv1width;
        osBuffer[k] = z;
        m_phase += phaseInc;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y - 1.f;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
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

    // float y;
    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      float phaseInc = freq[i] * m_freqMul / (float)oversample.getOversamplingRatio();

      float nextWidth = sc_clip(inWidth[i], 0.001, 0.999);
      float width = m_width;
      float invwidth = 2.f / width;
      float inv1width = 2.f / (1 - width);

      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        if (m_phase >= 1.f)
        {
          m_phase -= 1.f;
          width = m_width = sc_clip(nextWidth, 0.001, 0.999);
          invwidth = m_invwidth = 2.f / width;
          inv1width = m_inv1width = 2.f / (1.f - width);
        }

        float z = m_phase < width ? m_phase * invwidth : (1.f - m_phase) * inv1width;
        osBuffer[k] = z;
        m_phase += phaseInc;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y - 1.f;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
    }
  }

  void VarSawOS::next_kk(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float *inWidth = in(Width);

    float *outbuf = out(Out1);
    float nextWidth = sc_clip(in0(Width), 0.001, 0.999);
    float width = m_width;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      const float freq = slopedFreq.consume();


      float phaseInc = freq * m_freqMul / oversample.getOversamplingRatio();

      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        if (m_phase >= 1.f)
        {
          m_phase -= 1.f;
          width = m_width = sc_clip(nextWidth, 0.001, 0.999);
          invwidth = m_invwidth = 2.f / width;
          inv1width = m_inv1width = 2.f / (1.f - width);
        }

        float z = m_phase < width ? m_phase * invwidth : (1.f - m_phase) * inv1width;
        osBuffer[k] = z;
        m_phase += phaseInc;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y - 1.f;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
    }
    m_freq_past = slopedFreq.value;
  }

  void VarSawOS::next_ka(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
    const float *inWidth = in(Width);

    float *outbuf = out(Out1);
    float nextWidth = sc_clip(in0(Width), 0.001, 0.999);
    float width = m_width;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);

    int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
    if (osIndexIn != m_oversamplingIndex)
    {
      m_oversamplingIndex = osIndexIn;
      oversample.setOversamplingIndex(m_oversamplingIndex);
    }

    float *osBuffer = oversample.getOSBuffer();

    for (int i = 0; i < nSamples; ++i)
    {
      const float freq = slopedFreq.consume();

      float nextWidth = sc_clip(inWidth[i], 0.001, 0.999);
      float width = m_width;
      float invwidth = 2.f / width;
      float inv1width = 2.f / (1 - width);

      float phaseInc = freq * m_freqMul / oversample.getOversamplingRatio();

      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        if (m_phase >= 1.f)
        {
          m_phase -= 1.f;
          width = m_width = sc_clip(nextWidth, 0.001, 0.999);
          invwidth = m_invwidth = 2.f / width;
          inv1width = m_inv1width = 2.f / (1.f - width);
        }

        float z = m_phase < width ? m_phase * invwidth : (1.f - m_phase) * inv1width;
        osBuffer[k] = z;
        m_phase += phaseInc;
      }
      if (m_oversamplingIndex != 0)
      {
        float y = oversample.downsample();
        outbuf[i] = y - 1.f;
      }
      else
      {
        outbuf[i] = osBuffer[0];
      }
    }
    m_freq_past = slopedFreq.value;
  }
}

namespace SawPn
{

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

  SawPn::SawPn()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    m_counter = 0;

    if (inRate(0) == 2)
      mCalcFunc = make_calc_function<SawPn, &SawPn::next_a>();
    else
      mCalcFunc = make_calc_function<SawPn, &SawPn::next_k>();
    next_a((int)(m_freq / sampleRate()));
  }

  SawPn::~SawPn() {}

  void SawPn::next_a(int nSamples)
  {

    const float *freq = in(Freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float lastFreq = m_freq;

      m_freq = abs(freq[i]);
      m_freq = sc_max(m_freq, 0.0001);

      float lastPhase = m_phase;
      float lastPoly4 = m_poly4;
      float lastPoly2 = m_poly2;
      float p0n = sampleRate() / lastFreq;

      m_phase += lastFreq * m_freqMul;
      if (m_phase >= 1.f)
        m_phase -= 2.f;
      else if (m_phase <= -1.f)
        m_phase += 2.f;

      m_poly4 = m_phase * m_phase * (m_phase * m_phase - 2.0);
      m_poly4 = diff4_1.diff(m_poly4, p0n);
      m_poly4 = diff4_2.diff(m_poly4, p0n);
      m_poly4 = diff4_3.diff(m_poly4, p0n);
      m_poly4 = m_poly4 / 24.f;

      m_poly2 = m_phase * m_phase;
      m_poly2 = diff2_1.diff(m_poly2, p0n);
      m_poly2 = m_poly2 / 2.f;


      float cross = 1.f;

      if (lastFreq < 5000)
      {
        cross = 0.f;
      }
      else if (lastFreq > 10000)
      {
        cross = 1.f;
      }
      else
      {
        cross = (lastFreq - 5000) / 5000;
      }

      if(m_counter < 2)
        lastPoly2 = 0.f;
      if (m_counter < 4)
      {
        m_counter ++;
        lastPoly4 = 0.f;
      }
      //Print("%f %f \n", lastPoly2, lastPoly4);
      float out = (lastPoly2 * (1 - cross)) + (lastPoly4 * cross);

      // if(out>1.f||out<(-1.f))
      //   Print("%f %f %f \n", lastPoly2, lastPoly4, cross);

      outbuf[i] = out;
    }
  }

  void SawPn::next_k(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float lastFreq = m_freq;
      
      m_freq = abs(slopedFreq.consume());
      m_freq = sc_max(m_freq, 0.0001);

      float lastPhase = m_phase;
      float lastPoly4 = m_poly4;
      float lastPoly2 = m_poly2;
      float p0n = sampleRate() / lastFreq;

      m_phase += lastFreq * m_freqMul;
      if (m_phase >= 1.f)
        m_phase -= 2.f;
      else if (m_phase <= -1.f)
        m_phase += 2.f;

      m_poly4 = m_phase * m_phase * (m_phase * m_phase - 2.0);
      m_poly4 = diff4_1.diff(m_poly4, p0n);
      m_poly4 = diff4_2.diff(m_poly4, p0n);
      m_poly4 = diff4_3.diff(m_poly4, p0n);
      m_poly4 = m_poly4 / 24.f;

      m_poly2 = m_phase * m_phase;
      m_poly2 = diff2_1.diff(m_poly2, p0n);
      m_poly2 = m_poly2 / 2.f;

      float cross = 1.f;

      if (lastFreq < 5000)
      {
        cross = 0.f;
      }
      else if (lastFreq > 10000)
      {
        cross = 1.f;
      }
      else
      {
        cross = (lastFreq - 5000) / 5000;
      }

      if(m_counter < 2)
        lastPoly2 = 0.f;
      if (m_counter < 4)
      {
        m_counter ++;
        lastPoly4 = 0.f;
      }
      outbuf[i] = lastPoly2 * (1 - cross) + lastPoly4 * cross;
    }
  }
}

PluginLoad(SawOSUGens)
{
  ft = inTable;
  registerUnit<SawOS::SawOS>(ft, "SawOS", false);
  registerUnit<VarSawOS::VarSawOS>(ft, "VarSawOS", false);
  registerUnit<SawPn::SawPn>(ft, "SawPn", false);
}
