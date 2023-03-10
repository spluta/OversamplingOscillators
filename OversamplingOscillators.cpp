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
      mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_a>();
    else
      mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_k>();
    next_a(1);
  }

  VarSawOS::~VarSawOS() {}

  void VarSawOS::next_a(int nSamples)
  {

    const float *freq = in(Freq);

    float *outbuf = out(Out1);
    float nextDuty = sc_clip(in0(Duty), 0.001, 0.999);
    float duty = m_duty;
    float invduty = 2.f / duty;
    float inv1duty = 2.f / (1 - duty);

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
          duty = m_duty = sc_clip(nextDuty, 0.001, 0.999);
          invduty = m_invduty = 2.f / duty;
          inv1duty = m_inv1duty = 2.f / (1.f - duty);
        }

        float z = m_phase < duty ? m_phase * invduty : (1.f - m_phase) * inv1duty;
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

  void VarSawOS::next_k(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);

    float *outbuf = out(Out1);
    float nextDuty = sc_clip(in0(Duty), 0.001, 0.999);
    float duty = m_duty;
    float invduty = 2.f / duty;
    float inv1duty = 2.f / (1 - duty);

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
          duty = m_duty = sc_clip(nextDuty, 0.001, 0.999);
          invduty = m_invduty = 2.f / duty;
          inv1duty = m_inv1duty = 2.f / (1.f - duty);
        }

        float z = m_phase < duty ? m_phase * invduty : (1.f - m_phase) * inv1duty;
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

      m_counter = m_counter + 1;

      m_freq = abs(freq[i]);

      float lastPhase = m_phase;
      float lastPoly = m_poly;

      float p0n = sampleRate() / lastFreq;

      m_phase += lastFreq * m_freqMul;
      if (m_phase >= 1.f)
        m_phase -= 2.f;
      else if (m_phase <= -1.f)
        m_phase += 2.f;

      m_poly = m_phase * m_phase * (m_phase * m_phase - 2.0);

      m_poly = diff1.diff(m_poly, p0n);
      m_poly = diff2.diff(m_poly, p0n);
      m_poly = diff3.diff(m_poly, p0n);

      m_poly = m_poly / 24.f;

      float cross = 1.f;

      if (lastFreq < 100)
      {
        cross = 0.f;
      }
      else if (lastFreq > 500)
      {
        cross = 1.f;
      }
      else
      {
        cross = (lastFreq - 100) / 400;
      }

      if (m_counter < 4)
        lastPoly *= 0;

      outbuf[i] = lastPhase * (1 - cross) + lastPoly * cross;
    }
  }

  void SawPn::next_k(int nSamples)
  {

    SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float lastFreq = m_freq;
      m_counter = m_counter + 1;
      m_freq = abs(slopedFreq.consume());

      float lastPhase = m_phase;
      float lastPoly = m_poly;
      float p0n = sampleRate() / lastFreq;

      m_phase += lastFreq * m_freqMul;
      if (m_phase >= 1.f)
        m_phase -= 2.f;
      else if (m_phase <= -1.f)
        m_phase += 2.f;

      m_poly = m_phase * m_phase * (m_phase * m_phase - 2.0);

      m_poly = diff1.diff(m_poly, p0n);
      m_poly = diff2.diff(m_poly, p0n);
      m_poly = diff3.diff(m_poly, p0n);

      m_poly = m_poly / 24.f;
      float cross = 1.f;

      if (lastFreq < 100)
      {
        cross = 0.f;
      }
      else if (lastFreq > 500)
      {
        cross = 1.f;
      }
      else
      {
        cross = (lastFreq - 100) / 400;
      }

      if (m_counter < 4)
        lastPoly *= 0;

      outbuf[i] = lastPhase * (1 - cross) + lastPoly * cross;
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
