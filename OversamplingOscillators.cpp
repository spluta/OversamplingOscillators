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
  float SawOSNext::next(float freq, float phaseIn, float m_freqMul)
  {
    float phaseDiff = (phaseIn - m_lastPhase);
    m_lastPhase = phaseIn;

    m_phase += (phaseDiff);
    m_phase += freq * m_freqMul;
    if (m_phase >= 4.f)
      m_phase -= 8.f;
    else if (m_phase <= -4.f)
      m_phase += 8.f;

    float out = sc_wrap(m_phase, -1.f, 1.f);
    return out;
  }

  SawOS::SawOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<SawOS, &SawOS::next_aa>();
    next_aa(1);
  }
  SawOS::~SawOS() {}

  void SawOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *phase = in(Phase);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
        osBuffer[k] = saw.next(freq[i], phase[i], m_freqMul / oversample.getOversamplingRatio());
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
      table[i] = sin((float)i / 2048.f * pi) * (-1);
    table[4096] = 0.f;
  }

  float SinTable::lookup(float index)
  {
    int low = (int)floor(index);
    float dec = index - (float)low;
    int high = (int)ceil(index);
    float val = table[low] * (1.f - dec) + (table[high] * dec);
    // Print("%f %f \n", index, val);
    return val;
  }
  SinTable::~SinTable() {}

  SinOSNext::SinOSNext() {}

  SinOSNext::SinOSNext(float startPhase)
  {
    // m_lastPhase = startPhase;
    // m_phase = startPhase;
  }

  float SinOSNext::next(float freq, float phase, float m_freqMul)
  {
    float out = saw.next(freq, phase + m_phaseOffset, m_freqMul); // saw returns a value between -1 and 1
    m_val = sin((out + 1) * pi);
    return m_val;
    // return sinTable.lookup((out+1.f)*2048.f);
    //  for(int k = 0; k<overSamplingRatio; k++)
    //      osBuffer[k] = saw.next(freq, phase, m_freqMul/overSamplingRatio);
    //  for (int i2 = 0; i2 < overSamplingRatio; i2++)
    //    osBuffer[i2] = sinTable.lookup((osBuffer[i2]+1.f)*2048.f);
  }

  SinOSNext::~SinOSNext() {}

  SinOscOS::SinOscOS()
  {
    sample_rate = sampleRate();

    oversample.reset(sample_rate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<SinOscOS, &SinOscOS::next_aa>();
    next_aa(1);
  }
  SinOscOS::~SinOscOS() {}

  void SinOscOS::next_aa(int nSamples)
  {
    const float *freq = in(Freq);
    const float *phase = in(Phase);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
        osBuffer[k] = sine.next(freq[i], phase[i], m_freqMul / oversample.getOversamplingRatio());

      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
  }
}

namespace PMOscOS
{
  PMOscOS::PMOscOS()
  {
    sample_rate = sampleRate();

    oversample.reset(sample_rate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<PMOscOS, &PMOscOS::next_aa>();
    next_aa(1);
  }
  PMOscOS::~PMOscOS() {}

  void PMOscOS::next_aa(int nSamples)
  {

    // CarFreq, ModFreq, PMMul, PMModPhase
    const float *carFreq = in(CarFreq);
    const float *modFreq = in(ModFreq);
    const float *pm_mul = in(PMMul);
    const float *pm_phase = in(PMModPhase);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      float mod;
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      {
        mod = sine0.next(modFreq[i], pm_phase[i], m_freqMul / oversample.getOversamplingRatio());
        osBuffer[k] = sine1.next(carFreq[i], mod * pm_mul[i], m_freqMul / oversample.getOversamplingRatio());
      }
      if (m_oversamplingIndex != 0)
        out = oversample.downsample();
      else
        out = osBuffer[0];
      outbuf[i] = out;
    }
  }
}

namespace FM7bOS
{
  FM7bOS::FM7bOS()
  {
    sample_rate = (float)sampleRate();
    for (int i = 0; i < 6; i++)
      oversamples[i].reset(sample_rate);

    for (int i = 0; i < 6; i++)
      oversamples[i].setOversamplingIndex(m_oversamplingIndex);

    m_oversampleRatio = oversamples[0].getOversamplingRatio();
    m_freqMul = m_freqMul / (float)m_oversampleRatio;

    for (int k = 0; k < 6; k++)
    {
      osBuffers[k] = oversamples[k].getOSBuffer();
      m_vals[k] = 0;
    }

    mCalcFunc = make_calc_function<FM7bOS, &FM7bOS::next_aa>();
    next_aa(1);
  }
  FM7bOS::~FM7bOS() {}

  void FM7bOS::next_aa(int nSamples)
  {
    const float *freqs[6] = {in(ctl0), in(ctl3), in(ctl6), in(ctl9), in(ctl12), in(ctl15)};
    const float *amps[6] = {in(ctl2), in(ctl5), in(ctl8), in(ctl11), in(ctl14), in(ctl17)};

    const float *mods[6][6] = {{in(modNum0), in(modNum1), in(modNum2), in(modNum3), in(modNum4), in(modNum5)},
                               {in(modNum6), in(modNum7), in(modNum8), in(modNum9), in(modNum10), in(modNum11)},
                               {in(modNum12), in(modNum13), in(modNum14), in(modNum15), in(modNum16), in(modNum17)},
                               {in(modNum18), in(modNum19), in(modNum20), in(modNum21), in(modNum22), in(modNum23)},
                               {in(modNum24), in(modNum25), in(modNum26), in(modNum27), in(modNum28), in(modNum29)},
                               {in(modNum30), in(modNum31), in(modNum32), in(modNum33), in(modNum34), in(modNum35)}};

    float *outs[6] = {out(Out1), out(Out2), out(Out3), out(Out4), out(Out5), out(Out6)};

    float freqs2[6];

    int wavetypes[6];
    for (int m = 0; m < 6; m++)
      wavetypes[m]=(int)mods[m][m][0]; 

    for (int i = 0; i < nSamples; ++i)
    {
      float outSamps[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      for (int sawNum = 0; sawNum < 6; sawNum++)
      {
        freqs2[sawNum] = abs(freqs[sawNum][i]);

        
        float mod = 0;
        for (int m = 0; m < 6; m++)
          if (m != sawNum)
            mod = mod + (m_vals[m] * mods[sawNum][m][i]);
        mod = sc_clip(mod, -1.f, 1.f);

        freqs2[sawNum] = freqs2[sawNum] * powf(2.f, mod * 5.f);
        for (int k = 0; k < m_oversampleRatio; k++)
        {
          float out = saws[sawNum].next(freqs2[sawNum], m_phases[sawNum], m_freqMul) * amps[sawNum][i];
          switch (wavetypes[sawNum])
          {
          case 0:
            out = sin((out + 1) * pi);
            break;
          case 1:
            out = abs(out) * 2.f - 1.f;
            break;
          case 2:
            if (out >= 0)
              out = 1.f;
            else
              out = -1.f;
            break;
          default:
            out = out;
          }
          m_vals[sawNum] = out;
          osBuffers[sawNum][k] = out;
        }
      }

      for (int k = 0; k < 6; k++)
      {
        if (m_oversamplingIndex != 0)
          outSamps[k] = oversamples[k].downsample();
        else
          outSamps[k] = osBuffers[k][0];
        outs[k][i] = outSamps[k];
      }
    }
  }
}


namespace FM7aOS
{
  FM7aOS::FM7aOS()
  {
    sample_rate = (float)sampleRate();
    for (int i = 0; i < 6; i++)
      oversamples[i].reset(sample_rate);

    for (int i = 0; i < 6; i++)
      oversamples[i].setOversamplingIndex(m_oversamplingIndex);

    m_oversampleRatio = oversamples[0].getOversamplingRatio();
    m_freqMul = m_freqMul / (float)m_oversampleRatio;

    for (int k = 0; k < 6; k++)
    {
      osBuffers[k] = oversamples[k].getOSBuffer();
      m_vals[k] = 0;
    }

    mCalcFunc = make_calc_function<FM7aOS, &FM7aOS::next_aa>();
    next_aa(1);
  }
  FM7aOS::~FM7aOS() {}

  void FM7aOS::next_aa(int nSamples)
  {
    const float *freqs[6] = {in(ctl0), in(ctl3), in(ctl6), in(ctl9), in(ctl12), in(ctl15)};
    const float *amps[6] = {in(ctl2), in(ctl5), in(ctl8), in(ctl11), in(ctl14), in(ctl17)};

    const float *mods[6][6] = {{in(modNum0), in(modNum1), in(modNum2), in(modNum3), in(modNum4), in(modNum5)},
                               {in(modNum6), in(modNum7), in(modNum8), in(modNum9), in(modNum10), in(modNum11)},
                               {in(modNum12), in(modNum13), in(modNum14), in(modNum15), in(modNum16), in(modNum17)},
                               {in(modNum18), in(modNum19), in(modNum20), in(modNum21), in(modNum22), in(modNum23)},
                               {in(modNum24), in(modNum25), in(modNum26), in(modNum27), in(modNum28), in(modNum29)},
                               {in(modNum30), in(modNum31), in(modNum32), in(modNum33), in(modNum34), in(modNum35)}};

    float *outs[6] = {out(Out1), out(Out2), out(Out3), out(Out4), out(Out5), out(Out6)};

    float freqs2[6];

    int wavetypes[6];
    for (int m = 0; m < 6; m++)
      wavetypes[m]=(int)mods[m][m][0]; 


    for (int i = 0; i < nSamples; ++i)
    {
      float outSamps[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      for (int sawNum = 0; sawNum < 6; sawNum++)
      {
        freqs2[sawNum] = freqs[sawNum][i];

          float phaseMod = 0.f;
          for (int m = 0; m < 6; m++)
          {
            if (m != sawNum)
              freqs2[sawNum] = freqs2[sawNum] + (m_vals[m] * mods[sawNum][m][i]);
          }
    
          //osBuffers[sineNum][k] = sines[sineNum].next(freqs2[sineNum], m_phases[sineNum] + phaseMod, m_freqMul) * amps[sineNum][i]; // m_phases[sineNum]

        for (int k = 0; k < m_oversampleRatio; k++)
        {
          float out = saws[sawNum].next(freqs2[sawNum], m_phases[sawNum], m_freqMul) * amps[sawNum][i];
          switch (wavetypes[sawNum])
          {
          case 0:
            out = sin((out + 1) * pi);
            break;
          case 1:
            out = abs(out) * 2.f - 1.f;
            break;
          case 2:
            if (out >= 0)
              out = 1.f;
            else
              out = -1.f;
            break;
          default:
            out = out;
          }
          m_vals[sawNum] = out;
          osBuffers[sawNum][k] = out;
        }
      }

      for (int k = 0; k < 6; k++)
      {
        if (m_oversamplingIndex != 0)
          outSamps[k] = oversamples[k].downsample();
        else
          outSamps[k] = osBuffers[k][0];
        outs[k][i] = outSamps[k];
      }
    }
  }
}

namespace FM7OS
{
  FM7OS::FM7OS()
  {
    sample_rate = (float)sampleRate();
    for (int i = 0; i < 6; i++)
      oversamples[i].reset(sample_rate);

    // m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);

    for (int i = 0; i < 6; i++)
      oversamples[i].setOversamplingIndex(m_oversamplingIndex);

    m_oversampleRatio = oversamples[0].getOversamplingRatio();
    m_freqMul = m_freqMul / (float)m_oversampleRatio;

    for (int k = 0; k < 6; k++)
      osBuffers[k] = oversamples[k].getOSBuffer();

    mCalcFunc = make_calc_function<FM7OS, &FM7OS::next_aa>();
    next_aa(1);
  }
  FM7OS::~FM7OS() {}

  void FM7OS::next_aa(int nSamples)
  {
    const float *freqs[6] = {in(ctl0), in(ctl3), in(ctl6), in(ctl9), in(ctl12), in(ctl15)};
    const float *amps[6] = {in(ctl2), in(ctl5), in(ctl8), in(ctl11), in(ctl14), in(ctl17)};

    const float *mods[6][6] = {{in(modNum0), in(modNum1), in(modNum2), in(modNum3), in(modNum4), in(modNum5)},
                               {in(modNum6), in(modNum7), in(modNum8), in(modNum9), in(modNum10), in(modNum11)},
                               {in(modNum12), in(modNum13), in(modNum14), in(modNum15), in(modNum16), in(modNum17)},
                               {in(modNum18), in(modNum19), in(modNum20), in(modNum21), in(modNum22), in(modNum23)},
                               {in(modNum24), in(modNum25), in(modNum26), in(modNum27), in(modNum28), in(modNum29)},
                               {in(modNum30), in(modNum31), in(modNum32), in(modNum33), in(modNum34), in(modNum35)}};

    float *outs[6] = {out(Out1), out(Out2), out(Out3), out(Out4), out(Out5), out(Out6)};

    float freqs2[6];

    for (int i = 0; i < nSamples; ++i)
    {
      float outSamps[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      for (int k = 0; k < m_oversampleRatio; k++)
      {
        for (int sineNum = 0; sineNum < 6; sineNum++)
        {
          freqs2[sineNum] = freqs[sineNum][i];
          float phaseMod = 0.f;
          for (int m = 0; m < 6; m++)
          {
            if (m != sineNum)
              freqs2[sineNum] = freqs2[sineNum] + (sines[m].m_val * mods[sineNum][m][i] * amps[m][i]);
            else
              phaseMod = sines[m].m_val * mods[sineNum][m][i];
          }
          osBuffers[sineNum][k] = sines[sineNum].next(freqs2[sineNum], m_phases[sineNum] + phaseMod, m_freqMul) * amps[sineNum][i]; // m_phases[sineNum]
        }
      }

      for (int k = 0; k < 6; k++)
      {
        if (m_oversamplingIndex != 0)
          outSamps[k] = oversamples[k].downsample();
        else
          outSamps[k] = osBuffers[k][0];
        outs[k][i] = outSamps[k];
      }
    }
  }
}

namespace PM7OS
{
  PM7OS::PM7OS()
  {
    sample_rate = sampleRate();
    for (int i = 0; i < 6; i++)
      oversamples[i].reset(sample_rate);

    for (int i = 0; i < 6; i++)
      oversamples[i].setOversamplingIndex(m_oversamplingIndex);

    m_oversampleRatio = oversamples[0].getOversamplingRatio();
    m_freqMul = m_freqMul / (float)m_oversampleRatio;

    for (int k = 0; k < 6; k++)
      osBuffers[k] = oversamples[k].getOSBuffer();

    mCalcFunc = make_calc_function<PM7OS, &PM7OS::next_aa>();
    next_aa(1);
  }
  PM7OS::~PM7OS() {}

  void PM7OS::next_aa(int nSamples)
  {
    const float *freqs[6] = {in(ctl0), in(ctl3), in(ctl6), in(ctl9), in(ctl12), in(ctl15)};
    const float *amps[6] = {in(ctl2), in(ctl5), in(ctl8), in(ctl11), in(ctl14), in(ctl17)};

    const float *mods[6][6] = {{in(modNum0), in(modNum1), in(modNum2), in(modNum3), in(modNum4), in(modNum5)},
                               {in(modNum6), in(modNum7), in(modNum8), in(modNum9), in(modNum10), in(modNum11)},
                               {in(modNum12), in(modNum13), in(modNum14), in(modNum15), in(modNum16), in(modNum17)},
                               {in(modNum18), in(modNum19), in(modNum20), in(modNum21), in(modNum22), in(modNum23)},
                               {in(modNum24), in(modNum25), in(modNum26), in(modNum27), in(modNum28), in(modNum29)},
                               {in(modNum30), in(modNum31), in(modNum32), in(modNum33), in(modNum34), in(modNum35)}};

    float *outs[6] = {out(Out1), out(Out2), out(Out3), out(Out4), out(Out5), out(Out6)};

    for (int i = 0; i < nSamples; ++i)
    {
      float outSamps[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      for (int k = 0; k < m_oversampleRatio; k++)
      {

        for (int sineNum = 0; sineNum < 6; sineNum++)
        {
          float phaseMod = 0.f;
          for (int m = 0; m < 6; m++)
            phaseMod = phaseMod + (sines[m].m_val * mods[sineNum][m][i] * amps[m][i]);
          osBuffers[sineNum][k] = sines[sineNum].next(freqs[sineNum][i], m_phases[sineNum] + phaseMod, m_freqMul) * amps[sineNum][i];
        }
      }

      for (int k = 0; k < 6; k++)
      {

        if (m_oversamplingIndex != 0)
          outSamps[k] = oversamples[k].downsample();
        else
          outSamps[k] = osBuffers[k][0];

        outs[k][i] = outSamps[k];
      }
    }
  }
}

namespace TriOS
{
  TriOS::TriOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<TriOS, &TriOS::next_aa>();
  }
  TriOS::~TriOS() {}

  void TriOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *phase = in(Phase);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      float out;
      for (int k = 0; k < oversample.getOversamplingRatio(); k++)
        osBuffer[k] = saw.next(freq[i], phase[i] + 0.5, m_freqMul / oversample.getOversamplingRatio());
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

}

namespace VarSawOS
{

  VarSawOS::VarSawOS()
  {
    const float samplerate = sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<VarSawOS, &VarSawOS::next_aa>();
    next_aa(1);
  }

  VarSawOS::~VarSawOS() {}

  float VarSawOS::next(float freq, float phase, float width)
  {
    float out;
    float invwidth = 2.f / width;
    float inv1width = 2.f / (1 - width);

    for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      osBuffer[k] = saw.next(freq, phase - 0.5, m_freqMul / oversample.getOversamplingRatio());
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

  void VarSawOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *phase = in(Phase);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], phase[i], inWidth[i]);
    }
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

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<SquareOS, &SquareOS::next_aa>();
    next_aa(1);
  }

  SquareOS::~SquareOS() {}

  float SquareOS::next(float freq, float phase, float width)
  {
    float out;

    for (int k = 0; k < oversample.getOversamplingRatio(); k++)
      osBuffer[k] = saw.next(freq, phase, m_freqMul / oversample.getOversamplingRatio());
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

  void SquareOS::next_aa(int nSamples)
  {

    const float *freq = in(Freq);
    const float *phase = in(Phase);
    const float *inWidth = in(Width);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next(freq[i], phase[i], inWidth[i]);
    }
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

    if (out > 1.f)
      out = 1.f;
    else if (out < (-1.f))
      out = (-1.f);

    return out;
  }

  SawBL::SawBL()
  {

    m_counter = 0;

    mCalcFunc = make_calc_function<SawBL, &SawBL::next_a>();
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

}

namespace SquareBL
{

  SquareNext::SquareNext()
  {
  }

  SquareNext::~SquareNext() {}

  float SquareNext::next(float freq, float duty)
  {
    // freq = sc_clip(abs(freq), m_fmin*4, m_sampleRate/2);
    freq = abs(freq);
    float p0n = m_sampleRate / freq;
    float ddel = duty * p0n;
    float del = sc_max(0, sc_min(ddel, (float)delMax));

    float sawval = saw.next(freq, &m_phase, &m_counter, p0n, m_freqMul);

    delArray[delArrayCounter] = sawval;

    float dec = del - (long)del;
    float delaySig = (delArray[sc_wrap(delArrayCounter - (int)del, 0, delMax)] * (1 - dec)) +
                     (delArray[sc_wrap(delArrayCounter - ((int)del + 1), 0, delMax)] * dec);

    float out = (sawval - delaySig);

    delArrayCounter++;
    if (delArrayCounter >= (delMax + 1))
      delArrayCounter = 0;

    float lf_square = m_phase + ((0.5f - duty) * 2.f);
    if (lf_square >= 0)
      lf_square = 1.f;
    else
      lf_square = -1.f;
    lf_square -= ((0.5f - duty) * 2.f);

    if (freq < (m_fmin4))
    {
      if (freq < m_fmin2)
        out = lf_square;
      else
      {
        float mul = (freq - m_fmin2) / (m_fmin4 - m_fmin2);
        out = (out * mul) + (lf_square * (1 - mul));
      }
    }

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
    mCalcFunc = make_calc_function<SquareBL, &SquareBL::next_aa>();

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
}

namespace TriBL
{

  TriBL::TriBL()
  {
    square.setRatePhase(sampleRate(), m_phase);

    mCalcFunc = make_calc_function<TriBL, &TriBL::next_aa>();

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

PluginLoad(OversamplingOscillators)
{
  ft = inTable;
  registerUnit<SawOS::SawOS>(ft, "SawOS", false);
  registerUnit<SinOscOS::SinOscOS>(ft, "SinOscOS", false);
  registerUnit<TriOS::TriOS>(ft, "TriOS", false);
  registerUnit<VarSawOS::VarSawOS>(ft, "VarSawOS", false);
  registerUnit<SquareOS::SquareOS>(ft, "SquareOS", false);
  registerUnit<PMOscOS::PMOscOS>(ft, "PMOscOS", false);

  registerUnit<FM7OS::FM7OS>(ft, "FM7OS", false);
  registerUnit<FM7aOS::FM7aOS>(ft, "FM7aOS", false);
  registerUnit<FM7bOS::FM7bOS>(ft, "FM7bOS", false);
  registerUnit<PM7OS::PM7OS>(ft, "PM7OS", false);

  registerUnit<SawBL::SawBL>(ft, "SawBL", false);
  registerUnit<SquareBL::SquareBL>(ft, "SquareBL", false);
  registerUnit<TriBL::TriBL>(ft, "TriBL", false);
  registerUnit<ImpulseBL::ImpulseBL>(ft, "ImpulseBL", false);
}
