// OversamplingOscillators.cpp
// Sam Pluta
// uses the VariableOversampling class from Jatin Chowdhurry's ChowDSP library in the ported plugins library by Mads Kjeldgaard
// edit the c_cpp_properties.json and CMakeLists.txt files to correctly point to these files on your system

#include "OversamplingOscillators.hpp"
#include "SC_PlugIn.hpp"
#include "SC_PlugIn.h"

static InterfaceTable *ft;

namespace OSaw2 {

OSaw2::OSaw2() {
  const float samplerate = sampleRate();

  oversample.reset(samplerate);
  m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
  oversample.setOversamplingIndex(m_oversamplingIndex);
  
  if(inRate(0)==2)
    mCalcFunc = make_calc_function<OSaw2, &OSaw2::next_a>(); 
  else 
    mCalcFunc = make_calc_function<OSaw2, &OSaw2::next_k>(); 
  next_a(1);
}

OSaw2::~OSaw2() {}

void OSaw2::next_a(int nSamples) {

  const float *freq = in(Freq); 
  float *outbuf = out(Out1);

  int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
  if(osIndexIn!=m_oversamplingIndex){
    m_oversamplingIndex = osIndexIn;
    oversample.setOversamplingIndex(m_oversamplingIndex);
  }

  float *osBuffer = oversample.getOSBuffer();

  for (int i = 0; i < nSamples; ++i) {
    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
      float z = m_phase;
      m_phase += freq[i] * m_freqMul / (float)oversample.getOversamplingRatio();
      if (m_phase >= 1.f) 
        m_phase -= 2.f;
        else if (m_phase <= -1.f) m_phase += 2.f;
      osBuffer[k] = z;
    }
    if(m_oversamplingIndex!=0){
      float y = oversample.downsample();
      outbuf[i] = y;
    } else {outbuf[i] = osBuffer[0];}
  }
}

void OSaw2::next_k(int nSamples) {

  SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
  float *outbuf = out(Out1);

  int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
  if(osIndexIn!=m_oversamplingIndex){
    m_oversamplingIndex = osIndexIn;
    oversample.setOversamplingIndex(m_oversamplingIndex);
  }
  
  float *osBuffer = oversample.getOSBuffer();

  for (int i = 0; i < nSamples; ++i) {
    const float freq = slopedFreq.consume();
    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
      float z = m_phase;
      m_phase += freq * m_freqMul / (float)oversample.getOversamplingRatio();
      if (m_phase >= 1.f) 
        m_phase -= 2.f;
        else if (m_phase <= -1.f) m_phase += 2.f;
      osBuffer[k] = z;
    }
    if(m_oversamplingIndex!=0){
      float y = oversample.downsample();
      outbuf[i] = y;
    } else {outbuf[i] = osBuffer[0];}
  }
   m_freq_past = slopedFreq.value;
}
}

namespace OVarSaw2 {

OVarSaw2::OVarSaw2() {
  const float samplerate = sampleRate();

  oversample.reset(samplerate);
  m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
  oversample.setOversamplingIndex(m_oversamplingIndex);
  
  if(inRate(0)==2)
    mCalcFunc = make_calc_function<OVarSaw2, &OVarSaw2::next_a>(); 
  else 
    mCalcFunc = make_calc_function<OVarSaw2, &OVarSaw2::next_k>(); 
  next_a(1);
}

OVarSaw2::~OVarSaw2() {}

void OVarSaw2::next_a(int nSamples) {

  const float *freq = in(Freq); 

  float *outbuf = out(Out1);
  float nextDuty = sc_clip(in0(Duty), 0.001, 0.999);
  float duty = m_duty;
  float invduty = 2.f/duty;
  float inv1duty = 2.f/(1-duty);

  int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
  if(osIndexIn!=m_oversamplingIndex){
    m_oversamplingIndex = osIndexIn;
    oversample.setOversamplingIndex(m_oversamplingIndex);
  }

  //float y;
  float *osBuffer = oversample.getOSBuffer();

  for (int i = 0; i < nSamples; ++i) {
    float phaseInc = freq[i] * m_freqMul / (float)oversample.getOversamplingRatio();

    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
      if (m_phase >= 1.f) {
        m_phase -= 1.f;
        duty = m_duty = sc_clip(nextDuty, 0.001, 0.999);
        invduty = m_invduty = 2.f / duty;
        inv1duty = m_inv1duty = 2.f / (1.f - duty);
      }
       
      float z = m_phase < duty ? m_phase * invduty : (1.f - m_phase) * inv1duty;
      osBuffer[k] = z;
      m_phase+=phaseInc;
    }
    if(m_oversamplingIndex!=0){
      float y = oversample.downsample();
      outbuf[i] = y-1.f;
    } else {outbuf[i] = osBuffer[0];}
  }
}

void OVarSaw2::next_k(int nSamples) {

  SlopeSignal<float> slopedFreq = makeSlope(in0(Freq), m_freq_past);
  
  float *outbuf = out(Out1);
  float nextDuty = sc_clip(in0(Duty), 0.001, 0.999);
  float duty = m_duty;
  float invduty = 2.f/duty;
  float inv1duty = 2.f/(1-duty);

  int osIndexIn = sc_clip((int)in0(OverSample), 0, 4);
  if(osIndexIn!=m_oversamplingIndex){
    m_oversamplingIndex = osIndexIn;
    oversample.setOversamplingIndex(m_oversamplingIndex);
  }
  
  float *osBuffer = oversample.getOSBuffer();

  for (int i = 0; i < nSamples; ++i) {
    const float freq = slopedFreq.consume();
    float phaseInc = freq * m_freqMul / oversample.getOversamplingRatio();
    

    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
      if (m_phase >= 1.f) {
        m_phase -= 1.f;
        duty = m_duty = sc_clip(nextDuty, 0.001, 0.999);
        invduty = m_invduty = 2.f / duty;
        inv1duty = m_inv1duty = 2.f / (1.f - duty);
      }
       
      float z = m_phase < duty ? m_phase * invduty : (1.f - m_phase) * inv1duty;
      osBuffer[k] = z;
      m_phase+=phaseInc;
    }
    if(m_oversamplingIndex!=0){
      float y = oversample.downsample();
      outbuf[i] = y-1.f;
    } else {outbuf[i] = osBuffer[0];}
  }
  m_freq_past = slopedFreq.value;
}
}

PluginLoad(OSaw2UGens) {
  ft = inTable;
  registerUnit<OSaw2::OSaw2>(ft, "OSaw2", false);
  registerUnit<OVarSaw2::OVarSaw2>(ft, "OVarSaw2", false);
}
