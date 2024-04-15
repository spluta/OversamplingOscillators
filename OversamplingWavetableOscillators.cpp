#include "OversamplingWavetableOscillators.hpp"
#include "SC_PlugIn.hpp"
#include "SC_PlugIn.h"
#include "sergeWavetable.h"

static InterfaceTable *ft;

namespace Extras {
  SawOSNext::SawOSNext()
  {
  }

  SawOSNext::~SawOSNext() {}

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

  void SawOSNext::reset(float phaseIn)
  {
    m_lastPhase = phaseIn;
    m_phase = phaseIn;
  }

  ProcessFuncs::ProcessFuncs()
  {
  }

  ProcessFuncs::~ProcessFuncs() {}

  float ProcessFuncs::get_phase(const float* phase_buf_data, float phase, float phase_buf_divs, float phase_buf_loc, int phase_table_size, float phase_fmaxindex) {

    int iphase_buf_divs = (int)phase_buf_divs;
    float fphase_index = (phase*phase_fmaxindex);
    int phase_index = (int)fphase_index;
    float frac_phase_index = phase_index - phase_index;

    float ramp = 0.f;

    phase_buf_loc = phase_buf_loc*(float)(phase_buf_divs-1);
    int iphase_buf_loc = (int)phase_buf_loc;

    //if the phase_buf_loc is the last index, we only look at that index
    if(iphase_buf_loc == iphase_buf_divs-1)
    {
      
      int zero_index = (iphase_buf_loc*phase_table_size);
      int loc = phase_index+zero_index;

      //if it is the last index, we only look at that index
      if(phase_index==phase_table_size-1)
        ramp = phase_buf_data[loc];
      else
        //otherwise interpolate between the indexa dn the next index
        ramp = phase_buf_data[loc]*(1.0f-frac_phase_index) + phase_buf_data[loc+1]*frac_phase_index;
 
    } else {
      float frac_phase_loc = phase_buf_loc - iphase_buf_loc;

      int final_index = phase_index + (iphase_buf_loc*phase_table_size);
      int final_index2 = final_index + phase_table_size;

      float ramp1, ramp2;
      if (phase_index==phase_table_size-1)
      {
        ramp1 = phase_buf_data[final_index]*(1.0f-frac_phase_index) + phase_buf_data[iphase_buf_loc*phase_table_size]*frac_phase_index;
        ramp2 = phase_buf_data[final_index2]*(1.0f-frac_phase_index) + phase_buf_data[(iphase_buf_loc+1)*phase_table_size]*frac_phase_index;
      } else {
        ramp1 = phase_buf_data[final_index]*(1.0f-frac_phase_index) + phase_buf_data[final_index+1]*frac_phase_index;
        ramp2 = phase_buf_data[final_index2]*(1.0f-frac_phase_index) + phase_buf_data[final_index2+1]*frac_phase_index;
      }
      //interpolate between the phase of the two tables
      ramp = ramp1*(1.0f-frac_phase_loc) + ramp2*frac_phase_loc;
    }
    return ramp;
  };

  float ProcessFuncs::get_out(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc) {
    //now that we have the ramp, use it to get the value from the buffer

    float findex = (ramp*fmaxindex); 
    int index = (int)findex;
    float frac = findex - index;
    int ibuf_divs = (int)buf_divs;
    
    float out;

    buf_loc = buf_loc*(float)(buf_divs-1);
    int ibuf_loc = (int)buf_loc;

    if(ibuf_loc == ibuf_divs-1)
    {
      int zero_index = (ibuf_loc*each_table_size);
      int loc = (index+zero_index)*num_chans + chan_loc;

      if (index==each_table_size-1)
        out = buf_data[loc]*(1.0f-frac) + buf_data[zero_index*num_chans]*frac;
      else
        out = buf_data[loc]*(1.0f-frac) + buf_data[loc+num_chans]*frac;
    } else {

      float frac_loc = buf_loc - ibuf_loc;

      int final_index = index + (ibuf_loc*each_table_size);
      int final_index2 = final_index + each_table_size;

      final_index = final_index*num_chans + chan_loc;
      final_index2 = final_index2*num_chans + chan_loc;

      float out1, out2;
      if (index==each_table_size-1)
      {
        out1 = buf_data[final_index]*(1.0f-frac) + buf_data[ibuf_loc*each_table_size*num_chans + chan_loc]*frac;
        out2 = buf_data[final_index2]*(1.0f-frac) + buf_data[(ibuf_loc+1)*each_table_size*num_chans + chan_loc]*frac;
      } else {
        out1 = buf_data[final_index]*(1.0f-frac) + buf_data[final_index+num_chans]*frac;
        out2 = buf_data[final_index2]*(1.0f-frac) + buf_data[final_index2+num_chans]*frac;
      }

      out = out1*(1.0f-frac_loc) + out2*frac_loc;
    }

    return out;
  }
}

namespace BuchlaFoldOS
{

  BuchlaFoldOS::BuchlaFoldOS()
  {
    const float samplerate = (float) sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    m_oversampling_ratio = oversample.getOversamplingRatio();

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<BuchlaFoldOS, &BuchlaFoldOS::next_aa>();
    next_aa(1);
  }

  BuchlaFoldOS::~BuchlaFoldOS() {}

  float BuchlaFoldOS::buchla_cell(float sig, float sign, float thresh, float sig_mul1, float sign_mul, float sig_mul2) {
    if (std::abs(sig) > thresh) {
        return (sig * sig_mul1 - (sign * sign_mul)) * sig_mul2;
    } else {
        return 0.0f;
    }
  }

  float BuchlaFoldOS::next(float sig, float amp) {
    float out;
    if(amp>0.f){
      sig = sig * amp;
      float sign = (sig >= 0.0f) ? 1.0f : -1.0f;
      float v1 = buchla_cell(sig, sign, 0.6f, 0.8333f, 0.5f, 12.0f);
      float v2 = buchla_cell(sig, sign, 2.994f, 0.3768f, 1.1281f, 27.777f);
      float v3 = buchla_cell(sig, sign, 5.46f, 0.2829f, 1.5446f, 21.428f);
      float v4 = buchla_cell(sig, sign, 1.8f, 0.5743f, 1.0338f, 17.647f);
      float v5 = buchla_cell(sig, sign, 4.08f, 0.2673f, 1.0907f, 36.363f);
      float v6 = sig * 5.0f;
      out = ((v1 + v2 + v3)*(-1.f))+ v4 + v5 + v6;
      out = out / 5.0f;
    } else {
      out = 0.f;
    }
    
    return out;
  }

  float BuchlaFoldOS::next_os(float sig, float amp)
  {
    float out;
    
    oversample.upsample(sig);

    for (int k = 0; k < m_oversampling_ratio; k++){
      osBuffer[k] = next(osBuffer[k], amp);
    }
    if (m_oversamplingIndex != 0)
      out = oversample.downsample();
    else
      out = osBuffer[0];
    return out;
  }

  void BuchlaFoldOS::next_aa(int nSamples)
  {

    const float *sig = in(Sig);
    const float *amp = in(Amp);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next_os(sig[i], amp[i]);
    }
  }
} // namespace BuchlaFoldOS


namespace SergeFoldOS {

  SergeFoldOS::SergeFoldOS()
  {
    const float samplerate = (float) sampleRate();

    sergeWavetable = getSergeWavetable();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);
  m_oversampling_ratio = oversample.getOversamplingRatio();

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<SergeFoldOS, &SergeFoldOS::next_aa>();
    next_aa(1);
  } 
  SergeFoldOS::~SergeFoldOS() {}

  float SergeFoldOS::next(float sig, float amp) {
    double out = tanh(sig*amp);
    float findex = ((out*0.5+0.5)*(float)sergeWavetable.size());
    float frac = findex - (int)findex;
    int index = (int)findex;

    if (index < 0) {
      index = 0;
    }
    else if (index > sergeWavetable.size()-2) {
      index = sergeWavetable.size()-2;
    }

    out = sergeWavetable[index]*(1.0f-frac) + sergeWavetable[index+1]*frac;

    out = tanh(out);

    return out;
  }

  float SergeFoldOS::next_os(float sig, float amp)
  {
    float out;
    
    oversample.upsample(sig);

    for (int k = 0; k < m_oversampling_ratio; k++){
      osBuffer[k] = next(osBuffer[k], amp);
    }
    if (m_oversamplingIndex != 0)
      out = oversample.downsample();
    else
      out = osBuffer[0];
    return out;
  }

  void SergeFoldOS::next_aa(int nSamples)
  {

    const float *sig = in(Sig);
    const float *amp = in(Amp);
    float *outbuf = out(Out1);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = next_os(sig[i], amp[i]);
    }
  }
} // namespace SergeFoldOS

namespace BufUnit {
  BufUnit::BufUnit()
  {
    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;
  }

  BufUnit::~BufUnit() {}

  bool BufUnit::GetTable(World* world, float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData,
                               int& tableSize) {
    
    if (fbufnum < 0.f) {                                                                                               \
        fbufnum = 0.f;                                                                                                 \
    }    

    if (fbufnum != m_fbufnum) {
        uint32 bufnum = (uint32)fbufnum;
        if (bufnum >= world->mNumSndBufs) {
            uint32 localBufNum = bufnum - world->mNumSndBufs;
            Graph* parent = mParent;
            if (localBufNum <= parent->localBufNum)
                m_buf = parent->mLocalSndBufs + localBufNum;
            else {
                bufnum = 0;
                m_buf = world->mSndBufs + bufnum;
            }
        } else
            m_buf = world->mSndBufs + bufnum;

        m_fbufnum = fbufnum;
    }
    buf = m_buf;
    if (!buf) {
        ClearUnitOutputs(this, inNumSamples);
        return false;
    }

    bufData = buf->data;
    if (!bufData) {
        ClearUnitOutputs(this, inNumSamples);
        return false;
    }
    tableSize = buf->samples;
    return true;
}
} // namespace BufUnit

namespace ShaperOS {

  ShaperOS::ShaperOS()
  {
    const float samplerate = (float) sampleRate();

    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;

    buf_unit = BufUnit::BufUnit();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

    m_oversampling_ratio = oversample.getOversamplingRatio();

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<ShaperOS, &ShaperOS::next_aa>();
    next_aa(1);
  } 
  ShaperOS::~ShaperOS() {}

 float ShaperOS::Perform(const float* table0, float in, float fmaxindex) {

  float findex = ((in*0.5+0.5)*fmaxindex);
    float frac = findex - (int)findex;
    int index = (int)findex;

    index = sc_clip(index, 0.f, fmaxindex-1.f);

    float out = table0[index]*(1.0f-frac) + table0[index+1]*frac;

  return out;
}

  float ShaperOS::next_os(const float* table0, float in, float fmaxindex)
  {
    float out;
    
    oversample.upsample(in);

    for (int k = 0; k < m_oversampling_ratio; k++){
      osBuffer[k] = Perform(table0, osBuffer[k], fmaxindex);
    }
    if (m_oversamplingIndex != 0)
      out = oversample.downsample();
    else
      out = osBuffer[0];
    return out;
  }

  void ShaperOS::next_aa(int nSamples)
  {
    const float *sig = in(Sig);
    const float buf_num = in0(BufNum);
    float *outbuf = out(Out1);
    
    // get table
    const SndBuf* buf; const float* bufData; int tableSize;
    
    const bool verify_buf = buf_unit.GetTable(mWorld, buf_num, nSamples, buf, bufData, tableSize);

    if (!verify_buf){
        ClearUnitOutputs(this, nSamples);
        return;
    }

    const float* table0 = bufData;
    float fmaxindex = (float)(tableSize) - 1.f/(float)(tableSize);

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = ShaperOS::next_os(table0, sig[i], fmaxindex);
    }
  }
} // namespace ShaperOS

namespace ShaperOS2 {

  ShaperOS2::ShaperOS2()
  {
    const float samplerate = (float) sampleRate();


    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;

    m_last_phase = 0.f;
    m_last_buf_loc = 0.f;

    buf_unit = BufUnit::BufUnit();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);
    osBuffer = oversample.getOSBuffer();

    upsample_buf_loc.reset(samplerate);
    upsample_buf_loc.setOversamplingIndex(m_oversamplingIndex);
    upsample_buf_ptr = upsample_buf_loc.getOSBuffer();

    upsample_input.reset(samplerate);
    upsample_input.setOversamplingIndex(m_oversamplingIndex);
    upsample_input_ptr = upsample_input.getOSBuffer();

    m_oversampling_ratio = oversample.getOversamplingRatio();

    mCalcFunc = make_calc_function<ShaperOS2, &ShaperOS2::next_aa>();
    next_aa(1);
  } 
  ShaperOS2::~ShaperOS2() {}

//  float ShaperOS2::Perform(const float* table0, float input, float buf_divs, float fbuf_loc, int table_size, float fmaxindex) {

//     float findex = (input*fmaxindex);
//     int index = (int)findex;
//     float frac = findex - index;
//     int ibuf_divs = (int)buf_divs;
//     float out;

//     fbuf_loc = fbuf_loc*(buf_divs-1);
//     int ibuf_loc = (int)fbuf_loc;

//     if(ibuf_loc == ibuf_divs-1)
//     {
//       int loc = index+(ibuf_loc*table_size);
//       if (index==table_size-1)
//         out = table0[loc];
//       else
//         out = table0[loc]*(1.0f-frac) + (table0[loc+1]*frac);
//     } else {

//       float frac_loc = fbuf_loc - ibuf_loc;

//       int final_index = index + (ibuf_loc*table_size);
//       int final_index2 = final_index + table_size;

//       float out1, out2;
//       if (index==table_size-1)
//       {
//         out1 = table0[final_index]*(1.0f-frac) + table0[ibuf_loc*table_size]*frac;
//         out2 = table0[final_index2]*(1.0f-frac) + table0[(ibuf_loc+1)*table_size]*frac;
//       } else {
//         out1 = table0[final_index]*(1.0f-frac) + table0[final_index+1]*frac;
//         out2 = table0[final_index2]*(1.0f-frac) + table0[final_index2+1]*frac;
//       }

//       out = out1*(1.0f-frac_loc) + out2*frac_loc;
//     }

//   return out;
// }

  float ShaperOS2::next_os(const float* table0, float input, float buf_divs, float buf_loc, int each_table_size, float fmaxindex)
  {
    float out;
    
    input = (input+1.f)*0.5f;
    //input = sc_clip(input, 0, 1.0f);

    //Print("input: %f\n", input);
    
    float buf_loc1 = sc_clip(buf_loc, 0.f, 1.0f);

    upsample_input.upsample(input);
    upsample_buf_loc.upsample(buf_loc1);

    for (int k = 0; k < m_oversampling_ratio; k++){
      //osBuffer[k] = Perform(table0, upsample_input_ptr[k], buf_divs, upsample_buf_ptr[k], table_size, fmaxindex);
      float val = sc_clip(upsample_input_ptr[k],0.0,1.0);
      osBuffer[k] = process_funcs.get_out(table0, val, buf_divs, buf_loc, each_table_size, fmaxindex, 1, 0); //only one channel
    }

    out = oversample.downsample();

    return out;
  }

  void ShaperOS2::next_aa(int nSamples)
  {
    const float *input = in(In);
    float buf_divs = floor(in0(BufDivs));
    const float *buf_loc = in(BufLoc);
    const float buf_num = in0(BufNum);
    float *outbuf = out(Out1);

    if (buf_divs < 1.f)
      buf_divs = 1.f;
    
    // get table
    const SndBuf* buf; const float* bufData; int tableSize;
    
    const bool verify_buf = buf_unit.GetTable(mWorld, buf_num, nSamples, buf, bufData, tableSize);

    if (!verify_buf){
        ClearUnitOutputs(this, nSamples);
        return;
    }

    const float* table0 = bufData;

    int each_table_size = tableSize/buf_divs;
    float fmaxindex = (float)each_table_size - 1.f;
     

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = ShaperOS2::next_os(table0, input[i], buf_divs, buf_loc[i], each_table_size, fmaxindex);
    }
  }
} // namespace ShaperOS2

namespace OscOS {

  OscOS::OscOS()
  {
    const float samplerate = (float) sampleRate();


    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;

    m_last_phase = 0.f;
    m_last_buf_loc = 0.f;

    buf_unit = BufUnit::BufUnit();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);
    osBuffer = oversample.getOSBuffer();

    upsample_buf_loc.reset(samplerate);
    upsample_buf_loc.setOversamplingIndex(m_oversamplingIndex);
    upsample_buf = upsample_buf_loc.getOSBuffer();

    m_oversampling_ratio = oversample.getOversamplingRatio();

    mCalcFunc = make_calc_function<OscOS, &OscOS::next_aa>();
    next_aa(1);
  } 
  OscOS::~OscOS() {}

 float OscOS::Perform(const float* table0, float phase, float buf_divs, float fbuf_loc, int table_size, float fmaxindex) {

    float findex = (phase*fmaxindex);
    int index = (int)findex;
    float frac = findex - index;
    int ibuf_divs = (int)buf_divs;
    float out;
    //if(ibuf_divs > 1){

    fbuf_loc = fbuf_loc*(buf_divs-1);
    int ibuf_loc = (int)fbuf_loc;

      if(ibuf_loc == ibuf_divs-1)
      {
        int loc = index+(ibuf_loc*table_size);
        if (index==table_size-1)
          out = table0[loc]*(1.0f-frac) + table0[ibuf_loc*table_size]*frac;
        else
          out = table0[loc]*(1.0f-frac) + table0[loc+1]*frac;
      } else {

        float frac_loc = fbuf_loc - ibuf_loc;

        int final_index = index + (ibuf_loc*table_size);
        int final_index2 = final_index + table_size;

        float out1, out2;
        if (index==table_size-1)
        {
          out1 = table0[final_index]*(1.0f-frac) + table0[ibuf_loc*table_size]*frac;
          out2 = table0[final_index2]*(1.0f-frac) + table0[(ibuf_loc+1)*table_size]*frac;
        } else {
          out1 = table0[final_index]*(1.0f-frac) + table0[final_index+1]*frac;
          out2 = table0[final_index2]*(1.0f-frac) + table0[final_index2+1]*frac;
        }

        out = out1*(1.0f-frac_loc) + out2*frac_loc;
      }

  return out;
}

  float OscOS::next_os(const float* table0, float phase, float buf_divs, float buf_loc, int table_size, float fmaxindex)
  {
    float out;
    
    float phase1 = sc_clip(phase, 0.f, 1.0f);
    float buf_loc1 = sc_clip(buf_loc, 0.f, 1.0f);

    float phase_diff = (phase1 - m_last_phase);

    upsample_buf_loc.upsample(buf_loc1);

  //the phase_diff should not be more than 0.5 except when the phase crosses from 1 to 0 or vice versa
  //even at the nyquist frequency, the phase_diff should not be more than 0.5
  if(abs(phase_diff) > 0.5f){
    if (phase1<m_last_phase)
      phase_diff = (phase1 + 1.0f) - m_last_phase;
    else
      phase_diff = (phase1 - 1.0f) - m_last_phase;
    phase_diff = phase_diff/m_oversampling_ratio;
    for (int k = 0; k < m_oversampling_ratio; k++){
      m_last_phase += phase_diff;
      osBuffer[k] = Perform(table0, sc_wrap(m_last_phase, 0.f, 1.0f), buf_divs, upsample_buf[k], table_size, fmaxindex);
    }
  }  
  else {
    phase_diff = phase_diff/m_oversampling_ratio;
    for (int k = 0; k < m_oversampling_ratio; k++){
      osBuffer[k] = Perform(table0, m_last_phase+(k*phase_diff), buf_divs, upsample_buf[k], table_size, fmaxindex);
    }
  }

    m_last_phase = phase1;
    m_last_buf_loc = buf_loc1;

    if (m_oversamplingIndex != 0)
      out = oversample.downsample();
    else
      out = osBuffer[0];
    return out;
  }

  void OscOS::next_aa(int nSamples)
  {
    const float *phase = in(Phase);
    float buf_divs = floor(in0(BufDivs));
    const float *buf_loc = in(BufLoc);
    const float buf_num = in0(BufNum);
    float *outbuf = out(Out1);

    if (buf_divs < 1.f)
      buf_divs = 1.f;
    
    // get table
    const SndBuf* buf; const float* bufData; int tableSize;
    
    const bool verify_buf = buf_unit.GetTable(mWorld, buf_num, nSamples, buf, bufData, tableSize);

    if (!verify_buf){
        ClearUnitOutputs(this, nSamples);
        return;
    }

    const float* table0 = bufData;

    int each_table_size = tableSize/buf_divs;
    float fmaxindex = (float)each_table_size - 1.f/(float)(each_table_size);
     

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = OscOS::next_os(table0, phase[i], buf_divs, buf_loc[i], each_table_size, fmaxindex);
    }
  }
} // namespace OscOS


namespace OscOS2 {

  OscOS2::OscOS2()
  {
    const float samplerate = (float) sampleRate();


    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;

    buf_unit = BufUnit::BufUnit();
    phase_buf_unit = BufUnit::BufUnit();

    m_oversampling_index = sc_clip((int)in0(OverSample), 0, 4);

    oversample.reset(samplerate);
    oversample.setOversamplingIndex(m_oversampling_index);
    os_buffer = oversample.getOSBuffer();

    upsample_buf_loc.reset(samplerate);
    upsample_buf_loc.setOversamplingIndex(m_oversampling_index);
    os_buf_loc = upsample_buf_loc.getOSBuffer();

    upsample_phase_buf_loc.reset(samplerate);
    upsample_phase_buf_loc.setOversamplingIndex(m_oversampling_index);
    os_phase_buf_loc = upsample_phase_buf_loc.getOSBuffer();

    m_oversampling_ratio = oversample.getOversamplingRatio();

    mCalcFunc = make_calc_function<OscOS2, &OscOS2::next_aa>();
    next_aa(1);
  } 
  OscOS2::~OscOS2() {}

  

 float OscOS2::Perform(const float* buf_data, const float* phase_buf_data, float phase, float buf_divs, float buf_loc, float phase_buf_divs, float phase_buf_loc, int each_table_size, float fmaxindex, int phase_table_size, float phase_fmaxindex) {

    float ramp;

    //if the phase_buf is nil, just use the phase as the ramp
    if(phase_buf_divs < 0.f)
      ramp = phase;
    else
      ramp = process_funcs.get_phase(phase_buf_data, phase, phase_buf_divs, phase_buf_loc, phase_table_size, phase_fmaxindex);

    float out = process_funcs.get_out(buf_data, ramp, buf_divs, buf_loc, each_table_size, fmaxindex, 1, 0); //only one channel

    return out;
}

  float OscOS2::next_os(const float* buf_data, const float* phase_buf_data, const float freq, const float phase, float buf_divs, float buf_loc, float phase_buf_divs, float phase_buf_loc, int each_table_size, float fmaxindex, int phase_table_size, float phase_fmaxindex)
  {
    float out;

    buf_loc = sc_clip(buf_loc, 0.f, 1.0f);
    phase_buf_loc = sc_clip(phase_buf_loc, 0.f, 1.0f);

    upsample_buf_loc.upsample(buf_loc);
    upsample_phase_buf_loc.upsample(phase_buf_loc);

    for (int k = 0; k < m_oversampling_ratio; k++){
      float saw_phase = saw.next(freq, phase, m_freqMul / m_oversampling_ratio)*0.5+0.5;
      os_buffer[k] = Perform(buf_data, phase_buf_data, saw_phase, buf_divs, os_buf_loc[k], phase_buf_divs, os_phase_buf_loc[k], each_table_size, fmaxindex, phase_table_size, phase_fmaxindex);
    } 

    if (m_oversampling_index != 0)
      out = oversample.downsample();
    else
      out = os_buffer[0];
    return out;
  }

  void OscOS2::next_aa(int n_samples)
  {
    const float buf_num = in0(BufNum);
    const float phase_buf_num = in0(PhaseBuf);

    const float *freq = in(Freq);
    const float *phase = in(Phase);

    float buf_divs = floor(in0(BufDivs));
    const float *buf_loc = in(BufLoc);

    float phase_buf_divs = floor(in0(PhaseBufDivs));
    const float *phase_buf_loc = in(PhaseBufLoc);
    
    float *outbuf = out(Out1);

    if (buf_divs < 1.f)
      buf_divs = 1.f;
    if (phase_buf_divs < 1.f)
      phase_buf_divs = 1.f;
    
    // get table
    const SndBuf* buf; const float* buf_data; int table_size;
    const SndBuf* phase_buf; const float* phase_buf_data; int phase_table_size;
    
    bool verify_buf = buf_unit.GetTable(mWorld, buf_num, n_samples, buf, buf_data, table_size);
    
    if (!verify_buf){
        ClearUnitOutputs(this, n_samples);
        return;
    }

    verify_buf = phase_buf_unit.GetTable(mWorld, phase_buf_num, n_samples, phase_buf, phase_buf_data, phase_table_size);

    if (!verify_buf){
        ClearUnitOutputs(this, n_samples);
        return;
    }

    int each_table_size = table_size/buf_divs;
    float fmaxindex = (float)each_table_size - 1.f/(float)(each_table_size);
     
    int each_phase_table_size = phase_table_size/phase_buf_divs;
    int phase_fmaxindex = (float)each_phase_table_size - 1.f;
    

    //Print("table_size: %i, each_table_size: %i, each_phase_table_size: %i, buf_divs: %f\n", table_size, each_table_size, each_phase_table_size, buf_divs);

    for (int i = 0; i < n_samples; ++i)
    {
      outbuf[i] = OscOS2::next_os(buf_data, phase_buf_data, freq[i], phase[i], buf_divs, buf_loc[i], phase_buf_divs, phase_buf_loc[i], each_table_size, fmaxindex, each_phase_table_size, phase_fmaxindex);
    }
  }
} // namespace OscOS2


namespace OscOS3 {

OscOS3::OscOS3()
  {
    const float samplerate = (float) sampleRate();

    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;

    buf_unit = BufUnit::BufUnit();
    phase_buf_unit = BufUnit::BufUnit();

    m_sync_trig = 0.f;

    m_oversampling_index = sc_clip((int)in0(OverSample), 0, 4);

    oversample.reset(samplerate);
    oversample.setOversamplingIndex(m_oversampling_index);
    os_buffer = oversample.getOSBuffer();

    upsample_buf_loc.reset(samplerate);
    upsample_buf_loc.setOversamplingIndex(m_oversampling_index);
    os_buf_loc = upsample_buf_loc.getOSBuffer();

    upsample_chan_loc.reset(samplerate);
    upsample_chan_loc.setOversamplingIndex(m_oversampling_index);
    os_chan_loc = upsample_chan_loc.getOSBuffer();

    upsample_phase_buf_loc.reset(samplerate);
    upsample_phase_buf_loc.setOversamplingIndex(m_oversampling_index);
    os_phase_buf_loc = upsample_phase_buf_loc.getOSBuffer();

    m_oversampling_ratio = oversample.getOversamplingRatio();

    mCalcFunc = make_calc_function<OscOS3, &OscOS3::next_aa>();
    next_aa(1);
  } 

  OscOS3::~OscOS3() {}


//Perform(buf_data, phase_buf_data, saw_phase, buf_divs, upsample_buf_loc[k], phase_buf_divs, upsample_phase_buf_loc[k], table_size, fmaxindex);
 float OscOS3::Perform(const float* buf_data, const float* phase_buf_data, float phase, float buf_divs, float buf_loc, float num_chans, float chan_loc, float phase_buf_divs, float phase_buf_loc, int each_table_size, float fmaxindex, int phase_table_size, float phase_fmaxindex) 
 {

    float ramp;
    
    //if the phase_buf is nil, just use the phase as the ramp
    if(phase_buf_divs < 0.f)
      ramp = phase;
    else
      ramp = process_funcs.get_phase(phase_buf_data, phase, phase_buf_divs, phase_buf_loc, phase_table_size, phase_fmaxindex);

    float out = process_funcs.get_out(buf_data, ramp, buf_divs, buf_loc, each_table_size, fmaxindex, (int)num_chans, chan_loc); 

  return out;
}

float OscOS3::next_os(const float* buf_data, const float* phase_buf_data, const float freq, const float phase, float buf_divs, float buf_loc, float num_chans, float chan_loc, float phase_buf_divs, float phase_buf_loc, int each_table_size, float fmaxindex, int each_phase_table_size, float phase_fmaxindex)
  {
    float out;

    //Print("next_os\n");

    buf_loc = sc_clip(buf_loc, 0.f, 1.0f);
    chan_loc = sc_clip(chan_loc, 0.f, 1.0f-1.f/num_chans);
    phase_buf_loc = sc_clip(phase_buf_loc, 0.f, 1.0f);

    //Print("buf_loc: %f, chan_loc: %f, phase_buf_loc: %f\n", buf_loc, chan_loc, phase_buf_loc);

    //upsample the inputs

    upsample_buf_loc.upsample(buf_loc);
    upsample_chan_loc.upsample(chan_loc);
    upsample_phase_buf_loc.upsample(phase_buf_loc);


    for (int k = 0; k < m_oversampling_ratio; k++){
      float saw_phase = saw.next(freq, phase, m_freqMul / m_oversampling_ratio)*0.5+0.5;

      if (num_chans <= 1.f){
        os_buffer[k] = Perform(buf_data, phase_buf_data, saw_phase, buf_divs, os_buf_loc[k], num_chans, 0.f, phase_buf_divs, os_phase_buf_loc[k], each_table_size, fmaxindex, each_phase_table_size, phase_fmaxindex);
      } else {
        float full_chan_loc = os_chan_loc[k]*(num_chans-1.f);
        float low_chan_loc = floor(full_chan_loc);
        float high_chan_loc = ceil(full_chan_loc);
        float frac = full_chan_loc - low_chan_loc;

        float chan0 = Perform(buf_data, phase_buf_data, saw_phase, buf_divs, os_buf_loc[k], num_chans, low_chan_loc, phase_buf_divs, os_phase_buf_loc[k], each_table_size, fmaxindex, each_phase_table_size, phase_fmaxindex);

        float chan1 = Perform(buf_data, phase_buf_data, saw_phase, buf_divs, os_buf_loc[k], num_chans, high_chan_loc, phase_buf_divs, os_phase_buf_loc[k], each_table_size, fmaxindex, each_phase_table_size, phase_fmaxindex);

        os_buffer[k] = chan0*(1.0f-frac) + chan1*frac;
      }
    } 

    out = oversample.downsample();
    return out;
  }

  void OscOS3::next_aa(int n_samples)
  {
    const float buf_num = in0(BufNum);
    float phase_buf_num = in0(PhaseBuf);

    const float *freq = in(Freq);
    const float *phase = in(Phase);

    const float *sync_trig = in(SyncTrig);

    float buf_divs = floor(in0(BufDivs));
    const float *buf_loc = in(BufLoc);

    float num_chans = floor(in0(NumChans));
    const float *chan_loc = in(ChanLoc);

    float phase_buf_divs = floor(in0(PhaseBufDivs));
    const float *phase_buf_loc = in(PhaseBufLoc);
    
    float *outbuf = out(Out1);

    if (buf_divs < 1.f)
      buf_divs = 1.f;
    if (phase_buf_divs < 1.f)
      phase_buf_divs = 1.f;
    
    // get table
    const SndBuf* buf; const float* buf_data; int table_size;
    const SndBuf* phase_buf; const float* phase_buf_data; int phase_table_size;

    bool verify_buf = buf_unit.GetTable(mWorld, buf_num, n_samples, buf, buf_data, table_size);
    
    if (!verify_buf){
        //figure out what unit is
        ClearUnitOutputs(this, n_samples);
        return;
    }

    int each_phase_table_size = 0;
    float phase_fmaxindex = 0;

    if(phase_buf_num < 0.f){
      //Print("phase_buf_num is less than 0\n");
      phase_buf_num = -1.f;
      phase_buf_divs = -1.f;
    } else {
      //Print("phase_buf_num is greater than 0\n");
      verify_buf = phase_buf_unit.GetTable(mWorld, phase_buf_num, n_samples, phase_buf, phase_buf_data, phase_table_size);

      if (!verify_buf){
          phase_buf_divs = -1.f;
          return;
      } else {
        each_phase_table_size = phase_table_size/phase_buf_divs;
        phase_fmaxindex = (float)each_phase_table_size - 1.f;
      }
    }

    //Print("each_phase_table_size: %i, phase_buf_divs: %f\n", each_phase_table_size, phase_buf_divs);

    //fmaxindex and each_table_size are size of a single channel of a single table
    int each_table_size = table_size/(buf_divs*num_chans); 
    float fmaxindex = (float)each_table_size - 1.f/(float)each_table_size;

    for (int i = 0; i < n_samples; ++i)
    {
      if((m_sync_trig <= 0.f)&&(sync_trig[i] > 0.f)){
        saw.reset(phase[i]);
      };
      m_sync_trig = sync_trig[i];
      outbuf[i] = OscOS3::next_os(buf_data, phase_buf_data, freq[i], phase[i], buf_divs, buf_loc[i], num_chans, chan_loc[i], phase_buf_divs, phase_buf_loc[i], each_table_size, fmaxindex, each_phase_table_size, phase_fmaxindex);
    }
  }
} // namespace OscOS3

PluginLoad(OversamplingOscillators)
{
  ft = inTable;
  registerUnit<BuchlaFoldOS::BuchlaFoldOS>(ft, "BuchlaFoldOS", false);
  registerUnit<SergeFoldOS::SergeFoldOS>(ft, "SergeFoldOS", false);
  registerUnit<ShaperOS::ShaperOS>(ft, "ShaperOS", false);
  registerUnit<ShaperOS2::ShaperOS2>(ft, "ShaperOS2", false);
  registerUnit<OscOS::OscOS>(ft, "OscOS", false);
  registerUnit<OscOS2::OscOS2>(ft, "OscOS2", false);
  registerUnit<OscOS3::OscOS3>(ft, "OscOS3", false);
}

/*
void fillTables(WaveTableOsc *osc, double *freqWaveRe, double *freqWaveIm, int numSamples) {
    int idx;
    
    // zero DC offset and Nyquist
    freqWaveRe[0] = freqWaveIm[0] = 0.0;
    freqWaveRe[numSamples >> 1] = freqWaveIm[numSamples >> 1] = 0.0;
    
    // determine maxHarmonic, the highest non-zero harmonic in the wave
    int maxHarmonic = numSamples >> 1;
    const double minVal = 0.000001; // -120 dB
    while ((fabs(freqWaveRe[maxHarmonic]) + fabs(freqWaveIm[maxHarmonic]) < minVal)
        && maxHarmonic) --maxHarmonic;

    // calculate topFreq for the initial wavetable
    // maximum non-aliasing playback rate is 1 / (2 * maxHarmonic), but we allow
    // aliasing up to the point where the aliased harmonic would meet the next
    // octave table, which is an additional 1/3
    double topFreq = 2.0 / 3.0 / maxHarmonic;
    
    // for subsquent tables, double topFreq and remove upper half of harmonics
    double *ar = new double [numSamples];
    double *ai = new double [numSamples];
    double scale = 0.0;
    while (maxHarmonic) {
        // fill the table in with the needed harmonics
        for (idx = 0; idx < numSamples; idx++)
            ar[idx] = ai[idx] = 0.0;
        for (idx = 1; idx <= maxHarmonic; idx++) {
            ar[idx] = freqWaveRe[idx];
            ai[idx] = freqWaveIm[idx];
            ar[numSamples - idx] = freqWaveRe[numSamples - idx];
            ai[numSamples - idx] = freqWaveIm[numSamples - idx];
        }
        
        // make the wavetable
        scale = makeWaveTable(osc, numSamples, ar, ai, scale, topFreq);

        // prepare for next table
        topFreq *= 2;
        maxHarmonic >>= 1;
    }
}

*/