#include "OversamplingWavefolders.hpp"
#include "SC_PlugIn.hpp"
#include "SC_PlugIn.h"
#include "sergeWavetable.h"

static InterfaceTable *ft;

namespace BuchlaFoldOS
{

  BuchlaFoldOS::BuchlaFoldOS()
  {
    const float samplerate = (float) sampleRate();

    oversample.reset(samplerate);
    m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
    oversample.setOversamplingIndex(m_oversamplingIndex);

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

    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
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

    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
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

    osBuffer = oversample.getOSBuffer();

    mCalcFunc = make_calc_function<ShaperOS, &ShaperOS::next_aa>();
    next_aa(1);
  } 
  ShaperOS::~ShaperOS() {}

//   bool ShaperOS::GetTable(float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData,
//                                int& tableSize) {
    
//     if (fbufnum != m_fbufnum) {
//         uint32 bufnum = (uint32)fbufnum;
//         World* world = mWorld;
//         if (bufnum >= world->mNumSndBufs) {
//             uint32 localBufNum = bufnum - world->mNumSndBufs;
//             Graph* parent = mParent;
//             if (localBufNum <= parent->localBufNum)
//                 m_buf = parent->mLocalSndBufs + localBufNum;
//             else {
//                 bufnum = 0;
//                 m_buf = world->mSndBufs + bufnum;
//             }
//         } else
//             m_buf = world->mSndBufs + bufnum;

//         m_fbufnum = fbufnum;
//     }
//     buf = m_buf;
//     if (!buf) {
//         ClearUnitOutputs(this, inNumSamples);
//         return false;
//     }

//     bufData = buf->data;
//     if (!bufData) {
//         ClearUnitOutputs(this, inNumSamples);
//         return false;
//     }
//     tableSize = buf->samples;
//     return true;
// }

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

    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
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
    
    buf_unit.GetTable(mWorld, buf_num, nSamples, buf, bufData, tableSize);

    const float* table0 = bufData;
    float fmaxindex = (float)(tableSize) - 0.001;

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = ShaperOS::next_os(table0, sig[i], fmaxindex);
    }
  }
} // namespace ShaperOS

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
    if(ibuf_divs > 1){

      fbuf_loc = fbuf_loc*(buf_divs-1);
      int buf_loc = (int)fbuf_loc;

      if(buf_loc == buf_divs-1)
      {
        int loc = index+((buf_loc)*table_size);
        if (index==table_size-1)
          out = table0[loc]*(1.0f-frac) + table0[buf_loc*table_size]*frac;
        else
          out = table0[loc]*(1.0f-frac) + table0[loc+1]*frac;
      } else {

        float frac_loc = fbuf_loc - buf_loc;

        int final_index = index + (buf_loc*table_size);
        int final_index2 = final_index + table_size;

        float out1, out2;
        if (index==table_size-1)
        {
          out1 = table0[final_index]*(1.0f-frac) + table0[buf_loc*table_size]*frac;
          out2 = table0[final_index2]*(1.0f-frac) + table0[(buf_loc+1)*table_size]*frac;
        } else {
          out1 = table0[final_index]*(1.0f-frac) + table0[final_index+1]*frac;
          out2 = table0[final_index2]*(1.0f-frac) + table0[final_index2+1]*frac;
        }

        out = out1*(1.0f-frac_loc) + out2*frac_loc;
      }
    } else {
      out = table0[index]*(1.0f-frac) + table0[index+1]*frac;
    }

  return out;
}

  float OscOS::next_os(const float* table0, float phase, float buf_divs, float buf_loc, int table_size, float fmaxindex)
  {
    float out;
    
    float phase1 = sc_clip(phase, 0.f, 1.0f);
    float buf_loc1 = sc_clip(buf_loc, 0.f, 1.0f);

    float phase_diff = (phase1 - m_last_phase)/oversample.getOversamplingRatio();
    float loc_diff = (buf_loc1 - m_last_buf_loc)/oversample.getOversamplingRatio();

    
    //Print("m_last_phase: %f, phase1: %f, phase_diff: %f, m_last_buf_loc: %f, buf_loc1: %f, loc_diff: %f\n", m_last_phase, phase1, phase_diff, m_last_buf_loc, buf_loc1, loc_diff);
  //m_last_phase+(k*phase_diff)
    for (int k = 0; k < oversample.getOversamplingRatio(); k++){
      osBuffer[k] = Perform(table0, phase1, buf_divs, buf_loc1, table_size, fmaxindex);
    }

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
    
    buf_unit.GetTable(mWorld, buf_num, nSamples, buf, bufData, tableSize);

    const float* table0 = bufData;
    float fmaxindex = (float)(tableSize)/buf_divs - 0.001;
    int each_table_size = tableSize/buf_divs; 

    for (int i = 0; i < nSamples; ++i)
    {
      outbuf[i] = OscOS::next_os(table0, phase[i], buf_divs, buf_loc[i], each_table_size, fmaxindex);
    }
  }
} // namespace OscOS



PluginLoad(OversamplingOscillators)
{
  ft = inTable;
  registerUnit<BuchlaFoldOS::BuchlaFoldOS>(ft, "BuchlaFoldOS", false);
  registerUnit<SergeFoldOS::SergeFoldOS>(ft, "SergeFoldOS", false);
  registerUnit<ShaperOS::ShaperOS>(ft, "ShaperOS", false);
  registerUnit<OscOS::OscOS>(ft, "OscOS", false);
}