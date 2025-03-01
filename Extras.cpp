#include "Extras.hpp"
#include "wavetables.h"
#include "SC_PlugIn.hpp"

float quadratic_interpolation(float y0, float y1, float y2, float x) {
  // Calculate the coefficients of the quadratic polynomial
  float a = ((x - 1) * (x - 2)) * 0.5f * y0;
  float b = (x * (x - 2)) * (-1.0f) * y1;
  float c = (x * (x - 1)) * 0.5f * y2;

  // Return the estimated value
  return a + b + c;
}

namespace SincExtras {

  SincFuncs::SincFuncs(int ripples)
  {
    if (ripples==32){
      m_sinc_table32 = get_sinc_window32();
      m_sinc_table = m_sinc_table32.data();
      sinc_points32 = {0, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192, 9216, 10240, 11264, 12288, 13312, 14336, 15360, 16384, 17408, 18432, 19456, 20480, 21504, 22528, 23552, 24576, 25600, 26624, 27648, 28672, 29696, 30720, 31744};
      sinc_points = sinc_points32.data();
      sinc_len=32;
      sinc_half_len=16;
      sinc_table_size=32768;
    } else {
      m_sinc_table8 = get_sinc_window8();
      m_sinc_table = m_sinc_table8.data();
      sinc_len=8;
      sinc_half_len=4;
      sinc_table_size=8192;
      sinc_points8 = {0, 1024, 2048, 3072, 4096, 5120, 6144, 7168};
      sinc_points = sinc_points8.data();
    }

  }

  SincFuncs::~SincFuncs() {}

  float SincFuncs::get_spaced_out(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int spacing1, float sinc_crossfade, int num_chans, int chan_loc) {
    //now that we have the ramp, use it to get the value from the buffer

    int spacing2 = spacing1*2;

    buf_loc = buf_loc*(float)(buf_divs-1);
    int ibuf_loc = (int)buf_loc;

    float findex = (ramp*fmaxindex);
    int ibuf_divs = (int)buf_divs;
    float out;

    float frac_loc = buf_loc - ibuf_loc;

    float sinc1 = get_spaced_sinc_sum(buf_data, each_table_size, findex, spacing1, ibuf_loc, num_chans, chan_loc);
    float sinc2 = 0.f;
    float outA = 0.f;
    //we only need to calculate sinc2 if we aren't in the top octave
    if(spacing1<max_sinc_offset){
      sinc2 = get_spaced_sinc_sum(buf_data, each_table_size, findex, spacing2, ibuf_loc, num_chans, chan_loc);
      outA = sinc1*(1.0f-sinc_crossfade) + sinc2*sinc_crossfade;
    } else {
      outA = sinc1;
    }

    float outB = 0.f;
    if(ibuf_loc < ibuf_divs-1)
    {
      sinc1 = get_spaced_sinc_sum(buf_data, each_table_size, findex, spacing1, ibuf_loc+1, num_chans, chan_loc);
      if(spacing1<max_sinc_offset){
        sinc2 = get_spaced_sinc_sum(buf_data, each_table_size, findex, spacing2, ibuf_loc+1, num_chans, chan_loc);
        outB = sinc1*(1.0f-sinc_crossfade) + sinc2*sinc_crossfade;
      } else {
        outB = sinc1;
      }
    }

    out = outA*(1.0f-frac_loc) + outB*frac_loc;
    return out;
  }

  float SincFuncs::get_spaced_sinc_sum(const float* table, int table_size, float findex, int spacing, int ibuf_loc, int num_chans, int chan_loc){
    
    float sinc_sum=0.f;

    //zero_index is the index of the first sample of the small table inside the big table
    int zero_index = (ibuf_loc*table_size);
    int sinc_mult = max_sinc_offset/spacing;
    int index = (int)findex;
    float frac = findex - index;
    
    for(int sp=0; sp<sinc_len; sp++){

      //the exact point along the 1D table
      int loc_point=index+((sp-sinc_half_len)*spacing);
      loc_point = sc_wrap(loc_point, 0, table_size-1);
      //round to the nearest spacing
      int spaced_point = (loc_point/spacing)*spacing;

      //the offset from the exact point to the nearest spacing
      int sinc_offset2 = (loc_point-spaced_point);

      //add the zero index to shift to the correct table
      spaced_point = spaced_point+zero_index;

      //quadratic interpolation
      int sinc_indexA = (sinc_points[sp]-(sinc_offset2*sinc_mult));
      int sinc_indexB = sinc_indexA - 1;
      int sinc_indexC = sinc_indexA - 2;

      sinc_indexA = sc_wrap(sinc_indexA, 0, sinc_table_size-1);
      sinc_indexB = sc_wrap(sinc_indexB, 0, sinc_table_size-1);
      sinc_indexC = sc_wrap(sinc_indexC, 0, sinc_table_size-1);

      float sinc_val = quadratic_interpolation(m_sinc_table[sinc_indexA], m_sinc_table[sinc_indexB], m_sinc_table[sinc_indexC], frac);

      sinc_sum += sinc_val * table[spaced_point];
    }
    return sinc_sum;
  }
} // namespace SincExtras

namespace Extras {
    SawOSNext::SawOSNext()
    {
    }
  
    SawOSNext::~SawOSNext() {}
  
    float SawOSNext::next(float freq, float phaseIn, double m_freqMul)
    {
      float phaseDiff = (phaseIn - m_lastPhase);
      m_lastPhase = phaseIn;
  
      m_phase += (phaseDiff);
      m_phase += freq * m_freqMul;
      if (m_phase >= 4.f)
        m_phase -= 8.f;
      else if (m_phase <= -4.f)
        m_phase += 8.f;
  
      float out = sc_wrap(m_phase, -1.0, 1.0);
      return out;
    }
  
    void SawOSNext::reset(double phaseIn)
    {
      m_lastPhase = phaseIn;
      m_phase = phaseIn;
    }
  
    ProcessFuncs::ProcessFuncs()
    {}
  
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
  
    float ProcessFuncs::get_out_no_interp(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc) {
    
      buf_loc = buf_loc*(float)(buf_divs-1);
      int ibuf_loc = (int)buf_loc;
  
      float findex = (ramp*fmaxindex);
      int index = (int)findex;
      float frac = findex - index;
      int ibuf_divs = (int)buf_divs;
      float out;
  
      float frac_loc = buf_loc - ibuf_loc;
      
      if(ibuf_loc == ibuf_divs-1)
      {
        int zero_index = (ibuf_loc*each_table_size);
        int loc = (index+zero_index)*num_chans + chan_loc;
  
        if (index==each_table_size-1)
          out = buf_data[loc]*(1.0f-frac) + buf_data[zero_index*num_chans]*frac;
        else
          out = buf_data[loc]*(1.0f-frac) + buf_data[loc+num_chans]*frac;
      } else {
  
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

    float ProcessFuncs::get_out_quadratic(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc, int wrap_clip) {
      
      buf_loc = buf_loc*(float)(buf_divs-1);
      int ibuf_loc = (int)buf_loc;
  
      float findex = (ramp*fmaxindex);
      int index = (int)findex;
      float frac = findex - index;
      int ibuf_divs = (int)buf_divs;
      float out;
  
      float frac_loc = buf_loc - ibuf_loc;
      

      int zero_index = (ibuf_loc*each_table_size);
      int loc0 = (index+zero_index)*num_chans + chan_loc;
      int loc1 = 0;
      int loc2 = 0;
      if(wrap_clip==1)
      {
        loc1 = (sc_clip(index+1, 0, each_table_size-1)+zero_index)*num_chans + chan_loc;
        loc2 = (sc_clip(index+2, 0, each_table_size-1)+zero_index)*num_chans + chan_loc;
      } else {
        loc1 = (sc_wrap(index+1, 0, each_table_size-1)+zero_index)*num_chans + chan_loc;
        loc2 = (sc_wrap(index+2, 0, each_table_size-1)+zero_index)*num_chans + chan_loc;
  
      }

      float out1 = quadratic_interpolation(buf_data[loc0], buf_data[loc1], buf_data[loc2], frac);

      float out2 = 0.f;
      if(ibuf_loc < ibuf_divs-1)
      {
        int zero_index2 = ((ibuf_loc+1)*each_table_size);
        int loc0 = (index+zero_index2)*num_chans + chan_loc;
        int loc1 = 0;
        int loc2 = 0;
        if(wrap_clip==1)
        {
          loc1 = (sc_clip(index+1, 0, each_table_size-1)+zero_index2)*num_chans + chan_loc;
          loc2 = (sc_clip(index+2, 0, each_table_size-1)+zero_index2)*num_chans + chan_loc;
        } else {
          loc1 = (sc_wrap(index+1, 0, each_table_size-1)+zero_index2)*num_chans + chan_loc;
          loc2 = (sc_wrap(index+2, 0, each_table_size-1)+zero_index2)*num_chans + chan_loc;
    
        }

        out2 = quadratic_interpolation(buf_data[loc0], buf_data[loc1], buf_data[loc2], frac);
      }
    
      out = out1*(1.0f-frac_loc) + out2*frac_loc;
    return out;
  
    }
  } // end namespace Extras
