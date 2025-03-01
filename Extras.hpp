#include <iostream>

namespace SincExtras {
  class SincFuncs {
    public:
      SincFuncs(int ripples=8);
      ~SincFuncs();
      float get_spaced_out(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int spacing, float sinc_crossfade, int num_chans, int chan_loc);
      
      float get_spaced_sinc_sum(const float* table, int table_size, float index, int spacing, int ibuf_loc, int num_chans, int chan_loc);

      std::array<int, 8> sinc_points8;
      std::array<int, 32> sinc_points32;
      int* sinc_points;
      int sinc_len{8};
      int sinc_half_len{4};
      int sinc_table_size{8192};
      float fmax_sinc_offset{1024.f};
      int max_sinc_offset{1024};
      double* m_sinc_table;
      std::array<double, 8192>m_sinc_table8;
      std::array<double, 32768>m_sinc_table32;
      int counter{0};
      
    
    private:

  };
}


namespace Extras {

    class SawOSNext {
      public:
        SawOSNext();
        ~SawOSNext();
        float m_lastPhase {0.f};
        double m_phase {0.0};
        float next(float freq, float phaseIn, double m_freqMul);
        void reset(double phaseIn);
      private:
  
    };
  
    class ProcessFuncs {
      public:
        ProcessFuncs();
        ~ProcessFuncs();
        float get_phase(const float* phase_buf_data, float phase, float phase_buf_divs, float phase_buf_loc, int phase_table_size, float phase_fmaxindex);
        float get_out(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc);

        template<typename T>
        float get_out_quadratic(const T* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc, int wrap_clip);

        float get_out_no_interp(const float* buf_data, float ramp, float buf_divs, float buf_loc, int each_table_size, float fmaxindex, int num_chans, int chan_loc);
        float get_sinc_sum(const float* table, int table_size, int index, int ibuf_loc, int sinc_offset, int num_chans, int chan_loc);

        // std::array<int, 8> sinc_points;
        // int sinc_len{8};
        // int sinc_half_len{4};
        // int sinc_table_size{8192};
        //float fmax_sinc_offset{1024.f};
        // int max_sinc_offset{1024};
        // std::array<double, 8192>m_sinc_table;
      
      private:
  
    };
  }