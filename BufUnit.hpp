#include "SC_PlugIn.hpp"

namespace BufUnit {
  class BufUnit : public SCUnit {
    public:

    BufUnit();

    // Destructor
    ~BufUnit();
      float m_fbufnum;
      bool m_buf_failed{false};
      SndBuf* m_buf;
    
      bool GetTable(World* world, float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData, int& tableSize);
    private:
};
} // namespace BufUnit