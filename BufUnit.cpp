#include "BufUnit.hpp"
#include <limits>

namespace BufUnit {
  BufUnit::BufUnit()
  {
    m_fbufnum = std::numeric_limits<float>::quiet_NaN();
    m_buf = nullptr;
    m_buf_failed = false;
  }

  BufUnit::~BufUnit() {}

  bool BufUnit::GetTable(World* world, float fbufnum, int inNumSamples, const SndBuf*& buf, const float*& bufData,
                               int& tableSize) {
    
    if (fbufnum < 0.f) {                                                                      
        fbufnum = 0.f;                                                                                   
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
        return false;
    }

    bufData = buf->data;
    
    if (!bufData) {
      return false;
    }
    tableSize = buf->samples;
    return true;
}
} // namespace BufUnit