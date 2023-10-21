/*************************************************************************

    This project implements a complete(!) JPEG (10918-1 ITU.T-81) codec,
    plus a library that can be used to encode and decode JPEG streams. 
    It also implements ISO/IEC 18477 aka JPEG XT which is an extension
    towards intermediate, high-dynamic-range lossy and lossless coding
    of JPEG. In specific, it supports ISO/IEC 18477-3/-6/-7/-8 encoding.

    Copyright (C) 2012-2017 Thomas Richter, University of Stuttgart and
    Accusoft.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/
/*
**
** This class represents one row of quantized data of coefficients, i.e. one
** row of 8x8 blocks.
**
** $Id: blockrow.hpp,v 1.9 2014/09/30 08:33:16 thor Exp $
**
*/

#ifndef CODING_BLOCKROW_HPP
#define CODING_BLOCKROW_HPP

/// Includes
#include "tools/environment.hpp"
///

/// class BlockRow
// This class represents one row of coefficients, i.e. one
// row of 8x8 blocks.
template<class T>
class BlockRow : public JKeeper {
  //
public:
  // the actual block. We always keep 32 bit data. In the future, this class
  // might be templated.
  struct Block  {
    T m_Data[64];
  };
  //
private:
  //
  // The block array itself.
  struct Block       *m_pBlocks;
  //
  // The extend in number of blocks.
  ULONG               m_ulWidth;
  //
  // The next row in a row stack.
  class QuantizedRow *m_pNext;
  //
public:
  BlockRow(class Environ *env);
  //
  ~BlockRow(void);
  //
  // Allocate a row of data, sufficient to hold the indicated number of cofficients. Note that
  // it is still up to the caller to include the subsampling factors.
  void AllocateRow(ULONG coefficients);
  //
  // Return the n'th block.
  struct Block *BlockAt(ULONG pos) const
  {
    assert(pos < m_ulWidth);
    return m_pBlocks + pos;
  }
  //
  // Return the next row.
  class QuantizedRow *NextOf(void) const
  {
    return m_pNext;
  }
  //
  // Return as lvalue.
  class QuantizedRow *&NextOf(void)
  {
    return m_pNext;
  }
  //
  // Width of this row in blocks.
  ULONG WidthOf(void) const
  {
    return m_ulWidth;
  }
  //
  // Tag a row on this row such that the passed argument is below the current row.
  void TagOn(class QuantizedRow *below)
  {
    m_pNext = below;
  }
};

///
#endif
