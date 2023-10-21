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
** This class implements a simple deringing filter to avoid
** DCT artifacts (Gibbs Phenomenon).
**
** $Id: deringing.hpp,v 1.2 2016/10/28 13:58:54 thor Exp $
**
*/

#ifndef DCT_DERINGING_HPP
#define DCT_DERINGING_HPP

/// Includes
#include "tools/environment.hpp"
///

/// Forwards
class Frame;
class DCT;
///

/// class DeRinger
// This implements a deringing filter on top of a DCT
class DeRinger : public JKeeper {
  //
  // The DCT for the transformation.
  class DCT *m_pDCT;
  //
  // Maximum and minimum values for the current frame.
  LONG m_lMin;
  LONG m_lMax;
  LONG m_lDelta;
  //
#if ACCUSOFT_CODE
  // Run a simple Gaussian smoothing filter on the second argument, place
  // the result in the first argument, or copy from the original if mask is
  // not set.
  void Smooth(LONG target[64],const LONG src[64],const LONG mask[64]);
#endif  
  //
  //
public:
  DeRinger(class Frame *frame,class DCT *dct);
  //
  ~DeRinger(void);
  //
  // Remove Gibbs' pheonomenon artifacts from the given image block 
  // (non-DCT-transformed) by including overshooting
  // in the extreme image parts, or undershooting in the dark image regions.
  void DeRing(const LONG block[64],LONG dst[64],LONG dcshift);
  //
};
///

///
#endif
