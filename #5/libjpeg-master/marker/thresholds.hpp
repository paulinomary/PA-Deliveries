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
** This class contains the marker that defines the JPEG LS thresholds.
**
** $Id: thresholds.hpp,v 1.9 2014/09/30 08:33:17 thor Exp $
**
*/

#ifndef MARKER_THRESHOLDS_HPP
#define MARKER_THRESHOLDS_HPP

/// Includes
#include "tools/environment.hpp"
///

/// Forwards
class ByteStream;
///

/// Thresholds
class Thresholds : public JKeeper {
  //
  // Parameters of the JPEG LS coder as defined by C.2.4.1.1
  //
  // Maximum sample value.
  UWORD m_usMaxVal;
  //
  // Bucket threshold 1
  UWORD m_usT1;
  //
  // Bucket threshold 2
  UWORD m_usT2;
  //
  // Bucket threshold 3
  UWORD m_usT3;
  //
  // The statistics reset value.
  UWORD m_usReset;
  //
public:
  Thresholds(class Environ *env);
  //
  ~Thresholds(void);
  //
  // Write the marker contents to a LSE marker.
  void WriteMarker(class ByteStream *io);
  //
  // Parse the marker contents of a LSE marker.
  // marker length and ID are already parsed off.
  void ParseMarker(class ByteStream *io,UWORD len);
  //
  // Return the maximum sample value.
  UWORD MaxValOf(void) const
  {
    return m_usMaxVal;
  }
  //
  // Return the T1 value.
  UWORD T1Of(void) const
  {
    return m_usT1;
  }
  //
  // Return the T2 value.
  UWORD T2Of(void) const
  {
    return m_usT2;
  }
  //
  // Return the T3 value.
  UWORD T3Of(void) const
  {
    return m_usT3;
  }
  //
  // Return the Reset interval.
  UWORD ResetOf(void) const
  {
    return m_usReset;
  }
  //
  // Install the defaults for a given bits per pixel value
  // and the quality/near value.
  void InstallDefaults(UBYTE bpp,UWORD near);
};
///

///
#endif
