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
** This class contains and maintains the AC conditioning
** parameters.
**
** $Id: actemplate.hpp,v 1.8 2014/09/30 08:33:16 thor Exp $
**
*/

#ifndef CODING_ACTEMPLATE_HPP
#define CODING_ACTEMPLATE_HPP

/// Includes
#include "tools/environment.hpp"
///

/// Forwards
class ByteStream;
///

/// ACTemplate
// This class contains and maintains the AC conditioning
// parameters.
#if ACCUSOFT_CODE
class ACTemplate : public JKeeper {
  //
  // The lower threshold, also known as 'L' parameter in the specs.
  UBYTE m_ucLower;
  //
  // The upper threshold parameter, also known as 'U' parameter.
  UBYTE m_ucUpper;
  //
  // The block index that discriminates between lower and upper
  // block indices for AC coding.
  UBYTE m_ucBlockEnd;
  //
public:
  ACTemplate(class Environ *env);
  //
  ~ACTemplate(void);
  //
  // Parse off DC conditioning parameters.
  void ParseDCMarker(class ByteStream *io);
  //
  // Parse off an AC conditioning parameter.
  void ParseACMarker(class ByteStream *io);
  //
  // Just install the defaults found in the standard.
  void InitDefaults(void);
  //
  // Return the largest block index that still counts
  // as a lower index for AC coding. This is the kx parameter
  UBYTE BandDiscriminatorOf(void) const
  {
    return m_ucBlockEnd;
  }
  //
  // Return the L parameter.
  UBYTE LowerThresholdOf(void) const
  {
    return m_ucLower;
  }
  //
  // Return the U parameter.
  UBYTE UpperThresholdOf(void) const
  {
    return m_ucUpper;
  }
};
#endif
///

///
#endif
