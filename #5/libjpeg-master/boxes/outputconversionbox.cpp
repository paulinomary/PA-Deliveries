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
** This box defines the output process, namely whether data is casted 
** to float, whether a nonlinearity is applied, whether the data is
** clamped and several other options required as last processing
** step of the output conversion.
**
** It is a subbox of the merging specification box.
**
** $Id: outputconversionbox.cpp,v 1.3 2014/09/30 08:33:15 thor Exp $
**
*/

/// Includes
#include "boxes/box.hpp"
#include "boxes/outputconversionbox.hpp"
#include "io/bytestream.hpp"
#include "io/memorystream.hpp"
///

/// OutputConversionBox::CreateBoxContent
// Create the contents of the box. The size
// is computed automatically by the superbox.
bool OutputConversionBox::CreateBoxContent(class MemoryStream *target)
{
  UBYTE v;

  v  = m_ucExtraRangeBits << 4;
  if (m_bLossless)
    v |= 0x08;
  if (m_bCastToFloat)
    v |= 0x04;
  if (m_bEnableClamping)
    v |= 0x02;
  if (m_bEnableLookup)
    v |= 0x01;

  target->Put(v);

  if (m_bEnableLookup) {
    target->Put((m_ucOutputLookup[0] << 4) | (m_ucOutputLookup[1]));
    target->Put((m_ucOutputLookup[2] << 4) | (m_ucOutputLookup[3]));
  } else {
    target->Put(0);
    target->Put(0);
  }

  return true;
}
///

/// OutputConversionBox::ParseBoxContent
// Read the contents of this box and place it in the fields.
bool OutputConversionBox::ParseBoxContent(class ByteStream *stream,UQUAD boxsize)
{
  UBYTE v;

  if (boxsize != 3)
    JPG_THROW(MALFORMED_STREAM,"OutputConversionBox::ParseBoxContent","Malformed JPEG stream, "
              "Output Conversion box size is invalid");

  v = stream->Get();

  m_ucExtraRangeBits = v >> 4;
  if (m_ucExtraRangeBits > 8)
    JPG_THROW(MALFORMED_STREAM,"OutputConversionBox::ParseBoxContent","Malformed JPEG stream, "
              "bit depths cannot be larger than 16");

  m_bLossless       = (v & 0x08)?true:false;
  m_bCastToFloat    = (v & 0x04)?true:false;
  m_bEnableClamping = (v & 0x02)?true:false;
  m_bEnableLookup   = (v & 0x01)?true:false;

  if (m_bEnableLookup) {
    v = stream->Get();
    m_ucOutputLookup[0] = v >> 4;
    m_ucOutputLookup[1] = v & 0x0f;
    v = stream->Get();
    m_ucOutputLookup[2] = v >> 4;
    m_ucOutputLookup[3] = v & 0x0f;
  } else {
    if (stream->GetWord() != 0)
      JPG_THROW(MALFORMED_STREAM,"OutputConversionBox::ParseBoxContent","Malformed JPEG stream, "
                "output conversion is disabled, but lookup information is not zero");
  }

  return true;
}
///
