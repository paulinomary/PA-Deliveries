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
** This class represents multiple boxes that all contain a color
** transformation specification, or rather, a single index to a matrix.
**
** $Id: colortrafobox.hpp,v 1.4 2014/09/30 08:33:14 thor Exp $
**
*/

#ifndef BOXES_COLORTRAFOBOX_HPP
#define BOXES_COLORTRAFOBOX_HPP

/// Includes
#include "boxes/box.hpp"
#include "std/string.hpp"
///

/// class ColorTrafoBox
// This class represents multiple boxes that all contain a color
// transformation specification, or rather, a single index to a matrix.
class ColorTrafoBox : public Box {
  //
  // Index of the color transformation matrix to use.
  UBYTE m_ucTrafoIndex;
  //
  // Second level parsing stage: This is called from the first level
  // parser as soon as the data is complete. Must be implemented
  // by the concrete box. Returns true in case the contents is
  // parsed and the stream can go away.
  virtual bool ParseBoxContent(class ByteStream *stream,UQUAD boxsize);
  //
  // Second level creation stage: Write the box content into a temporary stream
  // from which the application markers can be created.
  // Returns whether the box content is already complete and the stream
  // can go away.
  virtual bool CreateBoxContent(class MemoryStream *target);
  //
public:
  enum {
    Base_Type          = MAKE_ID('L','T','R','F'), // base transformation
    Color_Type         = MAKE_ID('C','T','R','F'), // color transformation
    Residual_Type      = MAKE_ID('R','T','R','F'), // residual transformation
    ResidualColor_Type = MAKE_ID('D','T','R','F'), // residual color transformation
    Prescaling_Type    = MAKE_ID('S','T','R','F')  // prescaling transformation
  };
  //
  // Create a color transformation box. This also requires a type since there are
  // multiple boxes that all share the same syntax.
  ColorTrafoBox(class Environ *env,class Box *&boxlist,ULONG type)
    : Box(env,boxlist,type), m_ucTrafoIndex(0)
  { 
  }
  //
  virtual ~ColorTrafoBox(void)
  { }
  //
  // Return the index of the color transformation to be used.
  UBYTE TransformationIndexOf(void) const
  {
    return m_ucTrafoIndex;
  }
  //
  // Define the transformation index of component comp.
  void DefineTransformationIndex(UBYTE idx)
  {
    assert(idx <= 15);

    m_ucTrafoIndex = idx;
  }
};
///

///
#endif

    
