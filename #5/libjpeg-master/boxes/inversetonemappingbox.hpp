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
** This box keeps an inverse tone mapping curve, as required for the
** R and L transformations.
**
** $Id: inversetonemappingbox.hpp,v 1.16 2014/09/30 08:33:15 thor Exp $
**
*/

#ifndef BOXES_INVERSETONEMAPPINGBOX_HPP
#define BOXES_INVERSETONEMAPPINGBOX_HPP

/// Includes
#include "boxes/box.hpp"
#include "boxes/tonemapperbox.hpp"
///

/// class InverseToneMappingBox
class InverseToneMappingBox : public ToneMapperBox {
  //
  // The table itself, indexed by the decoded sample value.
  LONG  *m_plTable;
  // 
  // Inverse (encoding) tone mapping curve, if available.
  LONG  *m_plInverseMapping;
  // 
  // Number of additional residual bits bits. Actually, this is only used to
  // test whether the table fits to the data (since it is always
  // unscaled). This is called R_d in the standard.
  UBYTE  m_ucResidualBits;
  //
  // Second level parsing stage: This is called from the first level
  // parser as soon as the data is complete. Must be implemented
  // by the concrete box.
  virtual bool ParseBoxContent(class ByteStream *stream,UQUAD boxsize);
  //
  // Second level creation stage: Write the box content into a temporary stream
  // from which the application markers can be created.
  virtual bool CreateBoxContent(class MemoryStream *target);
  //
public:
  enum {
    Type = MAKE_ID('T','O','N','E')
  };
  //
  InverseToneMappingBox(class Environ *env,class Box *&boxlist)
    : ToneMapperBox(env,boxlist,Type), m_plTable(NULL), m_plInverseMapping(NULL)
  { }
  //
  virtual ~InverseToneMappingBox(void);
  //
  //
  // Return the size of the table in entries.
  ULONG EntriesOf(void) const
  {
    return m_ulTableEntries;
  }
  //
  // Return the table itself.
  const LONG *TableOf(void) const
  {
    return m_plTable;
  }
  //
  // Define the table from an external source.
  void DefineTable(UBYTE tableidx,const UWORD *table,ULONG size,UBYTE residualbits);
  //
  // Check whether the given table is identical to the table stored here, and thus
  // no new index is required (to save rate). Returns true if the two are equal.
  bool CompareTable(const UWORD *table,ULONG size,UBYTE residualbits) const;
  //
  // Return a table that maps inputs in the range 0..2^inputbits-1
  // to output bits in the range 0..2^oututbits-1.
  virtual const LONG *ScaledTableOf(UBYTE inputbits,UBYTE outputbits,UBYTE infract,UBYTE outfract);
  //
  // This is the floating point version of the above. It returns floating point sample
  // values instead of integer sample values.
  virtual const FLOAT *FloatTableOf(UBYTE,UBYTE,UBYTE,UBYTE)
  {
    // Currently, no floating point workflow for the integer tables.
    return NULL;
  }
  //
  // Return the inverse of the table, where the first argument is the number
  // of bits in the DCT domain (the output bits) and the second argument is
  // the number of bits in the spatial (image) domain, i.e. the argument
  // order is identical to that of the backwards table generated above.
  virtual const LONG *InverseScaledTableOf(UBYTE dctbits,UBYTE spatialbits,UBYTE dctfract,UBYTE spatialfract);
};
///

///
#endif
