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
** This class provides the necessary mechanisms for superbox parsing.
** It describes a box that contains sub-boxes and parses the sub-boxes.
**
** $Id: superbox.hpp,v 1.5 2014/09/30 08:33:15 thor Exp $
**
*/

#ifndef BOXES_SUPERBOX_HPP
#define BOXES_SUPERBOX_HPP

/// Includes
#include "interface/types.hpp"
#include "tools/environment.hpp"
#include "boxes/box.hpp"
///

/// Forwards
class NameSpace;
///

/// class SuperBox
// This is the base class for all superboxes.
// It decodes its contents as sub-boxes.
class SuperBox : public Box {
  //
  // The list of sub-boxes of this superbox. This is a singly linked list
  // of boxes queued by their m_pNext pointer.
  class Box *m_pSubBoxes;
  //
  // Parse the contents of the superbox as sub-boxes. This is implemented here
  // because the contents is already structured at the box-level. It
  // creates boxes from the box types, but leaves the actual box-parsing to the
  // correspondig super-box implementation.
  virtual bool ParseBoxContent(class ByteStream *stream,UQUAD boxsize);
  //
  // Write the super box content, namely all the boxes, into the output stream.
  // This also calls the superbox implementation to perform the operation.
  virtual bool CreateBoxContent(class MemoryStream *target);
  //
  // To be implemented by derived classes: Create a box for this superbox
  // in the given box list.
  virtual class Box *CreateBox(class Box *&boxlist,ULONG tbox) = 0;
  // 
  // Inform the superbox that the box is now created. Does nothing by default,
  // can be overloaded to sort the new box in.
  virtual void AcknowledgeBox(class Box *,ULONG)
  {
    // by default, nothing happens here.
  }
  //
protected:
  //
  // Register this box as primary namespace.
  void RegisterNameSpace(class NameSpace *names);
  //
public:
  //
  SuperBox(class Environ *env,class Box *&boxlist,ULONG boxtype);
  //
  virtual ~SuperBox(void);
  //
  // Create a new box as sub-box of this superbox.
  class Box *CreateBox(ULONG type);
};
///

///
#endif
