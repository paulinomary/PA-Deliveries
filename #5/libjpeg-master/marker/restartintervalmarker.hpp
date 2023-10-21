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
** This class keeps the restart interval size in MCUs.
**
** $Id: restartintervalmarker.hpp,v 1.7 2014/09/30 08:33:17 thor Exp $
**
*/

#ifndef MARKER_RESTARTINTERVALMARKER_HPP
#define MARKER_RESTARTINTERVALMARKER_HPP

/// Includes
#include "tools/environment.hpp"
///

/// Forwards
class ByteStream;
///

/// RestartIntervalMarker
// This class keeps the restart interval size in MCUs.
class RestartIntervalMarker : public JKeeper {
  //
  // The restart interval size in MCUs, or zero if restart markers
  // are disabled.
  UWORD     m_usRestartInterval;
  //
public:
  RestartIntervalMarker(class Environ *env);
  //
  ~RestartIntervalMarker(void)
  {
  }
  //
  // Install the defaults, namely the interval.
  void InstallDefaults(UWORD inter)
  {
    m_usRestartInterval = inter;
  }
  //
  // Return the currently active restart interval.
  UWORD RestartIntervalOf(void) const
  {
    return m_usRestartInterval;
  }
  //
  // Write the marker (without the marker id) to the stream.
  void WriteMarker(class ByteStream *io) const;
  //
  // Parse the marker from the stream.
  void ParseMarker(class ByteStream *io);
  //
};
///

///
#endif
