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
** This header provides the main function.
** This main function is only a demo, it is not part of the libjpeg code.
** It is here to serve as an entry point for the command line image
** compressor.
**
** $Id: main.hpp,v 1.9 2015/03/11 16:02:42 thor Exp $
**
*/

#ifndef CMD_MAIN_HPP
#define CMD_MAIN_HPP

/// Includes
#include "interface/types.hpp"
///

/// Prototypes
extern int main(int argc,char **argv);
extern void ParseSubsamplingFactors(UBYTE *sx,UBYTE *sy,const char *sub,int cnt);
///

///
#endif
