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
** This class contains and maintains the huffman code parsers.
**
** $Id: huffmantable.cpp,v 1.19 2015/09/17 11:20:35 thor Exp $
**
*/

/// Includes
#include "marker/huffmantable.hpp"
#include "io/bytestream.hpp"
#include "coding/huffmantemplate.hpp"
#include "coding/arithmetictemplate.hpp"
#include "std/string.hpp"
///

/// HuffmanTable::HuffmanTable
HuffmanTable::HuffmanTable(class Environ *env)
  : JKeeper(env)
{
  memset(m_pCoder,0,sizeof(m_pCoder));
}
///

/// HuffmanTable::~HuffmanTable
HuffmanTable::~HuffmanTable(void)
{
  int i;

  for(i = 0;i < 8;i++) {
    delete m_pCoder[i];
  }
}
///


/// HuffmanTable::isEmpty
// Check whether the tables are empty or not
// In such a case, do not write the tables.
bool HuffmanTable::isEmpty(void) const
{ 
  int i = 0;

  for(i = 0;i < 8;i++) {
    if (m_pCoder[i]) {
      return false;
    }
  }

  return true;
}
///

/// HuffmanTable::WriteMarker
// Write the currently defined huffman tables back to a stream.
void HuffmanTable::WriteMarker(class ByteStream *io)
{
  int i = 0;
  ULONG len = 2; // marker size itself.

  for(i = 0;i < 8;i++) {
    if (m_pCoder[i]) {
      len += 1; // Th and Tc field.
      len += m_pCoder[i]->MarkerOverhead();
    }
  }

  if (len > MAX_UWORD)
    JPG_THROW(OVERFLOW_PARAMETER,"HuffmanTable::WriteMarker","DHT marker overhead too large, Huffman tables too complex");

  io->PutWord(len);
  
  for(i = 0;i < 8;i++) {
    if (m_pCoder[i]) {
      UBYTE type = 0;
      if (i >= 4)
        type |= 0x10;    // is an AC table then.
      type |= i & 0x03;  // the huffman table identifier.
      io->Put(type);
      m_pCoder[i]->WriteMarker(io);
    }
  }
}
///  

/// HuffmanTable::ParseMarker
// Parse the marker contents of a DHT marker.
void HuffmanTable::ParseMarker(class ByteStream *io)
{
  LONG len = io->GetWord();

  if (len < 2)
    JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker","Huffman table length must be at least two bytes long");

  len -= 2; // remove the marker length.

  while(len > 0) {
    LONG  t = io->Get();
    UQUAD p = io->FilePosition();
    UQUAD q;
    
    if (t == ByteStream::EOF)
      JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker","Huffman table marker run out of data");
    len--;
    
    if ((t >> 4) > 1) {
      JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker","undefined Huffman table type");
      return;
    }
    if ((t & 0x0f) > 3) {
      JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker",
                "invalid Huffman table destination, must be between 0 and 3");
    }
    t = (t & 0x03) | ((t & 0xf0) >> 2);
    delete m_pCoder[t];m_pCoder[t] = NULL;
    m_pCoder[t] = new(m_pEnviron) HuffmanTemplate(m_pEnviron);
    m_pCoder[t]->ParseMarker(io);
    
    q = io->FilePosition();
    assert(q >= p);
    q -= p;

    if (q > ULONG(len))
      JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker","huffman table size corrupt");
    len -= q;
  }

  if (len)
    JPG_THROW(MALFORMED_STREAM,"HuffmanTable::ParseMarker","huffman table size is corrupt");
}
///

/// HuffmanTable::AdjustToStatistics
// Adjust all coders in here to the statistics collected before, i.e.
// find optimal codes.
void HuffmanTable::AdjustToStatistics(void)
{
  for(int i = 0;i < 8;i++) {
    if (m_pCoder[i])
      m_pCoder[i]->AdjustToStatistics();
  }
}
///

/// HuffmanTable::DCTemplateOf
// Get the template for the indicated DC table or NULL if it doesn't exist.
class HuffmanTemplate *HuffmanTable::DCTemplateOf(UBYTE idx,ScanType type,UBYTE depth,UBYTE hidden)
{
  assert(m_pCoder && idx < 4);
  
  if (m_pCoder[idx] == NULL) {
    m_pCoder[idx] = new(m_pEnviron) class HuffmanTemplate(m_pEnviron);
    // Provide a default that seems sensible. Everything else requires
    // measurement.
    if (idx == 0) {
      m_pCoder[idx]->InitDCLuminanceDefault(type,depth,hidden);
    } else {
      m_pCoder[idx]->InitDCChrominanceDefault(type,depth,hidden);
    }
  }
  
  return m_pCoder[idx];
}
///

/// HuffmanTable::ACTemplateOf
// Get the template for the indicated AC table or NULL if it doesn't exist.
class HuffmanTemplate *HuffmanTable::ACTemplateOf(UBYTE idx,ScanType type,UBYTE depth,UBYTE hidden)
{
  assert(m_pCoder && idx < 4);

  idx += 4;
  
  if (m_pCoder[idx] == NULL) {
    m_pCoder[idx] = new(m_pEnviron) class HuffmanTemplate(m_pEnviron);
    // Provide a default that seems sensible. Everything else requires
    // measurement.
    if (idx == 4) {
      m_pCoder[idx]->InitACLuminanceDefault(type,depth,hidden);
    } else {
      m_pCoder[idx]->InitACChrominanceDefault(type,depth,hidden);
    }
  }
  
  return m_pCoder[idx];
}
///
