// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-20.
//

#include "binaryio.h"
namespace bin_io
{

static const char b64str[]
   = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
     "abcdefghijklmnopqrstuvwxyz"
     "0123456789+/";

void WriteBase64(std::ostream &out, const void *bytes, size_t nbytes)
{
   const unsigned char *in = static_cast<const unsigned char *>(bytes);
   const unsigned char *end = in + nbytes;
   while (end - in >= 3)
   {
      out << b64str[in[0] >> 2];
      out << b64str[((in[0] & 0x03) << 4) | (in[1] >> 4)];
      out << b64str[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
      out << b64str[in[2] & 0x3f];
      in += 3;
   }
   if (end - in > 0) // Padding
   {
      out << b64str[in[0] >> 2];
      if (end - in == 1)
      {
         out << b64str[(in[0] & 0x03) << 4];
         out << '=';
      }
      else // end - in == 2
      {
         out << b64str[((in[0] & 0x03) << 4) | (in[1] >> 4)];
         out << b64str[(in[1] & 0x0f) << 2];
      }
      out << '=';
   }
}

static const unsigned char b64table[] =
{
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,62, 255,62, 255,63,
   52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 255,255,255,0,  255,255,
   255,0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
   15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 255,255,255,255,63,
   255,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
   41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255
};

void DecodeBase64(const char *src, size_t len, std::vector<char> &buf)
{
   const unsigned char *in = (const unsigned char *)src;
   buf.clear();
   size_t count = 0;
   for (size_t i=0; i<len; ++i) { if (b64table[in[i]] != 255) { ++count; } }
   if (count % 4 != 0) { return; }
   buf.resize(3*len/4);
   unsigned char *out = (unsigned char *)buf.data();
   count = 0;
   int pad = 0;
   unsigned char c[4];
   for (size_t i=0; i<len; ++i)
   {
      unsigned char t = b64table[in[i]];
      if (t == 255) { continue; }
      if (in[i] == '=') { ++pad; }
      c[count++] = t;
      if (count == 4)
      {
         *out++ = (c[0] << 2) | (c[1] >> 4);
         if (pad <= 1) { *out++ = (c[1] << 4) | (c[2] >> 2); }
         if (pad == 0) { *out++ = (c[2] << 6) | c[3]; }
         count = pad = 0;
      }
   }
   buf.resize(out - (unsigned char *)buf.data());
}

size_t NumBase64Chars(size_t nbytes) { return ((4*nbytes/3) + 3) & ~3; }

} // namespace mfem::bin_io