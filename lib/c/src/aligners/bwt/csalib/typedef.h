/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/

#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_

#include <inttypes.h>

#ifndef i32
typedef int32_t i32;
#endif
#ifndef u32
typedef uint32_t u32;
#endif
#ifndef i64
typedef int64_t i64;
#endif
#ifndef u64
typedef uint64_t u64;
#endif
#ifndef i16
typedef int16_t i16;
#endif
#ifndef u16
typedef uint16_t u16;
#endif


#ifndef uchar
typedef unsigned char uchar;
#endif

#ifndef ushort
typedef unsigned short ushort;
#endif

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif

#ifndef byte
typedef uint8_t byte;
#endif
#ifndef word
typedef uint16_t word;
#endif
#ifndef dword
typedef uint32_t dword;
#endif
#ifndef qword
typedef uint64_t qword;
#endif

#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif


#endif // _TYPEDEF_H_
