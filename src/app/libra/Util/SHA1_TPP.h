// -*- mode: c++ -*-


/*---------------------------------------------------------------------------
  File: SHA1.h

  Original Author: Dominik Reichl <Dominik.Reichl@tiscali.de>

  Description:
    A generic sha1 hash.
---------------------------------------------------------------------------*/


/*****************************************************************************

   This program is free software; you can redistribute it and/or modify  
   it under the terms of the GNU Library or "Lesser" General Public      
   License (LGPL) as published by the Free Software Foundation;          
   either version 2 of the License, or (at your option) any later        
   version.                                                              

*****************************************************************************/



/*
	100% free public domain implementation of the SHA-1
	algorithm by Dominik Reichl <Dominik.Reichl@tiscali.de>


	=== Test Vectors (from FIPS PUB 180-1) ===

	"abc"
		A9993E36 4706816A BA3E2571 7850C26C 9CD0D89D

	"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
		84983E44 1C3BD26E BAAE4AA1 F95129E5 E54670F1

	A million repetitions of "a"
		34AA973C D4C4DAA4 F61EEB2B DBAD2731 6534016F
*/


#ifndef _INCLUDED_SHA1_H_
#define _INCLUDED_SHA1_H_


#define MAX_FILE_READ_BUFFER 8000

class SHA1 {
 public:
	// Rotate x bits to the left
	#define ROL32(value, bits) (((value)<<(bits))|((value)>>(32-(bits))))
  
	#ifdef LITTLE_ENDIAN
	#define SHABLK0(i) (block->l[i] = (ROL32(block->l[i],24) & 0xFF00FF00)	\
			    | (ROL32(block->l[i],8) & 0x00FF00FF))
	#else // (big endian)
	#define SHABLK0(i) (block->l[i])
	#endif // LITTLE_ENDIAN

	#define SHABLK(i) (block->l[i&15] = ROL32(block->l[(i+13)&15] ^ block->l[(i+8)&15] \
						  ^ block->l[(i+2)&15] ^ block->l[i&15],1))

	// SHA-1 rounds
	#define R0(v,w,x,y,z,i) { z+=((w&(x^y))^y)+SHABLK0(i)+0x5A827999+ROL32(v,5); w=ROL32(w,30); }
	#define R1(v,w,x,y,z,i) { z+=((w&(x^y))^y)+SHABLK(i)+0x5A827999+ROL32(v,5); w=ROL32(w,30); }
	#define R2(v,w,x,y,z,i) { z+=(w^x^y)+SHABLK(i)+0x6ED9EBA1+ROL32(v,5); w=ROL32(w,30); }
	#define R3(v,w,x,y,z,i) { z+=(((w|x)&y)|(w&x))+SHABLK(i)+0x8F1BBCDC+ROL32(v,5); w=ROL32(w,30); }
	#define R4(v,w,x,y,z,i) { z+=(w^x^y)+SHABLK(i)+0xCA62C1D6+ROL32(v,5); w=ROL32(w,30); }

	typedef union {
		unsigned char c[64];
		unsigned int l[16];
	} SHA1_WORKSPACE_BLOCK;

	// Two different formats for ReportHash(...)
	enum { REPORT_HEX = 0, REPORT_DIGIT = 1 };

	// Constructor and Destructor
	SHA1();
	virtual ~SHA1();

	unsigned int m_state[5];
	unsigned int m_count[2];
	unsigned char m_buffer[64];
	unsigned char m_digest[20];

	void Reset();

	// Update the hash value
	void Update(unsigned char* data, unsigned int len);
	bool HashFile(char *szFileName);

	// Finalize hash and report
	void Final();
	void ReportHash(char *szReport, unsigned char uReportType = REPORT_HEX);
	void GetHash(unsigned char *uDest);

private:
	// Private SHA-1 transformation
	void Transform(unsigned int state[5], unsigned char buffer[64]);
};

#endif // header guards
