// -*- mode: c++ -*-


/*---------------------------------------------------------------------------
  File: SHA1_TPP.cpp

  Original Author: Dominik Reichl <Dominik.Reichl@tiscali.de>
  - Modified by Patrick Pedrioli
  
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



//#include "stdafx.h"

#include <stdio.h> // Needed for file access
#include <memory.h> // Needed for memset and memcpy
#include <string.h> // Needed for strcat and strcpy


#include "SHA1_TPP.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


SHA1::SHA1()
{
	Reset();
}

SHA1::~SHA1()
{
	Reset();
}


void SHA1::Reset()
{
	// SHA1 initialization constants
	m_state[0] = 0x67452301;
	m_state[1] = 0xEFCDAB89;
	m_state[2] = 0x98BADCFE;
	m_state[3] = 0x10325476;
	m_state[4] = 0xC3D2E1F0;

	m_count[0] = 0;
	m_count[1] = 0;
}

void SHA1::Transform(unsigned int state[5], unsigned char buffer[64])
{
	unsigned int a = 0, b = 0, c = 0, d = 0, e = 0;

	SHA1_WORKSPACE_BLOCK* block;
	static unsigned char workspace[64];
	block = (SHA1_WORKSPACE_BLOCK *)workspace;
	memcpy(block, buffer, 64);

	// Copy state[] to working vars
	a = state[0];
	b = state[1];
	c = state[2];
	d = state[3];
	e = state[4];

	// 4 rounds of 20 operations each. Loop unrolled.
	R0(a,b,c,d,e, 0); R0(e,a,b,c,d, 1); R0(d,e,a,b,c, 2); R0(c,d,e,a,b, 3);
	R0(b,c,d,e,a, 4); R0(a,b,c,d,e, 5); R0(e,a,b,c,d, 6); R0(d,e,a,b,c, 7);
	R0(c,d,e,a,b, 8); R0(b,c,d,e,a, 9); R0(a,b,c,d,e,10); R0(e,a,b,c,d,11);
	R0(d,e,a,b,c,12); R0(c,d,e,a,b,13); R0(b,c,d,e,a,14); R0(a,b,c,d,e,15);
	R1(e,a,b,c,d,16); R1(d,e,a,b,c,17); R1(c,d,e,a,b,18); R1(b,c,d,e,a,19);
	R2(a,b,c,d,e,20); R2(e,a,b,c,d,21); R2(d,e,a,b,c,22); R2(c,d,e,a,b,23);
	R2(b,c,d,e,a,24); R2(a,b,c,d,e,25); R2(e,a,b,c,d,26); R2(d,e,a,b,c,27);
	R2(c,d,e,a,b,28); R2(b,c,d,e,a,29); R2(a,b,c,d,e,30); R2(e,a,b,c,d,31);
	R2(d,e,a,b,c,32); R2(c,d,e,a,b,33); R2(b,c,d,e,a,34); R2(a,b,c,d,e,35);
	R2(e,a,b,c,d,36); R2(d,e,a,b,c,37); R2(c,d,e,a,b,38); R2(b,c,d,e,a,39);
	R3(a,b,c,d,e,40); R3(e,a,b,c,d,41); R3(d,e,a,b,c,42); R3(c,d,e,a,b,43);
	R3(b,c,d,e,a,44); R3(a,b,c,d,e,45); R3(e,a,b,c,d,46); R3(d,e,a,b,c,47);
	R3(c,d,e,a,b,48); R3(b,c,d,e,a,49); R3(a,b,c,d,e,50); R3(e,a,b,c,d,51);
	R3(d,e,a,b,c,52); R3(c,d,e,a,b,53); R3(b,c,d,e,a,54); R3(a,b,c,d,e,55);
	R3(e,a,b,c,d,56); R3(d,e,a,b,c,57); R3(c,d,e,a,b,58); R3(b,c,d,e,a,59);
	R4(a,b,c,d,e,60); R4(e,a,b,c,d,61); R4(d,e,a,b,c,62); R4(c,d,e,a,b,63);
	R4(b,c,d,e,a,64); R4(a,b,c,d,e,65); R4(e,a,b,c,d,66); R4(d,e,a,b,c,67);
	R4(c,d,e,a,b,68); R4(b,c,d,e,a,69); R4(a,b,c,d,e,70); R4(e,a,b,c,d,71);
	R4(d,e,a,b,c,72); R4(c,d,e,a,b,73); R4(b,c,d,e,a,74); R4(a,b,c,d,e,75);
	R4(e,a,b,c,d,76); R4(d,e,a,b,c,77); R4(c,d,e,a,b,78); R4(b,c,d,e,a,79);

	// Add the working vars back into state[]
	state[0] += a;
	state[1] += b;
	state[2] += c;
	state[3] += d;
	state[4] += e;

	// Wipe variables
	a = 0; b = 0; c = 0; d = 0; e = 0;
}

// Use this function to hash in binary data and strings
void SHA1::Update(unsigned char* data, unsigned int len)
{
	unsigned int i = 0, j = 0;

	j = (m_count[0] >> 3) & 63;

	if((m_count[0] += len << 3) < (len << 3)) m_count[1]++;

	m_count[1] += (len >> 29);

	if((j + len) > 63)
	{
		memcpy(&m_buffer[j], data, (i = 64 - j));
		Transform(m_state, m_buffer);

		for (; i+63 < len; i += 64)
		{
			Transform(m_state, &data[i]);
		}

		j = 0;
	}
	else i = 0;

	memcpy(&m_buffer[j], &data[i], len - i);
}

// Hash in file contents
bool SHA1::HashFile(char *szFileName)
{
	unsigned int ulFileSize = 0, ulRest = 0, ulBlocks = 0;
	unsigned int i = 0;
	unsigned char uData[MAX_FILE_READ_BUFFER];
	FILE *fIn = NULL;

	if((fIn = fopen(szFileName, "rb")) == NULL) return(false);


	fseek(fIn, 0, SEEK_END);
	ulFileSize = ftell(fIn);
	fseek(fIn, 0, SEEK_SET);

	ulRest = ulFileSize % MAX_FILE_READ_BUFFER;
	ulBlocks = ulFileSize / MAX_FILE_READ_BUFFER;

	for(i = 0; i < ulBlocks; i++)
	{
		size_t errCheckMe;	  
		errCheckMe = fread(uData, 1, MAX_FILE_READ_BUFFER, fIn);
		Update(uData, MAX_FILE_READ_BUFFER);
	}

	if(ulRest != 0)
	{
		size_t errCheckMe;	  
		errCheckMe = fread(uData, 1, ulRest, fIn);
		Update(uData, ulRest);
	}

	fclose(fIn);
	fIn = NULL;

	return(true);
}

void SHA1::Final()
{
	unsigned int i = 0, j = 0;
	unsigned char finalcount[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

	for (i = 0; i < 8; i++)
		finalcount[i] = (unsigned char)((m_count[(i >= 4 ? 0 : 1)]
			>> ((3 - (i & 3)) * 8) ) & 255); // Endian independent

	Update((unsigned char *)"\200", 1);

	while ((m_count[0] & 504) != 448)
		Update((unsigned char *)"\0", 1);

	Update(finalcount, 8); // Cause a SHA1Transform()

	for (i = 0; i < 20; i++)
	{
		m_digest[i] = (unsigned char)((m_state[i >> 2] >> ((3 - (i & 3)) * 8) ) & 255);
	}

	// Wipe variables for security reasons
	i = 0; j = 0;
	memset(m_buffer, 0, 64);
	memset(m_state, 0, 20);
	memset(m_count, 0, 8);
	memset(finalcount, 0, 8);

	Transform(m_state, m_buffer);
}

// Get the final hash as a pre-formatted string
void SHA1::ReportHash(char *szReport, unsigned char uReportType)
{
	unsigned char i = 0;
	char szTemp[4];

	if(uReportType == REPORT_HEX)
	{
		for(i = 0; i < 20; i++)
		{
			sprintf(szTemp, "%02x", m_digest[i]);
			strcat(szReport, szTemp);
		}
	}
	else if(uReportType == REPORT_DIGIT)
	{
		sprintf(szTemp, "%u", m_digest[0]);
		strcat(szReport, szTemp);

		for(i = 1; i < 20; i++)
		{
			sprintf(szTemp, " %u", m_digest[i]);
			strcat(szReport, szTemp);
		}
	}
	else strcpy(szReport, "Error: Unknown report type!");
}

// Get the raw message digest
void SHA1::GetHash(unsigned char *uDest)
{
	memcpy(uDest, m_digest, 20);
}
