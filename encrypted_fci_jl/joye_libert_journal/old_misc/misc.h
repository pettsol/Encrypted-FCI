#ifndef MISC_H
#define MISC_H

#include <gmp.h>
#include <NTL/GF2X.h>

void MyBytesFromGF2X(unsigned char* buffer, NTL::GF2X& p, int numbytes);
void MPZToZZ(NTL::ZZ *out, mpz_t in);
void ZZToMPZ(mpz_t out, NTL::ZZ *in);

#endif
