#include "misc.h"

#include <sstream>
// Map a NTL::GF2X polynomial representation of an element in GF(2^n)
// to an integer representation. Numbytes is 4 for unsigned int or 8 for
// unsigned long.
void MyBytesFromGF2X(unsigned char* buffer, NTL::GF2X& p, int numbytes) {
    int degree = NTL::deg(p);
    memset(buffer,0,numbytes);
    for(int i=0; i<=degree; i++) {
        if(NTL::IsOne(NTL::coeff(p,i))) {
            buffer[i/8] |= 1 << i%8;
        }
    }
}

void MPZToZZ(NTL::ZZ *out, mpz_t in)
{
	char *a = mpz_get_str(NULL, 10, in);

	*out = NTL::conv<NTL::ZZ>(a);
}

void ZZToMPZ(mpz_t out, NTL::ZZ *in)
{
	std::stringstream ss;
	ss << *in;

	mpz_set_str(out, ss.str().c_str(), 10);
}
