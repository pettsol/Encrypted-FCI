#include "misc.h"

#include <iostream>


int main()
{
	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_t q, recov_q;
	mpz_init(q);
	mpz_init(recov_q);

	mpz_urandomb(q, state, 1024);

	gmp_printf("GMP representation: %Zd\n", q);

	NTL::ZZ p;

	MPZToZZ(&p, q);
	
	std::cout << "NTL representation: " << p << std::endl;

	ZZToMPZ(recov_q, &p);

	gmp_printf("GMP representation: %Zd\n", recov_q);

}
