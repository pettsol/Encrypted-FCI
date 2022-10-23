#include "paillier.h"
#include <iostream>
#include <gmp.h>

int main()
{
	// We need a random state
	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_t summand1, summand2, output;
	mpz_init_set_ui(summand1, 10);
	mpz_init_set_ui(summand2, 20);
	mpz_init(output);

	// Define Paillier parameters
	mpz_t N, N2, lambda, mu, g;
	mpz_init(N);
	mpz_init(N2);
	mpz_init(lambda);
	mpz_init(mu);
	mpz_init(g);

	// Set keysize = 3072 bits
	uint32_t keysize = 3072;

	paillier_keygen(N, N2, g, lambda, mu, state, keysize);

	// Variables to hold encrypted data
	mpz_t c1, c2, c3;
	mpz_init(c1);
	mpz_init(c2);
	mpz_init(c3);

	paillier_encrypt(c1, state, summand1, g, N, N2);
	paillier_encrypt(c2, state, summand2, g, N, N2);

	mpz_mul(c3, c1, c2);
	mpz_mod(c3, c3, N2);

	paillier_decrypt(output, c3, lambda, mu, N, N2);

	// Print output to check for homomorphic property
	gmp_printf("Addition output = %Zd\n", output);

	// Check multiplication with plaintext constant
	mpz_t c4;
	mpz_init(c4);

	mpz_powm(c4, c3, summand1, N2);
	paillier_decrypt(output, c4, lambda, mu, N, N2);

	gmp_printf("Plaintext multiplication output = %Zd\n", output);
}
