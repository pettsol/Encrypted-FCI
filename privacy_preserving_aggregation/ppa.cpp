/////////////////////////////////////////////////
//
//
// Petter Solnør
////////////////////////////////////////////////

#include "ppa.h"
#include <gmp.h>
#include <iostream>

void ppa_setup(mpz_t N, mpz_t sk[], const uint32_t keysize, const uint32_t n_sensors)
{
	// We need to generate primes p, q of approximately the same size
	// to generate a safe modulus whose number of bits is given by keysize
	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_t p, q;
	mpz_init(p);
	mpz_init(q);

	// Set p to a random value with keysize/2 bits
	mpz_urandomb(p, state, keysize/2);
#ifdef DEBUG
	gmp_printf("Preparing to find primes\n");
#endif
	mpz_nextprime(p, p);
#ifdef DEBUG
	gmp_printf("Value for prime p: %Zd\n", p);
#endif

	// Set q to a random value with keysize/2 bits
	mpz_urandomb(q, state, keysize/2);
	mpz_nextprime(q, q);
#ifdef DEBUG
	gmp_printf("Value for prime q: %Zd\n", q);
#endif

	// Compute the modulus N
	mpz_mul(N, p, q);

#ifdef DEBUG
	gmp_printf("Value for factoring modulus N: %Zd\n", N);
#endif

	// Now we need to generate the user keys of bit-length keysize*2,
	// and the aggregation key sk0.
	mpz_t sk_hold;
	mpz_init(sk_hold);
	for (int i = 0; i < n_sensors; i++)
	{
		mpz_urandomb(sk_hold, state, keysize*2);
		mpz_set(sk[i+1], sk_hold);
		mpz_sub(sk[0], sk[0], sk_hold);
	}
}

void ppa_encrypt(mpz_t c, const mpz_t ski, const mpz_t x, const mpz_t t, const mpz_t N, const mpz_t N2)
{
	// We need to compute the hash of the time t
	mpz_t digest;
	mpz_init(digest);
	ppa_H(digest, t, N);
//	std::cout << "Taking the power in ppa-encrypt\n";
//	gmp_printf("Digest: %Zd\n", digest);
//	gmp_printf("ski: %Zd\n", ski);
//	gmp_printf("N2: %Zd\n", N2);
	// Take modular power of the hash digest with the user key ski
	mpz_powm(c, digest, ski, N2);
	// NB! TRY TO DO THIS BY MULTIPLYING? Take advantage of the binomial theorem
//	std::cout << "Done taking the power in ppa-encrypt\n";
	// Multiply with (1 + xN)
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(tmp, x, N);
	mpz_add_ui(tmp, tmp, 1);
	mpz_mul(c, c, tmp);
	mpz_mod(c, c, N2);

}

void ppa_aggregate_decrypt(mpz_t sum, const mpz_t sk0, const mpz_t t, const mpz_t c[], const mpz_t N, const mpz_t N2, const uint32_t n_sensors)
{
	// We need to compute the hash of the time t
	mpz_t digest;
	mpz_init(digest);
	ppa_H(digest, t, N);

//	gmp_printf("Decrypt digest: %Zd\n", digest);
	// Take modular power of the hash digest with the aggregation key sk0
	mpz_powm(sum, digest, sk0, N2);

	// Compute the product of the ciphertexts
	mpz_t tmp;
	mpz_init(tmp);

	mpz_mul(tmp, c[0], c[1]);
	mpz_mod(tmp, tmp, N2);

	for (int i = 2; i < n_sensors; i++)
	{
		mpz_mul(tmp, tmp, c[i]);
		mpz_mod(tmp, tmp, N2);
	}

//	gmp_printf("Product of ciphertexts: %Zd\n", tmp);
	// Multiply ciphertext product with the hash digest power
	mpz_mul(sum, sum, tmp);
	mpz_mod(sum, sum, N2);

	// Recover the plaintext sum by fast discrete logarithm
	mpz_sub_ui(sum, sum, 1);
//	gmp_printf("Sum: %Zd\n", sum);
	mpz_cdiv_q(sum, sum, N);
}

void ppa_H(mpz_t digest, const mpz_t input, const mpz_t N)
{
	// We first need to compute the SHA-256 digest of the input,
	// to we need to parse the input as a bitstring.
	size_t size = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
	
#ifdef DEBUG
//	std::cout << "Size = " << size << std::endl;
#endif
	
	uint8_t *bitstring = new uint8_t[size] {0};
	mpz_export(bitstring, &size, 1, 1, 0, 0, input);

	// Process the bitstring with SHA-256
	uint8_t sha_digest_array[256/8];
	sha256_process_message(sha_digest_array, bitstring, size/8);

	delete[] bitstring;

	// Import the SHA-256 digest back as as unsigned integer
	// in the range {0, ..., 2^{256} - 1}
	mpz_t sha_digest_integer;
	mpz_init(sha_digest_integer);
	mpz_import(sha_digest_integer, 32, 1, 1, 0, 0, sha_digest_array);

#ifdef DEBUG
//	gmp_printf("SHA-256 digest integer: %Zd\n", sha_digest_integer);
#endif
	
	// Map the sha_digest to an element of the multiplicative group
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(tmp, sha_digest_integer, N);
	
	mpz_set_ui(digest, 1);
	mpz_add(digest, digest, tmp);
}

