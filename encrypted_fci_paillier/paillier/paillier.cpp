///////////////////////////////////
// This implementation was place //
// in the public domain by       //
// Petter Solnør.                //
//                               //
// Petter Solnør - 23/10/2022    //
///////////////////////////////////

#include "paillier.h"

void paillier_keygen(mpz_t N, 
		     mpz_t N2,
		     mpz_t p2,
		     mpz_t q2,
		     mpz_t g,
		     mpz_t lambda,
		     mpz_t mu,
		     gmp_randstate_t state,
		     const uint32_t keysize)
{
	mpz_t p, q, div, gcd;
	mpz_init(p);
	mpz_init(q);
	mpz_init(div);
	mpz_init(gcd);

	uint8_t flag = 1;

	// Find p and q such that p does not divide q-1
	// and q does not divide p-1
	while (flag)
	{
		// Find random prime p of size keysize/2
		mpz_urandomb(p, state, keysize/2);
		mpz_nextprime(p, p);

		// Find random prime q of size keysize/2
		mpz_urandomb(q, state, keysize/2);
		mpz_nextprime(q, q);

		// Check that property is satisfied
		mpz_sub_ui(div, q, 1);
		mpz_gcd(gcd, p, div);

		if (!(mpz_cmp_ui(gcd, 0)))
		{
			// Not coprime, try again
			continue;
		}
		
		mpz_sub_ui(div, p, 1);
		mpz_gcd(gcd, q, div);

		if (!(mpz_cmp_ui(gcd, 0)))
		{
			// Not coprime, try again
			continue;
		}

		// Success!
		flag = 0;
	}

	// Compute p2 and q2 for CRT
	mpz_mul(p2, p, p);
	mpz_mul(q2, q, q);

	// Compute modulus N
	mpz_mul(N, p, q);

	// Compute N^{2}
	mpz_mul(N2, N, N);

	// Set g = N + 1 for efficiency
	mpz_add_ui(g, N, 1);

	// Set lambda = phi(N) = (p-1)(q-1)
	mpz_t tmp1, tmp2;
	mpz_init(tmp1);
	mpz_init(tmp2);

	mpz_sub_ui(tmp1, p, 1);
	mpz_sub_ui(tmp2, q, 1);
	mpz_mul(lambda, tmp1, tmp2);

	// Set mu = multiplicative inverse of lambda modulo N
	mpz_invert(mu, lambda, N);	
}

void paillier_encrypt(mpz_t c,
		      gmp_randstate_t state,
		      const mpz_t m,
		      const mpz_t g,
		      const mpz_t N,
		      const mpz_t N2)
{
	// Pick a random r where 0 < r < N holds
	mpz_t r;
	mpz_init(r);
	mpz_urandomm(r, state, N);

	// Encrypt
	mpz_t tmp;
	mpz_init(tmp);

	mpz_powm(tmp, r, N, N2);

	// Take advantage of g^{m} = (1 + N)^{m] = 1 + mN
	mpz_mul(c, m, N);
	mpz_add_ui(c, c, 1);

	// Complete encryption
	mpz_mul(c, c, tmp);
	mpz_mod(c, c, N2);
}

void paillier_decrypt(mpz_t m,
		      const mpz_t c,
		      const mpz_t lambda,
		      const mpz_t mu,
		      const mpz_t N,
		      const mpz_t p2,
		      const mpz_t q2)
{
	// We use the Chinese Remainder Theorem to accelerate
	// the decryption procedure. By the CRT, we can replace
	// one exponentiation with a large exponent and a large
	// modulus, with two exponentiations with smaller exponents
	// and smaller moduli.
	mpz_t tmp_p, tmp_q, x;
	mpz_init(tmp_p);
	mpz_init(tmp_q);
	mpz_init(x);

	mpz_powm(tmp_p, c, lambda, p2);
	mpz_powm(tmp_q, c, lambda, q2);
	mpz_sub(x, tmp_p, tmp_q);
	mpz_invert(tmp_p, q2, p2);
	mpz_mul(x, x, tmp_p);
	mpz_mod(x, x, p2);
	mpz_mul(x, x, q2);
	mpz_add(x, x, tmp_q);

	// Finish decryption
	L(m, x, N);
	mpz_mul(m, m, mu);
	mpz_mod(m, m, N);
}

void L(mpz_t output, const mpz_t input, const mpz_t N)
{
	mpz_sub_ui(output, input, 1);
	mpz_divexact(output, output, N);
}
