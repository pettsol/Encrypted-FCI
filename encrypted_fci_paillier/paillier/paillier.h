///////////////////////////////////
// This implementation was place //
// in the public domain by       //
// Petter Solnør.                //
//                               //
// Petter Solnør - 23/10/2022    //
///////////////////////////////////

#ifndef PAILLIER_H
#define PAILLIER_H

#include <gmp.h>
#include <stdint.h>

void paillier_keygen(mpz_t N,
		     mpz_t N2,
		     mpz_t p2,
		     mpz_t q2,
		     mpz_t g,
		     mpz_t lambda,
		     mpz_t mu,
		     gmp_randstate_t state,
		     const uint32_t keysize);

void paillier_encrypt(mpz_t c,
		      gmp_randstate_t state,
		      const mpz_t m,
		      const mpz_t g,
		      const mpz_t N,
		      const mpz_t N2);

void paillier_decrypt(mpz_t m,
		      const mpz_t c,
		      const mpz_t lambda,
		      const mpz_t mu,
		      const mpz_t N,
		      const mpz_t p2,
		      const mpz_t q2);

void L(mpz_t output,
       const mpz_t input,
       const mpz_t N);

#endif
