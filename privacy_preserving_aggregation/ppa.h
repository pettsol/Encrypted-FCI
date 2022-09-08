/////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////

#ifndef PPA_H
#define PPA_H

#include "SHA-256/sha-256.h"
#include <gmp.h>
#include <stdint.h>

void ppa_setup(mpz_t N, mpz_t sk[], const uint32_t keysize, const uint32_t n_sensors);

void ppa_encrypt(mpz_t c, const mpz_t ski, const mpz_t x, const mpz_t t, const mpz_t N, const mpz_t N2);

void ppa_aggregate_decrypt(mpz_t sum, const mpz_t sk0, const mpz_t t, const mpz_t c[], const mpz_t N, const mpz_t N2, const uint32_t n_sensors);

void ppa_H(mpz_t digest, const mpz_t input, const mpz_t N);

#endif
