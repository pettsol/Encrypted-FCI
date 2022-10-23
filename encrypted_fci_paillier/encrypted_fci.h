#include <gmp.h>

#include "paillier/paillier.h"

void efci_encrypt(
		mpz_t c_tr,
		mpz_t c_Px[],
		mpz_t c_P[],
		gmp_randstate_t state,
		const mpz_t trace,
		const mpz_t Px[],
		const mpz_t P[],
		const mpz_t N,
		const mpz_t N2,
		const mpz_t g,
		const uint32_t dim);

void efci_fusion(
		mpz_t sum_tr,
		mpz_t sum_Px[],
		mpz_t sum_P[],
		const mpz_t c_tr[],
		const mpz_t c_Px[],
		const mpz_t c_P[],
		const mpz_t N2,
		const uint32_t dim,
		const uint32_t n_sensors);

void efci_decrypt(
		double Px[],
		double P[],
		const mpz_t sum_tr,
		const mpz_t sum_Px[],
		const mpz_t sum_P[],
		const mpz_t lambda,
		const mpz_t mu,
		const mpz_t N,
		const mpz_t p2,
		const mpz_t q2,
		const mpz_t gamma,
		const uint32_t dim);

void rho(
                mpz_t out,
                const mpf_t in,
                const mpz_t gamma,
                const mpz_t ptspace);

void rho_inv(
                mpf_t out,
                const mpz_t in,
                const mpz_t gamma,
                const mpz_t ptspace);

