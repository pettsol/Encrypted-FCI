#include <gmp.h>

#include "joye_libert_journal/joye_libert.h"

void efci_encrypt(
		mpz_t c_tr,
		mpz_t c_Px[],
		mpz_t c_P[],
		gmp_randstate_t state,
		const mpz_t trace,
		const mpz_t Px[],
		const mpz_t P[],
		const mpz_t N,
		const mpz_t y,
		const uint32_t dim,
		const uint32_t msgsize);

void efci_fusion(
		mpz_t sum_tr,
		mpz_t sum_Px[],
		mpz_t sum_P[],
		const mpz_t c_tr[],
		const mpz_t c_Px[],
		const mpz_t c_P[],
		const mpz_t N,
		const uint32_t dim,
		const uint32_t n_sensors);

void efci_decrypt(
		double Px[],
		double P[],
		const mpz_t sum_tr,
		const mpz_t sum_Px[],
		const mpz_t sum_P[],
		const mpz_t p,
		const mpz_t y,
		const mpz_t gamma,
		const uint32_t dim,
		const uint32_t msgsize);

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

