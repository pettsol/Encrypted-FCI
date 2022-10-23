#include <gmp.h>

#include "paillier/paillier.h"
#include "encrypted_fci.h"

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
		const uint32_t dim)
{
	// We assume quantization and multiplication with scaling has been done beforehand
	
	// Encrypt the trace
	paillier_encrypt(c_tr, state, trace, g, N, N2);

	// Encrypt the vector
	for (uint32_t i = 0; i < dim; i++)
	{
		// Encrypt elements element-wise
		paillier_encrypt(c_Px[i], state, Px[i], g, N, N2);
	}

	// Encrypt the matrix
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		// Encrypt elements element-wise
		paillier_encrypt(c_P[i], state, P[i], g, N, N2);
	}
}

void efci_fusion(
		mpz_t sum_tr,
		mpz_t sum_Px[],
		mpz_t sum_P[],
		const mpz_t c_tr[],
		const mpz_t c_Px[],
		const mpz_t c_P[],
		const mpz_t N2,
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// Homomorphically add the traces
	mpz_set(sum_tr, c_tr[0]);
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		mpz_mul(sum_tr, sum_tr, c_tr[i]);
		mpz_mod(sum_tr, sum_tr, N2);
	}

	// Homomorphically add the vectors element-wise
	for (uint32_t i = 0; i < dim; i++)
	{
		mpz_set(sum_Px[i], c_Px[i]);
	}
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_mul(sum_Px[j], sum_Px[j], c_Px[i*dim + j]);
			mpz_mod(sum_Px[j], sum_Px[j], N2);
		}
	}

	// Homomorphically add the matrices element-wise
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		mpz_set(sum_P[i], c_P[i]);
	}
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		for (uint32_t j = 0; j < dim*dim; j++)
		{
			mpz_mul(sum_P[j], sum_P[j], c_P[i*dim*dim + j]);
			mpz_mod(sum_P[j], sum_P[j], N2);
		}
	}
}

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
		const uint32_t dim)
{
	mpf_t output;
	mpf_init(output);

	// Decrypt and decode the sum of traces
	mpz_t pt_traces;
	mpz_init(pt_traces);
	paillier_decrypt(pt_traces, sum_tr, lambda, mu, N, p2, q2);
	rho_inv(output, pt_traces, gamma, N);

	double traces = mpf_get_d(output);
	traces = 1/traces;

	// Decrypt and decode the fused vector
	mpz_t *pt_Px = new mpz_t[dim];

	for (uint32_t i = 0; i < dim; i++)
	{
		mpz_init(pt_Px[i]);
		paillier_decrypt(pt_Px[i], sum_Px[i], lambda, mu, N, p2, q2);
		rho_inv(output, pt_Px[i], gamma, N);
		Px[i] = mpf_get_d(output);
		Px[i] = Px[i] * traces;
	}

	// Decrypt and decode the fused matrix
	mpz_t *pt_P = new mpz_t[dim*dim];

	for (uint32_t i = 0; i < dim*dim; i++)
	{
		mpz_init(pt_P[i]);
		paillier_decrypt(pt_P[i], sum_P[i], lambda, mu, N, p2, q2);
		rho_inv(output, pt_P[i], gamma, N);
		P[i] = mpf_get_d(output);
		P[i] = P[i] * traces;
	}

	// Clean up
	delete[] pt_Px;
	delete[] pt_P;
}

void rho(mpz_t out, const mpf_t in, const mpz_t gamma, const mpz_t ptspace)
{
        mpf_t tmp;
        mpf_init(tmp);
        mpf_set_z(tmp, gamma);
        mpf_mul(tmp, in, tmp);

        mpz_set_f(out, tmp);
        mpz_mod(out, out, ptspace);
}

void rho_inv(mpf_t out, const mpz_t in, const mpz_t gamma, const mpz_t ptspace)
{
        mpz_t halfsize;
        mpz_init(halfsize);

        mpz_div_ui(halfsize, ptspace, 2);
        mpz_t test;
        mpz_init(test);

        mpz_sub(test, in, halfsize);

        mpf_t gamma_f;
        mpf_init(gamma_f);

        mpf_set_z(gamma_f, gamma);

        if (mpz_sgn(test) != -1)
        {
                // negative number
                mpz_sub(test, in, ptspace);
        } else {
                mpz_set(test, in);
        }
        mpf_set_z(out, test);
        mpf_div(out, out, gamma_f);
}

