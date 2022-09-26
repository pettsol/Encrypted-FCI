///////////////////////////////////////////
//
// Petter Solnør
///////////////////////////////////////////
#include "ppfci.h"
#include "joye_libert_journal/joye_libert.h"
#include <iostream>
#include <unistd.h>

void ppfci_setup(
		gmp_randstate_t rand_state,
		mpz_t y,
		mpz_t p,
		mpz_t N,
		mpz_t ppaN,
		mpz_t upk[],
		mpz_t usk[],
		mpz_t sk[],
		const uint32_t msgsize,
		const uint32_t keysize,
		const uint32_t n_sensors)
{
	// Generate ppa keys
	ppa_setup(ppaN, sk, keysize, n_sensors);

	// Generate mu-labhe master keys
	// mpk = (y, N)
	// msk = (y, p)
	mu_he_setup(N, y, p, msgsize, keysize);

	// Generate mu-labhe user keys
	for (int i = 0; i < n_sensors; i++)
	{
		mu_he_keygen(upk[i], usk[i], rand_state, N, y, msgsize);
	}
}

void ppfci_sensor_encrypt(
		mpz_t c_tr,
		mpz_t label,
		gmp_randstate_t rand_state,
		he_ct *C_tr,
		he_ct C_P[],
		he_ct C_Px[],
		const mpz_t Pm[],
		const mpz_t P[],
		const mpz_t Px[],
		const mpz_t usk,
		const mpz_t y,
		const mpz_t N,
		const mpz_t ski,
		const mpz_t timestep,
		const mpz_t ppaN,
		const mpz_t ppaN2,
		const uint32_t msgsize,
		const uint32_t dim, 
		const uint32_t n_sensors)
{
	// Compute the trace of P
	mpz_t trace;
	mpz_init_set_ui(trace, 0);
	for (uint32_t i = 0; i < dim; i++)
	{
		// Matrix is stored by rows:
		mpz_add(trace, trace, Pm[i*dim + i]);
	}

	// Do PPA encryption of the trace of P
	ppa_encrypt(c_tr, ski, trace, timestep, ppaN, ppaN2);

	// Do MU-LabHE encryption of the trace of P
	// Encrypt?
	mpz_t b;
	mpz_init(b);
	mu_he_encrypt(C_tr, rand_state, b, usk, trace, y, N, label, msgsize);

	mpz_add_ui(label, label, 1);

	// Do element-wise MU-LabHE encryption of P^{-1}
	uint32_t dim2 = dim*dim;
	for (uint32_t i = 0; i < dim2; i++)
	{
		// Encrypt
		mu_he_encrypt(&C_P[i], rand_state, b, usk, P[i], y, N, label, msgsize);
		mpz_add_ui(label, label, 1);
	}

	// Do element-wise MU-LabHE encryption of P^{-1}x
	for (uint32_t i = 0; i < dim; i++)
	{
		// Encrypt
		mu_he_encrypt(&C_Px[i], rand_state, b, usk, Px[i], y, N, label, msgsize);
		mpz_add_ui(label, label, 1);
	}

}

void ppfci_encrypted_fusion(
		mpz_t c_P0[],
		mpz_t c_P0x0[],
		mpz_t m_den,
		gmp_randstate_t rand_state,
	        const mpz_t sk0,
		const mpz_t timestep,
		const mpz_t y,
		const mpz_t N,
		const mpz_t ppaN,
		const mpz_t ppaN2,
		const mpz_t c_tr[],
		const he_ct C_tr[],
		const he_ct C_P[],
		const he_ct C_Px[],
		const uint32_t msgsize,
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// Compute the aggregated sum of traces from the encrypted traces
	// using the privacy-preserving aggregation scheme descrbied in
	// "A Scalable Scheme for Privacy-Preserving Aggregation of Time-Series Data"
	// by Marc Joye and Benoit Libert
	ppa_aggregate_decrypt(m_den, sk0, timestep, c_tr, ppaN, ppaN2, n_sensors);

	mpz_t ptspace;
	mpz_init(ptspace);
	mpz_ui_pow_ui(ptspace, 2, msgsize);
	
	mpz_mul_ui(m_den, m_den, n_sensors-1);
	mpz_mod(m_den, m_den, ptspace);

	// Find the inverse of m_den
	mpz_t gamma_tr;
	mpz_init(gamma_tr);
	mpz_ui_pow_ui(gamma_tr, 2, 10);
	mpf_t mpf_true_den;
	mpf_init(mpf_true_den);
	rho_inv(mpf_true_den, m_den, gamma_tr, ptspace);
	double true_den;
	true_den = mpf_get_d(mpf_true_den);
	double inv_den = 1/true_den;
	mpf_set_d(mpf_true_den, inv_den);

	mpz_t gamma;
	mpz_init(gamma);
	mpz_ui_pow_ui(gamma, 2, 20);
	rho(m_den, mpf_true_den, gamma, ptspace);
	
	// Compute the encrypted cumulative sum of traces
	he_ct sum;
	mu_he_eval_add(&sum, &C_tr[0], &C_tr[1], N, msgsize);
	for (uint32_t i = 2; i < n_sensors; i++)
	{
		mu_he_eval_add(&sum, &sum, &C_tr[i], N, msgsize);
	}

	// Compute the weights and multiply
	he_ct numerator, weight;
	uint32_t dim2 = dim*dim;
	mpz_t *c_P0_array = new mpz_t[n_sensors*dim2];
	mpz_t *c_P0x0_array = new mpz_t[n_sensors*dim];
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		// Compute encrypted weight for sensor system i
		mu_he_eval_sub(&numerator, &sum, &C_tr[i], N, msgsize);
		mu_he_eval_cmult(&weight, &numerator, m_den, N, msgsize);

		// Perform element-wise ciphertext multiplication 
		// between weight i and matrix P[i]
		for (uint32_t j = 0; j < dim2; j++)
		{
			mpz_init(c_P0_array[i*dim2 + j]);
			mu_he_eval_mult(c_P0_array[i*dim2 + j], rand_state, &weight, &C_P[i*dim2 + j], y, N, msgsize);
		}

		// Perform element-wise ciphertext multiplication
		// between weight i and vector P[i]x[i]
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init(c_P0x0_array[i*dim + j]);
			mu_he_eval_mult(c_P0x0_array[i*dim + j], rand_state, &weight, &C_Px[i*dim + j], y, N, msgsize);
		}
	}

	// Sum the weighted contributions of P0
	for (uint32_t i = 0; i < dim2; i++)
	{
		mu_he_eval_add(c_P0[i], c_P0_array[i], c_P0_array[dim2 + i], N);
		for (uint32_t j = 2; j < n_sensors; j++)
		{
			mu_he_eval_add(c_P0[i], c_P0[i], c_P0_array[j*dim2 + i], N);
		}
	}

	// Sum the weighted contributions of P0x0
	for (uint32_t i = 0; i < dim; i++)
	{
		mu_he_eval_add(c_P0x0[i], c_P0x0_array[i], c_P0x0_array[dim + i], N);
		for (uint32_t j = 2; j < n_sensors; j++)
		{
			mu_he_eval_add(c_P0x0[i], c_P0x0[i], c_P0x0_array[j*dim + i], N);
		}
	}

	delete[] c_P0_array;
	delete[] c_P0x0_array;
}

void ppfci_decrypt(
		mpz_t P0[],
		mpz_t P0x0[],
		mpz_t label,
		const mpz_t c_P0[],
		const mpz_t c_P0x0[],
		const mpz_t m_den,
		const mpz_t usk[],
		const mpz_t msk,
		const mpz_t y,
		const uint32_t msgsize,
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// NB! WE ONLY NEED TO RECOVER THE SECRET KEYS ONCE,
	//     AND, THEREFORE, ASSUME THAT THIS IS DONE BEFORE
	//     WE DECRYPT.
	// We have to evaluate the labeled program P(f, tau_i)
	// First we need to recover the secret keys.
	//mpz_t *usk = new mpz_t[n_sensors];
	//for (uint32_t i = 0; i < n_sensors; i++)
	//{
	//	mpz_init(usk[i]);
	//	joye_libert_decrypt(usk[i], upk[i], msk, y, msgsize);
	//}

	// We then need to compute the b's
	mpz_t label_counter, input;
	mpz_init_set(label_counter, label);
	mpz_init(input);

	uint32_t dim2 = dim*dim;
	// Declare array to hold b's corresponding to P, Px, and Trace
	mpz_t *b_P = new mpz_t[n_sensors*dim2];
	mpz_t *b_Px = new mpz_t[n_sensors*dim];
	mpz_t *b_trace = new mpz_t[n_sensors];

	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init(b_trace[i]);
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init(b_Px[j + i*dim]);
		}
		for (uint32_t j = 0; j < dim2; j++)
		{
			mpz_init(b_P[j + i*dim2]);
		}
	}

	uint8_t digest[32] = {0};

	// Compute b's for each sensor
	mpz_t ptspace;
	mpz_init(ptspace);
	mpz_ui_pow_ui(ptspace, 2, msgsize);
	
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_add(input, usk[i], label);
		mpz_mod(input, input, ptspace);

	        size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       		uint8_t *input_array = new uint8_t[size_array] {0};
		size_t size;
	        mpz_export(input_array, &size, 1, 1, 0, 0, input);
		sha256_process_message(digest, input_array, size_array);

		// We first compute the b corresponding the the trace // NEW
		mpz_import(b_trace[i], 32, 1, 1, 0, 0, digest);
		mpz_mod(b_trace[i], b_trace[i], ptspace);
	
		mpz_add_ui(label, label, 1);

		delete[] input_array;

		// For the P matrix encryptions
		for (uint32_t j = 0; j < dim2; j++)
		{
			// Add the label and the user secret key to produce the
			// input to the hash function. Then export it to a byte string
			// and compute the digest before importing the result again.
			mpz_add(input, usk[i], label);
			mpz_mod(input, input, ptspace);
	        	size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t *input_array = new uint8_t[size_array] {0};
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			size_t size;
			sha256_process_message(digest, input_array, size_array);
			
			// Compute b_P
			mpz_import(b_P[j + i*dim2], 32, 1, 1, 0, 0, digest);
			mpz_mod(b_P[j + i*dim2], b_P[j + i*dim2], ptspace);
			mpz_add_ui(label, label, 1);

			delete[] input_array;
		}

		// For the Px vector
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_add(input, usk[i], label);
			mpz_mod(input, input, ptspace);
	        	size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t *input_array = new uint8_t[size_array] {0};
			size_t size;
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			sha256_process_message(digest, input_array, size_array);
			
			// Compute b_Px
			mpz_import(b_Px[j + i*dim], 32, 1, 1, 0, 0, digest);
			mpz_mod(b_Px[j + i*dim], b_Px[j + i*dim], ptspace);
			mpz_add_ui(label, label, 1);

			delete[] input_array;
		}
	}

	// Now the b's have been computed. The next step is to evaluate the function f
	// over the b's.
	
	// First we compute the b associated with the sum of traces
	
	mpz_t b_sum;
	mpz_init_set(b_sum, b_trace[0]);

	for (uint32_t i = 1; i < n_sensors; i++)
	{
		mpz_add(b_sum, b_sum, b_trace[i]);
		mpz_mod(b_sum, b_sum, ptspace);
	}

	mpz_t *b_wP0 = new mpz_t[dim2];
	mpz_t *b_wP0x0 = new mpz_t[dim];

	for (int i = 0; i < dim2; i++)
	{
		mpz_init_set_ui(b_wP0[i], 0);
	}
	for (int i = 0; i < dim; i++)
	{
		mpz_init_set_ui(b_wP0x0[i], 0);
	}

	// Compute the b's corresponding to weights
	mpz_t b_weight;
	mpz_init(b_weight);

	mpz_t hold;
	mpz_init(hold);

	for (uint32_t i = 0; i < n_sensors; i++)
	{
		//mpz_sub(b_weight, b_sum, b[i*(1 + dim + dim2)]); // OLD
		mpz_sub(b_weight, b_sum, b_trace[i]); // NEW
		mpz_mod(b_weight, b_weight, ptspace);
		mpz_mul(b_weight, b_weight, m_den);
		mpz_mod(b_weight, b_weight, ptspace);

		// Compute the weighted b_P matrix of sensor i
		// and add to the fused b_P0 matrix?
		for (uint32_t j = 0; j < dim2; j++)
		{
			mpz_mul(hold, b_weight, b_P[j + i*dim2]);
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0[j], b_wP0[j], hold);
			mpz_mod(b_wP0[j], b_wP0[j], ptspace);
		}

		// Compute the weighted b_Px vector of sensor i
		// and add to the fused b_P0x0 vector?
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_mul(hold, b_weight, b_Px[j + i*dim]);
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0x0[j], b_wP0x0[j], hold);
			mpz_mod(b_wP0x0[j], b_wP0x0[j], ptspace);
		}
	}

	// Decrypt the P0 matrix
	for (uint32_t i = 0; i < dim2; i++)
	{
		mu_he_decrypt(P0[i], c_P0[i], b_wP0[i], msk, y, msgsize);
	}

	// Decrypt the P0x0 vector
	for (uint32_t i = 0; i < dim; i++)
	{
		mu_he_decrypt(P0x0[i], c_P0x0[i], b_wP0x0[i], msk, y, msgsize);
	}

	//delete[] usk;
	delete[] b_P;
	delete[] b_Px;
	delete[] b_trace;
	delete[] b_wP0;
	delete[] b_wP0x0;
}

void ppfci_normalize(
                double P0_matrix[],
                double P0x0_vector[],
                const mpz_t P0[],
                const mpz_t P0x0[],
                const mpz_t m_den,
                const mpz_t ptspace,
                const mpz_t gamma,
		const uint32_t dim)
{
	// We perform 2 multiplications per element,
	// so the cumulative scaling is given by
	// c_gamma = gamma*gamma*gamma
	mpz_t c_gamma;
	mpz_init(c_gamma);
	//mpz_mul(c_gamma, gamma, gamma);
	//mpz_mul(c_gamma, c_gamma, gamma);
	// 64 bit test
	mpz_ui_pow_ui(c_gamma, 2, 50);

	mpf_t tmp_float;
	mpf_init(tmp_float);
	// Inverse map matrix
	for (int i = 0; i < dim*dim; i++)
	{
		rho_inv(tmp_float, P0[i], c_gamma, ptspace);
		P0_matrix[i] = mpf_get_d(tmp_float);
	}

	// Inverse map vector
	for (int i = 0; i < dim; i++)
	{
		rho_inv(tmp_float, P0x0[i], c_gamma, ptspace);
		P0x0_vector[i] = mpf_get_d(tmp_float);
	}
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
