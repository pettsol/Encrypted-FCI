///////////////////////////////////////////
//
// Petter Solnør
///////////////////////////////////////////
#include "ppfci.h"
#include "joye_libert_journal/joye_libert.h"

void sensor_encrypt(
		mpz_t c_tr,
		mpz_t label,
		gmp_randstate_t rand_state,
		he_ct *C_tr,
		he_ct C_P[],
		he_ct C_Px[],
		const mpz_t Pm[];
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
	//! NB! Trace is computed from original matrix, the rest is done on the inverse matrix


	// Compute the trace of P
	mpz_t trace;
	mpz_init_set_ui(trace, 0);
	for (uint32_t i = 0; i < dim; i++)
	{
		// What format is the matrix stored in? If by rows:
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

void encrypted_fci(
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
	mpz_mul_ui(m_den, m_den, n_sensors-1);

	// If m_den is not invertible plaintext space, then 
	// add 1.
	mpz_t tmp;
	mpz_init(tmp);
	if (mpz_gcd_ui(NULL, m_den, 2) != 1)
	{
		mpz_add_ui(tmp, m_den, 1);
	}
	else
	{
		mpz_set(tmp, m_den);
	}
	
	// Compute the multiplicative inverse of tmp in
	// plaintext space.
	mpz_t ptspace;
	mpz_init(ptspace);
	mpz_ui_pow_ui(ptspace, 2, msgsize);
	mpz_invert(tmp, tmp, ptspace);

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
	mpz_t c_P0_array[n_sensors*dim2], c_P0x0_array[n_sensors*dim2];
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		// Compute encrypted weight for sensor system i
		mu_he_eval_sub(&numerator, &sum, &C_tr[i], N, msgsize);
		mu_he_eval_cmult(&weight, &numerator, tmp, N, msgsize);

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
			mpz_init(c_P0x0_array[i*dim2 + j]);
			mu_he_eval_mult(c_P0x0_array[i*dim2 + j], rand_state, &weight, &C_Px[i*dim2 + j], y, N, msgsize);
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
			mu_he_eval_add(c_P0x0[i], c_P0x0[i], c_P0x0_array[j*dim], N);
		}
	}
}

void decrypt(
		mpz_t P0[],
		mpz_t P0x0[],
		mpz_t label,
		const mpz_t c_P0[],
		const mpz_t c_P0x0[],
		const mpz_t m_den,
		const mpz_t upk[],
		const mpz_t msk,
		const mpz_t y,
		const uint32_t msgsize,
		const uint32_t dim,
		const uint32_t n_sensors)
{
	// We have to evaluate the labeled program P(f, tau_i)
	// First we need to recover the secret keys.
	mpz_t usk[n_sensors];
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init(usk[i]);
		joye_libert_decrypt(usk[i], upk[i], msk, y, msgsize);
	}

	// We then need to compute the b's
	mpz_t label_counter, input;
	mpz_init_set(label_counter, label);
	mpz_init(input);

	uint32_t dim2 = dim*dim;
	mpz_t b[n_sensors*(1 + dim + dim2)];

	uint8_t digest[32];

	// Compute b's for each sensor
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_set(label, label_counter);
		// How many labels are there for each sensor?
		mpz_init(b[i*(1 + dim + dim2)]);
		mpz_add(input, usk[i], label);
	        size_t size = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       		uint8_t input_array[size] = {0};
	        mpz_export(input_array, &size, 1, 1, 0, 0, input);
		sha256_process_message(digest, input_array, size);
		mpz_import(b[i*(1 + dim + dim2)], msgsize/32, 1, 4, 0, 0, digest);

		mpz_add_ui(label, label, 1);

		// For the P matrix encryptions
		for (uint32_t j = 0; j < dim2; j++)
		{
			//
			mpz_init(b[i*(1 + dim + dim2) + 1 + j]);

			// Add the label and the user secret key to produce the
			// input to the hash function. Then export it to a byte string
			// and compute the digest before importing the result again.
			mpz_add(input, usk[i], label);
	        	size_t size = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t input_array[size] = {0};
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			sha256_process_message(digest, input_array, size);
			mpz_import(b[i*(1 + dim + dim2)], msgsize/32, 1, 4, 0, 0, digest);

			mpz_add_ui(label, label, 1);
		}

		// For the Px vector
		for (uint32_t j = 0; j < dim; j++)
		{
			//
			mpz_init(b[i*(1 + dim + dim2) + 1 + dim2 + j]);

			mpz_add(input, usk[i], label);
	        	size_t size = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t input_array[size] = {0};
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			sha256_process_message(digest, input_array, size);
			mpz_import(b[i*(1 + dim + dim2)], msgsize/32, 1, 4, 0, 0, digest);

			mpz_add_ui(label, label, 1);
		}
	}

	// Now the b's have been computed. The next step is to evaluate the function f
	// over the b's.
	
	// First we compute the b associated with the sum of traces
	mpz_t b_sum, ptspace;
	mpz_init_set(b_sum, b[0]);
	mpz_init(ptspace);

	mpz_ui_pow_ui(ptspace, 2, msgsize);
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		mpz_add(b_sum, b_sum, b[i*(1 + dim + dim2)]);
		mpz_mod(b_sum, b_sum, ptspace);
	}

	// Check that m_den is invertible. If it is not, then add 1.
	mpz_t tmp;
	mpz_init(tmp);
	if (mpz_gcd_ui(NULL, m_den, 2) != 1)
	{
		mpz_add_ui(tmp, m_den, 1);
	}
	else
	{
		mpz_set(tmp, m_den);
	}
	mpz_invert(tmp, tmp, ptspace);

	mpz_t b_wP0[dim2];
	mpz_t b_wP0x0[dim];

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

	mpz_t b_P[n_sensors*dim2];
	mpz_t b_Px[n_sensors*dim];
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		mpz_sub(b_weight, b_sum, b[i*(1 + dim + dim2)]);
		mpz_mod(b_weight, b_weight, ptspace);
		mpz_mul(b_weight, b_weight, tmp);
		mpz_mod(b_weight, b_weight, ptspace);

		// Compute the weighted b_P matrix of sensor i
		// and add to the fused b_P0 matrix?
		for (uint32_t j = 0; j < dim2; j++)
		{
			//
			mpz_mul(hold, b_weight, b_P[1 + j + i*(1 + dim + dim2)]);
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0[j], b_wP0[j], hold);
			mpz_mod(b_wP0[j], b_wP0[j], ptspace);
		}

		// Compute the weighted b_Px vector of sensor i
		// and add to the fused b_P0x0 vector?
		for (uint32_t j = 0; j < dim; j++)
		{
			//
			mpz_mul(hold, b_weight, b_P[1 + dim2 + j + i*(1 + dim + dim2)]);
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0x0[j], b_wP0x0[j], hold);
			mpz_mod(b_wP0x0[j], b_wP0x0[j], ptspace);
		}
	}

	// Decrypt the P0 matrix
	for (uint32_t i = 1; i < dim2; i++)
	{
		mu_he_decrypt(P0[i], c_P0[i], b_wP0[i], msk, y, msgsize);
	}

	// Decrypt the P0x0 vector
	for (uint32_t i = 0; i < dim; i++)
	{
		mu_he_decrypt(P0x0[i], c_P0x0[i], b_wP0x0[i], msk, y, msgsize);
	}
}


