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
		// DEBUG - CHECK SIZE OF USK[i] and USK[i]
		//gmp_printf("usk[%u] = %Zd\n", i, usk[i]);
		//size_t test_size = (mpz_sizeinbase(usk[i], 2) + CHAR_BIT-1)/CHAR_BIT;
		//std::cout << "Size of usk[i]: " << test_size << std::endl;

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
		// What format is the matrix stored in? If by rows:
		mpz_add(trace, trace, Pm[i*dim + i]);
	}

	//gmp_printf("Trace: %Zd\n", trace);

	std::cout << "Before ppa_encrypt\n";
	// Do PPA encryption of the trace of P
	ppa_encrypt(c_tr, ski, trace, timestep, ppaN, ppaN2);

	// Do MU-LabHE encryption of the trace of P
	// Encrypt?
	mpz_t b;
	mpz_init(b);
	mu_he_encrypt(C_tr, rand_state, b, usk, trace, y, N, label, msgsize);

	//gmp_printf("Label encrypt trace: %Zd\n", label);
	//gmp_printf("Encrypt trace b: %Zd\n", b);
	
	mpz_add_ui(label, label, 1);

	std::cout << "Element-wise encryption\n";
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
	
	//gmp_printf("Recovered sum of trace: %Zd\n", m_den);
	
	mpz_mul_ui(m_den, m_den, n_sensors-1);

	//gmp_printf("(N-1)*Trace: %Zd\n", m_den);

	// DEBUG - Check m_den
	//gmp_printf("Encrypt m_den = %Zd\n", m_den);
	
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

	// DEBUG - Check m_den_inv
	//gmp_printf("Encrypt m_den_inv = %Zd\n", tmp);

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
		//std::cout << "CT mult between weight and matrix for sensor " << i << std::endl;
		for (uint32_t j = 0; j < dim2; j++)
		{
			mpz_init(c_P0_array[i*dim2 + j]);
			mu_he_eval_mult(c_P0_array[i*dim2 + j], rand_state, &weight, &C_P[i*dim2 + j], y, N, msgsize);
		}

		//std::cout << "CT mult between weight and matrix for sensor " << i << std::endl;
		// Perform element-wise ciphertext multiplication
		// between weight i and vector P[i]x[i]
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init(c_P0x0_array[i*dim + j]);
			mu_he_eval_mult(c_P0x0_array[i*dim + j], rand_state, &weight, &C_Px[i*dim + j], y, N, msgsize);
		}
	}

	// Sum the weighted contributions of P0
	std::cout << "Summing weighted contributions of P0\n";
	for (uint32_t i = 0; i < dim2; i++)
	{
		mu_he_eval_add(c_P0[i], c_P0_array[i], c_P0_array[dim2 + i], N);
		for (uint32_t j = 2; j < n_sensors; j++)
		{
			mu_he_eval_add(c_P0[i], c_P0[i], c_P0_array[j*dim2 + i], N);
		}
	}

	std::cout << "Summing weighted contributions of P0x0\n";
	// Sum the weighted contributions of P0x0
	for (uint32_t i = 0; i < dim; i++)
	{
		mu_he_eval_add(c_P0x0[i], c_P0x0_array[i], c_P0x0_array[dim + i], N);
		for (uint32_t j = 2; j < n_sensors; j++)
		{
			mu_he_eval_add(c_P0x0[i], c_P0x0[i], c_P0x0_array[j*dim + i], N);
		}
	}
	std::cout << "Finished iteration\n";
}

void ppfci_decrypt(
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

	std::cout << "Recovering usk" << std::endl;
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init(usk[i]);
		joye_libert_decrypt(usk[i], upk[i], msk, y, msgsize);
		// DEBUG - CHECK SIZE OF USK[i] and USK[i]
		//gmp_printf("usk[%u] = %Zd\n", i, usk[i]);
		//size_t test_size = (mpz_sizeinbase(usk[i], 2) + CHAR_BIT-1)/CHAR_BIT;
		//std::cout << "Size of usk[i]: " << test_size << std::endl;
	}

	// We then need to compute the b's
	mpz_t label_counter, input;
	mpz_init_set(label_counter, label);
	mpz_init(input);

	uint32_t dim2 = dim*dim;
	// Declare array to hold b's corresponding to P, Px, and Trace
	mpz_t b_P[n_sensors*dim2];
	mpz_t b_Px[n_sensors*dim];
	mpz_t b_trace[n_sensors];

	/*mpz_t b[n_sensors*(1 + dim + dim2)];
	//for (uint32_t i = 0; i < n_sensors*(1+dim+dim2); i++)
	//{
	//	mpz_init(b[i]);
	//} OLD */
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init(b_trace[i]);
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init(b_Px[j + i*dim]);
			//for (uint32_t k = 0; k < dim; k++)
			//{
			//	mpz_init(b_P[k + j*dim + i*dim2]);
			//}
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
	
	std::cout << "Computing the b's\n";
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		//mpz_set(label, label_counter);
		// How many labels are there for each sensor?
		//mpz_init(b[i*(1 + dim + dim2)]);
		mpz_add(input, usk[i], label);

		// DEBUG - CHECK SIZE OF USK[i]
		size_t test_size = (mpz_sizeinbase(usk[i], 2) + CHAR_BIT-1)/CHAR_BIT;
		std::cout << "Size of usk[i]: " << test_size << std::endl;
 

		mpz_mod(input, input, ptspace);
	        size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       		uint8_t input_array[size_array] = {0};
		size_t size;
	        mpz_export(input_array, &size, 1, 1, 0, 0, input);
		sha256_process_message(digest, input_array, size_array);
		//mpz_import(b[i*(1 + dim + dim2)], msgsize/32, 1, 4, 0, 0, digest); //OLD

		// We first compute the b corresponding the the trace // NEW
		mpz_import(b_trace[i], 32, 1, 1, 0, 0, digest);	
	
		std::cout << "Size_array: " << size_array << std::endl;
		gmp_printf("Label decrypt trace: %Zd\n", label);	
		gmp_printf("Decrypt trace b: %Zd\n", b_trace[i]);

		mpz_add_ui(label, label, 1);

		// For the P matrix encryptions
		std::cout << "Recovering P matrix b's for sensor " << i << std::endl;
		for (uint32_t j = 0; j < dim2; j++)
		{
			//
			//mpz_init(b[i*(1 + dim + dim2) + 1 + j]);

			// Add the label and the user secret key to produce the
			// input to the hash function. Then export it to a byte string
			// and compute the digest before importing the result again.
			mpz_add(input, usk[i], label);
	        	size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t input_array[size_array] = {0};
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			size_t size;
			sha256_process_message(digest, input_array, size_array);
			//mpz_import(b[i*(1 + dim + dim2) + 1 + j], msgsize/32, 1, 4, 0, 0, digest); // OLD
			
			// NEW compute b_P
			//mpz_import(b_P[j + i*dim2], msgsize/32, 1, 4, 0, 0, digest); 
			mpz_import(b_P[j + i*dim2], 32, 1, 1, 0, 0, digest);
			mpz_add_ui(label, label, 1);
		}

		// For the Px vector
		std::cout << "Recovering Px vector b's for sensor " << i << std::endl;
		for (uint32_t j = 0; j < dim; j++)
		{
			//
			//mpz_init(b[i*(1 + dim + dim2) + 1 + dim2 + j]);

			mpz_add(input, usk[i], label);
	        	size_t size_array = (mpz_sizeinbase(input, 2) + CHAR_BIT-1)/CHAR_BIT;
       			uint8_t input_array[size_array] = {0};
			size_t size;
	        	mpz_export(input_array, &size, 1, 1, 0, 0, input);
			sha256_process_message(digest, input_array, size_array);
			//mpz_import(b[i*(1 + dim + dim2) + 1 + dim2 + j], msgsize/32, 1, 4, 0, 0, digest); // OLD
			
			// NEW compute b_Px
			//mpz_import(b_Px[j + i*dim], msgsize/32, 1, 4, 0, 0, digest);
			mpz_import(b_Px[j + i*dim], 32, 1, 1, 0, 0, digest);

			mpz_add_ui(label, label, 1);
		}
	}

	// Now the b's have been computed. The next step is to evaluate the function f
	// over the b's.
	
	// First we compute the b associated with the sum of traces
	
	// DEBUG - Check m_den
	gmp_printf("Decrypt m_den = %Zd\n", m_den);
	
	mpz_t b_sum;
	mpz_init_set(b_sum, b_trace[0]);

	std::cout << "Computing sum of traces b's\n";
	for (uint32_t i = 1; i < n_sensors; i++)
	{
		//mpz_add(b_sum, b_sum, b[i*(1 + dim + dim2)]);
		mpz_add(b_sum, b_sum, b_trace[i]);
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

	// DEBUG - Check m_den_inv
	gmp_printf("Decrypt m_den_inv = %Zd\n", tmp);


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

	/*
	mpz_t b_P[n_sensors*dim2];
	mpz_t b_Px[n_sensors*dim];

	// Initialize b_P and b_Px
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		for (uint32_t j = 0; j < dim2; j++)
		{
			mpz_init(b_P[j + i*dim2]);
		}
		for (uint32_t j = 0; j < dim; j++)
		{
			mpz_init(b_Px[j + i*dim]);
		}
	}*/

	std::cout << "Computing b's corresponding to weights\n";
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		//mpz_sub(b_weight, b_sum, b[i*(1 + dim + dim2)]); // OLD
		mpz_sub(b_weight, b_sum, b_trace[i]); // NEW
		mpz_mod(b_weight, b_weight, ptspace);
		mpz_mul(b_weight, b_weight, tmp);
		mpz_mod(b_weight, b_weight, ptspace);

		// Compute the weighted b_P matrix of sensor i
		// and add to the fused b_P0 matrix?
		for (uint32_t j = 0; j < dim2; j++)
		{
			//
			//mpz_mul(hold, b_weight, b_P[1 + j + i*(1 + dim + dim2)]); // This should just be the b? Not bp
			//mpz_mul(hold, b_weight, b[1 + j + i*(1 + dim + dim2)]);  // OLD
			mpz_mul(hold, b_weight, b_P[j + i*dim2]); // NEW
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0[j], b_wP0[j], hold);
			mpz_mod(b_wP0[j], b_wP0[j], ptspace);
		}

		// Compute the weighted b_Px vector of sensor i
		// and add to the fused b_P0x0 vector?
		for (uint32_t j = 0; j < dim; j++)
		{
			//
			//mpz_mul(hold, b_weight, b_P[1 + dim2 + j + i*(1 + dim + dim2)]); // Same as above
			//mpz_mul(hold, b_weight, b[1 + dim2 + j + i*(1 + dim + dim2)]); // OLD
			mpz_mul(hold, b_weight, b_Px[j + i*dim]);
			mpz_mod(hold, hold, ptspace);
			mpz_add(b_wP0x0[j], b_wP0x0[j], hold);
			mpz_mod(b_wP0x0[j], b_wP0x0[j], ptspace);
		}
	}

	//sleep(5);

	// Decrypt the P0 matrix
	std::cout << "Decrypting P0 matrix\n";
	for (uint32_t i = 0; i < dim2; i++)
	{
		mu_he_decrypt(P0[i], c_P0[i], b_wP0[i], msk, y, msgsize);
	}

	// Decrypt the P0x0 vector
	std::cout << "Decrypting P0x0 vector\n";
	for (uint32_t i = 0; i < dim; i++)
	{
		//std::cout << "Element: " << i << std::endl;
		//sleep(2);
		mu_he_decrypt(P0x0[i], c_P0x0[i], b_wP0x0[i], msk, y, msgsize);
	}

	std::cout << "Done decrypting\n";
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
	mpz_mul(c_gamma, gamma, gamma);
	mpz_mul(c_gamma, c_gamma, gamma);

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
	
	// We need to check if we should normalize
	if (mpz_gcd_ui(NULL, m_den, 2) != 1)
	{
		mpz_t tmp;
		mpz_init(tmp);
		mpz_add_ui(tmp, m_den, 1);
		
		// We need to map m_den and tmp to doubles,
		// compute the ratio tmp/m_den, and multiply each
		// element in P0_matrix and P0x0_vector with the
		// ratio to normalize.
		double m_den_double, tmp_double;
		rho_inv(tmp_float, m_den, gamma, ptspace);
		m_den_double = mpf_get_d(tmp_float);
		rho_inv(tmp_float, tmp, gamma, ptspace);
		tmp_double = mpf_get_d(tmp_float);
		double ratio = tmp_double/m_den_double;

		// Normalize P0
		for (int i = 0; i < dim*dim; i++)
		{
			P0_matrix[i] = P0_matrix[i] * ratio;
		}

		// Normalize P0x0
		for (int i = 0; i < dim; i++)
		{
			P0x0_vector[i] = P0x0_vector[i] * ratio;
		}
	}
}

void rho(mpz_t out, const mpf_t in, const mpz_t gamma, const mpz_t ptspace)
{
        mpf_t tmp;
        mpf_init(tmp);
        mpf_set_z(tmp, gamma);
        mpf_mul(tmp, in, tmp);

//        mpz_t size;
//        mpz_init(size);
//        mpz_ui_pow_ui(size, 2, msgsize);

        mpz_set_f(out, tmp);
        mpz_mod(out, out, ptspace);
}

void rho_inv(mpf_t out, const mpz_t in, const mpz_t gamma, const mpz_t ptspace)
{
        
//        mpz_t size, halfsize;
//        mpz_init(size);
        mpz_t halfsize;
	mpz_init(halfsize);
//        mpz_ui_pow_ui(size, 2, msgsize);
//        mpz_ui_pow_ui(halfsize, 2, msgsize-1);
  
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
