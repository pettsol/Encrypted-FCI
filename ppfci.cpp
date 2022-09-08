///////////////////////////////////////////
//
//
///////////////////////////////////////////


void sensor_encrypt(
		mpz_t c_tr,
		he_ct *C_tr,
		he_ct C_P[],
		he_ct C_Px[],
		const mpz_t P[],
		const mpz_t Px[],
		const mpz_t ski,
		const mpz_t timestep,
		const mpz_t ppaN,
		const mpz_t ppaN2,
		const uint32_t dim, 
		const uint32_t n_sensors)
{
	// Compute the trace of P
	mpz_t trace;
	mpz_init_set_ui(trace, 0);
	for (int i = 0; i < dim; i++)
	{
		// What format is the matrix stored in? If by rows:
		mpz_add(trace, trace, P[i*dim + i]);
	}

	// Do PPA encryption of the trace of P
	ppa_encrypt(c_tr, ski, trace, timestep, ppaN, ppaN2);

	// Do MU-LabHE encryption of the trace of P
	// Encrypt?
	mu_he_encrypt(C_tr, rand_state, ?, ??, trace, y, N, label, msgsize);

	// Do element-wise MU-LabHE encryption of P^{-1}
	uint32_t dim2 = dim*dim;
	for (int i = 0; i < dim2; i++)
	{
		// Encrypt
	}

	// Do element-wise MU-LabHE encryption of P^{-1}x
	for (int i = 0; i < dim; i++)
	{
		// Encrypt
	}
}

void encrypted_fci(
		mpz_t c_P0[], 
		mpz_t c_P0x0[], 
		mpz_t m_den,
	        const mpz_t sk0,
		const mpz_t timestep,
		const ppaN,
		const ppaN2,	
		const mpz_t c_tr[], 
		const he_ct C_tr[], 
		const he_ct C_P[], 
		const he_ct C_Px[], 
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
	mpz_invert(tmp, tmp, msgsize);

	// Compute the encrypted cumulative sum of traces
	he_ct sum;
	mu_he_eval_add(&sum, C_tr[0], C_tr[1], N, msgsize);
	for (int i = 2; i < n_sensors; i++)
	{
		mu_he_eval_add(&sum, &sum, C_tr[i], N, msgsize);
	}

	// Compute the weights and multiply
	uint32_t dim2 = dim*dim;
	for (int i = 0; i < n_sensors; i++)
	{
		// Compute encrypted weight for sensor system i
		

		// Perform element-wise ciphertext multiplication 
		// between weight i and matrix P[i]
		for (int i = 0; i < dim*dim; i++)
		{
			??
		}

		// Perform element-wise ciphertext multiplication
		// between weight i and vector P[i]x[i]
		for (int i = 0; i < dim; i++)
		{

		}
	}

	// Sum the weighted contributions
}
