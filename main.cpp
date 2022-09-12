#include "ppfci.h"

#include <fstream>
#include <iostream>
#include <unistd.h>

int main()
{
	// We need to generate muLabHE keys, i.e., 1 msk and 1 mpk,
	// and N usk and upk.
	// We need to generate N ppa encryption keys
	uint32_t msgsize = 64;
	//uint32_t keysize = 2048;
	uint32_t keysize = 1024;
	uint32_t dim = 4;
	uint32_t n_sensors = 2;
	//uint32_t timesteps = 2;
	uint32_t timesteps = 1;

	gmp_randstate_t rand_state;
	gmp_randinit_mt(rand_state);

	mpz_t y, p, N, ppaN, upk[n_sensors], usk[n_sensors], sk[n_sensors+1];
	mpz_t ppaN2;
	mpz_init(ppaN2);
	
	mpz_init_set_ui(sk[0], 0);
	mpz_init(y);
	mpz_init(p);
	mpz_init(N);
	mpz_init(ppaN);
	for (uint32_t i = 0; i < n_sensors; i++)
	{
		mpz_init(upk[i]);
		mpz_init(usk[i]);
		mpz_init(sk[i+1]);
	}

	// *** RUN THE SETUP TO GENERATE THE APPROPRIATE KEYS *** //
	ppfci_setup(rand_state, y, p, N, ppaN, upk, usk, sk, msgsize, keysize, n_sensors);

	// Multiply before - it is convenient and saves a lot of time!
	mpz_mul(ppaN2, ppaN, ppaN);

	// We need to extract information from the files and store in 3 arrays	
	std::ifstream Pm_file, P_file, Px_file;
	Pm_file.open("datasets/n_sensors=2/covariance_0.csv");
	P_file.open("datasets/n_sensors=2/inv_covariance_0.csv");
	Px_file.open("datasets/n_sensors=2/inv_covariance_x_mean_0.csv");
	uint32_t count;
	count = dim*dim*n_sensors*timesteps;

	double Pm_double_array[count];
	double P_double_array[count];
	for (uint32_t i = 0; i < count; i++)
	{
		Pm_file >> Pm_double_array[i];
		P_file >> P_double_array[i];
	}

	count = count / dim;
	double Px_double_array[count];
	for (uint32_t i = 0; i < count; i++)
	{
		Px_file >> Px_double_array[i];
		//std::cout << Px_double_array[i] << std::endl;
	}

	// Map the arrays to unsigned integers
	mpz_t gamma, ptspace;
	mpz_init(gamma);
	mpz_init(ptspace);

	mpz_ui_pow_ui(ptspace, 2, msgsize);
	mpz_ui_pow_ui(gamma, 2, 10);

	mpz_t Pm[count*dim], P[count*dim], Px[count];

	mpf_t input;
	mpf_init(input);

	std::cout << "Mapping to integers\n";
	for (uint32_t i = 0; i < count*dim; i++)
	{
		mpz_init(Pm[i]);
		mpz_init(P[i]);

		mpf_set_d(input, Pm_double_array[i]);
		rho(Pm[i], input, gamma, ptspace);

		mpf_set_d(input, P_double_array[i]);
		rho(P[i], input, gamma, ptspace);
	}

	for (uint32_t i = 0; i < count; i++)
	{
		mpz_init(Px[i]);
		mpf_set_d(input, Px_double_array[i]);
		rho(Px[i], input, gamma, ptspace);
	}


	// For each timestep we need to encrypt the sensor output from each sensor,
	// fuse it, decrypt the result, and store in a new text file
	mpz_t timestep, label, label_start;
	mpz_init(timestep);
	mpz_init(label_start);
	mpz_init_set_ui(label, 1);

	// Variables holding the encrypted sensor data
	mpz_t c_trace[n_sensors];
	he_ct C_P[dim*dim*n_sensors];
	he_ct C_Px[dim*n_sensors];
	he_ct C_trace[n_sensors];

	// Variables holding the fused encrypted data and the sum of traces
	mpz_t c_P0[dim*dim];
	mpz_t c_P0x0[dim];
	mpz_t m_den;

	// Variables holding the decrypted fused data
	mpz_t P0[dim*dim];
	mpz_t P0x0[dim];

	// Variables to hold the final fused data
	double P0_matrix[dim*dim] = {0};
	double P0x0_vector[dim] = {0};

	// Initialize the arrays
	for (uint32_t i = 0; i < dim*dim; i++)
	{
		mpz_init(c_P0[i]);
		mpz_init(P0[i]);
	}

	for (uint32_t i = 0; i < dim; i++)
	{
		mpz_init(c_P0x0[i]);
		mpz_init(P0x0[i]);
	}
	
	// Initialize the sum of traces
	mpz_init(m_den);

	// We need to write the fused P0 and P0x0 to files
	std::ofstream P0_file("fused_inv_covariance_0.csv");
	std::ofstream P0x0_file("fused_inv_covariance_x_mean_0.csv");

	for (uint32_t t = 0; t < timesteps; t++)
	{
		// For each timestep, run 1 iteration
		std::cout << "Fusing timestep: " << t << std::endl;
		mpz_set(label_start, label);
		mpz_set_ui(timestep, t);

		// First encrypt the data from each sensor using the
		// sensor-specific secret key
		for (uint32_t i = 0; i < n_sensors; i++)
		{
			std::cout << "Encrypting sensor " << i << std::endl;
			mpz_init(c_trace[i]);
			ppfci_sensor_encrypt(c_trace[i], label, rand_state, &C_trace[i], 
					&C_P[i*(dim*dim)], &C_Px[i*dim], &Pm[i*(dim*dim) + t*(dim*dim*n_sensors)], 
					&P[i*(dim*dim) + t*(dim*dim*n_sensors)], &Px[i*dim + t*(dim*n_sensors)], 
					usk[i], y, N, sk[i+1], timestep, ppaN, ppaN2, msgsize, dim, n_sensors);
		}

		/*
		// DEBUG! Check that the sk0 is correct
		mpz_t test;
		mpz_init(test);
		mpz_add(test, sk[0], sk[1]);
		mpz_add(test, test, sk[2]);
		gmp_printf("Test: %Zd\n", test);

		// DEBUG! Check the expected sum of traces
		mpz_t true_sum, tmp_sum;
		mpz_init_set_ui(true_sum, 0);
		mpz_init(tmp_sum);
		for (uint32_t i = 0; i < n_sensors; i++)
		{
			mpz_set_ui(tmp_sum, 0);
			for (uint32_t j = 0; j < dim; j++)
			{
				mpz_add(tmp_sum, tmp_sum, Pm[j + j*dim + i*dim*dim]);
			}
			mpz_add(true_sum, true_sum, tmp_sum);
		}
		gmp_printf("Expected sum: %Zd\n", true_sum); */
		
		// Proceed by fusing the encrypted data
		std::cout << "Fusing timestep " << t << std::endl;
		ppfci_encrypted_fusion(c_P0, c_P0x0, m_den, rand_state, sk[0], timestep, y, N, ppaN,
					ppaN2, c_trace, C_trace, C_P, C_Px, msgsize, dim, n_sensors);

		/*// Check the sum of traces variable
		gmp_printf("sk[0] = %Zd\n", sk[0]);
		gmp_printf("Sum of traces: %Zd\n", m_den);*/

		// Decrypt the fused data
		std::cout << "Decrypting timestep " << t << std::endl;
		ppfci_decrypt(P0, P0x0, label_start, c_P0, c_P0x0, m_den, upk, p, y, msgsize, dim, n_sensors);

		// Map back to floating point numbers and normalize
		ppfci_normalize(P0_matrix, P0x0_vector, P0, P0x0, m_den, ptspace, gamma, dim);

		// Write the fused data to file
		// We write the matrix first
		for (uint32_t i = 0; i < dim*dim; i++)
		{
			std::cout << "Writing to file\n";
			P0_file << P0_matrix[i];
			P0_file << " ";
		}
		// Timesteps separated by newline
		P0_file << std::endl;

		// Then we write the vector
		for (uint32_t i = 0; i < dim; i++)
		{
			P0x0_file << P0x0_vector[i];
			P0x0_file << " ";
		}
		P0x0_file << std::endl;
	}

	// Close the files
	P0_file.close();
	P0x0_file.close();
}
