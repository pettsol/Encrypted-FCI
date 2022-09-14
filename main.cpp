#include "ppfci.h"

#include <fstream>
#include <iostream>
#include <unistd.h>

int main()
{
	// We need to generate muLabHE keys, i.e., 1 msk and 1 mpk,
	// and N usk and upk.
	// We need to generate N ppa encryption keys
	uint32_t msgsize = 128;
	uint32_t keysize = 1024;
	uint32_t dim = 4;
	//uint32_t n_sensors = 6;
	uint32_t n_datasets = 100;
	uint32_t timesteps = 150;

	uint32_t sensor_array[] = {2, 4, 6, 8, 10, 15};
	size_t total(sizeof(sensor_array) / sizeof(uint32_t));

	gmp_randstate_t rand_state;
	gmp_randinit_mt(rand_state);

	// WE ITERATE OVER ALL NUMBERS OF SENSORS
	for (uint32_t number = 0; number < total; number++)
	{

	uint32_t n_sensors = sensor_array[number];

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
	std::cout << "*** RUNNING PPFCI SETUP ***" << std::endl;
	ppfci_setup(rand_state, y, p, N, ppaN, upk, usk, sk, msgsize, keysize, n_sensors);

	// Multiply before - it is convenient and saves a lot of time!
	mpz_mul(ppaN2, ppaN, ppaN);

	// WE ITERATE THROUGH ALL DATASETS	
	for (uint32_t dataset = 0; dataset < n_datasets; dataset++)
	{

	// We need to extract information from the files and store in 3 arrays
	char buf[100];
        sprintf(buf, "datasets/n_sensors=%u/covariance_%u.csv", sensor_array[number], dataset);
	std::string covariance_string = buf;
	sprintf(buf, "datasets/n_sensors=%u/inv_covariance_%u.csv", sensor_array[number], dataset);
	std::string inv_covariance_string = buf;
	sprintf(buf, "datasets/n_sensors=%u/inv_covariance_x_mean_%u.csv", sensor_array[number], dataset);
	std::string inv_covariance_x_mean_string = buf;

	std::ifstream Pm_file, P_file, Px_file;
	Pm_file.open(covariance_string);
	P_file.open(inv_covariance_string);
	Px_file.open(inv_covariance_x_mean_string);
	
	// We need to write the fused P0 and P0x0 to files
	sprintf(buf, "output/n_sensors=%u/inv_covariance_%u.csv", sensor_array[number], dataset);
	std::string inv_covariance_string_out = buf;
	sprintf(buf, "output/n_sensors=%u/inv_covariance_x_mean_%u.csv", sensor_array[number], dataset);
	std::string inv_covariance_x_mean_string_out = buf;

	std::ofstream P0_file(inv_covariance_string_out);
	std::ofstream P0x0_file(inv_covariance_x_mean_string_out);

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
	}

	// Map the arrays to unsigned integers
	mpz_t gamma, ptspace;
	mpz_init(gamma);
	mpz_init(ptspace);

	mpz_ui_pow_ui(ptspace, 2, msgsize);
	mpz_ui_pow_ui(gamma, 2, 20);

	mpz_t Pm[count*dim], P[count*dim], Px[count];

	mpf_t input;
	mpf_init(input);

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

	for (uint32_t t = 0; t < timesteps; t++)
	{
		// For each timestep, run 1 iteration
		mpz_set(label_start, label);
		mpz_set_ui(timestep, t);

		// First encrypt the data from each sensor using the
		// sensor-specific secret key
		for (uint32_t i = 0; i < n_sensors; i++)
		{
			mpz_init(c_trace[i]);
			ppfci_sensor_encrypt(c_trace[i], label, rand_state, &C_trace[i], 
					&C_P[i*(dim*dim)], &C_Px[i*dim], &Pm[i*(dim*dim) + t*(dim*dim*n_sensors)], 
					&P[i*(dim*dim) + t*(dim*dim*n_sensors)], &Px[i*dim + t*(dim*n_sensors)], 
					usk[i], y, N, sk[i+1], timestep, ppaN, ppaN2, msgsize, dim, n_sensors);
		}

		// Proceed by fusing the encrypted data
		std::cout << "*** FUSING ENCRYPTED DATA FROM " << sensor_array[number]  << " SENSORS AT TIMESTEP " << t+1 << " IN DATASET " << dataset+1 << " ***" << std::endl;
		ppfci_encrypted_fusion(c_P0, c_P0x0, m_den, rand_state, sk[0], timestep, y, N, ppaN,
					ppaN2, c_trace, C_trace, C_P, C_Px, msgsize, dim, n_sensors);

		// Decrypt the fused data
		ppfci_decrypt(P0, P0x0, label_start, c_P0, c_P0x0, m_den, upk, p, y, msgsize, dim, n_sensors);

		// Map back to floating point numbers and normalize
		ppfci_normalize(P0_matrix, P0x0_vector, P0, P0x0, m_den, ptspace, gamma, dim);

		// Write the fused data to file
		// We write the matrix first
		for (uint32_t i = 0; i < dim*dim; i++)
		{
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
	}
	std::cout << " *** FUSION OF ALL DATASETS COMPLETE ***"  << std::endl;
}
