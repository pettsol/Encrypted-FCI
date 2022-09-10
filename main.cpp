#include "ppfci.h"

#include <fstream>
#include <iostream>

int main()
{
	// We need to generate muLabHE keys, i.e., 1 msk and 1 mpk,
	// and N usk and upk.
	// We need to generate N ppa encryption keys
	uint32_t msgsize = 64;
	uint32_t keysize = 2048;
	uint32_t dim = 4;
	uint32_t n_sensors = 2;
	uint32_t timesteps = 150;

	gmp_randstate_t rand_state;
	gmp_randinit_mt(rand_state);

	mpz_t y, p, N, ppaN, upk[n_sensors], usk[n_sensors], sk[n_sensors+1];
	mpz_init(sk[0]);
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

	ppfci_setup(rand_state, y, p, N, ppaN, upk, usk, sk, msgsize, keysize, n_sensors);

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
		std::cout << Px_double_array[i] << std::endl;
	}

	// Map the arrays to unsigned integers

	// For each timestep we need to encrypt the sensor output from each sensor,
	// fuse it, decrypt the result, and store in a new text file	
}
