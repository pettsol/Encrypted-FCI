////////////////////////////////////
//
//
////////////////////////////////////

#include "ppa.h"
#include <gmp.h>
#include <iostream>

int main()
{
	uint32_t n_sensors = 3;
	uint32_t keysize = 4096;

	std::cout << "Hello World\n";

#ifdef DEBUG
	gmp_printf("Declaring keys\n");
#endif

	mpz_t sk[n_sensors+1];
#ifdef DEBUG
	gmp_printf("Initializing keys\n");
#endif
	for (int i = 0; i < n_sensors+1; i++)
	{
		mpz_init(sk[i]);
	}

	mpz_t N;
	mpz_init(N);

	ppa_setup(N, sk, keysize, n_sensors);

	mpz_t N2;
	mpz_init(N2);
	mpz_mul(N2, N, N);

	// TEST ADDING 3 SENSORS
	mpz_t x1, x2, x3;
	mpz_init_set_ui(x1, 1);
	mpz_init_set_ui(x2, 2);
	mpz_init_set_ui(x3, 3);

	mpz_t x[n_sensors];
	mpz_t c[n_sensors];

	for (int i = 0; i < n_sensors; i++)
	{
		mpz_init_set_ui(x[i], i+100);
		mpz_init(c[i]);
	}

	mpz_t t;
	mpz_init_set_ui(t, 1);

	// Encrypt
	for (int i = 0; i < n_sensors; i++)
	{
		ppa_encrypt(c[i], sk[i+1], x[i], t, N, N2);
	}

#ifdef DEBUG
	std::cout << "Finished encrypting\n";
#endif

	std::cout << "Aggregating and decrypting\n";
	// Aggregate and decrypt
	mpz_t sum;
	mpz_init(sum);
	ppa_aggregate_decrypt(sum, sk[0], t, c, N, N2, n_sensors);

#ifdef DEBUG
	gmp_printf("The sum is: %Zd\n", sum);
#endif
}
