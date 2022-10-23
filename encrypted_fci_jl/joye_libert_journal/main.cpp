#include "joye_libert.h"

#include <gmp.h>
//#include <iostream>
//#include <chrono>
//#include <ctime>

int main()
{
	mpz_t N, y, p;

	mpz_init(N);
	mpz_init(y);
	mpz_init(p);

	uint32_t keysize = 3072;
	uint32_t msgsize = 64;

	joye_libert_keygen(N, y, p, msgsize, keysize);

	mpz_t m1, m2, c1, c2;

	mpz_init_set_ui(m1, 999);
	mpz_init_set_ui(m2, 150);

	mpz_init(c1);
	mpz_init(c2);

	gmp_randstate_t state;
	gmp_randinit_mt(state);

	gmp_printf("m1: %Zd\n", m1);
	gmp_printf("m2: %Zd\n", m2);

	joye_libert_encrypt(c1, state, m1, y, N, msgsize);
	joye_libert_encrypt(c2, state, m2, y, N, msgsize);

	mpz_mul(c1, c1, c2);
	mpz_mod(c1, c1, N);

	mpz_t recov_m;
	mpz_init(recov_m);

	joye_libert_decrypt(recov_m, c1, p, y, msgsize);

	gmp_printf("recov_m: %Zd\n", recov_m);

	// Encrypt
	//auto start_enc = std::chrono::system_clock::now();
	/*	
	for (int i = 0; i < 100000; i++)
	{
		mpz_set_ui(m, i);
		joye_libert_encrypt(c, state, m, y, N, msgsize);
	//auto end_enc = std::chrono::system_clock::now();

	//std::chrono::duration<double> elapsed_seconds_enc = end_enc - start_enc;

	//std::cout << "Encryption time: " << elapsed_seconds_enc.count() << "s\n";

		gmp_printf("c: %Zd\n", c);
	// Recover
		mpz_t recov_m;
		mpz_init(recov_m);

	//auto start_dec = std::chrono::system_clock::now();
		joye_libert_decrypt(recov_m, c, p, y, msgsize);
	//auto end_dec = std::chrono::system_clock::now();

	//std::chrono::duration<double> elapsed_seconds_dec = end_dec - start_dec;

	//std::cout << "Decryption time: " << elapsed_seconds_dec.count() << "s\n";

		gmp_printf("recov_m: %Zd\n", recov_m);

		if ( mpz_cmp(m, recov_m) != 0 )
		{
			exit(1);
		}
	}
	*/
}
