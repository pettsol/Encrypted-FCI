//////////////////////////////////
// This implementation has been //
// placed in the public domain  //
// by			        //
//			      	//
// Petter Solnoer - 1/11/2021 	//
//////////////////////////////////

#include "joye_libert.h"
#include "misc.h"

#include <iostream>
//#include <NTL/ZZ_pXFactoring.h>
//#include <NTL/ZZ_pEX.h>
//#include <NTL/GF2X.h>
//#include <NTL/GF2E.h>
//#include <NTL/GF2XFactoring.h>

void jl_key_gen(mpz_t N, mpz_t y, mpz_t k, mpz_t p)
{
	//** Generate primes p, q congruent to 1 mod 2^{k}. **
	// Start by initializing a prng.
	gmp_randstate_t state;
	gmp_randinit_mt(state);

	// Variable to hold q
	mpz_t q;
	mpz_init(q);
	// Variables to hold test primes
	mpz_t p_test, q_test;
	mpz_init(p_test);
	mpz_init(q_test);
	
	// Variable to hold pseudo-safeness test divisor
	// and pseudo_test assume 32-bit messages.
	mpz_t pseudo_safeness_divisor, pseudo_test_p, pseudo_test_q;
	mpz_init(pseudo_safeness_divisor);
	mpz_init(pseudo_test_p);
	mpz_init(pseudo_test_q);
	
	// Set divisor to 2^{32}
	mpz_ui_pow_ui(pseudo_safeness_divisor, 2, 32);

	// Set p to random value with 1024 bits
	mpz_urandomb(p_test, state, 1024);
	
	// Find the modulus
	mpz_t mod;
	mpz_init(mod);
	mpz_mod(mod, p_test, pseudo_safeness_divisor);

	mpz_sub(p_test, p_test, mod);
	mpz_add_ui(p_test, p_test, 1);

	// We now have a number that is congruent to 1, as desired.
	// Now we must check if it is prime!

	//mpz_nextprime(p_test, p_test);
	// Start by finding appropriate p
	mpz_t one;
	mpz_init_set_ui(one, 1);

	std::cout << "Find appropriate p\n";
	mpz_t t;
	mpz_init(t);
	while(true)
	{
		// check if prime
		//mpz_mod(t, p_test, pseudo_safeness_divisor);
		//gmp_printf("%Zd is the congruence.\n", t);
		//if (mpz_congruent_2exp_p(p_test, one, 32))
		if (mpz_probab_prime_p(p_test, 30))
		{
			// check for pseudo-safeness
			// condition: p_test - 1 / 2^{k} is also prime
			mpz_t tmp;
			mpz_init(tmp);
			mpz_sub_ui(tmp, p_test, 1);

			mpz_divexact(pseudo_test_p, tmp, pseudo_safeness_divisor);
			if (mpz_probab_prime_p(pseudo_test_p, 30))
			{
				// We have probably found a good value
				// for p
				gmp_printf("Good value for p: %Zd\n", p_test);
				mpz_set(p, p_test);
				break;
			}	
		}
		// Check next prime, probabilistic, chance of composite passing
		// is extremely small
		mpz_add(p_test, p_test, pseudo_safeness_divisor);
	}

	// Proceed to find appropriate q
	std::cout << "Find appropriate q\n";
	mpz_add(q_test, p_test, pseudo_safeness_divisor);
	//mpz_nextprime(q_test, p_test);
	while(true)
	{
		// check if prime
		if (mpz_probab_prime_p(q_test, 30))
		{
			// Check for pseudo-safeness
			// condition: q_test - 1 / 2^{k} is also prime
			mpz_t tmp;
			mpz_init(tmp);
			mpz_sub_ui(tmp, q_test, 1);

			mpz_cdiv_q(pseudo_test_q, tmp, pseudo_safeness_divisor);
			if (mpz_probab_prime_p(pseudo_test_q, 30))
			{
				// We have probably found a good value
				// for q
				gmp_printf("Good value for q: %Zd\n", q_test);
				mpz_set(q, q_test);
				break;
			}	
		}
		// Check next prime, probabilistic, chance of composite passing
		// is extremely small
		mpz_add(q_test, q_test, pseudo_safeness_divisor);
	}

	// Set N = p*q
	mpz_mul(N, p, q);

	// Find y, using algorithm 2 from Joye-Libert paper
	// Assume that the authors mean an arbitrary primitive 2^{th}
	// root of unity, i.e., a generator for the cyclic group under multiplication.
	
	// First select elements from the multiplicative groups Zp* and Zq*
	//
	// Select a generator for the multiplicative group 2^{k}
	// Iterate through primes and divide (?)
	//std::set<mpz_t> prime_factors;
	//std::set<mpz_t>::iterator it;
	// Array to hold prime factors
	// and number of prime factors
	/* ###################################### CONFERENCE PAPER SHIT #############
	mpz_t prime_factors[1000];
	int num_prime_factors = 0;

	mpz_t i, j;
	mpz_init(i);
	mpz_init(j);

	mpz_set_ui(i, 2);
	mpz_sub_ui(j, pseudo_safeness_divisor, 1);
	std::cout << "Find prime factors\n";
	while(true)
	{
		if(mpz_divisible_p(j, i))
		{
			gmp_printf("%Zd is a prime factor\n", i);
			mpz_set(prime_factors[num_prime_factors], i);
			num_prime_factors++;
			mpz_divexact(j, j, i);
			gmp_printf("j is now: %Zd\n", j);
			if(!mpz_cmp(j,one))
			{
				std::cout << "Breaking out\n";
				break;
			}
		} else {
			mpz_nextprime(i, i);
		}
	}

	std::cout << "Find a generator\n";


	// We need to 2^{k} is an extension field with characteristic 2.
	// So we need to use polynomials to represent the elements of
	// the multiplicative group. We use the NTL library by
	// Victor Shoup to do this.
	
	// We begin by defining the underlying field GF(2):
	//NTL::ZZ_p::init(NTL::ZZ(2));
	// Declare a polynomial and initialize with degree k=32
	
	NTL::GF2X irred_polynomial;
	

	
	NTL::BuildIrred(irred_polynomial, 32);
	NTL::GF2XModulus pol_mod(irred_polynomial);
	NTL::GF2E::init(irred_polynomial);
	// Declare and initialize a polynomial used
	// to search for generator.
	NTL::GF2X g_pol;
	NTL::GF2X f_pol;
	NTL::BuildIrred(f_pol, 30);

	// Variable to check if generator is found
	int generator_found = 0;

	mpz_t exponent, flag, g;
	mpz_init(exponent);
	mpz_init(flag);
	mpz_init(g);

	NTL::ZZ ntl_exponent;
	NTL::GF2X ntl_flag, g_pol_min;
	*/ // #### CONFERENCE PAPER SHIT

	/*mpz_t flag, exponent, g;
	mpz_init(flag);
	mpz_init(exponent);
	mpz_init_set_ui(g, 1);*/
	// j should contain 2^{k} - 1
	//
	/* #### CONFERENCE PAPER SHIT
	mpz_sub_ui(j, pseudo_safeness_divisor, 1);
	while(true)
	{
		generator_found = 1;
		// Pick a value in the set
		//mpz_add_ui(g, g, 1);
		g_pol = NTL::GF2X::zero();
		while (NTL::IsZero(g_pol))
		{
			// A primitive polynomial must be
			// irreducible
			NTL::BuildRandomIrred(g_pol, f_pol);
		}
		std::cout << "pol_mod: " << irred_polynomial << std::endl;
		std::cout << "g_pol: " << g_pol << std::endl;

		// Test if value is a generator
		for (int i = 0; i < num_prime_factors; ++i)
		{
			mpz_divexact(exponent, j, prime_factors[i]);
			gmp_printf("Exponent being tested: %Zd\n", exponent);

			// Convert exponent to MPZ::ZZ
			MPZToZZ(&ntl_exponent, exponent);

			//mpz_powm(flag, g, exponent, pseudo_safeness_divisor);
			NTL::PowerMod(ntl_flag, g_pol, ntl_exponent, pol_mod);

			std::cout << "Flag: " << ntl_flag << std::endl;

			unsigned char buffer[4];
			unsigned int integer = 0;		
			MyBytesFromGF2X(buffer, ntl_flag, 4);
			std::memcpy(&integer, buffer, 4);
			//ZZToMPZ(flag, &ntl_flag);
			
			std::cout << "Integer: " << integer << std::endl;

			//gmp_printf("Flag: %Zd\n", flag);
			//if (!mpz_cmp(flag, one))
			if(integer == 1)
			{
				// Not a generator
				//gmp_printf("%Zd is not a generator!\n", g);
				generator_found = 0;
				//exit(1);
				break;
			}
		}
		if(generator_found)
		{
			std::cout << "g_pol before: " << g_pol << std::endl;
			unsigned char buffer[4];
			unsigned int g_int = 0;
			MyBytesFromGF2X(buffer, g_pol, 4);
			std::memcpy(&g_int, buffer, 4);
			mpz_set_ui(g, g_int);
			break;
		}
	}

	gmp_printf("The generator found: %Zd\n", g);
	
	// We have found a generator! Now we complete algorithm 2 from Joye-Libert paper
	// Pick yp and yq at random from the multiplicative cyclic groups of Zp and Zq
	mpz_t yp, yq;
	mpz_init(yp);
	mpz_init(yq);

	while(true)
	{
		mpz_urandomm(yp, state, p);
		mpz_urandomm(yq, state, q);
		if ( mpz_cmpabs_ui(yp, 0) && mpz_cmpabs_ui(yq, 0) )
		{
			break;
		}
	}
	// Check legendre symbol, update using generator
	int legendre_p, legendre_q;

	gmp_printf("yp: %Zd\n", yp);
	gmp_printf("yq: %Zd\n", yq);

	int test_leg;

	legendre_p = mpz_legendre(yp, p);
	legendre_q = mpz_legendre(yq, q);
	std::cout << "Legendre yp/p pre-check: " << legendre_p << std::endl;
	std::cout << "Legendre yq/q pre-check: " << legendre_q << std::endl;

	std::cout << "Check legendre symbols\n";
	
	if (legendre_p == 1)
	{
		// Update yp using generator
		mpz_mul(yp, yp, g);
		mpz_mod(yp, yp, p);
		legendre_p = mpz_legendre(yp, p);
		std::cout << "Legendre yp/p post-check: " << legendre_p << std::endl;
	}
	if (legendre_q == 1)
	{
		// Update yq
		mpz_mul(yq, yq, g);
		mpz_mod(yq, yq, q);
		legendre_q = mpz_legendre(yq, q);
		std::cout << "Legendre yq/q post-check: " << legendre_q << std::endl;
	}

	std::cout << "Legendre symbols are good\n";
	// Set y
	mpz_t p_inv;
	mpz_init(p_inv);
	//mpz_invert(p_inv, p, q);
	mpz_invert(p_inv, p, N);

	mpz_sub(y, yq, yp);
	mpz_mod(y, y, q);
	mpz_mul(y, y, p_inv);
	//mpz_mod(y, y, q);
	mpz_mul(y, y, p);
	mpz_add(y, y, yp);

	gmp_printf("N: %Zd\n", N);
	gmp_printf("The y found: %Zd\n", y);

	// Check if desirable properties for y has been attained:
	// Check jacibo first:
	int jacobi_n;
	*/ // ################### CONFERENCE PAPER SHIT
	// Generate a random y in Z_{N}
	mpz_set_ui(y, 0);
	while ((mpz_legendre(y, p) != -1) & (mpz_legendre(y,q) != -1))
	{
		// Generate random y in Z_N
		mpz_urandomm(y, state, N);
	}
	int jacobi_n = mpz_jacobi(y, N);
	int legendre_p = mpz_legendre(y, p);
	int legendre_q = mpz_legendre(y, q);

	std::cout << "Jacobi: " << jacobi_n << std::endl;
	std::cout << "Legendre y/p final: " << legendre_p << std::endl;
	std::cout << "Legendre y/q final: " << legendre_q << std::endl;
}

void jl_encrypt(mpz_t *c, mpz_t *m, mpz_t *y)
{

}

void jl_decrypt(mpz_t *m, mpz_t *c, mpz_t *p)
{
}
