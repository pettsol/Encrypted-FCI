///////////////////////////////////////////
//
//
//////////////////////////////////////////

#ifndef PPFCI_H
#define PPFCI_H

#include "muLabeledHE/mu_labeled_he.h"
#include "privacy_preserving_aggregation/ppa.h"
#include <gmp.h>
#include <stdint.h>

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
                const uint32_t n_sensors);

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
		const uint32_t n_sensors);

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
		const uint32_t n_sensors);

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
		const uint32_t n_sensors);

void rho(
		mpz_t out,
		const mpf_t in,
		const mpz_t gamma,
		const mpz_t ptspace);

void rho_inv(
		mpf_t out,
		const mpz_t in,
		const mpz_t gamma,
		const mpz_t ptspace);

void normalize(void);

#endif
