///////////////////////////////////////////
//
//
//////////////////////////////////////////

#ifndef PPFCI_H
#define PPFCI_H

#include "muLabeledHE/mu_labeled_he.h"
#include <gmp.h>
#include <stdint.h>

void sensor_encrypt(
		mpz_t c_tr, 
		he_ct C_tr[], 
		he_ct C_P[], 
		he_ct C_Px[], 
		const mpz_t P[], 
		const mpz_t Px[], 
		const uint32_t dim, 
		const uint32_t n_sensors);

void encrypted_fci(
		mpz_t c_P0[], 
		mpz_t c_P0x0[], 
		mpz_t m_den, 
		const mpz_t c_tr[], 
		const he_ct C_tr[], 
		const he_ct C_P[], 
		const he_ct C_Px[], 
		const uint32_t dim, 
		const uint32_t n_sensors);

void decrypt(
		mpz_t P0[], 
		mpz_t P0x0[], 
		const mpz_t c_P0[], 
		const mpz_t c_x0[], 
		const uint32_t dim);

void normalize(?);

#endif
