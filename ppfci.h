///////////////////////////////////////////
//
//
//////////////////////////////////////////

#ifndef PPFCI_H
#define PPFCI_H

void sensor_encrypt(mpz_t c_tr, mpz_t C_tr, mpz_t C_P[], mpz_t C_Px[], const mpz_t P[], const mpz_t Px[], const uint32_t dim);

void encrypted_fci(mpz_t c_P0[], mpz_t c_P0x0[], mpz_t m_den, const mpz_t c_tr[],
	       	const mpz_t C_tr[], const mpz_t C_P[], const mpz_t C_Px[]);

void decrypt(mpz_t P0[], mpz_t P0x0[]);

void normalize(?);

#endif
