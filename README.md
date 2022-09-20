This repository contains an implementation of a privacy-preserving fast covariance intersection algorithm that fuses the state estimates and covariance matrices of N sensor systems in encrypted form. The scheme utilizes the labeled homomorphic encryption scheme by Manuel Barbosa, Dario Catalano, and Dario Fiore and the privacy-preserving aggregation scheme designed by Marc Joye and Benoit Libert.

The implementation and sample program can be compiled using the g++;

g++ -Werror=pedantic -std=c++11 main.cpp ppfci.cpp joye\_libert\_journal/joye\_libert.cpp muLabeledHE/mu\_labeled\_he.cpp privacy\_preserving\_aggregation/ppa.cpp privacy\_preserving\_aggregation/SHA-256/sha-256.cpp -o main -lgmp


NB! Compilation requires the GNU Multiple Precision Arithmetic Library to be installed. More information available here: https://gmplib.org/manual/Introduction-to-GMP

References:
Labeled homomorphic encryption: https://eprint.iacr.org/2017/326.pdf
Privacy-preserving aggregation: https://www.ifca.ai/fc13/proc/3-3.pdf
