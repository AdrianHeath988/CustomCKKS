#include <stdint.h>
#include "encoder.h"
#include <iostream>
#include <cmath>
#include "encryptor.h"

void printPolynomial(const Polynomial& p, const std::string& title) {
    std::cout << "----- " << title << " -----" << std::endl;
    std::cout << "[ ";
    for (size_t i = 0; i < p.size(); ++i) {
        std::cout << p[i] << (i == p.size() - 1 ? "" : ", ");
    }
    std::cout << " ]" << std::endl;
}

int main() {
    size_t poly_degree = 8;
    uint64_t scale = 128; 
    int64_t modulus = 8192;
    Encoder myEncoder(poly_degree, scale);
    std::vector<std::complex<double>> input;
    input.push_back(std::complex<double>(-5, 0));
    input.push_back(std::complex<double>(1, 0));

    std::cout << "----- Input -----" << std::endl;
    myEncoder.printVector(input);
    std::vector<std::complex<double>> complex_poly = myEncoder.encode(input);
    Polynomial plaintext_poly(poly_degree);
    for (size_t i = 0; i < complex_poly.size(); ++i) {
        plaintext_poly[i] = static_cast<long>(std::round(complex_poly[i].real()));
    }
    printPolynomial(plaintext_poly, "Plaintext Polynomial (Encoded & Rounded)");

    std::vector<std::complex<double>> output_no_encrypt(poly_degree);
    std::vector<std::complex<double>> pt_no_encrypt(poly_degree);
    for (size_t i = 0; i < pt_no_encrypt.size(); ++i) {
        pt_no_encrypt[i] = std::complex<double>(static_cast<double>(plaintext_poly[i]), 0.0);
    }
    output_no_encrypt = myEncoder.decode(pt_no_encrypt);
    std::cout << "----- Output (Decoded, not Encrypted) -----" << std::endl;
    myEncoder.printVector(output_no_encrypt);

    Encryptor myEncryptor(poly_degree, modulus);
    myEncryptor.generate_keys();
    
    Ciphertext encrypted_data = myEncryptor.encrypt(plaintext_poly);
    std::cout << "\n--- Polynomial has been encrypted. ---" << std::endl;
    Polynomial decrypted_poly = myEncryptor.decrypt(encrypted_data);
    std::cout << "--- Ciphertext has been decrypted. ---\n" << std::endl;
    printPolynomial(decrypted_poly, "Decrypted Polynomial");

    for (long& val : decrypted_poly) {
        if (val > modulus / 2) {
            val -= modulus;
        }
    }
    printPolynomial(decrypted_poly, "Decrypted Polynomial (Corrected for Sign)");


    std::vector<std::complex<double>> poly_to_decode(poly_degree);
    for (size_t i = 0; i < decrypted_poly.size(); ++i) {
        poly_to_decode[i] = std::complex<double>(static_cast<double>(decrypted_poly[i]), 0.0);
    }
    std::vector<std::complex<double>> output = myEncoder.decode(poly_to_decode);
    std::cout << "----- Output (Decoded) -----" << std::endl;
    myEncoder.printVector(output);

    return 0;
}
