#include "encryptor.h"
#include <random>
#include <stdexcept>
#include <numeric>

// Helper function for correct modular arithmetic on potentially negative numbers
long long custom_mod(long long value, long long modulus) {
    return (value % modulus + modulus) % modulus;
}

Encryptor::Encryptor(size_t poly_degree, int64_t modulus)
    : poly_degree_N(poly_degree), modulus_q(modulus) {
    // Check if N is a power of 2
    if (poly_degree_N == 0 || (poly_degree_N & (poly_degree_N - 1)) != 0) {
        throw std::invalid_argument("Polynomial degree N must be a power of 2.");
    }
}

void Encryptor::generate_keys() {
    //Generate the secret key 's' from a noise distribution
    secret_key_s = sample_noise_poly();

    // Generate the public key p = (b, a)
    // Sample a uniform polynomial 'a'
    Polynomial a = sample_uniform_poly();
    // Sample a noise polynomial 'e'
    Polynomial e = sample_noise_poly();

    // Calculate -a*s
    Polynomial a_s = poly_mult(a, secret_key_s);
    Polynomial neg_a_s = poly_neg(a_s);

    // Calculate b = -a*s + e
    Polynomial b = poly_add(neg_a_s, e);

    public_key_p = std::make_pair(b, a);
}

Ciphertext Encryptor::encrypt(const Polynomial& plain_poly) {
    // The public key is p = (b, a)
    const Polynomial& b = public_key_p.first;
    const Polynomial& a = public_key_p.second;

    // The ciphertext is c = (u + b, a)
    Polynomial c0 = poly_add(plain_poly, b);
    
    return std::make_pair(c0, a);
}

Polynomial Encryptor::decrypt(const Ciphertext& cipher_poly) {
    const Polynomial& c0 = cipher_poly.first;
    const Polynomial& c1 = cipher_poly.second;

    // Decryption is calculated as c0 + c1*s
    Polynomial c1_s = poly_mult(c1, secret_key_s);
    Polynomial decrypted_poly = poly_add(c0, c1_s);

    // The result is an approximation u + e
    return decrypted_poly;
}

const PublicKey& Encryptor::get_public_key() const {
    return public_key_p;
}


Polynomial Encryptor::sample_uniform_poly() const {
    // Set up a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> distrib(0, modulus_q - 1);

    Polynomial p(poly_degree_N);
    for (size_t i = 0; i < poly_degree_N; ++i) {
        p[i] = distrib(gen);
    }
    return p;
}

Polynomial Encryptor::sample_noise_poly() const {
    // Use a discrete Gaussian distribution for sampling small noise coefficients
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 3.2);

    Polynomial p(poly_degree_N);
    for (size_t i = 0; i < poly_degree_N; ++i) {
        p[i] = static_cast<int64_t>(std::round(d(gen)));
    }
    return p;
}

Polynomial Encryptor::poly_add(const Polynomial& a, const Polynomial& b) const {
    Polynomial result(poly_degree_N);
    for (size_t i = 0; i < poly_degree_N; ++i) {
        result[i] = custom_mod(a[i] + b[i], modulus_q);
    }
    return result;
}

Polynomial Encryptor::poly_neg(const Polynomial& a) const {
    Polynomial result(poly_degree_N);
    for (size_t i = 0; i < poly_degree_N; ++i) {
        result[i] = custom_mod(-a[i], modulus_q);
    }
    return result;
}

Polynomial Encryptor::poly_mult(const Polynomial& a, const Polynomial& b) const {
    // Implements polynomial multiplication in the ring Z_q[X] / (X^N + 1)
    // This is a simpler, but less efficient, O(N^2) implementation
    // A faster and professional implementation would be to use NTTT or FTT
    Polynomial result(poly_degree_N, 0);
    for (size_t i = 0; i < poly_degree_N; ++i) {
        for (size_t j = 0; j < poly_degree_N; ++j) {
            if (i + j < poly_degree_N) {
                // Standard term
                result[i + j] = custom_mod(result[i + j] + a[i] * b[j], modulus_q);
            } else {
                // Term wraps around due to X^N = -1
                // a_i * b_j * X^(i+j) = a_i * b_j * X^(i+j-N) * X^N = -a_i * b_j * X^(i+j-N)
                size_t wrapped_index = i + j - poly_degree_N;
                result[wrapped_index] = custom_mod(result[wrapped_index] - a[i] * b[j], modulus_q);
            }
        }
    }
    return result;
}