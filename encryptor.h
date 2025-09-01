#ifndef ENCRYPTOR_H
#define ENCRYPTOR_H

#include <vector>
#include <cstdint>
#include <utility>

// Type alias for representing a polynomial as a vector of its coefficients.
using Polynomial = std::vector<int64_t>;

// Type alias for the secret key, which is a single polynomial 's'.
using SecretKey = Polynomial;

// Type alias for the public key, a pair of polynomials (b, a) where b = -a*s + e.
using PublicKey = std::pair<Polynomial, Polynomial>;

// Type alias for a ciphertext, which is a pair of polynomials (c0, c1)
using Ciphertext = std::pair<Polynomial, Polynomial>;


class Encryptor {
public:
    /**
     * @brief Constructs an Encryptor instance.
     * @param poly_degree The degree 'N' of the polynomial ring Z_q[X] / (X^N + 1). [cite_start]This is a primary security parameter[cite: 51].
     * @param modulus The coefficient modulus 'q' for the polynomial ring.
     */
    Encryptor(size_t poly_degree, int64_t modulus);

    /**
     * Generates the secret and public keys.
     */
    void generate_keys();

    /**
     * @brief Encrypts a plaintext polynomial message.
     *
     * The encryption of a message polynomial 'u' is performed using the public key p = (b, a).
     * The resulting ciphertext is c = (u + b, a). The message 'u' is masked
     * @param plain_poly The plaintext polynomial 'u' to be encrypted.
     * @return A Ciphertext pair (c0, c1) representing the encrypted message.
     */
    Ciphertext encrypt(const Polynomial& plain_poly);

    /**
     * @brief Decrypts a ciphertext to retrieve the approximate original message.
     *
     * The decryption of a ciphertext c = (c0, c1) is performed using the secret key 's'.
     * The approximate message is recovered by computing c0 + c1*s.
     * This calculation results in u + e, 
     * @param cipher_poly The ciphertext to be decrypted.
     * @return A Polynomial that is an approximation of the original plaintext message.
     */
    Polynomial decrypt(const Ciphertext& cipher_poly);

    /**
     * @brief Retrieves the public key.
     *
     * @return A constant reference to the PublicKey.
     */
    const PublicKey& get_public_key() const;

private:
    size_t poly_degree_N;
    int64_t modulus_q;
    SecretKey secret_key_s;
    PublicKey public_key_p;

    // --- Private Helper Functions for Polynomial Arithmetic ---
    // These functions would be defined in a corresponding .cpp file to implement
    // the necessary modular arithmetic for polynomials in the Z_q[X] / (X^N + 1) ring.

    /**
     * @brief Samples a polynomial with coefficients chosen uniformly from Z_q.
     */
    Polynomial sample_uniform_poly() const;

    /**
     * @brief Samples a small polynomial from a discrete Gaussian distribution for noise.
     */
    Polynomial sample_noise_poly() const;

    /**
     * @brief Performs polynomial multiplication in the ring Z_q[X] / (X^N + 1).
     * This is more efficient than matrix multiplication in LWE, with a complexity
     * [cite_start]of O(N*log(N)) using transforms like the DFT[cite: 59].
     */
    Polynomial poly_mult(const Polynomial& a, const Polynomial& b) const;

    /**
     * @brief Performs polynomial addition in the ring.
     */
    Polynomial poly_add(const Polynomial& a, const Polynomial& b) const;
    
    /**
     * @brief Performs polynomial negation in the ring.
     */
    Polynomial poly_neg(const Polynomial& a) const;
};

#endif // ENCRYPTOR_H