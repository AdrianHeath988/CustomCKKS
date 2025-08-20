#ifndef ENCODER_H
#define ENCODER_H

#include <complex>
#include <vector>

class Encoder {
public:

    Encoder(int M);
    void printVector(std::vector<std::complex<double>> vec);
    void printMatrix(std::vector<std::vector<std::complex<double>>> matrix);
    std::vector<std::complex<double>> evaluate_sigma(std::vector<std::complex<double>> coeffs);
    std::vector<std::complex<double>> evaluate_sigma_inverse(std::vector<std::complex<double>> input_vector);

private:
    int M, N;
    std::complex<double> xi; // M-th root of unity: e^(2*i*pi/M)
    std::vector<std::vector<std::complex<double>>> compute_vandermonde();
    std::vector<std::vector<std::complex<double>>> matrix_inverse(std::vector<std::vector<std::complex<double>>> matrix);
    std::vector<std::complex<double>> dotProduct(std::vector<std::vector<std::complex<double>>> matrix, std::vector<std::complex<double>> vec);
    std::complex<double> evaluatePolynomial(std::vector<std::complex<double>> coeffs, std::complex<double> x);
};

#endif // ENCODER_H