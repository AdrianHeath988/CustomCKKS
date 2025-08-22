#ifndef ENCODER_H
#define ENCODER_H

#include <complex>
#include <vector>

class Encoder {
public:

    Encoder(int M);
    void printVector(std::vector<std::complex<double>> vec);
    void printMatrix(std::vector<std::vector<std::complex<double>>> matrix);
    std::vector<std::complex<double>> decode(std::vector<std::complex<double>> coeffs);
    std::vector<std::complex<double>> encode(std::vector<std::complex<double>> input_vector);
private:
    int M, N;
    std::complex<double> xi; // M-th root of unity: e^(2*i*pi/M)
    std::vector<std::complex<double>> evaluate_sigma(std::vector<std::complex<double>> coeffs);
    std::vector<std::complex<double>> evaluate_sigma_inverse(std::vector<std::complex<double>> input_vector);
    std::vector<std::vector<std::complex<double>>> compute_vandermonde();
    std::vector<std::vector<std::complex<double>>> matrix_inverse(std::vector<std::vector<std::complex<double>>> matrix);
    std::vector<std::complex<double>> dotProduct(std::vector<std::vector<std::complex<double>>> matrix, std::vector<std::complex<double>> vec);
    std::complex<double> vectorDot(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2);
    std::complex<double> evaluatePolynomial(std::vector<std::complex<double>> coeffs, std::complex<double> x);

    //for CKKS encoding
    std::vector<std::complex<double>> pi_function(std::vector<std::complex<double>>);
    std::vector<std::complex<double>> pi_inverse(std::vector<std::complex<double>>);
    std::vector<std::complex<double>> sigma_basis(std::vector<std::complex<double>>);

    //for coordinate-wise random rounding
    std::vector<std::complex<double>> coordinate_wise_random_rounding(std::vector<std::complex<double>>);   
    std::vector<std::complex<double>> sigma_R_discretization(std::vector<std::complex<double>>);
    std::vector<std::complex<double>> compute_basis_coordinates(std::vector<std::complex<double>>);
    std::vector<std::complex<double>> round_coordinates(std::vector<std::complex<double>>);

};

#endif // ENCODER_H