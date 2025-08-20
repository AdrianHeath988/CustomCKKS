#include "encoder.h"
#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

Encoder::Encoder(int M){
    this->M = M;
    this->N = M/2;
    std::complex<double> i (0, 1.0);
    std::complex<double> pi (M_PI, 0.0);
    std::complex<double> two (2.0, 0.0);
    std::complex<double> complex_M (M + 0.0, 0.0);
    std::complex<double> power = (i * pi * two) / complex_M;
    this->xi = std::exp(power);
    std::cout << "----- Initial Root Of Unity:  "<< this->xi << "----" << std::endl;

}


void Encoder::printVector(std::vector<std::complex<double>> vec) {
    for (const auto& val : vec) {
        std::cout << "(" << val.real() << "," << val.imag() << ")\t";
    }
    std::cout << std::endl;
}

void Encoder::printMatrix(std::vector<std::vector<std::complex<double>>> matrix) { 
    for (const auto& row : matrix) {
        for (const auto& val : row) {
            std::cout << "(" << val.real() << "," << val.imag() << ")\t";
        }
        std::cout << std::endl;
    }
}

/**
 * Decodes by evaluating polynomial on roots of unity
 */
std::vector<std::complex<double>> Encoder::evaluate_sigma(std::vector<std::complex<double>> coeffs){
    std::vector<std::complex<double>> toReturn(N);
    for(int i=0;i<N;i++){
        int power = 2 * i + 1;
        std::complex<double> root = pow(this->xi, power);
        toReturn[i] = this->evaluatePolynomial(coeffs, root);
    }


    return toReturn;
}

/**
 * Encodes by multiplying the inverse vandermonde matrix by the input vector
 */
std::vector<std::complex<double>> Encoder::evaluate_sigma_inverse(std::vector<std::complex<double>> input_vector){
    //Multiply A^-1 by input vector
    std::vector<std::vector<std::complex<double>>> A = this->compute_vandermonde();
    std::cout << "----- Computed A -----" << std::endl;
    this->printMatrix(A);
    std::vector<std::vector<std::complex<double>>> A_inverse = this->matrix_inverse(A);
    std::cout << "----- Computed A inverse -----" << std::endl;
    this->printMatrix(A_inverse);
    std::vector<std::complex<double>> poly_coeffs = this->dotProduct(A_inverse,input_vector);
    std::cout << "----- Computed Dot Product -----" << std::endl;
    this->printVector(poly_coeffs);
    return poly_coeffs;
}



/**
 * Computes vandermonde matrix from root of unity
 */
std::vector<std::vector<std::complex<double>>> Encoder::compute_vandermonde(){
    //The matrix must be NxN size and each row will be: [1, p, p^2, p^3, ...]
    //Where p is the root of unity of that row
    //Root of unity for row 0 is xi, then for each row i, root of unity = xi^(2 * i + 1)
    std::vector<std::vector<std::complex<double>>> toReturn(N, std::vector<std::complex<double>>(N));
    for(int i=0;i<N;i++){   //row increment
        int power = 2 * i + 1;
        std::complex<double> root = pow(this->xi, power);

        for(int j=0;j<N;j++){   //column increment
            std::complex<double> toAdd = pow(root, j);
            toReturn[i][j] = toAdd;
        }
    }
    return toReturn;
}
/**
 * Computes Matrix Inverse
 */
std::vector<std::vector<std::complex<double>>> Encoder::matrix_inverse(std::vector<std::vector<std::complex<double>>> matrix){
    size_t n = matrix.size();
    

    // Create an augmented matrix [A|I]
    std::vector<std::vector<std::complex<double>>> aug(n, std::vector<std::complex<double>>(2 * n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            aug[i][j] = matrix[i][j];
        }
        aug[i][i + n] = 1.0;
    }

    // Apply simplified Gauss-Jordan elimination
    for (size_t i = 0; i < n; ++i) {
        std::complex<double> pivot = aug[i][i];
        if (std::abs(pivot) == 0) throw std::runtime_error("Matrix is singular.");
        for (size_t j = i; j < 2 * n; ++j) {
            aug[i][j] /= pivot;
        }

        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                std::complex<double> factor = aug[j][i];
                for (size_t k = i; k < 2 * n; ++k) {
                    aug[j][k] -= factor * aug[i][k];
                }
            }
        }
    }
    // Extract the inverse
    std::vector<std::vector<std::complex<double>>> inverse(n, std::vector<std::complex<double>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inverse[i][j] = aug[i][j + n];
        }
    }

    return inverse;

}
std::vector<std::complex<double>> Encoder::dotProduct(std::vector<std::vector<std::complex<double>>> matrix, std::vector<std::complex<double>> vec) {
    size_t matrix_rows = matrix.size();
    if (matrix_rows == 0) {
        return {};
    }
    size_t matrix_cols = matrix[0].size();
    size_t vec_size = vec.size();

    if (matrix_cols != vec_size) {
        throw std::invalid_argument("Matrix columns must match vector size for multiplication.");
    }

    std::vector<std::complex<double>> result(matrix_rows, {0.0, 0.0});
    for (size_t i = 0; i < matrix_rows; ++i) {
        for (size_t j = 0; j < matrix_cols; ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}

std::complex<double> Encoder::evaluatePolynomial(std::vector<std::complex<double>> coeffs, std::complex<double> x) {

    if (coeffs.empty()) {
        return {0.0, 0.0};
    }

    // Horner's method
    // We iterate backwards from the highest-degree coefficient
    std::complex<double> result (0.0, 0.0);
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
    }

    return result;
}