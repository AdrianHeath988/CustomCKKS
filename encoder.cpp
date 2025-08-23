#include "encoder.h"
#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <chrono>

Encoder::Encoder(int M, int scale){
    this->M = M;
    this->scale = scale;
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
std::vector<std::complex<double>> Encoder::decode(std::vector<std::complex<double>> coeffs){
    std::vector<std::complex<double>> de_scaled = this->scale_vector(coeffs, (1.0/this->scale));
    std::vector<std::complex<double>> decoded = this->evaluate_sigma(de_scaled);
    std::vector<std::complex<double>> truncated = this->pi_function(decoded);
    return truncated;
}


std::vector<std::complex<double>> Encoder::encode(std::vector<std::complex<double>> input_vector){
    //truncate:
    std::vector<std::complex<double>> projected_to_H = this->pi_inverse(input_vector);
    std::vector<std::complex<double>> scaled = this->scale_vector(projected_to_H, (this->scale/1.0));
    std::vector<std::complex<double>> projected_sigma_basis = this->sigma_R_discretization(scaled);
    std::vector<std::complex<double>> encoded = this->evaluate_sigma_inverse(projected_sigma_basis);
    return encoded;
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
std::complex<double> Encoder::vectorDot(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2){
    std::complex<double> result = {0.0, 0.0};
    for (size_t i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * std::conj(vec2[i]);
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
/**
 * Take vector in Canonical Subspace (H) and output vector in C^(N/2)
 * Key is that Canonical Subspace is defined as 2nd half elements being complex conjugates of first half
 * therefore, we can just split in half 
 */
std::vector<std::complex<double>> Encoder::pi_function(std::vector<std::complex<double>> input){
    std::vector<std::complex<double>> toReturn (input.begin(), input.begin() + N/2);
    std::cout << "----- Computed Pi Function -----" << std::endl;
    this->printVector(toReturn);
    return toReturn;
}
/**
 * Doubnles vector side, adding complex conjugate in the 2nd half
 */
std::vector<std::complex<double>> Encoder::pi_inverse(std::vector<std::complex<double>> input){
    std::vector<std::complex<double>> toReturn = input;
    toReturn.resize(N); 
    for(int i = 0; i < N/2; i++){
        toReturn[N - 1 - i] = std::conj(input[i]);
    }
    
    std::cout << "----- Computed Inverse Pi Function -----" << std::endl;
    this->printVector(toReturn);
    return toReturn;
}

/**
 * Multiplies by scaling factor
 */
std::vector<std::complex<double>> Encoder::scale_vector(std::vector<std::complex<double>> input, double factor){
    std::vector<std::complex<double>> toReturn;
    for(int i=0;i<input.size();i++){
        toReturn.push_back(input[i] * static_cast<double>(factor));
    }
    return toReturn;
}

/**
 * Rounds coordinates randomly to (x) or (x) - 1, with probablility [1-x, x] -> [x, x-1]
 * Ex. if value is 0.2, theres a 80% chance it resolves to .2, and a 20% chance it resolves to -.8
 */
std::vector<std::complex<double>> Encoder::coordinate_wise_random_rounding(std::vector<std::complex<double>> input){
    std::vector<std::complex<double>> rounded = this->round_coordinates(input);
    static std::mt19937 generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<std::complex<double>> toReturn;
    toReturn.reserve(rounded.size());
    for (size_t i = 0; i < rounded.size(); ++i) {
        std::complex<double> c = rounded[i];
        std::complex<double> f_val;

        if (distribution(generator) < std::real(c)) {
            f_val = c - 1.0; // Chosen with probability 'c'.
        } else {
            f_val = c;       // Chosen with probability '1-c'.
        }
        std::complex<double> result = input[i] - f_val;
        result = {static_cast<int> (std::real(result)), std::imag(result)};
        toReturn.push_back((result)); 
    }
    return toReturn;

}
/**
 * Projects a vector on the lattice using coordinate-wise random rounding
 */
std::vector<std::complex<double>> Encoder::sigma_R_discretization(std::vector<std::complex<double>> input){
    std::vector<std::complex<double>> coords = this->compute_basis_coordinates(input);
    std::vector<std::complex<double>> rounded_coords = this->coordinate_wise_random_rounding(coords);
    std::vector<std::vector<std::complex<double>>> vandermonde_matrix = this->compute_vandermonde();
    // Dot prodduct between rounded_coords and vandermonde matrix
    std::vector<std::complex<double>> result = this->dotProduct(vandermonde_matrix, rounded_coords);
    std::cout << "----- Computed sigma_R_discretization -----" << std::endl;
    this->printVector(result);
    return result;
}
/**
 * Computes the coordinates of a vector in the new basis
 */
std::vector<std::complex<double>> Encoder::compute_basis_coordinates(std::vector<std::complex<double>> input){
    std::vector<std::vector<std::complex<double>>> vandermonde_matrix = this->compute_vandermonde();
    std::vector<std::complex<double>> toReturn;
    for(int i=0;i<vandermonde_matrix.size();i++){   
        //compute vdot(intput, ith column)/vdot(ith column, ith column)
        std::vector<std::complex<double>> ith_column;
        for(int j=0;j<vandermonde_matrix.size();j++){
            ith_column.push_back(vandermonde_matrix[j][i]);
        }
        std::cout << "----- Computed ith column -----" << std::endl;
        this->printVector(ith_column);
        std::complex<double> quotient = this->vectorDot(input, ith_column);
        std::complex<double> denominator = this->vectorDot(ith_column, ith_column);
        std::complex<double> toAdd = quotient/denominator;
        toReturn.push_back(toAdd);
    }
    std::cout << "----- Computed basis coordinates -----" << std::endl;
    this->printVector(toReturn);
    return toReturn;
}
/**
 * Returns coordinate - floor(coordinate) for each
 */
std::vector<std::complex<double>> Encoder::round_coordinates(std::vector<std::complex<double>> input){
    std::vector<std::complex<double>> toReturn;
    for(int i=0;i<input.size();i++){
        toReturn.push_back(input[i] - std::complex<double> (floor(std::real(input[i])), floor(std::imag(input[i]))));

    }
    std::cout << "----- Computed Rounded Coords -----" << std::endl;
    this->printVector(toReturn);
    return toReturn;
}

