#include <stdint.h>
#include "encoder.h"
#include <iostream>
int main(){
    Encoder myEncoder(8);
    std::vector<std::complex<double>> input = {1, 2, 3, 4};
    std::cout << "----- Input -----" << std::endl;
    myEncoder.printVector(input);
    std::vector<std::complex<double>> encoded = myEncoder.evaluate_sigma_inverse(input);
    std::cout << "----- Encoded -----" << std::endl;
    myEncoder.printVector(encoded);
    std::vector<std::complex<double>> decoded = myEncoder.evaluate_sigma(encoded);
    std::cout << "----- Output -----" << std::endl;
    myEncoder.printVector(decoded);

}