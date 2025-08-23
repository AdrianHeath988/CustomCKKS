#include <stdint.h>
#include "encoder.h"
#include <iostream>
int main(){
    Encoder myEncoder(8, 128);
    std::vector<std::complex<double>> input;
    input.push_back(std::complex<double>(-5, 0));
    input.push_back(std::complex<double>(1, 0));
    // input.push_back(std::complex<double>(1, 0));
    // input.push_back(std::complex<double>(0, 0));
    std::cout << "----- Input -----" << std::endl;
    myEncoder.printVector(input);
    std::vector<std::complex<double>> encoded = myEncoder.encode(input);
    std::cout << "----- Encoded -----" << std::endl;
    myEncoder.printVector(encoded);
    std::vector<std::complex<double>> decoded = myEncoder.decode(encoded);
    std::cout << "----- Output -----" << std::endl;
    myEncoder.printVector(decoded);

}