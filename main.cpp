#include <iostream>
#include <iomanip>
#include "genetic.h"

int main()
{
    double long result;
    genetic prima("x^(1/x)", true, 3, 0, 9, 1, 0.1, 100, 400);
    result=prima.mejor_individuo();
    std::cout<<std::setprecision(15)<<result<<std::endl;

    return 0;
}
