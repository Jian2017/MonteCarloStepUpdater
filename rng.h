#ifndef RNG_H
#define RNG_H

// include this header file to use dis(gen) as pseduo random numbers
// dis(gen)= rectangle distribution [0,1)
// distribution(gen)= either 0 or 1


#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0,1);
std::uniform_int_distribution<int> dis01(0,1);



#endif
