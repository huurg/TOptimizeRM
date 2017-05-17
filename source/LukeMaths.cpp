#include "LukeMaths.h"

#include <iostream>
using namespace std;

#include <cmath>

unsigned long long int LukeMaths::fact(int n) {
    unsigned long long int out = 1;

    for(int i = 2; i <= n; i++) out *= i;

    return out;
}

unsigned long long int LukeMaths::nCr(int n, int r) {
    unsigned long long int out = fact(n)/(fact(r)*fact(n-r));
    return out;
}
