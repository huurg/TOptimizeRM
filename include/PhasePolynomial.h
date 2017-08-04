#ifndef PHASE_POLYNOMIAL_HEADER
#define PHASE_POLYNOMIAL_HEADER

#include <iostream>
using namespace std;

#include <string>

class PhasePolynomial {
    private:
        int* a = NULL;
        int n = 0;
        int N = 0;
    public:
        PhasePolynomial(int in_n);
        ~PhasePolynomial();

        void print() const;

        int get_n() const;
        int get_N() const;
        int getNoTerms() const;

        void operator+=(const PhasePolynomial& inPP);
        void operator*=(const int in_I);
        void operator-=(const PhasePolynomial& inPP);

        // Getters
        int operator[](const int in_I) const;
        int operator[](const bool* in_x) const;
        int operator[](const string in_str) const;

        // Setters
        int& operator[](const int in_I);
        int& operator[](const bool* in_x);
        int& operator[](const string in_str);

        // Methods
        void mod8(); // Maps elements from Z to Z_8
        void mod2(); // Maps elements from Z or Z_8 to Z_2

};

#endif // PHASE_POLYNOMIAL_HEADER
