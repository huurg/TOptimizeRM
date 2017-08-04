#include "PhasePolynomial.h"

#include <iostream>
using namespace std;

#include <cmath>

#include "LukeBool.h"

PhasePolynomial::PhasePolynomial(int in_n) {
    n = in_n;
    N = (int)pow(2,n);
    if(n>0) {
        a = new int[N];
        for(int i = 0; i < N; i++) {
            a[i] = 0;
        }
    }
}

PhasePolynomial::~PhasePolynomial() {
    if(a) {
        delete [] a;
        a = NULL;
    }
}

void PhasePolynomial::print() const {
    /*for(int i = 1; i < N; i++) {
        cout << a[i] << " ";
    }
    cout << endl;*/
    cout << "f(x) =";
    bool first_term = 1;
    for(int i = 1; i < N; i++) {
        if(a[i]!=0) {
            cout << " ";
            if(a[i]>0) {
                if(!first_term) cout << "+ ";
                else first_term = 0;
            } else {
                cout << "- ";
                if(first_term) first_term = 0;
            }

            if(((int)fabs(a[i]))!=1) {
                cout << (int)fabs(a[i]);
            }

            cout << "(";
            bool x[n];
            LukeBool::IntToBoolVec(x,i,n);
            bool first_var = 1;
            for(int j = 0; j < n; j++) {
                if(x[j]) {
                    if(!first_var) cout << "^";
                    else first_var = 0;
                    cout << "x_" << (j+1);
                }
            }
            cout << ")";
        }
    }
    cout << endl;
}

int PhasePolynomial::getNoTerms() const {
    int out = 0;
    for(int i = 0; i < N; i++) {
        out += (a[i]>0);
    }
    return out;
}

int PhasePolynomial::operator[](const int in_I) const {
    int out = -1;

    if(a&&(in_I>=0)&&(in_I<N)) {
        out = a[in_I];
    }

    return out;
}

int PhasePolynomial::operator[](const bool* in_x) const {
    int out = -1;

    int in_I = LukeBool::BoolVecToInt(in_x, n);
    if(a&&(in_I>=0)&&(in_I<N)) {
        out = a[in_I];
    }

    return out;
}

int& PhasePolynomial::operator[](const int in_I) {
    if(a&&(in_I>=0)&&(in_I<N)) {
        return a[in_I];
    }
}

int& PhasePolynomial::operator[](const bool* in_x) {
    int in_I = LukeBool::BoolVecToInt(in_x, n);
    if(a&&(in_I>=0)&&(in_I<N)) {
        return a[in_I];
    }
}

int PhasePolynomial::get_n() const {
    return n;
}

int PhasePolynomial::get_N() const {
    return N;
}

void PhasePolynomial::operator+=(const PhasePolynomial& inPP) {
    if(n==inPP.n) {
        for(int i = 0; i < N; i++) {
            a[i] += inPP.a[i];
        }
    } else {
        cout << "Error! Phase polynomial dimensions must match." << endl;
    }
}

void PhasePolynomial::operator*=(const int in_I) {
    for(int i = 0; i < N; i++) {
        a[i] *= in_I;
    }
}

void PhasePolynomial::operator-=(const PhasePolynomial& inPP) {
    if(n==inPP.n) {
        for(int i = 0; i < N; i++) {
            a[i] -= inPP.a[i];
        }
    } else {
        cout << "Error! Phase polynomial dimensions must match." << endl;
    }
}

void PhasePolynomial::mod8() {
    for(int i = 0; i < N; i++) {
        while(a[i]<0) a[i] += 8;
        a[i] %= 8;
    }
}

void PhasePolynomial::mod2() {
    for(int i = 0; i < N; i++) {
        while(a[i]<0) a[i] += 2;
        a[i] %= 2;
    }
}
