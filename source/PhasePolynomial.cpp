#include "PhasePolynomial.h"

#include <iostream>
using namespace std;

#include <cmath>
#include <unordered_set>

#include "LukeBool.h"
#include "LukeConsoleOut.h"
using namespace LukeConsoleOut;


PhasePolynomial::PhasePolynomial(int in_n) {
    n = in_n;
    N = (int)pow(2,n);
}


PhasePolynomial::~PhasePolynomial() {
    ;
}

void PhasePolynomial::print() const {
    //warning("TODO.", "print", "PhasePolynomial");
    LOut() << "f(x) =";

    for(int t = 0; t < T(); t++) {
        if((m[t]<0)||(t>0)) {
            // Do print the sign
            if(m[t]<0) LOut() << " -";
            else LOut() << " +";
        }
        LOut() << " ";
        if(((int)fabs(m.at(t))!=1)) LOut() << (int)fabs(m[t]);
        LOut() << "(";
        bool x[n];
        LukeBool::IntToBoolVec(x,a[t],n);
        bool first_var = 1;
        for(int i = 0; i < n; i++) {
            if(x[i]) {
                if(first_var) {
                    first_var = 0;
                } else {
                    LOut() << "^";
                }
                LOut() << "x_" << (i+1);
            }
        }
        LOut() << ")";
    }
    LOut() << endl;
}

int PhasePolynomial::get_n() const {
    return n;
}

int PhasePolynomial::get_N() const {
    return N;
}

int PhasePolynomial::T() const {
    int out = 0;

    if(m.size()!=a.size()) {
        warning("m and a should be the same length.", "getNoTerms", "PhasePolynomial");
    }

    out = m.size();

    return out;
}

int PhasePolynomial::get_m_at(const int in_t) const {
    return m.at(in_t);
}

int PhasePolynomial::get_a_at(const int in_t) const {
    return a.at(in_t);
}

void PhasePolynomial::operator+=(const PhasePolynomial& inPP) {
    if(inPP.get_n()==n) {
        for(int t = 0; t < inPP.T(); t++) {
            operator[](inPP.a.at(t)) += inPP.m.at(t);
        }
    } else {
        error("Operands should have matching n.", "operator+=", "PhasePolynomial");
    }
}

void PhasePolynomial::operator*=(const int in_m) {
    //warning("TODO.", "operator*=", "PhasePolynomial");
    for(int t = 0; t < T(); t++) {
        m.at(t) *= in_m;
    }
}

void PhasePolynomial::operator-=(const PhasePolynomial& inPP) {
    //warning("TODO.", "operator-=", "PhasePolynomial");
    if(inPP.get_n()==n) {
        for(int t = 0; t < inPP.T(); t++) {
            operator[](inPP.a.at(t)) -= inPP.m.at(t);
        }
    } else {
        error("Operands should have matching n.", "operator-=", "PhasePolynomial");
    }
}

void PhasePolynomial::operator=(const PhasePolynomial& inPP) {
    //warning("TODO.", "operator=", "PhasePolynomial");
    if(inPP.get_n()==n) {
        m.clear();
        a.clear();
        for(int t = 0; t < inPP.T(); t++) {
            operator[](inPP.a.at(t)) = inPP.m.at(t);
        }
    } else {
        error("Operands should have matching n.", "operator=", "PhasePolynomial");
    }
}

void PhasePolynomial::operator%=(const int in_m) {
    //warning("TODO.", "operator%=", "PhasePolynomial");
    if(in_m>0) {
        for(int t = 0; t < T(); t++) {
            while(m.at(t)<0) m.at(t) += in_m;
            m.at(t) %= in_m;
        }
    } else {
        error("Modulus should be a positive, non-zero integer.", "operator%=", "PhasePolynomial");
    }
}

int PhasePolynomial::operator[](const int in_I) const {
    int out = 0;

    //warning("TODO.", "operator[]", "PhasePolynomial");
    bool found = 0;
    for(int t = 0; (!found)&&(t<T()); t++) {
        if(a.at(t)==in_I) {
            found = 1;
            out = m.at(t);
        }
    }

    return out;
}

int PhasePolynomial::operator[](const bool* in_x) const {
    int out = 0;

    //warning("TODO.", "operator[]", "PhasePolynomial");
    int this_I = LukeBool::BoolVecToInt(in_x,n);
    out = operator[](this_I);

    return out;
}

int PhasePolynomial::operator[](const string in_str) const {
    int out = 0;

    warning("TODO.", "operator[]", "PhasePolynomial");

    return out;
}

int& PhasePolynomial::operator[](const int in_I) {
    //warning("TODO.", "operator[]", "PhasePolynomial");

    for(int t = 0; t<T(); t++) {
        if(a.at(t)==in_I) {
            return m.at(t);
        }
    }
    m.push_back(0);
    a.push_back(in_I);
    return m.back();
}

int& PhasePolynomial::operator[](const bool* in_x) {
    //warning("TODO.", "operator[]", "PhasePolynomial");
    int this_I = LukeBool::BoolVecToInt(in_x,n);
    return operator[](this_I);
}

int& PhasePolynomial::operator[](const string in_str) {
    int out = 0;

    warning("TODO.", "operator[]", "PhasePolynomial");

    return out;
}

/*
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

PhasePolynomial::PhasePolynomial(const PhasePolynomial& in) {//*
    PhasePolynomial(in.n);
    for(int i = 0; i < in.N; i++) {
        a[i] = in.a[i];
    }
}

PhasePolynomial::~PhasePolynomial() {
    if(a) {
        delete [] a;
        a = NULL;
    }
}

void PhasePolynomial::print() const {
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

void PhasePolynomial::operator=(const PhasePolynomial& in) {
    for(int i = 0; i < N; i++) {
        a[i] = in.a[i];
    }
}
*/
