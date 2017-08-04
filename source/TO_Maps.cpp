#include "TO_Maps.h"

#include <iostream>
using namespace std;

#include "PhasePolynomial.h"
#include "SQC_Circuit.h"

#include "BoolMat.h"
#include "LukeBool.h"

PhasePolynomial TO_Maps::SQC_Circuit_to_PhasePolynomial(const SQC_Circuit& in) {
    int n = in.n;
    PhasePolynomial out(n);

    BoolMat E_mat(n,n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            E_mat.E(i,j) = (i==j);
        }
    }

    for(int t = 0; t < in.m; t++) {
        switch(in.operator_list[t][0]) {
            case SQC_OPERATOR_CNOT:
                {
                    int c = in.operator_list[t][2]-1;
                    int d = in.operator_list[t][1]-1;
                    BoolMat this_CNOT(n,n);
                    for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++) {
                            this_CNOT.E(i,j) = ((i==j)+(i==d)*(j==c))%2;
                        }
                    }
                    E_mat = (this_CNOT*E_mat);
                }
                break;
            case SQC_OPERATOR_Z:
            case SQC_OPERATOR_S:
            case SQC_OPERATOR_T:
                {
                    int q = in.operator_list[t][1]-1;
                    bool x[n];
                    for(int j = 0; j < n; j++) x[j] = E_mat.E(q,j);
                    int m = 0;
                    switch(in.operator_list[t][0]) {
                        case SQC_OPERATOR_Z:
                            m = 4;
                            break;
                        case SQC_OPERATOR_S:
                            m = 2;
                            break;
                        case SQC_OPERATOR_T:
                            m = 1;
                            break;
                    }
                    out[x] += m;
                }
                break;
            default:
                cout << "Warning in TO_Maps::PhasePolynomial_to_SQC_Circuit! Expected gates are {CNOT, Z, S, T} only." << endl;
                break;
        }
    }

    return out;
}

GateStringSparse TO_Maps::PhasePolynomial_to_GateStringSparse(const PhasePolynomial& in) {
    int n = in.get_n();
    GateStringSparse out(n);

    for(int i = 0; i < in.get_N(); i++) {
        if(in[i]%2) {
            out.set(i);
        }
    }

    return out;
}
