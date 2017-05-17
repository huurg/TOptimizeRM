#include "SQC_Circuit.h"

#include <iostream>
using namespace std;
#include <ostream>
#include <string>
#include <cstring>
#include <fstream>
#include <utility>
#include <cmath>

#include "Signature.h"
#include "BMSparse.h"
#include "LukeMat_GF2.h"
#include "Interface_SigBMS.h"
#include "GateSynthesisMatrix.h"
#include "GateStringSparse.h"
#include "GateSigInterface.h"
#include "Interface_BMSGSS.h"

SQC_Circuit::SQC_Circuit() {
    ;
}

void SQC_Circuit::Construct() {
    if(!operator_list) operator_list = new int*[max_m];
    for(int i = 0; i < max_m; i++){
        operator_list[i] = new int[n+1];
        for(int j = 0; j < (n+1); j++) {
            operator_list[i][j] = 0;
        }
    }
}

void SQC_Circuit::Destruct() {
    if(operator_list) {
        for(int i = 0; i < max_m; i++) {
            if(operator_list[i]) {
                    delete [] operator_list[i];
                operator_list[i] = NULL;
            }
        }
        delete [] operator_list;
        operator_list = NULL;
    }
}

void SQC_Circuit::Copy(const SQC_Circuit& in_C) {
    Destruct();
    n = in_C.n;
    p = in_C.p;
    max_m = in_C.max_m;
    m = 0;
    ancilla_mode = in_C.ancilla_mode;
    Construct();
    for(int i = 0; i < in_C.m; i++) {
        AddOperator(in_C.operator_list[i]);
    }
}

void SQC_Circuit::Print(ostream* in_OS) const {
    (*in_OS) << "n " << n << endl;
    (*in_OS) << "m " << max_m << endl;
    (*in_OS) << "p " << p << endl;
    for(int i = 0; i < m; i++) {
        string this_op_str;
        switch(operator_list[i][0]) {
        case SQC_OPERATOR_IDENTITY:
            this_op_str = SQC_OPSTRING_IDENTITY;
            break;
        case SQC_OPERATOR_HADAMARD:
            this_op_str = SQC_OPSTRING_HADAMARD;
            break;
        case SQC_OPERATOR_CNOT:
            this_op_str = SQC_OPSTRING_CNOT;
            break;
        case SQC_OPERATOR_T:
            this_op_str = SQC_OPSTRING_T;
            break;
        case SQC_OPERATOR_CS:
            this_op_str = SQC_OPSTRING_CS;
            break;
        case SQC_OPERATOR_CCZ:
            this_op_str = SQC_OPSTRING_CCZ;
            break;
        case SQC_OPERATOR_S:
            this_op_str = SQC_OPSTRING_S;
            break;
        case SQC_OPERATOR_Z:
            this_op_str = SQC_OPSTRING_Z;
            break;
        case SQC_OPERATOR_TOFFOLI:
            this_op_str = SQC_OPSTRING_TOFFOLI;
            break;
        case SQC_OPERATOR_TOFFOLI_4:
            this_op_str = SQC_OPSTRING_TOFFOLI_4;
            break;
        case SQC_OPERATOR_TOFFOLI_N:
            this_op_str = SQC_OPSTRING_TOFFOLI_N;
            break;
        }
        (*in_OS) << this_op_str;
        int j = 0;
        int this_qubit=0;
        while((this_qubit = operator_list[i][j+1])&&(j<n)) {
            (*in_OS) << " " << this_qubit;
            j++;
        }
        (*in_OS) << endl;
    }
    //cout << "Finish print." << endl;
}

void SQC_Circuit::AddOperator(const char* in_op_str) {
    // Resize if necessary
    if(m>=max_m) {
        Resize(max_m+1);
        cout << "Circuit expanded by 1" << endl;
    }
    if(m<max_m) {
        int linewidth = 2*(2*n+2);
        char* op_str = NULL;
        if(linewidth>0) {
            op_str = new char[linewidth];
            strcpy(op_str,in_op_str);
            char* this_tok = NULL;
            this_tok = strtok(op_str," ");
            if(this_tok) {
                if(!strcmp(this_tok,SQC_OPSTRING_IDENTITY)) {
                    operator_list[m][0] = SQC_OPERATOR_IDENTITY;
                } else if(!strcmp(this_tok,SQC_OPSTRING_HADAMARD)) {
                    operator_list[m][0] = SQC_OPERATOR_HADAMARD;
                } else if(!strcmp(this_tok,SQC_OPSTRING_CNOT)) {
                    operator_list[m][0] = SQC_OPERATOR_CNOT;
                } else if(!strcmp(this_tok,SQC_OPSTRING_T)) {
                    operator_list[m][0] = SQC_OPERATOR_T;
                } else if(!strcmp(this_tok,SQC_OPSTRING_CS)) {
                    operator_list[m][0] = SQC_OPERATOR_CS;
                } else if(!strcmp(this_tok,SQC_OPSTRING_CCZ)) {
                    operator_list[m][0] = SQC_OPERATOR_CCZ;
                } else if(!strcmp(this_tok,SQC_OPSTRING_S)) {
                    operator_list[m][0] = SQC_OPERATOR_S;
                } else if(!strcmp(this_tok,SQC_OPSTRING_Z)) {
                    operator_list[m][0] = SQC_OPERATOR_Z;
                } else if(!strcmp(this_tok, SQC_OPSTRING_TOFFOLI)) {
                    operator_list[m][0] = SQC_OPERATOR_TOFFOLI;
                } else if(!strcmp(this_tok, SQC_OPSTRING_TOFFOLI_4)) {
                    operator_list[m][0] = SQC_OPERATOR_TOFFOLI_4;
                } else if(!strcmp(this_tok, SQC_OPSTRING_TOFFOLI_N)) {
                    operator_list[m][0] = SQC_OPERATOR_TOFFOLI_N;
                }
                int i = 0;
                bool exit = false;
                while((!exit)&&(i<n)&&(this_tok = strtok(NULL," "))) {
                    //cout << this_tok << endl;
                    int this_qubit = atoi(this_tok);

                    if((this_qubit>0)&&(this_qubit<=n)) {
                        operator_list[m][i+1] = this_qubit;
                        i++;
                    } else {
                        exit = true;
                    }
                }

                bool success = true;
                switch(operator_list[m][0]) {
                    case SQC_OPERATOR_CNOT:
                    case SQC_OPERATOR_CS:
                        if(i!=2) success = false;
                        break;
                    case SQC_OPERATOR_CCZ:
                    case SQC_OPERATOR_TOFFOLI:
                        if(i!=3) success = false;
                        break;
                    case SQC_OPERATOR_TOFFOLI_4:
                        if(i!=5) success = false;
                        break;
                    case SQC_OPERATOR_TOFFOLI_N:
                        if((p>0)&&((i/2)==((double)i/2.0))) success = false;
                        break;
                }
                if(!success) cout << "WARNING! Incorrect number of qubits." << endl;
                m++;
                //cout << "op_str = " << op_str << endl;
                //cout << "n_args = " << i << endl;
            }
        }

        //cout << "QWER" << endl;
        if(op_str) delete [] op_str;
        //cout << "ASDF" << endl;
    }

}

void SQC_Circuit::Load(const char* in_filename) {
    //cout << "QWER" << endl;
    n = 0;
    max_m = 0;
    m = 0;
    ifstream my_file(in_filename);
    if(my_file.good()) {
        int linewidth = 1000;
        char* this_line = new char[linewidth];

        //cout << "asdf" << endl;
        // Get n (1st line)
        my_file.getline(this_line,linewidth);
        //cout << "Out 1: " << this_line << endl;
        char* this_tok = NULL;
        this_tok = strtok(this_line," ");
        this_tok = strtok(NULL," ");
        //cout << "Out 2: " << this_tok << endl;
        n = atoi(this_tok);
        if(n==0) cout << "Warning! n should take non-zero value. n = " << n << endl;

        //cout << "zxcv" << endl;
        // Get max_m (2nd line)
        my_file.getline(this_line,linewidth);
        //cout << "Out 2: " << this_line << endl;
        this_tok = strtok(this_line," ");
        this_tok = strtok(NULL," ");
        max_m = atoi(this_tok);
        if(max_m<=0) max_m = SQC_DEFAULT_MAX_M;

        // Get p (3rd line)
        my_file.getline(this_line,linewidth);
        //cout << "Out 3: " << this_line << endl;
        this_tok = strtok(this_line," ");
        this_tok = strtok(NULL," ");
        p = atoi(this_tok);
        if(p==-1) ancilla_mode = SQC_ANCILLA_MODE_PER_GATE;
        else if(p==-2) ancilla_mode = SQC_ANCILLA_MODE_PER_CIRCUIT;

        //cout << "tuyi" << endl;
        if(n&&max_m) {
            Construct();
            //cout << "ghjk" << endl;
            //my_file.getline(this_line,linewidth);
            while(/*(m<max_m)&&*/(!my_file.eof())) {
                my_file.getline(this_line,linewidth);
                //cout << "Line " << m << ": " << this_line << endl;
                if(strlen(this_line)) AddOperator(this_line);
            }
        }

        delete [] this_line;
        my_file.close();
    }
}

void SQC_Circuit::Save(const char* in_filename) const {
    ofstream my_file(in_filename);

    Print(&my_file);

    my_file.close();
}

void SQC_Circuit::Clear() {
    m = 0;
}

bool SQC_Circuit::GetPartition(SQC_Circuit* out_Hadamards, SQC_Circuit* out_CNOT_T) {
    cout << "QWER" << endl;
    bool out = false;
    if(operator_list&&out_Hadamards->operator_list&&out_CNOT_T->operator_list) {
        out = (bool)m;

        if(out) {
            // Pre-process Toffoli gates by converting them to CCZ conjugated by Hadamards on the target
            cout << "TYUI" << endl;
            ConvertFromToffoli();

            cout << "Out 1 n = " << n << endl;
            int init_m = m;
            bool* qubit_used = new bool[n];
            for(int i = 0; i < n; i++) qubit_used[i] = 0;
            bool** hadamard_met = new bool*[init_m];
            for(int i = 0; i < init_m; i++) {
                hadamard_met[i] = new bool[n];
                for(int j = 0; j < n; j++) {
                    hadamard_met[i][j] = 0;
                }
            }
            cout << "Out 2" << endl;
            for(int t = 0; t < m; t++) {
                //cout << "Out t " << t << endl;
                // If the operator is a Hadamard
                if(operator_list[t][0]==SQC_OPERATOR_HADAMARD) {
                    //cout << "Out Hadamard" << endl;
                    // If the qubit upon which it acts hasn't been used yet
                    if(!qubit_used[operator_list[t][1]-1]) {
                        //cout << "Out not used yet" << endl;
                        // Add it to out_Hadamards and remove it from the circuit
                        out_Hadamards->AddOperator(operator_list[t]);
                        //cout << "Out added" << endl;
                        DeleteOperator(t);
                        //cout << "Out deleted" << endl;
                        t--;
                    }
                } else {
                    // If non of the qubits upon which the operator acts have been acted upon by a Hadamard (except 'free' Hadamards)
                    //bool can_move_left = 0;
                    //cout << "Out Not Hadamard" << endl;
                    switch(operator_list[t][0]) {
                    case SQC_OPERATOR_IDENTITY:
                        break;
                    case SQC_OPERATOR_T:
                    case SQC_OPERATOR_S:
                    case SQC_OPERATOR_Z:
                        qubit_used[operator_list[t][1]-1]=1;
                        break;
                    case SQC_OPERATOR_CNOT:
                    case SQC_OPERATOR_CS:
                        qubit_used[operator_list[t][1]-1]=1;
                        qubit_used[operator_list[t][2]-1]=1;
                        break;
                    case SQC_OPERATOR_CCZ:
                        qubit_used[operator_list[t][1]-1]=1;
                        qubit_used[operator_list[t][2]-1]=1;
                        qubit_used[operator_list[t][3]-1]=1;
                        break;
                    }
                }
            }
            cout << "Out 3" << endl;
            for(int t = 0; t < m; t++) {
                // If the operator is a Hadamard
                if(operator_list[t][0]==SQC_OPERATOR_HADAMARD) {
                    for(int i = t; i < m; i++) {
                        hadamard_met[i][operator_list[t][1]-1]=1;
                    }
                } else {
                    // If non of the qubits upon which the operator acts have been acted upon by a Hadamard (except 'free' Hadamards)
                    bool can_move_left = 0;
                    switch(operator_list[t][0]) {
                    case SQC_OPERATOR_IDENTITY:
                    case SQC_OPERATOR_T:
                    case SQC_OPERATOR_S:
                    case SQC_OPERATOR_Z:
                        can_move_left = !hadamard_met[t][operator_list[t][1]-1];
                        break;
                    case SQC_OPERATOR_CNOT:
                    case SQC_OPERATOR_CS:
                        can_move_left = !hadamard_met[t][operator_list[t][1]-1]&&!hadamard_met[t][operator_list[t][2]-1];

                        break;
                    case SQC_OPERATOR_CCZ:
                        can_move_left = !hadamard_met[t][operator_list[t][1]-1]
                                     && !hadamard_met[t][operator_list[t][2]-1]
                                     && !hadamard_met[t][operator_list[t][3]-1];
                        break;
                    }
                    if(can_move_left) {
                        out_CNOT_T->AddOperator(operator_list[t]);
                        DeleteOperator(t);
                        t--;
                    }
                }
            }
            cout << "Out 4 init_m = " << init_m << endl;
            for(int i = 0; i < init_m; i++) {
                delete [] hadamard_met[i];
                hadamard_met[i] = NULL;
            }
            cout << "Out 4.1" << endl;
            delete [] hadamard_met;
            hadamard_met = NULL;
            cout << "Out 4.2 qubit_used = " << qubit_used << endl;
            delete [] qubit_used;
            qubit_used = NULL;
            cout << "Out 4.3" << endl;
            out = (bool)m;
            cout << "Out 5" << endl;
        }
    }
    return out;
}

void SQC_Circuit::DeleteOperator(int t) {
    if((t<m)&&(t>=0)) {
        for(int i = 0; i < (n+1); i++) operator_list[t][i] = 0;
        for(int i = t; i < (m-1); i++) {
            swap(operator_list[i],operator_list[i+1]);
        }
        m--;
    }
}

void SQC_Circuit::AddOperator(const SQC_Operator in_op) {
    // Resize if necessary
    if(m>=max_m) {
        Resize(max_m+1);
        cout << "Circuit expanded by 1" << endl;
    }
    if(m<max_m) {
        for(int i = 0; i < (n+1); i++) {
            operator_list[m][i] = in_op[i];
        }
        m++;
    }
}

void SQC_Circuit::DecompositionVW(SQC_Circuit* out_V, SQC_Circuit* out_W) const {
    if(operator_list&&out_V->operator_list&&out_W->operator_list) {
        SQC_Circuit circuit_Cp;
        circuit_Cp.n = n;
        circuit_Cp.max_m = max_m;
        circuit_Cp.Construct();

        // Create the CNOT-only version of C
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0]==SQC_OPERATOR_CNOT) {
                circuit_Cp.AddOperator(operator_list[i]);
            }
        }

        // Loop until no cancellations are found for all operators
        bool cancellations = true;
        while(cancellations) {
            cancellations = false;
            // For each CNOT operator until a cancellation is found
            for(int i = 0; (!cancellations)&&(i < (circuit_Cp.m-1)); i++) {
                SQC_Operator this_CNOT = circuit_Cp.operator_list[i];
                bool terminate_search = false;
                for(int ip = i+1; (!terminate_search)&&(!cancellations)&&(ip<circuit_Cp.m); ip++) {
                    SQC_Operator search_CNOT = circuit_Cp.operator_list[ip];
                    // If t=tp and c=cp, delete both operators, set cancellations to true
                    if((this_CNOT[1]==search_CNOT[1])&&(this_CNOT[2]==search_CNOT[2])) {
                        cancellations = true;
                        circuit_Cp.DeleteOperator(ip);
                        circuit_Cp.DeleteOperator(i);
                    }
                    // If t=cp or c=tp, terminate loop as can't cancel for this CNOT
                    else if((this_CNOT[1]==search_CNOT[2])||(this_CNOT[2]==search_CNOT[1])) {
                        terminate_search = true;
                    }
                }
            }
        }

        // Copy Cp to V
        for(int i = 0; i < circuit_Cp.m; i++) {
            out_V->AddOperator(circuit_Cp.operator_list[i]);
        }


        // Copy C to W
        for(int i = 0; i < m; i++) {
            out_W->AddOperator(operator_list[i]);
        }

        // Copy Cp^dagger to W
        for(int i = (circuit_Cp.m-1); i >= 0; i--) {
            out_W->AddOperator(circuit_Cp.operator_list[i]);
        }

        circuit_Cp.Destruct();
    } else {
        cout << "Error in SQC_Circuit::DecompositionVW. At least one SQC_Circuit has not been intialized." << endl;
    }
    // TODO remove redundant CNOTs from output W circuit.
}

BMSparse SQC_Circuit::toGateSynthesisMatrix() const {
    Signature out(n);

    BMSparse A(n,0);
    BMSparse P = BMSparse::eye(n,n);
    //cout << "P" << endl;
    //P.printFull();

    for(int t = 0; t < m; t++) {
        Signature temp_sig(n);
        switch(operator_list[t][0]) {
            case SQC_OPERATOR_CNOT:
                P = P.addRows(operator_list[t][1]-1, operator_list[t][2]-1);
                //cout << "P" << endl;
                //P.printFull();
                break;
            case SQC_OPERATOR_T:
                temp_sig.set(operator_list[t][1]);
                break;
            case SQC_OPERATOR_CS:
                temp_sig.set(operator_list[t][1], operator_list[t][2]);
                break;
            case SQC_OPERATOR_CCZ:
                temp_sig.set(operator_list[t][1], operator_list[t][2], operator_list[t][3]);
                break;
        }
        if(!temp_sig.isEmpty()) {
            GateStringSparse new_cols_GSS = GateSigInterface::SigToGSS(temp_sig);
            BMSparse new_cols = Interface_BMSGSS::GSSToBMS(new_cols_GSS);
            //cout << "new cols" << endl;
            //new_cols.printFull();
            new_cols = P.T()*new_cols;
            //cout << "new cols'" << endl;
            //new_cols.printFull();
            A = A && new_cols;
        }
    }

    //cout << "A" << endl;
    //A.printFull();

    int A_n = A.get_n();
    int A_m = A.get_m();
    bool** A_ba = LukeMat_GF2::construct(A_n, A_m);
    A.toBool(A_ba);
    //LukeMat_GF2::print(A_ba,A_n,A_m,"A_ba");
    int A_mp;
    GateSynthesisMatrix::cleanup(A_ba,A_n,A_m,A_mp);

    A.fromBool(A_ba,A_n,A_mp);

    LukeMat_GF2::destruct(A_ba,A_n,A_m);

    out = Interface_SigBMS::BMSToSig(A);

    //cout << "A'" << endl;
    //A.printFull();

    return A;
}

bool SQC_Circuit::NextSignature(Signature& outsig) {
    //cout << "QWER" << endl;
    outsig = Signature(n);
    SQC_Circuit hadamards, CNOT_Ts;
    hadamards.n = CNOT_Ts.n = n;
    hadamards.max_m = CNOT_Ts.max_m = max_m;
    hadamards.Construct();
    CNOT_Ts.Construct();
    cout << "xcvxcv" << endl;
    bool out = GetPartition(&hadamards, &CNOT_Ts);
    //cout << "Hadamards" << endl;
    //hadamards.Print();
    //cout << endl;
    //cout << "{CNOT,T}" << endl;
    //CNOT_Ts.Print();
    //cout << endl;
    cout << "asdf" << endl;
    if((bool)CNOT_Ts.m) {
        //cout << "TYUI" << endl;
        SQC_Circuit this_V, this_W;
        this_V.n = this_W.n = n;
        this_V.max_m = this_W.max_m = CNOT_Ts.max_m;
        this_V.Construct();
        this_W.Construct();
        CNOT_Ts.DecompositionVW(&this_V, &this_W);
        //cout << "V" << endl;
        //this_V.Print();
        //cout << endl;
        //cout << "W" << endl;
        //this_W.Print();
        //cout << endl;
        BMSparse this_W_BMS = this_W.toGateSynthesisMatrix();
        //cout << "W matrix" << endl;
        //this_W_BMS.printFull();
        outsig = Interface_SigBMS::BMSToSig(this_W_BMS);
        //cout << "Signature of W" << endl;
        //outsig.print();
        this_V.Destruct();
        this_W.Destruct();
    }
    //cout << "zxcv" << endl;
    hadamards.Destruct();
    CNOT_Ts.Destruct();
    return out;
}

void SQC_Circuit::ReplaceOperator(SQC_Circuit* in_new_ops, int t, int n_rep) {
    // Create temporary copy of original circuit
    SQC_Circuit temp;
    temp.n = n;
    temp.max_m = max_m;
    temp.Construct();
    for(int i = 0; i < m; i++) {
        temp.AddOperator(operator_list[i]);
    }
    // Clear this circuit
    Clear();
    // Resize this circuit if necessary
    if((m-n_rep+in_new_ops->m)>(max_m)) {
        Resize(m-n_rep+in_new_ops->m);
        cout << "Circuit expanded by " << (m-n_rep+in_new_ops->m-max_m) << endl;
    }

    // Copy original operators back into circuit up to t
    for(int i = 0; i < t; i++) {
        AddOperator(temp.operator_list[i]);
    }
    // Copy new operators from in_new_ops
    for(int i = 0; i < in_new_ops->m; i++) {
        AddOperator(in_new_ops->operator_list[i]);
    }
    // Copy rest from temp starting at t+n_rep
    for(int i = (t+n_rep); i < temp.m; i++) {
        AddOperator(temp.operator_list[i]);
    }
    temp.Destruct();
}

void SQC_Circuit::ConvertFromToffoli() {
    // Pre-convert to explicit ancilla representation
    if(ancilla_mode!=SQC_ANCILLA_MODE_MANUAL) {
        SQC_Circuit temp;
        temp.Copy(*this);
        AllocateAncillas(temp);
        temp.Destruct();
        cout << "Converted to explicit ancilla mode" << endl;
        Print();
    }

    // Convert Toffoli-N to Toffoli-3's
    if(toffoli_n_mode==SQC_TOFFOLI_N_MODE_TOFF3) {
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_N) {
                SQC_Circuit this_toffoli_N;
                int n_args = 0;
                bool found = 0;
                for(int a = 1; (!found)&&(a <= n); a++) {
                    if((operator_list[i][a]>0)&&(operator_list[i][a]<=n)) {
                        n_args++;
                    } else {
                        found = 1;
                    }
                }
                int N_toff = (n_args+3)/2;
                this_toffoli_N.n = n;
                this_toffoli_N.max_m = 2*(N_toff-3)+1;
                this_toffoli_N.p = p;
                this_toffoli_N.Construct();

                int* this_toff = new int[n+1];
                for(int c = 0; c < (n+1); c++) this_toff[c] = 0;
                for(int k = 0; k < (N_toff-3); k++) {
                    this_toff[0] = SQC_OPERATOR_TOFFOLI;
                    this_toff[1] = ((N_toff + 1 + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + 1 + k)];
                    this_toff[2] = ((N_toff + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + k)];
                    this_toff[3] = ((N_toff - 1 - k)>N_toff?(n-p):0) + operator_list[i][(N_toff - 1 - k)];
                    this_toffoli_N.AddOperator(this_toff);
                }
                this_toff[0] = SQC_OPERATOR_TOFFOLI;
                this_toff[1] = operator_list[i][1];
                this_toff[2] = operator_list[i][2];
                this_toff[3] = ((2*N_toff - 3)>N_toff?(n-p):0) + operator_list[i][(2*N_toff - 3)];
                this_toffoli_N.AddOperator(this_toff);
                for(int k = 0; k < (N_toff-3); k++) {
                    this_toff[0] = SQC_OPERATOR_TOFFOLI;
                    this_toff[1] = ((2*N_toff - 3 - k)>N_toff?(n-p):0) + operator_list[i][(2*N_toff - 3 - k)];
                    this_toff[2] = ((2*N_toff - 4 - k)>N_toff?(n-p):0) + operator_list[i][(2*N_toff - 4 - k)];
                    this_toff[3] = ((3 + k)>N_toff?(n-p):0) + operator_list[i][(3 + k)];
                    this_toffoli_N.AddOperator(this_toff);
                }
                delete [] this_toff;

                ReplaceOperator(&this_toffoli_N,i);
                this_toffoli_N.Destruct();
            }
        }
    } else if(toffoli_n_mode==SQC_TOFFOLI_N_MODE_JONES) {
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_N) {
                SQC_Circuit this_toffoli_N;
                int n_args = GetNArgs(i);
                int N_toff = (n_args+3)/2;

                this_toffoli_N.n = n;
                this_toffoli_N.max_m = 3*fmax(0,N_toff-3)+1;
                this_toffoli_N.p = p;
                this_toffoli_N.Construct();

                int* this_operator = new int[n+1];
                for(int c = 0; c < (n+1); c++) this_operator[c] = 0;

                if(N_toff>=3) {
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = ((N_toff + 1 + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + 1 + k)];
                        this_operator[2] = ((N_toff + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = ((N_toff - 1 - k)>N_toff?(n-p):0) + operator_list[i][(N_toff - 1 - k)];
                        this_toffoli_N.AddOperator(this_operator);
                        this_operator[0] = SQC_OPERATOR_CS;
                        this_operator[1] = ((N_toff - 1 - k)>N_toff?(n-p):0) + operator_list[i][(N_toff - 1 - k)];
                        this_operator[2] = ((N_toff + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = 0;
                        this_toffoli_N.AddOperator(this_operator);
                    }
                    this_operator[0] = SQC_OPERATOR_TOFFOLI;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_operator[3] = ((2*N_toff - 3)>N_toff?(n-p):0) + operator_list[i][(2*N_toff - 3)];
                    this_toffoli_N.AddOperator(this_operator);
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_HADAMARD;
                        this_operator[1] = ((2*N_toff - 3 - k)>N_toff?(n-p):0) + operator_list[i][(2*N_toff - 3 - k)];
                        this_operator[2] = 0;
                        this_operator[3] = 0;
                        this_toffoli_N.AddOperator(this_operator);
                    }
                } else if(N_toff == 2) {
                    this_operator[0] = SQC_OPERATOR_CNOT;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_toffoli_N.AddOperator(this_operator);
                } else if(N_toff == 1) {
                    this_operator[0] = SQC_OPERATOR_HADAMARD;
                    this_operator[1] = operator_list[i][1];
                    this_toffoli_N.AddOperator(this_operator);
                    this_operator[0] = SQC_OPERATOR_Z;
                    this_operator[1] = operator_list[i][1];
                    this_toffoli_N.AddOperator(this_operator);
                    this_operator[0] = SQC_OPERATOR_HADAMARD;
                    this_operator[1] = operator_list[i][1];
                    this_toffoli_N.AddOperator(this_operator);
                }
                delete [] this_operator;

                ReplaceOperator(&this_toffoli_N,i);
                this_toffoli_N.Destruct();
            }
        }
    }
    cout << "Converted Toffoli-N" << endl;
    Print();
    // Convert Toffoli-4 to Toffoli-3's
    for(int i = 0; i < m; i++) {
        if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_4) {
            SQC_Circuit this_toffoli_4;
            this_toffoli_4.n = n;
            this_toffoli_4.max_m = 3;
            this_toffoli_4.p = p;
            this_toffoli_4.Construct();
            this_toffoli_4.AddOperator("Toffoli 1 2 3");
            this_toffoli_4.AddOperator("Toffoli 1 2 3");
            this_toffoli_4.AddOperator("Toffoli 1 2 3");
            int ancilla_ind = operator_list[i][5]+n-p;
            this_toffoli_4.operator_list[0][1] = ancilla_ind;
            this_toffoli_4.operator_list[0][2] = operator_list[i][2];
            this_toffoli_4.operator_list[0][3] = operator_list[i][3];
            this_toffoli_4.operator_list[1][1] = operator_list[i][1];
            this_toffoli_4.operator_list[1][2] = operator_list[i][4];
            this_toffoli_4.operator_list[1][3] = ancilla_ind;
            this_toffoli_4.operator_list[2][1] = ancilla_ind;
            this_toffoli_4.operator_list[2][2] = operator_list[i][2];
            this_toffoli_4.operator_list[2][3] = operator_list[i][3];
            ReplaceOperator(&this_toffoli_4,i);
            this_toffoli_4.Destruct();
        }
    }
    cout << "Converted Toffoli-4" << endl;
    Print();
    // Convert Toffoli-3 to Hadamards and CCZ
    for(int i = 0; i < m; i++) {
        if(operator_list[i][0]==SQC_OPERATOR_TOFFOLI) {
            SQC_Circuit this_toffoli;
            this_toffoli.n = n;
            this_toffoli.max_m = 3;
            this_toffoli.Construct();
            this_toffoli.AddOperator("H 1");
            this_toffoli.AddOperator("CCZ 1 2 3");
            this_toffoli.AddOperator("H 1");
            this_toffoli.operator_list[0][1] = this_toffoli.operator_list[2][1] = operator_list[i][1];
            this_toffoli.operator_list[1][1] = operator_list[i][1];
            this_toffoli.operator_list[1][2] = operator_list[i][2];
            this_toffoli.operator_list[1][3] = operator_list[i][3];
            ReplaceOperator(&this_toffoli,i);
            this_toffoli.Destruct();
        }
    }
    // Cancel adjacent Hadamards
    cout << "Converted Toffoli to {H,CCZ}" << endl;
    Print();

    //cout << "Out of print." << endl;

    for(int i = 0; i < (m-1); i++) {
        //cout << "i = " << i << endl;
        if(operator_list[i][0]==SQC_OPERATOR_HADAMARD) {
            int acting_qubit = operator_list[i][1];
            bool qubit_used = 0;
            for(int j = (i+1); (!qubit_used)&&(j < m); j++) {
                switch(operator_list[j][0]) {
                    case SQC_OPERATOR_HADAMARD:
                        if((!qubit_used)&&(operator_list[j][1]==acting_qubit)) {
                            DeleteOperator(j);
                            DeleteOperator(i);
                            i--;
                        }
                        break;
                    case SQC_OPERATOR_S:
                    case SQC_OPERATOR_Z:
                    case SQC_OPERATOR_T:
                        if(operator_list[j][1]==acting_qubit) {
                            qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CS:
                    case SQC_OPERATOR_CNOT:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)) {
                            qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CCZ:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)||(operator_list[j][3]==acting_qubit)) {
                            qubit_used = 1;
                        }
                        break;
                }
            }
        }
    }
    cout << "Cancelled adjacent Hadamards" << endl;
    Print();
}

void SQC_Circuit::AllocateAncillas(const SQC_Circuit& in_C) {
    // Ensure this circuit is empty
    Destruct();

    // Check which ancilla mode
    switch(in_C.ancilla_mode) {
        case SQC_ANCILLA_MODE_MANUAL:
            // Simply copy this circuit into in_C
            Copy(in_C);
            break;
        case SQC_ANCILLA_MODE_PER_GATE:
            {
                // Determine the maximum N for each Toffoli-N gate in the circuit
                int N_ancillas = 0;
                int N_args = 3;
                for(int i = 0; i < in_C.m; i++) {
                    if(in_C.operator_list[i][0]==SQC_OPERATOR_TOFFOLI_N) {
                        bool found = 0;
                        int a = 0;
                        while((!found)&&(a<in_C.n)) {
                            if((in_C.operator_list[i][a+1]>0)&&(in_C.operator_list[i][a+1]<=in_C.n))
                            a++;
                            else
                            found = 1;
                        }
                        if(a>N_args) {
                            N_args = a;
                        }
                    }
                }
                N_ancillas = fmax(0,N_args-3);
                cout << "N_ancillas = " << N_ancillas << endl;
                n = in_C.n + N_ancillas;
                p = N_ancillas;
                m = 0;
                max_m = in_C.max_m;
                ancilla_mode = SQC_ANCILLA_MODE_MANUAL;
                Construct();
                // Go through each gate in in_C, adding ancilla indices to each ensure they are unique for each gate
                for(int i = 0; i < in_C.m; i++) {
                    if(in_C.operator_list[i][0]==SQC_OPERATOR_TOFFOLI_N) {
                        bool found = 0;
                        int a = 0;
                        while((!found)&&(a<in_C.n)) {
                            if((in_C.operator_list[i][a+1]>0)&&(in_C.operator_list[i][a+1]<=in_C.n))
                            a++;
                            else
                            found = 1;
                        }
                        int this_N_ancillas = fmax(0,a-3);
                        AddOperator(in_C.operator_list[i]);
                        for(int j = 1; j <= this_N_ancillas; j++) {
                            operator_list[i][a+j] = j;
                        }
                    }
                }
            }
            break;
        case SQC_ANCILLA_MODE_PER_CIRCUIT:
            {
                // Go through each gate in in_C, adding ancilla indices to each ensure they are unique for entire circuit
                int N_ancillas = 0;
                for(int i = 0; i < in_C.m; i++) {
                    if(in_C.operator_list[i][0]==SQC_OPERATOR_TOFFOLI_N) {
                        bool found = 0;
                        int a = 0;
                        while((!found)&&(a<in_C.n)) {
                            if((in_C.operator_list[i][a+1]>0)&&(in_C.operator_list[i][a+1]<=in_C.n))
                            a++;
                            else
                            found = 1;
                        }
                        N_ancillas += fmax(0,a-3);
                    }
                }
                cout << "N_ancillas = " << N_ancillas << endl;
            }
            break;
    }
}

int SQC_Circuit::GetNArgs(int i) const {
    bool found = 0;
    int a = 0;
    while((!found)&&(a<n)) {
        if((operator_list[i][a+1]>0)&&(operator_list[i][a+1]<=n))
        a++;
        else
        found = 1;
    }
    return a;
}

void SQC_Circuit::Resize(int in_max_m) {
    SQC_Circuit temp;
    temp.Copy(*this);
    Destruct();
    max_m = in_max_m;
    m = 0;
    Construct();
    for(int i = 0; i < temp.m; i++) {
        AddOperator(temp.operator_list[i]);
    }
    temp.Destruct();
}

void SQC_Circuit::LoadMaslovFile(const char* in_filename) {
    ifstream my_file(in_filename);
    if(my_file.good()) {
        Destruct();
        int linewidth = 1000;
        char* this_line = new char[linewidth];
        my_file.getline(this_line,linewidth);
        char* this_tok = NULL;
        this_tok = strtok(this_line,",");
        int this_n = 0;
        if(this_tok) this_n = 1;
        while(this_tok = strtok(NULL,",")) this_n++;
        cout << "n = " << this_n << endl;
        n = this_n;
        max_m = SQC_DEFAULT_MAX_M;
        m = 0;
        p = -1;
        ancilla_mode = SQC_ANCILLA_MODE_PER_GATE;
        Construct();
        my_file.getline(this_line,linewidth);
        my_file.getline(this_line,linewidth);
        my_file.getline(this_line,linewidth);
        my_file.getline(this_line,linewidth);
        int* this_operator = new int[n+1];
        this_operator[0] = SQC_OPERATOR_TOFFOLI_N;
        while(!my_file.eof()) {
            if(strlen(this_line)&&(this_line[0]=='t')) {
                char* this_args;
                this_args = strchr(this_line,' ');
                this_args+=1;
                cout << "this_args = " << this_args << endl;
                for(int i = 1; i < (n+1); i++) this_operator[i] = 0;
                char* this_arg = NULL;
                this_arg = strtok(this_args,",");
                if(this_arg) this_operator[1] = this_arg[0] - 'a' + 1;
                int c = 2;
                while((this_arg = strtok(NULL,","))&&(c<=n)) {
                    this_operator[c] = this_arg[0] - 'a' + 1;
                    c++;
                }
                AddOperator(this_operator);
            }
            my_file.getline(this_line,linewidth);
        }
        delete [] this_operator;

        delete [] this_line;
    }
    my_file.close();
}
