#include "SQC_Circuit.h"

#include <iostream>
using namespace std;
#include <ostream>
#include <string>
#include <cstring>
#include <fstream>
#include <utility>
#include <cmath>
#include <cstdlib>

#include "Signature.h"
#include "BMSparse.h"
#include "LukeMat_GF2.h"
#include "Interface_SigBMS.h"
#include "GateSynthesisMatrix.h"
#include "GateStringSparse.h"
#include "GateSigInterface.h"
#include "Interface_BMSGSS.h"
#include "LukeConsoleOut.h"
#include "Matrix.h"
#include "TO_Maps.h"
using namespace LukeConsoleOut;
#include "Utils.h"
using namespace Utils;


SQC_Circuit::SQC_Circuit() {
    ;
}

SQC_Circuit::SQC_Circuit(int in_n) {
    n = in_n;
    Construct();
}

SQC_Circuit::SQC_Circuit(int in_n, int in_d) {
    n = in_n;
    d = in_d;
    Construct();
}

SQC_Circuit::SQC_Circuit(int in_n, int in_d, int in_p) {
    n = in_n;
    d = in_d;
    p_hads = in_p;
    Construct();
}

SQC_Circuit::~SQC_Circuit() {
    Destruct();
}

SQC_Circuit::SQC_Circuit(const SQC_Circuit& in_C) {
    n = in_C.n;
    d = in_C.d;
    m = 0;
    max_m = in_C.max_m;
    hadamard_mode_max_ancillas = in_C.hadamard_mode_max_ancillas;
    p_hads = in_C.p_hads;
    Construct();
    for(int t = 0; t < in_C.m; t++) {
        AddOperator(in_C.operator_list[t]);
    }
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

void SQC_Circuit::Print(ostream* in_OS, int start_i, int print_n, bool with_t) const {
    (*in_OS) << "n " << n << endl;
    (*in_OS) << "n_x " << d << endl;
    (*in_OS) << "n_y " << (n-d-p_hads) << endl;
    (*in_OS) << "n_z " << (p_hads) << endl;
    (*in_OS) << "m " << max_m << endl;
    if(print_n<0) print_n = m; else print_n = (start_i+print_n);
    for(int i = start_i; i < fmin(m,print_n); i++) {
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
        case SQC_OPERATOR_X:
            this_op_str = SQC_OPSTRING_X;
            break;
        case SQC_OPERATOR_Y:
            this_op_str = SQC_OPSTRING_Y;
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
        case SQC_OPERATOR_PARTITION:
            this_op_str = SQC_OPSTRING_PARTITION;
            break;
        case SQC_OPERATOR_POST_0:
            this_op_str = SQC_OPSTRING_POST_0;
            break;
        case SQC_OPERATOR_CZ:
            this_op_str = SQC_OPSTRING_CZ;
            break;
        case SQC_OPERATOR_T_DAG:
            this_op_str = SQC_OPSTRING_T_DAG;
            break;
        case SQC_OPERATOR_CS_DAG:
            this_op_str = SQC_OPSTRING_CS_DAG;
            break;
        case SQC_OPERATOR_S_DAG:
            this_op_str = SQC_OPSTRING_S_DAG;
            break;
        }
        if(with_t) (*in_OS) << "[" << i << "] = ";
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
        Resize(max_m+SQC_CIRCUIT_EXPAND);
        //LOut() << "Circuit expanded by " << SQC_CIRCUIT_EXPAND << endl;
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
                } else if(!strcmp(this_tok, SQC_OPSTRING_PARTITION)) {
                    operator_list[m][0] = SQC_OPERATOR_PARTITION;
                } else if(!strcmp(this_tok, SQC_OPSTRING_X)) {
                    operator_list[m][0] = SQC_OPERATOR_X;
                } else if(!strcmp(this_tok, SQC_OPSTRING_Y)) {
                    operator_list[m][0] = SQC_OPERATOR_Y;
                } else if(!strcmp(this_tok, SQC_OPSTRING_POST_0)) {
                    operator_list[m][0] = SQC_OPERATOR_POST_0;
                } else if(!strcmp(this_tok, SQC_OPSTRING_CZ)) {
                    operator_list[m][0] = SQC_OPERATOR_CZ;
                } else if(!strcmp(this_tok, SQC_OPSTRING_T_DAG)) {
                    operator_list[m][0] = SQC_OPERATOR_T_DAG;
                } else if(!strcmp(this_tok, SQC_OPSTRING_CS_DAG)) {
                    operator_list[m][0] = SQC_OPERATOR_CS_DAG;
                } else if(!strcmp(this_tok, SQC_OPSTRING_S_DAG)) {
                    operator_list[m][0] = SQC_OPERATOR_S_DAG;
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
                    case SQC_OPERATOR_CZ:
                    case SQC_OPERATOR_CS_DAG:
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
                if(!success) LOut() << "WARNING! Incorrect number of qubits." << endl;
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
        if(n==0) LOut() << "Warning! n should take non-zero value. n = " << n << endl;

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

    Print(false,0,-1,&my_file);

    my_file.close();
}

void SQC_Circuit::Clear() {
    m = 0;
}

bool SQC_Circuit::GetPartition(SQC_Circuit* out_Hadamards, SQC_Circuit* out_CNOT_T) {
    //cout << "QWER" << endl;
    bool out = false;
    if(operator_list&&out_Hadamards->operator_list&&out_CNOT_T->operator_list) {
        out = (bool)m;
        //cout << "TYUI" << endl;

    if(out) {
        // Pre-process Toffoli gates by converting them to CCZ conjugated by Hadamards on the target
        //cout << "KHGFXC" << endl;
        ConvertFromToffoli();

        //cout << "Out 1 n = " << n << endl;
        int init_m = m;
        bool* qubit_used = new bool[n];
        for(int i = 0; i < n; i++) qubit_used[i] = 0;
        bool** hadamard_met = new bool*[init_m];
        bool* hadamard_used = new bool[n];
        for(int i = 0; i < n; i++) hadamard_used[i] = 0;
        for(int i = 0; i < init_m; i++) {
            hadamard_met[i] = new bool[n];
            for(int j = 0; j < n; j++) {
                hadamard_met[i][j] = 0;
            }
        }
        //cout << "Out 2" << endl;
        for(int t = 0; t < m; t++) {
            if((t==0)&&(operator_list[t][0]==SQC_OPERATOR_PARTITION)) {
                DeleteOperator(t);
                t--;
            } else {
                //cout << "Out t " << t << endl;
                // If the operator is a Hadamard
                if(operator_list[t][0]==SQC_OPERATOR_HADAMARD) {
                    //cout << "Out Hadamard" << endl;
                    // If the qubit upon which it acts hasn't been used yet
                    if(((operator_list[t][1]-1)<0)||((operator_list[t][1]-1)>(n-1))) {
                            LOut() << "Place 1: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][1]-1) << endl;
                            Print(&LOut(),t,1);
                            return 0;
                    }
                    if(qubit_used[(operator_list[t][1]-1)]==0) {
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
                        case SQC_OPERATOR_PARTITION:
                            {
                                if(((operator_list[t][1]-1)<0)||((operator_list[t][1]-1)>(n-1))) {
                                        LOut() << "Place 2: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][1]-1) << endl;
                                        Print(&LOut(),t,1);
                                        return 0;
                                }
                                int this_nargs = GetNArgs(t);
                                for(int is = 0; is < this_nargs; is++) qubit_used[(operator_list[t][is+1]-1)]=1;
                                //cout << "This argument(" << t << ") = " << (operator_list[t][1]-1) << endl;
                                //cout << "qwe = " << (qubit_used[(operator_list[t][1]-1)]=1) << endl;
                                //delete [] qubit_used;
                                //return 0;
                            }
                            //cout << "This argument(" << t << ") = " << (operator_list[t][1]-1) << endl;
                            //cout << "n = " << n << endl;
                            //Print();
                            break;
                        case SQC_OPERATOR_CNOT:
                        case SQC_OPERATOR_CS:
                            {
                                if(((operator_list[t][1]-1)<0)||((operator_list[t][1]-1)>(n-1))) {
                                        LOut() << "Place 3a: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][1]-1) << endl;
                                        Print(&LOut(),t,1);
                                        return 0;
                                }
                                if(((operator_list[t][2]-1)<0)||((operator_list[t][2]-1)>(n-1))) {
                                        LOut() << "Place 3b: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][2]-1) << endl;
                                        Print(&LOut(),t,1);
                                        return 0;
                                }
                                qubit_used[(operator_list[t][1]-1)]=1;
                                qubit_used[(operator_list[t][2]-1)]=1;
                            }

                            break;
                        case SQC_OPERATOR_CCZ:
                            {
                                if(((operator_list[t][1]-1)<0)||((operator_list[t][1]-1)>(n-1))) {
                                        LOut() << "Place 4a: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][1]-1) << endl;
                                        Print(&LOut(),t,1);
                                        return 0;
                                }
                                if(((operator_list[t][2]-1)<0)||((operator_list[t][2]-1)>(n-1))) {
                                        LOut() << "Place 4b: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][2]-1) << endl;
                                        Print(&LOut(),t,1);;
                                        return 0;
                                }
                                if(((operator_list[t][3]-1)<0)||((operator_list[t][3]-1)>(n-1))) {
                                        LOut() << "Place 4c: AAAAAAGHH!!! t = " << t << "; a = " << (operator_list[t][3]-1) << endl;
                                        Print(&LOut(),t,1);
                                        return 0;
                                }
                                qubit_used[(operator_list[t][1]-1)]=1;
                                qubit_used[(operator_list[t][2]-1)]=1;
                                qubit_used[(operator_list[t][3]-1)]=1;
                            }
                            break;
                        }
                    }
                }
            }
            //LOut() << "out_Hadamards" << endl;
            //out_Hadamards->Print();
            //delete [] qubit_used;
            //return 0;
            //for(int pk = 0; pk < n; pk++) {
            //    cout << "qubit_used[" << pk << "] = " << qubit_used[pk] << endl;
            //}
            //qubit_used = NULL;
            //cout << "Out 3" << endl;
            for(int t = 0; t < m; t++) {
                if((t==0)&&(operator_list[t][0]==SQC_OPERATOR_PARTITION)) {
                    DeleteOperator(t);
                    t--;
                } else {
                    // If the operator is a Hadamard
                    if((operator_list[t][0]==SQC_OPERATOR_HADAMARD)||(operator_list[t][0]==SQC_OPERATOR_PARTITION)) {
                        int this_nargs = GetNArgs(t);
                        bool old_hadamard_met = 1;
                        for(int is = 0; is < this_nargs; is++) old_hadamard_met *= hadamard_used[operator_list[t][is+1]-1];
                        for(int i = t; i < m; i++) {
                            hadamard_met[i][operator_list[t][1]-1]=1;
                        }

                        for(int is = 0; is < this_nargs; is++) hadamard_used[operator_list[t][is+1]-1] = 1;
                        //if((!old_hadamard_met)&&(operator_list[t][0]==SQC_OPERATOR_PARTITION)) {
                            //DeleteOperator(t);
                            //t--;
                        //}
                    } else {
                        // If non of the qubits upon which the operator acts have been acted upon by a Hadamard (except 'free' Hadamards)
                        bool can_move_left = 0;
                        switch(operator_list[t][0]) {
                            case SQC_OPERATOR_IDENTITY:
                            case SQC_OPERATOR_T:
                            case SQC_OPERATOR_S:
                            case SQC_OPERATOR_Z:
                                //can_move_left = !hadamard_met[t][operator_list[t][1]-1];
                                can_move_left = !hadamard_used[operator_list[t][1]-1];
                                break;
                            case SQC_OPERATOR_CNOT:
                            case SQC_OPERATOR_CS:
                                //can_move_left = !hadamard_met[t][operator_list[t][1]-1]&&!hadamard_met[t][operator_list[t][2]-1];
                                can_move_left = (!hadamard_used[operator_list[t][1]-1])&&(!hadamard_used[operator_list[t][2]-1]);
                                break;
                            case SQC_OPERATOR_CCZ:
                                /*can_move_left = !hadamard_met[t][operator_list[t][1]-1]
                                             && !hadamard_met[t][operator_list[t][2]-1]
                                             && !hadamard_met[t][operator_list[t][3]-1];*/
                                can_move_left = (!hadamard_used[operator_list[t][1]-1])&&(!hadamard_used[operator_list[t][2]-1])&&(!hadamard_used[operator_list[t][3]-1]);
                                break;
                        }
                        if(can_move_left) {
                            //cout << "Can partition." << endl;
                            //Print(&cout,t,1);
                            out_CNOT_T->AddOperator(operator_list[t]);
                            //LOut() << "out_CNOT_T" << endl;
                            //out_CNOT_T->Print();
                            DeleteOperator(t);
                            t--;
                        } else {
                            //cout << "Can't partition." << endl;
                            //Print(&cout,t,1);
                        }
                    }
                }
            }
            //cout << "Out 4 init_m = " << init_m << endl;
            for(int i = 0; i < init_m; i++) {
                delete [] hadamard_met[i];
                hadamard_met[i] = NULL;
            }
            //cout << "Out 4.1" << endl;
            delete [] hadamard_met;
            hadamard_met = NULL;
            //cout << "Out 4.2 qubit_used = " << qubit_used << endl;

            delete [] qubit_used;
            qubit_used = NULL;

            delete [] hadamard_used;

            //cout << "Out 4.3" << endl;
            out = (bool)m;
            //cout << "Out 5" << endl;
        }

    }
    return out;
}

bool SQC_Circuit::GetPartition(SQC_Circuit* out) {
    bool outb = 0;

    //cout << "JYGDSZF" << endl;
    ConvertFromToffoli();
    //cout << "vbmhgcyj" << endl;
    //LOut() << "This circuit" << endl;
    //Print();
    //LOut() << endl;
    //cout << "VBTRSTRG" << endl;
    ConvertHadamard(out,hadamard_mode_max_ancillas);
    //cout << "vbnxt" << endl;
    outb = (bool)m;

    return outb;
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

void SQC_Circuit::AddOperator(const SQC_Operator in_op, int len) {
    // Resize if necessary
    if(m>=max_m) {
        Resize(max_m+SQC_CIRCUIT_EXPAND);
        //LOut() << "Circuit expanded by " << SQC_CIRCUIT_EXPAND << endl;
    }
    if(m<max_m) {
        for(int i = 0; i < (n+1); i++) {
            if((len>0)&&(i>=len)) {
                operator_list[m][i] = 0;
            } else {
                operator_list[m][i] = in_op[i];
            }
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
        LOut() << "Error in SQC_Circuit::DecompositionVW. At least one SQC_Circuit has not been intialized." << endl;
    }
    // TODO remove redundant CNOTs from output W circuit.
}

BMSparse SQC_Circuit::toGateSynthesisMatrix() const {
    //cout << "tGSM start" << endl;
    Signature out(n);

    BMSparse A(n,0);
    BMSparse P = BMSparse::eye(n,n);
    //cout << "P" << endl;
    //P.printFull();

    for(int t = 0; t < m; t++) {
        //cout << "t = " << t << endl;
        //cout << "n = " << n << endl;
        //for(int c = 0; c < (n+1); c++) cout << (int)operator_list[t][c];
        //cout << endl;
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
            //temp_sig.print();
            GateStringSparse new_cols_GSS = GateSigInterface::SigToGSS(temp_sig);

            BMSparse new_cols = Interface_BMSGSS::GSSToBMS(new_cols_GSS);
            //cout << "new cols" << endl;
            //new_cols.printFull();
            new_cols = P.T()*new_cols;
            //cout << "new cols'" << endl;
            //new_cols.printFull();
            A = A && new_cols;
        }

        //cout << "end iter" << endl;
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
    //cout << "tGSM end" << endl;
    return A;
}

bool SQC_Circuit::NextSignature(Signature& outsig) {
    //cout << "QWER" << endl;
    outsig = Signature(n);
    bool out = 0;
    SQC_Circuit hadamards, CNOT_Ts;
    switch(hadamard_mode) {
        case SQC_HADAMARD_MODE_PARTITION:
            {
                hadamards.n = CNOT_Ts.n = n;
                hadamards.max_m = CNOT_Ts.max_m = max_m;
                hadamards.Construct();
                CNOT_Ts.Construct();
                //cout << "xcvxcv" << endl;
                out = GetPartition(&hadamards, &CNOT_Ts);
                //LOut() << "Hadamards" << endl;
                //hadamards.Print();
                //LOut() << endl;

                //cout << "zxcv" << endl;
                hadamards.Destruct();
            }
            break;
        case SQC_HADAMARD_MODE_MONTANARO: {
            {
                //cout << "SGHGSA" << endl;
                out = GetPartition(&CNOT_Ts);
                //cout << "jhjgh" << endl;
            }
            break;
        }
    }

    //LOut() << "This circuit" << endl;
    //Print();
    //LOut() << endl;
    //LOut() << "{CNOT,T}" << endl;
    //CNOT_Ts.Print(&cout,0,70);
    //LOut() << endl;
    if((bool)CNOT_Ts.m) {
        //cout << "TYUI" << endl;
        SQC_Circuit this_V, this_W;
        this_V.n = this_W.n = CNOT_Ts.n;
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
        //cout << "ZCXG" << endl;
    }
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
        //LOut() << "Circuit expanded by " << (m-n_rep+in_new_ops->m-max_m) << endl;
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
        LOut() << "Converted to explicit ancilla mode" << endl;
        //temp.Print();
        temp.Destruct();
        //return;
    }

    // Convert Toffoli-N to Toffoli-3's
    if(toffoli_n_mode==SQC_TOFFOLI_N_MODE_TOFF3) {
        //LOut() << "Yay" << endl;
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_N) {
                SQC_Circuit this_toffoli_N;
                int n_args = GetNArgs(i);
                int N_toff;
                N_toff = n_args;
                if(n_args>3) N_toff = (n_args+3)/2;

                this_toffoli_N.n = n;
                this_toffoli_N.max_m = 3*fmax(0,N_toff-3)+1;
                this_toffoli_N.p = p;
                this_toffoli_N.Construct();

                int* this_operator = new int[n+1];
                for(int c = 0; c < (n+1); c++) this_operator[c] = 0;

                if(N_toff>3) {
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = (n-p) + operator_list[i][(N_toff + 1 + k)];
                        this_operator[2] = (((N_toff + k)>N_toff)?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = operator_list[i][(N_toff - 1 - k)];
                        this_toffoli_N.AddOperator(this_operator);
                    }
                    this_operator[0] = SQC_OPERATOR_TOFFOLI;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_operator[3] = (n-p) + operator_list[i][(2*N_toff - 3)];
                    this_toffoli_N.AddOperator(this_operator);
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = (n-p) + operator_list[i][(2*N_toff - 3-k)];
                        this_operator[2] = (((2*N_toff - 4 - k)>N_toff)?(n-p):0) + operator_list[i][(2*N_toff - 4 - k)];
                        this_operator[3] = operator_list[i][(3 + k)];
                        this_toffoli_N.AddOperator(this_operator);
                    }
                } else if(N_toff == 3) {
                    this_operator[0] = SQC_OPERATOR_TOFFOLI;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_operator[3] = operator_list[i][3];
                    this_toffoli_N.AddOperator(this_operator);
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
                } /*else if(N_toff <=0) {
                    this_operator[0] = SQC_OPERATOR_IDENTITY;
                    this_toffoli_N.AddOperator(this_operator);
                }*/
                delete [] this_operator;

                ReplaceOperator(&this_toffoli_N,i);
                this_toffoli_N.Destruct();
            }
        }
    } else if(toffoli_n_mode==SQC_TOFFOLI_N_MODE_JONES) {
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_N) {
                SQC_Circuit this_toffoli_N;
                int n_args = GetNArgs(i);
                int N_toff;
                N_toff = n_args;
                if(n_args>3) N_toff = (n_args+3)/2;

                this_toffoli_N.n = n;
                this_toffoli_N.max_m = 3*fmax(0,N_toff-3)+1;
                this_toffoli_N.p = p;
                this_toffoli_N.Construct();

                int* this_operator = new int[n+1];
                for(int c = 0; c < (n+1); c++) this_operator[c] = 0;

                if(N_toff>=3) {
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = (n-p) + operator_list[i][(N_toff + 1 + k)];
                        this_operator[2] = ((N_toff + k)>N_toff?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = operator_list[i][(N_toff - 1 - k)];
                        this_toffoli_N.AddOperator(this_operator);
                        this_operator[0] = SQC_OPERATOR_CS;
                        this_operator[1] = operator_list[i][(N_toff - 1 - k)];
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
                        this_operator[0] = SQC_OPERATOR_PARTITION;
                        for(int c = 1; c <= n; c++) this_operator[c] = c;
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
                    if((this_operator[1]<=0)||(this_operator[1]>(this_toffoli_N.n))) {
                        LOut() << "CPlace 4: AAAAAAGHH!!!" << endl;
                        this_toffoli_N.Print();
                        return;
                    }
                    this_operator[0] = SQC_OPERATOR_Z;
                    this_operator[1] = operator_list[i][1];
                    this_toffoli_N.AddOperator(this_operator);
                    if((this_operator[1]<=0)||(this_operator[1]>(this_toffoli_N.n))) {
                        LOut() << "CPlace 5: AAAAAAGHH!!!" << endl;
                        this_toffoli_N.Print();
                        return;
                    }
                    this_operator[0] = SQC_OPERATOR_HADAMARD;
                    this_operator[1] = operator_list[i][1];
                    this_toffoli_N.AddOperator(this_operator);
                    if((this_operator[1]<=0)||(this_operator[1]>(this_toffoli_N.n))) {
                        LOut() << "CPlace 6: AAAAAAGHH!!!" << endl;
                        this_toffoli_N.Print();
                        return;
                    }
                }
                delete [] this_operator;

                ReplaceOperator(&this_toffoli_N,i);
                this_toffoli_N.Destruct();
            }
        }
    } else if(toffoli_n_mode==SQC_TOFFOLI_N_MODE_PSCJ) {
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_N) {
                SQC_Circuit this_toffoli_N;
                int n_args = GetNArgs(i);
                int N_toff;
                N_toff = n_args;
                if(n_args>3) N_toff = (n_args+3)/2;

                this_toffoli_N.n = n;
                this_toffoli_N.max_m = 3*fmax(0,N_toff-3)+1;
                this_toffoli_N.p = p;
                this_toffoli_N.Construct();

                int* this_operator = new int[n+1];
                for(int c = 0; c < (n+1); c++) this_operator[c] = 0;

                if(N_toff>3) {
                    for(int k = 0; k < (N_toff-3); k++) {
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = (n-p) + operator_list[i][(N_toff + 1 + k)];
                        this_operator[2] = (((N_toff + k)>N_toff)?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = operator_list[i][(N_toff - 1 - k)];
                        this_toffoli_N.AddOperator(this_operator);
                        this_operator[0] = SQC_OPERATOR_CS;
                        this_operator[1] = operator_list[i][(N_toff - 1 - k)];
                        this_operator[2] = (((N_toff + k)>N_toff)?(n-p):0) + operator_list[i][(N_toff + k)];
                        this_operator[3] = 0;
                        this_toffoli_N.AddOperator(this_operator);
                        if(hadamard_mode==SQC_HADAMARD_MODE_PARTITION) {
                            this_operator[0] = SQC_OPERATOR_PARTITION;
                            this_operator[1] = operator_list[i][(N_toff - 1 - k)];
                            this_operator[2] = (((N_toff + k)>N_toff)?(n-p):0) + operator_list[i][(N_toff + k)];
                            this_operator[3] = 0;
                            this_toffoli_N.AddOperator(this_operator);
                        }
                    }

                    /*{
                        this_operator[0] = SQC_OPERATOR_PARTITION;
                        int ih = 1;
                        for(int c = 0; c < (N_toff-3); c++) {
                            this_operator[ih] = operator_list[i][(N_toff - 1 - c)];
                            ih++;
                        }
                        for(int c = 0; c < (N_toff-3); c++) {
                            this_operator[ih] = (((N_toff + c)>N_toff)?(n-p):0) + operator_list[i][(N_toff + c)];
                            ih++;
                        }
                        this_toffoli_N.AddOperator(this_operator);
                        for(int c = 0; c < (n+1); c++) this_operator[c] = 0;
                    }*/

                    this_operator[0] = SQC_OPERATOR_TOFFOLI;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_operator[3] = (n-p) + operator_list[i][(2*N_toff - 3)];
                    this_toffoli_N.AddOperator(this_operator);

                    if(hadamard_mode==SQC_HADAMARD_MODE_PARTITION) {
                        for(int c = 0; c < (n+1); c++) this_operator[c] = 0;
                        this_operator[0] = SQC_OPERATOR_PARTITION;
                    }
                    /*{
                        int ih = 1;
                        for(int c = 0; c < (N_toff-3); c++) {
                            this_operator[ih] = operator_list[i][(3 + c)];
                            ih++;
                        }
                        for(int c = 0; c < (N_toff-3); c++) {
                            this_operator[ih] = (((2*N_toff - 4 - c)>N_toff)?(n-p):0) + operator_list[i][(2*N_toff - 4 - c)];
                            ih++;
                        }
                        this_toffoli_N.AddOperator(this_operator);
                        for(int c = 0; c < (n+1); c++) this_operator[c] = 0;
                    }*/

                    for(int k = 0; k < (N_toff-3); k++) {
                        if(hadamard_mode==SQC_HADAMARD_MODE_PARTITION) {
                            this_operator[0] = SQC_OPERATOR_PARTITION;
                            this_operator[1] = operator_list[i][(3 + k)];
                            this_operator[2] = (((2*N_toff - 4 - k)>N_toff)?(n-p):0) + operator_list[i][(2*N_toff - 4 - k)];
                            this_operator[3] = 0;
                            this_toffoli_N.AddOperator(this_operator);
                        }
                        this_operator[0] = SQC_OPERATOR_CS;
                        this_operator[1] = operator_list[i][(3 + k)];
                        this_operator[2] = (((2*N_toff - 4 - k)>N_toff)?(n-p):0) + operator_list[i][(2*N_toff - 4 - k)];
                        this_operator[3] = 0;
                        this_toffoli_N.AddOperator(this_operator);
                        this_operator[0] = SQC_OPERATOR_TOFFOLI;
                        this_operator[1] = (n-p) + operator_list[i][(2*N_toff - 3-k)];
                        this_operator[2] = (((2*N_toff - 4 - k)>N_toff)?(n-p):0) + operator_list[i][(2*N_toff - 4 - k)];
                        this_operator[3] = operator_list[i][(3 + k)];
                        this_toffoli_N.AddOperator(this_operator);
                    }
                } else if(N_toff == 3) {
                    this_operator[0] = SQC_OPERATOR_TOFFOLI;
                    this_operator[1] = operator_list[i][1];
                    this_operator[2] = operator_list[i][2];
                    this_operator[3] = operator_list[i][3];
                    this_toffoli_N.AddOperator(this_operator);
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
                } /*else if(N_toff <=0) {
                    this_operator[0] = SQC_OPERATOR_IDENTITY;
                    this_toffoli_N.AddOperator(this_operator);
                }*/
                delete [] this_operator;

                ReplaceOperator(&this_toffoli_N,i);
                this_toffoli_N.Destruct();
            }
        }
    }
    LOut() << "Converted Toffoli-N" << endl;
    //Print(&cout,0,27);
    //Print();
    // Convert Toffoli-4 to Toffoli-3's
    for(int i = 0; i < m; i++) {
        if(operator_list[i][0] == SQC_OPERATOR_TOFFOLI_4) {
            LOut() << "Warning! Toffoli-4 deprocated. Use Toffoli-4." << endl;
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
    //LOut() << "Converted Toffoli-4" << endl;
    //Print();
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
    LOut() << "Converted Toffoli to {H,CCZ}" << endl;
    //Print(&cout,0,70);
    //Print();

    //cout << "Out of print." << endl;
    CancelAdjacentHadamards();
    /*int n_deleted = 0;
    for(int i = 0; i < (m-1); i++) {
        //cout << "i = " << i << endl;
        if(operator_list[i][0]==SQC_OPERATOR_HADAMARD) {
            int acting_qubit = operator_list[i][1];
            bool this_qubit_used = 0;
            for(int j = (i+1); (!this_qubit_used)&&(j < m); j++) {
                switch(operator_list[j][0]) {
                    case SQC_OPERATOR_HADAMARD:
                        if((!this_qubit_used)&&(operator_list[j][1]==acting_qubit)) {
                            //cout << "Deleting operator " << j << endl;
                            //for(int c = 0; c < (GetNArgs(j)+1); c++) cout << operator_list[j][c] << " ";
                            //cout << endl;
                            //cout << "Deleting operator " << i << endl;
                            //for(int c = 0; c < (GetNArgs(i)+1); c++) cout << operator_list[i][c] << " ";
                            //cout << endl;
                            DeleteOperator(j);
                            DeleteOperator(i);
                            i--;
                            n_deleted++;
                            n_deleted++;
                            this_qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_S:
                    case SQC_OPERATOR_Z:
                    case SQC_OPERATOR_T:
                        if(operator_list[j][1]==acting_qubit) {
                            this_qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CS:
                    case SQC_OPERATOR_CNOT:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)) {
                            this_qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CCZ:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)||(operator_list[j][3]==acting_qubit)) {
                            this_qubit_used = 1;
                        }
                        break;
                }
            }
        }
    }*/
    //cout << "N deleted = " << n_deleted << endl;
    LOut() << "Cancelled adjacent Hadamards" << endl;
    //Print();
    //Print(&cout,0,70);
    //RemoveExternalHadamards();
    //LOut() << "Removed external Hadamards" << endl;
    //Print();
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
                //LOut() << "N_ancillas = " << N_ancillas << endl;
                n = in_C.n + N_ancillas;
                p = N_ancillas;
                m = 0;
                max_m = in_C.max_m;
                ancilla_mode = SQC_ANCILLA_MODE_MANUAL;
                Construct();
                // Go through each gate in in_C, adding ancilla indices to each ensure they are unique for each gate
                for(int i = 0; i < in_C.m; i++) {
                    if(in_C.operator_list[i][0]==SQC_OPERATOR_TOFFOLI_N) {
                        bool found = 0;/*
                        int a = 0;
                        while((!found)&&(a<in_C.n)) {
                            if((in_C.operator_list[i][a+1]>0)&&(in_C.operator_list[i][a+1]<=in_C.n))
                            a++;
                            else
                            found = 1;
                        }*/
                        int a = in_C.GetNArgs(i);
                        int this_N_ancillas = fmax(0,a-3);
                        int* this_operator = new int[n+1];
                        for(int c = 0; c < (n+1); c++) this_operator[c] = 0;
                        this_operator[0] = SQC_OPERATOR_TOFFOLI_N;
                        //AddOperator(in_C.operator_list[i]);
                        for(int j = 1; j <= a; j++) this_operator[j] = in_C.operator_list[i][j];
                        for(int j = 1; j <= this_N_ancillas; j++) {
                            this_operator[a+j] = j;
                        }
                        AddOperator(this_operator);
                        delete [] this_operator;
                    } else {
                        AddOperator(in_C.operator_list[i]);
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
                LOut() << "N_ancillas = " << N_ancillas << endl;
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
        while(!strstr(this_line,".v")) {
            my_file.getline(this_line,linewidth);
        }
        char* this_tok = NULL;
        this_tok = strtok(this_line,",");
        int this_n = 0;
        if(this_tok) this_n = 1;
        while(this_tok = strtok(NULL,",")) this_n++;
        //cout << "n = " << this_n << endl;
        n = this_n;
        max_m = SQC_DEFAULT_MAX_M;
        m = 0;
        p = -1;
        ancilla_mode = SQC_ANCILLA_MODE_PER_GATE;
        Construct();
        my_file.getline(this_line,linewidth);
        while(!strstr(this_line,"BEGIN")) {
            my_file.getline(this_line,linewidth);
        }
        my_file.getline(this_line,linewidth);
        int* this_operator = new int[n+1];
        this_operator[0] = SQC_OPERATOR_TOFFOLI_N;
        while(!my_file.eof()) {
            if(strlen(this_line)&&((this_line[0]=='t')||(this_line[0]=='T'))) {
                char* this_args;
                this_args = strchr(this_line,' ');
                this_args+=1;
                //cout << "this_args = " << this_args << endl;
                for(int i = 1; i < (n+1); i++) this_operator[i] = 0;
                int* these_args = new int[n];
                char* this_arg = NULL;
                this_arg = strtok(this_args,",");
                int c = 0;
                if(this_arg) {
                    these_args[0] = this_arg[0] - 'a' + 1;
                    c++;
                }
                while((this_arg = strtok(NULL,","))&&(c<=n)) {
                    these_args[c] = this_arg[0] - 'a' + 1;
                    c++;
                }
                for(int i = 0; i < c; i++) {
                    this_operator[i+1] = these_args[(i+c-1)%c];
                }
                AddOperator(this_operator);
                delete [] these_args;
            }
            my_file.getline(this_line,linewidth);
        }
        delete [] this_operator;

        delete [] this_line;
    }
    my_file.close();
}

SQC_Circuit* SQC_Circuit::LoadTFCFile(const char* inFilename) {
    SQC_Circuit* out = NULL;
    ifstream my_file(inFilename);
    if(my_file.good()) {
        int max_line = 1000;
        int max_qubits = 100;
        char this_line[max_line];
        // Search for 1st string .v <var 1> <var 2> ... <var n>
            // Assign each var string to a qubit index 1 to n
        int n = 0;
        string qubit_strings[max_qubits];
        string v_line = "";
        string i_line = "";
        string o_line = "";
        string c_line = "";
        int file_state = 0; // 0-before .v; 1-between .v and BEGIN; 2-between BEGIN and END; 3-After END before eof;
        while(!my_file.eof()) {
            my_file.getline(this_line,max_line);
            switch(file_state) {
                case 0: //before .v
                    {
                        if((strlen(this_line)>=2)&&(this_line[0]=='.')) {
                            switch(this_line[1]) {
                                case 'v':
                                    {
                                        v_line = this_line;
                                        if(v_line.length()) {
                                            char v_line_c[max_line];
                                            strcpy(v_line_c,v_line.c_str());
                                            char* this_tok = NULL;
                                            this_tok = strtok(v_line_c," ,\t");
                                            do {
                                                this_tok = strtok(NULL," ,\t");
                                                if(this_tok) {
                                                    qubit_strings[n] = this_tok;
                                                    n++;
                                                }
                                            } while(this_tok);
                                        }
                                        if(n) out = new SQC_Circuit(n,n);
                                        file_state = 1;
                                    }
                                    break;
                                case 'i':
                                    {

                                    }
                                    break;
                                case 'o':
                                    {

                                    }
                                    break;
                                case 'c':
                                    {

                                    }
                                    break;
                                default:
                                    {

                                    }
                                    break;
                            }
                        }
                    }
                    break;
                case 1: //between .v and BEGIN
                    {
                        if(!strcmp(this_line,"BEGIN")) {
                            file_state = 2;
                        }
                    }
                    break;
                case 2: // between BEGIN and END
                    {
                        if(!strcmp(this_line,"END")) {
                            file_state = 3;
                        } else {
                            char* this_tok = NULL;
                            this_tok = strtok(this_line," ,\t");
                            if(strlen(this_tok)) {
                                switch(this_tok[0]) {
                                    case 't':
                                        {
                                            this_tok++;
                                            int toff_n = 0;
                                            toff_n = atoi(this_tok);
                                            if(toff_n>0) {
                                                int this_gate[n+1];
                                                this_gate[0] = SQC_OPERATOR_TOFFOLI_N;
                                                for(int i = 0; i < n; i++) this_gate[i+1] = 0;
                                                int q_count = 0;
                                                while(this_tok=strtok(NULL," ,\t")) {
                                                    if(strlen(this_tok)) {
                                                        int this_q = 0;
                                                        for(int i = 0; (this_q==0)&&(i < n); i++) {
                                                            if(!qubit_strings[i].compare(this_tok)) {
                                                                this_q = (i+1);
                                                            }
                                                        }
                                                        if(this_q>0) {
                                                            this_gate[q_count+1] = this_q;
                                                            q_count++;
                                                        }
                                                    }
                                                }
                                                if(q_count==toff_n)
                                                    out->AddOperator(this_gate);
                                            }
                                        }
                                        break;
                                    case 'H':
                                        {
                                            int this_gate[n+1];
                                            for(int i = 0; i < (n+1); i++) this_gate[i] = 0;
                                            this_gate[0] = SQC_OPERATOR_HADAMARD;
                                            int this_q = 0;
                                            this_tok = strtok(NULL," ,\t");
                                            for(int i = 0; (this_q==0)&&(i < n); i++) {
                                                if(!qubit_strings[i].compare(this_tok)) {
                                                    this_q = (i+1);
                                                }
                                            }
                                            if(this_q!=0) {
                                                this_gate[1] = this_q;
                                                out->AddOperator(this_gate);
                                            }
                                        }
                                        break;
                                    case 'T':
                                        {
                                            int this_gate[n+1];
                                            for(int i = 0; i < (n+1); i++) this_gate[i] = 0;
                                            if((strlen(this_tok)>1)&&(this_tok[1]=='*')) {
                                                this_gate[0] = SQC_OPERATOR_T_DAG;
                                            } else {
                                                this_gate[0] = SQC_OPERATOR_T;
                                            }
                                            int this_q = 0;
                                            this_tok = strtok(NULL," ,\t");
                                            for(int i = 0; (this_q==0)&&(i < n); i++) {
                                                if(!qubit_strings[i].compare(this_tok)) {
                                                    this_q = (i+1);
                                                }
                                            }
                                            if(this_q!=0) {
                                                this_gate[1] = this_q;
                                                out->AddOperator(this_gate);
                                            }
                                        }
                                        break;
                                    case 'S':
                                        {
                                            int this_gate[n+1];
                                            for(int i = 0; i < (n+1); i++) this_gate[i] = 0;
                                            if((strlen(this_tok)>1)&&(this_tok[1]=='*')) {
                                                this_gate[0] = SQC_OPERATOR_S_DAG;
                                            } else {
                                                this_gate[0] = SQC_OPERATOR_S;
                                            }
                                            int this_q = 0;
                                            this_tok = strtok(NULL," ,\t");
                                            for(int i = 0; (this_q==0)&&(i < n); i++) {
                                                if(!qubit_strings[i].compare(this_tok)) {
                                                    this_q = (i+1);
                                                }
                                            }
                                            if(this_q!=0) {
                                                this_gate[1] = this_q;
                                                out->AddOperator(this_gate);
                                            }
                                        }
                                        break;
                                    case 'Z':
                                        {
                                            int this_gate[n+1];
                                            for(int i = 0; i < (n+1); i++) this_gate[i] = 0;
                                            this_gate[0] = SQC_OPERATOR_Z;
                                            int this_q = 0;
                                            this_tok = strtok(NULL," ,\t");
                                            for(int i = 0; (this_q==0)&&(i < n); i++) {
                                                if(!qubit_strings[i].compare(this_tok)) {
                                                    this_q = (i+1);
                                                }
                                            }
                                            if(this_q!=0) {
                                                this_gate[1] = this_q;
                                                out->AddOperator(this_gate);
                                            }
                                        }
                                        break;
                                }
                            }
                        }
                    }
                    break;
                case 3: // After END before eof;
                    {

                    }
                    break;
            }
        }

        // Search for .i line
        // Search for .o line
        // Search for .c line
        // Search for BEGIN line

    }
    my_file.close();
    return out;
}

/*
void SQC_Circuit::Print() const {
    LOut() << "Circuit:" << endl;
    //LOut_Pad++;
    int start_i = 0;
    int print_n = -1;
    bool with_t = true;
    LOut() << "n " << n << endl;
    LOut() << "m " << max_m << endl;
    LOut() << "p " << p << endl;
    if(print_n<0) print_n = m; else print_n = (start_i+print_n);
    for(int i = start_i; i < fmin(m,print_n); i++) {
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
        case SQC_OPERATOR_PARTITION:
            this_op_str = SQC_OPSTRING_PARTITION;
            break;
        case SQC_OPERATOR_X:
            this_op_str = SQC_OPSTRING_X;
            break;
        case SQC_OPERATOR_Y:
            this_op_str = SQC_OPSTRING_Y;
            break;
        }
        if(with_t) {
            LOut() << "[" << i << "] = ";
            cout << this_op_str;
        } else LOut() << this_op_str;
        int j = 0;
        int this_qubit=0;
        while((this_qubit = operator_list[i][j+1])&&(j<n)) {
            cout << " " << this_qubit;
            j++;
        }
        LOut() << endl;
        //LOut_Pad--;
    }
}
*/

void SQC_Circuit::ConvertHadamard(SQC_Circuit* out, int max_ancillas) {
    //cout << "Go" << endl;
    //cout << "max_ancillas = " << hadamard_mode_max_ancillas << endl;
    // Assumes circuit is composed only of {T, CS, CCZ, Clifford}
    // First determine the number of hadamards in the circuit
    if(out) {
        int n_hadamards = CountOperators(SQC_OPERATOR_HADAMARD);
        if(out->operator_list) {
            out->Destruct();
        }
        if(max_ancillas==-1) {
            max_ancillas = n_hadamards;
        }
        if(hadamard_mode_max_ancillas==-1) {
            hadamard_mode_max_ancillas = n_hadamards;
        }
        out->n = n+max_ancillas;
        out->p = p+max_ancillas;

        out->Construct();
        int hadamard_count = 0;
        int i = 0;
        bool exit = 0;
        int* this_operator = new int[out->n+1];

        while(!exit) {
            for(int c = 0; c < out->n+1; c++) this_operator[c] = 0;
            if(operator_list[0][0] == SQC_OPERATOR_HADAMARD) {
                if(hadamard_count>=max_ancillas) {
                    exit = 1;
                } else {
                    //cout << "Adding hadamard" << endl;
                    /*this_operator[0] = SQC_OPERATOR_CS;
                    this_operator[1] = operator_list[0][1];
                    this_operator[2] = (n+1+hadamard_count);
                    out->AddOperator(this_operator);*/
                    this_operator[0] = SQC_OPERATOR_CNOT;
                    this_operator[1] = operator_list[0][1];
                    this_operator[2] = (n+1+hadamard_count);
                    out->AddOperator(this_operator);
                    this_operator[0] = SQC_OPERATOR_CNOT;
                    this_operator[1] = (n+1+hadamard_count);
                    this_operator[2] = operator_list[0][1];
                    out->AddOperator(this_operator);
                    this_operator[0] = SQC_OPERATOR_CNOT;
                    this_operator[1] = operator_list[0][1];
                    this_operator[2] = (n+1+hadamard_count);
                    out->AddOperator(this_operator);
                    DeleteOperator(0);
                    hadamard_count++;
                }
            } else {
                //cout << "Adding other" << endl;
                for(int c = 0; c < (n+1); c++) {
                    //cout << operator_list[0][c];
                    this_operator[c] = operator_list[0][c];
                }
                //cout << endl;
                out->AddOperator(this_operator);
                DeleteOperator(0);
            }

            if(m==0) {
                exit = 1;
            }
        }
        delete [] this_operator;
    }
    //cout << "Stop" << endl;
}

int SQC_Circuit::CountOperators(SQC_Operator_Label in_op) const {
    int out = 0;
    if(in_op!=SQC_OPERATOR_N) {
        for(int i = 0; i < m; i++) {
            if(operator_list[i][0]==in_op) {
                out++;
            }
        }
    } else {
        out = m;
    }
    return out;
}

Matrix SQC_Circuit::toMatrix() const {
    Matrix out = Matrix::identity(pow(2,n));

    for(int i = 0; i < m; i++) {
        Matrix this_mat = Matrix::identity(n);
        switch(operator_list[i][0]) {
            case SQC_OPERATOR_IDENTITY:
                break;
            case SQC_OPERATOR_HADAMARD:
                {
                this_mat = Matrix::H(n,operator_list[i][1]);
                }
                break;
            case SQC_OPERATOR_CNOT:
                {
                this_mat = Matrix::CNOT(n,operator_list[i][2],operator_list[i][1]);
                }
                break;
            case SQC_OPERATOR_T:
                {
                this_mat = Matrix::R(n,operator_list[i][1],L_PI/4.0);
                }
                break;
            case SQC_OPERATOR_CS:
                {
                this_mat = Matrix::CS(n,operator_list[i][1],operator_list[i][2]);
                }
                break;
            case SQC_OPERATOR_CCZ:
                {
                this_mat = Matrix::CCZ(n,operator_list[i][1],operator_list[i][2],operator_list[i][3]);
                }
                break;
            case SQC_OPERATOR_S:
                {
                this_mat = Matrix::R(n,operator_list[i][1],L_PI/2.0);
                }
                break;
            case SQC_OPERATOR_Z:
                {
                this_mat = Matrix::Z(n,operator_list[i][1]);
                }
                break;
            case SQC_OPERATOR_X:
                {
                this_mat = Matrix::X(n,operator_list[i][1]);
                }
                break;
            case SQC_OPERATOR_Y:
                {
                this_mat = Matrix::Y(n,operator_list[i][1]);
                }
                break;
            case SQC_OPERATOR_TOFFOLI:
            case SQC_OPERATOR_TOFFOLI_N:
                {
                    int n_args = GetNArgs(i);
                    int* these_args = new int[n_args];
                    for(int c = 0; c < n_args; c++) these_args[c] = operator_list[i][c+1];
                    this_mat = Matrix::CNOT(n,these_args,n_args);
                    delete [] these_args;
                }
                break;
            default:
                break;
        }
        out = (this_mat*out);
    }

    return out;
}

void SQC_Circuit::RemoveExternalHadamards() {
    bool* operator_used = new bool[n];
    for(int i = 0; i < n; i++) operator_used[i] = 0;

    for(int t = 0; t < m; t++) {
        if(operator_list[t][0]==SQC_OPERATOR_HADAMARD) {
            if(!operator_used[operator_list[t][1]-1]) {
                DeleteOperator(t);
                t--;
            }
        } else {
            int this_n_args = GetNArgs(t);
            for(int i = 0; i < this_n_args; i++) operator_used[operator_list[t][i+1]-1] = 1;
        }
    }

    for(int i = 0; i < n; i++) operator_used[i] = 0;

    for(int t = (m-1); t >= 0; t--) {
        if(operator_list[t][0]==SQC_OPERATOR_HADAMARD) {
            if(!operator_used[operator_list[t][1]-1]) {
                DeleteOperator(t);
            }
        } else {
            int this_n_args = GetNArgs(t);
            for(int i = 0; i < this_n_args; i++) operator_used[operator_list[t][i+1]-1] = 1;
        }
    }

    delete [] operator_used;
}

bool SQC_Circuit::CancelAdjacentHadamards() {
    bool out = 0;
    int n_deleted = 0;
    for(int i = 0; i < (m-1); i++) {
        //cout << "i = " << i << endl;
        if(operator_list[i][0]==SQC_OPERATOR_HADAMARD) {
            int acting_qubit = operator_list[i][1];
            bool this_qubit_used = 0;
            for(int j = (i+1); (!this_qubit_used)&&(j < m); j++) {
                switch(operator_list[j][0]) {
                    case SQC_OPERATOR_HADAMARD:
                        if((!this_qubit_used)&&(operator_list[j][1]==acting_qubit)) {
                            //cout << "Deleting operator " << j << endl;
                            //for(int c = 0; c < (GetNArgs(j)+1); c++) cout << operator_list[j][c] << " ";
                            //cout << endl;
                            //cout << "Deleting operator " << i << endl;
                            //for(int c = 0; c < (GetNArgs(i)+1); c++) cout << operator_list[i][c] << " ";
                            //cout << endl;
                            DeleteOperator(j);
                            DeleteOperator(i);
                            i--;
                            n_deleted++;
                            n_deleted++;
                            this_qubit_used = 1;
                            out = 1;
                        }
                        break;
                    case SQC_OPERATOR_S:
                    case SQC_OPERATOR_Z:
                    case SQC_OPERATOR_T:
                        if(operator_list[j][1]==acting_qubit) {
                            this_qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CS:
                    case SQC_OPERATOR_CNOT:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)) {
                            this_qubit_used = 1;
                        }
                        break;
                    case SQC_OPERATOR_CCZ:
                        if((operator_list[j][1]==acting_qubit)||(operator_list[j][2]==acting_qubit)||(operator_list[j][3]==acting_qubit)) {
                            this_qubit_used = 1;
                        }
                        break;
                }
            }
        }
    }
    return out;
}

bool SQC_Circuit::CancelAdjacentTs() {
    bool out = 0;
    //cout << "FDSAGHAF" << endl;
    //bool* qubit_used = new bool[n];
    //for(int i = 0; i < n; i++) qubit_used[i] = 0;
    for(int i = 0; i < (m-1); i++) {
        //cout << "i = " << i << endl;
        switch(operator_list[i][0]) {
            case SQC_OPERATOR_T:
            case SQC_OPERATOR_CS:
            case SQC_OPERATOR_CCZ:
                {
                    //cout << "jmnhgtyr" << endl;
                    int this_nargs = GetNArgs(i);
                    int acting_qubit[3];
                    acting_qubit[0] = operator_list[i][1];
                    acting_qubit[1] = operator_list[i][2];
                    acting_qubit[2] = operator_list[i][3];
                    if(acting_qubit[2]>acting_qubit[1]) swap(acting_qubit[2],acting_qubit[1]);
                    if(acting_qubit[1]>acting_qubit[0]) swap(acting_qubit[1],acting_qubit[0]);
                    if(acting_qubit[2]>acting_qubit[1]) swap(acting_qubit[2],acting_qubit[1]);
                    //cout << "\tActing qubit = (" << acting_qubit[0] << ", " << acting_qubit[1] << ", " << acting_qubit[2] << ")" << endl;
                    bool this_qubit_used = 0;
                    for(int j = (i+1); (!this_qubit_used)&&(j < m); j++) {
                        //cout << "\t\tj = " << j << endl;
                        int jacting_qubit[3];
                        jacting_qubit[0] = operator_list[j][1];
                        jacting_qubit[1] = operator_list[j][2];
                        jacting_qubit[2] = operator_list[j][3];
                        if(jacting_qubit[2]>jacting_qubit[1]) swap(jacting_qubit[2],jacting_qubit[1]);
                        if(jacting_qubit[1]>jacting_qubit[0]) swap(jacting_qubit[1],jacting_qubit[0]);
                        if(jacting_qubit[2]>jacting_qubit[1]) swap(jacting_qubit[2],jacting_qubit[1]);
                        //cout << "\t\t\tjActing qubit = (" << jacting_qubit[0] << ", " << jacting_qubit[1] << ", " << jacting_qubit[2] << ")" << endl;
                        switch(operator_list[j][0]) {
                            case SQC_OPERATOR_T:
                            case SQC_OPERATOR_CS:
                            case SQC_OPERATOR_CCZ:
                                if((!this_qubit_used)&&(jacting_qubit[0]==acting_qubit[0])&&(jacting_qubit[1]==acting_qubit[1])&&(jacting_qubit[2]==acting_qubit[2])) {
                                    //cout << "Deleting operator " << j << endl;
                                    //for(int c = 0; c < (GetNArgs(j)+1); c++) cout << operator_list[j][c] << " ";
                                    //cout << endl;
                                    //cout << "Deleting operator " << i << endl;
                                    //for(int c = 0; c < (GetNArgs(i)+1); c++) cout << operator_list[i][c] << " ";
                                    //cout << endl;
                                    DeleteOperator(j);
                                    DeleteOperator(i);
                                    i--;
                                    //n_deleted++;
                                    //n_deleted++;
                                    this_qubit_used = 1;
                                    out = 1;
                                }
                                break;
                            case SQC_OPERATOR_HADAMARD:
                                for(int c = 0; c < this_nargs; c++) {
                                    if(operator_list[j][1]==acting_qubit[c]) {
                                        this_qubit_used = 1;
                                    }
                                }
                                break;
                            case SQC_OPERATOR_S:
                                break;
                            case SQC_OPERATOR_Z:
                                break;
                            case SQC_OPERATOR_CNOT:
                                for(int c = 0; c < this_nargs; c++) {
                                    if(operator_list[j][1]==acting_qubit[c]) {
                                        this_qubit_used = 1;
                                    }
                                }
                                break;
                        }
                    }
                }
                break;


        }
    }
    //delete [] qubit_used;
    return out;
}

void SQC_Circuit::AddAncilla(int anc_inc) {
    ;
}

void SQC_Circuit::PrintOperatorDistribution() const {
    int hist[SQC_OPERATOR_N];
    for(int i = 0; i < SQC_OPERATOR_N; i++) {
        hist[i] = 0;
    }
    for(int t = 0; t < m; t++) {
        hist[operator_list[t][0]]++;
    }
    cout << "Operator distribution:" << endl;
    cout << "\tid\tCount" << endl;
    for(int i = 0; i < SQC_OPERATOR_N; i++) {
        cout << "\t" << i << "\t" << hist[i] << endl;
    }
}

double VerifyOptimization(const SQC_Circuit& U1, const SQC_Circuit& U2) {
    double out = 0.0;

    int n = max(U1.n, U2.n);

    int N = pow(2,n);

    Matrix matrix_in = TO_Maps::SQC_Circuit_to_Matrix(U1,n);
    cout << "U1:" << endl;
    matrix_in.print();
    Matrix matrix_out = TO_Maps::SQC_Circuit_to_Matrix(U2,n);
    cout << "U2:" << endl;
    matrix_out.print();

    Matrix U2d = matrix_out.adjoint();
    Matrix U1U2d = matrix_in*U2d;

    Complex U1U2d_trace;
    for(int i = 0; i < N; i++) {
        U1U2d_trace += U1U2d.E(i,i);
    }

    Matrix temp = (U1U2d*N)*(Complex(1.0,0.0)/U1U2d_trace);
    Matrix temp2 = Matrix::identity(N);
    Matrix matrix_brackets = (temp2 - temp);

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            out -= matrix_brackets.E(i,j).magnitude();
        }
    }


    return out;
}

double VerifyOptimization2(const SQC_Circuit& U1, const SQC_Circuit& U2) {
    LOut() << "Begin circuit verification..." << endl;
    LOut_Pad++;
    LOut() << "^ = tensor product" << endl;
    LOut() << "` = Hermitian conjugate" << endl;
    LOut() << "F = fidelity" << endl;
    double out = 0.0;
    int N_d = 1;

    int n = U1.d;
    int N = pow(2,n);
    Matrix Uin = TO_Maps::SQC_Circuit_to_Matrix(U1);
    LOut() << "U_in:" << endl;
    Uin.print(); LOut() << endl;
    Matrix Uout = TO_Maps::SQC_Circuit_to_Matrix(U2);
    LOut() << "U_out:" << endl;
    Uout.print(); LOut() << endl;
    //Uout = (Matrix::M(2,2)*Uout*Matrix::M(2,2));
    //LOut() << "Uout:" << endl;
    //Uout.print();
    Matrix psi(N,1);
    srand(time(NULL));
    for(int d = 0; d < N_d; d++) {
        for(int i = 0; i < N; i++) psi.E(i,0,Utils::rand_d(),Utils::rand_d());
        psi = (psi*sqrt(1.0/((psi.adjoint()*psi).E(0,0).getReal())));
        //LOut() << "norm(psi) = " << ((psi.adjoint()*psi).E(0,0).getReal()) << endl;
        LOut() << "|psi> = random state" << endl;
        psi.print(); LOut() << endl;
        int n_anc = U2.n-n;
        int N_anc = pow(2,n_anc);
        Matrix anc(N_anc,1);
        anc.E(0,0,1.0);
        LOut() << "|anc> = hadamard ancillas" << endl;
        anc.print(); LOut() << endl;
        Matrix phi = (anc^psi);
        LOut() << "|phi> = |anc> ^ |psi>" << endl;
        phi.print(); LOut() << endl;

        phi = (Uout*phi);
        LOut() << "U_out |phi> = " << endl;
        phi.print(); LOut() << endl;

        Matrix rho_out = phi*phi.adjoint();
        LOut() << "rho_out = U_out |phi><phi| U_out'" << endl;
        rho_out.print(); LOut() << endl;

        Matrix rho = rho_out;
        for(int i = U2.n; i > (U2.n-U2.p_hads); i--) {
            LOut() << "Tracing qubit " << i << " of rho_out" << endl;
            rho = rho.partialTrace(i);
        }
        /*for(int i = U2.n; i > (U2.n-U2.p_hads); i--) {
            LOut() << "Tracing qubit " << i << " of rho_out" << endl;
            rho = rho.partialTrace(i);
        }*/

        rho = rho*(1.0/rho.trace().magnitude());
        LOut() << "rho_out` = Tr_{anc}(rho_out)" << endl;
        rho.print(); LOut() << endl;

        psi = (Uin*psi);
        LOut() << "U_in|psi>" << endl;
        psi.print(); LOut() << endl;

        LOut() << "rho_in = U_in |psi><psi| U_in'" << endl;
        (psi*psi.adjoint()).print(); LOut() << endl;

        Complex expt = (psi.adjoint()*rho*psi).E(0,0);
        LOut() << "<psi|rho_out`|psi> =" << endl;
        expt.print(); LOut() << endl;

        out += sqrt(expt.magnitude());
        cout << "F(U_in|psi>, rho_out`) = " << sqrt(expt.magnitude()) << endl;
    }
    out /= N_d;

    LOut_Pad--;

    LOut() << "End circuit verification." << endl;


    return out;
}
