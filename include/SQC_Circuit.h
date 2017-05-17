#ifndef HEADER_SQC_CIRCUIT
#define HEADER_SQC_CIRCUIT

#include <iostream>
using namespace std;
#include <ostream>
#include "Signature.h"
#include "BMSparse.h"

typedef int* SQC_Operator; // 0^th element = Operator type, rest of elements = qubit labels
typedef SQC_Operator* SQC_Operator_List;

const int SQC_DEFAULT_MAX_M = 1;

enum SQC_AncillaMode {
    SQC_ANCILLA_MODE_MANUAL,
    SQC_ANCILLA_MODE_PER_GATE,
    SQC_ANCILLA_MODE_PER_CIRCUIT
};

enum SQC_ToffoliNMode {
    SQC_TOFFOLI_N_MODE_TOFF3, // Concatenate Toffoli-3 gates
    SQC_TOFFOLI_N_MODE_JONES   // Add controls using Cody Jones method
};

enum SQC_Operator_Label {
        SQC_OPERATOR_IDENTITY,
        SQC_OPERATOR_HADAMARD,
        SQC_OPERATOR_CNOT,
        SQC_OPERATOR_T,
        SQC_OPERATOR_CS,
        SQC_OPERATOR_CCZ,
        SQC_OPERATOR_S,
        SQC_OPERATOR_Z,
        SQC_OPERATOR_TOFFOLI,
        SQC_OPERATOR_TOFFOLI_4,
        SQC_OPERATOR_TOFFOLI_N

// To add operator, update:
//  Print, AddOperator, GetPartition
};
constexpr static char SQC_OPSTRING_IDENTITY[] = "I";
constexpr static char SQC_OPSTRING_HADAMARD[] = "H";
constexpr static char SQC_OPSTRING_CNOT[] = "CNOT";
constexpr static char SQC_OPSTRING_T[] = "T";
constexpr static char SQC_OPSTRING_CS[] = "CS";
constexpr static char SQC_OPSTRING_CCZ[] = "CCZ";
constexpr static char SQC_OPSTRING_S[] = "S";
constexpr static char SQC_OPSTRING_Z[] = "Z";
constexpr static char SQC_OPSTRING_TOFFOLI[] = "Toffoli";
constexpr static char SQC_OPSTRING_TOFFOLI_4[] = "Toffoli-4";
constexpr static char SQC_OPSTRING_TOFFOLI_N[] = "t"; // N is given by number of arguments using formula N = (nargs+3)/2


struct SQC_Circuit {
    // Properties
    int n = 0; //Number of qubits (total)
    int p = 0; //Number of ancillas, ancilla indices: n-p+1 -> n
               //Number of non-ancilla qubits = n-p
               //If p=-1, then allocate ancillas automatically using SQC_ANCILLA_MODE_PER_GATE
               //If p=-2, then allocate ancillas automatically using SQC_ANCILLA_MODE_PER_CIRCUIT
    int m = 0; //Number of operators
    int max_m = SQC_DEFAULT_MAX_M; //Max number of operators
    SQC_Operator_List operator_list = NULL;
    SQC_AncillaMode ancilla_mode = SQC_ANCILLA_MODE_MANUAL;
    SQC_ToffoliNMode toffoli_n_mode = SQC_TOFFOLI_N_MODE_JONES;

    // Methods

    SQC_Circuit();
    void Construct();
    void Destruct();
    void Copy(const SQC_Circuit& in_C);

    void Print(ostream* in_OS = &cout) const;
    void Load(const char* in_filename);
    void LoadMaslovFile(const char* in_filename);
    void Save(const char* in_filename) const;
    void Clear();
    void Resize(int in_max_m);

    void AddOperator(const char* in_op_str);
    void AddOperator(const SQC_Operator in_op);
    void DeleteOperator(int t);

    bool GetPartition(SQC_Circuit* out_Hadamards, SQC_Circuit* out_CNOT_T);
    void DecompositionVW(SQC_Circuit* out_V, SQC_Circuit* out_W) const;
    BMSparse toGateSynthesisMatrix() const;
    bool NextSignature(Signature& outSig);

    void ReplaceOperator(SQC_Circuit* in_new_ops, int t, int n_rep=1);
    void ConvertFromToffoli();
    void AllocateAncillas(const SQC_Circuit& in_C); // Constructs new circuit with automatically allocated ancillas
    int GetNArgs(int i) const;
};

#endif // HEADER_SQC_CIRCUIT
