#ifndef UNIVERSAL_OPTIMIZE_HEADER

#include "SQC_Circuit.h"

// 0. Full UniversalOptimize algorithm
void UniversalOptimize();

// 1. Convert Y to XZ
void convert_Y_to_XZ();

// 2. Expand Toff_n to Toff_3
int N_Toffn_ancillas() const;
void expand_Toffn_to_Toff3(SQC_Circuit& out) const; // Assumes out already has enough ancillas

// 3. Convert X and Toff_3 to H, Z and CCZ
void convert_X_to_HZH();
void convert_Toff3_to_HCCZH();

// 4. Convert CS and CCZ to CNOT and T
void convert_CS_to_CNOT_T();
void convert_CCZ_to_CNOT_T();

// 5. Strip external Hadamards / 6. Partition internal circuit such that each partition contains no more than n_H Hadamards
void decompose_into_Hadamard_partitions(SQC_Circuit* inHs, SQC_Circuit* inPs); // U = H_NP P_NP ... H_1 P_1 H_0
int N_max_Hadamards() const;

// 7. Convert internal Hadamards of each partition to external Hadamards and postselection
void Hadamards_to_Gadgets(SQC_Circuit& H_1, SQC_Circuit& Pp, SQC_Circuit& H_2, SQC_Circuit& E); // P = E H_2 Pp H_1; assumes outputs all have correct number of ancillas

// 8. Split each partition into CNOT circuit * Diagonal C_3 circuit
void decompose_C3_to_CNOT_D3(SQC_Circuit& CNOT, SQC_Circuit& D3);

// 9. Optimize each partition's diagonal circuit w.r.t. T-count and map back to circuit decomposition.
void optimize_D3(GateStringSparse (*in_Decoder)(Signature& in_S));

#endif // UNIVERSAL_OPTIMIZE_HEADER
