#ifndef TO_DECODER_HEADER
#define TO_DECODER_HEADER

#include <iostream>
using namespace std;

#include "GateStringSparse.h"
#include "Signature.h"
#include "GateSigInterface.h"
#include "LukeBool.h"
#include "PhasePolynomial.h"
#include "Matrix.h"
#include "BMSparse.h"
#include "LukeInt.h"
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <string>
#include <sstream>

const double update_time = 5.0;

// Global variables
extern string g_output_filename;
extern string g_csv_filename;
extern string g_indvar_out;
extern ostringstream g_h_order_stream;
extern string g_best_random_h_order_f, g_best_random_h_order_nf;
extern bool g_error_report;
extern string g_algorithm;
extern int g_random_circuit_seed;
extern bool g_lempel_feedback;
extern int g_Reed_Muller_max;
extern int g_Hadamard_ancillas;
extern double g_matrix_precision;

namespace SYNTHESIS_ALGORITHM_TAG {
    const string DAFT_GUESS = "d";
    const string LEMPEL_LEAST_GREEDY = "ll";
    const string LEMPEL_GREEDY = "lg";
    const string LEMPEL_RANDOM = "lr";
    const string LEMPEL_X = "lx";
    const string LEMPEL_X_2 = "lx2";
    const string REED_MULLER = "rm";
    const string ALL_LEMPEL_SELECTORS = "als";
    const string ALL_LEMPEL = "all";
}
enum SYNTHESIS_ALGORITHM_ID {
    DAFT_GUESS,
    REED_MULLER,
    LEMPEL_LEAST_GREEDY_F,
    LEMPEL_GREEDY_F,
    LEMPEL_RANDOM_F,
    LEMPEL_LEAST_GREEDY_NF,
    LEMPEL_GREEDY_NF,
    LEMPEL_RANDOM_NF,
    LEMPEL_X,
    LEMPEL_X_2
};

typedef GateStringSparse (*TO_Decoder)(const Signature& in_S);
typedef int (*LempelSelector) (const Signature& inS);

GateStringSparse ReedMullerSynthesis(const Signature& inS);
GateStringSparse ReedMullerSynthesis2(const Signature& inS);

//LempelX synthesis functions
GateStringSparse LempelXSynthesis(const Signature& inS);
GateStringSparse LempelXSynthesis2(const Signature& inS);



//Lempel selectors
int LempelSelector_LeastGreedy(const Signature& inS);
int LempelSelector_Greedy(const Signature& inS);
int LempelSelector_Random(const Signature& inS);

//Main Lempel synthesis functions
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, LempelSelector lempelSelector, bool feedback, int psrmc=0, int col_order=0, int n_orders=1);
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, bool feedback, int n_rand = 1, int psrmc=0, int col_order=0, int n_orders=1);
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, int n_rand = 1, int psrmc=0, int col_order=0, int n_orders=1);
GateStringSparse LempelSynthesis(const Signature& inS);

void result_analysis(const Signature& inS, const GateStringSparse& inResult, ostream& inOS = cout);
bool synthesis_success(const Signature& inS, const GateStringSparse& inResult);
bool T_fast(const GateStringSparse& inGSS);

PhasePolynomial FullDecoderWrapper(const PhasePolynomial& in, TO_Decoder decoder = LempelXSynthesis);

Matrix RemoveRoundingError(const Matrix& in);

//Post-synthesis Reed-Muller corrections
GateStringSparse PSRMC_Lightweight_Search(const GateStringSparse& inGSS);
GateStringSparse PSRMC_Complex_Submat(const GateStringSparse& inGSS, int n_RM=6, int col_order = 0, int n_orders = 1);
int trapezoid_width(const BMSparse& inBMS, int height);

#endif // TO_DECODER_HEADER
