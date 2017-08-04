// TODO: Move all TOptimizeRM class1 <-> class2 maps here. E.g. GateStringSparse <-> BMSparse
//      - Use format: class2 Class1_to_Class2(const class1& in)

#ifndef TO_MAPS_HEADER
#define TO_MAPS_HEADER

#include "PhasePolynomial.h"
#include "SQC_Circuit.h"
#include "GateStringSparse.h"

namespace TO_Maps {
    // PhasePolynomial <-> SQC_Circuit
    PhasePolynomial SQC_Circuit_to_PhasePolynomial(const SQC_Circuit& in);
    SQC_Circuit PhasePolynomial_to_SQC_Circuit(const PhasePolynomial& in);

    // PhasePolynomial <-> GateStringSparse
    GateStringSparse PhasePolynomial_to_GateStringSparse(const PhasePolynomial& in);
    PhasePolynomial GateStringSparse_to_PhasePolynomial(const GateStringSparse& in);

}

#endif // TO_MAPS_HEADER
