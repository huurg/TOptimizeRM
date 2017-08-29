#include <iostream>
using namespace std;

#include "TO_Decoder.h"
#include "Signature.h"
#include "GateStringSparse.h"
#include "GateSigInterface.h"
#include "LukeMaths.h"
#include "LukeBool.h"
#include "BoolMat.h"
#include "BMSparse.h"
#include "Interface_SigBMS.h"
#include "Interface_BMSGSS.h"
#include "LukeInt.h"
#include "Bool_Signature.h"
#include "LukeMat_GF2.h"
#include "GateSynthesisMatrix.h"
#include "SQC_Circuit.h"
#include "LukeConsoleOut.h"
#include "PhasePolynomial.h"
#include "TO_Maps.h"
#include "WeightedPolynomial.h"

using namespace LukeConsoleOut;
//using namespace LinBox;
//using namespace std;

#include <climits>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <utility>

// TODO (Refactoring): Consolidate all 'interface' namespaces into a single 'converter' namespace.



GateStringSparse findNearOptimal(const Signature& inS, int maxRM=5);

void NullSpaceTest();
void Sig2BMSTest();
void compositionTest();
void hdecomtest();
void debuggingTest();
void n3benchmark(int n, int maxRM=5, int randgen = 0);
void n3benchmark_d(int n, int d, int maxRM=5, int randgen = 0);

//Universal T-Count finders
int UniversalTCount(SQC_Circuit* inC, int* out_daft = NULL, double* out_exec_time = NULL, int* n_parts = NULL, int* n_hadancs = NULL);

//Structured circuit generators
Signature CircuitGenerator(const string& inS);
Signature CircuitGenerator_Toffhash(int N_hash);
Signature CircuitGenerator_Toffoli(int N_toff);
Signature CircuitGenerator_RandomComplex(int n, int in_seed = 0);


//BMSparse decompTest(const BMSparse& inB);

//Tests
void test_PSRMC_1();
void test_Nullspace_t_vs_n();
void test_LX_exec_times();


//int main(int argc, char* argv[]);



void NullSpaceTest() {
    //BMSparse blah = BMSparse::eye(6,6);
    BMSparse blah(6,6);

    BMSparse thisL, thisQ, thisU, thisP;


    blah.E(2,1,1);
    //blah.E(1,0,1);
    blah.E(0,2,1);
    blah.E(0,1,1);
    //blah.E(1,3,1);
    blah.E(2,3,1);
    //blah.E(1,1,1);
    blah.E(3,0,1);
    blah.E(3,1,1);
    blah.E(5,4,1);
    blah.E(5,5,1);

    cout << "Blah" << endl;
    blah.printFull();

    blah.LQUPDecomposition(thisL, thisQ, thisU, thisP);

    cout << "L" << endl;
    thisL.printFull();
    cout << "Q" << endl;
    thisQ.printFull();
    cout << "U" << endl;
    thisU.printFull();
    cout << "P" << endl;
    thisP.printFull();

    BMSparse result = blah.nullspace();
    cout << "nullspace" << endl;
    result.printFull();

    cout << "nullspace is empty? " << result.isEmpty() << endl;

    // TODO: Carry on with RowEchelon
}

void Sig2BMSTest() {
    Signature blah = Signature::SigFromFile("circuits/test3.lsf");
    blah.print();

    BMSparse blahBMS = Interface_SigBMS::SigToBMS(blah);
    cout << "A" << endl;
    blahBMS.printFull();

    BMSparse thisB = blahBMS.elementaryFactorFast();
    cout << "B" << endl;
    thisB.printFull();

    cout << "B'" << endl;
    (thisB.T()).printFull();

    BMSparse result = thisB*thisB.T();
    cout << "A = B*B'" << endl;
    result.printFull();

    BMSparse nullspace_B = thisB.nullspace();
    cout << "nullspace" << endl;
    nullspace_B.printFull();

    BMSparse minfac = blahBMS.minimalFactor();
    cout << "Minimal factor" << endl;
    minfac.printFull();

    GateStringSparse thisGSS = Interface_BMSGSS::BMSToGSS(minfac);
    thisGSS.print();

    /*
    BMSparse blee = BMSparse::eye(6,6);
    blee.E(1,0,1);
    blee.E(2,0,1);
    blee.E(2,1,1);
    blee.E(3,1,1);
    blee.E(5,1,1);
    blee.E(4,3,1);
    blee.E(5,3,1);
    blee.E(5,4,1);
    blee.E(4,5,1);
    blee.E(5,5,0);
    cout << "Factor" << endl;
    blee.printFull();
    BMSparse result2 = blee*blee.T();
    cout << "A = B*B'" << endl;
    result2.printFull();
    */
}

void compositionTest() {
    BMSparse blah = BMSparse::eye(4,4);
    BMSparse blee(1,4);
    BMSparse blan(4,1);
    BMSparse bloo = BMSparse::eye(2,2);

    BMSparse composedLR = (blan&&blah)||blah.sub(2,0,1,-1);
    BMSparse composedUD = (blee||blah)&&composedLR;
    composedUD.sub(1,3,bloo);

    composedLR.printFull();
    composedUD.printFull();

    BMSparse emptyCol(2,0);
    BMSparse emptyRow(0,2);
    BMSparse result = emptyCol*emptyRow;
    result.printFull();
}

void hdecomtest() {
    Signature blah = Signature::SigFromFile("circuits/test3.lsf");
    GateStringSparse blah_GSS = GateSigInterface::SigToGSS(blah);
    cout << "Initial gate string" << endl;
    blah.print();
    Signature dimblah = blah.diminish(1);
    cout << "Diminished gate string" << endl;
    dimblah.print();
    Signature augblah = dimblah.augment(1);
    cout << "Augmented gate string" << endl;
    augblah.print();

    Signature this_alpha, this_A, this_B;

    blah.h_decomposition(1,this_alpha, this_A, this_B);

    BMSparse this_A_BMS = Interface_SigBMS::SigToBMS(this_A);

    BMSparse minfac = this_A_BMS.minimalFactor();

    GateStringSparse this_A_GSS = Interface_BMSGSS::BMSToGSS(minfac);

    GateStringSparse result = GateSigInterface::SigToGSS(this_alpha) + this_A_GSS.mult2xh(1) + GateSigInterface::SigToGSS(this_B);
    cout << "Improved  gate string" << endl;
    result.printString();
}

GateStringSparse findNearOptimal(const Signature& inS, int maxRM) {
    cout << endl;
    cout << "Circuit to optimize:" << endl;
    inS.print();
    cout << endl;

    GateStringSparse out(inS.get_n());

    Signature currentSig(inS.get_n());
    currentSig.assign(inS);
    bool exit = false;

    cout << "Optimization begin..." << endl << endl;
    clock_t start_time = clock();
    for(int h = 1; (!exit)&&(h <= (inS.get_n()-maxRM)); h++) {
        cout << "currentSig iteration h = " << h << endl;
        Signature this_alpha, this_A, this_B;
        currentSig.h_decomposition(1,this_alpha,this_A,this_B);
        BMSparse this_A_BMS = Interface_SigBMS::SigToBMS(this_A);
        BMSparse minfac = this_A_BMS.minimalFactor();
        GateStringSparse this_A_GSS = Interface_BMSGSS::BMSToGSS(minfac);
        GateStringSparse this_alpha_GSS = GateSigInterface::SigToGSS(this_alpha);
        for(int i = 1; i < h; i++) {
            GateStringSparse tempGSS = this_A_GSS.augment(1);
            GateStringSparse tempGSS2 = this_alpha_GSS.augment(1);
            this_A_GSS.assign(tempGSS);
            this_alpha_GSS.assign(tempGSS2);
        }

        {
            GateStringSparse tempout = out + this_alpha_GSS + this_A_GSS.mult2xh(h);
            out = tempout;
        }
        currentSig.assign(this_B);

        {
            Signature tempSig = currentSig.diminish(1);
            currentSig.assign(tempSig);
        }

        exit = currentSig.isEmpty();
    }
    if(!exit) {
        //cout << "currentSig before dim" << endl;
        //currentSig.print();
        /*
        for(int h = 1; h <= (inS.get_n()-maxRM); h++) {
            Signature tempsig = currentSig.diminish(1);
            currentSig.assign(tempsig);
        }
        */
        //cout << "currentSig after dim" << endl;
        //currentSig.print();
        GateStringSparse thisRM = ReedMullerSynthesis(currentSig);
        for(int h = 1; h <= (inS.get_n()-maxRM); h++) {
            GateStringSparse tempGSS = thisRM.augment(1);
            thisRM = tempGSS;
        }
        {
            GateStringSparse tempout = out + thisRM;
            out = tempout;
        }
    }
    cout << endl;
    clock_t total_duration = clock() - start_time;
    double total_duration_s = (double)total_duration/(double)CLOCKS_PER_SEC;
    cout << "...optimization done! Time taken = " << total_duration_s << "s" << endl << endl;

    cout << "Near optimal gate string:" << endl;
    out.print();
    cout << "Near optimal weight = " << out.weight(true) << endl;
    Signature result_sig = GateSigInterface::expandGSSTerm(out);
    cout << "Near optimal signature:" << endl;
    result_sig.print();
    bool successful = (inS + result_sig).isEmpty();
    GateStringSparse daft_string = GateSigInterface::SigToGSS(inS);
    cout << "Successful? " << (successful?"Yes":"No") << endl << endl;
    cout << "Daft gate string. Weight = " << daft_string.weight(true) << "." << endl;
    daft_string.print();
    int improvement = daft_string.weight(true) - out.weight(true);
    cout << "Weight reduced by " << improvement << endl << endl;
    GateStringSparse difference = daft_string + out;
    cout << "Difference:" << endl;
    difference.print();

    return out;
}



int main(int argc, char* argv[]) {
    g_output_filename.clear();
    g_csv_filename.clear();

    srand(time(NULL));

    // TODO: Have UniversalT/optimize produce SQC_Circuit decompositions
    //      DONE: Implement the rest of PhasePolynomial
    // TODO: Implement map from SQC_Circuit to unitary matrix.
    //      Modify UniversalT/optimize to check init matrix = final matrix

    /*SQC_Circuit* blah = SQC_Circuit::LoadTFCFile("test2.tfc");

    cout << "C_in:" << endl;
    blah->Print(); cout << endl;

    SQC_Circuit blee;
    blee.Load("test1.tfc");
    blee.d = 1;
    blee.p_hads = 1;

    cout << "C_out:" << endl;
    blee.Print(); cout << endl;

    Matrix U1 = TO_Maps::SQC_Circuit_to_Matrix(*blah);
    Matrix U2 = TO_Maps::SQC_Circuit_to_Matrix(blee);

    double result = VerifyOptimization2(*blah,blee);
    cout << "result = " << result << endl;*/

    /*Matrix this_mat = Matrix::CNOT(2,2,1);
    Matrix this_mat2 = Matrix::identity();
    int this_qargs[] = {1,2};
    Matrix this_mat3 = Matrix::CNOT(2,this_qargs,2);
    this_mat2 ^= Matrix::X();
    cout << "this_mat" <<endl;
    this_mat.print();
    cout << "this_mat2" <<endl;
    this_mat2.print();
    cout << "this_mat3" <<endl;
    this_mat3.print();
    Matrix this_q1(2,1);
    this_q1.E(1,0,1);
    Matrix this_q2(2,1);
    this_q2.E(0,0,1);
    cout << "this_q1" << endl;
    this_q1.print();
    cout << "this_q2" << endl;
    this_q2.print();
    Matrix this_state = this_q1;
    this_state^=this_q2;
    cout << "before" << endl;
    this_state.print();
    cout << "After" << endl;
    (this_mat*this_state).print();*/

    /*Matrix this_q1(2,1);
    this_q1.E(0,0,1.0/sqrt(2));
    this_q1.E(1,0,1.0/sqrt(2));
    Matrix this_q2(2,1);
    this_q2.E(0,0,1);
    cout << "q1" << endl;
    this_q1.print();
    cout << "q2" << endl;
    this_q2.print();
    cout << "psi = q2 ^ q1" << endl;
    Matrix psi = this_q2^this_q1;
    psi.print();
    cout << "rho = density(psi)" << endl;
    Matrix rho = psi*psi.adjoint();
    rho.print();
    cout << "partial_2(rho)" << endl;
    Matrix ptr = rho.partialTrace(2);
    ptr.print();*/


    /*SQC_Circuit* blah = SQC_Circuit::LoadTFCFile("test2.tfc");
    blah->Print();

    SQC_Circuit result = SQC_Circuit::UniversalOptimize(*blah,LempelXSynthesis);
    result.Print();
    cout << "Before:" << endl;
    blah->PrintOperatorDistribution();
    cout << "After:" << endl;
    result.PrintOperatorDistribution();

    cout << "Matrix before:" << endl;
    Matrix matrix_before = TO_Maps::SQC_Circuit_to_Matrix(*blah);
    matrix_before.print();
    cout << "Matrix after:" << endl;
    Matrix matrix_after = TO_Maps::SQC_Circuit_to_Matrix(result);
    matrix_after.print();

    cout << "V(before, after) = " << VerifyOptimization2(*blah,result) << endl;*/

    /* blah = SQC_Circuit::LoadTFCFile("test1.tfc");
    SQC_Circuit* blee = SQC_Circuit::LoadTFCFile("test2.tfc");
    cout << "V(blah, blee) = " << VerifyOptimization(*blah,*blee) << endl;*/


    //(TO_Maps::SQC_Circuit_to_Matrix(result)).print();

    //warning("TODO: Debug step 4 of UniversalOptimize.");

    //LOut() << "matrix precision = " << g_matrix_precision << endl;

    /*SQC_Circuit blah;
    blah.Load("test.sqc");
    blah.Print();

    SQC_Circuit result = SQC_Circuit::UniversalOptimize(blah,LempelXSynthesis);
    result.simplify();
    result.Print();*/

    /*SQC_Circuit result = SQC_Circuit::convert_Cliff3s(blah);
    result.Print();
    result.simplify();
    result.Print();

    SQC_Circuit** Ps = NULL;
    SQC_Circuit** Hs = NULL;
    int N_Ps, N_Hs;
    SQC_Circuit::decompose_into_Hadamard_partitions(result,Hs,N_Hs,Ps,N_Ps);
    cout << "N_Ps = " << N_Ps << endl;
    cout << "N_Hs = " << N_Hs << endl;
    for(int i = 0; i < N_Hs; i++) {
        cout << "Hadamard Partition " << i << ":" << endl;
        Hs[i]->Print();
    }
    for(int i = 0; i < N_Ps; i++) {
        cout << "Partition " << i << ":" << endl;
        Ps[i]->Print();
    }

    SQC_Circuit* myL = NULL;
    SQC_Circuit* myPp = NULL;
    SQC_Circuit* myR = NULL;
    SQC_Circuit::Hadamards_to_Gadgets(*Ps[0],myL,myPp,myR);
    myPp->Print();
    SQC_Circuit* myCNOT = NULL;
    SQC_Circuit* myD3 = NULL;
    SQC_Circuit::decompose_C3_to_CNOT_D3(*myPp, myCNOT, myD3);
    myD3->simplify();
    myD3->Print();
    SQC_Circuit blugh = SQC_Circuit::optimize_D3(*myD3,LempelXSynthesis);
    blugh.Print();*/

    /*
    PhasePolynomial blee = TO_Maps::SQC_Circuit_to_PhasePolynomial(blah);
    blee.print();
    PhasePolynomial blon = FullDecoderWrapper(blee,LempelXSynthesis);
    blon.print();
    WeightedPolynomial blee_wp = TO_Maps::PhasePolynomial_to_WeightedPolynomial(blee);
    WeightedPolynomial blon_wp = TO_Maps::PhasePolynomial_to_WeightedPolynomial(blon);
    blee_wp.print();
    blon_wp.print();
    */

    if(argc>1) {
        string this_command = argv[1];
        if(!this_command.compare("optimize")&&(argc>2)) {
            string circuit_filename = argv[2];
            int option_maxRM = 5;
            bool option_fileoutput = false;
            string option_filename;
            bool option_default_algorithm = true;
            string option_algorithm;
            bool option_feedback = true;
            int option_n_rand = 1;
            int option_psrmc = 0;
            int option_col_order = 0;
            int option_n_orders = 1;
            for(int i = 3; i < argc; i++) {
                string this_option = argv[i];
                if((this_option[0]=='-')&&((i+1)<argc)) {
                    string this_value = argv[i+1];
                    char this_option_char = this_option[1];
                    switch(this_option_char) {
                        case 'r':
                            option_maxRM = atoi(this_value.c_str());
                            break;
                        case 'o':
                            option_fileoutput = true;
                            option_filename = this_value;
                            break;
                        case 'a':
                            option_default_algorithm = false;
                            option_algorithm = this_value;
                            g_algorithm = option_algorithm;
                            break;
                        case 'f':
                            if((!this_value.compare("true"))||(!this_value.compare("1"))) {
                                option_feedback = true;
                            } else if((!this_value.compare("false"))||(!this_value.compare("0"))) {
                                option_feedback = false;
                            }
                            break;
                        case 'd':
                            option_n_rand = atoi(this_value.c_str());
                            break;
                        case 'c':
                            g_csv_filename = this_value;
                            break;
                        case 'p':
                            option_psrmc = atoi(this_value.c_str());
                            break;
                        case 'q':
                            option_col_order = atoi(this_value.c_str());
                            break;
                        case 's':
                            option_n_orders = atoi(this_value.c_str());
                            break;
                    }
                    i++;
                }
            }
            g_lempel_feedback = option_feedback;
            g_Reed_Muller_max = option_maxRM;
            g_output_filename = option_filename;
            Signature this_sig;
            char first_char = circuit_filename.c_str()[0];
            if(first_char==':') {
                cout << "Structured circuit." << endl;
                this_sig = CircuitGenerator(circuit_filename.substr(1));
            } else {
                cout << "Circuit from file: " << circuit_filename << endl;
                this_sig = Signature::SigFromFile(circuit_filename.c_str());
            }
            GateStringSparse out(this_sig.get_n());
            cout << "Input signature polynomial:" << endl;
            this_sig.print();
            cout << endl << "Synthesis algorithm: ";
            if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY)) {
                cout << "Lempel least-greedy" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_LeastGreedy,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::REED_MULLER)) {
                cout << "Reed-Muller" << endl;
                GateStringSparse tempout = ReedMullerSynthesis2(this_sig);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY)) {
                cout << "Lempel greedy" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_Greedy,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::DAFT_GUESS)) {
                cout << "Daft guess" << endl;
                GateStringSparse tempout = GateSigInterface::SigToGSS(this_sig);
                if(!g_csv_filename.empty()) {
                    ofstream my_csv(g_csv_filename.c_str(),iostream::app);
                    int n = this_sig.get_n();
                    my_csv << SYNTHESIS_ALGORITHM_ID::DAFT_GUESS << "," << n << "," << g_random_circuit_seed << "," << tempout.weight(true) << endl;
                    my_csv.close();
                }
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM)) {
                cout << "Lempel random" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_Random,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::ALL_LEMPEL_SELECTORS)) {
                cout << "All Lempel Selectors" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,option_feedback,option_n_rand,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(option_default_algorithm||!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::ALL_LEMPEL)) {
                cout << "All Lempel Selectors w/ && w/o feedback" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,option_n_rand,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X)) {
                cout << "LempelX" << endl;
                GateStringSparse tempout = LempelXSynthesis(this_sig);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X_2)) {
                cout << "LempelX.2" << endl;
                GateStringSparse tempout = LempelXSynthesis(this_sig);
                out.assign(tempout);
            }
            result_analysis(this_sig,out);
            cout << "Output phase polynomial:" << endl;
            out.print();

            // PSRMC Test
            /*
            BMSparse out_daft = Interface_BMSGSS::GSSToBMS(GateSigInterface::SigToGSS(this_sig));
            BMSparse out_BMS = Interface_BMSGSS::GSSToBMS(out);
            //BMSparse myL, myE, myU, myK;
            //out_BMS.LEUKTrapezoidal(myL,myE,myU,myK,6);
            cout << "in BMS" << endl;
            out_daft.printFull();
            cout << "4x1x2x3 = " << ((out_daft.row(0).bitwise_AND(out_daft.row(1)).bitwise_AND(out_daft.row(2))).sum()%2) << endl;
            cout << "out BMS" << endl;
            out_BMS.printFull();
            cout << "4x1x2x3 = " << ((out_BMS.row(0).bitwise_AND(out_BMS.row(1)).bitwise_AND(out_BMS.row(2))).sum()%2) << endl;
            /*
            cout << "out U order = 0" << endl;
            myU.printFull();
            out_BMS.LEUKTrapezoidal(myL,myE,myU,myK,6,1);
            cout << "out U order = 1" << endl;
            myU.printFull();
            BMSparse myl, myp, myu;
            out_BMS.LPUTrapezoidal(myl,myp,myu,6);
            myu.printFull();
            */
        } else if(!this_command.compare("test")&&(argc>=3)) {
            string test_name = argv[2];
            if(!test_name.compare("PSRMC_1")) {
                test_PSRMC_1();
            } else if(!test_name.compare("Nullspace_t_vs_n")) {
                test_Nullspace_t_vs_n();
            } else if(!test_name.compare("LX_exec_times")) {
                test_LX_exec_times();
            }
        } else if(!this_command.compare("generate")&&(argc>=3)) {
            cout << "Procedurally generate a circuit." << endl;
            string circuit_name = argv[2];
            Signature this_sig = CircuitGenerator(circuit_name);
            cout << "Output signature:" << endl;
            this_sig.print();
            if(argc>=4) {
                string out_filename = argv[3];
                cout << "Saving to file: " << out_filename << "..." << endl;
                this_sig.save(out_filename.c_str());
                cout << "Saved." << endl;
            }
        } else if(!this_command.compare("list-algorithms")&&(argc>=2)) {
            cout << "ID\tCLA\t\tAlgorithm" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::DAFT_GUESS << "\t" << SYNTHESIS_ALGORITHM_TAG::DAFT_GUESS << "\t\t" << "Daft guess" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::REED_MULLER << "\t" << SYNTHESIS_ALGORITHM_TAG::REED_MULLER << "\t\t" << "Reed Muller decoder" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY << "\t\t" << "Lempel, least greedy order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY << "\t\t" << "Lempel, greedy order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM << "\t\t" << "Lempel, random order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY << " -f false\t" << "Lempel, least greedy order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY << " -f false\t" << "Lempel, greedy order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM << " -f false\t" << "Lempel, random order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_X << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_X << "\t\t" << "LempelX" << endl;
        } else if(!this_command.compare("UniversalT")&&(argc>=3)) {
            cout << "Find T-Count for Universal Circuit" << endl;
            LOut_Pad++;
            string circuit_filename = argv[2];
            bool option_maslov = 0;
            SQC_ToffoliNMode option_toffmode = SQC_TOFFOLI_N_MODE_TOFF3;
            SQC_HadamardMode option_hadmode = SQC_HADAMARD_MODE_MONTANARO;
            int option_hadmode_max_anc = -1;
            for(int i = 3; i < argc; i++) {
                string this_option = argv[i];
                if((this_option[0]=='-')&&((i+1)<argc)) {
                    string this_value = argv[i+1];
                    char this_option_char = this_option[1];
                    switch(this_option_char) {
                        case 'm':
                            option_maslov = 1;//atoi(this_value.c_str());
                            break;
                        case 't':
                            if(!this_value.compare("PSCJ")) option_toffmode = SQC_TOFFOLI_N_MODE_PSCJ;
                            else if(!this_value.compare("Jones")) option_toffmode = SQC_TOFFOLI_N_MODE_JONES;
                            else if(!this_value.compare("Toff3")) option_toffmode = SQC_TOFFOLI_N_MODE_TOFF3;
                            break;
                        case 'h':
                            if(!this_value.compare("Part")) option_hadmode = SQC_HADAMARD_MODE_PARTITION;
                            else if(!this_value.compare("Mont")) option_hadmode = SQC_HADAMARD_MODE_MONTANARO;
                            break;
                        case 'a':
                            cout << "option_hadmode_max_anc = " << this_value << endl;
                            option_hadmode_max_anc = atoi(this_value.c_str());
                            break;
                    }
                    i++;
                }
            }
            SQC_Circuit this_circuit;
            if(option_maslov) {
                this_circuit.LoadMaslovFile(circuit_filename.c_str());
            } else {
                this_circuit.Load(circuit_filename.c_str());
            }
            this_circuit.Print();
            double exec_time = 0;
            int daft_t_count = 0;
            int n_parts;
            int n_hadancs = 0;
            this_circuit.toffoli_n_mode = option_toffmode;
            this_circuit.hadamard_mode = option_hadmode;
            this_circuit.hadamard_mode_max_ancillas = option_hadmode_max_anc;
            int result = UniversalTCount(&this_circuit, &daft_t_count, &exec_time,&n_parts,&n_hadancs);
            cout << "Final T-count (LempelX) = " << result << endl;
            cout << "Final T-count (Daft) = " << daft_t_count << endl;
            cout << "Final T-count reduced by " << (daft_t_count-result) << endl;
            cout << "Total execution time = " << exec_time << "s" << endl;
            cout << "No. partitions = " << n_parts << endl;
            cout << "No. hadamard ancillas = " << n_hadancs << endl;
            this_circuit.Destruct();
        } else if(!this_command.compare("UniversalOptimize")&&(argc>=3)) {
            cout << "Find T-Count optimal circuit decomposition for {Clifford,T}" << endl;
            SQC_Circuit* this_circuit = NULL;
            TO_Decoder this_decoder = LempelXSynthesis;
            string circuit_filename = argv[2];
            LOut() << "Input filename = " << circuit_filename << endl;
            string this_input_filetype = "sqc";
            bool this_verify = 0;

            for(int i = 3; i < argc; i++) {
                string this_option = argv[i];
                if((this_option[0]=='-')&&((i+1)<argc)) {
                    string this_value = argv[i+1];
                    char this_option_char = this_option[1];
                    switch(this_option_char) {
                        case 'f':
                            // filetype
                            this_input_filetype = this_value;
                            break;
                        case 'o':
                            // output filename
                            g_output_filename = this_value;
                            break;
                        case 'a':
                            // algorithm
                            g_algorithm = this_value;
                            break;
                        case 'v':
                            // verify?
                            this_verify = atoi(this_value.c_str());
                            break;
                        case 'h':
                            //no. had. ancs
                            g_Hadamard_ancillas = atoi(this_value.c_str());
                            break;
                    }
                    i++;
                }
            }

            // Load circuit depending on filetype
            if(!this_input_filetype.compare("sqc")) {
                this_circuit = new SQC_Circuit();
                this_circuit->Load(circuit_filename.c_str());
            } else if(!this_input_filetype.compare("tfc")) {
                this_circuit = SQC_Circuit::LoadTFCFile(circuit_filename.c_str());
            }

            // Determine algorithm to be used
            if(!g_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X)) {
                this_decoder = LempelXSynthesis;
            } /*else if(!g_algorithm.compare()) {
                REST OF SYNTHESIS ALGORITHMS GO HERE
            }*/

            if(this_circuit) {
                LOut() << "Input circuit:" << endl;
                this_circuit->Print();
                SQC_Circuit result = SQC_Circuit::UniversalOptimize(*this_circuit,this_decoder);
                LOut() << "Output circuit:" << endl;
                result.Print();
                LOut() << "Gate distributions:" << endl;
                LOut() << "Input:" << endl;
                this_circuit->PrintOperatorDistribution();
                LOut() << "Output:" << endl;
                result.PrintOperatorDistribution();
                if(this_verify) {
                    VerifyOptimization2(*this_circuit,result);
                }
            }

            // Delete circuits
            if(this_circuit) {
                delete this_circuit;
                this_circuit = NULL;
            }
        }
    }
    return 0;
}





Signature CircuitGenerator_Toffhash(int N_hash) {
    int n = 2*N_hash+1;
    Signature out(n);
    for(int t = 0; t < N_hash; t++){
        out.set(1,2*t+2,2*t+3);
    }
    ostringstream temp_ss;
    temp_ss << N_hash;
    g_indvar_out = temp_ss.str();
    return out;
}

Signature CircuitGenerator_Toffoli(int N_toff) {
    int n = 3*N_toff;
    Signature out(n);
    for(int t = 0; t < N_toff; t++){
        out.set(3*t+1,3*t+2,3*t+3);
    }
    ostringstream temp_ss;
    temp_ss << N_toff;
    g_indvar_out = temp_ss.str();
    return out;
}

Signature CircuitGenerator(const string& inS) {
    Signature out;

    char split_char = '-';
    int split_pos = inS.find(split_char);
    string circuit_name = inS.substr(0,split_pos);
    string circuit_args = inS.substr(split_pos+1);

    if(!circuit_name.compare("Toffhash")) {
        int N_hash = atoi(circuit_args.c_str());
        cout << "Toffhash circuit; N_hash = " << N_hash << endl;
        out = CircuitGenerator_Toffhash(N_hash);
    } else if(!circuit_name.compare("Toffoli")) {
        int N_toff = atoi(circuit_args.c_str());
        cout << "Toffoli circuit; N_toff = " << N_toff << endl;
        out = CircuitGenerator_Toffoli(N_toff);
    } else if(!circuit_name.compare("RandomComplex")) {
        int n;
        int this_seed = 0;
        int arg_split_pos = circuit_args.find(split_char);
        if(arg_split_pos!=circuit_args.npos) {
            string seed_string = circuit_args.substr(arg_split_pos+1);
            this_seed = atoi(seed_string.c_str());
        }
        n = atoi(circuit_args.substr(0,arg_split_pos).c_str());

        cout << "Random complex circuit; n = " << n;
        if(this_seed>0) {
            cout << "; seed = " << this_seed;
            g_random_circuit_seed = this_seed;
        }
        cout << endl;
        out = CircuitGenerator_RandomComplex(n,this_seed);
    }

    return out;
}

Signature CircuitGenerator_RandomComplex(int n, int in_seed) {
    Signature out;

    GateStringSparse complex_GSS = GateStringSparse::randomUpTo3Qu(n, in_seed);
    out = GateSigInterface::expandGSSTerm(complex_GSS);

    ostringstream temp_ss;
    temp_ss << n;
    g_indvar_out = temp_ss.str();

    return out;
}





/*
BMSparse decompTest(const BMSparse& inB) {
    BMSparse out;



    return out;
}*/

void test_PSRMC_1() {
    int n = 10;
    int n_RM = 5;
    Signature mySig = GateSigInterface::expandGSSTerm(GateStringSparse::randomUpTo3Qu(n));
    cout << "Input signature" << endl;
    mySig.print();

    GateStringSparse result = LempelSynthesis(mySig,n_RM,LempelSelector_Random,true);
    cout << "Output phase polynomial" << endl;
    result.print();
    result.printString();

    BMSparse result_BMS = Interface_BMSGSS::GSSToBMS(result);
    cout << "Synthesis matrix" << endl;
    result_BMS.printFull();

    {
        BMSparse myL, myE, myU, myK;
        result_BMS.LEUKTrapezoidal(myL,myE,myU,myK,n_RM,1);
        cout << "leUk random" << endl;
        myU.printFull();
        cout << "width = " << trapezoid_width(myU,n_RM) << endl;
        cout << "Nullspace" << endl;
        BMSparse myNull = myU.nullspace();
        myNull.printFull();
        (myU*myNull.col(0)).printFull();
    }
    {
        BMSparse myL, myE, myU, myK;
        result_BMS.LEUKTrapezoidal(myL,myE,myU,myK,n_RM,0);
        cout << "leUk many 1's" << endl;
        myU.printFull();
        cout << "width = " << trapezoid_width(myU,n_RM) << endl;
        cout << "Nullspace" << endl;
        BMSparse myNull = myU.nullspace();
        myNull.printFull();
        (myU*myNull.col(0)).printFull();
    }

}


void test_Nullspace_t_vs_n() {
    string out_filename = "./test_Nullspace_t_vs_n.csv";
    int start_n = 6;
    int end_n = 30;
    double max_time = 5*60.0;
    int max_d = 1000;
    double P = 0.5;
    ofstream myfile(out_filename.c_str(),iostream::app);
    cout << "Test: Nullspace_t_vs_n begin." << endl;
    for(int n = start_n; n <= end_n; n++) {
        cout << "n = " << n << endl;
        double total_exec_time = 0.0;
        int d = 0;
        while((total_exec_time<max_time)&&(d<max_d)) {
            cout << "Iteration d = " << d << " begin." << endl;
            BMSparse this_BMS = BMSparse::random(n,n,P);
            clock_t tic = clock();
            BMSparse this_nullspace = this_BMS.nullspace();
            clock_t toc = clock();
            double exec_time = (double)(toc-tic)/CLOCKS_PER_SEC;
            cout << "Executed in " << exec_time << "s" << endl;
            total_exec_time += exec_time;
            d++;
        }
        double mean_exec_time = total_exec_time/d;
        // Line format for csv file: <n>,<t>,<N_samples>
        myfile << n << "," << mean_exec_time << "," << d << endl;
    }
    cout << "Test: Nullspace_t_vs_n end." << endl;
    myfile.close();
}





int UniversalTCount(SQC_Circuit* inC, int* out_daft, double* out_exec_time, int* n_parts, int* n_hadancs) {
    //cout << "C" << endl;
    int out = 0;
    if(out_daft) *out_daft = 0;
    //inC->Print();
    //cout << endl;
    int tic = clock();
    bool exit = true;
    Signature this_sig;
    int this_part = 1;
    LOut() << "UniversalTCount: Algorithm begin." << endl;
    LOut_Pad++;
    inC->ConvertFromToffoli();
    inC->RemoveExternalHadamards();
    //LOut() << "Before" << endl;
    //inC->Print();
    {
        bool this_exit = 0;
        while(!this_exit) {
            //cout << "A" << endl;
            this_exit = 0;
            this_exit += inC->CancelAdjacentTs();
            this_exit += inC->CancelAdjacentHadamards();
            this_exit = !this_exit;
        }
    }

    //LOut() << "After" << endl;
    //inC->Print();
    LOut() << "UniversalTCount: Converted to elementary gate set." << endl;
    //inC->Print(&cout,0,65);
    LOut() << endl;
    //return out;

    do {
        LOut() << "Partition " << this_part << endl;
        //inC->Print();
        LOut_Pad++;
        int init_gates = inC->m;
        LOut() << "Remaining gates = " << init_gates << endl;
        //inC->Print();
        exit = !inC->NextSignature(this_sig);
        int reduces_gates = init_gates-inC->m;
        LOut() << "Converted" << endl;
        //inC->Print();
        //LOut() << "m = " << inC->m << endl;
        this_sig.print();
        GateStringSparse result = LempelXSynthesis(this_sig);
        GateStringSparse result_daft = GateSigInterface::SigToGSS(this_sig);
        int this_t_count = result.weight(true);
        int this_t_count_daft = result_daft.weight(true);
        out += this_t_count;
        *out_daft += this_t_count_daft;
        LOut_Pad--;
        LOut() << "T-count for partition " << this_part << " = " << this_t_count << endl;
        LOut() << "Total T-count so far = " << out << endl;
        LOut() << "T-gate per gate for this partition " << this_part << " = " << (double)this_t_count/(double)reduces_gates << endl;
        LOut() << "Hadamard ancillas for partition " << this_part << " = " << inC->hadamard_mode_max_ancillas << endl;
        LOut_Pad++;
        LOut() << "..." << endl;
        LOut_Pad--;
        /*
        if(!g_algorithm.compare(DAFT_GUESS)) {

        } else if(!g_algorithm.compare(LEMPEL_GREEDY)) {

        } else if(!g_algorithm.compare(LEMPEL_LEAST_GREEDY)) {

        } else if(!g_algorithm.compare(LEMPEL_RANDOM)) {

        } else if(!g_algorithm.compare(LEMPEL_X)) {

        } else if(!g_algorithm.compare(REED_MULLER)) {

        } else if(!g_algorithm.compare(ALL_LEMPEL_SELECTORS)) {

        } else if(!g_algorithm.compare(ALL_LEMPEL)) {

        }
        */
        this_part++;
    } while(!exit);
    if(n_parts) *n_parts = (this_part-1);
    if(n_hadancs) *n_hadancs = inC->hadamard_mode_max_ancillas;
    int toc = clock();
    double exec_time = ((double)toc-(double)tic)/CLOCKS_PER_SEC;
    if(out_exec_time) *out_exec_time = exec_time;
    LOut_Pad--;
    LOut() << "UniversalTCount: Algorithm end." << endl << endl;
    return out;
}

void test_LX_exec_times() {
    cout << "LempelX execution times." << endl;
    int d = 9;
    int n_start = 6;
    int n_end = 16;
    for(int i = 0; i < d; i++) {
        cout << "i = " << i << " of " << d << endl;
        for(int n = n_start; n <= n_end; n++) {
            cout << "n = " << n << endl;
            Signature this_sig = CircuitGenerator_RandomComplex(n);
            this_sig.print();

            clock_t tic = clock();
            GateStringSparse result = LempelXSynthesis(this_sig);
            clock_t toc = clock();
            double exec_time = (double)(toc-tic)/(double)CLOCKS_PER_SEC;
            ofstream outfile("LX_exec_times.csv", iostream::app);
            outfile << n << "," << exec_time << endl;
            outfile.close();
        }
    }
}
